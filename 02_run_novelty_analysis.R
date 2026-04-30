##################################################################
####   STAGE 2 — RUN NOVELTY ANALYSIS                         ####
####                                                          ####
####   Loops over all future time steps defined in config.R.  ####
####   For each step: runs inclusion test (hypervolume) and   ####
####   optionally ExDet. Saves rasters and a summary table.   ####
####                                                          ####
####   Run 01_build_hypervolumes.R first.                     ####
##################################################################

library(terra)
library(hypervolume)

source("00_config.R")
source("functions/data_utils.R")
source("functions/spatial_grid.R")
source("functions/hypervolume_fns.R")
source("functions/exdet_fns.R")

init_output_dirs(output_dir)

# ── Load pre-computed hypervolumes & scaling params ───────────
message("Loading hypervolumes …")
hvs <- load_hypervolumes(
  gcm_name      = gcm_name,
  include_space = include_space,
  time_col      = time_col,
  time_levels   = time_levels,
  output_dir    = output_dir
)

scaling_params <- readRDS(
  file.path(output_dir, paste0("scaling_params_", gcm_name, ".rds"))
)

space_params <- NULL
if (include_space) {
  space_params <- readRDS(
    file.path(output_dir, paste0("space_params_", gcm_name, ".rds"))
  )
}

# ── Load and prepare historical data (for ExDet) ─────────────
hist_data      <- NULL
time_col_exdet <- NULL
if (run_exdet) {
  message("Loading historical data for ExDet …")
  hist_data <- readRDS(hist_file)
  hist_data <- scale_data(hist_data, scaling_params, env_vars)
  names(hist_data)[names(hist_data) == lon_col] <- "x"
  names(hist_data)[names(hist_data) == lat_col] <- "y"
  if (!is.null(time_col) && time_col %in% names(hist_data)) {
    names(hist_data)[names(hist_data) == time_col] <- "time_group"
    time_col_exdet <- "time_group"
  }
}

# ── Resolve spatial grid ──────────────────────────────────────
coast_data <- NULL
cells      <- NULL
if (!is.null(spatial_grid_type)) {
  if (file.exists(coast_file)) coast_data <- readRDS(coast_file)

  if (spatial_grid_type == "auto" && !is.null(coast_data) && run_exdet) {
    hist_data <- ensure_dist_coast(hist_data, coast_data, "x", "y")
  }

  cells <- resolve_spatial_grid(
    spatial_grid_type = spatial_grid_type,
    lat_breaks  = if (exists("lat_breaks"))  lat_breaks  else NULL,
    dist_breaks = if (exists("dist_breaks")) dist_breaks else NULL,
    coast_file  = coast_file,
    region_file = if (exists("region_file")) region_file else NULL
  )
  message("Spatial grid: ", length(cells), " cell(s)")
}

# ── Bundle shared state passed to each time-step worker ───────
ctx <- list(
  gcm_name              = gcm_name,
  env_vars              = env_vars,
  include_space         = include_space,
  time_col_exdet        = time_col_exdet,
  output_dir            = output_dir,
  hvs                   = hvs,
  scaling_params        = scaling_params,
  space_params          = space_params,
  hist_data             = hist_data,
  coast_data            = coast_data,
  cells                 = cells,
  run_exdet             = run_exdet,
  save_inclusion_rasters = save_inclusion_rasters,
  save_exdet_rasters    = save_exdet_rasters,
  future_timesteps      = future_timesteps
)

# ── Per-time-step worker function ─────────────────────────────
run_one_step <- function(step_idx, ctx) {
  row      <- ctx$future_timesteps[step_idx, ]
  year_i   <- row$year
  time_lbl <- if ("time_label" %in% names(row)) row$time_label else NA
  rpath    <- row$raster_path

  step_label <- if (is.na(time_lbl)) {
    as.character(year_i)
  } else {
    paste0("y", year_i, "_t", time_lbl)
  }
  n_steps <- nrow(ctx$future_timesteps)
  message(sprintf("  [%d/%d] %s", step_idx, n_steps, step_label))

  if (!file.exists(rpath)) {
    warning("Raster file not found, skipping: ", rpath)
    return(NULL)
  }

  rstack_env <- load_and_scale_future(rpath, ctx$env_vars, ctx$scaling_params)

  rstack_space <- NULL
  if (ctx$include_space && !is.null(ctx$space_params)) {
    rstack_space <- add_space_layers(
      rstack_env, ctx$space_params$lon, ctx$space_params$lat
    )
  }

  # ── Inclusion tests ──────────────────────────────────────────
  incl <- run_inclusion_all_scales(
    rstack_env   = rstack_env,
    rstack_space = rstack_space,
    hvs          = ctx$hvs,
    time_label   = if (!is.na(time_lbl)) time_lbl else NULL
  )

  props <- lapply(
    incl,
    function(r) if (!is.null(r)) novelty_proportion(r) else NA_real_
  )

  if (ctx$save_inclusion_rasters) {
    incl_dir <- file.path(ctx$output_dir, "inclusion_rasters", ctx$gcm_name)
    dir.create(incl_dir, recursive = TRUE, showWarnings = FALSE)
    for (sc in names(incl)) {
      if (!is.null(incl[[sc]])) {
        terra::writeRaster(
          incl[[sc]],
          file.path(incl_dir, paste0("incl_", sc, "_", step_label, ".tif")),
          overwrite = TRUE
        )
      }
    }
  }

  # ── ExDet ────────────────────────────────────────────────────
  exdet_stats <- list()
  if (ctx$run_exdet && !is.null(ctx$hist_data)) {
    futr_df <- raster_to_df(rstack_env, ctx$env_vars)

    if (!is.null(ctx$cells) && !is.null(ctx$coast_data)) {
      futr_df <- ensure_dist_coast(futr_df, ctx$coast_data, "x", "y")
    }

    exdet_result <- run_exdet_all_scales(
      hist_df      = ctx$hist_data,
      futr_df      = futr_df,
      covars       = ctx$env_vars,
      time_col     = ctx$time_col_exdet,
      time_label   = if (!is.na(time_lbl)) time_lbl else NULL,
      cells        = ctx$cells,
      save_rasters = ctx$save_exdet_rasters
    )

    if (ctx$save_exdet_rasters) {
      exdet_dir <- file.path(ctx$output_dir, "exdet_rasters", ctx$gcm_name)
      dir.create(exdet_dir, recursive = TRUE, showWarnings = FALSE)
      for (sc in c("atap", "stap")) {
        r <- exdet_result[[sc]]$raster
        if (!is.null(r)) {
          terra::writeRaster(
            r,
            file.path(exdet_dir,
                      paste0("exdet_", sc, "_", step_label, ".tif")),
            overwrite = TRUE
          )
        }
      }
      for (sc in c("atsp", "stsp")) {
        r_mosaic <- mosaic_cell_rasters(exdet_result[[sc]])
        if (!is.null(r_mosaic)) {
          terra::writeRaster(
            r_mosaic,
            file.path(exdet_dir,
                      paste0("exdet_", sc, "_", step_label, ".tif")),
            overwrite = TRUE
          )
        }
      }
    }

    for (sc in c("atap", "stap")) {
      s <- exdet_result[[sc]]$stats
      if (!is.null(s)) {
        names(s) <- paste0("exdet_", sc, "_", names(s))
        exdet_stats <- c(exdet_stats, as.list(s))
      }
    }
    for (sc in c("atsp", "stsp")) {
      s <- summarize_spatial_exdet(exdet_result[[sc]])
      if (!is.null(s)) {
        names(s) <- paste0("exdet_", sc, "_", names(s))
        exdet_stats <- c(exdet_stats, as.list(s))
      }
    }
  }

  # ── Assemble summary row ──────────────────────────────────────
  base_row <- data.frame(
    gcm        = ctx$gcm_name,
    year       = year_i,
    time_label = if (is.na(time_lbl)) NA_character_ else as.character(time_lbl),
    prop_novel_atap = unlist(props$atap),
    prop_novel_stap = unlist(props$stap),
    prop_novel_atsp = unlist(props$atsp),
    prop_novel_stsp = unlist(props$stsp),
    stringsAsFactors = FALSE
  )
  if (length(exdet_stats) > 0)
    base_row <- cbind(base_row, as.data.frame(exdet_stats))

  base_row
}

# ── Sequential or parallel execution ─────────────────────────
n_steps <- nrow(future_timesteps)
message(sprintf("\nProcessing %d future time steps for GCM '%s' …\n",
                n_steps, gcm_name))

if (n_cores > 1) {
  if (!requireNamespace("parallel", quietly = TRUE))
    stop("n_cores > 1 requires the 'parallel' package")
  cl <- parallel::makeCluster(n_cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterEvalQ(cl, {
    library(terra); library(hypervolume)
    source("functions/data_utils.R")
    source("functions/spatial_grid.R")
    source("functions/hypervolume_fns.R")
    source("functions/exdet_fns.R")
  })
  summary_rows <- parallel::parLapply(
    cl, seq_len(n_steps), run_one_step, ctx = ctx
  )
} else {
  summary_rows <- lapply(seq_len(n_steps), run_one_step, ctx = ctx)
}

# ── Combine and save ──────────────────────────────────────────
summary_df <- do.call(rbind, Filter(Negate(is.null), summary_rows))

if (save_summary_table) {
  out_path <- file.path(
    output_dir, "summaries",
    paste0("novelty_summary_", gcm_name, ".rds")
  )
  saveRDS(summary_df, out_path)
  message("\nSummary saved to: ", out_path)
}

message("Stage 2 complete.")
