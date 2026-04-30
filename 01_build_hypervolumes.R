##################################################################
####   STAGE 1 — BUILD HISTORICAL HYPERVOLUMES                ####
####                                                          ####
####   Run once per GCM. Saves hypervolume .rds files to     ####
####   outputs/hypervolumes/                                  ####
##################################################################
#
# The historical hypervolumes form the baseline to determine future novelty
#
# Before running:
#   1. Edit config.R (gcm_name, env_vars, hist_file, time_col, …)
#   2. Ensure your historical data file exists at hist_file
#
# Historical data format:
#   An .rds file containing a data.frame with columns:
#     - one column per env_var  (required)
#     - lon, lat                (required if include_space = TRUE)
#     - time_col value          (required if time_col is not NULL)
#   One row per observation (e.g. one row per grid-cell per month per year).
#   No subsampling needed here — that is handled internally.

library(terra)
library(hypervolume)

source("00_config.R")
source("functions/data_utils.R")
source("functions/hypervolume_fns.R")
source("functions/spatial_grid.R")

init_output_dirs(output_dir)

# ── Load historical data ──────────────────────────────────────
message("Loading historical data: ", hist_file)
hist_data <- readRDS(hist_file)

# Basic check
missing_cols <- setdiff(env_vars, names(hist_data))
if (length(missing_cols) > 0)
  stop("Historical data missing columns: ", paste(missing_cols, collapse = ", "))

if (include_space) {
  if (!lon_col %in% names(hist_data))
    stop("lon_col '", lon_col, "' not found in historical data")
  if (!lat_col %in% names(hist_data))
    stop("lat_col '", lat_col, "' not found in historical data")
}

# ── Compute scaling parameters (saved for use in Stage 2) ─────
message("Computing scaling parameters …")
scaling_params <- compute_scaling_params(hist_data, env_vars)
saveRDS(scaling_params,
        file.path(output_dir, paste0("scaling_params_", gcm_name, ".rds")))

# ── Scale historical data ─────────────────────────────────────
hist_scaled <- scale_data(hist_data, scaling_params, env_vars)

# Scale lon and lat too (needed for "with-space" hypervolumes)
if (include_space) {
  lon_params <- list(mean = mean(hist_data[[lon_col]]),
                     sd   = sd(hist_data[[lon_col]]))
  lat_params <- list(mean = mean(hist_data[[lat_col]]),
                     sd   = sd(hist_data[[lat_col]]))
  hist_scaled[[lon_col]] <- (hist_data[[lon_col]] - lon_params$mean) / lon_params$sd
  hist_scaled[[lat_col]] <- (hist_data[[lat_col]] - lat_params$mean) / lat_params$sd
  # Save lon/lat scaling for use in Stage 2
  saveRDS(list(lon = lon_params, lat = lat_params),
          file.path(output_dir, paste0("space_params_", gcm_name, ".rds")))
}

# ── Build hypervolumes ────────────────────────────────────────
message("\nBuilding hypervolumes for GCM: ", gcm_name)
message("  env_vars      : ", paste(env_vars, collapse = ", "))
message("  include_space : ", include_space)
message("  time grouping : ",
        if (is.null(time_col)) "none (annual only)" else
          paste0(time_col, " — ", length(time_levels), " levels"))
message("  sample_frac   : ", hv_sample_fraction,
        "  (~", round(nrow(hist_scaled) * hv_sample_fraction), " rows)")

build_all_hypervolumes(
  hist_data     = hist_scaled,
  env_vars      = env_vars,
  lon_col       = lon_col,
  lat_col       = lat_col,
  include_space = include_space,
  time_col      = time_col,
  time_levels   = time_levels,
  method        = hv_method,
  sample_frac   = hv_sample_fraction,
  seed          = hv_seed,
  gcm_name      = gcm_name,
  output_dir    = output_dir
)

message("\nStage 1 complete. Hypervolumes saved to: ",
        file.path(output_dir, "hypervolumes"))
