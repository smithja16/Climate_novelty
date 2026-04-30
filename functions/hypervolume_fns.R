##################################################################
####   HYPERVOLUME FUNCTIONS                                   ####
####   Build hypervolumes and run inclusion tests              ####
##################################################################
#
# Key design: env_vars, time grouping, and spatial dimensions are
# always passed explicitly — nothing is hard-coded.

library(hypervolume)
# terra is used via terra:: prefix; run install_packages.R if needed.

# ── Build a single hypervolume ────────────────────────────────
#' @param data        data.frame; columns must include all vars in 'dims'
#' @param dims        character vector of column names to use as dimensions
#' @param method      "svm" | "box" | "gaussian"
#' @param sample_frac fraction of rows to use (1 = all)
#' @param seed        RNG seed for subsampling
#' @return hypervolume object
build_hv <- function(data, dims, method = "svm",
                     sample_frac = 1, seed = 123) {
  d <- data[, dims, drop = FALSE]
  if (sample_frac < 1) {
    set.seed(seed)
    d <- d[sample(nrow(d), floor(sample_frac * nrow(d))), ]
  }
  hypervolume::hypervolume(data = d, method = method, verbose = FALSE)
}

# ── Build all hypervolumes for one GCM ────────────────────────
#' Builds the four sets of hypervolumes needed for the four novelty scales.
#' Saves each as an .rds file in output_dir/hypervolumes/ and returns them
#' in a named list (invisible, since files may be large).
#'
#' @param hist_data     scaled historical data.frame
#' @param env_vars      env variable column names
#' @param lon_col/lat_col column names for spatial coordinates
#' @param include_space whether to build space-inclusive hypervolumes
#' @param time_col      column for temporal grouping (NULL = annual only)
#' @param time_levels   values of time_col to iterate over
#' @param method, sample_frac, seed  passed to build_hv()
#' @param gcm_name      label for filenames
#' @param output_dir    parent output directory
#' @return named list with elements atap, atsp (if include_space),
#'         stap (list by time level), stsp (list by time level, if include_space)
build_all_hypervolumes <- function(hist_data,
                                   env_vars,
                                   lon_col       = "lon",
                                   lat_col       = "lat",
                                   include_space = TRUE,
                                   time_col      = NULL,
                                   time_levels   = NULL,
                                   method        = "svm",
                                   sample_frac   = 0.05,
                                   seed          = 123,
                                   gcm_name      = "gcm",
                                   output_dir    = "outputs") {

  hv_dir <- file.path(output_dir, "hypervolumes")
  dir.create(hv_dir, recursive = TRUE, showWarnings = FALSE)

  dims_env   <- env_vars
  dims_space <- c(env_vars, lon_col, lat_col)

  hvs <- list()

  # ATAP — any time, any place (no spatial or temporal subsetting)
  message("Building ATAP hypervolume …")
  hvs$atap <- build_hv(hist_data, dims_env, method, sample_frac, seed)
  saveRDS(hvs$atap, file.path(hv_dir, paste0("hv_atap_", gcm_name, ".rds")))

  # ATSP — any time, same place (spatial dims included)
  if (include_space) {
    message("Building ATSP (with-space) hypervolume …")
    hvs$atsp <- build_hv(hist_data, dims_space, method, sample_frac, seed)
    saveRDS(hvs$atsp, file.path(hv_dir, paste0("hv_atsp_", gcm_name, ".rds")))
  }

  # STAP and STSP — one hypervolume per time level
  if (!is.null(time_col) && !is.null(time_levels)) {
    hvs$stap <- list()
    if (include_space) hvs$stsp <- list()

    for (tl in time_levels) {
      label <- as.character(tl)
      subset_t <- hist_data[hist_data[[time_col]] == tl, ]

      message(sprintf("  STAP time=%s …", label))
      hv_stap <- build_hv(subset_t, dims_env, method, sample_frac, seed)
      hvs$stap[[label]] <- hv_stap
      saveRDS(hv_stap, file.path(hv_dir,
                paste0("hv_stap_t", label, "_", gcm_name, ".rds")))

      if (include_space) {
        message(sprintf("  STSP time=%s …", label))
        hv_stsp <- build_hv(subset_t, dims_space, method, sample_frac, seed)
        hvs$stsp[[label]] <- hv_stsp
        saveRDS(hv_stsp, file.path(hv_dir,
                  paste0("hv_stsp_t", label, "_", gcm_name, ".rds")))
      }
    }
  }

  invisible(hvs)
}

# ── Load pre-built hypervolumes ───────────────────────────────
#' Loads the .rds files written by build_all_hypervolumes().
#' Call this at the top of 02_run_novelty_analysis.R.
load_hypervolumes <- function(gcm_name, include_space, time_col,
                               time_levels, output_dir) {
  hv_dir <- file.path(output_dir, "hypervolumes")
  hvs <- list()

  hvs$atap <- readRDS(file.path(hv_dir, paste0("hv_atap_", gcm_name, ".rds")))

  if (include_space)
    hvs$atsp <- readRDS(file.path(hv_dir, paste0("hv_atsp_", gcm_name, ".rds")))

  if (!is.null(time_col) && !is.null(time_levels)) {
    hvs$stap <- list()
    if (include_space) hvs$stsp <- list()
    for (tl in time_levels) {
      label <- as.character(tl)
      hvs$stap[[label]] <- readRDS(file.path(hv_dir,
        paste0("hv_stap_t", label, "_", gcm_name, ".rds")))
      if (include_space)
        hvs$stsp[[label]] <- readRDS(file.path(hv_dir,
          paste0("hv_stsp_t", label, "_", gcm_name, ".rds")))
    }
  }
  hvs
}

# ── Run inclusion test for a raster stack ────────────────────
#' Tests whether each pixel in a (scaled) future raster stack falls
#' inside the hypervolume.
#'
#' @param hv        hypervolume object
#' @param rstack    SpatRaster with layers named as the hypervolume's dimensions
#' @param hv_dims   dimension names the hypervolume was built with
#' @return SpatRaster: 1 = included (analog), 0 = excluded (novel), NA = land
run_inclusion <- function(hv, rstack, hv_dims) {
  # hypervolume_project() expects layers in the same order as hv was built
  rstack <- rstack[[hv_dims]]

  # hypervolume >= 3.1 accepts SpatRaster directly; fall back to raster if needed
  result <- tryCatch(
    hypervolume::hypervolume_project(hv, rasters = rstack,
                                     type = "inclusion",
                                     fast.or.accurate = "accurate"),
    error = function(e) {
      # fall back: convert to raster::RasterStack
      rs <- raster::stack(lapply(hv_dims,
              function(nm) raster::raster(rstack[[nm]])))
      hypervolume::hypervolume_project(hv, rasters = rs,
                                       type = "inclusion",
                                       fast.or.accurate = "accurate")
    }
  )
  # normalise back to SpatRaster
  if (!inherits(result, "SpatRaster"))
    result <- terra::rast(result)
  result
}

# ── Run all four scales for one future time step ──────────────
#' Returns a named list of inclusion SpatRasters: atap, atsp, stap, stsp.
#' Scales that cannot be computed (e.g. no spatial dims, no time grouping)
#' are NULL in the returned list.
#'
#' @param rstack_env   scaled SpatRaster with env_vars layers
#' @param rstack_space scaled SpatRaster with env_vars + lon + lat layers
#' @param hvs          output of load_hypervolumes()
#' @param time_label   value of the current time level (e.g. month number)
run_inclusion_all_scales <- function(rstack_env, rstack_space = NULL,
                                     hvs, time_label = NULL) {
  out <- list(atap = NULL, atsp = NULL, stap = NULL, stsp = NULL)

  # Extract dimension names stored in the hypervolume object
  hv_dims <- function(hv) colnames(hv@Data)

  out$atap <- run_inclusion(hvs$atap, rstack_env, hv_dims(hvs$atap))

  if (!is.null(hvs$atsp) && !is.null(rstack_space))
    out$atsp <- run_inclusion(hvs$atsp, rstack_space, hv_dims(hvs$atsp))

  if (!is.null(time_label) && !is.null(hvs$stap)) {
    label <- as.character(time_label)
    if (!is.null(hvs$stap[[label]]))
      out$stap <- run_inclusion(hvs$stap[[label]], rstack_env,
                                hv_dims(hvs$stap[[label]]))
    if (!is.null(hvs$stsp) && !is.null(hvs$stsp[[label]]) && !is.null(rstack_space))
      out$stsp <- run_inclusion(hvs$stsp[[label]], rstack_space,
                                hv_dims(hvs$stsp[[label]]))
  }

  out
}

# ── Proportion of pixels that are novel ──────────────────────
novelty_proportion <- function(inclusion_raster) {
  v <- terra::values(inclusion_raster, na.rm = FALSE)
  valid <- !is.na(v)
  if (sum(valid) == 0) return(NA_real_)
  sum(v[valid] == 0) / sum(valid)
}

# ── Add lon/lat layers to a scaled env raster stack ──────────
add_space_layers <- function(rstack, lon_params, lat_params) {
  r_lon <- rstack[[1]]; terra::values(r_lon) <- terra::crds(rstack)[, 1]
  r_lat <- rstack[[1]]; terra::values(r_lat) <- terra::crds(rstack)[, 2]
  r_lon <- (r_lon - lon_params$mean) / lon_params$sd
  r_lat <- (r_lat - lat_params$mean) / lat_params$sd
  names(r_lon) <- "lon"; names(r_lat) <- "lat"
  c(rstack, r_lon, r_lat)
}
