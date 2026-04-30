##################################################################
####   DATA UTILITIES                                          ####
####   Scaling, path helpers, raster loading                  ####
##################################################################

# terra is used via terra:: prefix throughout; run install_packages.R if needed.

# ── Compute and store historical scaling parameters ───────────
#' @param hist_data  data.frame with columns matching env_vars
#' @param env_vars   character vector of variable names
#' @return list(means, SDs) — used to scale all future data
compute_scaling_params <- function(hist_data, env_vars) {
  list(
    means = colMeans(hist_data[, env_vars, drop = FALSE], na.rm = TRUE),
    SDs   = apply(hist_data[, env_vars, drop = FALSE], 2, sd, na.rm = TRUE)
  )
}

# ── Scale a data.frame using pre-computed params ──────────────
#' @param data      data.frame (or subset) with env_var columns
#' @param params    output of compute_scaling_params()
#' @param env_vars  variables to scale (others left untouched)
#' @return data.frame with env_vars standardised
scale_data <- function(data, params, env_vars) {
  scaled <- scale(
    data[, env_vars, drop = FALSE],
    center = params$means,
    scale  = params$SDs
  )
  out <- data
  out[, env_vars] <- as.data.frame(scaled)
  out
}

# ── Scale a SpatRaster stack using pre-computed params ────────
#' @param rstack    SpatRaster with layers named as env_vars
#' @param params    output of compute_scaling_params()
#' @param env_vars  variables to scale (must match layer names)
#' @return SpatRaster with scaled layers
scale_raster <- function(rstack, params, env_vars) {
  for (v in env_vars) {
    rstack[[v]] <- (rstack[[v]] - params$means[v]) / params$SDs[v]
  }
  rstack
}

# ── Load a future raster stack and standardise ────────────────
#' @param path      file path to raster stack (.tif, .grd, etc.)
#' @param env_vars  expected layer names
#' @param params    scaling parameters
#' @return scaled SpatRaster
load_and_scale_future <- function(path, env_vars, params) {
  r <- terra::rast(path)
  missing <- setdiff(env_vars, names(r))
  if (length(missing) > 0) {
    stop("Future raster at ", path,
         " is missing layers: ", paste(missing, collapse = ", "))
  }
  r <- r[[env_vars]]           # reorder to match env_vars
  scale_raster(r, params, env_vars)
}

# ── Raster to data.frame (including coordinates) ─────────────
raster_to_df <- function(rstack, env_vars) {
  df <- as.data.frame(rstack, xy = TRUE, na.rm = TRUE)
  names(df)[1:2] <- c("x", "y")
  df
}

# ── Build future_timesteps table from a directory ─────────────
#' Convenience helper: generates the future_timesteps data.frame
#' from a directory of raster stacks following a naming pattern.
#'
#' @param dir       directory containing raster files
#' @param pattern   regex to match filenames; must capture named groups:
#'                  (?P<gcm>...), (?P<time>...), (?P<year>...)
#'                  OR supply a format_fn instead.
#' @param format_fn function(gcm, time, year) → filename (without dir)
#' @param gcm_name  GCM label
#' @param time_labels vector of time labels (e.g. 1:12)
#' @param years     vector of years (e.g. 2010:2100)
#' @param ext       file extension (default ".tif")
#' @return data.frame with columns year, time_label, raster_path
build_future_timesteps <- function(dir, format_fn, gcm_name,
                                   time_labels, years, ext = ".tif") {
  grid <- expand.grid(time_label = time_labels, year = years,
                      stringsAsFactors = FALSE)
  grid$raster_path <- file.path(
    dir,
    mapply(format_fn, gcm_name, grid$time_label, grid$year,
           USE.NAMES = FALSE)
  )
  grid
}

# ── Split a data.frame into per-cell subsets ─────────────────
#' Used by both exdet_fns.R and 02_run_novelty_analysis.R.
#' @param data   data.frame with x, y (and dist_coast for "auto" grids)
#' @param cells  output of resolve_spatial_grid() — named list of filter fns
#' @return named list of data.frame subsets, one per non-empty cell
split_by_cell <- function(data, cells) {
  out <- lapply(cells, function(fn) data[fn(data), , drop = FALSE])
  Filter(function(d) nrow(d) > 0, out)
}

# ── Subsample rows of a data.frame ────────────────────────────
subsample_data <- function(data, fraction, seed = 123) {
  if (fraction >= 1) return(data)
  set.seed(seed)
  data[sample(nrow(data), floor(fraction * nrow(data))), ]
}

# ── Ensure output sub-directories exist ───────────────────────
init_output_dirs <- function(output_dir) {
  dirs <- c(
    file.path(output_dir, "hypervolumes"),
    file.path(output_dir, "inclusion_rasters"),
    file.path(output_dir, "exdet_rasters"),
    file.path(output_dir, "summaries")
  )
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
}

# ── terra → raster conversion (for packages needing raster) ──
#' Some packages (dsmextra) still require {raster} objects.
#' This converts a SpatRaster to RasterStack without loading {raster}
#' as a direct dependency in the main workflow.
spatraster_to_rasterstack <- function(spatraster) {
  if (!requireNamespace("raster", quietly = TRUE))
    stop("Package 'raster' needed for ExDet; install it with install.packages('raster')")
  raster::stack(raster::raster(spatraster))   # single-layer fast path
  # multi-layer:
  raster::stack(lapply(names(spatraster), function(nm) raster::raster(spatraster[[nm]])))
}
