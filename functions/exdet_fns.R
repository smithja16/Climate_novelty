##################################################################
####   EXDET (EXTRAPOLATION DETECTION) FUNCTIONS               ####
####   Uses dsmextra::compute_extrapolation()                  ####
##################################################################
#
# ExDet quantifies three types of conditions in the future data:
#   Analogue     (ExDet in 0–1): within historical range and correlations
#   Univariate   (ExDet < 0)   : at least one variable outside historical range
#   Combinatorial(ExDet > 1)   : within variable ranges but novel combinations
#
# "Same place" analyses split the domain into spatial cells and run
# compute_extrapolation() separately for each cell, then mosaic results.

# terra is used via terra:: prefix; run install_packages.R if needed.
# split_by_cell() is defined in functions/data_utils.R — source that file first.

# ── Run ExDet for a single (historical, future) pair ─────────
#' @param hist_df   data.frame with x, y, and env_var columns (reference)
#' @param futr_df   data.frame with x, y, and env_var columns (target)
#' @param covars    character vector — variable names to include
#' @param crs_str   CRS string for sp::CRS() (default: WGS84 geographic)
#' @return list with: summary_stats (data.frame), exdet_raster (SpatRaster)
run_exdet_single <- function(hist_df, futr_df, covars,
                              crs_str = "+proj=longlat +datum=WGS84") {
  if (!requireNamespace("dsmextra", quietly = TRUE))
    stop("Package 'dsmextra' required; install from CRAN or GitHub.")
  if (!requireNamespace("sp", quietly = TRUE))
    stop("Package 'sp' required for dsmextra CRS handling.")

  roms_crs <- sp::CRS(crs_str)

  extrap <- tryCatch(
    dsmextra::compute_extrapolation(
      samples          = hist_df,
      covariate.names  = covars,
      prediction.grid  = futr_df,
      coordinate.system = roms_crs
    ),
    error = function(e) {
      warning("compute_extrapolation failed: ", conditionMessage(e))
      return(NULL)
    }
  )
  if (is.null(extrap)) return(NULL)

  # Extract summary percentages
  sx      <- summary(extrap)
  analog  <- sx$extrapolation$analogue.p
  univ    <- sx$extrapolation$univariate.p
  if (length(analog) == 0) analog <- 0
  if (length(univ)   == 0) univ   <- 0
  combo   <- 100 - analog - univ
  total   <- univ + combo

  # Per-covariate univariate contribution
  covar_contribs <- setNames(rep(0, length(covars)), covars)
  if (length(sx$mic) > 0) {
    mic_df <- sx$mic[[1]]  # univariate MIC table
    for (cv in covars) {
      pct <- mic_df$perc[mic_df$covariate == cv]
      if (length(pct) > 0) covar_contribs[cv] <- pct
    }
  }

  stats_df <- data.frame(
    perc_analog = analog,
    perc_univ   = univ,
    perc_combo  = combo,
    perc_total  = total,
    as.list(covar_contribs),
    check.names = FALSE
  )

  # Build ExDet raster
  exdet_vals <- extrap$data$all
  r <- terra::rast(
    data.frame(x = exdet_vals$x, y = exdet_vals$y, z = exdet_vals$ExDet),
    type = "xyz"
  )
  names(r) <- "ExDet"

  list(stats = stats_df, raster = r)
}

# ── Run all four scales for one future time step ──────────────
#' Handles ATAP, STAP (temporal subset), ATSP and STSP (spatial cells).
#' Returns a list of stats data.frames and (optionally) rasters.
#'
#' @param hist_df      full historical data.frame (env_vars + x, y, time_col)
#' @param futr_df      future data.frame for this time step (env_vars + x, y)
#' @param covars       env_vars to pass to ExDet
#' @param time_col     column name for temporal grouping in hist_df (or NULL)
#' @param time_label   value of time_col for this time step (or NULL)
#' @param cells        named list of filter functions (from spatial_grid.R); NULL = no SP
#' @param save_rasters logical; if TRUE, include rasters in output
#' @return named list: atap, stap, atsp (list by cell), stsp (list by cell)
run_exdet_all_scales <- function(hist_df, futr_df, covars,
                                  time_col    = NULL,
                                  time_label  = NULL,
                                  cells       = NULL,
                                  save_rasters = FALSE) {
  out <- list()

  # ATAP — all time, all place
  out$atap <- run_exdet_single(hist_df, futr_df, covars)
  if (!save_rasters && !is.null(out$atap)) out$atap$raster <- NULL

  # STAP — same time, all place
  if (!is.null(time_col) && !is.null(time_label)) {
    hist_t <- hist_df[hist_df[[time_col]] == time_label, ]
    out$stap <- run_exdet_single(hist_t, futr_df, covars)
    if (!save_rasters && !is.null(out$stap)) out$stap$raster <- NULL
  }

  # Spatial scales (ATSP + STSP)
  if (!is.null(cells)) {
    out$atsp <- list()
    out$stsp <- list()

    hist_cells <- split_by_cell(hist_df, cells)
    futr_cells <- split_by_cell(futr_df, cells)

    for (cell_id in names(futr_cells)) {
      futr_c <- futr_cells[[cell_id]]
      hist_c <- if (!is.null(hist_cells[[cell_id]])) hist_cells[[cell_id]] else hist_df

      # ATSP
      out$atsp[[cell_id]] <- run_exdet_single(hist_c, futr_c, covars)

      # STSP
      if (!is.null(time_col) && !is.null(time_label)) {
        hist_ct <- hist_c[hist_c[[time_col]] == time_label, ]
        if (nrow(hist_ct) > 0)
          out$stsp[[cell_id]] <- run_exdet_single(hist_ct, futr_c, covars)
      }
    }
  }

  out
}

# ── Mosaic ExDet rasters from spatial cells ───────────────────
#' Takes the per-cell ExDet raster list and mosaics back to full domain.
#' @param raster_list named list of SpatRasters (one per cell)
#' @return single SpatRaster (mean of overlapping cells; shouldn't overlap)
mosaic_cell_rasters <- function(raster_list) {
  rasters <- Filter(Negate(is.null), lapply(raster_list, `[[`, "raster"))
  if (length(rasters) == 0) return(NULL)
  if (length(rasters) == 1) return(rasters[[1]])
  do.call(terra::mosaic, c(rasters, list(fun = "mean")))
}

# ── Summarise ExDet results across cells (area-weighted) ──────
#' @param cell_results named list (one entry per cell) of ExDet results
#' @param cell_areas   named numeric vector of cell areas (km²);
#'                     if NULL, equal weighting is used
#' @return single data.frame row with weighted-average percentages
summarize_spatial_exdet <- function(cell_results, cell_areas = NULL) {
  valid <- Filter(Negate(is.null), cell_results)
  if (length(valid) == 0) return(NULL)

  stats_list <- lapply(valid, `[[`, "stats")
  ids        <- names(valid)

  if (is.null(cell_areas)) {
    weights <- rep(1 / length(ids), length(ids))
  } else {
    w <- cell_areas[ids]
    weights <- w / sum(w, na.rm = TRUE)
  }

  # weighted column means
  stat_mat <- do.call(rbind, lapply(stats_list, as.numeric))
  col_names <- names(stats_list[[1]])
  weighted_means <- colSums(stat_mat * weights, na.rm = TRUE)
  setNames(as.data.frame(t(weighted_means)), col_names)
}
