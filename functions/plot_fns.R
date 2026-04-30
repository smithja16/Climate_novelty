##################################################################
####   PLOTTING FUNCTIONS                                      ####
####   Novelty maps (Figure 3 style) and time series          ####
##################################################################

library(ggplot2)
# terra is used via terra:: prefix; run install_packages.R if needed.

# ── Classify pixels into novelty hierarchy ────────────────────
#' Applies the hierarchical novelty classification:
#'   1 = Analog       (not novel at any scale)
#'   2 = ATAP         (novel at broadest scale)
#'   3 = STAP or ATSP (novel only at temporal or spatial scale, not ATAP)
#'   4 = STSP         (novel only at most stringent scale)
#'
#' @param atap, stap, atsp, stsp  SpatRasters (0=novel, 1=analog); NULL if not computed
#' @return SpatRaster with integer values 1–4 and NA for missing data
classify_novelty <- function(atap, stap = NULL, atsp = NULL, stsp = NULL) {
  # start with analog everywhere
  val <- atap; terra::values(val) <- NA

  valid <- !is.na(terra::values(atap))
  v <- rep(1L, sum(valid))   # default: analog

  atap_v <- terra::values(atap)[valid]
  v[atap_v == 0] <- 2L

  if (!is.null(stap)) {
    stap_v <- terra::values(stap)[valid]
    v[atap_v == 1 & stap_v == 0] <- 3L
  }
  if (!is.null(atsp)) {
    atsp_v <- terra::values(atsp)[valid]
    v[atap_v == 1 & atsp_v == 0] <- 3L
  }
  if (!is.null(stap) && !is.null(atsp) && !is.null(stsp)) {
    stsp_v  <- terra::values(stsp)[valid]
    stap_v2 <- if (is.null(stap)) rep(1L, sum(valid)) else terra::values(stap)[valid]
    atsp_v2 <- if (is.null(atsp)) rep(1L, sum(valid)) else terra::values(atsp)[valid]
    v[atap_v == 1 & stap_v2 == 1 & atsp_v2 == 1 & stsp_v == 0] <- 4L
  }

  result_v              <- terra::values(val)
  result_v[valid]       <- v
  terra::values(val)    <- result_v
  val
}

# ── Multi-year modal classification ──────────────────────────
#' Averages a classification raster across multiple years by taking
#' the most frequent category (mode) per pixel — smooths inter-annual noise.
#'
#' @param raster_list list of SpatRasters from classify_novelty()
#' @return SpatRaster of modal category
modal_classification <- function(raster_list) {
  stacked <- terra::rast(raster_list)
  terra::modal(stacked, na.rm = TRUE)
}

# ── Figure-3 style map ────────────────────────────────────────
#' @param novelty_raster  output of classify_novelty() or modal_classification()
#' @param coast           SpatVector or sf of coastline (optional)
#' @param xlim, ylim      axis limits (default: full extent)
#' @param title           plot title
#' @param legend_labels   character vector of 4 category labels
#' @return ggplot object
plot_novelty_map <- function(novelty_raster,
                              coast         = NULL,
                              xlim          = NULL,
                              ylim          = NULL,
                              title         = "Climate novelty",
                              legend_labels = c("Analog",
                                                "Novel (any time & place)",
                                                "Novel (time OR place)",
                                                "Novel (same time & place)")) {
  # raster to data.frame
  df <- as.data.frame(novelty_raster, xy = TRUE, na.rm = TRUE)
  names(df)[1:3] <- c("lon", "lat", "category")
  df$category <- factor(df$category, levels = 1:4, labels = legend_labels)

  pal <- c("1" = "#2E8B57",   # analog: sea green
           "2" = "#CC3333",   # ATAP:   red
           "3" = "#E69F00",   # STAP|ATSP: orange
           "4" = "#F0E442")   # STSP only: yellow
  names(pal) <- legend_labels

  p <- ggplot() +
    geom_raster(data = df, aes(x = lon, y = lat, fill = category)) +
    scale_fill_manual(values = pal, name = "Novelty scale",
                      na.value = "grey80", drop = FALSE) +
    coord_fixed(ratio = 1) +
    labs(title = title, x = "Longitude", y = "Latitude") +
    theme_bw(base_size = 11)

  if (!is.null(xlim)) p <- p + xlim(xlim)
  if (!is.null(ylim)) p <- p + ylim(ylim)

  if (!is.null(coast)) {
    coast_df <- as.data.frame(terra::geom(terra::as.lines(
      if (!inherits(coast, "SpatVector")) terra::vect(coast) else coast
    )))[, c("x", "y", "geom")]
    p <- p + geom_path(data = coast_df, aes(x = x, y = y, group = geom),
                       colour = "black", linewidth = 0.3)
  }
  p
}

# ── Novelty time series (proportion novel through time) ───────
#' @param summary_df  data.frame with columns: year, time_label (optional),
#'                    and columns named e.g. prop_novel_atap, prop_novel_stap, …
#' @param scales      character vector of scales to plot
#' @param facet_by    column to facet by (e.g. "time_label" for monthly panels)
#' @return ggplot object
plot_novelty_timeseries <- function(summary_df,
                                     scales    = c("atap", "stap", "atsp", "stsp"),
                                     facet_by  = NULL,
                                     smooth    = TRUE) {
  # pivot longer
  prop_cols <- paste0("prop_novel_", scales)
  present   <- intersect(prop_cols, names(summary_df))
  if (length(present) == 0) stop("No prop_novel_* columns found in summary_df")

  long_df <- tidyr::pivot_longer(
    summary_df,
    cols        = all_of(present),
    names_to    = "scale",
    names_prefix = "prop_novel_",
    values_to   = "prop_novel"
  )
  long_df$scale <- toupper(long_df$scale)

  scale_pal <- c(ATAP = "#CC3333", STAP = "#E69F00",
                 ATSP = "#56B4E9", STSP = "#F0E442")

  p <- ggplot(long_df, aes(x = year, y = prop_novel * 100,
                            colour = scale)) +
    geom_line(alpha = 0.4) +
    scale_colour_manual(values = scale_pal, name = "Scale") +
    ylim(0, 100) +
    labs(x = NULL, y = "% area novel") +
    theme_bw(base_size = 11)

  if (smooth) {
    p <- p + geom_smooth(method = "gam", formula = y ~ s(x),
                         se = FALSE, linewidth = 1)
  }
  if (!is.null(facet_by) && facet_by %in% names(long_df)) {
    p <- p + facet_wrap(as.formula(paste("~", facet_by)))
  }
  p
}

# ── ExDet driver bar chart ────────────────────────────────────
#' Stack bar chart of which variable drives univariate extrapolation.
#' @param summary_df   data.frame; must have columns yr (or year), mo (or time_label),
#'                     and one column per covariate ending in _uni
#' @param covars       variable names (will look for {var}_uni columns)
#' @param facet_by     column to facet by (e.g. "mo")
plot_exdet_drivers <- function(summary_df, covars, facet_by = NULL,
                                palette = NULL) {
  uni_cols <- paste0(covars, "_uni")
  present  <- intersect(uni_cols, names(summary_df))

  yr_col <- intersect(c("year", "yr"), names(summary_df))[1]
  id_cols <- c(yr_col, if (!is.null(facet_by)) facet_by)
  long_df <- tidyr::pivot_longer(summary_df, cols = all_of(present),
                                  names_to = "variable",
                                  names_suffix = "_uni",
                                  values_to = "perc")
  long_df$variable <- sub("_uni$", "", long_df$variable)

  if (is.null(palette))
    palette <- scales::hue_pal()(length(present))

  p <- ggplot(long_df, aes(x = .data[[yr_col]], y = perc, fill = variable)) +
    geom_bar(stat = "identity", position = "stack", width = 1) +
    scale_fill_manual(values = setNames(palette, sub("_uni$", "", present)),
                      name = "Variable") +
    ylim(0, 100) +
    labs(x = NULL, y = "% univariate extrapolation") +
    theme_bw(base_size = 11)

  if (!is.null(facet_by) && facet_by %in% names(long_df))
    p <- p + facet_wrap(as.formula(paste("~", facet_by)))
  p
}
