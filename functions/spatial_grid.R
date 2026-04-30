##################################################################
####   SPATIAL GRID FUNCTIONS                                  ####
####   Define and apply regional subgrouping for              ####
####   "same place" novelty (ATSP / STSP)                     ####
##################################################################
#
# Three modes (set spatial_grid_type in config.R):
#   NULL   — no spatial subgrouping; ATAP + STAP only
#   "auto" — regular grid from lat_breaks + dist_breaks
#   "custom" — user-supplied sf/SpatVector polygon layer

# terra is used via terra:: prefix throughout; run install_packages.R if needed.

# ── Distance to coast ─────────────────────────────────────────
#' Compute minimum great-circle distance (km) from each (lon, lat)
#' point to the nearest coast point.
#'
#' @param lon        numeric vector of longitudes
#' @param lat        numeric vector of latitudes
#' @param coast_data data.frame/matrix with x/lon and y/lat columns,
#'                   or a SpatVector of coast points/lines
#' @return numeric vector of distances in km
compute_dist_coast <- function(lon, lat, coast_data) {
  if (inherits(coast_data, "SpatVector")) {
    coast_xy <- terra::crds(coast_data)
  } else {
    xcol <- intersect(c("x", "lon"), names(coast_data))[1]
    ycol <- intersect(c("y", "lat"), names(coast_data))[1]
    coast_xy <- as.matrix(coast_data[, c(xcol, ycol)])
  }

  pts_v   <- terra::vect(data.frame(x = lon, y = lat),
                         geom = c("x", "y"), crs = "EPSG:4326")
  coast_v <- terra::vect(
    data.frame(x = coast_xy[, 1], y = coast_xy[, 2]),
    geom = c("x", "y"), crs = "EPSG:4326"
  )
  nn <- terra::nearest(pts_v, coast_v)
  nn$distance / 1000   # metres → km
}

# ── Build an "auto" grid from lat + distance breaks ───────────
#' Creates a SpatVector of rectangular (in lat/dist space) regions.
#' Each region is labelled "r{row}c{col}" (row = lat band, col = dist band).
#'
#' @param lat_breaks  southern boundaries of lat bands (northern boundary
#'                    is inferred as the maximum lat in the domain + epsilon)
#' @param dist_breaks km-from-coast boundaries between columns
#'   (first column: 0 → dist_breaks[1]; last column: dist_breaks[end] → Inf)
#' @return list of named filter functions; each takes (data with x,y,dist_coast)
#'         and returns logical index vector
build_auto_grid <- function(lat_breaks, dist_breaks) {
  # build all combinations of lat band + dist band
  lat_lo  <- c(-Inf, lat_breaks)
  lat_hi  <- c(lat_breaks,  Inf)
  dist_lo <- c(0,    dist_breaks)
  dist_hi <- c(dist_breaks, Inf)

  cells <- list()
  for (ri in seq_along(lat_lo)) {
    for (ci in seq_along(dist_lo)) {
      label <- sprintf("r%dc%d", ri, ci)
      # capture by value (not reference) using local()
      cells[[label]] <- local({
        ll <- lat_lo[ri]
        lh <- lat_hi[ri]
        dl <- dist_lo[ci]
        dh <- dist_hi[ci]
        function(data) {
          data$y          >= ll &
            data$y          <  lh &
            data$dist_coast >= dl &
            data$dist_coast <  dh
        }
      })
    }
  }
  cells
}

# ── Build a "custom" grid from polygon features ───────────────
#' @param region_file path to polygon file readable by terra::vect()
#' @param id_col      name of the attribute that labels each region
#' @return named list of filter functions (same API as build_auto_grid)
build_custom_grid <- function(region_file, id_col = "region_id") {
  regions <- terra::vect(region_file)
  ids <- terra::values(regions)[[id_col]]
  structure(
    lapply(seq_along(ids), function(i) {
      poly_i <- regions[i, ]
      function(data) {
        pts <- terra::vect(data, geom = c("x", "y"), crs = "EPSG:4326")
        !is.na(terra::extract(poly_i, pts)[[id_col]])
      }
    }),
    names = as.character(ids)
  )
}

# ── Resolve spatial grid from config ─────────────────────────
#' Call once after sourcing config.R.
#' Returns NULL (no spatial split) or a named list of filter functions.
resolve_spatial_grid <- function(spatial_grid_type,
                                 lat_breaks    = NULL,
                                 dist_breaks   = NULL,
                                 coast_file    = NULL,
                                 region_file   = NULL,
                                 id_col        = "region_id") {
  if (is.null(spatial_grid_type)) return(NULL)

  if (spatial_grid_type == "auto") {
    if (is.null(lat_breaks) || is.null(dist_breaks))
      stop("spatial_grid_type='auto' requires lat_breaks and ",
           "dist_breaks in config.R")
    return(build_auto_grid(lat_breaks, dist_breaks))
  }

  if (spatial_grid_type == "custom") {
    if (is.null(region_file))
      stop("spatial_grid_type='custom' requires region_file in config.R")
    return(build_custom_grid(region_file, id_col))
  }

  stop("spatial_grid_type must be NULL, 'auto', or 'custom'")
}

# ── Add dist_coast to a data.frame (if using "auto" grid) ─────
#' Attaches a dist_coast column to data (computing if missing).
#' Skips computation if the column already exists.
ensure_dist_coast <- function(data, coast_data,
                              lon_col = "x", lat_col = "y") {
  if ("dist_coast" %in% names(data)) return(data)
  message("Computing distance-to-coast … (this may take a few minutes)")
  data$dist_coast <- compute_dist_coast(data[[lon_col]], data[[lat_col]],
                                        coast_data)
  data
}
