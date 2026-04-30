##################################################################
####   CLIMATE NOVELTY ANALYSIS — USER CONFIGURATION          ####
####   Edit this file, then source the numbered scripts       ####
##################################################################
#
# Workflow:
#   1. Edit this config.R            (once per GCM)
#   2. Run 01_build_hypervolumes.R   (once per GCM)
#   3. Run 02_run_novelty_analysis.R (once per GCM)
#   4. Run 03_plot_results.R
#
# Requirements:
#   terra (>= 1.6), hypervolume (>= 3.1), dsmextra (>= 1.1.5),
#   dplyr, ggplot2, sf


# ── Paths ─────────────────────────────────────────────────────
base_dir   <- "."                 # working directory (where this config lives)
output_dir <- file.path(base_dir, "outputs")   # all outputs written here
data_dir   <- file.path(base_dir, "data")      # data is saved here

# ── GCM / model run ───────────────────────────────────────────
# One analysis is run at a time; re-source config and re-run scripts for each GCM.
gcm_name <- "ipsl"   # label used in saved output filenames

# ── Environmental variables ───────────────────────────────────
# Must match column names in historical data AND layer names in future raster stacks.
env_vars <- c("sst", "ild", "oxygen")

# ── Historical reference data ─────────────────────────────────
# A data.frame (or tibble) with one row per observation.
# Required columns: all env_vars.
# Optional columns: lon, lat (if include_space = TRUE), time_col (see below).
hist_file <- file.path(data_dir, paste0("hist_data_", gcm_name, ".rds"))

# ── Spatial dimensions in hypervolume (ATSP / STSP) ───────────
# TRUE  = include lon and lat as hypervolume dimensions ("with Space" hypervolumes)
# FALSE = environmental dimensions only ("without Space" hypervolumes)
include_space <- TRUE
lon_col <- "lon"   # column name for longitude in hist data
lat_col <- "lat"   # column name for latitude  in hist data

# ── Temporal resolution ───────────────────────────────────────
# time_col: the column in hist data that identifies the time group (season, month, doy, …).
#   Set to NULL for a single annual hypervolume (no temporal stratification).
# time_levels: which values of time_col to build separate (STAP/STSP) hypervolumes for.
#   e.g. 1:12 for months, 1:365 for day-of-year, c("spring","summer") for seasons.
time_col    <- "month"    # column name in hist data; set NULL for annual only
time_levels <- 1:12       # unique time groups to fit separate hypervolumes for

# ── Future time steps ─────────────────────────────────────────
# A data.frame with (at minimum) columns:
#   year       : calendar year of the future time step
#   time_label : value matching time_levels (or NA if no temporal stratification)
#   raster_path: full path to the raster stack for this time step
#                (layers must be named with env_vars, in any order)
#
# Example: build from a folder of rasters named "sst_ipsl_7_2080.tif" etc.
#   See helper function `build_future_timesteps()` in functions/data_utils.R
#
future_timesteps <- local({
  expand.grid(time_label = 1:12, year = 2010:2100) |>
    transform(raster_path = file.path(
      data_dir, "future_rasters",
      paste0("stack_", gcm_name, "_", time_label, "_", year, ".tif")
    ))
})
# ↑ Replace with your actual file paths or use build_future_timesteps() helper.

# ── Spatial grid for "same place" novelty (ATSP / STSP) ───────
# NULL  = no spatial subgrouping; only ATAP and STAP scales are computed.
# "auto"= build a regular grid from lat_breaks + dist_breaks (requires coast data).
# sf/SpatVector = a polygon layer; each feature defines one spatial cell.
#   Use sf::read_sf() or terra::vect() to load your own region polygons.
#
spatial_grid_type <- "auto"   # NULL | "auto" | "custom"

# Settings for spatial_grid_type = "auto":
lat_breaks   <- c(34.5, 39.0, 43.5)      # southern latitude boundaries of rows
dist_breaks  <- c(250, 500, 750)          # km from coast; defines columns
coast_file   <- file.path(data_dir, "coast_data.rds")  # SpatVector or data.frame of coast points

# Settings for spatial_grid_type = "custom":
# region_file  <- file.path(data_dir, "my_regions.gpkg")  # polygon layer, one feature per region

# ── Hypervolume settings ──────────────────────────────────────
hv_method          <- "svm"   # "svm" | "box" | "gaussian"
hv_sample_fraction <- 0.05    # fraction of historical rows to use (1 = all; reduce for speed)
hv_seed            <- 123

# ── Inclusion test settings ───────────────────────────────────
inclusion_mode <- "accurate"  # "accurate" | "fast"

# ── ExDet settings ────────────────────────────────────────────
run_exdet <- FALSE   # FALSE to skip ExDet and use hypervolume only

# ── Output options ────────────────────────────────────────────
save_inclusion_rasters <- TRUE   # save per-timestep inclusion rasters (.tif)
save_exdet_rasters     <- TRUE   # save per-timestep ExDet rasters (.tif)
save_summary_table     <- TRUE   # save combined novelty summary as .rds

# ── Parallelisation (optional) ────────────────────────────────
# 1 = sequential (safe default). Set higher only if your system supports it.
# When measuring novelty for many time steps parallelisation is ~necessary.
n_cores <- 1
