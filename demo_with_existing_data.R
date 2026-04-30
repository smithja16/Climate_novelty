##################################################################
####   DEMO — Using simple example data                      ####
####                                                         ####
####   Runs a simplfied workflow using hist_data_example.rds ####
####   and three future raster pairs for IPSL July 2080      ####
####                                                         ####
####   Historical data has lon, lat, and mon columns, so     ####
####   all four novelty scales are demonstrated:             ####
####     ATAP — any time,  any place                         ####
####     STAP — same time, any place  (July reference)       ####
####     ATSP — any time,  same place (lon/lat included)     ####
####     STSP — same time, same place (July + lon/lat)       ####
####                                                         ####
####   This is a self-contained script — no config.R needed. ####
##################################################################

library(terra)
library(hypervolume)

# Portably set working directory to the folder containing this script,
# whether run interactively in RStudio or via Rscript from the command line.
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
} else {
  args     <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0)
    setwd(dirname(normalizePath(sub("^--file=", "", file_arg))))
}

source("functions/data_utils.R")
source("functions/hypervolume_fns.R")
source("functions/plot_fns.R")
source("functions/spatial_grid.R")

# ── Demo settings ─────────────────────────────────────────────
demo_data_dir      <- "demo_data"
output_dir         <- "outputs/demo"
gcm_name           <- "ipsl_demo"
env_vars           <- c("sst", "ild", "oxygen")
lon_col            <- "lon"
lat_col            <- "lat"
include_space      <- TRUE
time_col           <- "mon"
time_levels        <- 1:12
hv_method          <- "svm"
hv_sample_fraction <- 0.1   # kept low so demo runs quickly (~26 hypervolumes)
hv_seed            <- 123

init_output_dirs(output_dir)

# ── Load and scale historical data ────────────────────────────
message("Loading historical data ...")
hist_raw       <- readRDS(file.path(demo_data_dir, "hist_data_example.rds"))
scaling_params <- compute_scaling_params(hist_raw, env_vars)
hist_scaled    <- scale_data(hist_raw, scaling_params, env_vars)

# Scale lon and lat and save space_params for use when adding
# spatial layers to future rasters
lon_params <- list(mean = mean(hist_raw[[lon_col]]),
                   sd   = sd(hist_raw[[lon_col]]))
lat_params <- list(mean = mean(hist_raw[[lat_col]]),
                   sd   = sd(hist_raw[[lat_col]]))
hist_scaled[[lon_col]] <- (hist_raw[[lon_col]] - lon_params$mean) / lon_params$sd
hist_scaled[[lat_col]] <- (hist_raw[[lat_col]] - lat_params$mean) / lat_params$sd

space_params <- list(lon = lon_params, lat = lat_params)

saveRDS(scaling_params,
        file.path(output_dir, paste0("scaling_params_", gcm_name, ".rds")))
saveRDS(space_params,
        file.path(output_dir, paste0("space_params_", gcm_name, ".rds")))

# ── Build all hypervolumes ────────────────────────────────────
# Builds: ATAP (1), ATSP (1), STAP x12, STSP x12 = 26 total
message(sprintf("\nBuilding hypervolumes (sample_frac = %.2f) ...", hv_sample_fraction))
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
# ignore the 'consider removing some axes' warnings

# ── Load hypervolumes ─────────────────────────────────────────
hvs <- load_hypervolumes(
  gcm_name      = gcm_name,
  include_space = include_space,
  time_col      = time_col,
  time_levels   = time_levels,
  output_dir    = output_dir
)
message(sprintf("Loaded hypervolumes: %s", paste(names(hvs), collapse = ", ")))

# ── Load and scale future rasters (July 2080) ─────────────────
r_sst    <- terra::rast(file.path(demo_data_dir, "sst_ipsl_7_2080.grd"))
r_ild    <- terra::rast(file.path(demo_data_dir, "ild_ipsl_7_2080.grd"))
r_oxygen <- terra::rast(file.path(demo_data_dir, "oxygen_ipsl_7_2080.grd"))

rstack        <- c(r_sst, r_ild, r_oxygen)
names(rstack) <- env_vars
rstack_scaled <- scale_raster(rstack, scaling_params, env_vars)

# Add scaled lon/lat layers for ATSP and STSP tests
rstack_space  <- add_space_layers(rstack_scaled, space_params$lon, space_params$lat)

# ── Run inclusion tests — all four scales ─────────────────────
# time_label = 7 selects the July hypervolumes for STAP and STSP
message("\nRunning inclusion tests for July 2080 (all four scales) ...")
incl <- run_inclusion_all_scales(
  rstack_env   = rstack_scaled,
  rstack_space = rstack_space,
  hvs          = hvs,
  time_label   = 7
)

# ── Report results ────────────────────────────────────────────
scales <- c("atap", "stap", "atsp", "stsp")
labels <- c("ATAP (any time,  any place )",
            "STAP (same time, any place )",
            "ATSP (any time,  same place)",
            "STSP (same time, same place)")

message("\n── Proportion novel, IPSL July 2080 ──────────────────")
for (i in seq_along(scales)) {
  sc <- scales[i]
  r  <- incl[[sc]]
  if (!is.null(r)) {
    p <- novelty_proportion(r)
    message(sprintf("  %s : %5.1f%%", labels[i], p * 100))
  } else {
    message(sprintf("  %s : not computed", labels[i]))
  }
}

# ── Novelty classification map ────────────────────────────────
novelty_r <- classify_novelty(
  atap = incl$atap,
  stap = incl$stap,
  atsp = incl$atsp,
  stsp = incl$stsp
)
p <- plot_novelty_map(novelty_r,
                       title = "IPSL July 2080 — hierarchical novelty (demo)")
print(p)

message("\nDemo complete.")
