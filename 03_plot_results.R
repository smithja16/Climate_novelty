##################################################################
####   STAGE 3 — PLOT RESULTS                                 ####
####                                                          ####
####   Loads the summary table and saved rasters from         ####
####   Stage 2, then produces:                                ####
####     • Smith et al Figure 3-style novelty maps (per time  ####
####       step or multi-year modal maps)                     ####
####     • Time-series of proportion novel                    ####
####     • ExDet driver bar charts                            ####
##################################################################

source("00_config.R")
source("functions/data_utils.R")
source("functions/hypervolume_fns.R")
source("functions/plot_fns.R")

library(ggplot2)
library(tidyr)     # for pivot_longer in plot_novelty_timeseries

# ── Load summary table ────────────────────────────────────────
summary_path <- file.path(output_dir, "summaries",
                           paste0("novelty_summary_", gcm_name, ".rds"))
if (!file.exists(summary_path))
  stop("Summary file not found. Run 02_run_novelty_analysis.R first.\n", summary_path)

summary_df <- readRDS(summary_path)
message("Loaded summary: ", nrow(summary_df), " time steps")

# ── Optional: load coastline for map overlays ─────────────────
coast <- NULL
# coast <- terra::vect("path/to/coast.gpkg")   # load your coast here

# ── (A) Time-series: proportion novel through time ────────────
# All scales, faceted by time_label (e.g. month)
p_ts <- plot_novelty_timeseries(
  summary_df = summary_df,
  scales     = c("atap", "stap", "atsp", "stsp"),
  facet_by   = if ("time_label" %in% names(summary_df) &&
                   !all(is.na(summary_df$time_label))) "time_label" else NULL,
  smooth     = TRUE
)
print(p_ts)
ggsave(file.path(output_dir, paste0("timeseries_", gcm_name, ".pdf")),
       p_ts, width = 12, height = 8)

# ── (B) ExDet driver plot ─────────────────────────────────────
exdet_atap_cols <- grep("^exdet_atap_.*_uni$", names(summary_df), value = TRUE)
if (length(exdet_atap_cols) > 0) {
  covars_uni <- sub("^exdet_atap_", "", sub("_uni$", "", exdet_atap_cols))
  p_drivers <- plot_exdet_drivers(
    summary_df = setNames(summary_df,
                           sub("^exdet_atap_", "", names(summary_df))),
    covars     = covars_uni,
    facet_by   = if ("time_label" %in% names(summary_df)) "time_label" else NULL
  )
  print(p_drivers)
  ggsave(file.path(output_dir, paste0("exdet_drivers_", gcm_name, ".pdf")),
         p_drivers, width = 12, height = 8)
}

# ── (C) Novelty maps ─────────────────────────────────────────
# Choose which time steps to map — either a single snapshot or a
# multi-year modal map (averaged to reduce mesoscale noise).

map_time_label <- time_levels[1]   # e.g. month 1 — change as needed
map_years      <- seq(2071, 2100)  # decade to average over (or single year)

incl_dir <- file.path(output_dir, "inclusion_rasters", gcm_name)

if (dir.exists(incl_dir)) {
  # Collect rasters for chosen time label across years
  rasters_by_year <- lapply(map_years, function(yy) {
    step_label <- if (!is.na(map_time_label))
      paste0("y", yy, "_t", map_time_label) else as.character(yy)

    paths <- list(
      atap = file.path(incl_dir, paste0("incl_atap_", step_label, ".tif")),
      stap = file.path(incl_dir, paste0("incl_stap_", step_label, ".tif")),
      atsp = file.path(incl_dir, paste0("incl_atsp_", step_label, ".tif")),
      stsp = file.path(incl_dir, paste0("incl_stsp_", step_label, ".tif"))
    )
    paths <- Filter(file.exists, paths)
    if (length(paths) == 0) return(NULL)

    incl <- lapply(paths, terra::rast)
    classify_novelty(
      atap = incl$atap,
      stap = incl$stap,
      atsp = incl$atsp,
      stsp = incl$stsp
    )
  })
  rasters_by_year <- Filter(Negate(is.null), rasters_by_year)

  if (length(rasters_by_year) > 0) {
    # Modal map across years
    novelty_map <- if (length(rasters_by_year) == 1) {
      rasters_by_year[[1]]
    } else {
      modal_classification(rasters_by_year)
    }

    map_title <- sprintf("%s — time=%s, %d–%d modal",
                         gcm_name, map_time_label,
                         min(map_years), max(map_years))
    p_map <- plot_novelty_map(novelty_map, coast = coast, title = map_title)
    print(p_map)
    ggsave(
      file.path(output_dir,
                paste0("novelty_map_", gcm_name,
                       "_t", map_time_label,
                       "_", min(map_years), "-", max(map_years), ".pdf")),
      p_map, width = 7, height = 9
    )
  }
} else {
  message("No inclusion rasters found at: ", incl_dir,
          "\n  (set save_inclusion_rasters = TRUE in config.R to generate them)")
}

message("Stage 3 complete.")
