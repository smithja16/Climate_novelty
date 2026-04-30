##################################################################
####   PACKAGE INSTALLATION                                    ####
####   Run this once before using the refactored workflow      ####
##################################################################

pkgs_cran <- c(
  "terra",        # modern raster handling (replaces raster/rgdal/sp)
  "hypervolume",  # hypervolume novelty analysis
  "dsmextra",     # ExDet extrapolation detection
  "dplyr",
  "tidyr",
  "ggplot2",
  "sf",           # vector spatial data
  "mgcv",         # GAM smoothing for time series plots
  "scales"        # colour palettes in ggplot2
)

# raster is still needed internally by dsmextra and as a hypervolume fallback
pkgs_legacy <- c("raster", "sp")

to_install <- setdiff(c(pkgs_cran, pkgs_legacy),
                      rownames(installed.packages()))

if (length(to_install) > 0) {
  message("Installing: ", paste(to_install, collapse = ", "))
  install.packages(to_install)
} else {
  message("All packages already installed.")
}

# Verify terra loads correctly
if (!requireNamespace("terra", quietly = TRUE))
  stop("terra failed to install — check your R version (>= 4.1 required)")

packageVersion("terra")
