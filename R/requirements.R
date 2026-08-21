# requirements.R — install the R packages the pipeline uses
required_packages <- c(
  "here", "terra", "sf", "dplyr", "tidyr", "purrr", "readr", "stringr",
  "ggplot2", "patchwork", "scales", "tibble",
  "readxl", "TNRS",
  # No library() call names these two: the mammal-consumption RDS stores a raster-package
  # RasterLayer, so readRDS needs the S4 class defined before terra::rast() can coerce it.
  "raster", "sp"
)

installed <- rownames(installed.packages())
to_install <- setdiff(required_packages, installed)
if (length(to_install) > 0) install.packages(to_install, repos = "https://cloud.r-project.org")

cat("Base project packages installed.\n")
cat("IMPORTANT: CatBoost must be installed manually if unavailable through your R setup.\n")
cat("See: https://catboost.ai/docs/en/installation/r-installation\n")
