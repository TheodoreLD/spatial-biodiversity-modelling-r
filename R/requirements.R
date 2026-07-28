# requirements.R — install the R packages the pipeline uses
required_packages <- c(
  "here", "terra", "sf", "dplyr", "tidyr", "purrr", "readr", "stringr",
  "ggplot2", "viridis", "patchwork", "scales", "tibble",
  "readxl", "TNRS"   # step 0: read the species xlsx + resolve names (TNRS needs internet)
)
installed <- rownames(installed.packages())
to_install <- setdiff(required_packages, installed)
if (length(to_install) > 0) install.packages(to_install, repos = "https://cloud.r-project.org")

cat("Base project packages installed.\n")
cat("IMPORTANT: CatBoost must be installed manually if unavailable through your R setup.\n")
cat("See: https://catboost.ai/docs/en/installation/r-installation\n")
