# requirements.R  (revised) — install R packages for the pipeline + revision analyses
required_packages <- c(
  "here", "terra", "sf", "dplyr", "tidyr", "purrr", "readr",
  "ggplot2", "viridis", "forcats", "rnaturalearth", "blockCV",
  "ppcor"   # added: partial correlations for the environment-controlled clade models
)
installed <- rownames(installed.packages())
to_install <- setdiff(required_packages, installed)
if (length(to_install) > 0) install.packages(to_install, repos = "https://cloud.r-project.org")

cat("Base project packages installed.\n")
cat("IMPORTANT: CatBoost must be installed manually if unavailable through your R setup.\n")
cat("See: https://catboost.ai/docs/en/installation/r-installation\n")
