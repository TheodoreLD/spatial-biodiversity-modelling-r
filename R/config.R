# R/config.R
# Centralized project paths (revised: single consistent lowercase `outputs/` tree).

library(here)

create_project_dirs <- function() {
  dirs <- c(
    here('data', 'raw'),
    here('data', 'processed'),
    here('outputs'),
    here('outputs', 'biogeography', 'rasters'),
    here('outputs', 'biogeography', 'plots'),
    here('outputs', 'biogeography', 'tables'),
    here('outputs', 'ecoregion_mammal'),
    here('outputs', 'raster_model'),
    here('outputs', 'model1'),
    here('outputs', 'model2'),
    here('outputs', 'model3'),
    here('outputs', 'revision')   # environment-controlled clade models etc.
  )
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
}

require_file <- function(path, description = 'required file') {
  if (!file.exists(path)) stop(sprintf('Missing %s: %s', description, path), call. = FALSE)
  path
}

paths <- list(
  data_raw       = here('data', 'raw'),
  data_processed = here('data', 'processed'),
  outputs        = here('outputs'),
  biogeo_rasters = here('outputs', 'biogeography', 'rasters'),
  biogeo_plots   = here('outputs', 'biogeography', 'plots'),
  biogeo_tables  = here('outputs', 'biogeography', 'tables'),
  ecoregion_mammal = here('outputs', 'ecoregion_mammal'),
  raster_model   = here('outputs', 'raster_model'),
  model1         = here('outputs', 'model1'),
  model2         = here('outputs', 'model2'),
  model3         = here('outputs', 'model3'),
  revision       = here('outputs', 'revision')
)
