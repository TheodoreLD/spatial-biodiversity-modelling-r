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
    here('outputs', 'model2')
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

  # source species list (00 standardises its names) + derived inputs consumed by
  # the descriptive step (02) and the models (05/06); built by 01/04, or repoint
  # these at a Dryad download.
  species_xlsx    = here('data', 'spiny_trunk_species.xlsx'),
  predictor_stack = here('outputs', 'raster_model', 'predictor_stack.tif'),
  clade_stack     = here('outputs', 'raster_model', 'Clade_richness_stack.tif'),
  abs_richness    = here('outputs', 'biogeography', 'rasters', 'Absolute_Richness_Tot.tif'),
  rel_richness    = here('outputs', 'biogeography', 'rasters', 'Relative_Richness_Tot.tif'),
  ecoregions      = here('data', 'raw', 'Ecoregion', 'Ecoregions2017.shp'),
  spiny_ranges    = here('data', 'raw', 'Range_spiny_species', 'SpeciesTotal'),
  borders         = here('data', 'assets', 'world_borders.geojson')
)
