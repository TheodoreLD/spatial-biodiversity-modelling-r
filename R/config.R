# R/config.R
# Centralized project paths (single consistent lowercase `outputs/` tree).
# Every script sources this file; no path should be hardcoded outside it.

library(here)

paths <- list(
  # ---- directories -----------------------------------------------------------
  data_raw         = here('data', 'raw'),
  data_processed   = here('data', 'processed'),
  outputs          = here('outputs'),
  biogeo_rasters   = here('outputs', 'biogeography', 'rasters'),
  descriptive      = here('outputs', 'descriptive'),
  ecoregion_mammal = here('outputs', 'ecoregion_mammal'),
  raster_model     = here('outputs', 'raster_model'),
  model1           = here('outputs', 'model1'),
  model2           = here('outputs', 'model2'),

  # ---- source species list (00 standardises its names) -----------------------
  species_xlsx     = here('data', 'spiny_trunk_species.xlsx'),

  # ---- raw layers (all from Dryad; see data/README.md) -----------------------
  woody_richness   = here('data', 'raw', 'Map_global_plant_richness',
                          'sr_Ensemble_Prediction_7774_LongLat.RData'),
  env_preds_dir    = here('data', 'raw', 'Map_environmental_predictors'),
  consumption_rds  = here('data', 'raw', 'Map_mammal_consumption', 'Map.rds'),
  mammal_traits    = here('data', 'raw', 'Range_mammals', 'Traits', 'Trait_data.csv'),
  mammal_ranges    = here('data', 'raw', 'Range_mammals', 'Ranges', 'Present_natural'),
  ecoregions       = here('data', 'raw', 'Ecoregion', 'Ecoregions2017.shp'),
  spiny_ranges_dir = here('data', 'raw', 'Range_spiny_species'),
  spiny_ranges     = here('data', 'raw', 'Range_spiny_species', 'SpeciesTotal'),

  # ---- derived inputs (built by 01 / 04, or repoint at a Dryad download) -----
  predictor_stack  = here('outputs', 'raster_model', 'predictor_stack.tif'),
  clade_stack      = here('outputs', 'raster_model', 'Clade_richness_stack.tif'),
  abs_richness     = here('outputs', 'biogeography', 'rasters', 'Absolute_Richness_Tot.tif'),
  rel_richness     = here('outputs', 'biogeography', 'rasters', 'Relative_Richness_Tot.tif'),
  borders          = here('data', 'assets', 'world_borders.geojson')
)

create_project_dirs <- function() {
  dirs <- paths[c('data_raw', 'data_processed', 'outputs', 'biogeo_rasters', 'descriptive',
                  'ecoregion_mammal', 'raster_model', 'model1', 'model2')]
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
}

require_file <- function(path, description = 'required file') {
  if (!file.exists(path)) stop(sprintf('Missing %s: %s', description, path), call. = FALSE)
  path
}

# Read a raster stack and return its layers in `expected` order.
# Prefers the layer names stored in the file; falls back to positional assignment only
# when the layer count matches, and says so. Assigning names positionally without this
# check silently mislabels every predictor if the band order ever changes.
load_stack <- function(path, expected, description = 'raster stack') {
  r <- terra::rast(require_file(path, description))
  if (!all(expected %in% names(r))) {
    if (terra::nlyr(r) != length(expected))
      stop(sprintf('%s has %d layers named [%s]; expected %d: [%s]',
                   description, terra::nlyr(r), paste(names(r), collapse = ', '),
                   length(expected), paste(expected, collapse = ', ')), call. = FALSE)
    warning(sprintf('%s: layer names not preserved in the file; assigning by position (%s)',
                    description, paste(expected, collapse = ', ')), call. = FALSE)
    names(r) <- expected
  }
  r[[expected]]
}
