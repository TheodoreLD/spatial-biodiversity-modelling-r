# run_pipeline.R — run the whole all-R pipeline in order.
# Prereqs: source("R/requirements.R") once, install CatBoost for R (see README),
# and place the Dryad data under data/raw/ (and the derived rasters where noted).
#
# Order: name standardization -> biogeography (richness) -> DESCRIPTIVE (maps, biome/
# ecoregion distributions) -> mammal-ecoregions -> model data prep -> Model 1 -> Model 2.
# Descriptive precedes the modelling because the models (and the mammal-ecoregion step)
# build on the richness layers and richest ecoregions it characterises.
#
# Scripts 00, 01, 03, 04 use the here()-based paths in R/config.R; the descriptive and
# the models (02, 05, 06) read the derived rasters via the TRUNK_* environment variables.

# point the TRUNK_* vars at the derived rasters (Dryad output/, or your Phase-1 outputs/)
RM  <- Sys.getenv("TRUNK_RASTER_MODEL", "data/raw/output/Raster_model")
BIO <- Sys.getenv("TRUNK_BIOGEO_RAST",  "data/raw/output/Biogeography/rasters")
ECO <- Sys.getenv("TRUNK_ECO",          "data/raw/Ecoregion/Ecoregions2017.shp")
Sys.setenv(
  TRUNK_STACK   = file.path(RM,  "predictor_stack.tif"),
  TRUNK_CLADE   = file.path(RM,  "Clade_richness_stack.tif"),
  TRUNK_REL_DIR = BIO,
  TRUNK_ECO     = ECO,
  TRUNK_BORDERS = "data/assets/world_borders.geojson"
)

source("R/00_standardization.R")     # harmonise species names
source("R/01_biogeography.R")         # absolute + relative richness rasters
source("R/02_descriptive.R")          # richness maps + biome/ecoregion distributions + composite
source("R/03_mammal_ecoregions.R")    # megafauna of the richest ecoregion per continent
source("R/04_data_preparation.R")     # model-ready predictor + clade richness stacks
source("R/05_model1.R")               # Model 1 — environmental characterisation (ALE)
source("R/05_model1_sensitivity.R")   # Model 1 design justifications (resolution/transform/spatial-term)
source("R/05_model1_jackknife.R")     # Model 1 species-drop robustness (Reviewer 1); needs TRUNK_RANGES + TRUNK_ABS/REL
source("R/06_model2.R")               # Model 2 — mammal-clade associations
message("pipeline complete -> outputs/")
