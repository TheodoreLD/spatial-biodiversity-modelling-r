# run_pipeline.R — run the whole all-R pipeline in order.
#
# Prereqs: source("R/requirements.R") once, install CatBoost for R (see README), and
# place the Dryad data under data/raw/ (see data/README.md).
#
# All input/output paths live in R/config.R (the `paths` list, here()-based); each
# script sources it. Run this from the project root (or with the RStudio project open)
# so here() resolves correctly. Descriptive precedes the modelling because the models
# and the mammal-ecoregion step build on the richness layers it produces.

source("R/00_standardization.R")     # harmonise species names
source("R/01_biogeography.R")        # absolute + relative richness rasters
source("R/02_descriptive.R")         # richness maps + biome/ecoregion distributions + composite
source("R/03_mammal_ecoregions.R")   # megafauna of the richest ecoregion per continent
source("R/04_data_preparation.R")    # model-ready predictor + clade richness stacks
source("R/05_model1.R")              # Model 1 — environmental characterisation (ALE)
source("R/05_model1_sensitivity.R")  # Model 1 design justifications (resolution/transform/spatial-term)
source("R/05_model1_jackknife.R")    # Model 1 species-drop robustness (Reviewer 1)
source("R/06_model2.R")              # Model 2 — mammal-clade associations
message("pipeline complete -> outputs/")
