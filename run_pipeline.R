# run_pipeline.R — run the whole all-R pipeline in order.

source("R/00_standardization.R")     # harmonise species names
source("R/01_biogeography.R")        # absolute + relative richness rasters
source("R/02_descriptive.R")         # richness maps + biome/ecoregion distributions + composite
source("R/03_mammal_ecoregions.R")   # megafauna of the richest ecoregion per continent
source("R/04_data_preparation.R")    # model-ready predictor + clade richness stacks
source("R/05_model1.R")              # Model 1 — environmental characterisation
source("R/05_model1_sensitivity.R")  # Model 1 design justifications (resolution/transform/spatial-term)
source("R/05_model1_jackknife.R")    # Model 1 species-drop robustness
source("R/06_model2.R")              # Model 2 — mammal-clade associations
source("R/06_model2_jackknife.R")    # Model 2 — species-drop robustness of the clade ranking
message("pipeline complete -> outputs/")
