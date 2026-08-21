# run_pipeline.R — run the whole all-R pipeline in order.
#
# Run from the project root so here() resolves. Order matters: 03 reads the ecoregion
# ranking that 02 writes, and 05/06 read the stacks that 04 writes.
#
# Step 0 queries the TNRS web service, so it needs an internet connection. Its output
# is a nomenclature record that nothing downstream reads (01 matches range files by
# filename), so comment it out to run the analysis offline.

source("R/00_standardization.R")     # harmonise species names (needs internet; not read downstream)
source("R/01_biogeography.R")        # absolute + relative richness rasters
source("R/02_descriptive.R")         # richness maps + biome/ecoregion distributions + composite
source("R/03_mammal_ecoregions.R")   # megafauna of the richest ecoregion per continent
source("R/04_data_preparation.R")    # model-ready predictor + clade richness stacks
source("R/05_model1.R")              # Model 1 — fit, design justifications, robustness (see M1_RUN)
source("R/06_model2.R")              # Model 2 — mammal-clade associations
source("R/06_model2_jackknife.R")    # Model 2 — species-drop robustness of the clade ranking
source("R/06_model2_paper_figure.R") # Model 2 — assembles the combined manuscript figure
message("pipeline complete -> outputs/")
