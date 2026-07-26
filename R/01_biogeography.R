# =============================================================================
# 01_biogeography.R — absolute & relative spiny-trunk richness rasters.
# Rasterises the species ranges (total / prickly / thorny) and divides by global
# woody-plant richness. Maps and the biome/ecoregion distributions are produced
# downstream by 02_descriptive.R (with country borders and block-bootstrap CIs) —
# they are NOT duplicated here.
# =============================================================================
suppressPackageStartupMessages({ library(terra); library(dplyr); library(readr) })
source("R/config.R")
create_project_dirs()

# subset name -> species-range subfolder (unzip Range_spiny_species/TreeSpecies.zip here)
subsets <- c(Tot = "SpeciesTotal", Prickly = "SpeciesPrickly", Thorny = "SpeciesThorny")
predictions_path <- require_file(
  here::here("data", "raw", "Map_global_plant_richness", "sr_Ensemble_Prediction_7774_LongLat.RData"),
  "global woody-plant richness prediction")

# 0.1 degree WGS84 template + global woody-plant richness (object name: predictions_grid)
template <- rast(ext(-180, 180, -90, 90), res = 0.1, crs = "EPSG:4326"); values(template) <- 0
load(predictions_path)
nat_rast <- rasterize(vect(predictions_grid), template, field = "value", fun = "mean")

for (nm in names(subsets)) {
  message("=== Subset: ", nm, " ===")
  range_dir <- here::here("data", "raw", "Range_spiny_species", subsets[[nm]])
  if (!dir.exists(range_dir)) { warning("Missing range dir: ", range_dir, " (skipped)"); next }
  gpkg_list <- list.files(range_dir, "\\.gpkg$", full.names = TRUE)
  if (length(gpkg_list) == 0) { warning("No .gpkg in ", range_dir, " (skipped)"); next }

  # absolute richness = number of species ranges overlapping each 0.1 deg cell
  richness <- template
  for (gpkg in gpkg_list) {
    r <- rasterize(vect(gpkg), template, background = 0, touches = TRUE); r[r > 0] <- 1
    richness <- richness + r
  }
  writeRaster(richness, file.path(paths$biogeo_rasters, sprintf("Absolute_Richness_%s.tif", nm)), overwrite = TRUE)

  # relative richness = absolute / total woody-plant richness (Cai et al. 2023)
  nat_aligned <- resample(nat_rast, richness, method = "bilinear")
  richness[richness == 0] <- NA
  prop_rich <- richness / nat_aligned
  prop_rich[is.nan(prop_rich)] <- NA; prop_rich[prop_rich == 0] <- NA
  writeRaster(prop_rich, file.path(paths$biogeo_rasters, sprintf("Relative_Richness_%s.tif", nm)), overwrite = TRUE)
}
message("01_biogeography.R completed (richness rasters -> ", paths$biogeo_rasters, ")")
