# =============================================================================
# 01_biogeography.R — absolute & relative spiny-trunk richness rasters.
# Rasterises the species ranges (total / prickly / thorny) and divides by global
# vascular plant richness (Cai et al. 2023).
# =============================================================================
suppressPackageStartupMessages({ library(terra); library(dplyr); library(readr) })
source("R/config.R")
create_project_dirs()

# subset name -> species-range subfolder under paths$spiny_ranges_dir
subsets <- c(Tot = "SpeciesTotal", Prickly = "SpeciesPrickly", Thorny = "SpeciesThorny")
predictions_path <- require_file(paths$plant_richness,
                                 "global vascular plant richness prediction")

# 0.1 degree WGS84 template + global vascular plant richness (object: predictions_grid)
template <- rast(ext(-180, 180, -90, 90), res = 0.1, crs = "EPSG:4326"); values(template) <- 0
load(predictions_path)
nat_rast <- rasterize(vect(predictions_grid), template, field = "value", fun = "mean")
# Saved because Model 1's ratio check needs it as a predictor: the response is a ratio,
# and the only way to show the denominator is not generating the environmental signal is
# to model the numerator with the denominator entered as a covariate.
writeRaster(nat_rast, file.path(paths$biogeo_rasters, "Plant_Richness.tif"),
            overwrite = TRUE)

for (nm in names(subsets)) {
  message("=== Subset: ", nm, " ===")
  range_dir <- file.path(paths$spiny_ranges_dir, subsets[[nm]])
  if (!dir.exists(range_dir)) { warning("Missing range dir: ", range_dir, " (skipped)"); next }
  gpkg_list <- list.files(range_dir, "\\.gpkg$", full.names = TRUE)
  if (length(gpkg_list) == 0) { warning("No .gpkg in ", range_dir, " (skipped)"); next }

  # absolute richness = number of species ranges overlapping each 0.1 deg cell
  richness <- template
  for (gpkg in gpkg_list) {
    # touches = TRUE: a cell counts as occupied if the range polygon touches it at all,
    # so each range is inflated by up to one 0.1 deg cell along its edge. Applied
    # identically to all species, so it shifts absolute richness slightly upward at range
    # margins rather than favouring any species.
    r <- rasterize(vect(gpkg), template, background = 0, touches = TRUE); r[r > 0] <- 1
    richness <- richness + r
  }
  writeRaster(richness, file.path(paths$biogeo_rasters, sprintf("Absolute_Richness_%s.tif", nm)), overwrite = TRUE)

  # relative richness = absolute / total VASCULAR plant richness (Cai et al. 2023). The
  # denominator's native grain is 7,774 km2, so all fine structure in the ratio comes from
  # the numerator.
  # A no-op as things stand - nat_rast was rasterised onto `template` and `richness` is
  # built from it - kept only as a guard should either geometry change. The denominator is
  # NOT genuinely 0.1 deg: its native grain is 7,774 km2 and it was upsampled at the
  # rasterize() above, so the fine structure of the ratio comes entirely from the numerator.
  nat_aligned <- resample(nat_rast, richness, method = "bilinear")
  richness[richness == 0] <- NA
  prop_rich <- richness / nat_aligned
  prop_rich[!is.finite(prop_rich) | prop_rich <= 0] <- NA
  writeRaster(prop_rich, file.path(paths$biogeo_rasters, sprintf("Relative_Richness_%s.tif", nm)), overwrite = TRUE)
}
message("01_biogeography.R completed (richness rasters -> ", paths$biogeo_rasters, ")")
