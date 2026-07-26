# =============================================================================
# 01_biogeography_trunk.R  (revised)
# Absolute & relative spiny-trunk richness + biome/ecoregion summaries.
#
# FIXES vs original:
#  - single consistent lowercase output tree via config.R/here()
#    (was: "output/Biogeography/..." here vs "outputs/biogeography/..." in 03);
#  - loops over all three subsets (Tot / Prickly / Thorny) in ONE run, so it
#    always produces Relative_Richness_Tot.tif that 03 consumes
#    (was: hard-coded to "Prickly" with manual "Select" comments);
#  - data read from data/raw/... to match scripts 03/04.
# =============================================================================

suppressPackageStartupMessages({
  library(terra); library(sf); library(rnaturalearth)
  library(dplyr); library(ggplot2); library(viridis); library(readr)
})
source("R/config.R")
create_project_dirs()

# subset name -> species-range subfolder (unzip Range_spiny_species/TreeSpecies.zip here)
subsets <- c(Tot = "SpeciesTotal", Prickly = "SpeciesPrickly", Thorny = "SpeciesThorny")

predictions_path <- require_file(
  here::here("data", "raw", "Map_global_plant_richness", "sr_Ensemble_Prediction_7774_LongLat.RData"),
  "global woody-plant richness prediction")
ecoregions_path <- require_file(
  here::here("data", "raw", "Ecoregion", "Ecoregions2017.shp"), "WWF ecoregions shapefile")

# 0.1 degree WGS84 template
template <- rast(ext(-180, 180, -90, 90), res = 0.1, crs = "EPSG:4326")
values(template) <- 0

world <- ne_countries(scale = "medium", returnclass = "sf")

# Load global woody-plant richness once (object name: predictions_grid)
load(predictions_path)
nat_vect  <- vect(predictions_grid)
nat_rast  <- rasterize(nat_vect, template, field = "value", fun = "mean")

# Ecoregions once
eco <- st_read(ecoregions_path, quiet = TRUE) |> st_make_valid()

for (nm in names(subsets)) {
  message("=== Subset: ", nm, " ===")
  range_dir <- here::here("data", "raw", "Range_spiny_species", subsets[[nm]])
  if (!dir.exists(range_dir)) { warning("Missing range dir: ", range_dir, " (skipped)"); next }
  gpkg_list <- list.files(range_dir, "\\.gpkg$", full.names = TRUE)
  if (length(gpkg_list) == 0) { warning("No .gpkg in ", range_dir, " (skipped)"); next }

  # ---- Absolute richness ----
  richness <- template
  for (gpkg in gpkg_list) {
    r <- rasterize(vect(gpkg), template, background = 0, touches = TRUE)
    r[r > 0] <- 1
    richness <- richness + r
  }
  writeRaster(richness, file.path(paths$biogeo_rasters, sprintf("Absolute_Richness_%s.tif", nm)),
              overwrite = TRUE)

  df_abs <- as.data.frame(richness, xy = TRUE) |>
    setNames(c("lon", "lat", "richness")) |> filter(richness > 0)
  p_abs <- ggplot() +
    geom_tile(data = df_abs, aes(lon, lat, fill = richness)) +
    scale_fill_viridis_c(name = "Richness") +
    geom_sf(data = world, fill = NA, color = "black", linewidth = 0.3) +
    coord_sf() + labs(title = paste("Absolute richness -", nm)) + theme_minimal()
  ggsave(file.path(paths$biogeo_plots, sprintf("Absolute_Richness_%s.png", nm)),
         p_abs, width = 10, height = 6, dpi = 300)

  # ---- Relative richness ----
  nat_aligned <- resample(nat_rast, richness, method = "bilinear")
  richness[richness == 0] <- NA
  prop_rich <- richness / nat_aligned
  prop_rich[is.nan(prop_rich)] <- NA
  prop_rich[prop_rich == 0]    <- NA
  writeRaster(prop_rich, file.path(paths$biogeo_rasters, sprintf("Relative_Richness_%s.tif", nm)),
              overwrite = TRUE)

  df_prop <- as.data.frame(prop_rich, xy = TRUE) |>
    setNames(c("x", "y", "richness")) |> filter(!is.na(richness))
  p_prop <- ggplot() +
    geom_tile(data = df_prop, aes(x, y, fill = richness)) +
    geom_sf(data = world, fill = NA, color = "black", linewidth = 0.2) +
    scale_fill_viridis(name = "Prop. Richness", option = "D") +
    theme_minimal() + labs(title = paste("Relative richness -", nm)) +
    theme(axis.title = element_blank(), axis.text = element_blank(),
          axis.ticks = element_blank(), panel.grid = element_blank())
  ggsave(file.path(paths$biogeo_plots, sprintf("Relative_Richness_%s.png", nm)),
         p_prop, width = 10, height = 6, dpi = 300)

  # ---- Biome & ecoregion summaries ----
  eco_r <- st_transform(eco, crs(prop_rich))
  vals_df <- terra::extract(prop_rich, vect(eco_r), df = TRUE)
  names(vals_df)[2] <- "value"
  eco_meta <- eco_r |> st_drop_geometry() |> mutate(ID = row_number()) |>
    select(ID, ECO_NAME, BIOME_NAME)
  vals_df <- vals_df |> left_join(eco_meta, by = "ID") |>
    filter(!is.na(value), !is.na(BIOME_NAME), !is.na(ECO_NAME))

  bio_stats <- vals_df |> group_by(BIOME_NAME) |>
    summarise(mean_val = mean(value, na.rm = TRUE),
              median_val = median(value, na.rm = TRUE), n_cells = n(), .groups = "drop")
  write_csv(bio_stats, file.path(paths$biogeo_tables, sprintf("biome_summary_%s.csv", nm)))

  ecoregion_summary <- vals_df |> group_by(ID, ECO_NAME) |>
    summarise(mean_prop = mean(value, na.rm = TRUE), sd_prop = sd(value, na.rm = TRUE),
              n_cells = n(), .groups = "drop") |>
    mutate(se = sd_prop / sqrt(n_cells),
           ci_low = mean_prop - 1.96 * se, ci_high = mean_prop + 1.96 * se)
  write_csv(ecoregion_summary, file.path(paths$biogeo_tables, sprintf("ecoregion_summary_%s.csv", nm)))
}

message("01_biogeography_trunk.R completed.")
