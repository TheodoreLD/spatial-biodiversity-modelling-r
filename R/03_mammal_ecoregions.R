# =============================================================================
# 03_mammal_ecoregions.R
# Present-natural large-herbivore species list per richest ecoregion / continent.
#
# FIXES vs original:
#  - paths via config.R/here() and data/raw/... (consistent with 03/04);
#  - output directory is created before writing (was: write to a dir that may
#    not exist on a fresh clone -> error).
# =============================================================================

suppressPackageStartupMessages({
  library(sf); library(terra); library(dplyr)
  library(purrr); library(readr); library(tidyr)
})
source("R/config.R")
create_project_dirs()

traits_path <- require_file(
  here::here("data", "raw", "Range_mammals", "Traits", "Trait_data.csv"), "PHYLACINE traits")
ranges_nat  <- here::here("data", "raw", "Range_mammals", "Ranges", "Present_natural")
ecoregions_path <- require_file(
  here::here("data", "raw", "Ecoregion", "Ecoregions2017.shp"), "WWF ecoregions shapefile")

# ---- 1. Large terrestrial herbivorous mammals (>5 kg, plant-eating) ----------
traits <- read_csv(traits_path, show_col_types = FALSE) |>
  filter(Mass.g > 5000, Terrestrial == 1, Diet.Plant > 0,
         !Order.1.2 %in% c("Carnivora", "Dasyuromorphia")) |>
  mutate(raster_path = file.path(ranges_nat, paste0(Binomial.1.2, ".tif"))) |>
  filter(file.exists(raster_path))

# ---- 2. Target ecoregions (richest relative spiny-trunk richness per continent)
eco <- st_read(ecoregions_path, quiet = TRUE) |> st_make_valid()
target_ecos <- c("Beni savanna", "Cross-Niger transition forests",
                 "Tonle Sap freshwater swamp forests")
selected_ecos <- eco |> filter(ECO_NAME %in% target_ecos)

# ---- 3. Reproject to match species rasters -----------------------------------
ref_rast      <- rast(traits$raster_path[1])
selected_ecos <- st_transform(selected_ecos, crs = crs(ref_rast))

# ---- 4. Presence test --------------------------------------------------------
present_in_poly <- function(r_path, poly_sf) {
  r_vals <- terra::extract(rast(r_path), vect(poly_sf))[[2]]
  any(r_vals == 1, na.rm = TRUE)
}

# ---- 5. Species list per ecoregion -------------------------------------------
ecoregion_species <- selected_ecos |>
  dplyr::select(ECO_NAME, geometry) |>
  mutate(species = map(geometry, ~ {
    flags <- map_lgl(traits$raster_path, present_in_poly, poly_sf = .x)
    traits$Binomial.1.2[flags]
  })) |>
  st_set_geometry(NULL)

ecoregion_unnested <- ecoregion_species |>
  unnest(cols = species) |>
  rename(Ecoregion = ECO_NAME, Species = species)

write_csv(ecoregion_unnested, file.path(paths$ecoregion_mammal, "ecoregion_mammals.csv"))

# ---- 6. Order counts per region ----------------------------------------------
region_order_counts <- ecoregion_unnested |>
  distinct(Ecoregion, Species) |>
  left_join(traits |> dplyr::select(Binomial.1.2, Order.1.2),
            by = c("Species" = "Binomial.1.2")) |>
  count(Ecoregion, Order.1.2, name = "SpeciesCount") |>
  arrange(Ecoregion, desc(SpeciesCount))

write_csv(region_order_counts, file.path(paths$ecoregion_mammal, "ecoregion_order_counts.csv"))
message("03_mammal_ecoregions.R completed.")
