# =============================================================================
# 03_mammal_ecoregions.R
# Present-natural large-herbivore species list per richest ecoregion / continent.
# =============================================================================

suppressPackageStartupMessages({
  library(sf); library(terra); library(dplyr)
  library(purrr); library(readr); library(tidyr)
})
source("R/config.R")
create_project_dirs()

traits_path     <- require_file(paths$mammal_traits, "PHYLACINE traits")
ranges_nat      <- paths$mammal_ranges
ecoregions_path <- require_file(paths$ecoregions, "WWF ecoregions shapefile")

# ---- 1. Large terrestrial herbivorous mammals (>5 kg, plant-eating) ----------
traits <- read_csv(traits_path, show_col_types = FALSE) |>
  filter(Mass.g > 5000, Terrestrial == 1, Diet.Plant > 0,
         !Order.1.2 %in% c("Carnivora", "Dasyuromorphia")) |>
  mutate(raster_path = file.path(ranges_nat, paste0(Binomial.1.2, ".tif"))) |>
  filter(file.exists(raster_path))

# ---- 2. Target ecoregions (richest relative spiny-trunk richness per continent)
# Derived from the ranking that 02_descriptive.R writes, not hardcoded, so rerunning
# the descriptive step with different data updates this one automatically. Continents
# come from the Ecoregions2017 REALM field; Palearctic / Australasia / Oceania are not
# analysed, matching the three continents modelled in 06_model2.R.
eco <- st_read(ecoregions_path, quiet = TRUE) |> st_make_valid()
if (!"REALM" %in% names(eco)) stop("Ecoregions shapefile has no REALM field", call. = FALSE)

realm_to_continent <- function(r) {
  r <- tolower(as.character(r))
  ifelse(grepl("neotrop|nearctic", r), "Americas",
  ifelse(grepl("afrotrop", r),         "Africa",
  ifelse(grepl("indo.?malay", r),      "Asia", NA_character_)))
}

ranking_path <- require_file(
  file.path(paths$descriptive, "ecoregion_summary_total.csv"),
  "ecoregion ranking from 02_descriptive.R (run step 2 first)")

target_ecos <- read_csv(ranking_path, show_col_types = FALSE) |>
  left_join(eco |> st_drop_geometry() |> distinct(ECO_ID, REALM), by = "ECO_ID") |>
  mutate(continent = realm_to_continent(REALM)) |>
  filter(!is.na(continent)) |>
  group_by(continent) |>
  slice_max(mean, n = 1, with_ties = FALSE) |>
  ungroup() |>
  arrange(desc(mean))

if (nrow(target_ecos) == 0)
  stop("No ecoregion mapped to a modelled continent; check the REALM values: ",
       paste(unique(eco$REALM), collapse = ", "), call. = FALSE)

message("target ecoregions: ",
        paste(sprintf("%s (%s)", target_ecos$ECO_NAME, target_ecos$continent), collapse = " | "))
selected_ecos <- eco |> filter(ECO_NAME %in% target_ecos$ECO_NAME)

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
