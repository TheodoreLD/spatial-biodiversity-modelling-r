# =============================================================================
# 03_mammal_ecoregions.R
# Present-natural large-herbivore species list for the richest ecoregion of each modelled
# continental region. Only the Neotropic/Nearctic, Afrotropic and Indomalayan realms are
# mapped to a region (Americas, Africa, Asia); Palearctic, Australasian, Oceanian and
# Antarctic ecoregions are outside the study design and are dropped.
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
traits0 <- read_csv(traits_path, show_col_types = FALSE) |>
  filter(Mass.g > 5000, Terrestrial == 1, Diet.Plant > 0,
         !Order.1.2 %in% c("Carnivora", "Dasyuromorphia", "Afrosoricida")) |>
  mutate(raster_path = file.path(ranges_nat, paste0(Binomial.1.2, ".tif")))

# A partial range download would otherwise shorten every species list without a trace.
missing <- !file.exists(traits0$raster_path)
if (any(missing))
  warning(sum(missing), " of ", nrow(traits0),
          " eligible species have no natural-range raster and are dropped: ",
          paste(head(traits0$Binomial.1.2[missing], 5), collapse = ", "),
          if (sum(missing) > 5) ", ..." else "", call. = FALSE)
traits <- traits0[!missing, ]
if (nrow(traits) == 0)
  stop("No eligible species have range rasters under ", ranges_nat, call. = FALSE)
message(nrow(traits), " large herbivores with ranges")

# ---- 2. Target ecoregions (richest relative spiny-trunk richness per continent)
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
  # Picked on the point estimate alone. The ranking file already applies the >=10-cell
  # floor, but the bootstrap interval in the same file is not consulted: where the top two
  # ecoregions overlap, this choice is not resolved by the data.
  slice_max(mean, n = 1, with_ties = FALSE) |>
  ungroup() |>
  arrange(desc(mean))

if (nrow(target_ecos) == 0)
  stop("No ecoregion mapped to a modelled continent; check the REALM values: ",
       paste(unique(eco$REALM), collapse = ", "), call. = FALSE)

message("target ecoregions: ",
        paste(sprintf("%s (%s)", target_ecos$ECO_NAME, target_ecos$continent), collapse = " | "))

selected_ecos <- eco |>
  filter(ECO_ID %in% target_ecos$ECO_ID) |>
  left_join(target_ecos |> dplyr::select(ECO_ID, continent), by = "ECO_ID")

# ---- 3. Reproject to match species rasters -----------------------------------
ref_rast <- rast(traits$raster_path[1])
crs_ref  <- crs(ref_rast)
odd <- traits$raster_path[!vapply(traits$raster_path,
                                  function(p) identical(crs(rast(p)), crs_ref), logical(1))]
if (length(odd))
  stop(length(odd), " range rasters do not share the CRS of ",
       basename(traits$raster_path[1]), ", e.g. ", basename(odd[1]), call. = FALSE)
selected_ecos <- st_transform(selected_ecos, crs = crs_ref)

# ---- 4. Presence test --------------------------------------------------------
present_in_poly <- function(r_path, poly_sf) {
  vals <- terra::extract(rast(r_path), vect(poly_sf), ID = FALSE)[[1]]
  any(vals == 1, na.rm = TRUE)
}

# ---- 5. Species list per ecoregion -------------------------------------------
ecoregion_species <- map_dfr(seq_len(nrow(selected_ecos)), function(i) {
  poly  <- selected_ecos[i, ]
  flags <- map_lgl(traits$raster_path, present_in_poly, poly_sf = poly)
  tibble(ECO_ID = poly$ECO_ID, Ecoregion = poly$ECO_NAME, Continent = poly$continent,
         species = list(traits$Binomial.1.2[flags]))
})

ecoregion_unnested <- ecoregion_species |>
  unnest(cols = species) |>
  rename(Species = species)

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
