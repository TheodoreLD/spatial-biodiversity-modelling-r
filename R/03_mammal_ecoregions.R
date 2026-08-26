# =============================================================================
# 03_mammal_ecoregions.R
# Present-natural large-herbivore species list for the richest ecoregion of each modelled
# continental region. Only the Neotropic/Nearctic, Afrotropic and Indomalayan realms are
# mapped to a region (Americas, Africa, Asia); Palearctic, Australasian, Oceanian and
# Antarctic ecoregions are outside the study design and are dropped.
#
# Two ecoregions are reported for a continent whose richest one is too thinly replicated for
# 02 to attach a confidence interval to it (fewer than MIN_EFF_CI effective blocks). The
# leader is kept, because it is the answer to the question asked; the best CI-backed
# ecoregion is reported next to it, because the leader alone would present a species list
# resting on an estimate the descriptive step declines to bound. `Basis` distinguishes them
# and the effective block count travels with every row.
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

ranked <- read_csv(ranking_path, show_col_types = FALSE) |>
  left_join(eco |> st_drop_geometry() |> distinct(ECO_ID, REALM), by = "ECO_ID") |>
  mutate(continent = realm_to_continent(REALM)) |>
  filter(!is.na(continent))

# The leader is picked on the point estimate alone. The >=10-cell floor is already applied
# by 02, but the bootstrap interval in the same file is not consulted, so where the top two
# ecoregions overlap this choice is not resolved by the data.
leaders <- ranked |> group_by(continent) |> slice_max(mean, n = 1, with_ties = FALSE) |>
  ungroup() |> mutate(basis = "highest mean")

# Where the leader holds too few effective blocks for 02 to report an interval, the species
# list would rest on an estimate the descriptive step itself declines to bound. Report the
# best CI-backed ecoregion of that continent alongside it rather than silently substituting
# one for the other -- both lists go into the output, distinguished by `basis`.
backups <- ranked |> filter(continent %in% leaders$continent[is.na(leaders$lo)],
                            !is.na(lo), !ECO_ID %in% leaders$ECO_ID) |>
  group_by(continent) |> slice_max(mean, n = 1, with_ties = FALSE) |>
  ungroup() |> mutate(basis = "highest mean with a reportable CI")

target_ecos <- bind_rows(leaders, backups) |> arrange(continent, desc(mean))

for (i in which(is.na(leaders$lo)))
  warning(sprintf(paste("%s: the richest ecoregion (%s) rests on %.1f effective blocks,",
                        "too few for 02 to report an interval; %s is reported alongside it"),
                  leaders$continent[i], leaders$ECO_NAME[i], leaders$eff_blocks[i],
                  paste(backups$ECO_NAME[backups$continent == leaders$continent[i]],
                        collapse = ", ")), call. = FALSE)

if (nrow(target_ecos) == 0)
  stop("No ecoregion mapped to a modelled continent; check the REALM values: ",
       paste(unique(eco$REALM), collapse = ", "), call. = FALSE)

message("target ecoregions:")
for (i in seq_len(nrow(target_ecos))) with(target_ecos[i, ],
  message(sprintf("  %-9s %-42s mean %.2e  eff_blocks %4.1f  %s  [%s]",
                  continent, ECO_NAME, mean, eff_blocks,
                  ifelse(is.na(lo), "CI not reported     ", sprintf("CI %.2e-%.2e", lo, hi)),
                  basis)))

selected_ecos <- eco |>
  filter(ECO_ID %in% target_ecos$ECO_ID) |>
  left_join(target_ecos |> dplyr::select(ECO_ID, continent, basis, eff_blocks, lo, hi, n_cells),
            by = "ECO_ID")

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
# The reliability of the underlying richness estimate travels with every species row, so a
# list drawn from a 1.4-effective-block ecoregion cannot be read as if it were drawn from a
# well-replicated one.
ecoregion_species <- map_dfr(seq_len(nrow(selected_ecos)), function(i) {
  poly  <- selected_ecos[i, ]
  flags <- map_lgl(traits$raster_path, present_in_poly, poly_sf = poly)
  tibble(ECO_ID = poly$ECO_ID, Ecoregion = poly$ECO_NAME, Continent = poly$continent,
         Basis = poly$basis, EffBlocks = round(poly$eff_blocks, 1),
         CellsInEcoregion = poly$n_cells,
         RichnessCI = ifelse(is.na(poly$lo), NA_character_,
                             sprintf("%.2e-%.2e", poly$lo, poly$hi)),
         species = list(traits$Binomial.1.2[flags]))
})

ecoregion_unnested <- ecoregion_species |>
  unnest(cols = species) |>
  rename(Species = species)

write_csv(ecoregion_unnested, file.path(paths$ecoregion_mammal, "ecoregion_mammals.csv"))

# ---- 6. Order counts per region ----------------------------------------------
region_order_counts <- ecoregion_unnested |>
  distinct(Continent, Ecoregion, Basis, Species) |>
  left_join(traits |> dplyr::select(Binomial.1.2, Order.1.2),
            by = c("Species" = "Binomial.1.2")) |>
  count(Continent, Ecoregion, Basis, Order.1.2, name = "SpeciesCount") |>
  arrange(Continent, Ecoregion, desc(SpeciesCount))

write_csv(region_order_counts, file.path(paths$ecoregion_mammal, "ecoregion_order_counts.csv"))
message("03_mammal_ecoregions.R completed.")
