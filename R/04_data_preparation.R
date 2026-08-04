# =============================================================================
# 04_data_preparation.R
# Prepare predictor and mammal clade richness stacks for downstream modeling.
# Writes the canonical stacks (predictor_stack.tif / Clade_richness_stack.tif) that
# 05_model1.R and 06_model2.R read (paths defined in R/config.R).
# =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(readr)
})

source("R/config.R")
create_project_dirs()

# ---- Helpers -----------------------------------------------------------------
# Reproject onto the analysis grid, then discard anything outside the source's own
# range. Bilinear resampling is a convex combination, so a projected value CANNOT
# legitimately fall outside [min(source), max(source)] -- but terra >= 1.8 emits a
# large negative fill value (~ -2.7e37) wherever the interpolation window touches
# the source NA mask instead of propagating NA. Those are finite numbers, so
# `as.data.frame(..., na.rm = TRUE)` keeps them and they reach the models: without
# this clamp ~14k cells of vpi/cold/dry_mean carry -2.7e37 and the ALE grids,
# VIF and concurvity computed from them are meaningless.
align_to_template <- function(x, template, method = "bilinear") {
  out <- terra::project(x, template, method = method)
  lo <- suppressWarnings(min(terra::global(x, "min", na.rm = TRUE)[, 1], na.rm = TRUE))
  hi <- suppressWarnings(max(terra::global(x, "max", na.rm = TRUE)[, 1], na.rm = TRUE))
  if (is.finite(lo) && is.finite(hi)) out[out < lo | out > hi] <- NA
  out
}

load_large_herbivore_traits <- function(traits_csv) {
  readr::read_csv(traits_csv, show_col_types = FALSE) %>%
    filter(
      Mass.g > 5000,
      Terrestrial == 1,
      Diet.Plant > 0,
      !Order.1.2 %in% c("Carnivora", "Dasyuromorphia")
    ) %>%
    mutate(species_name = gsub("_", " ", Binomial.1.2))
}

# Every megafauna layer is built the same way: sum the species ranges in the source CRS,
# then bilinear-reproject once. Because bilinear resampling is linear it commutes with
# summation, so summing the clade layers reproduces total megafauna richness EXACTLY --
# Model 1's `mammal_richness` is derived that way below rather than computed separately.
build_clade_richness_stack <- function(large_mammals, ranges_dir, template, output_dir) {
  clades <- sort(unique(large_mammals$Order.1.2))
  all_tifs <- list.files(ranges_dir, "\\.tif$", full.names = TRUE, recursive = TRUE)
  tif_key <- tolower(basename(all_tifs))

  clade_layers <- vector("list", length(clades))
  names(clade_layers) <- clades

  for (clade in clades) {
    species <- large_mammals %>% filter(Order.1.2 == clade) %>% pull(species_name)
    idx <- match(tolower(paste0(gsub(" ", "_", species), ".tif")), tif_key)  # first hit per species
    sp_paths <- all_tifs[idx[!is.na(idx)]]
    message("Building clade richness layer: ", clade,
            " (", length(sp_paths), "/", length(species), " ranges found)")

    if (length(sp_paths) == 0) {
      clade_raster <- template
      terra::values(clade_raster) <- 0
    } else {
      # NB plain sum(): terra provides it through the S4 Summary group generic,
      # `terra::sum` is not an exported object.
      clade_raster <- align_to_template(
        sum(terra::rast(sp_paths), na.rm = TRUE), template, method = "bilinear")
    }

    names(clade_raster) <- clade
    clade_layers[[clade]] <- clade_raster
    terra::writeRaster(clade_raster, file.path(output_dir, paste0("richness_", clade, ".tif")), overwrite = TRUE)
  }

  # Absence is kept as 0. Writing it as NA would make "no species of this clade here"
  # indistinguishable from "no data" for anyone reading Clade_richness_stack.tif from
  # Dryad; 06_model2.R needs the zeros and would otherwise have to invent them back.
  terra::rast(clade_layers)
}

# ---- Main --------------------------------------------------------------------
message("Preparing predictor stack...")

trunk_richness <- terra::rast(require_file(paths$rel_richness, "relative trunk richness raster"))
names(trunk_richness) <- "trunkrichness"

vpi_stack <- terra::rast(require_file(file.path(paths$env_preds_dir, "vpi_stack.nc"), "VPI stack"))
cold_stack <- terra::rast(require_file(file.path(paths$env_preds_dir, "cold_stack.nc"), "cold-season deciduousness stack"))
dry_stack <- terra::rast(require_file(file.path(paths$env_preds_dir, "dry_stack.nc"), "dry-season deciduousness stack"))
height_nc <- terra::rast(require_file(file.path(paths$env_preds_dir, "height.nc"), "vegetation height raster"))

vpi_mean <- mean(align_to_template(vpi_stack, trunk_richness), na.rm = TRUE)
cold_mean <- mean(align_to_template(cold_stack, trunk_richness), na.rm = TRUE)
dry_mean <- mean(align_to_template(dry_stack, trunk_richness), na.rm = TRUE)
height <- align_to_template(height_nc, trunk_richness)

large_mammals <- load_large_herbivore_traits(require_file(paths$mammal_traits, "PHYLACINE traits"))

# Clade layers first: Model 1's total megafauna richness is their sum, so the two models
# share one definition of "megafauna richness" by construction.
message("Preparing mammal clade richness stack...")
clade_stack <- build_clade_richness_stack(
  large_mammals = large_mammals,
  ranges_dir    = paths$mammal_ranges,
  template      = trunk_richness,
  output_dir    = paths$raster_model
)
terra::writeRaster(clade_stack, paths$clade_stack, overwrite = TRUE)

mammal_richness <- sum(clade_stack, na.rm = TRUE)

consumption_list <- readRDS(require_file(paths$consumption_rds, "mammal consumption RDS"))
consumption <- terra::rast(consumption_list[["pn.consumption.maps"]][[1]])
consumption <- align_to_template(consumption, trunk_richness)

predictor_stack <- c(
  trunk_richness,
  vpi_mean,
  cold_mean,
  dry_mean,
  height,
  mammal_richness,
  consumption
)

names(predictor_stack) <- c(
  "trunkrichness",
  "vpi_mean",
  "cold_mean",
  "dry_mean",
  "height",
  "mammal_richness",
  "consumption"
)

# Guard: the models read this stack with na.rm = TRUE, which drops NA but happily
# keeps any finite fill value. Fail loudly here rather than train on nonsense.
guard <- as.data.frame(predictor_stack, na.rm = TRUE)
offending <- names(guard)[vapply(guard, function(v) any(abs(v) > 1e30), logical(1))]
if (length(offending))
  stop("predictor_stack carries fill-like values in: ", paste(offending, collapse = ", "),
       "\n  (check align_to_template's range clamp)", call. = FALSE)
message(sprintf("predictor stack OK: %d complete cells, all layers within source range", nrow(guard)))

terra::writeRaster(predictor_stack, paths$predictor_stack, overwrite = TRUE)

model_df <- as.data.frame(predictor_stack, xy = TRUE) %>%
  rename(longitude = x, latitude = y) %>%
  filter(!is.na(trunkrichness))
readr::write_csv(model_df, file.path(paths$raster_model, "predictor_table.csv"))

message("Data preparation completed successfully.")
