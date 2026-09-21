# =============================================================================
# 04_data_preparation.R
# Prepare predictor and mammal clade richness stacks for downstream modeling.
# Writes the canonical stacks (predictor_stack.tif / Clade_richness_stack.tif) that
# 05_model1.R and 06_model2.R read (paths defined in R/config.R).
#
# LAYOUT
#   1  settings       megafauna definition, fill-value limit, environmental inputs
#   2  helpers        reprojection, coverage reporting, trait loading, clade stack
#   3  template       relative trunk richness defines the grid everything is aligned to
#   4  environment    productivity and deciduousness stacks reduced to per-cell means
#   5  megafauna      clade richness stack, total richness, consumption
#   6  assemble       stack, validate against fill values, write raster + table
#
# =============================================================================

suppressPackageStartupMessages({
  library(terra)
  library(dplyr)
  library(readr)
})

source("R/config.R")
create_project_dirs()


# =============================================================================
# 1. SETTINGS
# =============================================================================

# MASS_MIN_G and EXCLUDE_ORDERS define what counts as a large herbivore; they live in
# config.R because 03 applies the same filter.

# Anything beyond this magnitude is a fill value, not data. 
FILL_LIMIT <- 1e30

# Multi-layer environmental stacks, each reduced to a per-cell mean over its layers.
ENV_STACKS <- c(vpi  = "vpi_stack.nc",
                cold = "cold_stack.nc",
                dry  = "dry_stack.nc")


# =============================================================================
# 2. HELPERS
# =============================================================================

# align_to_template() lives in config.R: 05's ratio check needs it too.

report_coverage <- function(x, nm){
  n <- sum(!is.na(x))
  q <- terra::global(n, c("min","mean"), na.rm = TRUE)
  message(sprintf("  %-10s layers per cell: min %d, mean %.1f of %d",
                  nm, as.integer(q[1,"min"]), q[1,"mean"], terra::nlyr(x)))
  invisible(NULL)
}

load_large_herbivore_traits <- function(traits_csv) {
  readr::read_csv(traits_csv, show_col_types = FALSE) %>%
    filter(
      Mass.g > MASS_MIN_G,
      Terrestrial == 1,
      Diet.Plant > 0,
      !Order.1.2 %in% EXCLUDE_ORDERS,
      !genus_of(Binomial.1.2) %in% EXCLUDE_GENERA
    ) %>%
    mutate(species_name = gsub("_", " ", Binomial.1.2))
}

build_clade_richness_stack <- function(large_mammals, ranges_dir, template, output_dir) {
  clades <- sort(unique(large_mammals$Order.1.2))
  if (!dir.exists(ranges_dir))
    stop("Mammal ranges directory does not exist: ", ranges_dir, call. = FALSE)
  all_tifs <- list.files(ranges_dir, "\\.tif$", full.names = TRUE, recursive = TRUE)
  if (!length(all_tifs))
    stop("No mammal range rasters found in: ", ranges_dir, call. = FALSE)
  tif_key <- tolower(basename(all_tifs))
  dup <- unique(tif_key[duplicated(tif_key)])
  if (length(dup))
    stop("Duplicate mammal raster filenames under ", ranges_dir, ": ",
         paste(head(dup, 5), collapse = ", "), call. = FALSE)

  clade_layers <- vector("list", length(clades))
  names(clade_layers) <- clades

  for (clade in clades) {
    species <- large_mammals %>% filter(Order.1.2 == clade) %>% pull(species_name)
    idx <- match(tolower(paste0(gsub(" ", "_", species), ".tif")), tif_key)
    sp_paths <- all_tifs[idx[!is.na(idx)]]
    message("Building clade richness layer: ", clade,
            " (", length(sp_paths), "/", length(species), " ranges found)")

    if (length(sp_paths) == 0)
      stop(clade, ": none of its ", length(species),
           " species has a range raster under ", ranges_dir, call. = FALSE)
    if (anyNA(idx))
      warning(clade, ": missing ranges for ", sum(is.na(idx)), " of ", length(species),
              " species: ", paste(head(species[is.na(idx)], 5), collapse = ", "),
              if (sum(is.na(idx)) > 5) ", ..." else "", call. = FALSE)


    clade_raster <- align_to_template(
      sum(terra::rast(sp_paths), na.rm = TRUE), template, method = "near")

    names(clade_raster) <- clade
    clade_layers[[clade]] <- clade_raster
    terra::writeRaster(clade_raster, file.path(output_dir, paste0("richness_", clade, ".tif")), overwrite = TRUE)
  }

  terra::rast(clade_layers)
}


# =============================================================================
# 3. TEMPLATE — relative trunk richness defines the analysis grid
# =============================================================================
message("Preparing predictor stack...")

trunk_richness <- terra::rast(require_file(paths$rel_richness, "relative trunk richness raster"))
if (terra::nlyr(trunk_richness) != 1)
  stop("relative trunk richness must contain exactly one layer", call. = FALSE)
names(trunk_richness) <- "trunkrichness"


# =============================================================================
# 4. ENVIRONMENT — productivity and deciduousness, reduced to per-cell means
# =============================================================================

vpi_stack  <- terra::rast(require_file(file.path(paths$env_preds_dir, ENV_STACKS[["vpi"]]),  "VPI stack"))
cold_stack <- terra::rast(require_file(file.path(paths$env_preds_dir, ENV_STACKS[["cold"]]), "cold-season deciduousness stack"))
dry_stack  <- terra::rast(require_file(file.path(paths$env_preds_dir, ENV_STACKS[["dry"]]),  "dry-season deciduousness stack"))
height_nc  <- terra::rast(require_file(file.path(paths$env_preds_dir, "height.nc"), "vegetation height raster"))

vpi_aligned  <- align_to_template(vpi_stack, trunk_richness);  report_coverage(vpi_aligned,  "vpi")
vpi_mean <- mean(vpi_aligned, na.rm = TRUE)
cold_aligned <- align_to_template(cold_stack, trunk_richness); report_coverage(cold_aligned, "cold")
cold_mean <- mean(cold_aligned, na.rm = TRUE)
dry_aligned  <- align_to_template(dry_stack, trunk_richness);  report_coverage(dry_aligned,  "dry")
dry_mean <- mean(dry_aligned, na.rm = TRUE)

# Height is already a single layer, so there is nothing to average over.
height <- align_to_template(height_nc, trunk_richness)


# =============================================================================
# 5. MEGAFAUNA — clade richness stack, total richness, consumption
# =============================================================================
large_mammals <- load_large_herbivore_traits(require_file(paths$mammal_traits, "PHYLACINE traits"))
if (anyDuplicated(large_mammals$Binomial.1.2))
  stop("PHYLACINE traits contain duplicate Binomial.1.2 rows", call. = FALSE)

message("Preparing mammal clade richness stack...")
clade_stack <- build_clade_richness_stack(
  large_mammals = large_mammals,
  ranges_dir    = paths$mammal_ranges,
  template      = trunk_richness,
  output_dir    = paths$raster_model
)
terra::writeRaster(clade_stack, paths$clade_stack, overwrite = TRUE)

valid_clades <- sum(!is.na(clade_stack))
mammal_richness <- sum(clade_stack, na.rm = TRUE)
mammal_richness[valid_clades == 0] <- NA

consumption_list <- readRDS(require_file(paths$consumption_rds, "mammal consumption RDS"))

consumption <- terra::rast(consumption_list[["pn.consumption.maps"]][["mean"]])
if (terra::nlyr(consumption) != 1)
  stop("consumption raster must contain exactly one layer", call. = FALSE)
consumption <- align_to_template(consumption, trunk_richness)


# =============================================================================
# 6. ASSEMBLE — stack, validate against fill values, write raster and table
# =============================================================================
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

rng <- terra::global(predictor_stack, c("min", "max"), na.rm = TRUE)
offending <- rownames(rng)[abs(rng$min) > FILL_LIMIT | abs(rng$max) > FILL_LIMIT]
if (length(offending))
  stop("predictor_stack carries fill-like values in: ", paste(offending, collapse = ", "),
       "\n  (check align_to_template's fill strip and range clamp)", call. = FALSE)
message("predictor stack OK: no fill-like values in any layer")

terra::writeRaster(predictor_stack, paths$predictor_stack, overwrite = TRUE)

model_df <- as.data.frame(predictor_stack, xy = TRUE) %>%
  rename(longitude = x, latitude = y) %>%
  filter(!is.na(trunkrichness))
message(sprintf("  %d cells with trunk richness written to predictor_table.csv", nrow(model_df)))
readr::write_csv(model_df, file.path(paths$raster_model, "predictor_table.csv"))

message("Data preparation completed successfully.")
