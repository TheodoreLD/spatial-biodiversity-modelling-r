# Data

## Included here
- `spiny_trunk_species.xlsx` — the compiled list of woody species bearing trunk spines,
  with columns: `Species`, `Family`, `Syndrome` (prickly / thorny), `Range`
  (presence of a mappable range), `Native` (continent of origin), `Source`, `Picture`.
  95 species compiled; 74 have a mappable range (46 prickly, 28 thorny) and are the
  ones every analysis uses.
- `assets/world_borders.geojson` — country outlines for the maps (exported from R's
  `maps` package).
- `processed/species_standardized.csv` — the WFO/WCVP + TNRS resolution of the 95 names,
  written by `R/00_standardization.R`. Documentation of nomenclature only: the spatial
  pipeline matches range files by filename and does not read this table.

## Archived on Dryad (not included here — large)
The raw and derived spatial layers are archived on Dryad:

> DOI: 10.5061/dryad.jdfn2z3pq  *(spiny-trunk trees dataset, rasters, and scripts)*

Unpack them under `data/raw/` (and `outputs/` for the derived stacks) so they match the
default entries in `R/config.R`:

| unpack to | contents | source |
|---|---|---|
| `data/raw/Range_spiny_species/` | per-species range polygons (GeoPackage), in the `SpeciesTotal/`, `SpeciesPrickly/` and `SpeciesThorny/` subfolders — unpack `TreeSpecies.zip` here | ranges from Boonman et al. 2024 |
| `data/raw/Range_mammals/` | PHYLACINE present-natural mammal ranges (`Ranges/Present_natural/`) + traits (`Traits/Trait_data.csv`) | Faurby et al. 2018, 2020 |
| `data/raw/Map_environmental_predictors/` | `vpi_stack.nc`, `cold_stack.nc`, `dry_stack.nc`, `height.nc` | Higgins et al. 2016; Simard et al. 2011 |
| `data/raw/Map_mammal_consumption/Map.rds` | present-natural megafauna consumption maps (`pn.consumption.maps`) | Pedersen et al. 2023 |
| `data/raw/Map_global_plant_richness/sr_Ensemble_Prediction_7774_LongLat.RData` | total woody-plant richness (`predictions_grid`) | Cai et al. 2023 |
| `data/raw/Ecoregion/Ecoregions2017.shp` | ecoregions + biomes + realms | Dinerstein et al. 2017 |
| `outputs/biogeography/rasters/` | absolute & relative richness (total / prickly / thorny), rebuilt by `01_biogeography.R` | this study |
| `outputs/raster_model/predictor_stack.tif` | model-ready predictor + response stack, rebuilt by `04_data_preparation.R` | this study |
| `outputs/raster_model/Clade_richness_stack.tif` | per-clade present-natural richness, rebuilt by `04_data_preparation.R` | this study |

The two `outputs/` entries are conveniences — `01_biogeography.R` and
`04_data_preparation.R` regenerate them from the raw layers.

## Pointing the scripts at the data
All input/output paths live in `R/config.R` (the `paths` list, `here()`-based) and every
script sources it. If your download sits elsewhere, edit the matching entries there —
no path is hardcoded in the analysis scripts.
