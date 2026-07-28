# Data

## Included here
- `spiny_trunk_species.xlsx` — the compiled list of woody species bearing trunk spines,
  with columns: `Species`, `Family`, `Syndrome` (prickly / thorny), `Range`
  (presence of a mappable range), `Native` (continent of origin), `Source`, `Picture`.

## Archived on Dryad (not included here — large)
The raw and derived spatial layers are archived on Dryad:

> DOI: 10.5061/dryad.jdfn2z3pq  *(spiny-trunk trees dataset, rasters, and scripts)*

They include:
- `Data/Range_spiny_species/` — per-species range polygons (GeoPackage).
- `Data/Range_mammals/` — PHYLACINE present-natural mammal ranges + traits.
- `Data/Map_environmental_predictors/` — vegetation/climate rasters (Higgins et al. 2016; Simard et al. 2011).
- `Data/Map_global_plant_richness/` — total woody-plant richness (Cai et al. 2023).
- `Data/Ecoregion/` — Ecoregions2017 (Dinerstein et al. 2017).
- `output/Raster_model/predictor_stack.tif` — the model-ready predictor + response stack.
- `output/Raster_model/Clade_richness_stack.tif` — per-clade present-natural richness.
- `output/Biogeography/rasters/` — absolute & relative richness (total / prickly / thorny).

## Pointing the scripts at the data
All input/output paths live in `R/config.R` (the `paths` list, `here()`-based) and every
script sources it. By default they resolve to `data/raw/` for the raw layers and
`outputs/` for the rasters that `01_biogeography.R` and `04_data_preparation.R` rebuild.
If your Dryad download sits elsewhere, edit the matching entries in `R/config.R`.
