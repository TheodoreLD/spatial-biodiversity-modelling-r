# Data

## Included here
- `spiny_trunk_species.xlsx` — the compiled list of woody species bearing trunk spines,
  with columns: `Species`, `Family`, `Syndrome` (prickly / thorny), `Range`
  (presence of a mappable range), `Native` (continent of origin), `Source`, `Picture`.
- `CHANGES.txt` — corrections applied to the submitted species list (genus-spelling fix,
  a missing continent, category standardisation) and the reconciled continent × syndrome
  cross-tab.

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
`run_pipeline.R` sets `TRUNK_*` environment variables that point the descriptive step
and the models at the derived rasters — `TRUNK_STACK` / `TRUNK_CLADE` (the
`Raster_model/` stacks), `TRUNK_REL_DIR` (the `Biogeography/rasters/` folder),
`TRUNK_ECO` (the Ecoregions2017 shapefile) and `TRUNK_BORDERS`
(`data/assets/world_borders.geojson`). Edit the paths at the top of `run_pipeline.R`
to match where your rasters live — either the Dryad `output/` layout or the
`outputs/` folders rebuilt by `01_biogeography.R` and `04_data_preparation.R`.
