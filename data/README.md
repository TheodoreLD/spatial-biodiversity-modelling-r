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
The Python analyses read the rasters through the `TRUNK_DATA` environment variable —
set it to the folder holding the Dryad `output/Raster_model/` rasters (and, for the
descriptive step, the `output/Biogeography/rasters/` + `Data/Ecoregion/` files).
