# Global ecology & evolutionary drivers of spiny-trunk trees — reproducibility archive (R)

A single **all-R** pipeline for the study of the global distribution, ecological
associations, and megafaunal-herbivory drivers of woody plants that bear spines on their
**trunks**. It compiles a dataset of spiny-trunk species, builds global richness layers,
and analyses them with gradient-boosted regression trees (**CatBoost**) under spatially
honest cross-validation.

Everything — data preparation, models, and figures — runs in R.

---

## Pipeline (`R/`)

Descriptive precedes the modelling — the models (and the mammal-ecoregion step) build on
the richness layers and the richest ecoregions that the descriptive step characterises.

| step | script | does | writes to |
|---|---|---|---|
| 0 | `00_standardization.R` | standardise species names against WFO + WCVP with synonym resolution via **TNRS** (needs internet) | `data/processed/` |
| 1 | `01_biogeography.R` | rasterise ranges ÷ woody richness → absolute + relative richness (total/prickly/thorny) | `outputs/biogeography/` |
| 2 | `02_descriptive.R` | richness maps (country borders, no axes), biome & ecoregion distributions (block-bootstrap CIs + rank-stability), Fig-3 composite | `outputs/descriptive/` |
| 3 | `03_mammal_ecoregions.R` | present-natural megafauna (>5 kg) of the richest ecoregion per continent | `outputs/ecoregion_mammal/` |
| 4 | `04_data_preparation.R` | model-ready predictor stack + per-clade richness stack | `outputs/raster_model/` |
| 5 | `05_model1.R` | **Model 1** — environmental characterisation via **ALE** (no importance ranking); collinearity + residual diagnostics | `outputs/model1/` |
| 6 | `06_model2.R` | **Model 2** — univariate mammal-clade associations: monotonic slope + unconstrained shape + OOB R² (both) + environment-adjusted **ALE** slope + variance-partitioning diagnostics | `outputs/model2/` |

Model 1 has two companion scripts (Reviewer responses / methods justification):
`05_model1_sensitivity.R` — resolution / transform / spatial-term sweeps (why 1°, why no
transform, why no spatial term); `05_model1_jackknife.R` — species-drop robustness (are
the ALE shapes and CV-R² stable when 30% of the 74 species are removed?).

`R/config.R` holds `here()`-based paths; `R/requirements.R` lists packages.

## The three "conditions"
Every descriptive and Model-2 analysis runs for **total** (all 74 mapped species),
**prickly** (46, phellogen-derived, Neotropical-centred), and **thorny** (28, bud-derived,
Afro-Asian-centred). The response is **relative richness** = spiny-trunk species ÷ total
woody-plant richness (partly controls for uneven botanical sampling).

## Method notes
- **1° analysis unit** (aggregated from 0.1°): defensible *n* (4,858 cells) instead of
  ~404k pseudo-replicates from 74 species ranges.
- **Spatial** block cross-validation / block bootstrap for all performance and CIs — never
  cell-level resampling (neighbouring cells are strongly autocorrelated).
- **Model 1** characterises *where* richness is high (ALE response shapes), not which
  predictor is "most important"; the two megafauna metrics are read as one axis
  (concurvity ≈ 0.85) and their redundancy is disclosed, not ranked.
- **Model 2** fits each clade independently (univariate), so clade–clade correlation is
  not a confounder; environmental confounding is handled by the env-adjusted ALE slope
  and variance partitioning.
- Residual short-range spatial autocorrelation persists in both models and is disclosed
  as a limitation (an explicit spatial term would erase the ecological signal being
  measured).

## Reproducing
1. **Install R packages** — `source("R/requirements.R")`.
2. **Install CatBoost for R** (not on CRAN): from the official release binary, e.g.
   `install.packages("https://github.com/catboost/catboost/releases/download/v1.2.10/catboost-R-windows-x86_64-1.2.10.tgz", repos=NULL)`
   (pick the file for your OS). This archive was built and verified with **CatBoost 1.2.10**.
3. **Get the heavy data** from **Dryad** (see `data/README.md` for the DOI): raw ranges,
   PHYLACINE mammal ranges, environmental rasters, the Ecoregions2017 shapefile, and the
   derived richness/predictor stacks. They are not in git (size).
4. **Run** `run_pipeline.R` (or the `R/` scripts in order). Every path lives in
   `R/config.R` (the `paths` list) and each script sources it; if your data isn't under
   the default `data/raw/` and `outputs/` locations, edit the entries there.

Country borders for the maps are in `data/assets/world_borders.geojson` (exported from R's
`maps` package). Committed outputs are figures (`.png`) and tables (`.csv`); rasters
(`.tif`) are regenerated / on Dryad.

## Citation & licence
Cite the associated article (in revision, *Global Ecology and Biogeography*).
Data: Dryad (DOI in `data/README.md`). Add a licence file before public release.
