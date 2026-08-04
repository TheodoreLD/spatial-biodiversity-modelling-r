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
| 0 | `00_standardization.R` | standardise species names against WFO + WCVP with synonym resolution via **TNRS** (needs internet). Documentation only — the spatial pipeline matches range files by filename and does not read its output | `data/processed/` |
| 1 | `01_biogeography.R` | rasterise ranges ÷ woody richness → absolute + relative richness (total/prickly/thorny) | `outputs/biogeography/` |
| 2 | `02_descriptive.R` | richness maps (country borders, no axes), biome & ecoregion distributions (block-bootstrap CIs + rank-stability), Fig-3 composite | `outputs/descriptive/` |
| 3 | `03_mammal_ecoregions.R` | present-natural megafauna (>5 kg) of the richest ecoregion per continent | `outputs/ecoregion_mammal/` |
| 4 | `04_data_preparation.R` | model-ready predictor stack + per-clade richness stack | `outputs/raster_model/` |
| 5 | `05_model1.R` | **Model 1** — environmental characterisation via **ALE** (no importance ranking); collinearity + residual diagnostics | `outputs/model1/` |
| 6 | `06_model2.R` | **Model 2** — mammal-clade associations, each clade fitted **with and without** the environmental covariates: unconstrained shape + OOB R², positive-monotonic ranking slope + OOB R², plus variance-partitioning diagnostics | `outputs/model2/` |

Each model has companion scripts (reviewer responses / methods justification):

| script | question it answers |
|---|---|
| `05_model1_sensitivity.R` | why 1°, why no transform, why no spatial term |
| `05_model1_jackknife.R` | are the ALE shapes and CV-R² stable when 30% of the 74 species are dropped? |
| `06_model2_jackknife.R` | does the **clade ranking** survive dropping 10/20/30% of the species? |

`R/config.R` holds `here()`-based paths; `R/requirements.R` lists packages;
`ENVIRONMENT.md` records the versions the pipeline was verified with — read its `terra`
note before changing `04_data_preparation.R`.

## The three "conditions"
Every descriptive and Model-2 analysis runs for **total** (all 74 mapped species),
**prickly** (46, phellogen-derived, Neotropical-centred), and **thorny** (28, bud-derived,
Afro-Asian-centred). The response is **relative richness** = spiny-trunk species ÷ total
woody-plant richness (partly controls for uneven botanical sampling).

## Method notes
- **1° analysis unit** (aggregated from 0.1°): defensible *n* (4,858 cells) instead of
  ~404k pseudo-replicates from 74 species ranges.
- **The response is conditioned on presence.** Cells with zero spiny-trunk species are
  set to `NA` in `01_biogeography.R`, so both models describe variation *among occupied
  cells*. They do not model presence/absence, and absences are not treated as zeros.
- **Spatial** block cross-validation / block bootstrap for all performance and CIs — never
  cell-level resampling (neighbouring cells are strongly autocorrelated).
- **Model 1** characterises *where* richness is high (ALE response shapes), not which
  predictor is "most important"; the two megafauna metrics are read as one axis
  (concurvity ≈ 0.85) and their redundancy is disclosed, not ranked.
- **Model 2** fits each clade independently, so clade–clade correlation is not a
  confounder, and fits every clade **twice** — alone, and alongside the four
  environmental covariates. The two settings share the same per-iteration block-bootstrap
  seeds, so the comparison is *paired*: the change in shape, R² and slope is the
  environment's contribution, not resampling noise. Continents are split at the
  Suez/Red Sea line, so East Africa is African (`cont_of`).
- **The ranking slope is one-sided.** The monotonic model is constrained non-decreasing,
  so the slope is ≥ 0 by construction and ranks the *strength of positive association*.
  It is not a test of direction, and `slope_lo > 0` should not be read as one — the
  unconstrained response shape is what shows whether a relationship is positive.
- **Two slope columns, two jobs.** Each continent is modelled on its own cells, so both
  the predictor and the response are standardised *within* continent — 1 SD of
  Proboscidea richness is 0.23 species in Africa but 0.71 in Asia. Use `slope` to rank
  clades **within** a continent, and `slope_std` (= slope ÷ that continent's response SD,
  hence unit-free) to compare the **same clade between** continents. `slope_env_std` is
  the environment-adjusted equivalent.
- **How firm is the ranking?** Two independent checks, because the order turned out to be
  the fragile part: `model2_rank_stability.csv` re-ranks the clades inside every spatial
  block-bootstrap resample, and `06_model2_jackknife.R` re-ranks them after dropping
  10/20/30% of the 74 species. Agreement between the two is the evidence that an ordering
  is real; disagreement means it is not identifiable.
- Residual short-range spatial autocorrelation persists in both models and is disclosed
  as a limitation. Adding an explicit lon/lat term *raises* CV R² (see
  `spatial_term_test.csv`) but roughly halves the permutation importance of both
  predictor groups — it absorbs the signal being characterised, which is why Model 1
  omits it. The better fit is reported rather than hidden.
- **One definition of megafauna richness across both models**: every clade layer is the
  sum of its species ranges in the source CRS, bilinear-reprojected once, and Model 1's
  `mammal_richness` is the sum of those clade layers — so the two models cannot drift
  apart (`04_data_preparation.R`).

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
   the default `data/raw/` and `outputs/` locations, edit the entries there. Steps must
   run in order: 03 reads the ecoregion ranking that 02 writes, and 05/06 read the
   stacks that 04 writes.

Run from the project root (or with the RStudio project open) so `here()` resolves.
Country borders for the maps are in `data/assets/world_borders.geojson` (exported from R's
`maps` package). Committed outputs are figures (`.png`) and summary tables (`.csv`);
rasters (`.tif`) and the large per-cell `predictor_table.csv` are regenerated / on Dryad.

## Headline numbers (this archive)
Reported so the manuscript and the archive cannot drift apart:

| quantity | value | where |
|---|---|---|
| Model 1 pooled out-of-fold R² | **0.426** | `model1/performance_summary.csv` |
| — mean of per-fold R² (weaker statistic) | 0.367 (SD 0.162) | same |
| — random-fold reference | 0.509 (SD 0.037) | same |
| Analysis unit | 1° → 4,806 cells (from 398,474 at 0.1°) | `model1/resolution_sweep.csv` |
| Block size (residual range ~675 km) | 7° | `model1/residual_correlogram.csv` |
| Residual Moran's I at 75 km | 0.735 | same |
| Species-drop stability (30% dropped) | CV-R² 0.377 vs 0.367 full | `model1/jackknife.csv` |
| Clade models, cells per continent | Americas 2,062 · Africa 1,414 · Asia 1,277 | `model2/` |

Prefer **pooled OOF R²** to the mean of per-fold R²: per-fold R² divides by each fold's own
variance, so a fold covering a homogeneous region scores badly even when its RMSE is low
(here RMSE varies ~1.2× across folds while sd(y) varies ~1.4×). The fold spread is mostly a
denominator artefact, not regional model failure.

## Citation & licence
Cite the associated article (in revision, *Global Ecology and Biogeography*).
Data: Dryad (DOI in `data/README.md`). Add a licence file before public release.
