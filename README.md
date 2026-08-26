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
| 2 | `02_descriptive.R` | richness maps (country borders, no axes), biome distributions (area-weighted IQR) and ecoregion means (block-bootstrap CIs), Fig-3 composite | `outputs/descriptive/` |
| 3 | `03_mammal_ecoregions.R` | present-natural megafauna (>5 kg) of the richest ecoregion per continent, plus the best CI-backed ecoregion where the richest one is too thinly replicated for `02` to bound | `outputs/ecoregion_mammal/` |
| 4 | `04_data_preparation.R` | model-ready predictor stack + per-clade richness stack | `outputs/raster_model/` |
| 5 | `05_model1.R` | **Model 1** — environmental characterisation via **ALE** (no importance ranking), plus its design justifications (resolution, transform) and species-drop robustness, all in one file; sections selectable with `M1_RUN` | `outputs/model1/` |
| 6 | `06_model2.R` | **Model 2** — mammal-clade associations, each clade fitted **with and without** the environmental covariates: unconstrained shape (direction) + OOB R², constrained shape + `slope_std` ranking + OOB R², plus unconstrained variance-partitioning diagnostics | `outputs/model2/` |

Each model has companion scripts (reviewer responses / methods justification):

| script | question it answers |
|---|---|
| `06_model2_jackknife.R` | does the **slope ranking** survive dropping 10/20/30% of the species (30 replicates each, both settings)? Shapes and R² are not refit |
| `06_model2_paper_figure.R` | assembles the combined manuscript figure (jackknife + unconstrained shapes); fits nothing, runs in seconds |

`R/config.R` holds `here()`-based paths; `R/requirements.R` lists packages;
`ENVIRONMENT.md` records the versions the pipeline was verified with — read its `terra`
note before changing `04_data_preparation.R`.

## The three "conditions"
Every descriptive and Model-2 analysis runs for **total** (all 74 mapped species),
**prickly** (46, phellogen-derived, Neotropical-centred), and **thorny** (28, bud-derived,
Afro-Asian-centred). The response is **relative richness** = spiny-trunk species ÷ total
woody-plant richness (partly controls for uneven botanical sampling).

## Method notes
- **1° analysis unit** (aggregated from 0.1°): defensible *n* (4,806 cells) instead of
  ~399k pseudo-replicates from 74 species ranges. That is the reason for the choice; the
  resolution sweep is a check, not the justification — it cannot separate 1° (0.357 ± 0.032)
  from 2° (0.373 ± 0.011), and only rules out 0.1° (0.289 ± 0.018) firmly.
- **Tuning does not discriminate.** All 12 configurations fall within one standard error of
  the minimum in every fold allocation, so the 1-SE rule returns the simplest corner of the
  grid by construction rather than by selection. Hyperparameters are also chosen on the same
  cells that produce the reported R² (no nested CV); with a grid this flat the leak is small,
  but the estimate is post-selection.
- **The response is conditioned on presence.** Cells with zero spiny-trunk species are
  set to `NA` in `01_biogeography.R`, so both models describe variation *among occupied
  cells*. They do not model presence/absence, and absences are not treated as zeros.
- **Spatial** block cross-validation / block bootstrap for all performance and CIs — never
  cell-level resampling (neighbouring cells are strongly autocorrelated). Model 1 derives
  **8°** from its residual autocorrelation range (~675 km, converted at the latitude where
  most cells sit rather than a flat 111 km/deg). Model 2 still uses a **hardcoded 7°**
  (`06_model2.R:61`, `:266`) that predates that derivation — the two have not been
  reconciled. The descriptive summaries use 2°, from the much shorter within-ecoregion
  range (125–175 km, `descriptive/within_ecoregion_correlogram.csv`).
- **Model 1** characterises *where* richness is high (ALE response shapes), not which
  predictor is "most important"; the two megafauna metrics are read as one axis
  (concurvity 0.84 each) and their redundancy is disclosed, not ranked. The entanglement is
  not specific to them: all six predictors score 0.67–0.84, with productivity and vegetation
  height both at 0.77. Note also that concurvity is measured under *random* 3-fold CV, not
  spatial blocks — deliberately, since blocking would confound predictor overlap with
  geography, but it means the values are optimistic as absolute numbers and should be read
  as a ranking of entanglement rather than as calibrated R².
- **Model 2** fits each clade independently, so clade–clade correlation is not a
  confounder, and fits every clade **twice** — alone, and alongside the four
  environmental covariates. The two settings share the same per-iteration block-bootstrap
  seeds, so the comparison is *paired*: the change in shape, R² and effect size is the
  environment's contribution, not resampling noise. Continents are split at the
  Suez/Red Sea line, so East Africa is African (`cont_of`).
- **The ranking metric is `slope_std`**, the standardised slope of the constrained
  response: the least-squares slope of the fitted curve against z-scored clade richness,
  ÷ sd(response). Its units are **SDs of relative spiny-trunk richness per SD of clade
  richness** — standardised on both axes, which is what makes it comparable across clades
  (whose richness scales differ) and across continents (whose response scales differ). It
  is fitted at the distinct values of a 21-point quantile grid, so it is evaluated across
  the range the data span; where a clade takes few distinct richness values these collapse
  to fewer points, each contributing equally. sd(response) is the full-sample value, so it
  is a per-continent constant and ranking by the raw slope and by `slope_std` cannot
  disagree. `slope_env` is the environment-adjusted equivalent.
- **`slope_std` is one-sided.** The monotonic model is constrained non-decreasing, so the
  slope is ≥ 0 by construction and ranks the *strength of positive association*. It is not
  a test of direction: the direction of a response is read from the unconstrained curves
  (`model2_curves*.csv`, `FIG_M2_shape_unconstrained_*`), which are free to fall. A clade
  that fits flat
  (`slope_std` ≈ 0, e.g. Diprotodontia) is reported as such — no positive association is a
  result, not a reason for exclusion — with `n_present` alongside so sparse occupancy stays
  visible.
- **Variance partitioning is deliberately unconstrained.** It asks how much variance a
  clade explains uniquely versus jointly with environment; constraining the fit would
  measure the constraint as much as the data.
- **How firm is a slope?** Two independent perturbations. The spatial block bootstrap
  (200 resamples) gives the confidence interval on `slope_std`, reported for both settings.
  `06_model2_jackknife.R` then drops 10/20/30% of the 74 mapped species (30 replicates per
  fraction), re-derives relative richness and refits the slope — again for both settings.
  The first asks whether the result depends on *where* you sample, the second on *which
  species* were mapped; the bootstrap cannot see the latter.
- **Neither perturbation measures confounding.** How much of a clade's association is
  separable from environment is read from the variance partition
  (`model2_diagnostics.csv`, `FIG_M2_r2_partition_*`): `unique_clade` is what the clade
  adds on top of the environmental model, and it is negative where adding the clade makes
  the out-of-fold fit worse.
- Residual short-range spatial autocorrelation persists in both models and is disclosed
  as a limitation. It is handled by resampling whole spatial blocks rather than by adding
  a coordinate term to the model.
- **One definition of megafauna richness across both models**: every clade layer is the
  sum of its species ranges in the source CRS, reprojected once by **nearest neighbour**
  (these are counts of 0/1 presences and PHYLACINE is coarser than the 0.1° template, so
  bilinear would turn a count into a fraction), and Model 1's `mammal_richness` is the sum
  of those clade layers — so the two models cannot drift apart (`04_data_preparation.R`).

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
| Model 1 pooled out-of-fold R² | **0.374** (SD 0.025 over 10 fold allocations, range 0.331–0.400) | `model1/performance_summary.csv` |
| — training R² / overfit gap | 0.564 / 0.191 | same |
| — RMSE | 0.00148 | same |
| — mean of per-fold R² (weaker statistic) | 0.327 (SD 0.157); fold 2.5–97.5% −0.096 to 0.562 | same |
| Analysis unit | 1° → 4,806 cells (from 398,497 at 0.1°) | `model1/resolution_sweep.csv` |
| Model 1 block size | **8°**, from a 675 km residual range at 91.7 km per degree of longitude (\|lat\| p90 = 34.5) | `model1/blocksize_recheck.csv` |
| — confirmatory CV at the block the final residuals imply | 9° (825 km) gives 0.380, against 0.374 at 8° | `model1/blocksize_confirmatory.csv` |
| Residual Moran's I at 75 km | 0.744 (SD 0.027 over 10 allocations × 5 subsamples); below 0.1 beyond 825 km | `model1/residual_correlogram.csv` |
| Model 1 calibration slope | 0.996 (SD 0.075 over 10 allocations, range 0.914–1.161) | `model1/calibration.csv` |
| Model 1 final fit | 804 trees, chosen by early stopping and refitted on all 4,806 cells | `model1/performance_summary.csv` |
| Species-drop stability (30% dropped) | 0.352 ± 0.078 over 300 replicates vs 0.388 with the full pool (52 of 74 species kept; pooled OOF, 8° blocks, one fixed fold allocation) | `model1/jackknife.csv` |
| Model 2 block size | **7°**, hardcoded rather than derived — it predates Model 1's derivation and has not been reconciled | `R/06_model2.R:61` |
| Descriptive block bootstrap | 2°; 2,000 resamples per ecoregion. Biomes are **not** bootstrapped: their bar is an area-weighted IQR across cells, not a CI | `descriptive/ecoregion_summary_*.csv` |
| Within-ecoregion range (block size source) | 175 km spiny and prickly, 125 km thorny → 2° | `descriptive/within_ecoregion_correlogram_summary.csv` |
| Ecoregions keeping a CI | 302 of 462 (total); the rest hold fewer than 3 effective 2° blocks | `descriptive/ecoregion_summary_total.csv` |
| Richest ecoregion per continent | Beni savanna (4.9 effective blocks) · Cross-Niger transition forests (1.4, **no CI** — Sudd flooded grasslands, 6.1, reported alongside) · Tonle Sap freshwater swamp forests (3.7) | `ecoregion_mammal/ecoregion_mammals.csv` |
| — megafauna species in each | Beni 40 · Cross-Niger 17 vs Sudd 34 · Tonle Sap 19 | same |
| Clade models, cells per continent | Americas 2,062 · Africa 1,414 · Asia 1,277 | `model2/` |

Prefer **pooled OOF R²** to the mean of per-fold R²: per-fold R² divides by each fold's own
variance, so a fold covering a homogeneous region scores badly even when its RMSE is low
(here RMSE varies ~1.2× across folds while sd(y) varies ~1.4×). The fold spread is mostly a
denominator artefact, not regional model failure.

## Citation & licence
Cite the associated article (in revision, *Global Ecology and Biogeography*).
Data: Dryad (DOI in `data/README.md`). Add a licence file before public release.
