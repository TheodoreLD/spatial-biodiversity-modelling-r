# Environment of record

The full pipeline was run end to end and verified in this environment.
Package versions matter here: see the note on `terra` below.

- **R**: R version 4.5.2 (2025-10-31 ucrt)
- **Platform**: x86_64-w64-mingw32

| package | version |
|---|---|
| `here` | 1.0.2 |
| `terra` | 1.8.93 |
| `sf` | 1.1.0 |
| `dplyr` | 1.2.0 |
| `tidyr` | 1.3.2 |
| `purrr` | 1.2.1 |
| `readr` | 2.2.0 |
| `stringr` | 1.6.0 |
| `ggplot2` | 4.0.2 |
| `patchwork` | 1.3.2 |
| `scales` | 1.4.0 |
| `tibble` | 3.3.1 |
| `readxl` | 1.4.5 |
| `TNRS` | 0.3.6 |
| `catboost` | 1.2.10 |
| `raster` | 3.6.32 |
| `sp` | 2.2.1 |

`raster` and `sp` appear in no `library()` call. They are required all the same: the
mammal-consumption RDS (`data/raw/Map_mammal_consumption/Map.rds`) stores a `RasterLayer`,
so `readRDS` needs that S4 class defined before `terra::rast()` can coerce it in `04`.

## terra version sensitivity (important)

`04_data_preparation.R` reprojects the vegetation netCDFs onto the analysis grid with
`terra::project(..., method = "bilinear")`. On terra >= 1.8 this emits a large negative
fill value (~ -2.7e37) wherever the interpolation window touches the source's NA mask,
instead of propagating NA. Those are finite numbers, so `as.data.frame(..., na.rm = TRUE)`
keeps them and they reach the models silently.

`align_to_template()` therefore clamps each reprojected layer to that layer's own range
(bilinear resampling is a convex combination, so it cannot legitimately leave that range),
and `04` refuses to write a stack that still carries fill-like values. Do not remove either
guard: without them ~14k cells of `vpi_mean` / `cold_mean` / `dry_mean` carry -2.7e37 and
the ALE grids, VIF and concurvity computed from them are meaningless.

Two details of that guard matter. The clamp is **per layer**, not one global range across
the stack -- a shared range would clip `height` to the range of `vpi_mean`. And it is
skipped entirely for `method = "near"`, which copies source values rather than interpolating
between them, so a nearest-neighbour layer can legitimately hold any value the source holds.
Clade richness is resampled with `near` for exactly that reason: it is a count, and bilinear
interpolation of counts invents non-integer richness at every cell it touches.

## Reproducibility notes

- Steps 1-3 are deterministic and reproduce byte-identical outputs.
- **Step 5 is too**: two consecutive full runs of `05_model1.R` on one stack produced all
  20 outputs byte-identical, PNGs included, and reading plus aggregating the stack agrees
  to 17 significant figures across separate R sessions. Model 1 is bit-reproducible given
  its input.
- **Step 4 is not.** Running `04_data_preparation.R` twice on the same machine, same R,
  same inputs, produces predictor stacks that differ. Confirmed by checksumming the
  in-memory `predictor_stack` in the session that built it, so no file is written or read
  in between and no sync race can explain it: across one pair `dry_mean` differed by 7,449
  NA cells and 2,382 in its sum, while `trunkrichness`, `cold_mean`, `height`,
  `mammal_richness` and `consumption` were identical. In another pair it was `cold_mean`
  that moved, by 3,632 cells. Across runs it
  is a different vegetation layer each time (`vpi_mean`, `cold_mean` or `dry_mean`), always
  at coverage edges, and the affected layer gains or loses a few thousand NA cells. The
  cause is the same terra behaviour documented above: `project()` treats the source NA mask
  differently depending on whether the raster is held in memory or file-backed, and terra
  makes that choice from free memory at run time. `trunkrichness`, `height`,
  `mammal_richness`, `consumption` and the whole clade stack are stable.
  The instability has two scales. Two `project()` calls on the same object in the same
  session differ in the last float32 bit (max 5.96e-08 = 2^-24) in a few dozen cells; that
  is always present and harmless. Separately, and less often, one layer shifts its NA
  boundary by a few thousand cells, which is the effect that matters.
  Consequence: rebuilding the stack moves Model 1's pooled OOF R2 by ~0.004, well inside
  the 0.025 spread across fold allocations, so no conclusion depends on it -- but the
  archive cannot promise a byte-identical `predictor_stack.tif`. Use the Dryad copy if you
  need to reproduce the exact numbers in the README. Forcing `terraOptions(todisk = TRUE)`
  before step 4 would likely pin the behaviour to one branch; this has not been tested.
- A separate hazard, worth knowing if you work in a synced folder: this project lives under
  OneDrive, and a 48 MB `predictor_stack.tif` read while a sync is still in flight yields a
  run that does not reproduce. One Model 1 run was lost that way before three later sessions
  agreed with each other. Let the sync settle before starting step 5, or work outside the
  synced tree.
- Models 1 and 2 seed every resampling step, so reruns reproduce within floating point.
- CatBoost is not on CRAN; install the release binary matching the version above.
- The 1 deg cell count depends on the terra version through NA-boundary handling
  (4,806 cells here). Small differences are expected across terra releases.

_Generated 2026-08-04; package list and `align_to_template` notes revised 2026-08-21._
