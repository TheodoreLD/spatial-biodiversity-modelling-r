# Repository review — 2026-11-16

Review of repository content at commit `e1f19337a08273da40ec2598cc3b19474e045c6e`,
with a documentation and figure-presentation update. This is a static archive review,
not an independent rerun or scientific validation of the models.

## Verified

- The archive contains 47 PNG figures and 47 output CSVs. The output index links all 94 files.
- The three figures featured in the main README were opened and visually inspected:
  descriptive total-condition composite, Model 1 ALE and Model 2 total-condition composite.
- All literal `source()` file references in the 10 R files resolve within the repository.
- Model 1's performance table supports the rounded headline values:
  pooled OOF R² 0.372, SD 0.019, range 0.338–0.395, training R² 0.565,
  mean per-fold R² 0.324 and RMSE 0.00148. The calibration table supports slope 0.999.
- The resolution table contains 4,806 cells at 1°, compared with 398,474 at 0.1°.
- Model 2 coverage has 90 rows, with 56 estimable combinations. The environment-results
  table has 56 rows; 3 clade-increment intervals exclude zero and 44 shared components
  are negative.
- The three Model 2 R² components sum to the joint-model R² within floating-point
  precision (maximum absolute discrepancy approximately 1.1e-16).
- Input paths, pipeline ordering and the package list are documented. Heavy spatial
  inputs are explicitly excluded from Git and must be supplied externally.

## Documentation corrections

- Added main figures, explanatory captions and links to their scripts and result tables.
- Added a complete output index and clickable headline evidence links.
- Distinguished a full rebuild from modelling with the archived stacks, highlighting the
  existing step-4 reproducibility limitation.
- Clarified that fitting clades separately does not identify independent effects among
  correlated clades.
- Removed reference to an RStudio project file that is not present in the repository.
- Corrected the count of derived-output entries in the data guide from two to three.

## Open items

1. **Confirm the richness denominator.** The README and script 01 describe vascular
   plant richness, while the data guide describes the Cai et al. layer as woody-plant
   richness. Confirm the source dataset's taxonomic scope and align those descriptions.
   No denominator or numerical result was changed during this review.
2. **Verify the external archive.** The DOI listed in the data guide could not be
   accessed by the review tool. This does not establish that it is invalid.
   Confirm that a fresh download provides every raw layer and the exact derived stacks
   listed in `R/config.R`; record file checksums for reproducibility.
3. **Retain the archived stacks.** `ENVIRONMENT.md` documents varying NA boundaries when
   step 4 rebuilds environmental rasters. Exact headline-number reproduction depends on
   using the archived inputs. The suggested disk-backed workaround has not been tested.
4. **Pin the runtime if distributing a reproducible release.** Package versions are
   documented, but there is no dependency lockfile. `R/requirements.R` installs missing
   packages without enforcing those versions.
5. **Complete publication metadata when available.** The citation currently refers to
   a submitted article without a full bibliographic entry or article DOI.

## Validation limits

The large raw layers and derived raster stacks were not available in this review,
and the execution environment did not permit launching R. R syntax, installation,
full pipeline execution, regenerated figures and exact numerical reproduction were
therefore not tested. Statistical model correctness and every manuscript claim are
not certified by these file and summary-table checks.
