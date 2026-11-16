# Figures and result tables

These are committed outputs from the analysis archive. Start with the
[main figures](../README.md#main-figures), then use this index for diagnostics,
robustness checks and the prickly/thorny variants.

- **Descriptive:** global richness, biome distributions and ecoregion summaries.
- **Model 1:** environmental ALE responses, collinearity, spatial validation and species-drop sensitivity.
- **Model 2:** clade associations, environmental adjustment, R² partitions and sensitivity analyses.
- **Ecoregion mammals:** species lists and order counts for the selected ecoregions.

Rasters and the full predictor table are external/generated inputs; see
[data requirements](../data/README.md) and [environment notes](../ENVIRONMENT.md).
PNG figures are previews; use the accompanying CSVs for numerical values.

## descriptive

### Figures

- [FIG_biome_prickly.png](descriptive/FIG_biome_prickly.png)
- [FIG_biome_thorny.png](descriptive/FIG_biome_thorny.png)
- [FIG_biome_total.png](descriptive/FIG_biome_total.png)
- [FIG_composite_prickly.png](descriptive/FIG_composite_prickly.png)
- [FIG_composite_thorny.png](descriptive/FIG_composite_thorny.png)
- [FIG_composite_total.png](descriptive/FIG_composite_total.png)
- [FIG_ecoregion_prickly.png](descriptive/FIG_ecoregion_prickly.png)
- [FIG_ecoregion_thorny.png](descriptive/FIG_ecoregion_thorny.png)
- [FIG_ecoregion_total.png](descriptive/FIG_ecoregion_total.png)
- [MAP_absolute_prickly.png](descriptive/MAP_absolute_prickly.png)
- [MAP_absolute_richness.png](descriptive/MAP_absolute_richness.png)
- [MAP_absolute_thorny.png](descriptive/MAP_absolute_thorny.png)
- [MAP_relative_prickly.png](descriptive/MAP_relative_prickly.png)
- [MAP_relative_richness.png](descriptive/MAP_relative_richness.png)
- [MAP_relative_thorny.png](descriptive/MAP_relative_thorny.png)

### Tables

- [biome_summary_prickly.csv](descriptive/biome_summary_prickly.csv)
- [biome_summary_thorny.csv](descriptive/biome_summary_thorny.csv)
- [biome_summary_total.csv](descriptive/biome_summary_total.csv)
- [ecoregion_summary_prickly.csv](descriptive/ecoregion_summary_prickly.csv)
- [ecoregion_summary_thorny.csv](descriptive/ecoregion_summary_thorny.csv)
- [ecoregion_summary_total.csv](descriptive/ecoregion_summary_total.csv)
- [within_ecoregion_correlogram.csv](descriptive/within_ecoregion_correlogram.csv)
- [within_ecoregion_correlogram_summary.csv](descriptive/within_ecoregion_correlogram_summary.csv)

## model1

### Figures

- [M1_ALE.png](model1/M1_ALE.png)
- [M1_collinearity.png](model1/M1_collinearity.png)
- [M1_diagnostics.png](model1/M1_diagnostics.png)
- [M1_jackknife_ALE.png](model1/M1_jackknife_ALE.png)

### Tables

- [ale_curves.csv](model1/ale_curves.csv)
- [blocksize_confirmatory.csv](model1/blocksize_confirmatory.csv)
- [blocksize_recheck.csv](model1/blocksize_recheck.csv)
- [calibration.csv](model1/calibration.csv)
- [input_summary.csv](model1/input_summary.csv)
- [jackknife.csv](model1/jackknife.csv)
- [model1_oof.csv](model1/model1_oof.csv)
- [performance_summary.csv](model1/performance_summary.csv)
- [predictor_concurvity.csv](model1/predictor_concurvity.csv)
- [predictor_correlation.csv](model1/predictor_correlation.csv)
- [ratio_check_ale.csv](model1/ratio_check_ale.csv)
- [ratio_check_performance.csv](model1/ratio_check_performance.csv)
- [ratio_check_shapes.csv](model1/ratio_check_shapes.csv)
- [residual_correlogram.csv](model1/residual_correlogram.csv)
- [resolution_sweep.csv](model1/resolution_sweep.csv)
- [species_drop_ale_envelope.csv](model1/species_drop_ale_envelope.csv)
- [species_drop_ale_stability.csv](model1/species_drop_ale_stability.csv)
- [transform_test.csv](model1/transform_test.csv)
- [tuning_grid.csv](model1/tuning_grid.csv)

## model2

### Figures

- [FIG_M2_composite.png](model2/FIG_M2_composite.png)
- [FIG_M2_composite_prickly.png](model2/FIG_M2_composite_prickly.png)
- [FIG_M2_composite_thorny.png](model2/FIG_M2_composite_thorny.png)
- [FIG_M2_confounding_prickly.png](model2/FIG_M2_confounding_prickly.png)
- [FIG_M2_confounding_thorny.png](model2/FIG_M2_confounding_thorny.png)
- [FIG_M2_confounding_total.png](model2/FIG_M2_confounding_total.png)
- [FIG_M2_jackknife.png](model2/FIG_M2_jackknife.png)
- [FIG_M2_r2_prickly.png](model2/FIG_M2_r2_prickly.png)
- [FIG_M2_r2_thorny.png](model2/FIG_M2_r2_thorny.png)
- [FIG_M2_r2_total.png](model2/FIG_M2_r2_total.png)
- [FIG_M2_scoring_alone_prickly.png](model2/FIG_M2_scoring_alone_prickly.png)
- [FIG_M2_scoring_alone_thorny.png](model2/FIG_M2_scoring_alone_thorny.png)
- [FIG_M2_scoring_alone_total.png](model2/FIG_M2_scoring_alone_total.png)
- [FIG_M2_scoring_env_prickly.png](model2/FIG_M2_scoring_env_prickly.png)
- [FIG_M2_scoring_env_thorny.png](model2/FIG_M2_scoring_env_thorny.png)
- [FIG_M2_scoring_env_total.png](model2/FIG_M2_scoring_env_total.png)
- [FIG_M2_shape_alone_mono_prickly.png](model2/FIG_M2_shape_alone_mono_prickly.png)
- [FIG_M2_shape_alone_mono_thorny.png](model2/FIG_M2_shape_alone_mono_thorny.png)
- [FIG_M2_shape_alone_mono_total.png](model2/FIG_M2_shape_alone_mono_total.png)
- [FIG_M2_shape_alone_unc_prickly.png](model2/FIG_M2_shape_alone_unc_prickly.png)
- [FIG_M2_shape_alone_unc_thorny.png](model2/FIG_M2_shape_alone_unc_thorny.png)
- [FIG_M2_shape_alone_unc_total.png](model2/FIG_M2_shape_alone_unc_total.png)
- [FIG_M2_shape_env_mono_prickly.png](model2/FIG_M2_shape_env_mono_prickly.png)
- [FIG_M2_shape_env_mono_thorny.png](model2/FIG_M2_shape_env_mono_thorny.png)
- [FIG_M2_shape_env_mono_total.png](model2/FIG_M2_shape_env_mono_total.png)
- [FIG_M2_shape_env_unc_prickly.png](model2/FIG_M2_shape_env_unc_prickly.png)
- [FIG_M2_shape_env_unc_thorny.png](model2/FIG_M2_shape_env_unc_thorny.png)
- [FIG_M2_shape_env_unc_total.png](model2/FIG_M2_shape_env_unc_total.png)

### Tables

- [TAB_M2_confounding.csv](model2/TAB_M2_confounding.csv)
- [TAB_M2_coverage.csv](model2/TAB_M2_coverage.csv)
- [TAB_M2_jackknife.csv](model2/TAB_M2_jackknife.csv)
- [TAB_M2_r2.csv](model2/TAB_M2_r2.csv)
- [TAB_M2_scoring.csv](model2/TAB_M2_scoring.csv)
- [TAB_M2_shapes.csv](model2/TAB_M2_shapes.csv)
- [TAB_S_model2_autocorrelation.csv](model2/TAB_S_model2_autocorrelation.csv)
- [TAB_S_model2_settings.csv](model2/TAB_S_model2_settings.csv)
- [TAB_S_model2_strata.csv](model2/TAB_S_model2_strata.csv)
- [TAB_S_model2_summary.csv](model2/TAB_S_model2_summary.csv)
- [model2_curves.csv](model2/model2_curves.csv)
- [model2_curves_constrained.csv](model2/model2_curves_constrained.csv)
- [model2_curves_env.csv](model2/model2_curves_env.csv)
- [model2_curves_env_constrained.csv](model2/model2_curves_env_constrained.csv)
- [model2_env_results.csv](model2/model2_env_results.csv)
- [model2_jackknife_ranks.csv](model2/model2_jackknife_ranks.csv)
- [model2_jackknife_summary.csv](model2/model2_jackknife_summary.csv)
- [model2_results.csv](model2/model2_results.csv)

## ecoregion_mammal

### Tables

- [ecoregion_mammals.csv](ecoregion_mammal/ecoregion_mammals.csv)
- [ecoregion_order_counts.csv](ecoregion_mammal/ecoregion_order_counts.csv)

