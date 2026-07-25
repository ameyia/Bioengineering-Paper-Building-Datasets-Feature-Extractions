# Script Reference

All commands should be run from the repository root.

## Dataset construction

| Script | Purpose | Main output |
|---|---|---|
| `src/data/build_proxy_integrated_no_overlap_dataset.py` | Combines original positives, proxy candidates, and negatives; removes duplicates/overlaps; extracts features | `data/processed/proxy/` |

## Feature extraction

| Script | Purpose | Notes |
|---|---|---|
| `src/features/extract_training_features.py` | Rebuilds all three baseline feature tables | Uses raw FASTA/TSV inputs |
| `src/features/extract_features.py` | General tunable model/feature utility from the earlier workflow | Retained for experiment history |

## Modeling

| Script | Purpose |
|---|---|
| `src/modeling/train_thrombin_models.py` | Tunes and evaluates one supplied feature table |
| `src/modeling/train_and_validate_all_models_by_dataset.py` | Runs models across feature sets and predicts S3 peptides |
| `src/modeling/train_and_validate_all_models_s3_s5.py` | Earlier combined S3/S5 workflow |
| `src/modeling/summarize_all_results.py` | Consolidates model outputs |
| `src/modeling/summarize_all_results_s3_s5.py` | Consolidates S3/S5 outputs |

## Validation

| Script | Purpose |
|---|---|
| `src/validation/external_validate_combined_no_overlap_10seed.py` | Current repeated internal/external validation; use for both baseline and proxy tables |
| `src/validation/external_validate_s3.py` | Earlier single-model S3 validation |

## Interpretation

| Script | Purpose |
|---|---|
| `src/interpretation/run_shap_all_models.py` | Original SHAP workflow |
| `src/interpretation/run_shap_all_models_high_quality.py` | Higher-quality SHAP export attempt |
| `src/interpretation/run_shap_all_models_high_quality_fixed.py` | Corrected high-quality SHAP workflow; prefer this over the preceding version |
| `src/interpretation/stable_shap_repeated_splits.py` | Measures feature stability across repeated splits |
| `src/interpretation/summarize_shap_outputs.py` | Consolidates SHAP outputs |
| `src/interpretation/shap_pcp_aac*.py` | Focused scripts for selected feature/model combinations |

## Visualization

| Script | Purpose |
|---|---|
| `src/visualization/make_all_confusion_matrices.py` | Produces confusion matrices across datasets, models, and seeds |
| `src/visualization/make_confusion_matrix.py` | Produces one selected confusion matrix |
| `src/visualization/build_pca_comparison_plots.py` | Creates feature-set PCA comparisons |
| `src/visualization/build_3d_pca_plot.py` | Creates interactive 3D PCA output |
| `src/visualization/make_internal_external_validation_graphs.R` | Main internal/external comparison figure set |
| `src/visualization/make_scientific_seed_validation_figure.R` | Seed-level scientific validation figure |
| `src/visualization/make_direct_proxy_validation_comparison.R` | Direct shared-set proxy comparison |
| `src/visualization/make_three_model_proxy_comparison.R` | Selected-model comparison |
| `src/visualization/make_shap_best_model_table.R` | Presentation/manuscript table for SHAP model choices |

