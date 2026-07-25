# File and Folder Guide

## Data inputs and processed feature tables

| Location | Purpose |
|---|---|
| `TableS1PositiveProteins.fasta` | Source positive peptide sequences |
| `filtered_negatives.fasta` | Filtered negative sequences |
| `final_negatives.tsv` | Final negative-sequence table |
| `pcp_aac.tsv` | Base PCP and AAC features |
| `pcp_aac_ctd.tsv` | Base PCP, AAC, and CTD features |
| `pcp_aac_ctd_dpc.tsv` | Base PCP, AAC, CTD, and DPC features |
| `s1_proxy_no_overlap_*` | Cleaned proxy-integrated inputs, features, and audit |

## Feature and dataset construction code

| Location | Purpose |
|---|---|
| `extract_features.py` | General peptide feature extraction |
| `extract_training_features.py` | Training-table feature extraction |
| `build_proxy_integrated_no_overlap_dataset.py` | Builds cleaned proxy dataset |

## Modeling code

| Location | Purpose |
|---|---|
| `train_thrombin_models.py` | Main single-dataset modeling workflow |
| `train_and_validate_all_models_by_dataset.py` | Compares models by feature set |
| `train_and_validate_all_models_s3_s5.py` | Trains and validates S3/S5 workflow |
| `external_validate_combined_no_overlap_10seed.py` | Repeated external validation |

## Current validation results

| Location | Status |
|---|---|
| `external_validation_no_overlap_rerun/` | Current cleaned non-proxy results |
| `external_validation_proxy_no_overlap/` | Current cleaned proxy results |
| `outputs/validation_graphs/` | Current comparison figures |
| `external_validation_no_overlap/` | Earlier non-proxy run; retain for history |
| `external_validation_no_overlap_Test_notActualResults/` | Test/non-final; do not report |

## Model interpretation

| Location | Purpose |
|---|---|
| `shap_outputs/` | Original SHAP outputs |
| `shap_outputs_high_quality/` | Higher-quality SHAP outputs |
| `stable_shap_results/` | Repeated-split stability results |
| `shap_pcp_aac*/` | Feature-set-specific SHAP results |

## Other results and figures

| Location | Purpose |
|---|---|
| `summary_outputs/` | Consolidated internal/external summaries |
| `confusion_matrices/` | Confusion matrices and supporting tables |
| `pca_plots/` | Interactive PCA visualizations |
| `results_by_dataset/` | Model results separated by feature set |
| `saved_models/` | Serialized trained models |

## Presentations

| File | Purpose |
|---|---|
| `new_peptide_external_validation_summary.pptx` | New-peptide external validation |
| `proxy_peptide_integration_validation_summary.pptx` | Proxy analysis summary |
| `external_validation_proxy_vs_nonproxy_slide.pptx` | Direct comparison slide |

The adjacent `.inspect.ndjson` files are machine-generated inspection metadata, not
primary scientific results.

## Third-party software

`Code-folder/`, `cd-hit-auxtools/`, `psi-cd-hit/`, `usecases/`, and `doc/` contain
CD-HIT source, executables, utilities, documentation, and examples. These materials
are retained for provenance but should be clearly separated from project-authored
analysis when uploaded.

## Files that can be omitted from a polished shared folder

- `.DS_Store`
- `tmp/`
- `.git/` unless version-control history is intentionally being transferred
- compiled `.o` files when working executables/source are already retained
- `external_validation_no_overlap_Test_notActualResults/` from the primary results
  view (retain it only in an archive/history section)

