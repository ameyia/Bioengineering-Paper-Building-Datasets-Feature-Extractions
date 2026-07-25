# Reorganization Map

This file helps readers translate paths found in older commits, notebooks,
presentations, or messages into the current layout.

| Previous location | Current location |
|---|---|
| Root-level `pcp_aac*.tsv` | `data/processed/baseline/` |
| Root-level `s1_proxy_no_overlap_*` | `data/processed/proxy/` |
| Root FASTA and negative TSV inputs | `data/raw/` |
| `external_predictions/`, `s3_predictions/` | `data/external/` |
| Root Python/R scripts | `src/` subfolders by purpose |
| `external_validation_no_overlap_rerun/` | `results/baseline/external_validation/` |
| `external_validation_proxy_no_overlap/` | `results/proxy/external_validation/` |
| `external_validation_no_overlap/` | `results/archive/nonproxy_earlier_run/` |
| `external_validation_no_overlap_Test_notActualResults/` | `results/archive/nonproxy_test_not_final/` |
| `outputs/validation_graphs/` | `results/comparisons/validation_graphs/` |
| `summary_outputs/` | `results/baseline/model_summaries/` |
| `results_by_dataset/` | `results/baseline/results_by_dataset/` |
| `saved_models/` and root `.pkl` files | `results/models/` |
| Root SHAP output folders | `results/interpretability/` |
| `confusion_matrices/`, `pca_plots/` | `results/figures/` |
| Root `.pptx` files | `presentations/` |
| `Code-folder/` and CD-HIT utilities | `third_party/` |
| `Modeling_workflow.md` | `docs/Modeling_workflow.md` |
| Temporary PDF renders/compiler caches | Removed; reproducible temporary files |

Git may display these as delete/add pairs until the reorganization is staged.
During staging or commit, Git's rename detection should recognize files whose
contents did not change.

