# Model interpretation results

This folder contains several generations of SHAP analysis:

- `shap_outputs/`: original outputs.
- `shap_outputs_high_quality/`: higher-resolution/repeated outputs.
- `stable_shap_results/`: feature stability across repeated splits.
- `shap_pcp_aac*`: feature-set-specific summaries.
- `shap_single_outputs/`: selected single model/feature analyses.

Prefer the corrected script
`src/interpretation/run_shap_all_models_high_quality_fixed.py` when regenerating
high-quality SHAP results.

