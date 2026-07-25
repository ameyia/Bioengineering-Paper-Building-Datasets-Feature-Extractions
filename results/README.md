# Results

- `baseline/`: results produced without proxy peptides.
- `proxy/`: results produced with proxy-integrated training.
- `comparisons/`: direct, shared-set comparisons between the two conditions.
- `archive/`: older or explicitly non-final runs retained for provenance.
- `figures/`: confusion matrices and PCA outputs.
- `interpretability/`: SHAP outputs and repeated-split feature stability.
- `models/`: serialized estimators.

Use the current external-validation folders under `baseline/` and `proxy/` for
reporting. Never treat `archive/nonproxy_test_not_final/` as a final result.

