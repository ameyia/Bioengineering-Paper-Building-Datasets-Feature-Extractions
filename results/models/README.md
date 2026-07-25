# Saved models

The three feature-set folders contain serialized estimators for each classifier:

- `pcp_aac/`
- `pcp_aac_ctd/`
- `pcp_aac_ctd_dpc/`

`best_model.pkl` and `xgboost_pcp_aac_ctd.pkl` are earlier standalone exports.
Serialized models can depend on the Python and library versions used during
training. When exact reproducibility matters, retrain from the feature tables and
save environment versions alongside the new model.

