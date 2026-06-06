# Thrombin activity classification workflow

This repository now includes `train_thrombin_models.py`, which implements the modeling procedure you described:

1. **Stratified split:** 60% train, 20% validation, 20% test.
2. **Models:** SVM (linear + RBF), Logistic Regression, Random Forest, kNN, and XGBoost (if installed).
3. **Imbalance handling:** `class_weight='balanced'` for sklearn classifiers; `scale_pos_weight` for XGBoost.
4. **Tuning:** `RandomizedSearchCV` with **5-fold CV** and MCC scoring on the **combined train+validation** set.
5. **Evaluation:** Compare CV MCC vs out-of-sample test MCC, plus test Accuracy and F1.

## Run

```bash
python train_thrombin_models.py --input pcp_aac_ctd_dpc.tsv --n-iter 30
```

You can run the same command for each feature set:

```bash
python train_thrombin_models.py --input pcp_aac.tsv --output-json results_pcp_aac.json
python train_thrombin_models.py --input pcp_aac_ctd.tsv --output-json results_pcp_aac_ctd.json
python train_thrombin_models.py --input pcp_aac_ctd_dpc.tsv --output-json results_pcp_aac_ctd_dpc.json
```

## Notes

- The script expects a `label` column with `0/1` classes.
- The script excludes `sequence` from model features.
- Hyperparameter ranges are broad defaults and can be narrowed for speed.