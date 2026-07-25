# Machine-Learning-Assisted Antithrombotic Peptide Research

This repository contains the complete computational workflow for comparing
antithrombotic-peptide classifiers trained:

1. without proxy peptides (`baseline`), and
2. with 30 candidate proxy peptides (`proxy`).

Both conditions now live in the same repository and use the same analysis code.
Keeping them side by side makes the comparison reproducible and prevents the two
Git branches from silently drifting apart.

## Research question

Can numerical features derived from peptide sequences help identify
antithrombotic peptides and prioritize candidates for later experimental testing?

This project is a computational screening study. A positive model prediction does
not demonstrate biological activity or clinical effectiveness.

## Quick orientation

```text
configs/        Human-readable settings for baseline and proxy runs
data/           Raw inputs, processed feature tables, and external predictions
docs/           Mentor notes, workflow history, data dictionary, and script guide
handoff/        Dated mentor handoff package
presentations/  Final presentation files and inspection metadata
results/        Baseline, proxy, comparison, model, figure, and SHAP outputs
src/            Project-authored Python and R code grouped by purpose
third_party/    CD-HIT source, executables, documentation, and examples
```

Each major folder contains its own README.

## Which results are current?

Use:

- `results/baseline/external_validation/`
- `results/proxy/external_validation/`
- `results/comparisons/validation_graphs/`

Do not report results from:

- `results/archive/nonproxy_test_not_final/`

That folder is retained only to document experiment history.

## Main results

### Baseline, without proxy peptides

- Best mean internal test MCC: **0.849**, DPC + SVM-RBF.
- Best mean external positive detection: **78.2%**, DPC + kNN.
- External detection was evaluated across 10 training seeds on 11 known-positive,
  non-overlapping external peptides.

### Proxy-integrated training

- Best mean internal test MCC shown in the project summary: **0.740**,
  CTD + SVM-RBF.
- DPC + kNN mean external positive detection: **78.6%** in the proxy-validation
  folder.

External positive detection is sensitivity/recall on known positives, not overall
accuracy. The external sets are small, so the results are preliminary.

## Reproduce the main workflows

Run commands from the repository root.

### 1. Install Python dependencies

```bash
python -m pip install -r requirements.txt
```

R figure scripts use base R and do not currently require an R package lockfile.

### 2. Rebuild baseline features

```bash
python src/features/extract_training_features.py
```

### 3. Rebuild proxy-integrated features

```bash
python src/data/build_proxy_integrated_no_overlap_dataset.py
```

### 4. Run current baseline external validation

```bash
python src/validation/external_validate_combined_no_overlap_10seed.py \
  --tables \
  data/processed/baseline/pcp_aac.tsv \
  data/processed/baseline/pcp_aac_ctd.tsv \
  data/processed/baseline/pcp_aac_ctd_dpc.tsv \
  --outdir results/baseline/external_validation
```

### 5. Run current proxy external validation

```bash
python src/validation/external_validate_combined_no_overlap_10seed.py \
  --tables \
  data/processed/proxy/s1_proxy_no_overlap_pcp_aac.tsv \
  data/processed/proxy/s1_proxy_no_overlap_pcp_aac_ctd.tsv \
  data/processed/proxy/s1_proxy_no_overlap_pcp_aac_ctd_dpc.tsv \
  --outdir results/proxy/external_validation
```

### 6. Rebuild comparison figures

```bash
Rscript src/visualization/make_internal_external_validation_graphs.R
Rscript src/visualization/make_scientific_seed_validation_figure.R
Rscript src/visualization/make_direct_proxy_validation_comparison.R
Rscript src/visualization/make_three_model_proxy_comparison.R
Rscript src/visualization/make_shap_best_model_table.R
```

## Important limitations

- The training dataset is imbalanced.
- External-validation sets contain only known positives, so specificity and
  full external accuracy cannot be calculated.
- There are only 11 non-overlapping positives in the current baseline external set.
- Proxy integration changes the composition of the positive class.
- Model outputs must be experimentally validated.
- Deep-learning numbers reported in the summer presentation are not reproducible
  from this repository because the training code, pair file, logs, and checkpoints
  are absent.

## Documentation for a mentor or new researcher

Read these files in order:

1. `README.md`
2. `docs/EXPERIMENT_NOTES.md`
3. `docs/DATA_DICTIONARY.md`
4. `docs/SCRIPT_REFERENCE.md`
5. `docs/BRANCH_HISTORY.md`
6. `docs/REORGANIZATION_MAP.md`
7. `docs/FILE_GUIDE.md`
