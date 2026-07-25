# Experiment Notes

## 1. Dataset assembly

Positive antithrombotic/thrombin-related peptide sequences and negative peptide
sequences were collected and cleaned. CD-HIT utilities and overlap checks were used
to reduce sequence redundancy and prevent obvious overlap between training and
external-validation peptides.

Relevant files:

- `TableS1PositiveProteins.fasta`
- `filtered_negatives.fasta`
- `final_negatives.tsv`
- `Code-folder/`
- `results/baseline/external_validation/combined_external_validation_manifest.tsv`

The manifest is important because it records which proposed external peptides were
kept or removed and why.

## 2. Feature extraction

Peptide sequences were converted into numerical feature tables:

- PCP + AAC: 25 predictor features.
- PCP + AAC + CTD: 172 predictor features.
- PCP + AAC + CTD + DPC: 572 predictor features.

The two non-feature columns are `sequence` and `label`, which explains why the TSV
files contain 27, 174, and 574 total columns.

Relevant scripts and outputs:

- `extract_features.py`
- `extract_training_features.py`
- `data/processed/baseline/pcp_aac.tsv`
- `data/processed/baseline/pcp_aac_ctd.tsv`
- `data/processed/baseline/pcp_aac_ctd_dpc.tsv`

## 3. Baseline model comparison

The following classifiers were compared:

- linear SVM
- RBF SVM
- logistic regression
- random forest
- k-nearest neighbors
- XGBoost

The work used stratified splitting, imbalance handling, hyperparameter tuning, and
metrics including MCC, accuracy, F1, precision, and recall. MCC was emphasized
because the dataset is imbalanced.

Relevant scripts and outputs:

- `train_thrombin_models.py`
- `train_and_validate_all_models_by_dataset.py`
- `results/baseline/results_by_dataset/`
- `results/baseline/model_summaries/`
- `model_comparison_results.json`

## 4. Repeated internal and external validation

Models were rerun across 10 random seeds. For each seed, the trained model was
evaluated internally and then applied to known-positive external peptides excluded
from training.

Current non-proxy results:

- `results/baseline/external_validation/internal_metrics_by_seed.tsv`
- `results/baseline/external_validation/internal_metrics_summary_ci.tsv`
- `results/baseline/external_validation/combined_external_predictions_by_seed.tsv`
- `results/baseline/external_validation/combined_external_summary_ci.tsv`

Example interpretation:

- DPC + SVM-RBF mean internal test MCC: 0.849.
- DPC + kNN mean external positive detection: 0.782, displayed as 78%.

The 78% is the mean across 10 seeds on 11 known positives. Across all runs, DPC +
kNN produced 86 positive calls out of 110 peptide/run opportunities.

Limitation: because all external examples are positive, this analysis measures
positive detection/sensitivity only. It does not measure external specificity,
false-positive rate, or complete external accuracy.

## 5. Proxy-peptide integration

Thirty proxy candidates were considered to broaden the positive class. The build
script performed overlap filtering and created cleaned positive, negative, feature,
and audit files.

Relevant files:

- `build_proxy_integrated_no_overlap_dataset.py`
- `s1_proxy_no_overlap_input_manifest.tsv`
- `s1_proxy_no_overlap_overlap_audit.tsv`
- `s1_proxy_no_overlap_positives.fasta`
- `s1_proxy_no_overlap_negatives.tsv`
- `s1_proxy_no_overlap_pcp_aac.tsv`
- `s1_proxy_no_overlap_pcp_aac_ctd.tsv`
- `s1_proxy_no_overlap_pcp_aac_ctd_dpc.tsv`
- `results/proxy/external_validation/`

Observed result:

- Best proxy-integrated internal result shown in the presentation: CTD + SVM-RBF,
  mean test MCC 0.740.
- DPC + kNN had a mean external positive-detection rate of 0.786 in the proxy
  external folder.

Interpretation:

The proxy set introduced greater sequence/mechanistic diversity. A possible
explanation is that rich DPC models had learned a narrower positive class, while
simpler AAC representations accommodated the added diversity somewhat better.
This is an interpretation, not a mechanism proven by the current experiments.

## 6. Interpretability

SHAP analyses were run for multiple models, datasets, and repeated splits to examine
which features influenced predictions and which features were stable.

Relevant scripts and folders:

- `run_shap_all_models.py`
- `run_shap_all_models_high_quality.py`
- `run_shap_all_models_high_quality_fixed.py`
- `stable_shap_repeated_splits.py`
- `results/interpretability/shap_outputs/`
- `results/interpretability/shap_outputs_high_quality/`
- `results/interpretability/stable_shap_results/`
- `results/interpretability/shap_pcp_aac/`
- `results/interpretability/shap_pcp_aac_ctd/`
- `results/interpretability/shap_pcp_aac_ctd_dpc/`

The `*_fixed.py` name suggests it supersedes the earlier high-quality script, but
both are retained to preserve experiment history.

## 7. Figures and presentations

Validation figures are stored in `results/comparisons/validation_graphs/`. PNG files are useful
for slides, while PDF files preserve vector quality for manuscripts.

Presentations currently stored in the repository summarize:

- new-peptide external validation;
- proxy-peptide integration;
- direct proxy versus non-proxy validation.

Presentation-building code was intentionally removed from the working repository at
the user's request. Existing presentation files and scientific figure-generation
scripts were retained.

## 8. Deep-learning feasibility work

Slides report an early sequence-to-sequence reconstruction experiment:

- 2,420 training pairs;
- 9.4x faster after batching;
- 31.2% held-out exact reconstruction accuracy at 300 epochs;
- decreasing reconstruction loss without a clear plateau.

The underlying training artifacts are absent from this repository. To make this
work reproducible, add the model script, pair-generation method, train/validation
split, random seed, dependency versions, metric log, and checkpoint.

## 9. Recommended next steps

1. Define a principled proxy weighting or positive-unlabeled strategy.
2. Compare proxy and non-proxy models on exactly the same external set.
3. Expand external validation to include both positive and negative peptides.
4. Report confidence intervals and seed-level distributions alongside means.
5. Save complete deep-learning inputs, code, checkpoints, and logs.
6. Select candidates for experimental testing only after checking applicability
   domain, sequence similarity, and uncertainty.

