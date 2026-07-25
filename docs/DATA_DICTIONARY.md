# Data Dictionary

## Shared identifiers and outcomes

| Field | Meaning |
|---|---|
| `sequence` | Amino-acid sequence using one-letter codes |
| `label` | `1` for positive and `0` for negative |
| `record_id` | Identifier assigned during dataset construction |
| `source` | Origin of the peptide record |
| `status` | Whether an audited record was kept or removed |
| `reason` | Reason for removal, usually overlap or duplication |

## Physicochemical and composition features

| Feature group | Meaning |
|---|---|
| PCP | Length, molecular weight, aromaticity, pI, and instability |
| AAC | Percentage composition for the 20 amino acids |
| CTD | Composition, transition, and distribution descriptors |
| DPC | Frequencies of the 400 possible adjacent amino-acid pairs |

The feature tables contain:

| Table | Predictors | Total columns |
|---|---:|---:|
| PCP + AAC | 25 | 27 |
| PCP + AAC + CTD | 172 | 174 |
| PCP + AAC + CTD + DPC | 572 | 574 |

Total columns include `sequence` and `label`.

## Validation fields

| Field | Meaning |
|---|---|
| `seed` | Random seed used for a repeated split/training run |
| `test_mcc` | Matthews correlation coefficient on the internal test split |
| `test_accuracy` | Fraction correct on the internal test split |
| `test_f1` | Harmonic mean of precision and recall |
| `test_precision` | Positive predictive value |
| `test_recall` | Sensitivity on the internal test split |
| `predicted_label` | Model class prediction |
| `positive_probability` | Estimated probability of the positive class |
| `external_positive_fraction` | Fraction of known external positives detected |

## Current row counts

- Baseline feature file: 88 positives and 435 negatives.
- Cleaned proxy feature file: 91 positives and 421 negatives.
- Current baseline external set: 11 known positives after overlap filtering.

Counts reflect the files currently stored in this repository and may differ from
older, uncleaned counts shown in early presentation drafts.

