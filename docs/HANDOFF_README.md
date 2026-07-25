# Antithrombotic Peptide Machine-Learning Project

## Start here

This repository contains the summer research workflow for building peptide datasets,
extracting sequence features, training machine-learning classifiers, evaluating them
on held-out and external peptides, testing proxy-peptide integration, and generating
interpretability and validation figures.

The main scientific question was:

> Can peptide-sequence features and machine-learning models help prioritize
> antithrombotic peptides for eventual use in vascular-material coatings?

The work is computational and preliminary. It prioritizes candidates for later
experimental testing; it does not establish that a peptide is clinically effective.

## Main conclusions

1. Internal test performance was strong for several high-dimensional feature/model
   combinations. In the cleaned non-proxy analysis, DPC + SVM-RBF had a mean test
   MCC of 0.849 across 10 seeds.
2. The model with the best internal MCC was not the best at detecting the separate
   external positives. DPC + kNN had the highest mean external positive-detection
   rate, 78.2%, in the non-proxy rerun.
3. Adding proxy peptides broadened the positive class and generally reduced the
   strongest internal MCC values. Simpler AAC models improved slightly in some
   comparisons, while specialized DPC models declined more.
4. External results should be treated as preliminary because the external sets are
   small and contain positive examples only. “External positive detection” is
   sensitivity/recall on known positives, not overall accuracy.

## Recommended results to use

Use these folders for the current presentation and handoff:

- `results/baseline/external_validation/`: cleaned non-proxy, 10-seed results.
- `results/proxy/external_validation/`: cleaned proxy-integrated results.
- `results/comparisons/validation_graphs/`: comparison figures generated from the validation
  summaries.
- `s1_proxy_no_overlap_*`: proxy-integrated inputs and overlap audit.

Do not use `results/archive/nonproxy_test_not_final/` for scientific
reporting. Its name marks it as a test/non-final output.

## Important dataset counts

- `data/processed/baseline/pcp_aac.tsv`, the older/base feature file, contains 88 positive and 435 negative
  rows after the header.
- `s1_proxy_no_overlap_pcp_aac.tsv`, the cleaned proxy-integrated feature file,
  contains 91 positive and 421 negative rows.
- The current non-proxy rerun evaluates 11 non-overlapping external positive
  peptides.
- The proxy external-validation folder contains 14 external positives. When
  comparing proxy and non-proxy performance directly, use the shared external set
  specified by the comparison scripts rather than comparing percentages with
  different denominators.

## Reproducibility warning

The presentation mentions a sequence-to-sequence deep-learning experiment with
2,420 training pairs, a 9.4x batching speedup, and 31.2% held-out reconstruction
accuracy. The corresponding training script, input-pair file, checkpoint, and
metric log are not present in this repository. Those values cannot be independently
reproduced from the current files and should be added separately if available.

## Documentation

- `EXPERIMENT_NOTES.md`: what was tried, in chronological/logical order.
- `FILE_GUIDE.md`: which files and folders contain inputs, code, and results.
- `Modeling_workflow.md`: original modeling notes and example commands.

## Suggested handoff reading order

1. Read this file.
2. Read `EXPERIMENT_NOTES.md`.
3. Review `results/comparisons/validation_graphs/`.
4. Inspect the two current external-validation folders.
5. Use `FILE_GUIDE.md` to locate scripts or detailed result tables.

