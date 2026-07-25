# Proxy-integrated dataset

These files were created by
`src/data/build_proxy_integrated_no_overlap_dataset.py`.

- `s1_proxy_no_overlap_input_manifest.tsv`: every proposed input and its source.
- `s1_proxy_no_overlap_overlap_audit.tsv`: records kept or removed and why.
- `s1_proxy_no_overlap_positives.fasta`: cleaned positive sequences.
- `s1_proxy_no_overlap_negatives.tsv`: cleaned negatives.
- `s1_proxy_no_overlap_pcp_aac*.tsv`: feature tables used for proxy training.

The current cleaned table contains 91 positive and 421 negative rows. Thirty proxy
records were retained, while overlapping original records were also removed; this
is why the positive count is not simply the baseline count plus 30.

