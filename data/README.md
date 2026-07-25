# Data

## `raw/`

Source positive sequences, negative sequences, and the proxy-peptide table used to
construct processed datasets. These files should not be overwritten by analysis.

## `processed/baseline/`

Feature tables for training without proxy peptides.

## `processed/proxy/`

The cleaned proxy-integrated positive and negative records, feature tables, input
manifest, and overlap audit.

## `external/`

Earlier S3 prediction tables and per-model external prediction outputs. The current
10-seed validation tables are under `results/`.

See `docs/DATA_DICTIONARY.md` for row counts, feature counts, and field meanings.

