# Git Branch History

## Current branches

- `main`: points to the original non-proxy baseline commit.
- `bioengineering-no-proxy-baseline`: historical baseline branch at the same commit
  as `main`.
- `proxy-peptide-integration`: contains the proxy integration and the unified
  organization documented here.

## Recommended future use

After this reorganization is reviewed:

1. Preserve `bioengineering-no-proxy-baseline` as a read-only historical reference.
2. Merge the organized proxy branch into `main`.
3. Use `main` for both conditions, controlled by the paths in `configs/`.
4. Avoid maintaining separate copies of shared modeling code on two branches.
5. Create short-lived feature branches only for new analyses.

This design makes the baseline/proxy comparison visible in one checkout and reduces
the risk that bug fixes are applied to only one condition.

