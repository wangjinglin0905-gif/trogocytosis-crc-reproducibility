# Changelog

## v1.0.0 — 2026-09-02

- Packaged the frozen V8.1 analysis, processed derivatives and final Figures 1–7.
- Replaced local absolute paths with command-line inputs and portable relative outputs.
- Added an independent 15-item core-statistics check and public-release audit.
- Corrected the matched-null equal-tail two-sided P-value label without changing the directional test or conclusion.
- Refit the TCGA model with categorical stage and recorded Schoenfeld proportional-hazards tests.
- Added matched-null balance QC, adjusted SCD–VPS72 replication and independent GSE132465 composition replication.
- Added a BFGS optimization fallback with a gradient check for compatibility with current NumPy/SciPy releases; reproduced the frozen results.
- Fixed the Figure 6 plotting seed at 20260901 so that jitter locations are deterministic; statistical results are unchanged.
- Excluded third-party raw data and added accession-, version- and hash-level provenance.

