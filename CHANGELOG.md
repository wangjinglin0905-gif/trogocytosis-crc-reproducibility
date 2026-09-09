# Changelog

## v1.1.0 — 2026-09-09

- Added the approved submission-display Figures 1–7 and their portable plotting sources and hash/source maps; preserved the v1.0.0 exports as historical assets.
- Figure 1 uses full display ranges, Figure 2 separates P-value labels from intervals, Figure 4 retains the categorical-stage model and top-positioned continuous-Cox labels, and Figures 3/6/7 retain the approved panel/label layout.
- No frozen numerical analysis, threshold, cohort definition, gene set or conclusion was changed. Statistical reruns were not performed for this display/documentation release.
- Clarified the separation between manuscript V8.6, frozen numerical basis V8.1 and software release v1.1.0; corrected the README funding-role wording from review to study.
- Private manuscript files and screening reports are not distributed. MSigDB GMT and other third-party raw data remain excluded.
- New Zenodo version pending; the old version-specific DOI continues to identify v1.0.0 only.

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
