# Trogocytosis–CRC reproducibility package

This repository contains the analysis code, processed derivatives, source-data tables, quality-control records and final figures for the manuscript:

> **Functional screening identifies an SCD dependency while matched-null analyses expose transcriptomic limits of trogocytosis inference in colorectal cancer**

The package corresponds to manuscript analysis version **V8.1** and repository release **v1.0.0** (2 September 2026).

Canonical code repository: <https://github.com/wangjinglin0905-gif/trogocytosis-crc-reproducibility>

## Scientific scope

The analyses address two deliberately separate questions:

1. Which literature-derived trogocytosis-panel genes are fitness dependencies in colorectal-cancer (CRC) organoids and cell lines?
2. How far can bulk and single-cell RNA data support inference about trogocytosis?

The package reproduces the reported SCD dependency, negative/attenuated RNA-proxy findings, matched-null calibration, survival analyses, SCD–VPS72 replication test and an independent full-cell CRC single-cell replication. It does **not** treat an RNA module, a double-positive transcript profile or an isolated epithelial transcript in a leukocyte as a directly observed trogocytosis event. Event-level confirmation requires surface-protein or imaging measurements with lineage or genotype tracing.

## Repository contents

| Path | Contents |
|---|---|
| `analysis/` | Final V8/V8.1 processed results and resampling distributions |
| `baseline/v7/analysis/` | Frozen V7 processed derivatives used by the portable recalculation scripts |
| `config/` | Frozen analysis contract and parameters |
| `figures/final/` | Final Figures 1–7 in PNG, TIFF, PDF and SVG |
| `figures/source_data/` | Plot-level source-data tables for corrected/extended figures |
| `methods/` | Analysis contracts and argument/terminology ledger |
| `provenance/` | External-source manifest, figure-to-source map and release file manifest |
| `qa/` | Independent numerical, layout and structural quality-control records |
| `scripts/analysis/` | Portable primary and independent recalculation scripts |
| `scripts/figures/` | Figure-generation scripts |
| `scripts/upstream/` | Raw-data acquisition/preparation and upstream reanalysis scripts |
| `scripts/validation/` | Public-release audit, path sanitizer and archive builder |

## Data availability and redistribution boundary

No controlled or identifiable participant data are included. The study uses public, de-identified datasets. Large third-party source files are **not redistributed** because their original repositories provide the authoritative copies and may impose their own terms. `provenance/external_sources.tsv` records the accession/DOI, version, retrieval date, expected file name and SHA-256 where available. The repository contains only processed derivatives needed for numerical checking and figure regeneration.

See [DATA_AVAILABILITY.md](DATA_AVAILABILITY.md) for the full source and licensing boundary.

## Quick numerical verification

Python 3.12 and R 4.6 were used for the frozen release. Install the pinned Python packages and the listed R packages from `environment/` before running the checks.

```bash
python -m pip install -r environment/python-requirements.txt
```

Run from the repository root:

```bash
python scripts/analysis/recalculate_v7_key_results.py --root baseline/v7 --outdir qa/recomputed/v7_python
Rscript scripts/analysis/recalculate_v7_key_results.R baseline/v7 qa/recomputed/v7_r
Rscript scripts/analysis/05_independent_core_statistics_check.R . qa/recomputed/v8_core
python scripts/analysis/04_gse178341_matching_balance_qc.py --analysis-dir analysis/matched_null_gse178341 --outdir qa/recomputed/matched_balance
Rscript scripts/analysis/recalculate_v8_1_corrections.R . qa/recomputed/v8_1_corrections
```

Convenience wrappers are also provided:

```powershell
powershell -ExecutionPolicy Bypass -File scripts/reproduce_core.ps1
```

```bash
bash scripts/reproduce_core.sh
```

The independent core-statistics check is expected to report **15/15 PASS**. The frozen V8.1 correction audit should reproduce:

- matched-null one-sided empirical P = 0.492651 and equal-tail two-sided P = 0.985301;
- TCGA age- and categorical-stage-adjusted HR = 1.010619 (95% CI 0.839653–1.216396), P = 0.911053;
- Schoenfeld module P = 0.873385 and global P = 0.701304.

## Figure regeneration

Figures 1–6 from the frozen V7 derivatives:

```bash
Rscript scripts/figures/make_v7_figures.R baseline/v7 qa/recomputed/figures_v7
```

Corrected/extended Figures 3, 4, 6 and 7:

```bash
Rscript scripts/figures/make_v8_1_changed_figures.R . analysis/v8_1_corrections qa/recomputed/figures_v8_1
```

Figure 6 uses a fixed plotting seed (`20260901`) so that jittered point positions are reproducible. The statistical values are not generated from plot jitter.

## Full raw-data reanalysis

Scripts under `scripts/upstream/` and `scripts/analysis/01_*` through `03_*` reconstruct the processed derivatives from the external public files. These runs require the source files listed in `provenance/external_sources.tsv`; users should download them from the original repositories and verify the recorded hashes before analysis. Each script exposes its required paths through command-line arguments (`--help`).

The GSE132465 library-size derivative included in `analysis/gse132465_replication/input_derivatives/` is a deterministic, non-expression summary extracted from the public matrix. Its extraction script and SHA-256 manifest are included.

## Release audit and archive

Before a release, run:

```bash
python scripts/validation/audit_public_release.py --root . --out qa/public_release_audit.json
python scripts/validation/build_release_archive.py --root . --version 1.0.0
```

The audit checks required metadata, panel formats, unexpectedly large files, local absolute paths, common credential patterns, symlinks and excluded raw source files. The archive builder writes a file-level SHA-256 manifest and a deterministic ZIP under `release/`.

## Citation

Use the metadata in [CITATION.cff](CITATION.cff). The immutable archive for release `v1.0.0` is available at the version-specific DOI [10.5281/zenodo.22239612](https://doi.org/10.5281/zenodo.22239612). The concept DOI [10.5281/zenodo.22239611](https://doi.org/10.5281/zenodo.22239611) resolves to the latest archived version. Cite the version-specific DOI when reproducing the analyses reported here; use the concept DOI when referring to the software package across versions.

## Licensing

- Code in `scripts/` is released under the MIT License; see [LICENSE](LICENSE).
- Original documentation, derived tables and figures in this repository are released under CC BY 4.0; see [LICENSE_DATA.md](LICENSE_DATA.md).
- Third-party source datasets are not relicensed or redistributed. Their original terms continue to apply; see [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).

## Contact and funding

Correspondence: Juan Yang, MD, PhD — `yj63yj63@163.com`.

This work was supported by the National Natural Science Foundation of China (No. 82360467) and the Science and Technology Fund Project of the Guizhou Provincial Science and Technology Program (Nos. QKH JC-ZK [2023]-358 and QKH JC-ZK [2023]-348). The funders had no role in the review design, interpretation, manuscript preparation or decision to submit the manuscript.
