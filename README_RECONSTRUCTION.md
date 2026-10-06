# Reconstruction release v1.2.0 (manuscript V09.3)

This release supports **Cellular composition and the interpretation of candidate RNA readouts of trogocytosis in colorectal cancer**. It is a reproducibility archive, not a journal publication or a preprint of the manuscript.

Version DOI: https://doi.org/10.5281/zenodo.23183934. Series DOI: https://doi.org/10.5281/zenodo.22239611. Earlier tags/releases are retained without replacement.

## Current entry points

- `supplement/v09_3/`: Table S1 evidence annotation and S2–S10 complete supplementary data.
- `methods/Supplementary_methods_v09_3.md`: analysis definitions and correction families.
- `analysis/reconstruction_v09/`: full corrected four-anchor associations, lineage comparisons, two-pair DepMap replication checks, bootstrap/permutation outputs and figure sources.
- `figures/reconstruction_v09_3/`: the author-approved four main and four supplementary figures, each in PNG, TIFF and PDF. These supersede display numbering, not historical files.
- `scripts/reconstruction/`: public path adapters for the actual analysis and plotting algorithms. Run 03 then 04 to obtain the approved Figure 1 layout; 03 preserves five earlier approved panels. Regenerated figures go to `qa/recomputed/`, never overwrite approved exports.
- `provenance/reconstruction_v09_inputs.json` and `qa/reconstruction_v09_gsea_run.json`: input hashes and actual prior numerical run settings. File names are relative resource identifiers, not redistributed raw inputs.

## What changed

Compared with v1.1.0, this release adds corrected organoid partial-correlation inference (df=80; self-gene excluded), separate single-anchor and four-anchor BH families, independent SCD–VPS72 and ATF3–SMYD2 checks, four-target lineage context and a separately labelled 2020-GMT Hallmark exploration. Table S1 now annotates all 39 entries while distinguishing experimental evidence from context; 34 occur in the organoid matrix and the RNA module remains eight genes. FADS2 is outside the 39-gene panel. The manuscript V09.2-to-V09.3 increment only revises source annotation and archive links; no numerical analysis or figure was changed in that increment.

## Reproduction and limits

Use the existing `environment/` specifications and `scripts/analysis/` baseline entry points for the frozen earlier analyses. The reconstruction analysis used Python 3.12.14, NumPy 2.4.6 and the SciPy/Pandas versions recorded in the input manifest. Set `TROGO_RAW_INPUTS` to a local read-only input directory containing `organoid/supplementary_table_6_revision.csv`, `depmap_26Q1/Model.csv` and `depmap_26Q1/CRISPRGeneEffect.csv`, with the exact recorded hashes. Then run `python scripts/reconstruction/01_reanalyse.py`. An optional `TROGO_OUTPUT` selects a fresh output directory. Obtain the source inputs from their original providers; see the existing provenance inventory and data-availability document.

For Hallmark exploration, set `TROGO_HALLMARK_GMT` to the legally obtained 2020-labelled GMT with SHA256 `4275592957a1587652092bb398cf77216fde5b8daa2aedaa0e016f7d10bbdb81`. Run `python scripts/reconstruction/02_hallmark.py`; this uses the supplied full ranking by default, with optional `TROGO_RANKING` and `TROGO_GSEA_OUTPUT`. MSigDB terms apply and GMT files are not redistributed. This exact historical GMT is not asserted equivalent to a current official collection. If it cannot be obtained, the exact Hallmark calculation cannot be reproduced from the archive alone; the complete 200-row output remains available.

No raw single-cell or whole-cohort model rerun was performed for this release. The public path adapters were syntax-checked and their computational operations compared with the actual local-run scripts; their full execution was not repeated. Input hashes and complete outputs preserve traceability but do not establish a new independent replication of the release. RNA measurements are not direct observations of trogocytosis, and absence of corrected significance is not evidence of biological absence.

## Rights and privacy

MIT applies to code; CC BY 4.0 to original documentation, derived tables and figures, subject to third-party source terms. Raw external matrices, MSigDB GMT files, unpublished manuscript files, private review reports and credentials are not included. Supplementary methods and evidence-table descriptions are supporting documentation, not the full manuscript. Existing baseline directories are explicitly historical; current table/figure numbering is under the entry points above.
