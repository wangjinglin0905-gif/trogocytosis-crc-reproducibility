# Current submission-display figures

The authoritative approved exports are `figures/frontiers_v8_6/Figure_1` through `Figure_7`, each with `.png`, `.tiff`, `.pdf` and `.svg`. They were copied byte-for-byte from the approved local submission assets, not redrawn for the V8.6 wording revision. The older `figures/final/` files remain historical v1.0.0 assets. The assets TSV records hashes for the current exports, portable plotting sources and supplemental display tables.

## Source map and regeneration

Run commands from the repository root using the existing R dependencies. Use new output directories; do not overwrite approved exports. Reproduced vector metadata or jitter can depend on the graphics environment; compare source values and labels, not only image hashes.

| Figures | Portable script | Inputs | Generated stems |
|---|---|---|---|
| 1 | `scripts/figures/frontiers_v8_6/make_figure_1.R` | `baseline/v7/analysis/organoid/` | `Figure_1` |
| 2, 5 | `scripts/figures/frontiers_v8_6/make_figures_2_5.R` | `baseline/v7/analysis/depmap_26Q1/`, `baseline/v7/analysis/gse146771/` | `Fig2_DepMap_validation_v7`, `Fig5_leukocyte_candidates_v7` |
| 3, 6, 7 | `scripts/figures/frontiers_v8_6/make_figures_3_6_7.R` | `baseline/v7/analysis/`, `analysis/v8_1_corrections/`, `analysis/matched_null_gse178341/`, `analysis/depmap_vps72_replication/`, `analysis/gse132465_replication/` | `Fig3_SCD_codependency_v8_1`, `Fig6_GSE178341_patient_level_v8_1`, `Fig7_composition_controls_and_vps72_replication_v8_1` |
| 4 | `scripts/figures/frontiers_v8_6/make_figure_4.R` | `baseline/v7/analysis/gse39582_recalculation/`, `baseline/v7/analysis/tcga_cbioportal_592/`, `analysis/v8_1_corrections/` | `Fig4_bulk_survival_v8_1` |

```bash
Rscript scripts/figures/frontiers_v8_6/make_figure_1.R baseline/v7 qa/recomputed/frontiers_figure1
Rscript scripts/figures/frontiers_v8_6/make_figures_2_5.R baseline/v7 qa/recomputed/frontiers_figures2_5
Rscript scripts/figures/frontiers_v8_6/make_figures_3_6_7.R . analysis/v8_1_corrections qa/recomputed/frontiers_figures3_6_7
Rscript scripts/figures/frontiers_v8_6/make_figure_4.R . analysis/v8_1_corrections qa/recomputed/frontiers_figure4
```

The two multi-figure scripts were mechanically restricted to the relevant panels from existing display scripts; superseded figure blocks were omitted. The single-figure scripts were retained without algorithm changes. Script-internal historical stem/version names identify numerical provenance, not additional results. Existing approved figures were visually checked and compared with their legends and source records for this release, but R re-execution was not performed.

Additional full-range Figure 1 source tables and the frozen Figure 4 forest table are in `figures/source_data/frontiers_v8_6/`. Their earlier subgroup test/display checks are historical evidence, not new statistical tests. Figure 2C's CI/P distinction, Figure 4's continuous-Cox rather than log-rank annotations, and Figure 7's one-sided versus equal-tail P definitions remain unchanged.
