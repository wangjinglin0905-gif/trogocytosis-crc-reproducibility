# V8 argument and terminology ledger

## One-sentence argument

In colorectal cancer, SCD is a reproducible tumour-cell fitness dependency across organoid and two-dimensional CRISPR screens, whereas prespecified matched-null and independent full-cell single-cell analyses indicate that the frozen mixed-marker RNA module principally reports immune composition; neither the RNA analyses nor the non-replicated SCD-VPS72 association establishes a trogocytosis event or mechanism.

## Claim boundary

- Supported: SCD is a CRC fitness dependency in the analysed screening platforms.
- Supported: the frozen eight-gene RNA score is strongly coupled to T-cell abundance and is not more T-cell-correlated than matched immune-lineage-specific random modules in GSE178341.
- Compatible cross-cohort evidence: GSE132465 reproduces the positive whole-tumour association and point attenuation after epithelial restriction, but the attenuation confidence interval includes zero.
- Not supported: SCD is a regulator of trogocytosis; VPS72 is a cross-platform SCD co-dependency; the RNA score detects protein transfer; the score predicts treatment response or prognosis.

## Locked terminology

| Canonical term | First-use definition | Variants avoided | Decision |
|---|---|---|---|
| trogocytosis event | contact-dependent acquisition of donor-cell membrane material or membrane-associated protein by a recipient cell | event-like cell, inferred event | Reserve “event” for event-resolving assays |
| frozen eight-gene module | CD4, PTPRC, CTLA4, PDCD1, HAVCR2, VSIR, LAG3 and CD38 | 8-gene signature, trogocytosis signature | “Module”, not a validated signature |
| RNA-defined candidate | a cell meeting an RNA filter without proof of membrane transfer | high-confidence event | Never call an event |
| tumour-wide pseudobulk | patient-level aggregation across all author-annotated tumour cells | whole tumour bulk | Use consistently for single-cell aggregation |
| epithelial pseudobulk | patient-level aggregation restricted to author-annotated epithelial cells | epithelial-only score | Use consistently |
| T-cell fraction | author-annotated T cells divided by all analysed tumour cells for each patient | T-cell abundance | Use fraction for the numeric variable |
| matched-null module | eight genes matched target-by-target for mean expression, cell detection and TNK/ILC specificity | random signature | Use “module” and state matching variables |
| QC/common-essentiality-adjusted partial Spearman correlation | correlation of residualized within-CRC gene-effect ranks after adjustment for global median, MAD and five common-essential PCs | corrected co-dependency | Use the full definition at first use |
| compatible replication | prespecified primary association and point attenuation criteria met, but attenuation CI crosses zero | validation, strong replication | Use only for GSE132465 result |
| non-replication | prespecified direction/significance/CI gate not met | negative validation | Use for DepMap SCD-VPS72 |
| SCD dependency | loss of fitness after SCD knockout | trogocytosis dependency | Never infer trogocytosis mechanism from fitness alone |
| VSIR | V-set immunoregulatory receptor | C10orf54 | Map legacy C10orf54 to VSIR and audit the mapping |

