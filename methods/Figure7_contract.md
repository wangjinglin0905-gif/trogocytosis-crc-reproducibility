# Figure 7 contract

- **Core conclusion:** A mixed immune/acquired-marker RNA module tracks T-cell composition rather than a specific membrane-transfer event, while the exploratory SCD-VPS72 organoid association fails independent two-dimensional replication.
- **Evidence chain:**
  - panel a: the observed GSE178341 module-T-cell correlation lies near the centre of expression/detection/specificity-matched null distributions;
  - panel b: matching-feature balance and a closest-10% sensitivity show that this conclusion is not rescued by tighter module-level balance;
  - panel c: DepMap raw and QC/common-essentiality-adjusted SCD-VPS72 estimates are negative but imprecise and cross zero;
  - panel d: GSE132465 reproduces a strong tumour-wide correlation that attenuates and becomes non-significant after epithelial restriction.
- **Archetype:** quantitative grid, with panels a and d as the principal evidence panels.
- **Backend:** R only for plotting, assembly, export and visual QA.
- **Analysis units:** patients for GSE178341 and GSE132465; cell lines for DepMap.
- **Statistics:** 10,000 matched-null modules; 10,000 bootstraps/permutations; effect sizes and 95% intervals reported; no cell-level P values.
- **Export:** 180 mm wide, white background, editable SVG/PDF, 600-dpi TIFF and 300-dpi PNG preview; lowercase panel labels.
- **Review risks:** imperfect module-mean matching, post-hoc nature of the closest-10% sensitivity, wide attenuation CI in n=23, and non-comparability of organoid and DepMap effect estimates must be explicit.

