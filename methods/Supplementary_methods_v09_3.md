# Supplementary methods and table guide

## Purpose and evidence limits

These materials accompany the manuscript “Cellular composition and the interpretation of candidate RNA readouts of trogocytosis in colorectal cancer”. The main patient analyses evaluate a fixed eight-gene module and a separate rare-transcript definition. They do not use a transfer-positive gold standard. Functional-screen analyses measure cancer-cell fitness. Negative associations, failed replication and post-review additions are retained.

## Gene definitions and temporal status

The fixed RNA module is CD4, PTPRC, CTLA4, PDCD1, HAVCR2, VSIR, LAG3 and CD38. These are not all experimentally validated recipient-cell RNA indicators of transfer. Table S1 lists the 39-member historical exploratory panel and its coverage, using descriptive evidence classes. Those classes do not reconstruct missing contemporaneous reasons for selection. The panel was not systematically assembled or prospectively registered. FADS2 is outside the original panel. CH25H, ATF3 and FADS2 functional comparisons were added after peer review; the resulting ATF3–SMYD2 association was selected after inspecting that exploratory screen.

GSE178341 T-lineage fractions combine author midpoint categories TCD4, TCD8, Tgd and TZBTB16 and divide their count by all tumour cells in each PID. Author-annotated exhausted fractions use TNK/ILC labels containing CXCL13 or PDCD1 and divide by T-lineage cells. The separate exhaustion-expression set is PDCD1, CTLA4, HAVCR2, LAG3, TIGIT, TOX and CXCL13. Patient-level inference uses 62 PIDs, not 64 PatientTypeID specimen units. Official epithelial CMScaller labels are retained; the 13 unclassified cases in the FDR-filtered sensitivity are not silently reassigned.

The seven-test group-comparison family contains three MMR comparisons (T-cell fraction, annotation-based exhaustion fraction and exhaustion-expression score) and four CMS comparisons (all-cell and epithelial scores, each using raw CMS calls and FDR-filtered calls). The separate four-test pathway family crosses all-cell/epithelial scores with TGF-β/EMT. These families are not combined with the post-review four-anchor functional screen.

## Matched-null calculations

The null generator compares the eight-gene score with 10,000 other eight-gene modules in the same 62 patients. It excludes all original panel genes, mitochondrial/ribosomal symbols, insufficiently detected or invariant genes, and ambiguous duplicate symbols. Its eligible pool has 21,904 genes. Three matching features measure mean expression, detection and TNK/ILC-versus-other specificity. Robustly standardised Euclidean distances define neighbourhoods; one gene per target is drawn without replacement within a module. The primary neighbourhood is k=50, with k=25 and k=100 sensitivity analyses. The random generator seed is 20260901.

If B=10,000 and r_obs is the frozen-module correlation, the upper-tail empirical P is (1 + count[r_b ≥ r_obs])/(B+1), using the implementation's numerical tie handling. The descriptive equal-tail two-sided P is min(1, 2 min[P_upper, P_lower]); it is centred on the empirical matched distribution. It is not the probability of a non-zero correlation. Thus P_upper≈0.493 and P_two≈0.985 are compatible with the observed coefficient near the 50.8th null percentile.

The complete-gene calibration has r_obs=0.7672181511. The archived contextual score gives 0.7665130568. Their patient score rank correlation is 0.9933016041. Both observed and random modules in the calibration use the same complete-gene CPM calculation; neither statistic is substituted into the other pipeline. This distinction explains why contextual plots and calibration can display the same rounded coefficient but differ in unrounded asymptotic P values.

Matching is assessed at the module-mean level rather than assumed perfect from nearest-neighbour selection. The frozen mean expression, detection and specificity lie at null percentiles 97.31, 32.94 and 90.86. Restricting to the nearest 10% of modules by the standardised distance between module-mean features is explicitly post hoc. Its upper-tail P=0.571 is a sensitivity result, not the primary test. Neither test quantifies causal mediation by cell composition.

## Independent compartment analysis

GSE132465 has 23 patients, 47,285 tumour cells, 17,469 epithelial cells and 16,739 T cells. Counts are summed by patient for all tumour cells and separately for the epithelial compartment. Expression is log2(CPM+0.25), with each gene standardised between patients within its compartment before averaging the eight genes. C10orf54 is mapped to VSIR; no other gene substitution is made. All eight genes vary in both compartments. The 10,000 patient bootstraps and permutations are inherited from the frozen analysis. The attenuation statistic is |ρ_all|−|ρ_epithelial|. Its interval, rather than different significance labels for its components, determines whether attenuation is established. The 95% interval −0.041 to 0.759 includes zero.

## CEACAM5 candidate definition

The epithelial-background markers are CDX2, AGR2, CLDN4, CLDN3, KRT7, KRT8, KRT18, KRT19, EPCAM, MUC1, CEACAM6, KRT20, MUC13, TFF3 and CDH1. The background score is log2(1 + the sum of their TPM values), not the arithmetic mean of individually log-transformed values. The primary candidate definition requires CEACAM5 TPM>1 and background score ≤ the 99.5th percentile in peripheral-blood leukocytes. Sensitivities use CEACAM5 >0, >1 and >2 and blood percentiles 99.0, 99.5 and 99.9. This is an operational RNA filter, not an ambient-RNA correction or a validated transfer-event detector.

Eight patients have paired normal and tumour specimens. Rates are calculated within each patient and tissue, then compared by enumerating all 2^8 sign flips of paired differences. The observed mean difference is −7.48×10−4 per leukocyte; its displayed sign-flip P=0.125 is patient-level. Tissue-level Wilson intervals in Figure 2 describe observed candidate proportions and do not account for all between-patient heterogeneity. They are not confidence intervals for biological event prevalence.

## TCGA score scaling and proportional hazards diagnostics

The frozen TCGA implementation standardises each gene and the composite score using sample standard deviations (ddof=1) in the expression cohort before excluding patients with unavailable or non-positive follow-up. Complete-case covariate filtering retains this score scale. The categorical-stage model uses R survival 3.8-6, Efron ties and stages I–IV as a factor. Schoenfeld-residual diagnostics use cox.zph with transform="rank". The module-term P is 0.873385 and the joint five-degree-of-freedom GLOBAL P is 0.701304. A minimum P across coefficient-specific tests is not the joint global diagnostic. Both the coefficient and the GLOBAL test are retained in the reproducibility record.

## Organoid partial correlations and multiple testing

The original selected-library table contains one record per each of 85 models. The same choices are used for every target; no model or library is selected to improve an association. The complete matrix has 16,562 genes. Ranks are residualised against an intercept and three ranked covariates: per-model genome-wide median LFC, AUC-ROC and AUC-PR. The design has rank four. For residual-rank correlation r, the approximate statistic is t=r sqrt[80/(1−r²)], evaluated against a t distribution with df=80. The parametric approximation is reported as such, not as an exact rank permutation test.

The anchor is removed before correction. Single-anchor families contain 16,561 tests; the pooled sensitivity contains 66,244 tests over four anchors. This distinction changes the SCD–VPS72 decision at q<0.05: q_single=0.042882 and q_pooled=0.085764. ATF3–SMYD2 has q_single=0.006661 and q_pooled=0.026642. A separate regression-t calculation was checked against the partial-correlation formula for a selected pair in each anchor analysis.

The SCD Hallmark rankings are unchanged by this correction, so the original full enrichment table is retained. It uses the original Hallmark 2023.2 GMT, 15–500 represented genes and 5,000 permutations with seed 20260831. The separate four-anchor exploratory extension is described below and reported in Table S10. Table S8 gives the full SCD results, including non-significant sets.

## DepMap adjustment and replication

The 26Q1 analysis selects records labelled Cell Line, giving 1,208 lines with 18,531 unique gene symbols. There are 63 CRC and 1,145 non-CRC lines. Common-essential genes are selected using all lines, requiring ≥90% finite values, median effect ≤−0.5 and effect <−0.5 in at least 80% of all lines. The tested genes are excluded. Both pair checks retain 1,179 genes. Missing common-essential values are median-imputed; standardisation uses sample standard deviations across all lines. Singular-value decomposition of this matrix gives the first five PC scores. This is not a CRC-only PCA.

Within CRC, each target is ranked and residualised against seven standardised covariates: global median effect, global median absolute deviation and five common-essential PC scores. Their residual correlation is the adjusted estimate. For each of 10,000 bootstrap samples, CRC lines are resampled, target ranks recomputed and the regression refitted. The global common-essential definition and PCA basis are held fixed; the interval is conditional on those choices. The two-sided empirical P compares the absolute observed coefficient with 10,000 correlations after permuting one residual vector and adds one to numerator and denominator. Seeds and resampling outputs are saved.

This procedure reproduced the original SCD–VPS72 result and gave ATF3–SMYD2 ρ=0.1563, 95% CI −0.1276–0.4209, P=0.2302. Lack of replication does not establish biological absence or a difference between platforms. Different model sets, culture conditions and perturbation libraries remain alternative explanations.

## Supplementary tables

- Table S1: 39 panel members, descriptive inclusion role, evidence limit and screen/module membership.
- Table S2: four-target CRC-versus-other-lineage effects, intervals, raw P and four-target BH q.
- Tables S3–S6: full SCD, CH25H, ATF3 and FADS2 organoid partner results; one row per non-self gene.
- Table S7: both DepMap pair checks, raw and adjusted associations, adjusted intervals and permutation P.
- Table S8: full original SCD Hallmark results for raw and QC-adjusted rankings.
- Table S9: combined four-anchor organoid results, including pooled BH q; compressed CSV to preserve all 66,244 rows.

- Table S10: all 200 pathway tests from the exploratory four-anchor Hallmark analysis, including both BH correction families.

Table definitions and CSV column names distinguish single-anchor and pooled correction. A descriptive proportion below −1 is never substituted for the organoid study's official depletion call. Source-data provenance is retained with the analysis scripts and audit records. The current supplementary tables and supporting code are included in release v1.2.0 (https://doi.org/10.5281/zenodo.23183934); the earlier v1.1.0 archive is retained as a historical baseline.


## Exploratory four-anchor Hallmark extension (Table S10)

The post-review extension uses the cached Hallmark 2020-labelled GMT for all four anchors, separate from the original Hallmark 2023.2 SCD analysis in Table S8. Version equivalence is not assumed. The local file contains 50 sets and is not redistributed. Its SHA-256 is 4275592957a1587652092bb398cf77216fde5b8daa2aedaa0e016f7d10bbdb81.

For each anchor, the complete adjusted rank-correlation list was sorted descending, with gene symbol as the deterministic tie-breaker. The anchor was removed, leaving 16,561 unique finite-valued genes. All 50 sets met the 15–500 represented-gene criterion. The weighted running sum uses absolute adjusted correlation, exponent 1, for hits and uniform decrements for misses. Positive/negative ES is the largest absolute excursion. NES divides observed ES by the mean absolute null ES of the same sign. There were 5,000 uniformly sampled gene-membership sets without replacement per pathway (seed 42, reset for each anchor), not patient-label permutations. Same-sign nominal P=(extreme count+1)/(same-sign null count+1). BH correction covers 50 tests per anchor and, as a sensitivity, all 200 tests. These are BH-adjusted nominal P values, not Broad GSEA pooled-NES FDR estimates.

The sparse running-sum implementation was checked against a full ranked-vector calculation for all 200 observed scores. All scores agreed to 1e-11. Full nominal P values, NES, represented-set sizes, null denominators, extreme counts, within-anchor BH and pooled BH are in Table S10. Finite-permutation variation affects exact P values; the chosen run was not repeated to seek significance. Minimum within-anchor adjusted P values were 0.495 (SCD), 0.416 (CH25H), 0.257 (ATF3), and 0.221 (FADS2). None met 0.05. This exploratory gene-set analysis does not test equivalence, measure pathway activation, or directly observe membrane transfer.
