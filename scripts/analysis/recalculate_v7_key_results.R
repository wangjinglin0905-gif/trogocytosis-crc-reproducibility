args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: recalculate_v7_key_results.R <root> <outdir>")

root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
outdir <- args[[2]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(jsonlite))

read_tsv <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
as_num <- function(x) as.numeric(x)

cox_summary <- function(fit, module_name = "module_z") {
  s <- summary(fit)
  coef_row <- s$coefficients[module_name, ]
  ci_row <- s$conf.int[module_name, ]
  zph <- cox.zph(fit, transform = "rank")$table
  list(
    n = unname(fit$n),
    events = unname(fit$nevent),
    hr = unname(ci_row[["exp(coef)"]]),
    ci95_low = unname(ci_row[["lower .95"]]),
    ci95_high = unname(ci_row[["upper .95"]]),
    p = unname(coef_row[["Pr(>|z|)"]]),
    module_ph_p_rank = unname(zph[module_name, "p"]),
    global_ph_p_rank = unname(zph["GLOBAL", "p"])
  )
}

results <- list()

# Organoid anchors and subgroup tests
org <- file.path(root, "analysis", "organoid")
panel <- read.csv(file.path(org, "organoid_panel_dependency_summary.csv"), check.names = FALSE)
screen <- read.csv(gzfile(file.path(org, "organoid_crc_best_library_screen.csv.gz")), check.names = FALSE)
ann <- read.csv(file.path(org, "organoid_model_annotations.csv"), check.names = FALSE)
screen <- merge(screen, ann[, c("sample_ID", "individual_ID", "msStatus", "CMS_prediction", "KRAS_class")],
                by = "sample_ID", all.x = TRUE, sort = FALSE)

results$organoid_anchors <- lapply(c("SCD", "EGFR"), function(gene) {
  values <- screen$LFC[screen$gene == gene]
  stored <- panel[panel$gene == gene, ]
  list(gene = gene, n = sum(!is.na(values)), median_lfc_recomputed = median(values, na.rm = TRUE),
       median_lfc_stored = stored$median_LFC[[1]], official_pct_depleted = stored$pct_depleted[[1]])
})

cms_order <- c("CMS1", "CMS2", "CMS3", "CMS4")
scd <- screen[screen$gene == "SCD", ]
cms_groups <- lapply(cms_order, function(cms) scd$LFC[scd$CMS_prediction == cms & !is.na(scd$CMS_prediction)])
results$organoid_scd_cms <- list(
  counts = sapply(cms_groups, length),
  medians = sapply(cms_groups, median, na.rm = TRUE),
  kruskal_p = kruskal.test(cms_groups)$p.value
)

msi_rows <- lapply(sort(unique(panel$gene)), function(gene) {
  subset <- screen[screen$gene == gene, ]
  msi <- subset$LFC[subset$msStatus == "MSI"]
  mss <- subset$LFC[subset$msStatus == "MSS"]
  data.frame(gene = gene, n_msi = sum(!is.na(msi)), n_mss = sum(!is.na(mss)),
             p = wilcox.test(msi, mss, exact = FALSE, correct = FALSE)$p.value)
})
msi_table <- do.call(rbind, msi_rows)
msi_table$q <- p.adjust(msi_table$p, method = "BH")
results$organoid_msi <- list(
  n_msi = sum(ann$msStatus == "MSI", na.rm = TRUE),
  n_mss = sum(ann$msStatus == "MSS", na.rm = TRUE),
  genes = lapply(c("SCD", "EGFR"), function(gene) as.list(msi_table[msi_table$gene == gene, ])),
  n_fdr_lt_005 = sum(msi_table$q < 0.05)
)

# DepMap cross-platform calculations
dep <- file.path(root, "analysis", "depmap_26Q1")
cross <- read.csv(file.path(dep, "depmap_organoid_cross_platform_panel.csv"), check.names = FALSE)
full_cor <- cor.test(cross$median_organoid, cross$median_crc_2d, method = "spearman", exact = FALSE)
legacy <- cross[cross$gene != "CD274", ]
legacy_cor <- cor.test(legacy$median_organoid, legacy$median_crc_2d, method = "spearman", exact = FALSE)
context <- read.csv(file.path(dep, "depmap_target_context_summary.csv"), check.names = FALSE)
anchor_context <- lapply(c("SCD", "EGFR"), function(gene) {
  row <- context[context$gene == gene, ]
  as.list(row[, c("gene", "median_crc_2d", "pct_crc_lt_minus_1", "median_difference_crc_minus_other",
                  "median_difference_ci95_lo", "median_difference_ci95_hi", "mannwhitney_p_crc_vs_other")])
})
results$cross_platform <- list(
  n_full = nrow(cross), rho_full = unname(full_cor$estimate), p_full = full_cor$p.value,
  n_legacy = nrow(legacy), rho_legacy = unname(legacy_cor$estimate), p_legacy = legacy_cor$p.value,
  anchors = anchor_context
)

# SCD co-dependency and Hallmark sensitivity outputs
codep <- file.path(root, "analysis", "scd_codependency_qc")
rank <- read.csv(file.path(codep, "scd_codependency_qc_adjusted.csv"), check.names = FALSE)
gsea <- read.csv(file.path(codep, "scd_hallmark_raw_vs_qc.csv"), check.names = FALSE)
target_rows <- rank[rank$gene %in% c("VPS72", "CDX1", "WLS"), ]
results$scd_codependency <- list(
  n_genes = nrow(rank),
  n_gene_fdr_lt_005 = sum(rank$FDR_BH < 0.05),
  target_rows = split(target_rows, target_rows$gene),
  n_hallmark_fdr_lt_005 = sum(gsea$FDR_BH < 0.05),
  dna_repair = split(gsea[gsea$pathway == "HALLMARK_DNA_REPAIR", ], gsea$ranking[gsea$pathway == "HALLMARK_DNA_REPAIR"]),
  cholesterol_homeostasis = split(gsea[gsea$pathway == "HALLMARK_CHOLESTEROL_HOMEOSTASIS", ], gsea$ranking[gsea$pathway == "HALLMARK_CHOLESTEROL_HOMEOSTASIS"])
)

# Survival recalculation from frozen patient-level tables
gse <- read_tsv(file.path(root, "analysis", "gse39582_recalculation", "gse39582_module_clinical.tsv"))
gse_univ_data <- gse[gse$rfs_months > 0 & complete.cases(gse[, c("rfs_months", "rfs_event", "module_z")]), ]
gse_adj_vars <- c("rfs_months", "rfs_event", "module_z", "age", "mmrd", "chemotherapy_yes", "stage_2", "stage_3", "stage_4")
gse_adj_data <- gse[gse$rfs_months > 0 & complete.cases(gse[, gse_adj_vars]), ]
gse_univ_fit <- coxph(Surv(rfs_months, rfs_event) ~ module_z, data = gse_univ_data, ties = "efron", x = TRUE)
gse_adj_fit <- coxph(Surv(rfs_months, rfs_event) ~ module_z + age + mmrd + chemotherapy_yes + stage_2 + stage_3 + stage_4,
                     data = gse_adj_data, ties = "efron", x = TRUE)
results$gse39582 <- list(rfs_univariable = cox_summary(gse_univ_fit), rfs_adjusted = cox_summary(gse_adj_fit), os_analysable = FALSE)

tcga <- read_tsv(file.path(root, "analysis", "tcga_cbioportal_592", "tcga_cbioportal_module_clinical.tsv"))
tcga$stage_factor <- factor(tcga$stage_num, levels = c(1, 2, 3, 4))
tcga_univ_data <- tcga[tcga$os_months > 0 & complete.cases(tcga[, c("os_months", "os_event", "module_z")]), ]
tcga_adj_data <- tcga[tcga$os_months > 0 & complete.cases(tcga[, c("os_months", "os_event", "module_z", "age", "stage_factor")]), ]
tcga_univ_fit <- coxph(Surv(os_months, os_event) ~ module_z, data = tcga_univ_data, ties = "efron", x = TRUE)
tcga_adj_fit <- coxph(Surv(os_months, os_event) ~ module_z + age + stage_factor, data = tcga_adj_data, ties = "efron", x = TRUE)
results$tcga <- list(expression_n = sum(!is.na(tcga$module_z)), os_univariable = cox_summary(tcga_univ_fit),
                     os_age_stage_adjusted = cox_summary(tcga_adj_fit))

# Patient-paired candidate-cell recalculation
leuko <- file.path(root, "analysis", "gse146771")
rates <- read.csv(file.path(leuko, "gse146771_candidate_rates_by_patient_tissue_lineage.csv"), check.names = FALSE)
all_rates <- rates[rates$lineage == "All leukocytes", ]
wide <- reshape(all_rates[, c("Sample", "Tissue", "rate")], idvar = "Sample", timevar = "Tissue", direction = "wide")
paired <- wide[complete.cases(wide[, c("rate.N", "rate.T")]), ]
differences <- paired$rate.T - paired$rate.N
sign_matrix <- expand.grid(rep(list(c(-1, 1)), length(differences)))
permuted <- apply(sign_matrix, 1, function(s) abs(mean(differences * as.numeric(s))))
sign_p <- mean(permuted >= abs(mean(differences)) - 1e-15)
tissue <- read.csv(file.path(leuko, "gse146771_candidate_tissue_descriptive.csv"), check.names = FALSE)
sensitivity <- read.csv(file.path(leuko, "gse146771_threshold_sensitivity.csv"), check.names = FALSE)
results$gse146771 <- list(
  paired_patients = nrow(paired), mean_difference_T_minus_N = mean(differences), exact_sign_flip_p = sign_p,
  tissue_counts = split(tissue, tissue$tissue),
  threshold_candidate_range = lapply(split(sensitivity$total_candidates, sensitivity$CEACAM5_TPM_threshold), range)
)

# Patient-level GSE178341 recalculation
g178 <- read_tsv(file.path(root, "analysis", "gse178341_recalculation", "gse178341_patient_scores.tsv"))
cms4 <- g178$prediction_raw == "CMS4"
cor_pairs <- list(
  module_all_vs_TGFb_all = c("module_all", "TGFb_all"),
  module_all_vs_EMT_all = c("module_all", "EMT_all"),
  module_Epi_vs_TGFb_Epi = c("module_Epi", "TGFb_Epi"),
  module_Epi_vs_EMT_Epi = c("module_Epi", "EMT_Epi"),
  module_all_vs_T_cell_fraction_all = c("module_all", "T_cell_fraction_all"),
  module_Epi_vs_T_cell_fraction_all = c("module_Epi", "T_cell_fraction_all")
)
correlations <- lapply(cor_pairs, function(pair) {
  test <- cor.test(g178[[pair[[1]]]], g178[[pair[[2]]]], method = "spearman", exact = FALSE)
  list(rho = unname(test$estimate), p = test$p.value)
})
mmrd <- g178$MMRStatus == "MMRd"
mmr_columns <- c("T_cell_fraction_all", "exhausted_label_fraction_T", "exhaustion_expression_TNKILC")
mmr_tests <- lapply(mmr_columns, function(column) {
  test <- wilcox.test(g178[[column]][mmrd], g178[[column]][!mmrd], exact = FALSE, correct = FALSE)
  list(variable = column, median_mmrd = median(g178[[column]][mmrd], na.rm = TRUE),
       median_mmrp = median(g178[[column]][!mmrd], na.rm = TRUE), p = test$p.value)
})
fdr_classified <- !is.na(g178$prediction_FDR05) & g178$prediction_FDR05 != "Unclassified"
cms_p_values <- c(
  wilcox.test(g178$module_all[cms4], g178$module_all[!cms4], exact = FALSE, correct = FALSE)$p.value,
  wilcox.test(g178$module_all[g178$prediction_FDR05 == "CMS4"],
              g178$module_all[fdr_classified & g178$prediction_FDR05 != "CMS4"], exact = FALSE, correct = FALSE)$p.value,
  wilcox.test(g178$module_Epi[cms4], g178$module_Epi[!cms4], exact = FALSE, correct = FALSE)$p.value,
  wilcox.test(g178$module_Epi[g178$prediction_FDR05 == "CMS4"],
              g178$module_Epi[fdr_classified & g178$prediction_FDR05 != "CMS4"], exact = FALSE, correct = FALSE)$p.value
)
group_q <- p.adjust(c(sapply(mmr_tests, function(x) x$p), cms_p_values), method = "BH")
for (i in seq_along(mmr_tests)) mmr_tests[[i]]$q <- group_q[[i]]
results$gse178341 <- list(
  patients = length(unique(g178$PID)), mmrd = sum(mmrd), mmrp = sum(!mmrd),
  cms_raw_counts = as.list(table(g178$prediction_raw)),
  cms4_whole_p = cms_p_values[[1]], cms4_whole_q = group_q[[4]],
  cms4_whole_fdr_call_p = cms_p_values[[2]], cms4_whole_fdr_call_q = group_q[[5]],
  cms4_epithelial_p = cms_p_values[[3]], cms4_epithelial_q = group_q[[6]],
  cms4_epithelial_fdr_call_p = cms_p_values[[4]], cms4_epithelial_fdr_call_q = group_q[[7]],
  correlations = correlations, mmr_tests = mmr_tests
)

output <- file.path(outdir, "v7_key_results_recalculated.json")
write_json(results, output, auto_unbox = TRUE, pretty = TRUE, digits = 16, na = "null")
cat(output, "\n")
cat(toJSON(results, auto_unbox = TRUE, pretty = TRUE, digits = 8, na = "null"), "\n")
