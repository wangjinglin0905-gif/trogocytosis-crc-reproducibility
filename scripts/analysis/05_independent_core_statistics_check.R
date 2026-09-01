options(stringsAsFactors = FALSE, warn = 1)

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args) >= 1L) args[[1]] else ".",
                      winslash = "/", mustWork = TRUE)
qa_dir <- if (length(args) >= 2L) args[[2]] else file.path(root, "qa", "recomputed")
dir.create(qa_dir, recursive = TRUE, showWarnings = FALSE)

rows <- list()
add_check <- function(analysis, statistic, recomputed, reported,
                      tolerance = 1e-10) {
  difference <- abs(as.numeric(recomputed) - as.numeric(reported))
  rows[[length(rows) + 1L]] <<- data.frame(
    analysis = analysis,
    statistic = statistic,
    recomputed_in_R = as.numeric(recomputed),
    reported_by_primary_script = as.numeric(reported),
    absolute_difference = difference,
    tolerance = tolerance,
    pass = is.finite(difference) && difference <= tolerance
  )
}

# 1. GSE178341 primary k=50 matched-null.
mn_dir <- file.path(root, "analysis", "matched_null_gse178341")
mn_pat <- read.delim(file.path(mn_dir, "gse178341_patient_scores_matched_null.tsv"),
                     check.names = FALSE)
mn_sum <- read.delim(file.path(mn_dir, "gse178341_matched_null_summary.tsv"),
                     check.names = FALSE)
mn_sum <- mn_sum[mn_sum$k == 50, ]
mn_null <- read.delim(gzfile(file.path(mn_dir, "gse178341_matched_null_rhos.tsv.gz")),
                      check.names = FALSE)
mn_null <- mn_null$rho[mn_null$k == 50]
mn_rho <- cor(mn_pat$frozen_module_score_complete_gene_cpm,
              mn_pat$t_cell_fraction, method = "spearman")
mn_p_one <- (1 + sum(mn_null >= mn_rho)) / (1 + length(mn_null))
add_check("GSE178341 matched null", "observed Spearman rho", mn_rho,
          mn_sum$observed_rho)
add_check("GSE178341 matched null", "null median", median(mn_null),
          mn_sum$null_median)
add_check("GSE178341 matched null", "empirical one-sided P", mn_p_one,
          mn_sum$empirical_p_one_sided)

# 2. DepMap SCD-VPS72 raw and QC/common-essential-adjusted estimates.
dm_dir <- file.path(root, "analysis", "depmap_vps72_replication")
dm_lines <- read.delim(file.path(dm_dir, "depmap_scd_vps72_crc_lines.tsv"),
                       check.names = FALSE)
dm_res <- read.delim(file.path(dm_dir, "depmap_scd_vps72_results.tsv"),
                     check.names = FALSE)
dm_boot <- read.delim(gzfile(file.path(dm_dir, "depmap_scd_vps72_bootstrap.tsv.gz")),
                      check.names = FALSE)
dm_perm <- read.delim(gzfile(file.path(dm_dir, "depmap_scd_vps72_permutation.tsv.gz")),
                      check.names = FALSE)
raw_row <- dm_res[dm_res$model == "raw Spearman", ]
adj_row <- dm_res[grepl("partial Spearman", dm_res$model), ]
raw_rho <- cor(dm_lines$SCD_gene_effect, dm_lines$VPS72_gene_effect,
               method = "spearman")
adj_rho <- cor(dm_lines$SCD_rank_residual, dm_lines$VPS72_rank_residual,
               method = "pearson")
adj_p <- (1 + sum(abs(dm_perm$permuted_rho) >= abs(adj_rho))) /
  (1 + nrow(dm_perm))
adj_ci <- quantile(dm_boot$adjusted_partial_spearman, c(0.025, 0.975),
                   type = 7, na.rm = TRUE)
add_check("DepMap 26Q1", "raw Spearman rho", raw_rho, raw_row$rho)
add_check("DepMap 26Q1", "adjusted partial Spearman", adj_rho, adj_row$rho)
add_check("DepMap 26Q1", "adjusted permutation P", adj_p, adj_row$p)
add_check("DepMap 26Q1", "adjusted bootstrap CI low", adj_ci[1],
          adj_row$ci95_low)
add_check("DepMap 26Q1", "adjusted bootstrap CI high", adj_ci[2],
          adj_row$ci95_high)

# 3. GSE132465 full-cell replication.
sc_dir <- file.path(root, "analysis", "gse132465_replication")
sc_pat <- read.delim(file.path(sc_dir, "gse132465_patient_scores.tsv"),
                     check.names = FALSE)
sc_sum <- read.delim(file.path(sc_dir, "gse132465_composition_summary.tsv"),
                     check.names = FALSE)
sc_boot <- read.delim(gzfile(file.path(sc_dir, "gse132465_bootstrap.tsv.gz")),
                      check.names = FALSE)
sc_perm <- read.delim(gzfile(file.path(sc_dir, "gse132465_permutation.tsv.gz")),
                      check.names = FALSE)
rho_whole <- cor(sc_pat$whole_tumour_module_score, sc_pat$t_cell_fraction,
                 method = "spearman")
rho_epi <- cor(sc_pat$epithelial_module_score, sc_pat$t_cell_fraction,
               method = "spearman")
# A 1e-12 tolerance preserves rank-identical permutation ties after TSV
# serialisation across Python and R floating-point implementations.
p_whole <- (1 + sum(abs(sc_perm$rho_whole) >= abs(rho_whole) - 1e-12)) /
  (1 + nrow(sc_perm))
p_epi <- (1 + sum(abs(sc_perm$rho_epithelial) >= abs(rho_epi) - 1e-12)) /
  (1 + nrow(sc_perm))
att <- abs(rho_whole) - abs(rho_epi)
att_ci <- quantile(sc_boot$attenuation_abs_rho, c(0.025, 0.975),
                   type = 7, na.rm = TRUE)
add_check("GSE132465", "whole-tumour Spearman rho", rho_whole,
          sc_sum$rho_whole)
add_check("GSE132465", "whole-tumour permutation P", p_whole,
          sc_sum$permutation_p_whole)
add_check("GSE132465", "epithelial Spearman rho", rho_epi,
          sc_sum$rho_epithelial)
add_check("GSE132465", "epithelial permutation P", p_epi,
          sc_sum$permutation_p_epithelial)
add_check("GSE132465", "absolute-rho attenuation", att,
          sc_sum$attenuation_abs_rho)
add_check("GSE132465", "attenuation bootstrap CI low", att_ci[1],
          sc_sum$attenuation_bootstrap_ci_low)
add_check("GSE132465", "attenuation bootstrap CI high", att_ci[2],
          sc_sum$attenuation_bootstrap_ci_high)

checks <- do.call(rbind, rows)
write.table(checks, file.path(qa_dir, "independent_core_statistics_check.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
stopifnot(all(checks$pass))

sink(file.path(qa_dir, "independent_core_statistics_check.log"))
cat("Independent R recalculation of V8 core statistics\n")
cat("Run UTC:", format(Sys.time(), tz = "UTC", usetz = TRUE), "\n\n")
print(checks, row.names = FALSE)
cat("\nSession information:\n")
print(sessionInfo())
sink()
print(checks, row.names = FALSE)
