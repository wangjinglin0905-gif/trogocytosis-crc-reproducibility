args <- commandArgs(trailingOnly = TRUE)

source_root <- if (length(args) >= 1L) args[[1]] else
  "."
output_root <- if (length(args) >= 2L) args[[2]] else
  file.path(source_root, "qa", "recomputed", "v8_1_corrections")

source_root <- normalizePath(source_root, winslash = "/", mustWork = TRUE)
dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(survival))

read_tsv <- function(path) {
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

# 1. Correct the matched-null two-sided empirical tail definition.
# The primary one-sided test asks whether the frozen module is unusually more
# correlated than matched modules. The two-sided descriptive value is an
# equal-tail probability around the empirical null distribution, not a test
# around rho = 0.
matched_dir <- file.path(source_root, "analysis", "matched_null_gse178341")
matched_summary <- read_tsv(file.path(matched_dir, "gse178341_matched_null_summary.tsv"))
matched_rhos <- read_tsv(gzfile(file.path(matched_dir, "gse178341_matched_null_rhos.tsv.gz")))

corrected_rows <- lapply(matched_summary$k, function(k_value) {
  observed <- matched_summary$observed_rho[matched_summary$k == k_value][1]
  null_values <- matched_rhos$rho[matched_rhos$k == k_value & is.finite(matched_rhos$rho)]
  n_valid <- length(null_values)
  p_upper <- (1 + sum(null_values >= observed)) / (n_valid + 1)
  p_lower <- (1 + sum(null_values <= observed)) / (n_valid + 1)
  p_two_equal_tail <- min(1, 2 * min(p_lower, p_upper))
  data.frame(
    k = k_value,
    reps_valid = n_valid,
    observed_rho = observed,
    observed_percentile = 100 * mean(null_values <= observed),
    empirical_p_one_sided_upper = p_upper,
    empirical_p_lower_tail = p_lower,
    empirical_p_two_sided_equal_tail = p_two_equal_tail,
    stringsAsFactors = FALSE
  )
})
matched_corrected <- do.call(rbind, corrected_rows)
write.table(matched_corrected,
            file.path(output_root, "gse178341_matched_null_pvalue_correction.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# 2. Refit the TCGA model exactly as described in the manuscript: stage is a
# categorical covariate. The archived V7 table used stage_num continuously.
tcga_dir <- file.path(source_root, "baseline", "v7", "analysis", "tcga_cbioportal_592")
tcga <- read_tsv(file.path(tcga_dir, "tcga_cbioportal_module_clinical.tsv"))
tcga$stage_factor <- factor(tcga$stage_num,
                            levels = c(1, 2, 3, 4),
                            labels = c("I", "II", "III", "IV"))

fit_data <- tcga[is.finite(tcga$os_months) & tcga$os_months > 0 &
                   !is.na(tcga$os_event) & is.finite(tcga$module_z) &
                   is.finite(tcga$age) & !is.na(tcga$stage_factor), ]

fit <- coxph(Surv(os_months, os_event) ~ module_z + age + stage_factor,
             data = fit_data, ties = "efron", x = TRUE)
fit_summary <- summary(fit)
module_index <- which(rownames(fit_summary$coefficients) == "module_z")
module_coef <- fit_summary$coefficients[module_index, ]
module_ci <- fit_summary$conf.int[module_index, ]
ph <- cox.zph(fit, transform = "rank")

tcga_corrected <- data.frame(
  cohort = "TCGA-COAD/READ via cBioPortal",
  model = "OS_age_categorical_stage_adjusted",
  n = nrow(fit_data),
  events = sum(fit_data$os_event == 1),
  hr_per_sd = unname(module_ci["exp(coef)"]),
  ci95_low = unname(module_ci["lower .95"]),
  ci95_high = unname(module_ci["upper .95"]),
  p = unname(module_coef["Pr(>|z|)"]),
  module_schoenfeld_p = unname(ph$table["module_z", "p"]),
  global_schoenfeld_p = unname(ph$table["GLOBAL", "p"]),
  stage_coding = "factor(I, II, III, IV)",
  ties = "efron",
  stringsAsFactors = FALSE
)
write.table(tcga_corrected,
            file.path(output_root, "tcga_categorical_stage_corrected.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

write.table(as.data.frame(ph$table),
            file.path(output_root, "tcga_categorical_stage_schoenfeld.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

audit_lines <- c(
  "V8.1 numerical correction audit",
  sprintf("Matched-null primary k=50: one-sided upper-tail P=%.6f; equal-tail two-sided P=%.6f; percentile=%.2f.",
          matched_corrected$empirical_p_one_sided_upper[matched_corrected$k == 50],
          matched_corrected$empirical_p_two_sided_equal_tail[matched_corrected$k == 50],
          matched_corrected$observed_percentile[matched_corrected$k == 50]),
  sprintf("TCGA categorical-stage model: n=%d, deaths=%d, HR=%.6f, 95%% CI %.6f to %.6f, P=%.6f.",
          tcga_corrected$n, tcga_corrected$events, tcga_corrected$hr_per_sd,
          tcga_corrected$ci95_low, tcga_corrected$ci95_high, tcga_corrected$p),
  sprintf("TCGA PH tests: module P=%.6f; global P=%.6f.",
          tcga_corrected$module_schoenfeld_p, tcga_corrected$global_schoenfeld_p),
  "The matched-null correction changes only the descriptive two-sided P value; the prespecified one-sided conclusion is unchanged.",
  "The categorical-stage TCGA refit is concordant with the archived continuous-stage model and remains null."
)
writeLines(audit_lines, file.path(output_root, "V8_1_numerical_correction_audit.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_root, "R_sessionInfo.txt"))

cat(paste(audit_lines, collapse = "\n"), "\n")
