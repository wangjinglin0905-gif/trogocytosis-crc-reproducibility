args <- commandArgs(trailingOnly = TRUE)
source_root <- if (length(args) >= 1L) args[[1]] else
  "."
correction_root <- if (length(args) >= 2L) args[[2]] else
  file.path(source_root, "analysis", "v8_1_corrections")
outdir <- if (length(args) >= 3L) args[[3]] else
  file.path(source_root, "qa", "recomputed", "figures_v8_1")

source_root <- normalizePath(source_root, winslash = "/", mustWork = TRUE)
correction_root <- normalizePath(correction_root, winslash = "/", mustWork = TRUE)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
source_dir <- file.path(outdir, "source_data")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(scales))
suppressPackageStartupMessages(library(grid))

blue <- "#3C5488"
light_blue <- "#91B7D6"
red <- "#E64B35"
orange <- "#F39B7F"
purple <- "#7E57C2"
grey <- "#737373"
light_grey <- "#D9D9D9"
black <- "#222222"

theme_v8_1 <- function(base_size = 7.5) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1.1,
                                hjust = 0, margin = margin(b = 2)),
      plot.subtitle = element_text(size = base_size - 0.45, colour = grey,
                                   hjust = 0, margin = margin(b = 3)),
      plot.caption = element_text(size = base_size - 1.2, colour = grey,
                                  hjust = 0, margin = margin(t = 3)),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.45, colour = black),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size),
      legend.title = element_blank(),
      legend.text = element_text(size = base_size - 0.7),
      legend.key.height = unit(3.5, "mm"),
      plot.tag = element_text(face = "bold", size = base_size + 2),
      plot.tag.position = c(0, 1),
      plot.margin = margin(6, 7, 6, 7)
    )
}

export_plot <- function(plot, stem, width_mm, height_mm) {
  ggsave(file.path(outdir, paste0(stem, ".png")), plot,
         width = width_mm, height = height_mm, units = "mm",
         dpi = 300, device = ragg::agg_png, background = "white")
  ggsave(file.path(outdir, paste0(stem, ".tiff")), plot,
         width = width_mm, height = height_mm, units = "mm",
         dpi = 600, device = "tiff", compression = "lzw", bg = "white")
  ggsave(file.path(outdir, paste0(stem, ".pdf")), plot,
         width = width_mm / 25.4, height = height_mm / 25.4, units = "in",
         device = cairo_pdf, bg = "white")
  svg(file.path(outdir, paste0(stem, ".svg")),
      width = width_mm / 25.4, height = height_mm / 25.4,
      family = "Arial", onefile = TRUE, bg = "white")
  print(plot)
  dev.off()
}

read_tsv <- function(path) {
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}

v7_analysis <- file.path(source_root, "baseline", "v7", "analysis")

# Figure 3: decouple the lower panel's long pathway labels from the upper-row
# column geometry. This lets panels a and b share the upper row equally and
# prevents panel a's title from encroaching on panel b.
codep <- file.path(v7_analysis, "scd_codependency_qc")
rank <- read.csv(file.path(codep, "scd_codependency_qc_adjusted.csv"),
                 check.names = FALSE)
gsea <- read.csv(file.path(codep, "scd_hallmark_raw_vs_qc.csv"),
                 check.names = FALSE)
rank$significant <- rank$FDR_BH < 0.05
rank$label <- ifelse(rank$gene %in% c("VPS72", "CDX1", "WLS", "MPIG6B"),
                     rank$gene, NA)

p3a <- ggplot(rank, aes(raw_spearman_rho, qc_adjusted_rank_correlation)) +
  geom_abline(slope = 1, intercept = 0, colour = light_grey,
              linetype = "dashed", linewidth = 0.4) +
  geom_point(colour = light_blue, alpha = 0.42, size = 0.55) +
  geom_point(data = rank[!is.na(rank$label), ],
             aes(colour = significant), size = 2) +
  geom_text_repel(data = rank[!is.na(rank$label), ],
                  aes(label = label, colour = significant),
                  seed = 20260901, size = 2.6, max.overlaps = Inf,
                  min.segment.length = 0) +
  scale_colour_manual(values = c(`TRUE` = red, `FALSE` = blue),
                      guide = "none") +
  coord_equal() +
  labs(title = "Screen-quality adjustment",
       x = "Raw Spearman rho with SCD",
       y = "QC-adjusted rank correlation") +
  theme_v8_1()

p3b <- ggplot(rank, aes(qc_adjusted_rank_correlation,
                        -log10(qc_adjusted_p))) +
  geom_vline(xintercept = 0, colour = light_grey, linewidth = 0.35) +
  geom_point(aes(colour = significant), alpha = 0.55, size = 0.65) +
  geom_text_repel(data = rank[rank$significant, ], aes(label = gene),
                  colour = red, seed = 20260901, size = 2.8,
                  fontface = "bold", max.overlaps = Inf) +
  scale_colour_manual(values = c(`TRUE` = red, `FALSE` = light_blue),
                      guide = "none") +
  labs(title = "Genome-wide co-dependency",
       subtitle = "VPS72 was the only gene with BH-FDR<0.05",
       x = "QC-adjusted rank correlation with SCD", y = "-log10 P") +
  theme_v8_1()

pathways <- c(
  "HALLMARK_DNA_REPAIR", "HALLMARK_CHOLESTEROL_HOMEOSTASIS",
  "HALLMARK_NOTCH_SIGNALING", "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_HYPOXIA", "HALLMARK_ADIPOGENESIS",
  "HALLMARK_KRAS_SIGNALING_UP", "HALLMARK_APICAL_JUNCTION"
)
gsea_sel <- gsea[gsea$pathway %in% pathways, ]
gsea_sel$pathway_label <- tools::toTitleCase(
  gsub("_", " ", sub("HALLMARK_", "", gsea_sel$pathway))
)
gsea_sel$ranking_label <- ifelse(
  gsea_sel$ranking == "raw_spearman",
  "Raw Spearman ranking", "QC-adjusted ranking"
)
wide_gsea <- reshape(
  gsea_sel[, c("pathway_label", "ranking_label", "NES")],
  idvar = "pathway_label", timevar = "ranking_label", direction = "wide"
)
gsea_sel$pathway_label <- factor(
  gsea_sel$pathway_label,
  levels = rev(tools::toTitleCase(
    gsub("_", " ", sub("HALLMARK_", "", pathways))
  ))
)
p3c <- ggplot() +
  geom_segment(
    data = wide_gsea,
    aes(y = pathway_label, yend = pathway_label,
        x = `NES.Raw Spearman ranking`, xend = `NES.QC-adjusted ranking`),
    colour = light_grey, linewidth = 0.7
  ) +
  geom_vline(xintercept = 0, colour = grey, linetype = "dotted",
             linewidth = 0.45) +
  geom_point(data = gsea_sel,
             aes(NES, pathway_label, colour = ranking_label,
                 shape = ranking_label), size = 2.3) +
  scale_colour_manual(values = c(
    "Raw Spearman ranking" = blue,
    "QC-adjusted ranking" = red
  )) +
  scale_shape_manual(values = c(
    "Raw Spearman ranking" = 16,
    "QC-adjusted ranking" = 15
  )) +
  labs(
    title = "Hallmark enrichment: no pathway passed FDR<0.05",
    subtitle = paste0(
      "DNA repair FDR 0.154 raw / 0.491 adjusted; cholesterol homeostasis ",
      "FDR 0.746 raw / 0.890 adjusted"
    ),
    x = "Normalised enrichment score", y = NULL
  ) +
  theme_v8_1() + theme(legend.position = "bottom")

top3 <- (p3a | p3b) + plot_layout(widths = c(1, 1))
fig3 <- (top3 / patchwork::free(p3c, side = "l")) +
  plot_layout(heights = c(1, 1.05)) +
  plot_annotation(tag_levels = "a")
export_plot(fig3, "Fig3_SCD_codependency_v8_1", 180, 160)
write.csv(rank, file.path(source_dir, "Fig3ab_SCD_codependency_v8_1.csv"),
          row.names = FALSE)
write.csv(gsea_sel, file.path(source_dir, "Fig3c_Hallmark_enrichment_v8_1.csv"),
          row.names = FALSE)

# Figure 4: replace only the adjusted TCGA estimate with the categorical-stage
# refit and keep all frozen display-only Kaplan-Meier curves unchanged.
gse_dir <- file.path(v7_analysis, "gse39582_recalculation")
tcga_dir <- file.path(v7_analysis, "tcga_cbioportal_592")
gse_cox <- read_tsv(file.path(gse_dir, "gse39582_cox_results.tsv"))
tcga_cox <- read_tsv(file.path(tcga_dir, "tcga_cbioportal_cox_results.tsv"))
tcga_corrected <- read_tsv(file.path(correction_root, "tcga_categorical_stage_corrected.tsv"))

forest <- rbind(
  data.frame(label = "GSE39582 RFS, univariable", cohort = "GSE39582",
             model_type = "Univariable", gse_cox[1, c("n", "events", "hr_per_sd", "ci95_low", "ci95_high", "p")]),
  data.frame(label = "GSE39582 RFS, adjusted", cohort = "GSE39582",
             model_type = "Adjusted", gse_cox[2, c("n", "events", "hr_per_sd", "ci95_low", "ci95_high", "p")]),
  data.frame(label = "TCGA OS, univariable", cohort = "TCGA",
             model_type = "Univariable",
             tcga_cox[tcga_cox$model == "OS_univariable",
                      c("n", "events", "hr_per_sd", "ci95_low", "ci95_high", "p")]),
  data.frame(label = "TCGA OS, age +\ncategorical stage adjusted", cohort = "TCGA",
             model_type = "Adjusted",
             tcga_corrected[, c("n", "events", "hr_per_sd", "ci95_low", "ci95_high", "p")])
)
forest$label <- factor(forest$label, levels = rev(forest$label))
forest$stat_label <- sprintf("HR %.3f (%.3f-%.3f)\nP=%.3f",
                             forest$hr_per_sd, forest$ci95_low,
                             forest$ci95_high, forest$p)
p4a <- ggplot(forest, aes(hr_per_sd, label, colour = cohort, shape = model_type)) +
  geom_vline(xintercept = 1, colour = grey, linetype = "dotted", linewidth = 0.45) +
  geom_errorbar(aes(xmin = ci95_low, xmax = ci95_high), orientation = "y",
                width = 0.11, linewidth = 0.7) +
  geom_point(size = 2.4) +
  geom_text(aes(x = 1.285, label = stat_label), hjust = 0,
            colour = black, size = 2.15, lineheight = 0.95) +
  scale_colour_manual(values = c(GSE39582 = blue, TCGA = red)) +
  scale_shape_manual(values = c(Univariable = 16, Adjusted = 15)) +
  coord_cartesian(xlim = c(0.70, 1.82), clip = "off") +
  labs(title = "Frozen eight-gene module",
       subtitle = sprintf("Categorical-stage TCGA refit remains null (P=%.3f)",
                          tcga_corrected$p),
       x = "Hazard ratio per module SD (95% CI)", y = NULL) +
  guides(colour = guide_legend(order = 1, nrow = 1),
         shape = guide_legend(order = 2, nrow = 1)) +
  theme_v8_1() + theme(legend.position = "bottom", legend.box = "vertical")

km_plot <- function(curve, title, stat_row) {
  cox_label <- sprintf("Continuous Cox\nHR %.3f (95%% CI %.3f-%.3f)\nP=%.3f",
                       stat_row$hr_per_sd, stat_row$ci95_low,
                       stat_row$ci95_high, stat_row$p)
  ggplot(curve, aes(time_months, survival, colour = group, fill = group)) +
    geom_ribbon(aes(ymin = ci95_low, ymax = ci95_high), alpha = 0.13, colour = NA) +
    geom_step(linewidth = 0.75) +
    annotate("label", x = max(curve$time_months) * 0.97, y = 0.87,
             label = cox_label, hjust = 1, vjust = 0.5,
             size = 2.05, lineheight = 0.95, colour = grey,
             fill = alpha("white", 0.88), linewidth = 0.18) +
    scale_colour_manual(values = c(Low = blue, High = red)) +
    scale_fill_manual(values = c(Low = blue, High = red)) +
    coord_cartesian(ylim = c(0, 1.02), expand = FALSE) +
    labs(title = title, subtitle = "Median split for display only",
         x = "Months", y = "Survival probability") +
    theme_v8_1() + theme(legend.position = "bottom")
}
gse_km <- read_tsv(file.path(gse_dir, "gse39582_km_curve.tsv"))
tcga_km <- read_tsv(file.path(tcga_dir, "tcga_cbioportal_km_curve.tsv"))
p4b <- km_plot(gse_km, "GSE39582 relapse-free survival", gse_cox[1, ])
p4c <- km_plot(tcga_km, "TCGA overall survival",
               tcga_cox[tcga_cox$model == "OS_univariable", ])
fig4 <- (p4a | (p4b / p4c)) + plot_layout(widths = c(1.18, 1)) +
  plot_annotation(tag_levels = "a")
export_plot(fig4, "Fig4_bulk_survival_v8_1", 180, 158)
write.csv(forest, file.path(source_dir, "Fig4a_survival_forest_v8_1.csv"), row.names = FALSE)

# Figure 6: preserve numbers; repair the cropped panel-d title and increase
# spacing in the lower row.
set.seed(20260901)
g178_dir <- file.path(v7_analysis, "gse178341_recalculation")
scores <- read_tsv(file.path(g178_dir, "gse178341_patient_scores.tsv"))
associations <- read_tsv(file.path(g178_dir, "gse178341_pathway_associations.tsv"))
composition <- read_tsv(file.path(g178_dir, "gse178341_composition_correlations.tsv"))
scores$prediction_raw <- factor(scores$prediction_raw, levels = c("CMS1", "CMS2", "CMS3", "CMS4"))
cms_labels <- c("CMS1\n(n=17)", "CMS2\n(n=23)", "CMS3\n(n=14)", "CMS4\n(n=8)")
cms_plot <- function(column, title, p_value) {
  ggplot(scores, aes(prediction_raw, .data[[column]], colour = prediction_raw == "CMS4")) +
    geom_boxplot(width = 0.58, outlier.shape = NA, fill = "white", linewidth = 0.5) +
    geom_jitter(width = 0.10, size = 1.25, shape = 21, fill = "white", stroke = 0.4) +
    annotate("text", x = 4, y = 2.35,
             label = sprintf("CMS4 vs CMS1-3: P=%.3f", p_value),
             hjust = 1, colour = grey, size = 2.5) +
    scale_colour_manual(values = c(`TRUE` = red, `FALSE` = blue), guide = "none") +
    scale_x_discrete(labels = cms_labels) +
    coord_cartesian(ylim = c(-2.45, 2.55)) +
    labs(title = title, x = NULL, y = "Frozen eight-gene module score") +
    theme_v8_1()
}
p6a <- cms_plot("module_all", "Whole-tumour pseudobulk module", 0.4506366)
p6b <- cms_plot("module_Epi", "Epithelial pseudobulk module", 0.7183457)

primary <- associations[associations$adjustment == "none", ]
primary$label <- paste(ifelse(primary$compartment == "all", "Whole tumour", "Epithelial"),
                       ifelse(grepl("TGFb", primary$analysis), "vs TGF-beta", "vs EMT"))
corr1 <- data.frame(label = primary$label, rho = primary$rho,
                    lo = primary$rho_ci95_low_bootstrap,
                    hi = primary$rho_ci95_high_bootstrap, type = "Pathway")
comp <- composition[composition$variable_2 == "T_cell_fraction_all", ]
corr2 <- data.frame(label = ifelse(comp$variable_1 == "module_all",
                                   "Whole tumour vs T-cell fraction",
                                   "Epithelial vs T-cell fraction"),
                    rho = comp$spearman_rho, lo = NA, hi = NA, type = "Composition")
corr <- rbind(corr1, corr2)
corr$label <- factor(corr$label, levels = rev(corr$label))
p6c <- ggplot(corr, aes(rho, label, colour = grepl("Epithelial", label), shape = type)) +
  geom_vline(xintercept = 0, colour = grey, linetype = "dotted", linewidth = 0.45) +
  geom_errorbar(data = corr[!is.na(corr$lo), ], aes(xmin = lo, xmax = hi),
                orientation = "y", width = 0.10, linewidth = 0.6) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = c(`TRUE` = red, `FALSE` = blue), guide = "none") +
  scale_shape_manual(values = c(Pathway = 16, Composition = 15), guide = "none") +
  scale_y_discrete(labels = function(value) sub(" vs ", "\nvs ", value,
                                                fixed = TRUE)) +
  coord_cartesian(xlim = c(-0.2, 0.92)) +
  labs(title = "Pathway and composition correlations",
       x = "Spearman rho (bootstrap 95% CI for pathways)", y = NULL) +
  theme_v8_1() +
  theme(axis.text.y = element_text(lineheight = 0.88))

mmr_long <- rbind(
  data.frame(MMRStatus = scores$MMRStatus, metric = "T-cell fraction",
             value = scores$T_cell_fraction_all),
  data.frame(MMRStatus = scores$MMRStatus, metric = "Exhausted-label fraction",
             value = scores$exhausted_label_fraction_T)
)
mmr_long$metric <- factor(mmr_long$metric,
                          levels = c("T-cell fraction", "Exhausted-label fraction"))
p6d <- ggplot(mmr_long, aes(metric, value, colour = MMRStatus)) +
  geom_boxplot(position = position_dodge(width = 0.68), width = 0.55,
               outlier.shape = NA, fill = "white", linewidth = 0.5) +
  geom_point(position = position_jitterdodge(jitter.width = 0.09,
                                             dodge.width = 0.68),
             shape = 21, fill = "white", size = 1.05, stroke = 0.35) +
  annotate("text", x = 1, y = 0.96, label = "FDR=0.0139", colour = grey, size = 2.5) +
  annotate("text", x = 2, y = 0.96, label = "FDR=4.0e-6", colour = grey, size = 2.5) +
  scale_colour_manual(values = c(MMRp = grey, MMRd = purple),
                      breaks = c("MMRp", "MMRd"),
                      labels = c(MMRp = "MMRp (n=28)", MMRd = "MMRd (n=34)")) +
  coord_cartesian(ylim = c(-0.03, 1.02)) +
  labs(title = "MMR-associated immune context", x = NULL, y = "Fraction") +
  theme_v8_1() +
  theme(legend.position = "top", plot.title = element_text(margin = margin(b = 5)))

top6 <- (p6a | p6b) + plot_layout(widths = c(1, 1))
bottom6 <- (patchwork::free(p6c, side = "l") | p6d) +
  plot_layout(widths = c(1, 1))
fig6 <- (top6 / bottom6) +
  plot_layout(heights = c(1, 1.08)) +
  plot_annotation(tag_levels = "a")
export_plot(fig6, "Fig6_GSE178341_patient_level_v8_1", 180, 162)

# Figure 7: correct the empirical two-sided definition in the source table,
# explicitly label the primary P value as one-sided, repair annotation overlap,
# and use an en dash in the SCD-VPS72 title.
mn_dir <- file.path(source_root, "analysis", "matched_null_gse178341")
null <- read_tsv(gzfile(file.path(mn_dir, "gse178341_matched_null_rhos.tsv.gz")))
null <- null[null$k == 50, ]
conditioned <- read_tsv(gzfile(file.path(mn_dir, "gse178341_balance_conditioned_modules.tsv.gz")))
mn_summary <- read_tsv(file.path(mn_dir, "gse178341_matched_null_summary.tsv"))
mn_summary <- mn_summary[mn_summary$k == 50, ]
mn_corrected <- read_tsv(file.path(correction_root, "gse178341_matched_null_pvalue_correction.tsv"))
mn_corrected <- mn_corrected[mn_corrected$k == 50, ]
conditioned_summary <- read_tsv(file.path(mn_dir, "gse178341_balance_conditioned_sensitivity.tsv"))
observed <- mn_summary$observed_rho

p7a <- ggplot(null, aes(rho)) +
  geom_density(aes(fill = "Prespecified k=50 null"), alpha = 0.42,
               colour = blue, linewidth = 0.55) +
  geom_density(data = conditioned, aes(rho, colour = "Closest 10% sensitivity"),
               inherit.aes = FALSE, linewidth = 0.7, adjust = 1.05) +
  geom_vline(xintercept = observed, colour = red, linewidth = 0.8) +
  annotate("label", x = min(null$rho) + 0.010, y = 0.72,
           label = sprintf("frozen\nrho=%.3f", observed),
           colour = red, hjust = 0, vjust = 0.5, size = 2.15,
           lineheight = 0.92, fontface = "bold", linewidth = 0,
           label.padding = unit(0.65, "mm"),
           fill = alpha("white", 0.92)) +
  annotate("label", x = min(null$rho) + 0.010, y = Inf,
           label = sprintf(paste0("full-null median %.3f\n",
                                  "one-sided empirical P=%.3f\n",
                                  "equal-tail two-sided P=%.3f\n",
                                  "closest-10%% sensitivity P=%.3f"),
                           mn_summary$null_median,
                           mn_corrected$empirical_p_one_sided_upper,
                           mn_corrected$empirical_p_two_sided_equal_tail,
                           conditioned_summary$empirical_p_one_sided),
           hjust = 0, vjust = 1.04, size = 1.95, lineheight = 0.94,
           colour = grey, linewidth = 0,
           label.padding = unit(0.65, "mm"),
           fill = alpha("white", 0.92)) +
  scale_fill_manual(values = c("Prespecified k=50 null" = light_blue)) +
  scale_colour_manual(values = c("Closest 10% sensitivity" = orange)) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.14))) +
  labs(title = "Matched-null calibration in GSE178341",
       x = "Spearman rho with patient T-cell fraction", y = "Density") +
  theme_v8_1() + theme(legend.position = "bottom")

balance_features <- read_tsv(gzfile(file.path(mn_dir, "gse178341_matching_balance_module_features.tsv.gz")))
balance_summary <- read_tsv(file.path(mn_dir, "gse178341_matching_balance_module_means.tsv"))
feature_labels <- c(mean_expression_log_cp10k = "Mean\nexpression",
                    cell_detection_rate = "Cell\ndetection",
                    tnkilc_vs_other_log2_cpm_ratio = "TNK/ILC\nspecificity")
balance_long <- do.call(rbind, lapply(names(feature_labels), function(feature) {
  values <- balance_features[[feature]]
  data.frame(feature = feature_labels[[feature]],
             standardized_null_value = (values - mean(values)) / sd(values))
}))
balance_points <- data.frame(
  feature = unname(feature_labels[balance_summary$feature]),
  standardized_null_value = balance_summary$standardized_difference_vs_null_mean,
  percentile = balance_summary$observed_percentile
)
balance_long$feature <- factor(balance_long$feature, levels = unname(feature_labels))
balance_points$feature <- factor(balance_points$feature, levels = unname(feature_labels))
balance_points$label_x <- c(1.23, 2.23, 2.77)
balance_points$label_hjust <- c(0, 0, 1)
p7b <- ggplot(balance_long, aes(feature, standardized_null_value)) +
  geom_violin(fill = light_blue, colour = blue, alpha = 0.42,
              linewidth = 0.45, width = 0.62, trim = TRUE) +
  geom_hline(yintercept = 0, colour = light_grey, linewidth = 0.4) +
  geom_point(data = balance_points, aes(feature, standardized_null_value),
             inherit.aes = FALSE, shape = 23, size = 2.8, stroke = 0.55,
             fill = red, colour = red) +
  geom_label(data = balance_points,
             aes(label_x, standardized_null_value,
                 label = sprintf("%.1f%%", percentile),
                 hjust = label_hjust),
             inherit.aes = FALSE, colour = red, fill = alpha("white", 0.92),
             size = 2.15, fontface = "bold", linewidth = 0,
             label.padding = unit(0.45, "mm")) +
  coord_cartesian(ylim = c(-3.2, 3.15)) +
  labs(title = "Module-mean matching balance", x = NULL,
       y = "Standardized value relative to matched null") +
  theme_v8_1()

dm_dir <- file.path(source_root, "analysis", "depmap_vps72_replication")
dm <- read_tsv(file.path(dm_dir, "depmap_scd_vps72_results.tsv"))
dm <- dm[grepl("DepMap", dm$platform), ]
dm$label <- factor(c("Raw Spearman", "QC/common-essentiality adjusted"),
                   levels = rev(c("Raw Spearman", "QC/common-essentiality adjusted")))
p7c <- ggplot(dm, aes(rho, label)) +
  geom_vline(xintercept = 0, colour = grey, linetype = "dotted", linewidth = 0.5) +
  geom_vline(xintercept = -0.492595177156636, colour = light_grey,
             linetype = "dashed", linewidth = 0.55) +
  geom_errorbar(aes(xmin = ci95_low, xmax = ci95_high, colour = label),
                orientation = "y", width = 0.12, linewidth = 0.65) +
  geom_point(aes(colour = label, shape = label), size = 2.7) +
  annotate("text", x = -0.492595177156636, y = 2.36,
           label = "organoid discovery estimate\n(different platform)",
           colour = grey, hjust = 0, vjust = 0.5, size = 2.1) +
  geom_text(aes(x = 0.22,
                label = ifelse(grepl("partial", model),
                               sprintf("rho=%.3f; permutation P=%.3f", rho, p),
                               sprintf("rho=%.3f; P=%.3f", rho, p))),
            hjust = 0, colour = grey, size = 2.2) +
  scale_colour_manual(values = c("Raw Spearman" = blue,
                                 "QC/common-essentiality adjusted" = red),
                      guide = "none") +
  scale_shape_manual(values = c("Raw Spearman" = 16,
                                "QC/common-essentiality adjusted" = 15),
                     guide = "none") +
  coord_cartesian(xlim = c(-0.58, 0.48), clip = "off") +
  labs(title = paste0("SCD", "\u2013", "VPS72 does not replicate in DepMap 2D CRC"),
       subtitle = "n=63 cell lines; 95% bootstrap intervals",
       x = "Correlation of gene-effect profiles", y = NULL) +
  theme_v8_1()

sc_dir <- file.path(source_root, "analysis", "gse132465_replication")
sc <- read_tsv(file.path(sc_dir, "gse132465_patient_scores.tsv"))
sc_summary <- read_tsv(file.path(sc_dir, "gse132465_composition_summary.tsv"))
sc_long <- rbind(
  data.frame(patient = sc$patient, t_cell_fraction = sc$t_cell_fraction,
             score = sc$whole_tumour_module_score, compartment = "Tumour-wide pseudobulk"),
  data.frame(patient = sc$patient, t_cell_fraction = sc$t_cell_fraction,
             score = sc$epithelial_module_score, compartment = "Epithelial pseudobulk")
)
sc_long$compartment <- factor(sc_long$compartment,
                              levels = c("Tumour-wide pseudobulk", "Epithelial pseudobulk"))
sc_annotations <- data.frame(
  compartment = factor(c("Tumour-wide pseudobulk", "Epithelial pseudobulk"),
                       levels = levels(sc_long$compartment)),
  x = rep(min(sc$t_cell_fraction) + 0.012, 2),
  y = rep(max(sc_long$score) + 0.14, 2),
  label = c(sprintf("rho=%.3f\npermutation P=%.4f", sc_summary$rho_whole,
                    sc_summary$permutation_p_whole),
            sprintf("rho=%.3f\npermutation P=%.3f", sc_summary$rho_epithelial,
                    sc_summary$permutation_p_epithelial))
)
p7d <- ggplot(sc_long, aes(t_cell_fraction, score, colour = compartment)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.65, alpha = 0.8) +
  geom_point(size = 1.65, alpha = 0.88) +
  geom_label(data = sc_annotations, aes(x, y, label = label), inherit.aes = FALSE,
             hjust = 0, vjust = 1, colour = black, fill = alpha("white", 0.92),
             size = 2.15, lineheight = 0.94, linewidth = 0,
             label.padding = unit(0.55, "mm")) +
  facet_wrap(~compartment, nrow = 1) +
  scale_colour_manual(values = c("Tumour-wide pseudobulk" = blue,
                                 "Epithelial pseudobulk" = grey), guide = "none") +
  scale_x_continuous(labels = percent_format(accuracy = 1)) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.12))) +
  labs(title = "Independent full-cell CRC replication (GSE132465)",
       subtitle = sprintf("n=23 patients; delta |rho|=%.3f (bootstrap 95%% CI %.3f to %.3f)",
                          sc_summary$attenuation_abs_rho,
                          sc_summary$attenuation_bootstrap_ci_low,
                          sc_summary$attenuation_bootstrap_ci_high),
       x = "Author-annotated T-cell fraction", y = "Frozen module score") +
  theme_v8_1()

fig7 <- ((p7a | p7b) / patchwork::free(p7c, side = "l") / p7d) +
  plot_layout(heights = c(1.0, 0.72, 1.05)) +
  plot_annotation(tag_levels = "a")
export_plot(fig7, "Fig7_composition_controls_and_vps72_replication_v8_1", 180, 210)

write.csv(null, file.path(source_dir, "Fig7a_matched_null_k50_v8_1.csv"), row.names = FALSE)
write.csv(conditioned, file.path(source_dir, "Fig7a_closest10_sensitivity_v8_1.csv"), row.names = FALSE)
write.csv(balance_long, file.path(source_dir, "Fig7b_matching_balance_null_v8_1.csv"), row.names = FALSE)
write.csv(balance_points, file.path(source_dir, "Fig7b_matching_balance_frozen_v8_1.csv"), row.names = FALSE)
write.csv(dm, file.path(source_dir, "Fig7c_depmap_scd_vps72_v8_1.csv"), row.names = FALSE)
write.csv(sc_long, file.path(source_dir, "Fig7d_gse132465_patient_scores_v8_1.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(outdir, "V8_1_changed_figures_R_sessionInfo.txt"))

cat("Saved corrected Figure 3, Figure 4, Figure 6 and Figure 7 sets to",
    outdir, "\n")
