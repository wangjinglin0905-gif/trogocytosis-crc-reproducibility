args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: make_v7_figures.R <root> <outdir>")

root <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
outdir <- args[[2]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(scales))

blue <- "#3C5488"
light_blue <- "#91B7D6"
red <- "#E64B35"
orange <- "#F39B7F"
green <- "#00A087"
purple <- "#7E57C2"
grey <- "#7A7A7A"
light_grey <- "#D9D9D9"
black <- "#222222"

theme_v7 <- function(base_size = 8.5) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      plot.title = element_text(face = "bold", size = base_size + 1.2, hjust = 0,
                                margin = margin(l = 9, b = 2)),
      plot.subtitle = element_text(size = base_size - 0.4, colour = grey, hjust = 0),
      plot.caption = element_text(size = base_size - 1.2, colour = grey, hjust = 0),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.5, colour = black),
      legend.title = element_blank(),
      legend.key.height = unit(3.5, "mm"),
      plot.tag = element_text(face = "bold", size = base_size + 2),
      plot.tag.position = c(0, 1),
      plot.margin = margin(5, 7, 5, 7)
    )
}

export_plot <- function(plot, stem, width_mm, height_mm) {
  png_path <- file.path(outdir, paste0(stem, ".png"))
  tiff_path <- file.path(outdir, paste0(stem, ".tiff"))
  pdf_path <- file.path(outdir, paste0(stem, ".pdf"))
  svg_path <- file.path(outdir, paste0(stem, ".svg"))
  ggsave(png_path, plot, width = width_mm, height = height_mm, units = "mm", dpi = 300,
         device = ragg::agg_png, background = "white")
  ggsave(tiff_path, plot, width = width_mm, height = height_mm, units = "mm", dpi = 600,
         device = "tiff", compression = "lzw", bg = "white")
  ggsave(pdf_path, plot, width = width_mm / 25.4, height = height_mm / 25.4,
         units = "in", device = cairo_pdf, bg = "white")
  svg(svg_path, width = width_mm / 25.4, height = height_mm / 25.4, family = "Arial",
      onefile = TRUE, bg = "white")
  print(plot)
  dev.off()
}

read_tsv <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)

# Figure 1 ------------------------------------------------------------------
org <- file.path(root, "analysis", "organoid")
summary_org <- read.csv(file.path(org, "organoid_panel_dependency_summary.csv"), check.names = FALSE)
screen <- read.csv(gzfile(file.path(org, "organoid_crc_best_library_screen.csv.gz")), check.names = FALSE)
ann <- read.csv(file.path(org, "organoid_model_annotations.csv"), check.names = FALSE)
screen <- merge(screen, ann[, c("sample_ID", "msStatus", "CMS_prediction", "KRAS_class")],
                by = "sample_ID", all.x = TRUE, sort = FALSE)
msi <- read.csv(file.path(org, "organoid_MSI_full_panel_tests.csv"), check.names = FALSE)

summary_org <- summary_org[order(summary_org$median_LFC), ]
summary_org$gene_factor <- factor(summary_org$gene, levels = summary_org$gene)
summary_org$highlight <- ifelse(summary_org$gene %in% c("SCD", "EGFR"), "anchor", "other")
p1a <- ggplot(summary_org, aes(y = gene_factor, x = median_LFC, colour = highlight)) +
  geom_segment(aes(x = 0, xend = median_LFC, yend = gene_factor), linewidth = 0.7) +
  geom_point(size = 1.8) +
  geom_vline(xintercept = 0, colour = light_grey, linewidth = 0.35) +
  scale_colour_manual(values = c(anchor = red, other = blue), guide = "none") +
  labs(title = "Screen-covered panel genes", x = "Median knockout log fold-change", y = NULL) +
  theme_v7(7.4)

scd_cms <- screen[screen$gene == "SCD" & !is.na(screen$CMS_prediction) & screen$CMS_prediction != "", ]
scd_cms$CMS_prediction <- factor(scd_cms$CMS_prediction, levels = c("CMS1", "CMS2", "CMS3", "CMS4"))
cms_counts <- as.data.frame(table(scd_cms$CMS_prediction))
p1b <- ggplot(scd_cms, aes(CMS_prediction, LFC)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, colour = blue, fill = "white", linewidth = 0.45) +
  geom_jitter(width = 0.10, size = 1.25, shape = 21, fill = "white", colour = blue, stroke = 0.4) +
  geom_text(data = cms_counts, aes(x = Var1, y = -1.78, label = paste0("n=", Freq)),
            inherit.aes = FALSE, colour = grey, size = 2.4) +
  annotate("text", x = 4, y = 0.30, label = "Kruskal-Wallis P=0.376", hjust = 1, colour = grey, size = 2.6) +
  coord_cartesian(ylim = c(-1.82, 0.38)) +
  labs(title = "SCD dependency by biobank CMS label", x = NULL, y = "SCD knockout log fold-change") +
  theme_v7()

msi_plot <- screen[screen$gene %in% c("SCD", "EGFR") & screen$msStatus %in% c("MSS", "MSI"), ]
msi_plot$gene <- factor(msi_plot$gene, levels = c("SCD", "EGFR"))
msi_plot$msStatus <- factor(msi_plot$msStatus, levels = c("MSS", "MSI"))
msi_labels <- msi[msi$gene %in% c("SCD", "EGFR"), ]
msi_labels$gene <- factor(msi_labels$gene, levels = c("SCD", "EGFR"))
msi_labels$label <- sprintf("P=%.3f; q=%.3f", msi_labels$p_two_sided, msi_labels$FDR_BH_across_screen_covered_panel)
p1c <- ggplot(msi_plot, aes(msStatus, LFC, colour = msStatus)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, fill = "white", linewidth = 0.45) +
  geom_jitter(width = 0.10, size = 0.95, shape = 21, fill = "white", stroke = 0.35) +
  geom_text(data = msi_labels, aes(x = 1.5, y = 0.34, label = label), inherit.aes = FALSE,
            colour = grey, size = 2.35) +
  facet_wrap(~gene, nrow = 1) +
  scale_colour_manual(values = c(MSS = blue, MSI = red), guide = "none") +
  coord_cartesian(ylim = c(-2.32, 0.45)) +
  labs(title = "Microsatellite-status sensitivity", x = NULL, y = "Knockout log fold-change") +
  theme_v7()

fig1 <- (p1a | (p1b / p1c)) + plot_layout(widths = c(1.08, 1), heights = c(1, 1)) +
  plot_annotation(tag_levels = "a", theme = theme(plot.tag = element_text(face = "bold", size = 11, family = "Arial")))
export_plot(fig1, "Fig1_organoid_dependencies_v7", 180, 198)

# Figure 2 ------------------------------------------------------------------
dep <- file.path(root, "analysis", "depmap_26Q1")
cross <- read.csv(file.path(dep, "depmap_organoid_cross_platform_panel.csv"), check.names = FALSE)
anchors <- read.csv(file.path(dep, "depmap_organoid_anchor_comparison.csv"), check.names = FALSE)
context <- read.csv(file.path(dep, "depmap_target_context_summary.csv"), check.names = FALSE)
cross$highlight <- ifelse(cross$gene %in% c("SCD", "EGFR"), "anchor",
                          ifelse(cross$gene == "CD274", "panel-sensitive", "other"))
label_cross <- cross[cross$highlight != "other", ]
p2a <- ggplot(cross, aes(median_organoid, median_crc_2d)) +
  geom_hline(yintercept = 0, colour = light_grey, linewidth = 0.35) +
  geom_vline(xintercept = 0, colour = light_grey, linewidth = 0.35) +
  geom_point(aes(fill = highlight, colour = highlight, size = highlight), shape = 21, stroke = 0.45) +
  geom_text_repel(data = label_cross, aes(label = gene, colour = highlight), size = 2.7,
                  seed = 20260901, box.padding = 0.28, max.overlaps = Inf, min.segment.length = 0) +
  scale_fill_manual(values = c(anchor = red, `panel-sensitive` = green, other = "white"), guide = "none") +
  scale_colour_manual(values = c(anchor = red, `panel-sensitive` = green, other = blue), guide = "none") +
  scale_size_manual(values = c(anchor = 2.8, `panel-sensitive` = 2.8, other = 1.7), guide = "none") +
  labs(title = "Rank-based cross-platform comparison",
       subtitle = "34 genes: rho=0.351, P=0.0418; sensitivity definition excluding CD274: 33 genes, rho=0.299, P=0.0912",
       x = "CRC organoid median knockout log fold-change",
       y = "DepMap 26Q1 CRC median Chronos gene effect") +
  theme_v7()

anchor_long <- rbind(
  data.frame(gene = anchors$gene, platform = "Organoid official depleted", rate = anchors$pct_officially_depleted),
  data.frame(gene = anchors$gene, platform = "DepMap gene effect < -1", rate = anchors$pct_crc_lt_minus_1)
)
anchor_long <- anchor_long[anchor_long$gene %in% c("SCD", "EGFR"), ]
anchor_long$gene <- factor(anchor_long$gene, levels = c("SCD", "EGFR"))
p2b <- ggplot(anchor_long, aes(gene, rate, fill = platform)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.64) +
  geom_text(aes(label = sprintf("%.1f%%", rate)), position = position_dodge(width = 0.72),
            vjust = -0.35, size = 2.5) +
  scale_fill_manual(values = c("Organoid official depleted" = blue, "DepMap gene effect < -1" = orange)) +
  coord_cartesian(ylim = c(0, 100)) +
  labs(title = "Gene-level anchor replication", x = NULL, y = "Models meeting platform criterion (%)") +
  theme_v7() + theme(legend.position = "bottom")

context2 <- context[context$gene %in% c("SCD", "EGFR"), ]
context2$gene <- factor(context2$gene, levels = c("EGFR", "SCD"))
context2$label <- sprintf("P=%.2g", context2$mannwhitney_p_crc_vs_other)
p2c <- ggplot(context2, aes(median_difference_crc_minus_other, gene, colour = gene)) +
  geom_vline(xintercept = 0, colour = grey, linetype = "dotted", linewidth = 0.45) +
  geom_errorbar(aes(xmin = median_difference_ci95_lo, xmax = median_difference_ci95_hi),
                orientation = "y", width = 0.12, linewidth = 0.65) +
  geom_point(size = 2.5) +
  geom_text(aes(x = 0.04, label = label), hjust = 1, colour = grey, size = 2.4) +
  scale_colour_manual(values = c(EGFR = blue, SCD = red), guide = "none") +
  coord_cartesian(xlim = c(-0.85, 0.06)) +
  labs(title = "CRC-specific dependency shift", x = "Median gene-effect difference\n(CRC minus other tumour lines)", y = NULL) +
  theme_v7()

fig2 <- (p2a / (p2b | p2c)) + plot_layout(heights = c(1.08, 0.92)) +
  plot_annotation(tag_levels = "a", theme = theme(plot.tag = element_text(face = "bold", size = 11, family = "Arial")))
export_plot(fig2, "Fig2_DepMap_validation_v7", 180, 166)

# Figure 3 ------------------------------------------------------------------
codep <- file.path(root, "analysis", "scd_codependency_qc")
rank <- read.csv(file.path(codep, "scd_codependency_qc_adjusted.csv"), check.names = FALSE)
gsea <- read.csv(file.path(codep, "scd_hallmark_raw_vs_qc.csv"), check.names = FALSE)
rank$significant <- rank$FDR_BH < 0.05
rank$label <- ifelse(rank$gene %in% c("VPS72", "CDX1", "WLS", "MPIG6B"), rank$gene, NA)
p3a <- ggplot(rank, aes(raw_spearman_rho, qc_adjusted_rank_correlation)) +
  geom_abline(slope = 1, intercept = 0, colour = light_grey, linetype = "dashed", linewidth = 0.4) +
  geom_point(colour = light_blue, alpha = 0.42, size = 0.55) +
  geom_point(data = rank[!is.na(rank$label), ], aes(colour = significant), size = 2) +
  geom_text_repel(data = rank[!is.na(rank$label), ], aes(label = label, colour = significant),
                  seed = 20260901, size = 2.6, max.overlaps = Inf, min.segment.length = 0) +
  scale_colour_manual(values = c(`TRUE` = red, `FALSE` = blue), guide = "none") +
  coord_equal() +
  labs(title = "Effect of screen-quality adjustment", x = "Raw Spearman rho with SCD", y = "QC-adjusted rank correlation") +
  theme_v7()

p3b <- ggplot(rank, aes(qc_adjusted_rank_correlation, -log10(qc_adjusted_p))) +
  geom_vline(xintercept = 0, colour = light_grey, linewidth = 0.35) +
  geom_point(aes(colour = significant), alpha = 0.55, size = 0.65) +
  geom_text_repel(data = rank[rank$significant, ], aes(label = gene), colour = red,
                  seed = 20260901, size = 2.8, fontface = "bold", max.overlaps = Inf) +
  scale_colour_manual(values = c(`TRUE` = red, `FALSE` = light_blue), guide = "none") +
  labs(title = "Genome-wide co-dependency", subtitle = "VPS72 was the only gene with BH-FDR<0.05",
       x = "QC-adjusted rank correlation with SCD", y = "-log10 P") +
  theme_v7()

pathways <- c("HALLMARK_DNA_REPAIR", "HALLMARK_CHOLESTEROL_HOMEOSTASIS", "HALLMARK_NOTCH_SIGNALING",
              "HALLMARK_G2M_CHECKPOINT", "HALLMARK_HYPOXIA", "HALLMARK_ADIPOGENESIS",
              "HALLMARK_KRAS_SIGNALING_UP", "HALLMARK_APICAL_JUNCTION")
gsea_sel <- gsea[gsea$pathway %in% pathways, ]
gsea_sel$pathway_label <- tools::toTitleCase(gsub("_", " ", sub("HALLMARK_", "", gsea_sel$pathway)))
gsea_sel$ranking_label <- ifelse(gsea_sel$ranking == "raw_spearman", "Raw Spearman ranking", "QC-adjusted ranking")
wide_gsea <- reshape(gsea_sel[, c("pathway_label", "ranking_label", "NES")], idvar = "pathway_label", timevar = "ranking_label", direction = "wide")
gsea_sel$pathway_label <- factor(gsea_sel$pathway_label, levels = rev(tools::toTitleCase(gsub("_", " ", sub("HALLMARK_", "", pathways)))))
p3c <- ggplot() +
  geom_segment(data = wide_gsea, aes(y = pathway_label, yend = pathway_label,
                                     x = `NES.Raw Spearman ranking`, xend = `NES.QC-adjusted ranking`),
               colour = light_grey, linewidth = 0.7) +
  geom_vline(xintercept = 0, colour = grey, linetype = "dotted", linewidth = 0.45) +
  geom_point(data = gsea_sel, aes(NES, pathway_label, colour = ranking_label, shape = ranking_label), size = 2.3) +
  scale_colour_manual(values = c("Raw Spearman ranking" = blue, "QC-adjusted ranking" = red)) +
  scale_shape_manual(values = c("Raw Spearman ranking" = 16, "QC-adjusted ranking" = 15)) +
  labs(title = "Hallmark enrichment: no pathway passed FDR<0.05",
       subtitle = "DNA repair FDR 0.154 raw / 0.491 adjusted; cholesterol homeostasis FDR 0.746 raw / 0.890 adjusted",
       x = "Normalised enrichment score", y = NULL) +
  theme_v7() + theme(legend.position = "bottom")

fig3 <- ((p3a | p3b) / p3c) + plot_layout(heights = c(1, 1.05)) +
  plot_annotation(tag_levels = "a", theme = theme(plot.tag = element_text(face = "bold", size = 11, family = "Arial")))
export_plot(fig3, "Fig3_SCD_codependency_v7", 180, 160)

# Figure 4 ------------------------------------------------------------------
gse_dir <- file.path(root, "analysis", "gse39582_recalculation")
tcga_dir <- file.path(root, "analysis", "tcga_cbioportal_592")
gse_cox <- read_tsv(file.path(gse_dir, "gse39582_cox_results.tsv"))
tcga_cox <- read_tsv(file.path(tcga_dir, "tcga_cbioportal_cox_results.tsv"))
forest_cols <- c("model", "n", "events", "hr_per_sd", "ci95_low", "ci95_high", "p")
forest_row <- function(df, index, label, cohort) {
  out <- df[index, forest_cols, drop = FALSE]
  out$label <- label
  out$cohort <- cohort
  out
}
forest <- rbind(
  forest_row(gse_cox, 1, "GSE39582 RFS, univariable", "GSE39582"),
  forest_row(gse_cox, 2, "GSE39582 RFS, adjusted", "GSE39582"),
  forest_row(tcga_cox, tcga_cox$model == "OS_univariable", "TCGA OS, univariable", "TCGA"),
  forest_row(tcga_cox, tcga_cox$model == "OS_age_stage_adjusted", "TCGA OS, age + stage adjusted", "TCGA")
)
forest$label <- factor(forest$label, levels = rev(forest$label))
forest$cohort <- ifelse(grepl("GSE39582", forest$label), "GSE39582", "TCGA")
forest$model_type <- ifelse(grepl("univariable", forest$label), "Univariable", "Adjusted")
p4a <- ggplot(forest, aes(hr_per_sd, label, colour = cohort, shape = model_type)) +
  geom_vline(xintercept = 1, colour = grey, linetype = "dotted", linewidth = 0.45) +
  geom_errorbar(aes(xmin = ci95_low, xmax = ci95_high),
                orientation = "y", width = 0.11, linewidth = 0.7) +
  geom_point(size = 2.4) +
  scale_colour_manual(values = c(GSE39582 = blue, TCGA = red)) +
  scale_shape_manual(values = c(Univariable = 16, Adjusted = 15)) +
  coord_cartesian(xlim = c(0.70, 1.49)) +
  labs(title = "Frozen eight-gene module", subtitle = "Proportional-hazards tests: all P>0.05",
       x = "Hazard ratio per module SD (95% CI)", y = NULL) +
  theme_v7() + theme(legend.position = "bottom")

km_plot <- function(curve, title) {
  ggplot(curve, aes(time_months, survival, colour = group, fill = group)) +
    geom_ribbon(aes(ymin = ci95_low, ymax = ci95_high), alpha = 0.13, colour = NA) +
    geom_step(linewidth = 0.75) +
    scale_colour_manual(values = c(Low = blue, High = red)) +
    scale_fill_manual(values = c(Low = blue, High = red)) +
    coord_cartesian(ylim = c(0, 1.02), expand = FALSE) +
    labs(title = title, subtitle = "Median split for display only", x = "Months", y = "Survival probability") +
    theme_v7() + theme(legend.position = "bottom")
}
gse_km <- read_tsv(file.path(gse_dir, "gse39582_km_curve.tsv"))
tcga_km <- read_tsv(file.path(tcga_dir, "tcga_cbioportal_km_curve.tsv"))
p4b <- km_plot(gse_km, "GSE39582 relapse-free survival")
p4c <- km_plot(tcga_km, "TCGA overall survival")
fig4 <- (p4a | (p4b / p4c)) + plot_layout(widths = c(1.12, 1), heights = c(1, 1)) +
  plot_annotation(tag_levels = "a", theme = theme(plot.tag = element_text(face = "bold", size = 11, family = "Arial")))
export_plot(fig4, "Fig4_bulk_survival_v7", 180, 158)

# Figure 5 ------------------------------------------------------------------
leuko <- file.path(root, "analysis", "gse146771")
rates <- read.csv(file.path(leuko, "gse146771_candidate_rates_by_patient_tissue_lineage.csv"), check.names = FALSE)
all_rates <- rates[rates$lineage == "All leukocytes", ]
wide <- reshape(all_rates[, c("Sample", "Tissue", "rate")], idvar = "Sample", timevar = "Tissue", direction = "wide")
paired <- wide[complete.cases(wide[, c("rate.N", "rate.T")]), ]
paired_long <- rbind(data.frame(Sample = paired$Sample, tissue = "Adjacent normal", rate = paired$rate.N),
                     data.frame(Sample = paired$Sample, tissue = "Tumour", rate = paired$rate.T))
paired_long$tissue <- factor(paired_long$tissue, levels = c("Adjacent normal", "Tumour"))
p5a <- ggplot(paired_long, aes(tissue, rate, group = Sample)) +
  geom_line(colour = light_grey, linewidth = 0.55) +
  geom_point(aes(colour = tissue), size = 2) +
  scale_colour_manual(values = c("Adjacent normal" = blue, "Tumour" = red), guide = "none") +
  annotate("text", x = 2, y = max(paired_long$rate) * 1.03,
           label = "8 paired patients\nMean T-N = -7.48e-4\nExact sign-flip P=0.125", hjust = 1, vjust = 1,
           colour = grey, size = 2.5) +
  labs(title = "Patient-paired candidate rates", x = NULL, y = "Candidate rate per leukocyte") +
  theme_v7()

tissue <- read.csv(file.path(leuko, "gse146771_candidate_tissue_descriptive.csv"), check.names = FALSE)
tissue$label <- factor(tissue$tissue, levels = c("T", "N", "P"), labels = c("Tumour", "Adjacent normal", "Peripheral blood"))
tissue$count_label <- sprintf("%d/%s", tissue$candidates, comma(tissue$n_cells))
p5b <- ggplot(tissue, aes(rate * 100, label)) +
  geom_errorbar(aes(xmin = wilson95_lo * 100, xmax = wilson95_hi * 100),
                orientation = "y", width = 0.12, colour = blue, linewidth = 0.65) +
  geom_point(colour = blue, size = 2.3) +
  geom_text(aes(x = 0.235, label = count_label), hjust = 1, colour = grey, size = 2.4) +
  coord_cartesian(xlim = c(-0.01, 0.24)) +
  labs(title = "Tissue distribution", x = "Candidate rate (%) with Wilson 95% CI", y = NULL) +
  theme_v7()

sensitivity <- read.csv(file.path(leuko, "gse146771_threshold_sensitivity.csv"), check.names = FALSE)
sensitivity$ceacam_label <- factor(paste0(">", sensitivity$CEACAM5_TPM_threshold), levels = c(">2", ">1", ">0"))
sensitivity$epi_label <- factor(sprintf("%.1f", 100 * sensitivity$negative_control_epi_quantile), levels = c("99.0", "99.5", "99.9"))
p5c <- ggplot(sensitivity, aes(epi_label, ceacam_label, fill = total_candidates)) +
  geom_tile(colour = "white", linewidth = 0.55) +
  geom_text(aes(label = total_candidates, colour = total_candidates > 30), fontface = "bold", size = 3) +
  scale_fill_gradient(low = "#EFF6FC", high = "#084D96", guide = "none") +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = black), guide = "none") +
  labs(title = "Threshold sensitivity",
       x = "Peripheral-blood epithelial-score percentile", y = "CEACAM5 TPM threshold") +
  theme_v7()

fig5 <- (p5a | (p5b / p5c)) + plot_layout(widths = c(1.08, 1)) +
  plot_annotation(tag_levels = "a", theme = theme(plot.tag = element_text(face = "bold", size = 11, family = "Arial")))
export_plot(fig5, "Fig5_leukocyte_candidates_v7", 180, 158)

# Figure 6 ------------------------------------------------------------------
g178_dir <- file.path(root, "analysis", "gse178341_recalculation")
scores <- read_tsv(file.path(g178_dir, "gse178341_patient_scores.tsv"))
associations <- read_tsv(file.path(g178_dir, "gse178341_pathway_associations.tsv"))
groups <- read_tsv(file.path(g178_dir, "gse178341_group_comparisons.tsv"))
composition <- read_tsv(file.path(g178_dir, "gse178341_composition_correlations.tsv"))
scores$prediction_raw <- factor(scores$prediction_raw, levels = c("CMS1", "CMS2", "CMS3", "CMS4"))
cms_labels <- c("CMS1\n(n=17)", "CMS2\n(n=23)", "CMS3\n(n=14)", "CMS4\n(n=8)")
cms_plot <- function(column, title, p_value) {
  ggplot(scores, aes(prediction_raw, .data[[column]], colour = prediction_raw == "CMS4")) +
    geom_boxplot(width = 0.58, outlier.shape = NA, fill = "white", linewidth = 0.5) +
    geom_jitter(width = 0.10, size = 1.25, shape = 21, fill = "white", stroke = 0.4) +
    annotate("text", x = 4, y = 2.35, label = sprintf("CMS4 vs CMS1-3: P=%.3f", p_value), hjust = 1, colour = grey, size = 2.5) +
    scale_colour_manual(values = c(`TRUE` = red, `FALSE` = blue), guide = "none") +
    scale_x_discrete(labels = cms_labels) +
    coord_cartesian(ylim = c(-2.45, 2.55)) +
    labs(title = title, x = NULL, y = "Frozen eight-gene module score") +
    theme_v7()
}
p6a <- cms_plot("module_all", "Whole-tumour pseudobulk module", 0.4506366)
p6b <- cms_plot("module_Epi", "Epithelial pseudobulk module", 0.7183457)

primary <- associations[associations$adjustment == "none", ]
primary$label <- paste(ifelse(primary$compartment == "all", "Whole tumour", "Epithelial"),
                       ifelse(grepl("TGFb", primary$analysis), "vs TGF-beta", "vs EMT"))
corr1 <- data.frame(label = primary$label, rho = primary$rho, lo = primary$rho_ci95_low_bootstrap,
                    hi = primary$rho_ci95_high_bootstrap, type = "Pathway")
comp <- composition[composition$variable_2 == "T_cell_fraction_all", ]
corr2 <- data.frame(label = ifelse(comp$variable_1 == "module_all", "Whole tumour vs T-cell fraction", "Epithelial vs T-cell fraction"),
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
  coord_cartesian(xlim = c(-0.2, 0.92)) +
  labs(title = "Pathway and composition correlations", x = "Spearman rho (bootstrap 95% CI for pathways)", y = NULL) +
  theme_v7()

mmr_long <- rbind(
  data.frame(MMRStatus = scores$MMRStatus, metric = "T-cell fraction", value = scores$T_cell_fraction_all),
  data.frame(MMRStatus = scores$MMRStatus, metric = "Exhausted-label fraction", value = scores$exhausted_label_fraction_T)
)
mmr_long$metric <- factor(mmr_long$metric, levels = c("T-cell fraction", "Exhausted-label fraction"))
p6d <- ggplot(mmr_long, aes(metric, value, colour = MMRStatus)) +
  geom_boxplot(position = position_dodge(width = 0.68), width = 0.55, outlier.shape = NA, fill = "white", linewidth = 0.5) +
  geom_point(position = position_jitterdodge(jitter.width = 0.09, dodge.width = 0.68), shape = 21, fill = "white", size = 1.05, stroke = 0.35) +
  annotate("text", x = 1, y = 0.96, label = "FDR=0.0139", colour = grey, size = 2.5) +
  annotate("text", x = 2, y = 0.96, label = "FDR=4.0e-6", colour = grey, size = 2.5) +
  scale_colour_manual(values = c(MMRp = grey, MMRd = purple), labels = c("MMRp (n=28)", "MMRd (n=34)")) +
  coord_cartesian(ylim = c(-0.03, 1.02)) +
  labs(title = "MMR-associated immune contexture", x = NULL, y = "Fraction") +
  theme_v7() + theme(legend.position = "top")

fig6 <- ((p6a | p6b) / (p6c | p6d)) +
  plot_annotation(tag_levels = "a", theme = theme(plot.tag = element_text(face = "bold", size = 11, family = "Arial")))
export_plot(fig6, "Fig6_GSE178341_patient_level_v7", 180, 156)

cat("Saved six V7 figure sets to", outdir, "\n")
