args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: make_figures_2_5.R <root> <outdir>")

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
  geom_text(aes(x = -0.03, label = label), hjust = 1, colour = grey, size = 2.4,
            position = position_nudge(y = 0.14)) +
  scale_colour_manual(values = c(EGFR = blue, SCD = red), guide = "none") +
  coord_cartesian(xlim = c(-0.85, 0.06)) +
  labs(title = "CRC-specific dependency shift", x = "Median gene-effect difference\n(CRC minus other tumour lines)", y = NULL) +
  theme_v7()

fig2 <- (p2a / (p2b | p2c)) + plot_layout(heights = c(1.08, 0.92)) +
  plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(face = "bold", size = 11, family = "Arial")))
export_plot(fig2, "Fig2_DepMap_validation_v7", 180, 166)

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
  plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(face = "bold", size = 11, family = "Arial")))
export_plot(fig5, "Fig5_leukocyte_candidates_v7", 180, 158)


cat("Saved current Figures 2 and 5.\n")
