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

# Figure 1: full-range redraw; no observation is removed.
set.seed(20260907)
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
cms_range <- range(scd_cms$LFC)
cms_span <- diff(cms_range)
cms_ticks <- setNames(paste0(cms_counts$Var1, "\n(n=", cms_counts$Freq, ")"), cms_counts$Var1)
p1b <- ggplot(scd_cms, aes(CMS_prediction, LFC)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, colour = blue, fill = "white", linewidth = 0.45) +
  geom_jitter(width = 0.10, height = 0, size = 1.25, shape = 21, fill = "white", colour = blue, stroke = 0.4) +
  annotate("text", x = 4, y = max(cms_range) + 0.12 * cms_span, label = "Kruskal-Wallis P=0.376", hjust = 1, colour = grey, size = 2.6) +
  scale_x_discrete(labels = cms_ticks) +
  coord_cartesian(ylim = c(min(cms_range) - 0.03 * cms_span, max(cms_range) + 0.23 * cms_span)) +
  labs(title = "SCD dependency by biobank CMS label", x = NULL, y = "SCD knockout log fold-change") +
  theme_v7()

msi_plot <- screen[screen$gene %in% c("SCD", "EGFR") & screen$msStatus %in% c("MSS", "MSI"), ]
msi_plot$gene <- factor(msi_plot$gene, levels = c("SCD", "EGFR"))
msi_plot$msStatus <- factor(msi_plot$msStatus, levels = c("MSS", "MSI"))
msi_labels <- msi[msi$gene %in% c("SCD", "EGFR"), ]
msi_labels$gene <- factor(msi_labels$gene, levels = c("SCD", "EGFR"))
msi_labels$label <- sprintf("P=%.3f; q=%.3f", msi_labels$p_two_sided, msi_labels$FDR_BH_across_screen_covered_panel)
msi_range <- range(msi_plot$LFC)
msi_span <- diff(msi_range)
p1c <- ggplot(msi_plot, aes(msStatus, LFC, colour = msStatus)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, fill = "white", linewidth = 0.45) +
  geom_jitter(width = 0.10, height = 0, size = 0.95, shape = 21, fill = "white", stroke = 0.35) +
  geom_text(data = msi_labels, aes(x = 1.5, y = max(msi_range) + 0.10 * msi_span, label = label), inherit.aes = FALSE,
            colour = grey, size = 2.35) +
  facet_wrap(~gene, nrow = 1) +
  scale_colour_manual(values = c(MSS = blue, MSI = red), guide = "none") +
  coord_cartesian(ylim = c(min(msi_range) - 0.03 * msi_span, max(msi_range) + 0.20 * msi_span)) +
  labs(title = "Microsatellite-status sensitivity", x = NULL, y = "Knockout log fold-change") +
  theme_v7()

fig1 <- (p1a | (p1b / p1c)) + plot_layout(widths = c(1.08, 1), heights = c(1, 1)) +
  plot_annotation(tag_levels = "A", theme = theme(plot.tag = element_text(face = "bold", size = 11, family = "Arial")))
export_plot(fig1, "Figure_1", 180, 198)


write.csv(scd_cms, file.path(outdir, "Figure_1B_source_data.csv"), row.names = FALSE)
write.csv(msi_plot, file.path(outdir, "Figure_1C_source_data.csv"), row.names = FALSE)
stopifnot(nrow(scd_cms) == 52, nrow(msi_plot) == 170)
yr_b <- ggplot_build(p1b)$layout$panel_params[[1]]$y.range
yr_c <- ggplot_build(p1c)$layout$panel_params[[1]]$y.range
visibility <- data.frame(
  panel=c("B", "C"), n=c(nrow(scd_cms), nrow(msi_plot)),
  data_min=c(min(scd_cms$LFC), min(msi_plot$LFC)),
  data_max=c(max(scd_cms$LFC), max(msi_plot$LFC)),
  old_visible_min=c(-1.93, -2.4585), old_visible_max=c(0.49, 0.5885),
  old_clipped_n=c(sum(scd_cms$LFC < -1.93 | scd_cms$LFC > 0.49),
                  sum(msi_plot$LFC < -2.4585 | msi_plot$LFC > 0.5885)),
  new_visible_min=c(yr_b[1], yr_c[1]), new_visible_max=c(yr_b[2], yr_c[2]),
  new_clipped_n=c(sum(scd_cms$LFC < yr_b[1] | scd_cms$LFC > yr_b[2]),
                  sum(msi_plot$LFC < yr_c[1] | msi_plot$LFC > yr_c[2])))
stopifnot(all(visibility$new_clipped_n == 0))
write.csv(visibility, file.path(outdir,"Figure_1_visibility_audit.csv"), row.names=FALSE)
retests <- data.frame(test="SCD_CMS_Kruskal_Wallis", p_recomputed=kruskal.test(LFC~CMS_prediction,scd_cms)$p.value, p_displayed=0.376)
for (g in c("SCD", "EGFR")) {
  d <- msi_plot[msi_plot$gene==g,]
  p <- wilcox.test(LFC~msStatus,d,exact=FALSE)$p.value
  stored <- msi[msi$gene==g,"p_two_sided"]
  retests <- rbind(retests,data.frame(test=paste0(g,"_MSI_MSS"),p_recomputed=p,p_displayed=stored))
}
write.csv(retests,file.path(outdir,"Figure_1_subgroup_retests.csv"),row.names=FALSE)
print(visibility)
print(retests)
