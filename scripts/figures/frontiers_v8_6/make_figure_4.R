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
       subtitle = sprintf("TCGA model with categorical stage (P=%.3f)",
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
  plot_annotation(tag_levels = "A")
export_plot(fig4, "Fig4_bulk_survival_v8_1", 180, 158)
write.csv(forest, file.path(source_dir, "Fig4a_survival_forest_v8_1.csv"), row.names = FALSE)
