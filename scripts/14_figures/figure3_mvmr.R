#!/usr/bin/env Rscript
# Figure 3: multivariable MR forest plot.

source(file.path(Sys.getenv("RAC_ROOT", "."), "config", "paths.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

root <- RAC_ROOT
in_file <- file.path(root, "results/mvmr_axis/metabolic_axis_mvmr_results.tsv")
out_dir <- file.path(root, "results/mvmr_axis/figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

res <- fread(in_file)
res[, exposure := factor(exposure, levels = rev(c("Total_C", "Total_TG", "ApoA1", "GlycA")))]
res[, sig := fifelse(fdr_within_outcome < 0.05, "FDR < 0.05", "Not FDR < 0.05")]
res[, lo := beta - 1.96 * se]
res[, hi := beta + 1.96 * se]

p <- ggplot(res, aes(x = beta, y = exposure, color = sig)) +
  geom_vline(xintercept = 0, linewidth = 0.35, color = "grey55") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0, linewidth = 0.7) +
  geom_point(size = 2.6) +
  facet_wrap(~ outcome, scales = "free_x") +
  scale_color_manual(values = c("FDR < 0.05" = "#1F7A6D", "Not FDR < 0.05" = "#707070")) +
  labs(
    x = "Multivariable MR estimate (beta per SD metabolite)",
    y = NULL,
    color = NULL,
    title = "Core Metabolic Axis MVMR",
    subtitle = "Mutual adjustment for cholesterol, triglycerides, ApoA1, and GlycA"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "grey94", color = "grey80"),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold")
  )

out_file <- file.path(out_dir, "metabolic_axis_mvmr_forest.png")
ggsave(out_file, p, width = 7.2, height = 4.2, dpi = 600)
message("Wrote ", out_file)
