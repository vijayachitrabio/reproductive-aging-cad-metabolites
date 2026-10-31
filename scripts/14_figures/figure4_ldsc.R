#!/usr/bin/env Rscript
# Figure 4: LDSC genetic correlations of prioritized metabolites with ANM and CAD.

source(file.path(Sys.getenv("RAC_ROOT", "."), "config", "paths.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

root <- RAC_ROOT
out_dir <- file.path(root, "results/ldsc/figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

d <- fread(file.path(root, "results/ldsc/ldsc_rg_summary_all_pairs.tsv"))
d <- d[pair_type != "ANM_vs_CAD"]
d[, outcome := fifelse(pair_type == "CAD_vs_metabolite", "CAD", "ANM")]
d[, sig := fifelse(fdr_across_all_37 < 0.05, "FDR < 0.05", "not significant")]

# order metabolites by CAD rg magnitude for a stable, readable layout
ord <- d[outcome == "CAD"][order(rg), trait2]
d[, trait2 := factor(trait2, levels = ord)]

p <- ggplot(d, aes(x = rg, y = trait2, color = outcome, alpha = sig)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_errorbarh(aes(xmin = rg - 1.96 * se, xmax = rg + 1.96 * se),
                 height = 0.2, position = position_dodge(width = 0.6)) +
  geom_point(size = 2.6, position = position_dodge(width = 0.6)) +
  scale_color_manual(values = c("CAD" = "#B22222", "ANM" = "#1F6FB2")) +
  scale_alpha_manual(values = c("FDR < 0.05" = 1, "not significant" = 0.35)) +
  labs(
    x = expression(paste("Genome-wide genetic correlation (", r[g], ") with outcome")),
    y = NULL, color = "Outcome", alpha = NULL,
    title = "Genome-wide genetic correlation (LDSC) between 18 prioritized metabolites\nand ANM / CAD"
  ) +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 11), legend.position = "bottom")

out_file <- file.path(out_dir, "ldsc_rg_anm_cad_metabolites.png")
ggsave(out_file, p, width = 7.5, height = 6, dpi = 320)
message("Wrote ", out_file)
