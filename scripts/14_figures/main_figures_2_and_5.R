#!/usr/bin/env Rscript
# Figures 2 and 5.
# Figure 2: two-sample MR screens of 56 metabolites against ANM and CAD.
# Figure 5: STRING network drawn from the saved edge list and module table.
# Inputs are the result tables; nothing is re-estimated.

source(file.path(Sys.getenv("RAC_ROOT", "."), "config", "paths.R"))

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(ggrepel); library(patchwork)
  library(igraph); library(ggraph); library(tidygraph)
})
root <- RAC_ROOT
out  <- file.path(root, "results", "figures")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
DPI <- 600
okabe <- c(red = "#D55E00", orange = "#E69F00", blue = "#0072B2", purple = "#CC79A7", grey = "#9A9A9A", teal = "#009E73")

# ------------------------------------------------------------------ Figure 2
anm <- fread(file.path(root, "results/metabolome_to_anm_mr_results.tsv"))[method == "Inverse variance weighted"]
cad <- fread(file.path(root, "results/cad_all_metabolites/metabolome_to_cad_ivw_global_and_targeted_fdr.tsv"))[method == "Inverse variance weighted"]
ev  <- fread(file.path(root, "results/final_integrated_evidence/metabolite_integrated_evidence.tsv"), na.strings = "")
tier <- fread(file.path(root, "results/cad_all_metabolites/metabolic_paper_evidence_table.tsv"))[, .(metabolite, evidence_tier)]
flag <- ev[!is.na(direction_sensitivity_flag), metabolite]
stopifnot(nrow(anm) == 56, nrow(cad) == 56)
lab <- function(x) gsub("_", " ", x)

anm[, cat := fifelse(fdr < 0.05, "FDR < 0.05", "Not significant")]
anm[, cat := factor(cat, levels = c("FDR < 0.05", "Not significant"))]
thrA <- anm[fdr < 0.05, max(pval)]

cad <- merge(cad, tier, by.x = "exposure", by.y = "metabolite", all.x = TRUE)
cad[, cat := fcase(
  exposure %in% flag & shared_anm_cad_global_fdr %in% TRUE, "Shared, direction-sensitive",
  shared_anm_cad_global_fdr %in% TRUE, "Shared: FDR < 0.05 for ANM and CAD",
  evidence_tier == "tier3_nominal_cad_followup", "ANM FDR < 0.05, CAD nominal only",
  cad_global_fdr < 0.05, "CAD FDR < 0.05 only",
  default = "Not significant")]
lv <- c("Shared: FDR < 0.05 for ANM and CAD", "Shared, direction-sensitive", "CAD FDR < 0.05 only", "ANM FDR < 0.05, CAD nominal only", "Not significant")
cad[, cat := factor(cat, levels = lv)]
cat("Figure 2B category counts:\n"); print(cad[, .N, by = cat][order(cat)])
cat("Figure 2A: ANM FDR<0.05 =", sum(anm$fdr < 0.05), "\n")
thrB <- cad[cad_global_fdr < 0.05, max(pval)]

th <- theme_classic(base_size = 9) + theme(legend.position = "bottom", legend.title = element_blank(),
        legend.text = element_text(size = 7.5), plot.title = element_text(face = "bold", size = 9.5),
        plot.tag = element_text(face = "bold", size = 12))
pA <- ggplot(anm, aes(b, -log10(pval), colour = cat)) +
  geom_hline(yintercept = -log10(thrA), linetype = "dashed", colour = "grey40", linewidth = 0.35) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_point(size = 1.9, alpha = 0.9) +
  geom_text_repel(data = anm[fdr < 0.05], aes(label = lab(exposure)), size = 2.35, max.overlaps = Inf,
                  segment.size = 0.2, min.segment.length = 0.1, box.padding = 0.25, show.legend = FALSE) +
  scale_colour_manual(values = c("FDR < 0.05" = okabe[["red"]], "Not significant" = okabe[["grey"]])) +
  labs(title = "Metabolites and age at natural menopause (ANM)",
       x = "IVW estimate: ANM change per 1-SD higher metabolite (units of the ANM GWAS)", y = expression(-log[10](italic(P)))) + th
pB <- ggplot(cad, aes(b, -log10(pval), colour = cat)) +
  geom_hline(yintercept = -log10(thrB), linetype = "dashed", colour = "grey40", linewidth = 0.35) +
  geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_point(size = 1.9, alpha = 0.9) +
  geom_text_repel(data = cad[cat != "Not significant"], aes(label = lab(exposure)), size = 2.35, max.overlaps = Inf,
                  segment.size = 0.2, min.segment.length = 0.1, box.padding = 0.25, show.legend = FALSE) +
  scale_colour_manual(values = setNames(c(okabe[["red"]], okabe[["orange"]], okabe[["blue"]], okabe[["purple"]], okabe[["grey"]]), lv)) +
  guides(colour = guide_legend(nrow = 3, override.aes = list(size = 2.5))) +
  labs(title = "Metabolites and coronary artery disease (CAD)",
       x = "IVW estimate: log odds of CAD per 1-SD higher metabolite", y = expression(-log[10](italic(P)))) + th
fig2 <- (pA / pB) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 12))
ggsave(file.path(out, "Figure_2_MR_screens.png"), fig2, width = 6.69, height = 8.6, dpi = DPI, bg = "white")

# ------------------------------------------------------------------ Figure 5
nd <- file.path(root, "results/ppi_network")
edges <- fread(file.path(nd, "string_network_edges.tsv"))
mem   <- fread(file.path(nd, "network_module_membership.tsv"))
stopifnot(nrow(mem) == 22, nrow(edges) == 43)
mem[, type := fifelse(gene == "FADS1", "FADS1", fifelse(seed %in% TRUE, "Seed gene", "Added STRING interactor"))]
mem[, type := factor(type, levels = c("FADS1", "Seed gene", "Added STRING interactor"))]
fills <- c("FADS1" = "#C7465B", "Seed gene" = "#2E5E8C", "Added STRING interactor" = "#D3D9DF")
txt   <- c("FADS1" = "white", "Seed gene" = "white", "Added STRING interactor" = "#1F2933")
mod_name <- c(`1` = "Fatty-acid biosynthesis / elongation", `2` = "AMH receptor signalling",
              `3` = "Lipoxygenase (ALOX15)", `4` = "Follicular gonadotropin signalling",
              `5` = "Cyclooxygenase", `6` = "Lipoxygenase accessory (ALOX5)")
cross <- merge(merge(edges, mem[, .(gene, mA = module)], by.x = "from", by.y = "gene"), mem[, .(gene, mB = module)], by.x = "to", by.y = "gene")
cat("Inter-module edges in the baseline network:", sum(cross$mA != cross$mB), "\n")
draw_module <- function(m, layout = "kk", pad = 0.18, lab_size = 2.5) {
  nodes <- mem[module == m, .(name = gene, type)]
  ed <- edges[from %in% nodes$name & to %in% nodes$name, .(from, to)]
  g <- tbl_graph(nodes = nodes, edges = ed, directed = FALSE)
  set.seed(7)
  ggraph(g, layout = layout) +
    geom_edge_link(colour = "#9AA5B1", width = 0.45) +
    geom_node_label(aes(label = name, fill = type, colour = type), size = lab_size, label.size = 0,
                    label.padding = unit(pad, "lines"), label.r = unit(0.25, "lines"), fontface = "bold") +
    scale_fill_manual(values = fills, drop = FALSE) + scale_colour_manual(values = txt, drop = FALSE) +
    scale_x_continuous(expand = expansion(mult = 0.18)) + scale_y_continuous(expand = expansion(mult = 0.22)) +
    theme_void(base_size = 8) + theme(legend.position = "none",
      plot.title = element_text(size = 7.5, face = "bold", hjust = 0.5), plot.margin = margin(2, 2, 2, 2)) +
    labs(title = sprintf("%s (module %s; %d proteins, %d edges)", mod_name[[as.character(m)]], m, nrow(nodes), nrow(ed)))
}
pa <- draw_module(1, layout = "kk", pad = 0.22, lab_size = 2.7) +
  labs(title = "Module 1: fatty-acid biosynthesis/elongation\n(10 proteins, 35 edges)") +
  theme(plot.title = element_text(size = 8, face = "bold", hjust = 0.5))
nodesB <- mem[module != 1, .(name = gene, type, module)]
nodesB[, mlab := factor(sprintf("Module %s: %s", module, mod_name[as.character(module)]),
                        levels = sprintf("Module %s: %s", 2:6, mod_name[as.character(2:6)]))]
edB <- edges[from %in% nodesB$name & to %in% nodesB$name, .(from, to)]
gB <- tbl_graph(nodes = nodesB, edges = edB, directed = FALSE)
set.seed(11)
pb <- ggraph(gB, layout = "kk") +
  geom_edge_link(colour = "#9AA5B1", width = 0.45) +
  geom_node_label(aes(label = name, fill = type, colour = type), size = 2.5, label.padding = unit(0.2, "lines"),
                  label.r = unit(0.25, "lines"), fontface = "bold") +
  facet_nodes(~mlab, scales = "free", ncol = 2, labeller = label_wrap_gen(24)) +
  scale_fill_manual(values = fills, drop = FALSE) + scale_colour_manual(values = txt, drop = FALSE) +
  scale_x_continuous(expand = expansion(mult = 0.45)) + scale_y_continuous(expand = expansion(mult = 0.45)) +
  theme_void(base_size = 8) + theme(legend.position = "none", strip.text = element_text(size = 7, face = "bold", margin = margin(2, 0, 2, 0)),
    panel.border = element_rect(colour = "#D5DDE5", fill = NA, linewidth = 0.4), panel.spacing = unit(3, "pt"),
    plot.title = element_text(size = 8, face = "bold", hjust = 0.5)) +
  labs(title = "Modules 2-6 (no edge connects them to module 1)")
leg <- ggplot(data.frame(x = 1:3, y = 1, t = factor(levels(mem$type), levels = levels(mem$type))), aes(x, y, fill = t)) +
  geom_point(shape = 22, size = 4, colour = "white") + scale_fill_manual(values = fills, name = NULL) +
  theme_void(base_size = 8) + theme(legend.position = "bottom", legend.text = element_text(size = 8))
legend_only <- cowplot::get_legend(leg)
fig5 <- (pa | pb) + plot_layout(widths = c(1.1, 1)) + plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(face = "bold", size = 12))
fig5 <- cowplot::plot_grid(fig5, legend_only, ncol = 1, rel_heights = c(1, 0.06))
ggsave(file.path(out, "Figure_5_STRING_network.png"), fig5, width = 6.69, height = 5.2, dpi = DPI, bg = "white")
message("Wrote Figures 2 and 5 to ", out)
