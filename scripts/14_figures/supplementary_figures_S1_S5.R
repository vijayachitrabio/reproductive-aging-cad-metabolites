#!/usr/bin/env Rscript
# Supplementary Figures S1-S5.
#   - S1: reverse MR, ANM and CAD -> prioritized metabolites
#   - S2: evidence matrix for the 18 prioritized metabolites, grouped by Tier 1/2/3
#   - S3: metabolite-CAD colocalization
#   - S4: metabolite-ovarian trait colocalization
#   - S5: metabolite -> AMH MR

source(file.path(Sys.getenv("RAC_ROOT", "."), "config", "paths.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(ggnewscale)
})

root <- RAC_ROOT
out_dir <- file.path(root, "results", "figures", "supplementary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
library(patchwork)
DPI <- 500

tier_of <- function(x) fifelse(grepl("^tier1", x), "Tier 1", fifelse(grepl("^tier2", x), "Tier 2", "Tier 3"))

# ---- S1: reverse MR, ANM and CAD -> prioritized metabolites -----------
rev_panel <- function(f, expo, ttl) {
  res <- fread(f)
  res[, `:=`(lo = reverse_beta - 1.96 * reverse_se, hi = reverse_beta + 1.96 * reverse_se)]
  res[, status := factor(fifelse(reverse_fdr < 0.05, "FDR < 0.05", fifelse(reverse_p < 0.05, "Nominal P < 0.05", "Not significant")),
                         levels = c("FDR < 0.05", "Nominal P < 0.05", "Not significant"))]
  res[, metabolite := factor(gsub("_", " ", metabolite), levels = rev(gsub("_", " ", res[order(reverse_beta)]$metabolite)))]
  ggplot(res, aes(reverse_beta, metabolite, color = status)) +
    geom_vline(xintercept = 0, linewidth = 0.35, color = "grey55") +
    geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0, linewidth = 0.55) +
    geom_point(size = 2) +
    scale_color_manual(values = c("FDR < 0.05" = "#1F7A6D", "Nominal P < 0.05" = "#B7791F", "Not significant" = "#707070"), drop = FALSE) +
    labs(x = "IVW estimate (95% CI), metabolite SD", y = NULL, color = NULL, title = ttl) +
    theme_classic(base_size = 9.5) + theme(legend.position = "bottom", plot.title = element_text(face = "bold", size = 10))
}
pA <- rev_panel(file.path(root, "results/reverse_mr_anm_to_metabolites/anm_to_prioritized_metabolites_reverse_summary.tsv"), "ANM",
                "Exposure: ANM (245 instruments)")
pB <- rev_panel(file.path(root, "results/reverse_mr_cad_to_metabolites/cad_to_prioritized_metabolites_reverse_summary.tsv"), "CAD",
                "Exposure: CAD (39 instruments)")
p1 <- (pA | pB) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") & theme(legend.position = "bottom", plot.tag = element_text(face = "bold"))
ggsave(file.path(out_dir, "Supplementary_Figure_S1_reverse_MR.png"), p1, width = 10, height = 5.4, dpi = DPI, bg = "white")

# ---- S5: prioritized metabolites -> AMH --------------------------------------------------
amh <- fread(file.path(root, "results/metabolome_to_amh/prioritized_metabolite_to_amh_summary.tsv"))
amh[, `:=`(lo = amh_beta - 1.96 * amh_se, hi = amh_beta + 1.96 * amh_se)]
amh[, status := fifelse(amh_fdr_all_metabolites < 0.05, "FDR < 0.05", "Not FDR < 0.05")]
amh[, tier_label := tier_of(evidence_tier)]
amh[, metabolite := factor(metabolite, levels = rev(amh[order(amh_beta)]$metabolite))]
p5 <- ggplot(amh, aes(amh_beta, metabolite, color = status, shape = tier_label)) +
  geom_vline(xintercept = 0, linewidth = 0.35, color = "grey55") +
  geom_errorbar(aes(xmin = lo, xmax = hi), orientation = "y", width = 0, linewidth = 0.55) +
  geom_point(size = 2.5) +
  scale_color_manual(values = c("FDR < 0.05" = "#1F7A6D", "Not FDR < 0.05" = "#707070")) +
  labs(x = "MR estimate (IVW, 95% CI): genetically predicted metabolite -> AMH",
       y = NULL, color = NULL, shape = NULL,
       title = "Secondary MR of the 18 prioritized metabolites on AMH",
       subtitle = "FDR computed across all 56 metabolites; marker shape shows the evidence tier") +
  theme_classic(base_size = 10.5) +
  theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
ggsave(file.path(out_dir, "Supplementary_Figure_S5_AMH_MR.png"), p5, width = 7.4, height = 5.4, dpi = DPI, bg = "white")

# ---- S4: colocalization of metabolite loci with ANM and AMH -----------------------------
co <- fread(file.path(root, "results/coloc_metabolite_ovarian/metabolite_ovarian_coloc_pilot_summary.tsv"))
co <- co[status == "ok" & PP.H4 >= 0.01]
co[, lead_pos := as.integer(sub(".*_", "", locus_id))]
co[, label := paste0(metabolite, " | chr", chr, ":", format(lead_pos, big.mark = ",", trim = TRUE))]
co[, label := factor(label, levels = rev(unique(co[order(PP.H4)]$label)))]
sup_lab <- c("weak_H4_lt_0.25" = "Weak (PP.H4 < 0.25)", "suggestive_H4_0.25_0.5" = "Suggestive (0.25-0.50)",
             "moderate_H4_ge_0.5" = "Moderate (0.50-0.79)", "strong_H4_ge_0.8" = "Strong (>= 0.80)")
co[, support := factor(sup_lab[coloc_support], levels = sup_lab)]
p4 <- ggplot(co, aes(PP.H4, label, fill = support)) +
  geom_col(width = 0.72) +
  geom_vline(xintercept = c(0.5, 0.8), linetype = "dashed", linewidth = 0.35, color = "grey50") +
  facet_wrap(~ outcome, scales = "free_y") +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = setNames(c("#8A8A8A", "#B7791F", "#1F7A6D", "#0B4F6C"), sup_lab), drop = TRUE) +
  labs(x = "Colocalization posterior probability of a shared variant (PP.H4)", y = NULL, fill = NULL,
       title = "Colocalization of metabolite loci with ovarian-reserve traits",
       subtitle = "Top two metabolite-associated loci per prioritized metabolite; loci with PP.H4 >= 0.01 shown") +
  theme_classic(base_size = 9.5) +
  theme(legend.position = "bottom", strip.background = element_rect(fill = "grey94", color = "grey80"),
        strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold"))
ggsave(file.path(out_dir, "Supplementary_Figure_S4_ovarian_colocalization.png"), p4, width = 8.2, height = 4.8, dpi = DPI, bg = "white")

# ---- S3: metabolite-CAD colocalization ------------------------
cc <- fread(file.path(root, "results/coloc_metabolite_cad/metabolite_cad_coloc_summary.tsv"))[status == "ok"]
cc[, label := paste0(gsub("_", " ", metabolite), " | chr", chr, ":", lead_snp)]
cc[, label := factor(label, levels = rev(unique(cc[order(-PP.H4)]$label)))]
lab3 <- c("weak_H4_lt_0.25" = "Weak (< 0.25)", "suggestive_H4_0.25_0.5" = "Suggestive (0.25-0.49)",
          "moderate_H4_ge_0.5" = "Moderate (0.50-0.79)", "strong_H4_ge_0.8" = "Strong (>= 0.80)")
cc[, support := factor(lab3[coloc_support], levels = lab3)]
p3 <- ggplot(cc, aes(PP.H4, label, fill = support)) + geom_col(width = 0.72) +
  geom_vline(xintercept = c(0.5, 0.8), linetype = "dashed", linewidth = 0.35, color = "grey50") +
  scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = setNames(c("#8A8A8A", "#B7791F", "#1F7A6D", "#0B4F6C"), lab3), drop = FALSE) +
  labs(x = "Posterior probability of a shared variant (PP.H4, coloc.abf)", y = NULL, fill = NULL,
       title = sprintf("Metabolite-CAD colocalization: %d of %d locus pairs with PP.H4 >= 0.80", sum(cc$PP.H4 >= 0.8), nrow(cc))) +
  theme_classic(base_size = 9) + theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
ggsave(file.path(out_dir, "Supplementary_Figure_S3_CAD_colocalization.png"), p3, width = 7.5, height = 7.2, dpi = DPI, bg = "white")

# ---- S2: integrated evidence matrix, all 18 prioritized metabolites ---------------------
ev <- fread(file.path(root, "results/final_integrated_evidence/metabolite_integrated_evidence.tsv"))
ev <- ev[evidence_tier != "not_prioritized"]
stopifnot(nrow(ev) == 18)
mv <- fread(file.path(root, "results/mvmr_axis/metabolic_axis_mvmr_results.tsv"))
ev[, tier := tier_of(evidence_tier)]
setorder(ev, tier, CAD_FDR)
lv <- rev(ev$metabolite)

cap <- 4
mr <- rbind(
  ev[, .(metabolite, tier, layer = "ANM MR", v = -log10(ANM_FDR) * sign(ANM_beta), sig = ANM_FDR < 0.05)],
  ev[, .(metabolite, tier, layer = "CAD MR", v = -log10(CAD_FDR) * sign(CAD_beta), sig = CAD_FDR < 0.05)],
  ev[, .(metabolite, tier, layer = "AMH MR", v = -log10(AMH_FDR) * sign(AMH_beta), sig = AMH_FDR < 0.05)]
)
mr[, v := pmin(pmax(v, -cap), cap)]
cl <- rbind(
  ev[, .(metabolite, tier, layer = "ANM/AMH coloc", v = ovarian_coloc_best_PPH4)],
  ev[, .(metabolite, tier, layer = "CAD coloc", v = CAD_coloc_best_PPH4)]
)
mv_status <- function(out) {
  s <- mv[outcome == out, .(metabolite = exposure, st = fifelse(fdr_within_outcome < 0.05, "Independent (FDR < 0.05)", "Not independent"))]
  x <- merge(ev[, .(metabolite, tier)], s, by = "metabolite", all.x = TRUE)
  x[is.na(st), st := "Not in MVMR model"]
  x[, layer := paste0("MVMR (", out, ")")]
  x
}
ct <- rbind(
  ev[, .(metabolite, tier, layer = "Reverse MR",
         st = fifelse(!is.na(reverse_MR_FDR) & reverse_MR_FDR < 0.05, "Reverse association (FDR < 0.05)", "No reverse association"))],
  mv_status("ANM"), mv_status("CAD")
)
layers <- c("ANM MR", "CAD MR", "AMH MR", "Reverse MR", "MVMR (ANM)", "MVMR (CAD)", "ANM/AMH coloc", "CAD coloc")
for (d in list(mr, cl, ct)) {
  d[, metabolite := factor(metabolite, levels = lv)]
  d[, layer := factor(layer, levels = layers)]
}

p2 <- ggplot() +
  geom_tile(data = mr, aes(layer, metabolite, fill = v), color = "white", linewidth = 0.6) +
  geom_text(data = mr[sig == TRUE], aes(layer, metabolite, label = "*"), size = 4, vjust = 0.75) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0, limits = c(-cap, cap),
                       name = paste0("MR: signed -log10 FDR\n(capped at +/-", cap, ")")) +
  new_scale_fill() +
  geom_tile(data = ct, aes(layer, metabolite, fill = st), color = "white", linewidth = 0.6) +
  scale_fill_manual(values = c("No reverse association" = "#C7E9C0", "Reverse association (FDR < 0.05)" = "#E6550D",
                               "Independent (FDR < 0.05)" = "#31A354", "Not independent" = "#E0E0E0",
                               "Not in MVMR model" = "#FAFAFA"), name = "Reverse MR / MVMR") +
  new_scale_fill() +
  geom_tile(data = cl, aes(layer, metabolite, fill = v), color = "white", linewidth = 0.6) +
  scale_fill_gradient(low = "#F7FBFF", high = "#08306B", limits = c(0, 1), na.value = "#D9D9D9",
                      name = "Colocalization:\nbest PP.H4\n(grey = not tested)") +
  geom_text(data = rbind(ct[st == "Not in MVMR model", .(tier, layer, metabolite)], cl[is.na(v), .(tier, layer, metabolite)]),
            aes(layer, metabolite, label = "n/a"), size = 2.6, color = "grey45") +
  facet_grid(tier ~ ., scales = "free_y", space = "free_y") +
  labs(x = NULL, y = NULL, title = "Integrated evidence for the 18 prioritized metabolites",
       subtitle = "* MR FDR < 0.05; n/a = not tested or not in the MVMR model. Rows grouped by manuscript evidence tier and ordered by CAD FDR.") +
  theme_minimal(base_size = 10) +
  theme(axis.text.x = element_text(angle = 40, hjust = 1, face = "bold"), axis.text.y = element_text(face = "bold"),
        strip.text.y = element_text(angle = 0, face = "bold"), panel.grid = element_blank(),
        plot.title = element_text(face = "bold"), legend.title = element_text(size = 8.5), legend.text = element_text(size = 8))
ggsave(file.path(out_dir, "Supplementary_Figure_S2_evidence_matrix.png"), p2, width = 9.5, height = 7, dpi = DPI, bg = "white")

message("Wrote S1, S2, S4, S5 to ", out_dir)
