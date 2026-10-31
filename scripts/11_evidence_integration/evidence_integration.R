# Integrate MR, sensitivity, colocalization and LDSC evidence per metabolite.

source(file.path(Sys.getenv("RAC_ROOT", "."), "config", "paths.R"))

library(data.table)
library(ggplot2)

# Set paths
base_dir <- RAC_ROOT
setwd(base_dir)

out_dir <- "results/final_integrated_evidence"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(paste0(out_dir, "/figures"), showWarnings = FALSE, recursive = TRUE)

# 1. Load core evidence table
cat("Loading MR base evidence...\n")
base_dt <- fread("results/cad_all_metabolites/metabolic_paper_evidence_table.tsv")
setnames(base_dt, "anm_beta", "ANM_beta")
setnames(base_dt, "anm_fdr", "ANM_FDR")
setnames(base_dt, "cad_beta", "CAD_beta")
setnames(base_dt, "cad_global_fdr", "CAD_FDR")

# 2. Load AMH sensitivity
cat("Loading AMH sensitivity...\n")
if(file.exists("results/metabolome_to_amh/metabolite_to_amh_summary_with_evidence.tsv")) {
  amh_dt <- fread("results/metabolome_to_amh/metabolite_to_amh_summary_with_evidence.tsv")
  base_dt <- merge(base_dt, amh_dt[, .(metabolite, AMH_beta = amh_beta, AMH_FDR = amh_fdr_all_metabolites)], by = "metabolite", all.x = TRUE)
} else {
  base_dt[, `:=`(AMH_beta = NA_real_, AMH_FDR = NA_real_)]
}

# 3. Load Reverse MR (both directions: ANM->metabolite and CAD->metabolite)
cat("Loading Reverse MR (both directions)...\n")
if(file.exists("results/reverse_mr_anm_to_metabolites/anm_to_prioritized_metabolites_reverse_summary.tsv")) {
  rev_dt <- fread("results/reverse_mr_anm_to_metabolites/anm_to_prioritized_metabolites_reverse_summary.tsv")
  base_dt <- merge(base_dt, rev_dt[, .(metabolite, reverse_ANM_MR_FDR = reverse_fdr)], by = "metabolite", all.x = TRUE)
} else {
  base_dt[, reverse_ANM_MR_FDR := NA_real_]
}
if(file.exists("results/reverse_mr_cad_to_metabolites/cad_to_prioritized_metabolites_reverse_summary.tsv")) {
  rev_cad_dt <- fread("results/reverse_mr_cad_to_metabolites/cad_to_prioritized_metabolites_reverse_summary.tsv")
  base_dt <- merge(base_dt, rev_cad_dt[, .(metabolite, reverse_CAD_MR_FDR = reverse_fdr)], by = "metabolite", all.x = TRUE)
} else {
  base_dt[, reverse_CAD_MR_FDR := NA_real_]
}
# reverse_MR_FDR kept for backward compatibility with the tier logic below: worst (most significant) of the two directions
base_dt[, reverse_MR_FDR := pmin(reverse_ANM_MR_FDR, reverse_CAD_MR_FDR, na.rm = TRUE)]

# 4. Load MVMR
cat("Loading MVMR...\n")
if(file.exists("results/mvmr_axis/metabolic_axis_mvmr_results.tsv")) {
  mvmr_dt <- fread("results/mvmr_axis/metabolic_axis_mvmr_results.tsv")
  mvmr_sig <- mvmr_dt[fdr_within_outcome < 0.05, .(MVMR_support = paste(unique(outcome), collapse=", ")), by = .(metabolite = exposure)]
  base_dt <- merge(base_dt, mvmr_sig, by = "metabolite", all.x = TRUE)
} else {
  base_dt[, MVMR_support := NA_character_]
}

# 5. Load Ovarian Coloc
cat("Loading Ovarian Coloc...\n")
if(file.exists("results/coloc_metabolite_ovarian/metabolite_ovarian_coloc_pilot_summary.tsv")) {
  ov_coloc <- fread("results/coloc_metabolite_ovarian/metabolite_ovarian_coloc_pilot_summary.tsv")
  ov_best <- ov_coloc[, .SD[which.max(PP.H4)], by = metabolite]
  base_dt <- merge(base_dt, ov_best[, .(metabolite, ovarian_coloc_best_trait = outcome, ovarian_coloc_best_PPH4 = PP.H4)], by = "metabolite", all.x = TRUE)
} else {
  base_dt[, `:=`(ovarian_coloc_best_trait = NA_character_, ovarian_coloc_best_PPH4 = NA_real_)]
}

# 6. Load CAD Coloc
cat("Loading CAD Coloc...\n")
if(file.exists("results/coloc_metabolite_cad/metabolite_cad_coloc_summary.tsv")) {
  cad_coloc <- fread("results/coloc_metabolite_cad/metabolite_cad_coloc_summary.tsv")
  cad_best <- cad_coloc[, .SD[which.max(PP.H4)], by = metabolite]
  base_dt <- merge(base_dt, cad_best[, .(metabolite, CAD_coloc_best_PPH4 = PP.H4)], by = "metabolite", all.x = TRUE)
} else {
  base_dt[, CAD_coloc_best_PPH4 := NA_real_]
}

# Add module/class definitions from the 12-module Nightingale-panel assignment
# (02_forward_mr/assign_metabolite_modules.R output).
cat("Assigning modules (from assign_metabolite_modules.R output)...\n")
module_map <- fread("results/final_integrated_evidence/metabolite_module_assignments.tsv")[, .(metabolite, module_class = module)]
base_dt <- merge(base_dt, module_map, by = "metabolite", all.x = TRUE)

# Merge MR-PRESSO pleiotropy/outlier-corrected estimates (metabolite -> CAD direction)
cat("Loading MR-PRESSO results...\n")
if (file.exists("results/mr_presso/mr_presso_metabolite_to_CAD.tsv")) {
  presso_dt <- fread("results/mr_presso/mr_presso_metabolite_to_CAD.tsv")
  base_dt <- merge(
    base_dt,
    presso_dt[, .(
      metabolite = group,
      CAD_presso_global_test_p = global_test_p,
      CAD_presso_n_outliers = n_outliers,
      CAD_outlier_corrected_beta = outlier_corrected_beta
    )],
    by = "metabolite", all.x = TRUE
  )
} else {
  base_dt[, `:=`(CAD_presso_global_test_p = NA_real_, CAD_presso_n_outliers = NA_integer_, CAD_outlier_corrected_beta = NA_real_)]
}

# Direction-robustness flag: LA and GlycA showed instrument-set-dependent CAD
# direction (single strongest instrument disagrees with the pooled IVW estimate,
# and MR-PRESSO outlier-correction does not resolve it toward the single-SNP
# direction). LA is supported by a published MR study using the same CAD GWAS;
# GlycA is flagged only.
base_dt[, direction_sensitivity_flag := NA_character_]
base_dt[metabolite == "LA", direction_sensitivity_flag :=
  "CAUTION: full-instrument CAD direction (+) contradicts the single strongest instrument (FADS locus rs174564, protective direction), matching a published same-CAD-GWAS MR study. Report both estimates; do not present as confident Tier-1."]
base_dt[metabolite == "GlycA", direction_sensitivity_flag :=
  "CAUTION: full-instrument CAD direction (+) contradicts the single strongest instrument (rs77303550, protective direction). Unlike LA, not externally validated -- open flag, not resolved. Report both estimates."]

# Create Priority Tiers
cat("Assigning priority tiers...\n")
base_dt[, manuscript_priority := "Uncategorized"]
base_dt[, interpretation_label := ""]

for (i in 1:nrow(base_dt)) {
  row <- base_dt[i]
  
  has_shared_mr <- !is.na(row$ANM_FDR) && row$ANM_FDR < 0.05 && !is.na(row$CAD_FDR) && row$CAD_FDR < 0.05
  has_amh_support <- is.na(row$AMH_FDR) || row$AMH_FDR < 0.05
  has_coloc <- (!is.na(row$ovarian_coloc_best_PPH4) && row$ovarian_coloc_best_PPH4 >= 0.50) || 
               (!is.na(row$CAD_coloc_best_PPH4) && row$CAD_coloc_best_PPH4 >= 0.50)
  has_reverse <- !is.na(row$reverse_MR_FDR) && row$reverse_MR_FDR < 0.05
  
  if (has_shared_mr && has_amh_support && has_coloc && !has_reverse) {
    base_dt[i, manuscript_priority := "Tier A"]
    base_dt[i, interpretation_label := "Strong shared lipid/metabolic genetic architecture across cardiometabolic risk and ovarian reserve"]
  } else if (has_shared_mr && !has_reverse) {
    base_dt[i, manuscript_priority := "Tier B"]
    base_dt[i, interpretation_label := "Supportive MR evidence for shared genetic architecture; colocalization lacking or untested"]
  } else {
    base_dt[i, manuscript_priority := "Tier C"]
    if (has_reverse) {
      base_dt[i, interpretation_label := "Evidence confounded by bidirectional/reverse effects"]
    } else {
      base_dt[i, interpretation_label := "Secondary finding / Nominal significance"]
    }
  }
}

# Sort by priority and then by CAD FDR
setorder(base_dt, manuscript_priority, CAD_FDR)

# Select and order final columns
final_cols <- c(
  "metabolite", "manuscript_priority", "interpretation_label", "evidence_tier", "module_class",
  "ANM_beta", "ANM_FDR", "CAD_beta", "CAD_FDR", "AMH_beta", "AMH_FDR",
  "reverse_ANM_MR_FDR", "reverse_CAD_MR_FDR", "reverse_MR_FDR", "MVMR_support",
  "CAD_presso_global_test_p", "CAD_presso_n_outliers", "CAD_outlier_corrected_beta", "direction_sensitivity_flag",
  "ovarian_coloc_best_trait", "ovarian_coloc_best_PPH4", "CAD_coloc_best_PPH4"
)
out_dt <- base_dt[, ..final_cols]

fwrite(out_dt, paste0(out_dir, "/metabolite_integrated_evidence.tsv"), sep="\t", quote=FALSE)
cat("Wrote evidence matrix to", paste0(out_dir, "/metabolite_integrated_evidence.tsv"), "\n")

# --- Heatmap Generation ---
cat("Generating Heatmap...\n")
plot_dt <- base_dt[manuscript_priority %in% c("Tier A", "Tier B")]
if (nrow(plot_dt) == 0) plot_dt <- base_dt[evidence_tier %like% "tier1|tier2"]

heatmap_data <- data.table(
  metabolite = plot_dt$metabolite,
  priority = plot_dt$manuscript_priority,
  `ANM MR` = -log10(plot_dt$ANM_FDR) * sign(plot_dt$ANM_beta),
  `CAD MR` = -log10(plot_dt$CAD_FDR) * sign(plot_dt$CAD_beta),
  `AMH MR` = -log10(plot_dt$AMH_FDR) * sign(plot_dt$AMH_beta),
  `Rev MR Absent` = as.numeric(is.na(plot_dt$reverse_MR_FDR) | plot_dt$reverse_MR_FDR > 0.05),
  `MVMR` = as.numeric(!is.na(plot_dt$MVMR_support)),
  `Ov Coloc` = plot_dt$ovarian_coloc_best_PPH4,
  `CAD Coloc` = plot_dt$CAD_coloc_best_PPH4
)

# Replace NAs with 0 for plotting
heatmap_data[is.na(`AMH MR`), `AMH MR` := 0]
heatmap_data[is.na(`Ov Coloc`), `Ov Coloc` := 0]
heatmap_data[is.na(`CAD Coloc`), `CAD Coloc` := 0]

# Melt for ggplot
melt_dt <- melt(heatmap_data, id.vars=c("metabolite", "priority"), variable.name="Evidence_Layer", value.name="Score")
melt_dt[, metabolite := factor(metabolite, levels=rev(plot_dt$metabolite))]

# Define standard max ranges
melt_dt[, Scaled_Score := Score]
# For MR, cap at -5 to 5. For Coloc, cap at 0 to 1. For boolean, 0 to 1.
melt_dt[Evidence_Layer %in% c("ANM MR", "CAD MR", "AMH MR"), Scaled_Score := pmin(pmax(Score, -5), 5)] 
# All layers share one diverging colour scale (scale_fill_gradient2).

p <- ggplot(melt_dt, aes(x=Evidence_Layer, y=metabolite, fill=Scaled_Score)) +
  geom_tile(color="white", linewidth=0.5) +
  facet_grid(priority ~ ., scales="free_y", space="free_y") +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle=45, hjust=1, vjust=1, face="bold"),
    axis.text.y = element_text(face="bold"),
    strip.text.y = element_text(angle=0, face="bold", size=10),
    panel.grid = element_blank()
  ) +
  labs(title="Integrated Evidence Matrix", x="", y="", fill="Strength\n(Dir. -log10FDR\nor PP.H4)")

ggsave(paste0(out_dir, "/figures/integrated_evidence_heatmap.png"), p, width=8, height=6, bg="white")
cat("Wrote heatmap to", paste0(out_dir, "/figures/integrated_evidence_heatmap.png"), "\n")

# Dump top prioritized list for report
cat("\nRanked Priority List:\n")
print(out_dt[, .N, by=manuscript_priority])
