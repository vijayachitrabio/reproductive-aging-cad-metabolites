#!/usr/bin/env Rscript
# Assign the 56 NMR metabolites to biochemical modules following standard
# Nightingale NMR panel groupings (not a data-driven clustering), used as an
# interpretive scaffold for the prioritized metabolites (Table 1).

source(file.path(Sys.getenv("RAC_ROOT", "."), "config", "paths.R"))

suppressPackageStartupMessages(library(data.table))

root <- RAC_ROOT
evidence_file <- file.path(root, "results/cad_all_metabolites/metabolic_paper_evidence_table.tsv")
out_dir <- file.path(root, "results/final_integrated_evidence")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

module_map <- fread(text = "metabolite\tmodule
Total_C\tCholesterol_lipoprotein_C
Total_FC\tCholesterol_lipoprotein_C
Remnant_C\tCholesterol_lipoprotein_C
VLDL_C\tCholesterol_lipoprotein_C
LDL_C\tCholesterol_lipoprotein_C
HDL_C\tCholesterol_lipoprotein_C
Clinical_LDL_C\tCholesterol_lipoprotein_C
non_HDL_C\tCholesterol_lipoprotein_C
ApoA1\tApolipoprotein
ApoB\tApolipoprotein
ApoB_by_ApoA1\tApolipoprotein
Total_TG\tTriglyceride_rich_lipoprotein
TG_by_PG\tTriglyceride_rich_lipoprotein
VLDL_size\tLipoprotein_particle_size
LDL_size\tLipoprotein_particle_size
HDL_size\tLipoprotein_particle_size
Total_P\tLipoprotein_particle_size
Total_FA\tFatty_acid
SFA\tFatty_acid
MUFA\tFatty_acid
PUFA\tFatty_acid
PUFA_by_MUFA\tFatty_acid
Omega_3\tFatty_acid
Omega_6\tFatty_acid
Omega_6_by_Omega_3\tFatty_acid
LA\tFatty_acid
DHA\tFatty_acid
Unsaturation\tFatty_acid
Total_PL\tPhospholipid
Phosphatidylc\tPhospholipid
Sphingomyelins\tPhospholipid
Phosphoglyc\tPhospholipid
Cholines\tPhospholipid
Total_L\tGeneral_lipid
Total_CE\tGeneral_lipid
Glucose\tGlycolysis_related
Lactate\tGlycolysis_related
Pyruvate\tGlycolysis_related
Citrate\tGlycolysis_related
Acetate\tKetone_body
Acetone\tKetone_body
Acetoacetate\tKetone_body
bOHbutyrate\tKetone_body
Ala\tAmino_acid
Gln\tAmino_acid
Gly\tAmino_acid
His\tAmino_acid
Ile\tAmino_acid
Leu\tAmino_acid
Val\tAmino_acid
Phe\tAmino_acid
Tyr\tAmino_acid
Total_BCAA\tAmino_acid
GlycA\tInflammation
Albumin\tOther_clinical
Creatinine\tOther_clinical
")

message("Loading metabolic paper evidence table")
evidence <- fread(evidence_file)

missing_module <- setdiff(evidence$metabolite, module_map$metabolite)
if (length(missing_module) > 0) {
  warning("No module assignment for: ", paste(missing_module, collapse = ", "))
}

merged <- merge(evidence, module_map, by = "metabolite", all.x = TRUE)
setorder(merged, module, -evidence_tier)
fwrite(merged, file.path(out_dir, "metabolite_module_assignments.tsv"), sep = "\t")

message("Building module-level summary")
merged[, anm_fdr_sig := !is.na(anm_fdr) & anm_fdr < 0.05]
merged[, cad_global_fdr_sig := !is.na(cad_global_fdr) & cad_global_fdr < 0.05]
merged[, prioritized := evidence_tier != "not_prioritized"]

module_summary <- merged[
  ,
  .(
    n_metabolites = .N,
    n_anm_fdr_sig = sum(anm_fdr_sig, na.rm = TRUE),
    n_cad_global_fdr_sig = sum(cad_global_fdr_sig, na.rm = TRUE),
    n_prioritized_tier1_2_3 = sum(prioritized, na.rm = TRUE),
    prioritized_metabolites = paste(sort(metabolite[prioritized]), collapse = ", ")
  ),
  by = module
]
setorder(module_summary, -n_prioritized_tier1_2_3)
fwrite(module_summary, file.path(out_dir, "metabolite_module_summary.tsv"), sep = "\t")

n_modules_with_hits <- module_summary[n_prioritized_tier1_2_3 > 0, .N]
message("Done.")
message("Prioritized (tier 1-3) metabolites span ", n_modules_with_hits, " of ",
        nrow(module_summary), " biochemical modules — see metabolite_module_summary.tsv")
message("Interpretation: within-module hits should be read as one correlated signal, not independent replications.")
print(module_summary[n_prioritized_tier1_2_3 > 0])
