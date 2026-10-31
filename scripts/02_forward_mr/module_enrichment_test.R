#!/usr/bin/env Rscript
# Biochemical-module enrichment test (Supplementary Table S2).
#
# The Nightingale panel is a small, pre-annotated, named biochemical panel, so
# enrichment is tested as a formal over-representation test (Fisher's exact /
# hypergeometric) of the prioritized (18) vs. full (56) metabolite set across the
# modules assigned in assign_metabolite_modules.R, each module mapped onto its
# corresponding KEGG/SMPDB pathway for interpretability.

source(file.path(Sys.getenv("RAC_ROOT", "."), "config", "paths.R"))

suppressPackageStartupMessages(library(data.table))

root <- RAC_ROOT
out_dir <- file.path(root, "results", "pathway_enrichment")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

d <- fread(file.path(root, "results/final_integrated_evidence/metabolite_module_assignments.tsv"))

# Map each internal module to its corresponding real pathway resource entry.
# KEGG map numbers are pathway-map IDs (map00xxx); SMPDB IDs are the Small
# Molecule Pathway Database's super-pathway groupings. These are the standard
# public pathway identifiers.
pathway_map <- fread(text = "module\tkegg_pathway\tsmpdb_superpathway
Cholesterol_lipoprotein_C\tmap04979 (Cholesterol metabolism)\tSMP0000004 (Cholesterol Metabolism)
Apolipoprotein\tmap04979 (Cholesterol metabolism)\tSMP0000004 (Cholesterol Metabolism)
Triglyceride_rich_lipoprotein\tmap00561 (Glycerolipid metabolism)\tSMP0000009 (Triacylglycerol Biosynthesis)
Lipoprotein_particle_size\tmap04979 (Cholesterol metabolism)\tSMP0000004 (Cholesterol Metabolism)
Fatty_acid\tmap01212 (Fatty acid metabolism)\tSMP0000456 (Fatty Acid Metabolism)
Phospholipid\tmap00564 (Glycerophospholipid metabolism)\tSMP0000015 (Glycerophospholipid Metabolism)
General_lipid\tmap00071 (Fatty acid degradation) / map04979\tSMP0000004 (Cholesterol Metabolism)
Glycolysis_related\tmap00010 (Glycolysis / Gluconeogenesis)\tSMP0000040 (Glycolysis)
Ketone_body\tmap00072 (Synthesis and degradation of ketone bodies)\tSMP0000185 (Ketone Body Metabolism)
Amino_acid\tmap01230 (Biosynthesis of amino acids)\tSMP0000021 (multiple amino acid pathways)
Inflammation\tNot a KEGG metabolic pathway (acute-phase glycoprotein signal)\tNot in SMPDB metabolic set
Other_clinical\tNot a metabolic pathway (routine clinical chemistry markers)\tNot in SMPDB metabolic set
")

d <- merge(d, pathway_map, by = "module", all.x = TRUE)
d[, prioritized := evidence_tier != "not_prioritized"]

n_total <- nrow(d)
n_prioritized_total <- sum(d$prioritized)

modules <- unique(d$module)
rows <- list()
for (m in modules) {
  n_in_module <- sum(d$module == m)
  n_prioritized_in_module <- sum(d$module == m & d$prioritized)
  # 2x2 table: prioritized/not-prioritized x in-module/not-in-module
  tbl <- matrix(c(
    n_prioritized_in_module, n_prioritized_total - n_prioritized_in_module,
    n_in_module - n_prioritized_in_module,
    n_total - n_in_module - (n_prioritized_total - n_prioritized_in_module)
  ), nrow = 2)
  ft <- fisher.test(tbl)
  rows[[m]] <- data.table(
    module = m,
    kegg_pathway = pathway_map[module == m, kegg_pathway],
    smpdb_superpathway = pathway_map[module == m, smpdb_superpathway],
    n_in_module = n_in_module,
    n_prioritized_in_module = n_prioritized_in_module,
    pct_of_module_prioritized = round(100 * n_prioritized_in_module / n_in_module, 1),
    odds_ratio = unname(ft$estimate),
    p_value = ft$p.value
  )
}

out <- rbindlist(rows)
out[, fdr := p.adjust(p_value, method = "BH")]
setorder(out, p_value)

fwrite(out, file.path(out_dir, "pathway_module_enrichment_test.tsv"), sep = "\t")
cat("Wrote", file.path(out_dir, "pathway_module_enrichment_test.tsv"), "\n\n")
print(out[, .(module, kegg_pathway, n_in_module, n_prioritized_in_module, pct_of_module_prioritized, odds_ratio, p_value, fdr)])
