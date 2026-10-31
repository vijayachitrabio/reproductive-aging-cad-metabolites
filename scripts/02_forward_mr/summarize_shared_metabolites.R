#!/usr/bin/env Rscript
# Summarize metabolites associated with both ANM and CAD.

source(file.path(Sys.getenv("RAC_ROOT", "."), "config", "paths.R"))

suppressPackageStartupMessages({
  library(data.table)
})

root <- RAC_ROOT
out_dir <- file.path(root, "results")

anm <- fread(file.path(out_dir, "metabolome_to_anm_mr_results.tsv"))
cad <- fread(file.path(out_dir, "metabolome_to_cad_mr_results.tsv"))

anm_ivw <- anm[method == "Inverse variance weighted",
  .(
    metabolite = exposure,
    anm_beta = b,
    anm_se = se,
    anm_p = pval,
    anm_fdr = fdr,
    anm_nsnp = nsnp
  )
]

cad_ivw <- cad[method == "Inverse variance weighted",
  .(
    metabolite = exposure,
    cad_beta = b,
    cad_se = se,
    cad_p = pval,
    cad_targeted_fdr = fdr,
    cad_nsnp = nsnp
  )
]

shared <- merge(anm_ivw, cad_ivw, by = "metabolite", all = FALSE)
shared[, `:=`(
  anm_fdr_sig = anm_fdr < 0.05,
  cad_targeted_fdr_sig = cad_targeted_fdr < 0.05,
  same_beta_direction = sign(anm_beta) == sign(cad_beta)
)]
shared[, evidence_class := fifelse(
  anm_fdr_sig & cad_targeted_fdr_sig,
  "shared_ivw_fdr_targeted_cad_followup",
  "anm_fdr_only_or_cad_suggestive"
)]

setorder(shared, cad_targeted_fdr, anm_fdr)
fwrite(shared, file.path(out_dir, "anm_cad_shared_metabolite_summary.tsv"), sep = "\t")

cat("ANM IVW metabolites tested:", nrow(anm_ivw), "\n")
cat("ANM IVW FDR<0.05:", sum(anm_ivw$anm_fdr < 0.05), "\n")
cat("CAD follow-up metabolites tested:", nrow(cad_ivw), "\n")
cat("Shared ANM + CAD targeted FDR<0.05:", sum(shared$anm_fdr_sig & shared$cad_targeted_fdr_sig), "\n")
cat("Summary written to results/anm_cad_shared_metabolite_summary.tsv\n")
