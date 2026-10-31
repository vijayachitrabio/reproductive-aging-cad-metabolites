#!/usr/bin/env Rscript
# Assign evidence tiers (Tier 1-3) to the prioritized metabolites.

source(file.path(Sys.getenv("RAC_ROOT", "."), "config", "paths.R"))

suppressPackageStartupMessages({
  library(data.table)
})

root <- RAC_ROOT
res_dir <- file.path(root, "results", "cad_all_metabolites")
ivw_file <- file.path(res_dir, "metabolome_to_cad_ivw_global_and_targeted_fdr.tsv")
qc_file <- file.path(res_dir, "metabolite_cad_instrument_qc.tsv")
pleio_file <- file.path(res_dir, "metabolome_to_cad_egger_intercepts.tsv")
het_file <- file.path(res_dir, "metabolome_to_cad_heterogeneity.tsv")

if (!file.exists(ivw_file)) {
  stop("Run scripts/02_forward_mr/metabolites_to_cad_mr_all.R first.")
}

ivw <- fread(ivw_file)
qc <- fread(qc_file)
pleio <- fread(pleio_file)
het <- fread(het_file)

evidence <- ivw[
  ,
  .(
    metabolite = exposure,
    cad_beta = b,
    cad_se = se,
    cad_p = pval,
    cad_global_fdr,
    cad_anm_prioritized_fdr,
    cad_nsnp = nsnp,
    anm_beta,
    anm_p,
    anm_fdr,
    anm_nsnp,
    shared_anm_cad_global_fdr,
    shared_anm_cad_targeted_fdr
  )
]
evidence <- merge(evidence, qc, by = "metabolite", all.x = TRUE)

if (nrow(pleio) > 0) {
  pleio_small <- pleio[, .(metabolite = exposure, egger_intercept, egger_intercept_p = pval)]
  evidence <- merge(evidence, pleio_small, by = "metabolite", all.x = TRUE)
}

het_ivw <- het[method == "Inverse variance weighted", .(metabolite = exposure, ivw_Q = Q, ivw_Q_p = Q_pval)]
evidence <- merge(evidence, het_ivw, by = "metabolite", all.x = TRUE)

evidence[, evidence_tier := fifelse(
  shared_anm_cad_global_fdr == TRUE & n_mr_keep >= 3 & (is.na(egger_intercept_p) | egger_intercept_p >= 0.05),
  "tier1_shared_global_fdr_no_egger_flag",
  fifelse(
    shared_anm_cad_targeted_fdr == TRUE & n_mr_keep >= 3,
    "tier2_shared_targeted_fdr_needs_sensitivity",
    fifelse(
      !is.na(anm_fdr) & anm_fdr < 0.05 & cad_p < 0.05,
      "tier3_nominal_cad_followup",
      "not_prioritized"
    )
  )
)]

setorder(evidence, evidence_tier, cad_global_fdr, cad_p)
fwrite(evidence, file.path(res_dir, "metabolic_paper_evidence_table.tsv"), sep = "\t")

summary <- evidence[
  ,
  .N,
  by = evidence_tier
][order(evidence_tier)]
fwrite(summary, file.path(res_dir, "metabolic_paper_evidence_tier_counts.tsv"), sep = "\t")

cat("Evidence table written to results/cad_all_metabolites/metabolic_paper_evidence_table.tsv\n")
print(summary)
