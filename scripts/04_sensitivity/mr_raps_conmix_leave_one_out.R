#!/usr/bin/env Rscript
# MR-RAPS, contamination-mixture and leave-one-out analyses for the 18
# prioritized metabolites against ANM and CAD (Supplementary Table S3b).
# Uses the cached harmonised SNP sets; no new API calls.

suppressPackageStartupMessages({
  library(data.table)
  library(MendelianRandomization)
  library(mr.raps)
  library(TwoSampleMR)
})

root <- normalizePath(getwd(), mustWork = TRUE)
out_dir <- file.path(root, "results", "mr_robustness")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

evidence <- fread(file.path(root, "results/final_integrated_evidence/metabolite_integrated_evidence.tsv"))
prioritized <- evidence[evidence_tier != "not_prioritized", metabolite]
cat("Prioritized metabolites (n=", length(prioritized), "):\n", paste(prioritized, collapse=", "), "\n\n", sep="")

anm <- fread(file.path(root, "results/sensitivity_directionality/metabolite_anm_harmonised.tsv.gz"))
cad <- fread(file.path(root, "results/cad_all_metabolites/metabolite_cad_harmonised.tsv.gz"))

run_robustness_one <- function(dat, exposure_name, outcome_label) {
  sub <- dat[exposure == exposure_name & mr_keep == TRUE]
  sub <- sub[!is.na(beta.exposure) & !is.na(beta.outcome) & !is.na(se.exposure) & !is.na(se.outcome)]
  n <- nrow(sub)
  if (n < 3) {
    return(data.table(metabolite = exposure_name, outcome = outcome_label, n_snps = n,
                       conmix_beta = NA_real_, conmix_ci_lower = NA_real_, conmix_ci_upper = NA_real_,
                       raps_beta = NA_real_, raps_se = NA_real_, raps_p = NA_real_,
                       loo_min_p = NA_real_, loo_max_p = NA_real_, loo_flag = "too_few_snps"))
  }

  # Contamination mixture (MendelianRandomization package)
  conmix_res <- tryCatch({
    mrinput <- mr_input(bx = sub$beta.exposure, bxse = sub$se.exposure,
                         by = sub$beta.outcome, byse = sub$se.outcome, snps = sub$SNP)
    mr_conmix(mrinput, alpha = 0.05)
  }, error = function(e) NULL)

  conmix_beta <- if (!is.null(conmix_res)) conmix_res@Estimate else NA_real_
  conmix_ci <- if (!is.null(conmix_res) && length(conmix_res@CIRange) > 0) {
    c(min(conmix_res@CIRange), max(conmix_res@CIRange))
  } else c(NA_real_, NA_real_)

  # MR-RAPS
  raps_res <- tryCatch({
    mr.raps::mr.raps.overdispersed.robust(
      b_exp = sub$beta.exposure, b_out = sub$beta.outcome,
      se_exp = sub$se.exposure, se_out = sub$se.outcome
    )
  }, error = function(e) tryCatch({
    mr.raps::mr.raps(b_exp = sub$beta.exposure, b_out = sub$beta.outcome,
                      se_exp = sub$se.exposure, se_out = sub$se.outcome)
  }, error = function(e2) NULL))

  raps_beta <- if (!is.null(raps_res)) raps_res$beta.hat else NA_real_
  raps_se <- if (!is.null(raps_res)) raps_res$beta.se else NA_real_
  raps_p <- if (!is.null(raps_res) && !is.na(raps_se) && raps_se > 0) {
    2 * pnorm(-abs(raps_beta / raps_se))
  } else NA_real_

  # Leave-one-out (IVW), using TwoSampleMR's built-in function on a minimal reconstructed dataset
  loo_res <- tryCatch({
    tsmr_dat <- sub
    tsmr_dat$id.exposure <- exposure_name
    tsmr_dat$id.outcome <- outcome_label
    tsmr_dat$outcome <- outcome_label
    tsmr_dat$mr_keep <- TRUE
    mr_leaveoneout(tsmr_dat)
  }, error = function(e) NULL)

  loo_p_range <- if (!is.null(loo_res) && nrow(loo_res) > 0) {
    range(loo_res$p[!is.na(loo_res$p)])
  } else c(NA_real_, NA_real_)

  full_ivw_p <- tryCatch({
    if (!is.null(loo_res)) loo_res$p[loo_res$SNP == "All"] else NA_real_
  }, error = function(e) NA_real_)
  loo_flag <- if (!is.null(loo_res) && length(full_ivw_p) == 1 && !is.na(full_ivw_p)) {
    single_snp_p <- loo_res$p[loo_res$SNP != "All"]
    if (any(!is.na(single_snp_p)) && (full_ivw_p < 0.05) != all(single_snp_p < 0.05, na.rm = TRUE)) {
      "sensitive_to_single_snp_removal"
    } else "stable"
  } else "not_computed"

  data.table(
    metabolite = exposure_name, outcome = outcome_label, n_snps = n,
    conmix_beta = conmix_beta, conmix_ci_lower = conmix_ci[1], conmix_ci_upper = conmix_ci[2],
    raps_beta = raps_beta, raps_se = raps_se, raps_p = raps_p,
    loo_min_p = loo_p_range[1], loo_max_p = loo_p_range[2], loo_flag = loo_flag
  )
}

results_list <- list()
for (m in prioritized) {
  cat("Running robustness for:", m, "\n")
  results_list[[paste0(m, "_ANM")]] <- run_robustness_one(anm, m, "ANM")
  results_list[[paste0(m, "_CAD")]] <- run_robustness_one(cad, m, "CAD")
}

out <- rbindlist(results_list, fill = TRUE)
out <- merge(out, evidence[, .(metabolite, evidence_tier, manuscript_priority, module_class)], by = "metabolite", all.x = TRUE)
setcolorder(out, c("metabolite", "outcome", "evidence_tier", "manuscript_priority", "module_class"))

fwrite(out, file.path(out_dir, "mr_robustness_raps_conmix_loo.tsv"), sep = "\t")
cat("\nWrote", file.path(out_dir, "mr_robustness_raps_conmix_loo.tsv"), "\n")
print(out[, .(metabolite, outcome, n_snps, conmix_beta, raps_beta, raps_p, loo_flag)])
