#!/usr/bin/env Rscript
# Steiger directionality and MR-Egger/weighted-median sensitivity analyses for
# the prioritized metabolites.

source(file.path(Sys.getenv("RAC_ROOT", "."), "config", "paths.R"))

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
})

root <- RAC_ROOT
out_dir <- file.path(root, "results", "sensitivity_directionality")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mqtl_file <- Sys.getenv(
  "MQTL_INSTRUMENTS",
  file.path(INSTRUMENT_DIR, "data/mqtl/mqtl_gwas_instruments.csv")
)
anm_file <- Sys.getenv(
  "ANM_SUMSTATS",
  file.path(REPRO_GWAS_DIR, "data/raw/ANM.txt.gz")
)
cad_harmonised_file <- file.path(root, "results", "cad_all_metabolites", "metabolite_cad_harmonised.tsv.gz")

cad_ncase <- as.numeric(Sys.getenv("CAD_NCASE", "60801"))
cad_ncontrol <- as.numeric(Sys.getenv("CAD_NCONTROL", "123504"))
cad_prevalence <- as.numeric(Sys.getenv("CAD_PREVALENCE", as.character(cad_ncase / (cad_ncase + cad_ncontrol))))

message("Loading mQTL instruments: ", mqtl_file)
mqtl <- fread(mqtl_file)
required <- c("metabolite", "SNP", "chr", "pos", "beta", "se", "eaf", "effect_allele", "other_allele", "pval", "F_stat", "n")
missing_cols <- setdiff(required, names(mqtl))
if (length(missing_cols) > 0) {
  stop("Missing required mQTL columns: ", paste(missing_cols, collapse = ", "))
}
mqtl[, SNP_chrpos := paste0("chr", chr, ":", pos)]

message("Preparing metabolite -> ANM harmonised data")
exp_anm <- format_data(
  as.data.frame(mqtl),
  type = "exposure",
  phenotype_col = "metabolite",
  snp_col = "SNP_chrpos",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  samplesize_col = "n"
)

anm <- fread(anm_file)
anm <- anm[SNP %in% unique(exp_anm$SNP)]
out_anm <- format_data(
  as.data.frame(anm),
  type = "outcome",
  snps = unique(exp_anm$SNP),
  snp_col = "SNP",
  beta_col = "Effect",
  se_col = "SE",
  eaf_col = "EAF",
  effect_allele_col = "Effect_Allele",
  other_allele_col = "Other_Allele",
  pval_col = "Pval",
  samplesize_col = "N"
)
out_anm$outcome <- "Age at Natural Menopause (ANM)"

harm_anm <- as.data.table(harmonise_data(exp_anm, out_anm))
mqtl_chr_key <- mqtl[, .(
  exposure = metabolite,
  SNP = SNP_chrpos,
  exposure_rsid = SNP,
  exposure_F_stat = F_stat,
  exposure_n = n
)]
harm_anm <- merge(harm_anm, mqtl_chr_key, by = c("exposure", "SNP"), all.x = TRUE)
fwrite(harm_anm, file.path(out_dir, "metabolite_anm_harmonised.tsv.gz"), sep = "\t")

message("Running metabolite -> ANM MR sensitivity outputs")
mr_anm <- as.data.table(mr(
  harm_anm,
  method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")
))
ivw_anm <- mr_anm[method == "Inverse variance weighted"]
ivw_anm[, anm_fdr := p.adjust(pval, method = "fdr")]
mr_anm <- merge(mr_anm, ivw_anm[, .(id.exposure, anm_fdr)], by = "id.exposure", all.x = TRUE)
pleio_anm <- as.data.table(mr_pleiotropy_test(harm_anm))
het_anm <- as.data.table(mr_heterogeneity(harm_anm))
qc_anm <- harm_anm[, .(
  n_harmonised = .N,
  n_mr_keep = sum(mr_keep %in% TRUE),
  n_palindromic = sum(palindromic %in% TRUE, na.rm = TRUE),
  n_ambiguous = sum(ambiguous %in% TRUE, na.rm = TRUE),
  mean_F = mean(exposure_F_stat, na.rm = TRUE),
  min_F = min(exposure_F_stat, na.rm = TRUE),
  median_F = median(exposure_F_stat, na.rm = TRUE)
), by = exposure]

fwrite(mr_anm, file.path(out_dir, "metabolite_to_anm_all_methods.tsv"), sep = "\t")
fwrite(pleio_anm, file.path(out_dir, "metabolite_to_anm_egger_intercepts.tsv"), sep = "\t")
fwrite(het_anm, file.path(out_dir, "metabolite_to_anm_heterogeneity.tsv"), sep = "\t")
fwrite(qc_anm, file.path(out_dir, "metabolite_to_anm_instrument_qc.tsv"), sep = "\t")

add_quant_r <- function(dat, side) {
  beta_col <- paste0("beta.", side)
  se_col <- paste0("se.", side)
  n_col <- paste0("samplesize.", side)
  r_col <- paste0("r.", side)
  dat[[r_col]] <- get_r_from_bsen(dat[[beta_col]], dat[[se_col]], dat[[n_col]])
  dat
}

message("Running Steiger directionality for metabolite -> ANM")
steiger_anm_dat <- copy(harm_anm[mr_keep == TRUE])
steiger_anm_dat <- add_quant_r(steiger_anm_dat, "exposure")
steiger_anm_dat <- add_quant_r(steiger_anm_dat, "outcome")
steiger_anm <- as.data.table(directionality_test(as.data.frame(steiger_anm_dat)))
fwrite(steiger_anm, file.path(out_dir, "steiger_metabolite_to_anm.tsv"), sep = "\t")

if (file.exists(cad_harmonised_file)) {
  message("Preparing cached metabolite -> CAD harmonised data for Steiger")
  harm_cad <- fread(cad_harmonised_file)
  mqtl_rsid_key <- mqtl[, .(
    exposure = metabolite,
    SNP,
    samplesize.exposure.from_mqtl = n,
    exposure_F_stat = F_stat
  )]
  harm_cad <- merge(harm_cad, mqtl_rsid_key, by = c("exposure", "SNP"), all.x = TRUE)
  if (!"samplesize.exposure" %in% names(harm_cad)) {
    harm_cad[, samplesize.exposure := samplesize.exposure.from_mqtl]
  }
  harm_cad[is.na(samplesize.exposure), samplesize.exposure := samplesize.exposure.from_mqtl]
  harm_cad[, r.exposure := get_r_from_bsen(beta.exposure, se.exposure, samplesize.exposure)]
  harm_cad[, r.outcome := get_r_from_lor(
    beta.outcome,
    eaf.outcome,
    cad_ncase,
    cad_ncontrol,
    cad_prevalence
  )]
  fwrite(harm_cad, file.path(out_dir, "metabolite_cad_harmonised_with_r.tsv.gz"), sep = "\t")

  steiger_cad <- as.data.table(directionality_test(as.data.frame(harm_cad[mr_keep == TRUE])))
  fwrite(steiger_cad, file.path(out_dir, "steiger_metabolite_to_cad.tsv"), sep = "\t")
} else {
  warning("CAD harmonised file not found; skipping CAD Steiger: ", cad_harmonised_file)
}

message("Building combined directionality summary")
anm_ivw <- ivw_anm[, .(
  metabolite = exposure,
  anm_beta = b,
  anm_p = pval,
  anm_fdr,
  anm_nsnp = nsnp
)]
anm_dir <- steiger_anm[, .(
  metabolite = exposure,
  anm_r2_exposure = snp_r2.exposure,
  anm_r2_outcome = snp_r2.outcome,
  anm_steiger_correct = correct_causal_direction,
  anm_steiger_p = steiger_pval
)]
combined <- merge(anm_ivw, anm_dir, by = "metabolite", all.x = TRUE)

cad_ivw_file <- file.path(root, "results", "cad_all_metabolites", "metabolome_to_cad_ivw_global_and_targeted_fdr.tsv")
if (file.exists(cad_ivw_file)) {
  cad_ivw <- fread(cad_ivw_file)
  cad_part <- cad_ivw[, .(
    metabolite = exposure,
    cad_beta = b,
    cad_p = pval,
    cad_global_fdr,
    cad_anm_prioritized_fdr,
    cad_nsnp = nsnp
  )]
  combined <- merge(combined, cad_part, by = "metabolite", all.x = TRUE)
}
steiger_cad_file <- file.path(out_dir, "steiger_metabolite_to_cad.tsv")
if (file.exists(steiger_cad_file)) {
  cad_dir <- fread(steiger_cad_file)[, .(
    metabolite = exposure,
    cad_r2_exposure = snp_r2.exposure,
    cad_r2_outcome = snp_r2.outcome,
    cad_steiger_correct = correct_causal_direction,
    cad_steiger_p = steiger_pval
  )]
  combined <- merge(combined, cad_dir, by = "metabolite", all.x = TRUE)
}

combined[, directionality_class := fifelse(
  anm_fdr < 0.05 & cad_global_fdr < 0.05 & anm_steiger_correct == TRUE & cad_steiger_correct == TRUE,
  "shared_fdr_and_steiger_supported",
  fifelse(
    anm_fdr < 0.05 & cad_global_fdr < 0.05,
    "shared_fdr_directionality_caution",
    "not_shared_fdr"
  )
)]
setorder(combined, directionality_class, cad_global_fdr, anm_fdr)
fwrite(combined, file.path(out_dir, "metabolic_directionality_summary.tsv"), sep = "\t")

summary_counts <- combined[, .N, by = directionality_class][order(directionality_class)]
fwrite(summary_counts, file.path(out_dir, "metabolic_directionality_summary_counts.tsv"), sep = "\t")

message("Done. Outputs written to: ", out_dir)
print(summary_counts)
