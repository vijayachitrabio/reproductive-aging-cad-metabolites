#!/usr/bin/env Rscript
# Two-sample MR of the prioritized metabolites on anti-Mullerian hormone (AMH).

source(file.path(Sys.getenv("RAC_ROOT", "."), "config", "paths.R"))

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
})

root <- RAC_ROOT
atlas_root <- REPRO_GWAS_DIR

mqtl_file <- Sys.getenv(
  "MQTL_INSTRUMENTS",
  file.path(INSTRUMENT_DIR, "data/mqtl/mqtl_gwas_instruments.csv")
)
amh_file <- Sys.getenv("AMH_SUMSTATS", file.path(atlas_root, "data/raw/AMH.txt.gz"))
amh_n <- as.numeric(Sys.getenv("AMH_N", "9668"))

out_dir <- file.path(root, "results/metabolome_to_amh")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading mQTL instruments")
mqtl <- fread(mqtl_file)
required <- c("metabolite", "SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "pval", "F_stat", "n")
missing_cols <- setdiff(required, names(mqtl))
if (length(missing_cols) > 0) {
  stop("Missing required mQTL columns: ", paste(missing_cols, collapse = ", "))
}
mqtl <- mqtl[!is.na(SNP) & SNP != ""]

exp_dat <- format_data(
  as.data.frame(mqtl),
  type = "exposure",
  phenotype_col = "metabolite",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  samplesize_col = "n"
)

message("Loading AMH outcome summary statistics")
amh <- fread(
  amh_file,
  select = c(
    "rs_id", "effect_allele", "other_allele", "beta", "standard_error",
    "effect_allele_frequency", "p_value"
  )
)
amh <- amh[rs_id %in% unique(exp_dat$SNP)]
amh[, `:=`(
  effect_allele = toupper(effect_allele),
  other_allele = toupper(other_allele),
  p_value = as.numeric(p_value)
)]

message("AMH overlap rows: ", nrow(amh))
out_dat <- format_data(
  as.data.frame(amh),
  type = "outcome",
  snps = unique(exp_dat$SNP),
  snp_col = "rs_id",
  beta_col = "beta",
  se_col = "standard_error",
  eaf_col = "effect_allele_frequency",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "p_value"
)
out_dat$outcome <- "Circulating AMH"
out_dat$id.outcome <- "AMH_GCST90428625"
out_dat$samplesize.outcome <- amh_n

message("Harmonising metabolite instruments with AMH")
harm <- as.data.table(harmonise_data(exp_dat, out_dat))
mqtl_key <- mqtl[, .(
  exposure = metabolite,
  SNP,
  exposure_F_stat = F_stat,
  exposure_n = n
)]
harm <- merge(harm, mqtl_key, by = c("exposure", "SNP"), all.x = TRUE)
fwrite(harm, file.path(out_dir, "metabolite_amh_harmonised.tsv.gz"), sep = "\t")

message("Running metabolite -> AMH MR")
mr_res <- as.data.table(mr(
  harm,
  method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")
))

ivw <- mr_res[method == "Inverse variance weighted"]
ivw[, amh_fdr_all_metabolites := p.adjust(pval, method = "fdr")]
mr_res <- merge(
  mr_res,
  ivw[, .(id.exposure, amh_fdr_all_metabolites)],
  by = "id.exposure",
  all.x = TRUE
)

pleio <- as.data.table(mr_pleiotropy_test(harm))
het <- as.data.table(mr_heterogeneity(harm))
qc <- harm[
  ,
  .(
    n_harmonised = .N,
    n_mr_keep = sum(mr_keep %in% TRUE),
    n_palindromic = sum(palindromic %in% TRUE, na.rm = TRUE),
    n_ambiguous = sum(ambiguous %in% TRUE, na.rm = TRUE),
    mean_F = mean(exposure_F_stat, na.rm = TRUE),
    min_F = min(exposure_F_stat, na.rm = TRUE),
    median_F = median(exposure_F_stat, na.rm = TRUE)
  ),
  by = exposure
]

message("Running Steiger directionality for metabolite -> AMH")
steiger_dat <- copy(harm[mr_keep == TRUE])
steiger_dat[, r.exposure := get_r_from_bsen(beta.exposure, se.exposure, samplesize.exposure)]
steiger_dat[, r.outcome := get_r_from_bsen(beta.outcome, se.outcome, samplesize.outcome)]
steiger <- as.data.table(directionality_test(as.data.frame(steiger_dat)))

evidence_file <- file.path(root, "results/cad_all_metabolites/metabolic_paper_evidence_table.tsv")
if (file.exists(evidence_file)) {
  evidence <- fread(evidence_file)
  summary <- merge(
    ivw[, .(
      metabolite = exposure,
      amh_beta = b,
      amh_se = se,
      amh_p = pval,
      amh_fdr_all_metabolites,
      amh_nsnp = nsnp
    )],
    evidence[, .(
      metabolite,
      evidence_tier,
      anm_beta,
      anm_fdr,
      cad_beta,
      cad_global_fdr
    )],
    by = "metabolite",
    all.x = TRUE
  )
  summary <- merge(
    summary,
    steiger[, .(
      metabolite = exposure,
      amh_steiger_correct = correct_causal_direction,
      amh_steiger_p = steiger_pval
    )],
    by = "metabolite",
    all.x = TRUE
  )
  summary[, prioritized := !is.na(evidence_tier) & evidence_tier != "not_prioritized"]
  setorder(summary, -prioritized, amh_fdr_all_metabolites, amh_p)
  fwrite(summary, file.path(out_dir, "metabolite_to_amh_summary_with_evidence.tsv"), sep = "\t")

  prioritized_summary <- summary[prioritized == TRUE]
  fwrite(prioritized_summary, file.path(out_dir, "prioritized_metabolite_to_amh_summary.tsv"), sep = "\t")
}

fwrite(mr_res, file.path(out_dir, "metabolite_to_amh_all_methods.tsv"), sep = "\t")
fwrite(ivw, file.path(out_dir, "metabolite_to_amh_ivw.tsv"), sep = "\t")
fwrite(pleio, file.path(out_dir, "metabolite_to_amh_egger_intercepts.tsv"), sep = "\t")
fwrite(het, file.path(out_dir, "metabolite_to_amh_heterogeneity.tsv"), sep = "\t")
fwrite(qc, file.path(out_dir, "metabolite_to_amh_instrument_qc.tsv"), sep = "\t")
fwrite(steiger, file.path(out_dir, "metabolite_to_amh_steiger.tsv"), sep = "\t")

message("Done. Outputs written to ", out_dir)
print(ivw[order(amh_fdr_all_metabolites)][1:20])
if (exists("prioritized_summary")) {
  message("Prioritized metabolite AMH summary")
  print(prioritized_summary[order(amh_fdr_all_metabolites)])
}
