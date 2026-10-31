#!/usr/bin/env Rscript
# Two-sample MR on CAD for the metabolites significant in the ANM screen.

source(file.path(Sys.getenv("RAC_ROOT", "."), "config", "paths.R"))

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
  library(dplyr)
})

if (Sys.getenv("OPENGWAS_JWT") == "") {
  stop("Set OPENGWAS_JWT in the environment before querying OpenGWAS.")
}

root <- RAC_ROOT
out_dir <- file.path(root, "results")

# Load significant metabolites from previous step
cat("Loading previous ANM MR results...\n")
res_anm <- fread(file.path(out_dir, "metabolome_to_anm_mr_results.tsv"))
sig_mets <- unique(res_anm[res_anm$fdr < 0.05 & res_anm$method == "Inverse variance weighted"]$exposure)
cat(sprintf("Found %d significant metabolites to test against CAD.\n", length(sig_mets)))

cat("Loading Local mQTL Exposure Instruments...\n")
mqtl_file <- file.path(INSTRUMENT_DIR, "data/mqtl/mqtl_gwas_instruments.csv")
mqtl_dat <- fread(mqtl_file)

# Filter mQTLs to just the significant ones
mqtl_dat <- mqtl_dat[mqtl_dat$metabolite %in% sig_mets, ]

cat("Querying OpenGWAS for Coronary Artery Disease (ieu-a-7)...\n")
# Using actual rsIDs for OpenGWAS query instead of chr:pos since OpenGWAS uses rsIDs
exp_dat_rsid <- format_data(
    as.data.frame(mqtl_dat),
    type = "exposure",
    phenotype_col = "metabolite",
    snp_col = "SNP",
    beta_col = "beta",
    se_col = "se",
    eaf_col = "eaf",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval"
)

out_dat <- extract_outcome_data(snps = exp_dat_rsid$SNP, outcomes = "ieu-a-7")
out_dat$outcome <- "Coronary Artery Disease"

cat("Harmonising Data...\n")
dat <- harmonise_data(exposure_dat = exp_dat_rsid, outcome_dat = out_dat)

cat("Running MR...\n")
res <- mr(dat, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))

ivw_res <- res[res$method == "Inverse variance weighted", ]
ivw_res$fdr <- p.adjust(ivw_res$pval, method="fdr")

out_file <- file.path(out_dir, "metabolome_to_cad_mr_results.tsv")
fwrite(ivw_res, out_file, sep="\t")
cat(sprintf("\nFound %d metabolites with significant causal effects on CAD (FDR < 0.05)!\n", sum(ivw_res$fdr < 0.05)))
print(ivw_res[ivw_res$fdr < 0.05, c("exposure", "b", "pval", "fdr")])
