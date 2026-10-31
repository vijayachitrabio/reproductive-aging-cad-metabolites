#!/usr/bin/env Rscript
# Two-sample MR of 56 NMR metabolites on age at natural menopause (ANM).

source(file.path(Sys.getenv("RAC_ROOT", "."), "config", "paths.R"))

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
  library(dplyr)
})

root <- RAC_ROOT
out_dir <- file.path(root, "results")
if(!dir.exists(out_dir)) dir.create(out_dir, recursive=TRUE)

cat("1. Loading Local mQTL Exposure Instruments...\n")
mqtl_file <- file.path(INSTRUMENT_DIR, "data/mqtl/mqtl_gwas_instruments.csv")
mqtl_dat <- fread(mqtl_file)

# Create chr:pos SNP identifier to match ANM formatting
mqtl_dat$SNP_chrpos <- paste0("chr", mqtl_dat$chr, ":", mqtl_dat$pos)

# Format for TwoSampleMR
exp_dat <- format_data(
    as.data.frame(mqtl_dat),
    type = "exposure",
    phenotype_col = "metabolite",
    snp_col = "SNP_chrpos",
    beta_col = "beta",
    se_col = "se",
    eaf_col = "eaf",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "pval"
)
cat(sprintf("Loaded %d instruments across %d unique metabolites.\n", nrow(exp_dat), length(unique(exp_dat$exposure))))

cat("2. Loading Local ANM Outcome Data...\n")
anm_file <- file.path(REPRO_GWAS_DIR, "data/raw/ANM.txt.gz")
anm_dat_full <- fread(anm_file)

# Filter ANM data to only the SNPs present in the exposure to save memory/time
target_snps <- unique(exp_dat$SNP)
anm_dat_filtered <- anm_dat_full[anm_dat_full$SNP %in% target_snps, ]

out_dat <- format_data(
    as.data.frame(anm_dat_filtered),
    type = "outcome",
    snps = target_snps,
    snp_col = "SNP",
    beta_col = "Effect",
    se_col = "SE",
    eaf_col = "EAF",
    effect_allele_col = "Effect_Allele",
    other_allele_col = "Other_Allele",
    pval_col = "Pval"
)
out_dat$outcome <- "Age at Natural Menopause (ANM)"
cat(sprintf("Found %d out of %d required SNPs in ANM outcome data.\n", nrow(out_dat), length(target_snps)))

cat("3. Harmonising Data...\n")
dat <- harmonise_data(exposure_dat = exp_dat, outcome_dat = out_dat)

cat("4. Running High-Throughput MR...\n")
res <- mr(dat, method_list=c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))

cat("5. Applying FDR Correction...\n")
# Calculate FDR P-values for IVW (main analysis)
ivw_res <- res[res$method == "Inverse variance weighted", ]
ivw_res$fdr <- p.adjust(ivw_res$pval, method="fdr")

# Merge FDR back
res <- merge(res, ivw_res[, c("id.exposure", "fdr")], by="id.exposure", all.x=TRUE)

out_file <- file.path(out_dir, "metabolome_to_anm_mr_results.tsv")
fwrite(res, out_file, sep="\t")
cat(sprintf("MR Complete! Results saved to %s\n", out_file))

# Output a quick summary of significant hits (FDR < 0.05)
sig_hits <- ivw_res[ivw_res$fdr < 0.05, ]
if(nrow(sig_hits) > 0) {
    cat(sprintf("\nFound %d significantly causal metabolites!\n", nrow(sig_hits)))
    print(sig_hits[, c("exposure", "b", "pval", "fdr")])
} else {
    cat("\nNo metabolies showed significant causal effects after FDR correction.\n")
}
