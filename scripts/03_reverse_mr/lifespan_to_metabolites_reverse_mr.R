#!/usr/bin/env Rscript
# Reverse MR: parental lifespan -> 18 prioritized metabolites, using the same
# methods as for CAD -> metabolite reverse MR. Instruments are clumped locally
# with the 1000 Genomes EUR panel (r2 = 0.001, 10,000-kb window).

source(file.path(Sys.getenv("RAC_ROOT", "."), "config", "paths.R"))

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(ieugwasr)
  library(data.table)
})

root <- RAC_ROOT
atlas_root <- REPRO_GWAS_DIR
plink_bin <- file.path(atlas_root, "tools/plink/plink")
plink_ref <- file.path(atlas_root, "data/lava_ref/g1000_eur")
out_dir <- file.path(root, "results/lifespan_mr")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

PRIORITIZED_18 <- c("ApoA1","Cholines","DHA","GlycA","LA","MUFA","Omega_3","Omega_6",
                     "Omega_6_by_Omega_3","PUFA","SFA","TG_by_PG","Total_C","Total_CE",
                     "Total_FA","Total_L","Total_PL","Total_TG")

message("Loading and clumping parental lifespan GWS SNPs")
sig <- fread(file.path(root, "data/raw_parental_lifespan/lifespan_gws_snps.tsv"))
clump_input <- data.frame(rsid = sig$rsid, pval = sig$p)
clumped <- ld_clump(
  clump_input, clump_kb = 10000, clump_r2 = 0.001, clump_p = 5e-8,
  plink_bin = plink_bin, bfile = plink_ref
)
message("Independent lifespan instruments after clumping: ", nrow(clumped))

exp_df <- sig[rsid %in% clumped$rsid]
exp_df <- exp_df[, .(SNP = rsid, effect_allele = toupper(a1), other_allele = toupper(a0),
                      eaf = freq1, beta = beta1, se = se, pval = p, n = n)]
exp_dat <- format_data(as.data.frame(exp_df), type = "exposure", snp_col = "SNP",
                        beta_col = "beta", se_col = "se", eaf_col = "eaf",
                        effect_allele_col = "effect_allele", other_allele_col = "other_allele",
                        pval_col = "pval", samplesize_col = "n")
exp_dat$exposure <- "ParentalLifespan"

message("Loading full metabolite GWAS at these instrument positions (local regenie files)")
mqtl_dir <- file.path(INSTRUMENT_DIR, "data/mqtl/mqtl_full_gwas")

all_res <- list()
for (met in PRIORITIZED_18) {
  full_file <- file.path(mqtl_dir, paste0(met, "_full_regenie.tsv.gz"))
  if (!file.exists(full_file)) { message("MISSING: ", met); next }
  x <- fread(full_file, select = c("ID", "ALLELE0", "ALLELE1", "A1FREQ", "N", "BETA", "SE", "LOG10P"))
  x <- x[ID %in% exp_df$SNP]
  if (nrow(x) == 0) next
  x[, p := 10^(-as.numeric(LOG10P))]
  out_df <- x[, .(SNP = ID, effect_allele = toupper(ALLELE1), other_allele = toupper(ALLELE0),
                   eaf = A1FREQ, beta = BETA, se = SE, pval = p, n = N)]
  out_dat <- format_data(as.data.frame(out_df), type = "outcome", snp_col = "SNP",
                          beta_col = "beta", se_col = "se", eaf_col = "eaf",
                          effect_allele_col = "effect_allele", other_allele_col = "other_allele",
                          pval_col = "pval", samplesize_col = "n")
  out_dat$outcome <- met
  harm <- tryCatch(harmonise_data(exp_dat, out_dat, action = 2), error = function(e) NULL)
  if (is.null(harm) || sum(harm$mr_keep) < 3) { message(met, ": too few harmonised SNPs, skipping"); next }
  res <- tryCatch(mr(harm, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")),
                  error = function(e) NULL)
  if (!is.null(res)) all_res[[met]] <- as.data.table(res)
}

combined <- rbindlist(all_res, fill = TRUE)
ivw <- combined[method == "Inverse variance weighted"]
ivw[, fdr := p.adjust(pval, method = "BH")]
setorder(ivw, pval)

fwrite(combined, file.path(out_dir, "lifespan_to_metabolite_reverse_mr_all_methods.tsv"), sep = "\t")
fwrite(ivw, file.path(out_dir, "lifespan_to_metabolite_reverse_mr_ivw_fdr.tsv"), sep = "\t")

message("\n=== Reverse MR: Parental Lifespan -> metabolites, IVW, FDR-corrected ===")
print(ivw[, .(outcome, nsnp, b, se, pval, fdr)])
message("\nSignificant at FDR<0.05: ", sum(ivw$fdr < 0.05, na.rm = TRUE), " of ", nrow(ivw))
