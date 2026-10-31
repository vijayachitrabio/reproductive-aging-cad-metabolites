#!/usr/bin/env Rscript
# Forward MR: 18 prioritized metabolites -> parental lifespan (Timmers et al. 2019),
# using the same methods as for metabolite -> CAD. Outcome data are read from the
# local parental-lifespan summary statistics.

source(file.path(Sys.getenv("RAC_ROOT", "."), "config", "paths.R"))

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
})

root <- RAC_ROOT
out_dir <- file.path(root, "results/lifespan_mr")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

PRIORITIZED_18 <- c("ApoA1","Cholines","DHA","GlycA","LA","MUFA","Omega_3","Omega_6",
                     "Omega_6_by_Omega_3","PUFA","SFA","TG_by_PG","Total_C","Total_CE",
                     "Total_FA","Total_L","Total_PL","Total_TG")

message("Loading metabolite instruments")
mqtl <- fread(file.path(INSTRUMENT_DIR, "data/mqtl/mqtl_gwas_instruments.csv"))
mqtl <- mqtl[metabolite %in% PRIORITIZED_18]
message("Instrument SNPs across 18 metabolites: ", nrow(mqtl), " (", length(unique(mqtl$SNP)), " unique)")

message("Loading parental lifespan outcome data (local file, rsID-matched)")
lifespan <- fread(
  file.path(root, "data/raw_parental_lifespan/lifegen_no_snpid.tsv.gz"),
  select = c("rsid", "a1", "a0", "n", "freq1", "beta1", "se", "p")
)
setnames(lifespan, c("SNP", "effect_allele", "other_allele", "n", "eaf", "beta", "se", "pval"))
lifespan <- lifespan[SNP %in% mqtl$SNP]
message("Parental lifespan rows matched by rsID: ", nrow(lifespan))

exp_dat <- format_data(
  as.data.frame(mqtl), type = "exposure",
  phenotype_col = "metabolite", snp_col = "SNP",
  beta_col = "beta", se_col = "se", eaf_col = "eaf",
  effect_allele_col = "effect_allele", other_allele_col = "other_allele",
  pval_col = "pval", samplesize_col = "n"
)
out_dat <- format_data(
  as.data.frame(lifespan), type = "outcome",
  snp_col = "SNP", beta_col = "beta", se_col = "se", eaf_col = "eaf",
  effect_allele_col = "effect_allele", other_allele_col = "other_allele",
  pval_col = "pval", samplesize_col = "n"
)
out_dat$outcome <- "ParentalLifespan"

harm <- harmonise_data(exp_dat, out_dat, action = 2)
fwrite(harm, file.path(out_dir, "metabolite_lifespan_harmonised.tsv.gz"), sep = "\t")

message("Running MR (IVW, MR-Egger, weighted median) per metabolite")
res <- mr(harm, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
res <- as.data.table(res)
res[, fdr := p.adjust(pval, method = "BH"), by = method]
setorder(res, method, pval)

fwrite(res, file.path(out_dir, "metabolite_to_lifespan_mr_all_methods.tsv"), sep = "\t")

ivw <- res[method == "Inverse variance weighted"]
message("\n=== IVW results, FDR-corrected across 18 metabolites ===")
print(ivw[, .(exposure, nsnp, b, se, pval, fdr)])
message("\nSignificant at FDR<0.05: ", sum(ivw$fdr < 0.05, na.rm = TRUE), " of ", nrow(ivw))
