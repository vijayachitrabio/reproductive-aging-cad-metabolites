#!/usr/bin/env Rscript
# Total-effect MR of ANM on parental lifespan, using the same methods as for
# ANM -> CAD, joined by chr:pos against the local parental-lifespan summary
# statistics (both files carry chromosome and position).

source(file.path(Sys.getenv("RAC_ROOT", "."), "config", "paths.R"))

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
})

root <- RAC_ROOT
atlas_root <- REPRO_GWAS_DIR
out_dir <- file.path(root, "results/lifespan_mr")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

anm_file <- file.path(atlas_root, "data/raw/ANM.txt.gz")
anm_lead_file <- file.path(atlas_root, "results/phase4_loci/ANM.lead_loci.tsv")

message("Loading ANM lead loci")
lead <- fread(anm_lead_file)
message("ANM lead loci: ", nrow(lead))
# lead_pos and the position embedded in lead_snp differ in this file -- lead_snp's own embedded position is the true top-SNP
# position; lead_pos appears to reflect a different locus-level metric. Parse
# position from lead_snp (format "chrN:POS"), not from lead_pos.
lead[, lead_snp_pos := as.integer(sub("^chr[0-9XY]+:", "", lead_snp))]

message("Loading raw ANM summary statistics at lead loci")
anm <- fread(
  anm_file,
  select = c("SNP", "CHR", "POS", "Effect_Allele", "Other_Allele", "EAF", "Effect", "SE", "Pval", "N")
)
anm[, `:=`(CHR = as.integer(CHR), POS = as.integer(POS))]
anm_lead <- anm[paste0(CHR, ":", POS) %in% paste0(lead$chrom, ":", lead$lead_snp_pos)]
anm_lead <- anm_lead[!duplicated(paste0(CHR, ":", POS))]
message("ANM lead-locus SNPs found in raw ANM file: ", nrow(anm_lead), " of ", nrow(lead))

message("Loading parental lifespan outcome data at the same positions")
lifespan <- fread(
  file.path(root, "data/raw_parental_lifespan/lifegen_no_snpid.tsv.gz"),
  select = c("rsid", "chr", "pos", "a1", "a0", "n", "freq1", "beta1", "se", "p")
)
lifespan[, `:=`(chr = as.integer(chr), pos = as.integer(pos))]
lifespan_lead <- lifespan[paste0(chr, ":", pos) %in% paste0(anm_lead$CHR, ":", anm_lead$POS)]
lifespan_lead <- lifespan_lead[!duplicated(paste0(chr, ":", pos))]
message("Parental lifespan rows matched by chr:pos: ", nrow(lifespan_lead))

exp_df <- anm_lead[, .(SNP = paste0(CHR, ":", POS), effect_allele = toupper(Effect_Allele),
                        other_allele = toupper(Other_Allele), eaf = EAF, beta = Effect,
                        se = SE, pval = as.numeric(Pval), n = N)]
out_df <- lifespan_lead[, .(SNP = paste0(chr, ":", pos), effect_allele = toupper(a1),
                             other_allele = toupper(a0), eaf = freq1, beta = beta1,
                             se = se, pval = p, n = n)]

exp_dat <- format_data(as.data.frame(exp_df), type = "exposure", snp_col = "SNP",
                        beta_col = "beta", se_col = "se", eaf_col = "eaf",
                        effect_allele_col = "effect_allele", other_allele_col = "other_allele",
                        pval_col = "pval", samplesize_col = "n")
exp_dat$exposure <- "ANM"
out_dat <- format_data(as.data.frame(out_df), type = "outcome", snp_col = "SNP",
                        beta_col = "beta", se_col = "se", eaf_col = "eaf",
                        effect_allele_col = "effect_allele", other_allele_col = "other_allele",
                        pval_col = "pval", samplesize_col = "n")
out_dat$outcome <- "ParentalLifespan"

harm <- harmonise_data(exp_dat, out_dat, action = 2)
fwrite(harm, file.path(out_dir, "anm_to_lifespan_harmonised.tsv.gz"), sep = "\t")
message("Harmonised instruments retained: ", sum(harm$mr_keep))

res <- mr(harm, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median"))
res <- as.data.table(res)
fwrite(res, file.path(out_dir, "anm_to_lifespan_total_effect_mr.tsv"), sep = "\t")

het <- tryCatch(mr_heterogeneity(harm), error = function(e) NULL)
if (!is.null(het)) fwrite(as.data.table(het), file.path(out_dir, "anm_to_lifespan_heterogeneity.tsv"), sep = "\t")
pleio <- tryCatch(mr_pleiotropy_test(harm), error = function(e) NULL)
if (!is.null(pleio)) fwrite(as.data.table(pleio), file.path(out_dir, "anm_to_lifespan_egger_intercept.tsv"), sep = "\t")

message("\n=== ANM -> Parental Lifespan, total-effect MR ===")
print(res)
if (!is.null(het)) { message("\n=== Heterogeneity ==="); print(het) }
if (!is.null(pleio)) { message("\n=== Egger intercept ==="); print(pleio) }
