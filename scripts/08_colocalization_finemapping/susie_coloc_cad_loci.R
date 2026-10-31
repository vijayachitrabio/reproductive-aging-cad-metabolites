#!/usr/bin/env Rscript
# SuSiE colocalization fine-mapping at the CAD-side loci with PP.H4 >= 0.8 in the
# single-variant coloc.abf analysis (results/coloc_metabolite_cad/metabolite_cad_coloc_summary.tsv).
# The 16 metabolite-locus rows collapse to 3 genomic windows (metabolites at the
# same locus share one LD structure):
#   - APOE/APOC region, chr19:44,912,079-45,913,233 (Total_C, Total_CE,
#     Total_L, LA, Omega_6, Total_PL)
#   - APOA1/APOA5 cluster, chr11:116,148,917-117,148,917 (LA, MUFA, Omega_6,
#     PUFA, SFA, TG_by_PG, Total_FA, Total_TG)
#   - LDLR locus, chr19:10,698,502-11,698,502 (Total_C, Total_CE)
# CAD summary statistics are read from the local file
# data/dense_sumstats/CAD_ieu-a-7.tsv.gz (Table 2, Supplementary Table S4).

source(file.path(Sys.getenv("RAC_ROOT", "."), "config", "paths.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
  library(susieR)
})

root <- RAC_ROOT
atlas_root <- REPRO_GWAS_DIR
mqtl_dir <- file.path(INSTRUMENT_DIR, "data/mqtl")
full_gwas_dir <- file.path(mqtl_dir, "mqtl_full_gwas")

cad_file <- file.path(root, "data/dense_sumstats/CAD_ieu-a-7.tsv.gz")
cad_s <- 60801 / 184305 # case proportion, CARDIoGRAMplusC4D (ieu-a-7), as in the coloc.abf analysis

plink_prefix <- file.path(atlas_root, "data/lava_ref/g1000_eur")
plink_bin <- file.path(atlas_root, "tools/plink/plink")

out_dir <- file.path(root, "results", "susie_coloc")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

loci <- list(
  list(name = "APOE_APOC_chr19", chr = 19, start = 44912079, end = 45913233,
       metabolite = c("Total_C", "Total_CE", "Total_L", "LA", "Omega_6", "Total_PL")),
  list(name = "APOA1_APOA5_chr11", chr = 11, start = 116148917, end = 117148917,
       metabolite = c("LA", "MUFA", "Omega_6", "PUFA", "SFA", "TG_by_PG", "Total_FA", "Total_TG")),
  list(name = "LDLR_chr19", chr = 19, start = 10698502, end = 11698502,
       metabolite = c("Total_C", "Total_CE"))
)

message("Loading CAD dense summary statistics")
cad <- fread(cad_file, select = c("rsid", "effect_allele", "other_allele", "eaf", "beta", "se", "pval", "n"))
cad[, maf := pmin(eaf, 1 - eaf)]
cad <- cad[!is.na(beta) & !is.na(se) & se > 0 & maf > 0 & maf < 1]
setnames(cad, "rsid", "SNP")

get_ld_matrix <- function(snp_ids, tmp_prefix) {
  snp_file <- paste0(tmp_prefix, "_snps.txt")
  writeLines(snp_ids, snp_file)
  cmd <- sprintf(
    "%s --bfile %s --extract %s --make-bed --r square --out %s --silent",
    shQuote(plink_bin), shQuote(plink_prefix), shQuote(snp_file), shQuote(tmp_prefix)
  )
  system_result <- system(cmd, intern = TRUE, ignore.stderr = FALSE)
  ld_file <- paste0(tmp_prefix, ".ld")
  if (!file.exists(ld_file)) {
    message("PLINK output:\n", paste(system_result, collapse = "\n"))
    return(NULL)
  }
  ld <- as.matrix(fread(ld_file, header = FALSE))
  bim <- fread(paste0(tmp_prefix, ".bim"), header = FALSE)
  colnames(bim) <- c("chr", "SNP", "cm", "pos", "A1", "A2")
  rownames(ld) <- colnames(ld) <- bim$SNP
  ld
}

run_susie_coloc_one <- function(locus, met_name) {
  message("=== ", locus$name, " / ", met_name, " vs CAD ===")
  full_file <- file.path(full_gwas_dir, paste0(met_name, "_full_regenie.tsv.gz"))
  if (!file.exists(full_file)) {
    warning("Missing regenie file for ", met_name)
    return(data.table(locus = locus$name, metabolite = met_name, outcome = "CAD", status = "missing_metabolite_gwas"))
  }
  met <- fread(full_file, select = c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "N", "BETA", "SE", "LOG10P"))
  met <- met[CHROM == locus$chr & GENPOS >= locus$start & GENPOS <= locus$end]
  met[, `:=`(pos = as.integer(GENPOS), maf = pmin(A1FREQ, 1 - A1FREQ), p = 10^(-as.numeric(LOG10P)))]
  met <- met[!is.na(BETA) & !is.na(SE) & SE > 0 & maf > 0 & maf < 1]

  common <- merge(
    met[, .(SNP = ID, pos, beta_met = BETA, se_met = SE, maf_met = maf, n_met = N)],
    cad[, .(SNP, beta_cad = beta, se_cad = se, maf_cad = maf, n_cad = n)],
    by = "SNP"
  )
  common <- unique(common, by = "SNP")
  message("Common SNPs before LD: ", nrow(common))
  if (nrow(common) < 50) {
    warning("Too few common SNPs (", nrow(common), ") for ", locus$name, "/", met_name)
    return(data.table(locus = locus$name, metabolite = met_name, outcome = "CAD",
                       n_common_snps = nrow(common), status = "too_few_snps"))
  }

  tmp_prefix <- file.path(tempdir(), paste0("susie_cad_", locus$name, "_", met_name))
  ld <- get_ld_matrix(common$SNP, tmp_prefix)
  if (is.null(ld)) {
    warning("LD matrix extraction failed for ", locus$name, "/", met_name)
    return(data.table(locus = locus$name, metabolite = met_name, outcome = "CAD",
                       n_common_snps = nrow(common), status = "ld_extraction_failed"))
  }
  common <- common[SNP %in% rownames(ld)]
  common <- common[match(rownames(ld), SNP)]
  common <- common[!is.na(SNP)]

  n_met <- median(common$n_met, na.rm = TRUE)
  n_cad <- median(common$n_cad, na.rm = TRUE)

  d1 <- list(beta = common$beta_met, varbeta = common$se_met^2, N = n_met, MAF = common$maf_met,
             type = "quant", snp = common$SNP, LD = ld, position = common$pos)
  d2 <- list(beta = common$beta_cad, varbeta = common$se_cad^2, N = n_cad, MAF = common$maf_cad,
             type = "cc", s = cad_s, snp = common$SNP, LD = ld, position = common$pos)

  susie1 <- tryCatch(runsusie(d1, max_iter = 500), error = function(e) {message("SuSiE fit 1 failed: ", conditionMessage(e)); NULL})
  susie2 <- tryCatch(runsusie(d2, max_iter = 500), error = function(e) {message("SuSiE fit 2 failed: ", conditionMessage(e)); NULL})

  if (is.null(susie1) || is.null(susie2)) {
    abf <- tryCatch(coloc.abf(d1[c("beta","varbeta","N","MAF","type","snp")], d2[c("beta","varbeta","N","MAF","type","s","snp")]),
                     error = function(e) NULL)
    pph4_abf <- if (!is.null(abf)) abf$summary["PP.H4.abf"] else NA_real_
    return(data.table(locus = locus$name, metabolite = met_name, outcome = "CAD",
                       n_common_snps = nrow(common), status = "susie_fit_failed_fallback_abf",
                       coloc_abf_pph4 = pph4_abf))
  }

  res <- tryCatch(coloc.susie(susie1, susie2), error = function(e) {message("coloc.susie failed: ", conditionMessage(e)); NULL})
  if (is.null(res) || is.null(res$summary) || nrow(res$summary) == 0) {
    return(data.table(locus = locus$name, metabolite = met_name, outcome = "CAD",
                       n_common_snps = nrow(common), status = "no_credible_sets_or_coloc_failed",
                       n_credible_sets_trait1 = length(susie1$sets$cs), n_credible_sets_trait2 = length(susie2$sets$cs)))
  }

  summ <- as.data.table(res$summary)
  summ[, `:=`(locus = locus$name, metabolite = met_name, outcome = "CAD",
              n_common_snps = nrow(common), status = "ok")]
  summ
}

all_results <- list()
for (locus in loci) {
  for (met_name in locus$metabolite) {
    key <- paste(locus$name, met_name, sep = "__")
    all_results[[key]] <- tryCatch(
      run_susie_coloc_one(locus, met_name),
      error = function(e) {
        message("FAILED ", key, ": ", conditionMessage(e))
        data.table(locus = locus$name, metabolite = met_name, outcome = "CAD", status = paste("error:", conditionMessage(e)))
      }
    )
  }
}

final <- rbindlist(all_results, fill = TRUE)
fwrite(final, file.path(out_dir, "susie_coloc_cad_loci_results.tsv"), sep = "\t")
message("\nWrote ", file.path(out_dir, "susie_coloc_cad_loci_results.tsv"))
print(final)
