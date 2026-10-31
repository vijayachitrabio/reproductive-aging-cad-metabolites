#!/usr/bin/env Rscript
# SuSiE colocalization (multiple-causal-variant fine-mapping) at two loci from the
# single-variant coloc.abf analysis:
#   - chr11 FADS region: Omega_3 vs ANM (coloc.abf PP.H4=0.707, moderate)
#   - chr19 APOE region: Total_C/Total_CE/Total_L/LA/Omega_6/Total_PL vs AMH
#     (coloc.abf PP.H4~0.458-0.459, suggestive)
#
# coloc.abf assumes a single causal variant per trait per locus; SuSiE allows
# multiple independent signals. The same full regenie metabolite GWAS as the
# coloc.abf analysis is used, so PP.H4 estimates are directly comparable.

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

anm_file <- file.path(atlas_root, "data/raw/ANM.txt.gz")
amh_file <- file.path(atlas_root, "data/raw/AMH.txt.gz")
plink_prefix <- file.path(atlas_root, "data/lava_ref/g1000_eur")
plink_bin <- file.path(atlas_root, "tools/plink/plink")

out_dir <- file.path(root, "results", "susie_coloc")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# The two loci (GRCh37 coordinates), from the coloc.abf analysis
loci <- list(
  list(name = "FADS_chr11_Omega3_ANM", chr = 11, start = 61088305, end = 62088305,
       metabolite = "Omega_3", outcome = "ANM"),
  list(name = "APOE_chr19_lipids_AMH", chr = 19, start = 44912079, end = 45913233,
       metabolite = c("Total_C", "Total_CE", "Total_L", "LA", "Omega_6", "Total_PL"),
       outcome = "AMH")
)

message("Loading ANM and AMH")
anm <- fread(anm_file, select = c("SNP", "CHR", "POS", "Effect_Allele", "Other_Allele", "EAF", "Effect", "SE", "Pval", "N"))
anm[, `:=`(chr = as.integer(CHR), pos = as.integer(POS), maf = pmin(EAF, 1 - EAF))]
anm <- anm[!is.na(chr) & !is.na(pos) & !is.na(Effect) & !is.na(SE) & SE > 0 & maf > 0 & maf < 1]

amh <- fread(amh_file, select = c("chromosome", "base_pair_location", "rs_id", "effect_allele", "other_allele", "beta", "standard_error", "effect_allele_frequency", "p_value"))
amh[, `:=`(chr = as.integer(chromosome), pos = as.integer(base_pair_location), maf = pmin(effect_allele_frequency, 1 - effect_allele_frequency))]
amh <- amh[!is.na(chr) & !is.na(pos) & !is.na(beta) & !is.na(standard_error) & standard_error > 0 & maf > 0 & maf < 1]

get_ld_matrix <- function(snp_ids, chr, tmp_prefix) {
  snp_file <- paste0(tmp_prefix, "_snps.txt")
  writeLines(snp_ids, snp_file)
  cmd <- sprintf(
    "%s --bfile %s --extract %s --make-bed --r square --out %s --silent",
    shQuote(plink_bin), shQuote(plink_prefix), shQuote(snp_file), shQuote(tmp_prefix)
  )
  system_result <- system(cmd, intern = TRUE, ignore.stderr = FALSE)
  ld_file <- paste0(tmp_prefix, ".ld")
  bim_file <- paste0(tmp_prefix, ".bim")
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

run_susie_coloc_one <- function(locus, met_name, outcome_dt, outcome_name) {
  message("=== ", locus$name, " / ", met_name, " vs ", outcome_name, " ===")
  full_file <- file.path(full_gwas_dir, paste0(met_name, "_full_regenie.tsv.gz"))
  if (!file.exists(full_file)) {
    warning("Missing regenie file for ", met_name)
    return(NULL)
  }
  met <- fread(full_file, select = c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "N", "BETA", "SE", "LOG10P"))
  met <- met[CHROM == locus$chr & GENPOS >= locus$start & GENPOS <= locus$end]
  met[, `:=`(pos = as.integer(GENPOS), maf = pmin(A1FREQ, 1 - A1FREQ), p = 10^(-as.numeric(LOG10P)))]
  met <- met[!is.na(BETA) & !is.na(SE) & SE > 0 & maf > 0 & maf < 1]

  out_sub <- outcome_dt[chr == locus$chr & pos >= locus$start & pos <= locus$end]

  common <- merge(
    met[, .(SNP = ID, pos, beta_met = BETA, se_met = SE, maf_met = maf, n_met = N)],
    out_sub[, .(pos, beta_out = if ("Effect" %in% names(out_sub)) Effect else beta,
                se_out = if ("SE" %in% names(out_sub)) SE else standard_error,
                n_out = if ("N" %in% names(out_sub)) N else NA_real_,
                maf_out = maf)],
    by = "pos"
  )
  common <- unique(common, by = "pos")
  message("Common SNPs before LD: ", nrow(common))
  if (nrow(common) < 50) {
    warning("Too few common SNPs (", nrow(common), ") for ", locus$name, "/", met_name)
    return(data.table(locus = locus$name, metabolite = met_name, outcome = outcome_name,
                       n_common_snps = nrow(common), status = "too_few_snps"))
  }

  tmp_prefix <- file.path(tempdir(), paste0("susie_", locus$name, "_", met_name))
  ld <- get_ld_matrix(common$SNP, locus$chr, tmp_prefix)
  if (is.null(ld)) {
    warning("LD matrix extraction failed for ", locus$name, "/", met_name)
    return(data.table(locus = locus$name, metabolite = met_name, outcome = outcome_name,
                       n_common_snps = nrow(common), status = "ld_extraction_failed"))
  }
  common <- common[SNP %in% rownames(ld)]
  common <- common[match(rownames(ld), SNP)]
  common <- common[!is.na(SNP)]

  n_met <- median(common$n_met, na.rm = TRUE)
  n_out <- if (all(is.na(common$n_out))) {
    if (outcome_name == "ANM") 201323 else 9668
  } else median(common$n_out, na.rm = TRUE)

  d1 <- list(beta = common$beta_met, varbeta = common$se_met^2, N = n_met, MAF = common$maf_met,
             type = "quant", snp = common$SNP, LD = ld, position = common$pos)
  d2 <- list(beta = common$beta_out, varbeta = common$se_out^2, N = n_out, MAF = common$maf_out,
             type = "quant", snp = common$SNP, LD = ld, position = common$pos)

  susie1 <- tryCatch(runsusie(d1, max_iter = 500), error = function(e) {message("SuSiE fit 1 failed: ", conditionMessage(e)); NULL})
  susie2 <- tryCatch(runsusie(d2, max_iter = 500), error = function(e) {message("SuSiE fit 2 failed: ", conditionMessage(e)); NULL})

  if (is.null(susie1) || is.null(susie2)) {
    # Fall back to standard coloc.abf on the same window for a direct comparison point
    abf <- tryCatch(coloc.abf(d1[c("beta","varbeta","N","MAF","type","snp")], d2[c("beta","varbeta","N","MAF","type","snp")]),
                     error = function(e) NULL)
    pph4_abf <- if (!is.null(abf)) abf$summary["PP.H4.abf"] else NA_real_
    return(data.table(locus = locus$name, metabolite = met_name, outcome = outcome_name,
                       n_common_snps = nrow(common), status = "susie_fit_failed_fallback_abf",
                       coloc_abf_pph4 = pph4_abf))
  }

  res <- tryCatch(coloc.susie(susie1, susie2), error = function(e) {message("coloc.susie failed: ", conditionMessage(e)); NULL})
  if (is.null(res) || is.null(res$summary) || nrow(res$summary) == 0) {
    return(data.table(locus = locus$name, metabolite = met_name, outcome = outcome_name,
                       n_common_snps = nrow(common), status = "no_credible_sets_or_coloc_failed",
                       n_credible_sets_trait1 = length(susie1$sets$cs), n_credible_sets_trait2 = length(susie2$sets$cs)))
  }

  summ <- as.data.table(res$summary)
  summ[, `:=`(locus = locus$name, metabolite = met_name, outcome = outcome_name,
              n_common_snps = nrow(common), status = "ok")]
  summ
}

all_results <- list()
for (locus in loci) {
  outcome_dt <- if (locus$outcome == "ANM") anm else amh
  for (met_name in locus$metabolite) {
    key <- paste(locus$name, met_name, sep = "__")
    all_results[[key]] <- tryCatch(
      run_susie_coloc_one(locus, met_name, outcome_dt, locus$outcome),
      error = function(e) {
        message("FAILED ", key, ": ", conditionMessage(e))
        data.table(locus = locus$name, metabolite = met_name, outcome = locus$outcome, status = paste("error:", conditionMessage(e)))
      }
    )
  }
}

final <- rbindlist(all_results, fill = TRUE)
fwrite(final, file.path(out_dir, "susie_coloc_fads_apoe_results.tsv"), sep = "\t")
message("\nWrote ", file.path(out_dir, "susie_coloc_fads_apoe_results.tsv"))
print(final)
