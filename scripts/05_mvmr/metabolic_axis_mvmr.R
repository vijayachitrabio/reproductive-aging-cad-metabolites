#!/usr/bin/env Rscript
# Multivariable MR of the four-trait metabolic axis (total cholesterol, total
# triglycerides, ApoA1, GlycA) on ANM and CAD (Figure 3, Supplementary Table S21).

source(file.path(Sys.getenv("RAC_ROOT", "."), "config", "paths.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(MVMR)
})

root <- RAC_ROOT
mqtl_dir <- file.path(INSTRUMENT_DIR, "data/mqtl")
full_gwas_dir <- file.path(mqtl_dir, "mqtl_full_gwas")
instrument_file <- file.path(mqtl_dir, "mqtl_gwas_instruments.csv")
anm_file <- file.path(REPRO_GWAS_DIR, "data/raw/ANM.txt.gz")
cad_file <- file.path(root, "results/cad_all_metabolites/cad_outcome_extracted.tsv.gz")

out_dir <- file.path(root, "results/mvmr_axis")
extract_dir <- file.path(out_dir, "extracted_metabolite_gwas")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)

modules <- list(
  lipid_inflammation_core = c("Total_C", "Total_TG", "ApoA1", "GlycA")
)

message("Loading mQTL instrument table")
inst <- fread(instrument_file)
inst[, snp_chrpos := paste0("chr", chr, ":", pos)]

needed_metabolites <- unique(unlist(modules))
missing_metabolites <- setdiff(needed_metabolites, unique(inst$metabolite))
if (length(missing_metabolites) > 0) {
  stop("Missing metabolites in instrument file: ", paste(missing_metabolites, collapse = ", "))
}

selected_inst <- inst[metabolite %in% needed_metabolites]
target_rsid <- unique(selected_inst$SNP)
target_chrpos <- unique(selected_inst$snp_chrpos)

is_palindromic <- function(a1, a2) {
  pair <- paste0(toupper(a1), toupper(a2))
  pair %in% c("AT", "TA", "CG", "GC")
}

align_to_outcome <- function(beta, ea, oa, out_ea, out_oa) {
  ea <- toupper(ea)
  oa <- toupper(oa)
  out_ea <- toupper(out_ea)
  out_oa <- toupper(out_oa)
  same <- ea == out_ea & oa == out_oa
  flip <- ea == out_oa & oa == out_ea
  out <- rep(NA_real_, length(beta))
  out[same] <- beta[same]
  out[flip] <- -beta[flip]
  out
}

extract_metabolite <- function(metabolite) {
  cache <- file.path(extract_dir, paste0(metabolite, "_union_instrument_extract.tsv.gz"))
  if (file.exists(cache)) {
    message("Using cached extract: ", basename(cache))
    return(fread(cache))
  }

  full_file <- file.path(full_gwas_dir, paste0(metabolite, "_full_regenie.tsv.gz"))
  if (!file.exists(full_file)) {
    stop("Full mQTL GWAS file not found for ", metabolite, ": ", full_file)
  }

  message("Extracting union instruments from ", basename(full_file))
  gwas <- fread(
    full_file,
    select = c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "N", "BETA", "SE", "LOG10P")
  )
  gwas[, snp_chrpos := paste0("chr", CHROM, ":", GENPOS)]
  gwas <- gwas[ID %in% target_rsid | snp_chrpos %in% target_chrpos]
  gwas[, metabolite := metabolite]
  setnames(
    gwas,
    old = c("ID", "ALLELE1", "ALLELE0", "A1FREQ", "N", "BETA", "SE", "LOG10P"),
    new = c("SNP", "effect_allele", "other_allele", "eaf", "n", "beta", "se", "log10p")
  )
  gwas <- gwas[
    ,
    .(metabolite, SNP, snp_chrpos, effect_allele, other_allele, eaf, n, beta, se, log10p)
  ]
  fwrite(gwas, cache, sep = "\t")
  gwas
}

message("Preparing exposure extracts")
exposure_long <- rbindlist(lapply(needed_metabolites, extract_metabolite), fill = TRUE)

deduplicate_exposures <- function(dt, key_col) {
  before <- nrow(dt)
  dup_keys <- dt[, .N, by = c("metabolite", key_col)][N > 1]
  if (nrow(dup_keys) > 0) {
    dt <- dt[!dup_keys, on = c("metabolite", key_col)]
  }
  list(data = dt, dropped = before - nrow(dt), duplicate_keys = nrow(dup_keys))
}

make_mvmr_input <- function(module_name, exposures, outcome_name, outcome, key_col,
                            outcome_beta_col, outcome_se_col, outcome_ea_col,
                            outcome_oa_col, outcome_eaf_col = NULL) {
  exp_dt <- exposure_long[metabolite %in% exposures]
  dedup <- deduplicate_exposures(copy(exp_dt), key_col)
  exp_dt <- dedup$data

  outcome_keep <- unique(
    outcome[
      ,
      c(key_col, outcome_beta_col, outcome_se_col, outcome_ea_col, outcome_oa_col, outcome_eaf_col),
      with = FALSE
    ],
    by = key_col
  )
  setnames(
    outcome_keep,
    old = c(outcome_beta_col, outcome_se_col, outcome_ea_col, outcome_oa_col),
    new = c("betaYG", "sebetaYG", "out_effect_allele", "out_other_allele")
  )
  if (!is.null(outcome_eaf_col)) setnames(outcome_keep, outcome_eaf_col, "out_eaf")

  joined <- merge(exp_dt, outcome_keep, by = key_col)
  joined[, beta_aligned := align_to_outcome(
    beta, effect_allele, other_allele, out_effect_allele, out_other_allele
  )]
  joined[, palindromic := is_palindromic(effect_allele, other_allele)]
  joined[, ambiguous_palindromic := palindromic & eaf > 0.42 & eaf < 0.58]
  if ("out_eaf" %in% names(joined)) {
    joined[, ambiguous_palindromic := ambiguous_palindromic | (palindromic & out_eaf > 0.42 & out_eaf < 0.58)]
  }

  aligned <- joined[!is.na(beta_aligned) & !ambiguous_palindromic]
  bx <- dcast(aligned, as.formula(paste(key_col, "+ betaYG + sebetaYG ~ metabolite")), value.var = "beta_aligned")
  sx <- dcast(aligned, as.formula(paste(key_col, "+ betaYG + sebetaYG ~ metabolite")), value.var = "se")
  setnames(sx, exposures, paste0("se_", exposures), skip_absent = TRUE)

  wide <- merge(bx, sx, by = c(key_col, "betaYG", "sebetaYG"))
  needed_cols <- c(exposures, paste0("se_", exposures))
  wide <- wide[complete.cases(wide[, ..needed_cols])]

  qc <- data.table(
    module = module_name,
    outcome = outcome_name,
    exposures = paste(exposures, collapse = ";"),
    n_union_instruments = uniqueN(selected_inst[metabolite %in% exposures][[key_col]]),
    n_outcome_overlap = uniqueN(joined[[key_col]]),
    n_after_allele_alignment = uniqueN(aligned[[key_col]]),
    n_complete_mvmr_snps = nrow(wide),
    duplicate_exposure_rows_dropped = dedup$dropped,
    duplicate_exposure_keys_dropped = dedup$duplicate_keys,
    ambiguous_palindromic_rows_dropped = nrow(joined[!is.na(beta_aligned) & ambiguous_palindromic]),
    allele_mismatch_rows_dropped = nrow(joined[is.na(beta_aligned)])
  )

  list(wide = wide, qc = qc)
}

run_one_mvmr <- function(module_name, exposures, outcome_name, input) {
  wide <- input$wide
  qc <- input$qc
  if (nrow(wide) <= length(exposures) + 2) {
    return(list(
      results = data.table(),
      strength = data.table(),
      pleiotropy = data.table(),
      qc = qc[, status := "skipped_too_few_complete_snps"]
    ))
  }

  bx <- as.matrix(wide[, ..exposures])
  sx <- as.matrix(wide[, paste0("se_", exposures), with = FALSE])
  r_input <- format_mvmr(
    BXGs = bx,
    BYG = wide$betaYG,
    seBXGs = sx,
    seBYG = wide$sebetaYG,
    RSID = wide[[1]]
  )

  message("Running MVMR: ", module_name, " -> ", outcome_name, " (", nrow(wide), " SNPs)")
  fit <- suppressWarnings(ivw_mvmr(r_input))
  strength <- suppressWarnings(strength_mvmr(r_input))
  pleio <- suppressWarnings(pleiotropy_mvmr(r_input))

  res <- as.data.table(fit, keep.rownames = "mvmr_exposure")
  res[, exposure := exposures]
  setnames(
    res,
    old = intersect(names(res), c("Estimate", "Std. Error", "t value", "Pr(>|t|)")),
    new = c("beta", "se", "t_stat", "pval")[seq_along(intersect(names(res), c("Estimate", "Std. Error", "t value", "Pr(>|t|)")))]
  )
  res[, `:=`(module = module_name, outcome = outcome_name, n_snps = nrow(wide))]
  setcolorder(res, c("module", "outcome", "exposure", "n_snps"))

  strength_dt <- as.data.table(strength)
  strength_dt <- melt(strength_dt, measure.vars = names(strength_dt), variable.name = "mvmr_exposure", value.name = "conditional_F")
  strength_dt[, exposure := exposures]
  strength_dt[, `:=`(module = module_name, outcome = outcome_name)]
  setcolorder(strength_dt, c("module", "outcome", "exposure", "conditional_F"))

  pleio_dt <- data.table(
    module = module_name,
    outcome = outcome_name,
    q_stat = as.numeric(pleio$Qstat),
    q_pval = as.numeric(pleio$Qpval)
  )

  cor_mat <- cor(bx, use = "pairwise.complete.obs")
  qc[, `:=`(
    status = "ok",
    exposure_beta_condition_number = kappa(cor_mat, exact = TRUE),
    max_abs_exposure_beta_correlation = max(abs(cor_mat[upper.tri(cor_mat)]), na.rm = TRUE)
  )]

  fwrite(wide, file.path(out_dir, paste0(module_name, "_", outcome_name, "_mvmr_input.tsv.gz")), sep = "\t")

  list(results = res, strength = strength_dt, pleiotropy = pleio_dt, qc = qc)
}

message("Loading outcomes")
anm <- fread(anm_file, select = c("SNP", "Effect", "SE", "Effect_Allele", "Other_Allele", "EAF", "Pval", "N"))
anm <- anm[SNP %in% target_chrpos]
setnames(anm, "SNP", "snp_chrpos")

cad <- fread(
  cad_file,
  select = c(
    "SNP", "beta.outcome", "se.outcome", "effect_allele.outcome",
    "other_allele.outcome", "eaf.outcome", "pval.outcome", "samplesize.outcome"
  )
)
cad <- cad[SNP %in% target_rsid]

all_results <- list()
all_strength <- list()
all_pleio <- list()
all_qc <- list()

for (module_name in names(modules)) {
  exposures <- modules[[module_name]]

  anm_input <- make_mvmr_input(
    module_name, exposures, "ANM", anm, "snp_chrpos",
    "Effect", "SE", "Effect_Allele", "Other_Allele", "EAF"
  )
  cad_input <- make_mvmr_input(
    module_name, exposures, "CAD", cad, "SNP",
    "beta.outcome", "se.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome"
  )

  anm_run <- run_one_mvmr(module_name, exposures, "ANM", anm_input)
  cad_run <- run_one_mvmr(module_name, exposures, "CAD", cad_input)

  all_results[[paste(module_name, "ANM")]] <- anm_run$results
  all_results[[paste(module_name, "CAD")]] <- cad_run$results
  all_strength[[paste(module_name, "ANM")]] <- anm_run$strength
  all_strength[[paste(module_name, "CAD")]] <- cad_run$strength
  all_pleio[[paste(module_name, "ANM")]] <- anm_run$pleiotropy
  all_pleio[[paste(module_name, "CAD")]] <- cad_run$pleiotropy
  all_qc[[paste(module_name, "ANM")]] <- anm_run$qc
  all_qc[[paste(module_name, "CAD")]] <- cad_run$qc
}

mvmr_results <- rbindlist(all_results, fill = TRUE)
mvmr_strength <- rbindlist(all_strength, fill = TRUE)
mvmr_pleio <- rbindlist(all_pleio, fill = TRUE)
mvmr_qc <- rbindlist(all_qc, fill = TRUE)

if (nrow(mvmr_results) > 0) {
  mvmr_results[, fdr_within_outcome := p.adjust(pval, method = "fdr"), by = outcome]
}

fwrite(mvmr_results, file.path(out_dir, "metabolic_axis_mvmr_results.tsv"), sep = "\t")
fwrite(mvmr_strength, file.path(out_dir, "metabolic_axis_mvmr_strength.tsv"), sep = "\t")
fwrite(mvmr_pleio, file.path(out_dir, "metabolic_axis_mvmr_pleiotropy.tsv"), sep = "\t")
fwrite(mvmr_qc, file.path(out_dir, "metabolic_axis_mvmr_qc.tsv"), sep = "\t")

message("Wrote MVMR outputs to ", out_dir)
print(mvmr_qc)
print(mvmr_results)
print(mvmr_strength)
print(mvmr_pleio)
