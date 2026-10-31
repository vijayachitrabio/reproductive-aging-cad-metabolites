#!/usr/bin/env Rscript
# Single-variant colocalization (coloc.abf) of prioritized metabolites with ANM and AMH.

source(file.path(Sys.getenv("RAC_ROOT", "."), "config", "paths.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
})

root <- RAC_ROOT
atlas_root <- REPRO_GWAS_DIR
mqtl_dir <- file.path(INSTRUMENT_DIR, "data/mqtl")
full_gwas_dir <- file.path(mqtl_dir, "mqtl_full_gwas")

instrument_file <- file.path(mqtl_dir, "mqtl_gwas_instruments.csv")
evidence_file <- file.path(root, "results/cad_all_metabolites/metabolic_paper_evidence_table.tsv")
anm_file <- file.path(atlas_root, "data/raw/ANM.txt.gz")
amh_file <- file.path(atlas_root, "data/raw/AMH.txt.gz")

out_dir <- file.path(root, "results/coloc_metabolite_ovarian")
extract_dir <- file.path(out_dir, "locus_extracts")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)

window_bp <- as.integer(Sys.getenv("COLOC_WINDOW_BP", "500000"))
top_loci_per_metabolite <- as.integer(Sys.getenv("COLOC_TOP_LOCI_PER_METABOLITE", "2"))
amh_n <- as.numeric(Sys.getenv("AMH_N", "9668"))
min_snps <- as.integer(Sys.getenv("COLOC_MIN_SNPS", "50"))

requested_metabolites <- Sys.getenv("COLOC_METABOLITES", "")

message("Loading prioritized metabolite evidence")
evidence <- fread(evidence_file)
if (requested_metabolites != "") {
  coloc_metabolites <- trimws(strsplit(requested_metabolites, ",")[[1]])
} else {
  coloc_metabolites <- evidence[
    evidence_tier != "not_prioritized" & cad_global_fdr < 0.05,
    unique(metabolite)
  ]
}

message("Coloc metabolites: ", paste(coloc_metabolites, collapse = ", "))

message("Loading instrument table")
inst <- fread(instrument_file)
inst <- inst[metabolite %in% coloc_metabolites]

select_lead_loci <- function(dt, top_n, window) {
  dt <- copy(dt)
  dt[, p_order := fifelse(is.na(pval) | pval <= 0, .Machine$double.xmin, pval)]
  setorder(dt, p_order, -F_stat)
  out <- list()
  for (i in seq_len(nrow(dt))) {
    candidate <- dt[i]
    if (length(out) > 0) {
      existing <- rbindlist(out)
      overlaps <- existing$chr == candidate$chr &
        abs(existing$pos - candidate$pos) <= (2 * window)
      if (any(overlaps)) next
    }
    candidate[, locus_rank := length(out) + 1L]
    out[[length(out) + 1L]] <- candidate
    if (length(out) >= top_n) break
  }
  rbindlist(out, fill = TRUE)
}

loci <- inst[, select_lead_loci(.SD, top_loci_per_metabolite, window_bp), by = metabolite]
loci[, `:=`(
  locus_id = paste0(metabolite, "_L", locus_rank, "_chr", chr, "_", pos),
  start = pmax(1L, as.integer(pos) - window_bp),
  end = as.integer(pos) + window_bp
)]
fwrite(
  loci[, .(metabolite, locus_id, locus_rank, lead_snp = SNP, chr, pos, start, end, lead_p = pval, lead_beta = beta, lead_F = F_stat)],
  file.path(out_dir, "coloc_pilot_loci.tsv"),
  sep = "\t"
)

message("Loading ANM and AMH outcome summary statistics")
anm <- fread(
  anm_file,
  select = c("SNP", "CHR", "POS", "Effect_Allele", "Other_Allele", "EAF", "Effect", "SE", "Pval", "N")
)
anm[, `:=`(
  chr = as.integer(CHR),
  pos = as.integer(POS),
  variant_key = paste0("chr", CHR, ":", POS),
  p = as.numeric(Pval),
  maf = pmin(EAF, 1 - EAF)
)]
anm <- anm[!is.na(chr) & !is.na(pos) & !is.na(Effect) & !is.na(SE) & SE > 0 & maf > 0 & maf < 1]
setkey(anm, chr, pos)

amh <- fread(
  amh_file,
  select = c("chromosome", "base_pair_location", "rs_id", "effect_allele", "other_allele", "beta", "standard_error", "effect_allele_frequency", "p_value")
)
amh[, `:=`(
  chr = as.integer(chromosome),
  pos = as.integer(base_pair_location),
  variant_key = paste0("chr", chromosome, ":", base_pair_location),
  p = as.numeric(p_value),
  maf = pmin(effect_allele_frequency, 1 - effect_allele_frequency)
)]
amh <- amh[!is.na(chr) & !is.na(pos) & !is.na(beta) & !is.na(standard_error) & standard_error > 0 & maf > 0 & maf < 1]
setkey(amh, chr, pos)

extract_metabolite_window <- function(met_name) {
  full_file <- file.path(full_gwas_dir, paste0(met_name, "_full_regenie.tsv.gz"))
  if (!file.exists(full_file)) {
    stop("Missing full metabolite GWAS: ", full_file)
  }
  cache <- file.path(extract_dir, paste0(met_name, "_coloc_windows.tsv.gz"))
  if (file.exists(cache)) {
    message("Using cached metabolite window extract: ", basename(cache))
    return(fread(cache))
  }

  met_loci <- loci[metabolite == met_name]
  chr_keep <- unique(met_loci$chr)

  message("Extracting coloc windows from ", basename(full_file))
  x <- fread(
    full_file,
    select = c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "N", "BETA", "SE", "LOG10P")
  )
  x[, `:=`(
    chr = as.integer(CHROM),
    pos = as.integer(GENPOS),
    variant_key = paste0("chr", CHROM, ":", GENPOS),
    p = 10^(-as.numeric(LOG10P)),
    maf = pmin(A1FREQ, 1 - A1FREQ)
  )]
  x <- x[chr %in% chr_keep]

  pieces <- lapply(seq_len(nrow(met_loci)), function(i) {
    loc <- met_loci[i]
    y <- x[chr == loc$chr & pos >= loc$start & pos <= loc$end]
    if (nrow(y) == 0) return(data.table())
    y[, `:=`(
      metabolite = loc$metabolite,
      locus_id = loc$locus_id,
      lead_snp = loc$SNP,
      window_start = loc$start,
      window_end = loc$end
    )]
    y
  })
  out <- rbindlist(pieces, fill = TRUE)
  out <- out[!is.na(BETA) & !is.na(SE) & SE > 0 & maf > 0 & maf < 1]
  fwrite(out, cache, sep = "\t")
  out
}

make_dataset <- function(dt, beta_col, se_col, maf_col, n_col, snp_col = "variant_key", pos_col = "position", n_override = NULL) {
  dt <- unique(dt, by = snp_col)
  n_value <- if (is.null(n_override)) median(dt[[n_col]], na.rm = TRUE) else n_override
  list(
    beta = dt[[beta_col]],
    varbeta = dt[[se_col]]^2,
    MAF = dt[[maf_col]],
    snp = dt[[snp_col]],
    position = dt[[pos_col]],
    N = n_value,
    type = "quant"
  )
}

run_coloc_pair <- function(met_dt, out_dt, outcome_name, loc) {
  common <- merge(
    met_dt[, .(variant_key, position = pos, beta_met = BETA, se_met = SE, maf_met = maf, n_met = N, p_met = p, rsid_met = ID)],
    out_dt,
    by = "variant_key"
  )
  common <- unique(common, by = "variant_key")

  base <- data.table(
    metabolite = loc$metabolite,
    locus_id = loc$locus_id,
    outcome = outcome_name,
    lead_snp = loc$SNP,
    chr = loc$chr,
    locus_start = loc$start,
    locus_end = loc$end,
    n_common_snps = nrow(common),
    min_metabolite_p = if (nrow(common)) min(common$p_met, na.rm = TRUE) else NA_real_
  )

  if (outcome_name == "ANM") {
    base[, min_outcome_p := if (nrow(common)) min(common$p, na.rm = TRUE) else NA_real_]
    if (nrow(common) < min_snps) return(base[, status := "skipped_too_few_snps"])
    d1 <- make_dataset(common, "beta_met", "se_met", "maf_met", "n_met")
    d2 <- make_dataset(common, "Effect", "SE", "maf", "N")
  } else {
    base[, min_outcome_p := if (nrow(common)) min(common$p, na.rm = TRUE) else NA_real_]
    if (nrow(common) < min_snps) return(base[, status := "skipped_too_few_snps"])
    d1 <- make_dataset(common, "beta_met", "se_met", "maf_met", "n_met")
    d2 <- make_dataset(common, "beta", "standard_error", "maf", "n_met", n_override = amh_n)
  }

  fit <- tryCatch(
    {
      quiet <- capture.output(ans <- coloc.abf(d1, d2))
      ans
    },
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    base[, `:=`(status = "error", error = fit$message)]
    return(base)
  }

  snp_post <- as.data.table(fit$results)
  snp_post[, `:=`(
    metabolite = loc$metabolite,
    locus_id = loc$locus_id,
    outcome = outcome_name
  )]
  fwrite(
    snp_post,
    file.path(extract_dir, paste0(loc$locus_id, "_", outcome_name, "_snp_posteriors.tsv.gz")),
    sep = "\t"
  )

  summ <- as.list(fit$summary)
  base[, `:=`(
    status = "ok",
    nsnps_coloc = as.numeric(summ[["nsnps"]]),
    PP.H0 = as.numeric(summ[["PP.H0.abf"]]),
    PP.H1 = as.numeric(summ[["PP.H1.abf"]]),
    PP.H2 = as.numeric(summ[["PP.H2.abf"]]),
    PP.H3 = as.numeric(summ[["PP.H3.abf"]]),
    PP.H4 = as.numeric(summ[["PP.H4.abf"]])
  )]
  fwrite(
    common,
    file.path(extract_dir, paste0(loc$locus_id, "_", outcome_name, "_common_variants.tsv.gz")),
    sep = "\t"
  )
  base
}

results <- list()
for (met in unique(loci$metabolite)) {
  met_extract <- extract_metabolite_window(met)
  for (i in seq_len(nrow(loci[metabolite == met]))) {
    loc <- loci[metabolite == met][i]
    met_dt <- met_extract[locus_id == loc$locus_id]

    anm_dt <- anm[chr == loc$chr & pos >= loc$start & pos <= loc$end]
    amh_dt <- amh[chr == loc$chr & pos >= loc$start & pos <= loc$end]

    results[[length(results) + 1L]] <- run_coloc_pair(met_dt, anm_dt, "ANM", loc)
    results[[length(results) + 1L]] <- run_coloc_pair(met_dt, amh_dt, "AMH", loc)
  }
}

summary <- rbindlist(results, fill = TRUE)
summary <- merge(
  summary,
  evidence[, .(metabolite, evidence_tier, anm_fdr, cad_global_fdr)],
  by = "metabolite",
  all.x = TRUE
)
summary[, coloc_support := fifelse(
  status == "ok" & PP.H4 >= 0.8,
  "strong_H4_ge_0.8",
  fifelse(
    status == "ok" & PP.H4 >= 0.5,
    "moderate_H4_ge_0.5",
    fifelse(
      status == "ok" & PP.H4 >= 0.25,
      "suggestive_H4_0.25_0.5",
      fifelse(status == "ok", "weak_H4_lt_0.25", status)
    )
  )
)]
setorder(summary, outcome, -PP.H4, min_outcome_p)

fwrite(summary, file.path(out_dir, "metabolite_ovarian_coloc_pilot_summary.tsv"), sep = "\t")
fwrite(summary[, .N, by = .(outcome, coloc_support)][order(outcome, coloc_support)], file.path(out_dir, "metabolite_ovarian_coloc_pilot_counts.tsv"), sep = "\t")

message("Done. Wrote coloc outputs to ", out_dir)
print(summary[, .(metabolite, locus_id, outcome, n_common_snps, min_metabolite_p, min_outcome_p, PP.H3, PP.H4, coloc_support)][order(outcome, -PP.H4)][1:30])
