#!/usr/bin/env Rscript
# Single-variant colocalization (coloc.abf) of prioritized metabolites with CAD.

source(file.path(Sys.getenv("RAC_ROOT", "."), "config", "paths.R"))

suppressPackageStartupMessages({
  library(data.table)
  library(coloc)
  library(ieugwasr)
})

if (Sys.getenv("OPENGWAS_JWT") == "") {
  stop("Set OPENGWAS_JWT in the environment before running this script.")
}

root <- RAC_ROOT
in_dir <- file.path(root, "results/coloc_metabolite_ovarian")
out_dir <- file.path(root, "results/coloc_metabolite_cad")
extract_dir <- file.path(out_dir, "locus_extracts")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)

cad_id <- "ieu-a-7"
cad_s <- 60801 / 184305 # case proportion for CARDIoGRAMplusC4D

# Load loci from the metabolite-ovarian trait colocalization
loci_file <- file.path(in_dir, "coloc_pilot_loci.tsv")
loci <- fread(loci_file)

results <- list()

for (i in seq_len(nrow(loci))) {
  loc <- loci[i]
  met_name <- loc$metabolite
  locus_id <- loc$locus_id
  
  cat(sprintf("\nProcessing %s: %s (chr%s:%s-%s)\n", met_name, locus_id, loc$chr, loc$start, loc$end))
  
  # Load cached metabolite data from the metabolite-ovarian trait colocalization
  met_cache <- file.path(in_dir, "locus_extracts", paste0(met_name, "_coloc_windows.tsv.gz"))
  if (!file.exists(met_cache)) {
    cat("  Missing metabolite cache for", met_name, "- skipping\n")
    next
  }
  met_all <- fread(met_cache)
  met_dt <- met_all[locus_id == loc$locus_id]
  if (nrow(met_dt) == 0) {
    cat("  No SNPs found for locus in metabolite cache - skipping\n")
    next
  }
  
  # Ensure matching keys
  met_dt[, variant_key := paste0("chr", chr, ":", pos)]
  
  # Query OpenGWAS for CAD regional data
  region_str <- sprintf("%s:%s-%s", loc$chr, loc$start, loc$end)
  cad_dt <- tryCatch({
    res <- associations(region_str, id = cad_id)
    as.data.table(res)
  }, error = function(e) {
    cat("  OpenGWAS query failed:", e$message, "\n")
    data.table()
  })
  
  if (nrow(cad_dt) == 0) {
    cat("  No CAD SNPs retrieved - skipping\n")
    next
  }
  
  cad_dt[, variant_key := paste0("chr", chr, ":", position)]
  
  # Merge
  common <- merge(
    met_dt[, .(variant_key, pos, beta_met = BETA, se_met = SE, maf_met = maf, n_met = N, p_met = p)],
    cad_dt[, .(variant_key, beta_cad = beta, se_cad = se, maf_cad = eaf, p_cad = p, n_cad = n)],
    by = "variant_key"
  )
  
  common <- unique(common, by = "variant_key")
  # Ensure MAF is actually MAF
  common[, maf_cad := pmin(maf_cad, 1 - maf_cad)]
  common <- common[!is.na(beta_met) & !is.na(beta_cad) & !is.na(se_met) & !is.na(se_cad) & se_met > 0 & se_cad > 0 & maf_met > 0 & maf_cad > 0]
  
  base <- data.table(
    metabolite = loc$metabolite,
    locus_id = loc$locus_id,
    outcome = "CAD",
    lead_snp = loc$lead_snp,
    chr = loc$chr,
    locus_start = loc$start,
    locus_end = loc$end,
    n_common_snps = nrow(common),
    min_metabolite_p = if(nrow(common)) min(common$p_met) else NA_real_,
    min_outcome_p = if(nrow(common)) min(common$p_cad) else NA_real_
  )
  
  if (nrow(common) < 50) {
    cat("  Skipping: too few overlapping SNPs (", nrow(common), ")\n")
    base[, status := "skipped_too_few_snps"]
    results[[i]] <- base
    next
  }
  
  cat(sprintf("  Found %d common SNPs. Running COLOC...\n", nrow(common)))
  
  d1 <- list(
    beta = common$beta_met,
    varbeta = common$se_met^2,
    MAF = common$maf_met,
    snp = common$variant_key,
    position = common$pos,
    N = median(common$n_met, na.rm=TRUE),
    type = "quant"
  )
  
  d2 <- list(
    beta = common$beta_cad,
    varbeta = common$se_cad^2,
    MAF = common$maf_cad,
    snp = common$variant_key,
    position = common$pos,
    N = median(common$n_cad, na.rm=TRUE),
    s = cad_s,
    type = "cc"
  )
  
  fit <- tryCatch(
    {
      quiet <- capture.output(ans <- coloc.abf(d1, d2))
      ans
    },
    error = function(e) e
  )
  
  if (inherits(fit, "error")) {
    cat("  COLOC error:", fit$message, "\n")
    base[, `:=`(status = "error", error = fit$message)]
    results[[i]] <- base
    next
  }
  
  # Save posteriors
  snp_post <- as.data.table(fit$results)
  snp_post[, `:=`(metabolite = met_name, locus_id = locus_id, outcome = "CAD")]
  fwrite(snp_post, file.path(extract_dir, paste0(locus_id, "_CAD_snp_posteriors.tsv.gz")), sep="\t")
  
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
  
  results[[i]] <- base
}

summary <- rbindlist(results, fill = TRUE)
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

setorder(summary, -PP.H4, min_outcome_p)

fwrite(summary, file.path(out_dir, "metabolite_cad_coloc_summary.tsv"), sep = "\t")
counts <- summary[, .N, by = .(outcome, coloc_support)][order(outcome, coloc_support)]
fwrite(counts, file.path(out_dir, "metabolite_cad_coloc_counts.tsv"), sep = "\t")

cat("\nDone! Wrote summary to", file.path(out_dir, "metabolite_cad_coloc_summary.tsv"), "\n")
print(summary[, .(metabolite, locus_id, n_common_snps, min_outcome_p, PP.H3, PP.H4, coloc_support)][order(-PP.H4)][1:10])
