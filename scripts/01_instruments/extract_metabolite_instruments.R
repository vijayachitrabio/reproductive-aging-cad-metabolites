#!/usr/bin/env Rscript
# Extract metabolite instruments from the full metabolite GWAS summary statistics.
#
# Method:
#   1. Genome-wide significance: p < 5e-8 (LOG10P >= 7.301)
#   2. Proximity clumping: within each 500-kb window keep the variant with the
#      smallest p-value (distance-based pruning without an LD reference panel)
#   3. F-statistic filter: F = beta^2 / se^2 > 10
#   4. Output formatted for TwoSampleMR

set.seed(42)

if (!require("data.table", quietly=TRUE)) install.packages("data.table", repos="https://cloud.r-project.org")
library(data.table)

source(file.path(Sys.getenv("RAC_ROOT", "."), "config", "paths.R"))
project_dir <- INSTRUMENT_DIR
gwas_dir    <- file.path(project_dir, "data", "mqtl", "mqtl_full_gwas")
out_dir     <- file.path(project_dir, "data", "mqtl")

P_THRESHOLD  <- 5e-8   # genome-wide significance
LOG10P_FLOOR <- -log10(P_THRESHOLD)   # 7.301
WINDOW_KB    <- 500    # clumping window (kb)
F_THRESHOLD  <- 10     # minimum F-statistic

gwas_files <- list.files(gwas_dir, pattern="full_regenie\\.tsv\\.gz$", full.names=TRUE)
cat(sprintf("Found %d full GWAS files.\n\n", length(gwas_files)))

all_instruments <- list()
summary_log     <- list()

for (f in gwas_files) {
  metabolite <- sub("_full_regenie\\.tsv\\.gz$", "", basename(f))
  cat(sprintf("[%s] Loading...\n", metabolite))

  # Read only needed columns; LOG10P column 13
  dat <- tryCatch(
    fread(f, select=c("CHROM","GENPOS","ID","ALLELE0","ALLELE1","A1FREQ","N","BETA","SE","LOG10P")),
    error=function(e) { cat("  ERROR:", conditionMessage(e), "\n"); NULL }
  )
  if (is.null(dat) || nrow(dat)==0) {
    summary_log[[metabolite]] <- data.table(metabolite=metabolite, n_gws=0, n_clumped=0, n_f_pass=0)
    next
  }

  # Remove rows with missing values
  dat <- dat[!is.na(BETA) & !is.na(SE) & SE > 0 & !is.na(LOG10P)]

  # GWS filter
  dat_gws <- dat[LOG10P >= LOG10P_FLOOR]
  n_gws <- nrow(dat_gws)
  cat(sprintf("  GWS (p<5e-8): %d hits\n", n_gws))

  if (n_gws == 0) {
    summary_log[[metabolite]] <- data.table(metabolite=metabolite, n_gws=0, n_clumped=0, n_f_pass=0)
    next
  }

  dat_gws[, CHROM := as.integer(CHROM)]

  # Proximity clumping: sort by significance, then greedily keep variants
  # that are >500kb from every already-kept variant on the same chromosome.
  # This keeps the most significant SNP per local locus while preserving
  # independent loci both upstream and downstream of the lead SNP.
  setorder(dat_gws, -LOG10P)
  kept <- list()
  for (chr in unique(dat_gws$CHROM)) {
    sub <- dat_gws[CHROM == chr]
    keep_idx <- logical(nrow(sub))
    kept_pos <- numeric()
    for (i in seq_len(nrow(sub))) {
      if (length(kept_pos) == 0 ||
          all(abs(sub$GENPOS[i] - kept_pos) > (WINDOW_KB * 1000))) {
        keep_idx[i] <- TRUE
        kept_pos <- c(kept_pos, sub$GENPOS[i])
      }
    }
    kept[[length(kept)+1]] <- sub[keep_idx]
  }
  clumped <- rbindlist(kept)
  setorder(clumped, CHROM, GENPOS)
  n_clumped <- nrow(clumped)
  cat(sprintf("  After 500kb clumping: %d independent loci\n", n_clumped))

  # F-statistic filter
  clumped[, F_stat := BETA^2 / SE^2]
  clumped_f <- clumped[F_stat > F_THRESHOLD]
  n_f_pass <- nrow(clumped_f)
  cat(sprintf("  After F>%d filter: %d instruments\n", F_THRESHOLD, n_f_pass))

  summary_log[[metabolite]] <- data.table(metabolite=metabolite, n_gws=n_gws, n_clumped=n_clumped, n_f_pass=n_f_pass)

  if (n_f_pass == 0) next

  # Format for TwoSampleMR
  clumped_f[, `:=`(
    metabolite     = metabolite,
    SNP            = ID,
    chr            = CHROM,
    pos            = GENPOS,
    effect_allele  = ALLELE1,  # ALLELE1 = effect allele in regenie
    other_allele   = ALLELE0,
    eaf            = A1FREQ,
    beta           = BETA,
    se             = SE,
    pval           = 10^(-LOG10P),
    n              = N
  )]

  all_instruments[[metabolite]] <- clumped_f[, .(metabolite, SNP, chr, pos,
                                                  effect_allele, other_allele,
                                                  eaf, beta, se, pval, F_stat, n)]
}

# Save
cat("\n--- Summary ---\n")
summary_dt <- rbindlist(summary_log, fill=TRUE)
print(summary_dt[order(-n_f_pass)])
fwrite(summary_dt, file.path(out_dir, "gwas_instrument_extraction_log.csv"))

if (length(all_instruments) > 0) {
  final <- rbindlist(all_instruments, fill=TRUE)

  cat(sprintf("\nTotal instruments: %d across %d metabolites\n",
              nrow(final), length(all_instruments)))
  cat(sprintf("Median instruments per metabolite: %.0f\n",
              median(summary_dt[n_f_pass>0, n_f_pass])))

  # Save
  fwrite(final, file.path(out_dir, "mqtl_gwas_instruments.csv"))
  saveRDS(final, file.path(out_dir, "mqtl_gwas_instruments.rds"))
  cat(sprintf("Saved to: %s\n", file.path(out_dir, "mqtl_gwas_instruments.csv")))
} else {
  cat("No instruments extracted.\n")
}

cat("\nDone.\n")
sessionInfo()
