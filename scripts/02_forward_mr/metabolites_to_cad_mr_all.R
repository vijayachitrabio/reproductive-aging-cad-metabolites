#!/usr/bin/env Rscript
# Two-sample MR of all 56 NMR metabolites on coronary artery disease (CAD), with
# metabolome-wide and targeted FDR correction.

source(file.path(Sys.getenv("RAC_ROOT", "."), "config", "paths.R"))

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(ieugwasr)
  library(data.table)
})

if (Sys.getenv("OPENGWAS_JWT") == "") {
  stop("Set OPENGWAS_JWT in the environment before querying OpenGWAS.")
}

root <- RAC_ROOT
out_dir <- file.path(root, "results", "cad_all_metabolites")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

mqtl_file <- Sys.getenv(
  "MQTL_INSTRUMENTS",
  file.path(INSTRUMENT_DIR, "data/mqtl/mqtl_gwas_instruments.csv")
)
cad_id <- Sys.getenv("CAD_GWAS_ID", "ieu-a-7")
chunk_size <- as.integer(Sys.getenv("OPENGWAS_CHUNK_SIZE", "500"))

message("Loading mQTL instruments: ", mqtl_file)
mqtl <- fread(mqtl_file)
required <- c("metabolite", "SNP", "beta", "se", "eaf", "effect_allele", "other_allele", "pval")
missing_cols <- setdiff(required, names(mqtl))
if (length(missing_cols) > 0) {
  stop("Missing required mQTL columns: ", paste(missing_cols, collapse = ", "))
}

exposure <- format_data(
  as.data.frame(mqtl),
  type = "exposure",
  phenotype_col = "metabolite",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval"
)

snps <- unique(exposure$SNP)
message("Exposure instruments: ", nrow(exposure), " rows across ", length(unique(exposure$exposure)), " metabolites")
message("Unique SNPs requested from CAD GWAS: ", length(snps))

cache_file <- file.path(out_dir, paste0("cad_outcome_", cad_id, "_snps.rds"))
if (file.exists(cache_file)) {
  message("Using cached CAD outcome data: ", cache_file)
  outcome <- readRDS(cache_file)
} else {
  chunks <- split(snps, ceiling(seq_along(snps) / chunk_size))
  outcome_list <- vector("list", length(chunks))
  for (i in seq_along(chunks)) {
    message("Querying CAD chunk ", i, "/", length(chunks), " (", length(chunks[[i]]), " SNPs)")
    outcome_list[[i]] <- tryCatch(
      extract_outcome_data(snps = chunks[[i]], outcomes = cad_id),
      error = function(e) {
        warning("Chunk ", i, " failed: ", conditionMessage(e))
        NULL
      }
    )
    Sys.sleep(0.5)
  }
  outcome <- rbindlist(outcome_list, fill = TRUE)
  if (nrow(outcome) == 0) {
    stop("No CAD outcome data returned from OpenGWAS.")
  }
  saveRDS(outcome, cache_file)
}
outcome$outcome <- "Coronary Artery Disease"
fwrite(as.data.table(outcome), file.path(out_dir, "cad_outcome_extracted.tsv.gz"), sep = "\t")

message("Harmonising mQTL instruments against CAD")
harm <- harmonise_data(exposure_dat = exposure, outcome_dat = outcome)
fwrite(as.data.table(harm), file.path(out_dir, "metabolite_cad_harmonised.tsv.gz"), sep = "\t")

instrument_qc <- as.data.table(harm)[
  ,
  .(
    n_harmonised = .N,
    n_mr_keep = sum(mr_keep %in% TRUE),
    n_palindromic = sum(palindromic %in% TRUE, na.rm = TRUE),
    n_ambiguous = sum(ambiguous %in% TRUE, na.rm = TRUE)
  ),
  by = .(metabolite = exposure)
]
setorder(instrument_qc, metabolite)
fwrite(instrument_qc, file.path(out_dir, "metabolite_cad_instrument_qc.tsv"), sep = "\t")

message("Running CAD MR for all metabolites")
mr_res <- as.data.table(mr(
  harm,
  method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")
))

ivw <- mr_res[method == "Inverse variance weighted"]
ivw[, cad_global_fdr := p.adjust(pval, method = "fdr")]

anm_file <- file.path(root, "results", "metabolome_to_anm_mr_results.tsv")
if (file.exists(anm_file)) {
  anm <- fread(anm_file)
  anm_ivw <- anm[method == "Inverse variance weighted",
    .(exposure, anm_beta = b, anm_p = pval, anm_fdr = fdr, anm_nsnp = nsnp)
  ]
  ivw <- merge(ivw, anm_ivw, by = "exposure", all.x = TRUE)
  ivw[, anm_fdr_sig := anm_fdr < 0.05]
  ivw[, cad_anm_prioritized_fdr := NA_real_]
  ivw[anm_fdr_sig == TRUE, cad_anm_prioritized_fdr := p.adjust(pval, method = "fdr")]
  ivw[, shared_anm_cad_global_fdr := anm_fdr_sig == TRUE & cad_global_fdr < 0.05]
  ivw[, shared_anm_cad_targeted_fdr := anm_fdr_sig == TRUE & cad_anm_prioritized_fdr < 0.05]
}

mr_res <- merge(
  mr_res,
  ivw[, .(id.exposure, cad_global_fdr, cad_anm_prioritized_fdr, anm_beta, anm_p, anm_fdr, anm_nsnp)],
  by = "id.exposure",
  all.x = TRUE
)

pleio <- as.data.table(mr_pleiotropy_test(harm))
het <- as.data.table(mr_heterogeneity(harm))

setorder(mr_res, exposure, method)
setorder(ivw, cad_global_fdr, pval)
setorder(pleio, exposure)
setorder(het, exposure, method)

fwrite(mr_res, file.path(out_dir, "metabolome_to_cad_all_methods.tsv"), sep = "\t")
fwrite(ivw, file.path(out_dir, "metabolome_to_cad_ivw_global_and_targeted_fdr.tsv"), sep = "\t")
fwrite(pleio, file.path(out_dir, "metabolome_to_cad_egger_intercepts.tsv"), sep = "\t")
fwrite(het, file.path(out_dir, "metabolome_to_cad_heterogeneity.tsv"), sep = "\t")

message("Done.")
message("IVW metabolites tested: ", nrow(ivw))
message("CAD global FDR<0.05: ", sum(ivw$cad_global_fdr < 0.05, na.rm = TRUE))
if ("shared_anm_cad_targeted_fdr" %in% names(ivw)) {
  message("Shared ANM-FDR and CAD targeted-FDR<0.05: ", sum(ivw$shared_anm_cad_targeted_fdr, na.rm = TRUE))
}
message("Outputs written to: ", out_dir)
