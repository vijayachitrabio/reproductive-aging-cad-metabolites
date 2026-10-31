#!/usr/bin/env Rscript
# MR-PRESSO horizontal-pleiotropy/outlier analysis on the harmonised forward and
# reverse MR datasets (Supplementary Table S3a).

source(file.path(Sys.getenv("RAC_ROOT", "."), "config", "paths.R"))

suppressPackageStartupMessages({
  library(MRPRESSO)
  library(data.table)
  library(parallel)
})

root <- RAC_ROOT
out_dir <- file.path(root, "results/mr_presso")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# NbDistribution=1000 is a runtime/precision tradeoff: MRPRESSO's own outlier-test
# precision warning recommends NbDistribution >= n_snps/SignifThreshold (e.g. ~3900
# for a 196-SNP group at SignifThreshold=0.05), which benchmarked at ~90s/group;
# at 1000 draws the largest groups (>100 SNPs) have somewhat under-powered outlier
# detection but the global pleiotropy test (this script's primary use) is unaffected.
# This is noted in the output, not silently accepted.
nb_distribution <- as.integer(Sys.getenv("MRPRESSO_NB_DISTRIBUTION", "1000"))
min_snps <- 4L
n_cores <- max(1L, min(10L, detectCores() - 2L))

datasets <- list(
  list(
    label = "metabolite_to_ANM",
    file = file.path(root, "results/sensitivity_directionality/metabolite_anm_harmonised.tsv.gz"),
    group_col = "exposure"
  ),
  list(
    label = "metabolite_to_CAD",
    file = file.path(root, "results/cad_all_metabolites/metabolite_cad_harmonised.tsv.gz"),
    group_col = "exposure"
  ),
  list(
    label = "metabolite_to_AMH",
    file = file.path(root, "results/metabolome_to_amh/metabolite_amh_harmonised.tsv.gz"),
    group_col = "exposure"
  ),
  list(
    label = "ANM_to_metabolite",
    file = file.path(root, "results/reverse_mr_anm_to_metabolites/anm_to_prioritized_metabolites_harmonised.tsv.gz"),
    group_col = "outcome"
  ),
  list(
    label = "CAD_to_metabolite",
    file = file.path(root, "results/reverse_mr_cad_to_metabolites/cad_to_prioritized_metabolites_harmonised.tsv.gz"),
    group_col = "outcome"
  )
)

run_presso_one <- function(dat, label, group_name) {
  n <- nrow(dat)
  base <- data.table(
    dataset = label,
    group = group_name,
    n_snp = n
  )
  if (n < min_snps) {
    base[, `:=`(
      status = "skipped_too_few_snps",
      global_rss_obs = NA_real_,
      global_test_p = NA_real_,
      global_test_p_below_threshold = NA,
      n_outliers = NA_integer_,
      raw_beta = NA_real_,
      raw_se = NA_real_,
      outlier_corrected_beta = NA_real_,
      outlier_corrected_se = NA_real_,
      outlier_corrected_p = NA_real_,
      distortion_test_p = NA_real_
    )]
    return(base)
  }

  res <- tryCatch(
    mr_presso(
      BetaOutcome = "beta.outcome",
      BetaExposure = "beta.exposure",
      SdOutcome = "se.outcome",
      SdExposure = "se.exposure",
      OUTLIERtest = TRUE,
      DISTORTIONtest = TRUE,
      data = as.data.frame(dat),
      NbDistribution = nb_distribution,
      SignifThreshold = 0.05
    ),
    error = function(e) {
      message("  PRESSO failed for ", label, "/", group_name, ": ", conditionMessage(e))
      NULL
    }
  )

  if (is.null(res)) {
    base[, `:=`(
      status = "presso_error",
      global_rss_obs = NA_real_,
      global_test_p = NA_real_,
      global_test_p_below_threshold = NA,
      n_outliers = NA_integer_,
      raw_beta = NA_real_,
      raw_se = NA_real_,
      outlier_corrected_beta = NA_real_,
      outlier_corrected_se = NA_real_,
      outlier_corrected_p = NA_real_,
      distortion_test_p = NA_real_
    )]
    return(base)
  }

  main <- res[["Main MR results"]]
  raw_row <- main[main$`MR Analysis` == "Raw", ]
  oc_row <- main[main$`MR Analysis` == "Outlier-corrected", ]

  global_test <- res[["MR-PRESSO results"]][["Global Test"]]
  outlier_idx <- res[["MR-PRESSO results"]][["Distortion Test"]][["Outliers Indices"]]
  n_outliers <- if (is.null(outlier_idx) || identical(outlier_idx, "No significant outliers")) 0L else length(outlier_idx)
  distortion_p <- res[["MR-PRESSO results"]][["Distortion Test"]][["Pvalue"]]

  # MRPRESSO returns Global Test Pvalue as a character string like "<0.001"
  # (not NA) whenever the observed RSS exceeds every value in the simulated
  # null distribution, i.e. p is below the precision floor of NbDistribution.
  # This is the STRONGEST possible pleiotropy signal, not missing data --
  # as.numeric() on "<0.001" silently returns NA, which would misrepresent
  # exactly the groups with the clearest evidence of horizontal pleiotropy.
  raw_global_p <- global_test$Pvalue
  global_p_below_threshold <- is.character(raw_global_p) && grepl("^<", raw_global_p)
  global_test_p_numeric <- if (global_p_below_threshold) {
    as.numeric(sub("^<", "", raw_global_p))
  } else {
    as.numeric(raw_global_p)
  }

  base[, `:=`(
    status = "ok",
    global_rss_obs = as.numeric(global_test$RSSobs),
    global_test_p = global_test_p_numeric,
    global_test_p_below_threshold = global_p_below_threshold,
    n_outliers = n_outliers,
    raw_beta = as.numeric(raw_row$`Causal Estimate`),
    raw_se = as.numeric(raw_row$Sd),
    outlier_corrected_beta = if (nrow(oc_row) > 0) as.numeric(oc_row$`Causal Estimate`) else NA_real_,
    outlier_corrected_se = if (nrow(oc_row) > 0) as.numeric(oc_row$Sd) else NA_real_,
    outlier_corrected_p = if (nrow(oc_row) > 0) as.numeric(oc_row$`P-value`) else NA_real_,
    distortion_test_p = if (is.null(distortion_p)) NA_real_ else as.numeric(distortion_p)
  )]
  base
}

all_results <- list()

for (ds in datasets) {
  if (!file.exists(ds$file)) {
    message("Skipping ", ds$label, ": file not found (", ds$file, ")")
    next
  }
  message("Loading ", ds$label, " from ", ds$file)
  dat <- fread(ds$file)
  if (!"mr_keep" %in% names(dat)) {
    dat[, mr_keep := mr_keep.exposure %in% TRUE & mr_keep.outcome %in% TRUE]
  }
  dat <- dat[mr_keep %in% TRUE]
  group_col <- ds$group_col
  groups <- unique(dat[[group_col]])
  message("  ", length(groups), " groups (by ", group_col, "), running on ", n_cores, " cores")

  t_start <- Sys.time()
  group_list <- mclapply(seq_along(groups), function(i) {
    g <- groups[i]
    sub <- dat[get(group_col) == g]
    res <- run_presso_one(sub, ds$label, g)
    message("    [", ds$label, "] ", i, "/", length(groups), " done: ", g, " (n_snp=", nrow(sub), ")")
    res
  }, mc.cores = n_cores, mc.preschedule = FALSE)
  message("  ", ds$label, " finished in ", round(as.numeric(Sys.time() - t_start, units = "mins"), 1), " min")
  group_results <- rbindlist(group_list, fill = TRUE)
  all_results[[ds$label]] <- group_results
  fwrite(group_results, file.path(out_dir, paste0("mr_presso_", ds$label, ".tsv")), sep = "\t")
  message("  Wrote ", nrow(group_results), " rows for ", ds$label)
}

combined <- rbindlist(all_results, fill = TRUE)
setorder(combined, dataset, global_test_p)
fwrite(combined, file.path(out_dir, "mr_presso_all_datasets_summary.tsv"), sep = "\t")

message("Done. Combined MR-PRESSO summary written to ", file.path(out_dir, "mr_presso_all_datasets_summary.tsv"))
message("Datasets/groups with global test P<0.05 (evidence of horizontal pleiotropy/outliers): ",
        sum(combined$global_test_p < 0.05, na.rm = TRUE), " of ", sum(!is.na(combined$global_test_p)))
message("Groups with >=1 detected outlier: ", sum(combined$n_outliers > 0, na.rm = TRUE))
