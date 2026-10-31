#!/usr/bin/env Rscript
# Reverse MR: CAD -> 18 prioritized metabolites.

source(file.path(Sys.getenv("RAC_ROOT", "."), "config", "paths.R"))

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
})

if (Sys.getenv("OPENGWAS_JWT") == "") {
  stop("Set OPENGWAS_JWT in the environment before querying OpenGWAS.")
}

root <- RAC_ROOT
mqtl_full_dir <- file.path(INSTRUMENT_DIR, "data/mqtl/mqtl_full_gwas")

cad_id <- Sys.getenv("CAD_GWAS_ID", "ieu-a-7")
evidence_file <- file.path(root, "results/cad_all_metabolites/metabolic_paper_evidence_table.tsv")

out_dir <- file.path(root, "results/reverse_mr_cad_to_metabolites")
extract_dir <- file.path(out_dir, "extracted_metabolite_outcomes")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(extract_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading prioritized metabolites")
evidence <- fread(evidence_file)
priority <- evidence[evidence_tier != "not_prioritized"]
priority_metabolites <- unique(priority$metabolite)
message("Prioritized metabolites: ", paste(priority_metabolites, collapse = ", "))

message("Extracting genome-wide significant CAD instruments from ", cad_id)
inst_cache <- file.path(out_dir, "cad_instruments.tsv")
if (file.exists(inst_cache)) {
  message("Using cached CAD instruments: ", inst_cache)
  exp_cad <- fread(inst_cache)
} else {
  exp_cad <- as.data.table(extract_instruments(outcomes = cad_id))
  fwrite(exp_cad, inst_cache, sep = "\t")
}
exp_cad <- as.data.frame(exp_cad)
if (nrow(exp_cad) == 0) {
  stop("No CAD instruments returned from OpenGWAS for ", cad_id)
}
exp_cad$exposure <- "Coronary Artery Disease (CAD)"
exp_cad$id.exposure <- "CAD_reverse_lead_loci"

target_snps <- unique(exp_cad$SNP)
message("CAD lead instruments: ", length(target_snps))
fwrite(as.data.table(exp_cad), file.path(out_dir, "cad_reverse_exposure_instruments.tsv"), sep = "\t")

extract_metabolite_outcome <- function(metabolite) {
  cache <- file.path(extract_dir, paste0(metabolite, "_cad_lead_outcome_extract.tsv.gz"))
  if (file.exists(cache)) {
    message("Using cached outcome extract: ", basename(cache))
    return(fread(cache))
  }

  full_file <- file.path(mqtl_full_dir, paste0(metabolite, "_full_regenie.tsv.gz"))
  if (!file.exists(full_file)) {
    warning("Full metabolite GWAS not found for ", metabolite)
    return(data.table())
  }

  message("Extracting CAD lead SNP outcomes from ", basename(full_file))
  x <- fread(
    full_file,
    select = c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "N", "BETA", "SE", "LOG10P")
  )
  x[, SNP := paste0("chr", CHROM, ":", GENPOS)]
  x <- x[ID %in% target_snps | SNP %in% target_snps]
  if (nrow(x) == 0) return(data.table())
  x[, `:=`(
    metabolite = metabolite,
    pval = 10^(-as.numeric(LOG10P)),
    ALLELE1 = toupper(ALLELE1),
    ALLELE0 = toupper(ALLELE0)
  )]
  setorder(x, ID, -LOG10P, SE)
  x <- unique(x, by = "ID")
  out <- x[
    ,
    .(
      metabolite,
      SNP = ID,
      rsid = ID,
      effect_allele = ALLELE1,
      other_allele = ALLELE0,
      eaf = A1FREQ,
      beta = BETA,
      se = SE,
      pval,
      samplesize = N
    )
  ]
  fwrite(out, cache, sep = "\t")
  out
}

message("Preparing metabolite outcome extracts")
outcome_long <- rbindlist(lapply(priority_metabolites, extract_metabolite_outcome), fill = TRUE)
if (nrow(outcome_long) == 0) {
  stop("No metabolite outcome rows extracted for CAD lead SNPs.")
}
fwrite(outcome_long, file.path(out_dir, "prioritized_metabolite_outcome_extracts.tsv.gz"), sep = "\t")

out_met <- format_data(
  as.data.frame(outcome_long),
  type = "outcome",
  phenotype_col = "metabolite",
  snp_col = "SNP",
  beta_col = "beta",
  se_col = "se",
  eaf_col = "eaf",
  effect_allele_col = "effect_allele",
  other_allele_col = "other_allele",
  pval_col = "pval",
  samplesize_col = "samplesize"
)

message("Harmonising CAD instruments against metabolite outcomes")
harm <- as.data.table(harmonise_data(exp_cad, out_met))
fwrite(harm, file.path(out_dir, "cad_to_prioritized_metabolites_harmonised.tsv.gz"), sep = "\t")

message("Running reverse MR")
mr_res <- as.data.table(mr(
  harm,
  method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")
))

ivw <- mr_res[method == "Inverse variance weighted"]
ivw[, reverse_cad_to_metabolite_fdr := p.adjust(pval, method = "fdr")]
mr_res <- merge(
  mr_res,
  ivw[, .(id.outcome, reverse_cad_to_metabolite_fdr)],
  by = "id.outcome",
  all.x = TRUE
)

pleio <- as.data.table(mr_pleiotropy_test(harm))
het <- as.data.table(mr_heterogeneity(harm))

message("Running reverse Steiger directionality")
steiger_dat <- copy(harm[mr_keep == TRUE])
steiger_dat[, r.exposure := get_r_from_lor(
  lor = beta.exposure,
  af = eaf.exposure,
  ncase = 60801,
  ncontrol = 123504,
  prevalence = 0.06,
  model = "logit",
  correction = FALSE
)]
steiger_dat[, r.outcome := get_r_from_bsen(beta.outcome, se.outcome, samplesize.outcome)]
steiger <- as.data.table(directionality_test(as.data.frame(steiger_dat)))

qc <- harm[
  ,
  .(
    n_harmonised = .N,
    n_mr_keep = sum(mr_keep %in% TRUE),
    n_palindromic = sum(palindromic %in% TRUE, na.rm = TRUE),
    n_ambiguous = sum(ambiguous %in% TRUE, na.rm = TRUE)
  ),
  by = .(outcome, id.outcome)
]

summary <- merge(
  ivw[
    ,
    .(
      metabolite = outcome,
      reverse_beta = b,
      reverse_se = se,
      reverse_p = pval,
      reverse_fdr = reverse_cad_to_metabolite_fdr,
      reverse_nsnp = nsnp
    )
  ],
  priority[
    ,
    .(
      metabolite,
      evidence_tier,
      forward_metabolite_to_anm_beta = anm_beta,
      forward_metabolite_to_anm_fdr = anm_fdr,
      forward_metabolite_to_cad_beta = cad_beta,
      forward_metabolite_to_cad_fdr = cad_global_fdr
    )
  ],
  by = "metabolite",
  all.x = TRUE
)
summary <- merge(
  summary,
  steiger[
    ,
    .(
      metabolite = outcome,
      reverse_steiger_correct = correct_causal_direction,
      reverse_steiger_p = steiger_pval
    )
  ],
  by = "metabolite",
  all.x = TRUE
)
summary[, bidirectional_flag := reverse_fdr < 0.05 & reverse_steiger_correct == TRUE]
setorder(summary, reverse_fdr, reverse_p)

fwrite(mr_res, file.path(out_dir, "cad_to_prioritized_metabolites_all_methods.tsv"), sep = "\t")
fwrite(ivw, file.path(out_dir, "cad_to_prioritized_metabolites_ivw.tsv"), sep = "\t")
fwrite(pleio, file.path(out_dir, "cad_to_prioritized_metabolites_egger_intercepts.tsv"), sep = "\t")
fwrite(het, file.path(out_dir, "cad_to_prioritized_metabolites_heterogeneity.tsv"), sep = "\t")
fwrite(steiger, file.path(out_dir, "cad_to_prioritized_metabolites_steiger.tsv"), sep = "\t")
fwrite(qc, file.path(out_dir, "cad_to_prioritized_metabolites_qc.tsv"), sep = "\t")
fwrite(summary, file.path(out_dir, "cad_to_prioritized_metabolites_reverse_summary.tsv"), sep = "\t")

message("Done. Reverse MR outputs written to ", out_dir)
print(summary)
