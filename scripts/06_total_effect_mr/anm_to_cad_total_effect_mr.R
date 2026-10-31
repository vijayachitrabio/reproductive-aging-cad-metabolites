#!/usr/bin/env Rscript
# Total-effect MR of ANM on CAD (Supplementary Table S22).
# Related estimates used in the mediation interpretation:
#   - ANM -> metabolite effects: results/reverse_mr_anm_to_metabolites/anm_to_prioritized_metabolites_reverse_summary.tsv
#   - metabolite -> CAD (MVMR-adjusted) effects: results/mvmr_axis/metabolic_axis_mvmr_results.tsv

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
atlas_root <- REPRO_GWAS_DIR
ld_dir <- file.path(atlas_root, "data/ldsc_ref/eur_w_ld_chr")

anm_file <- file.path(atlas_root, "data/raw/ANM.txt.gz")
anm_lead_file <- file.path(atlas_root, "results/phase4_loci/ANM.lead_loci.tsv")
cad_id <- Sys.getenv("CAD_GWAS_ID", "ieu-a-7")

out_dir <- file.path(root, "results/mediation_anm_cad")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading ANM lead loci")
lead <- fread(anm_lead_file)
message("ANM lead loci: ", nrow(lead))

message("Loading raw ANM summary statistics at lead loci")
anm <- fread(
  anm_file,
  select = c("SNP", "CHR", "POS", "Effect_Allele", "Other_Allele", "EAF", "Effect", "SE", "Pval", "N")
)
anm <- anm[SNP %in% lead$lead_snp]
anm[, `:=`(
  Effect_Allele = toupper(Effect_Allele),
  Other_Allele = toupper(Other_Allele),
  Pval = as.numeric(Pval)
)]
if (nrow(anm) == 0) stop("No ANM lead instruments found in raw ANM summary statistics.")

message("Mapping ANM lead loci chr:pos to HapMap3 rsIDs (needed for OpenGWAS CAD lookup)")
hm3_cache <- file.path(out_dir, "anm_lead_loci_rsid_map.tsv")
if (file.exists(hm3_cache)) {
  rsid_map <- fread(hm3_cache)
} else {
  chr_needed <- unique(anm$CHR)
  map_list <- vector("list", length(chr_needed))
  for (i in seq_along(chr_needed)) {
    chrom <- chr_needed[i]
    ld_path <- file.path(ld_dir, paste0(chrom, ".l2.ldscore.gz"))
    if (!file.exists(ld_path)) next
    ld <- fread(ld_path, select = c("CHR", "BP", "SNP"))
    ld <- ld[BP %in% anm[CHR == chrom]$POS]
    if (nrow(ld) > 0) map_list[[i]] <- ld[, .(CHR, BP, rsid = SNP)]
  }
  rsid_map <- rbindlist(map_list, fill = TRUE)
  fwrite(rsid_map, hm3_cache, sep = "\t")
}
message("Mapped ", nrow(rsid_map), " / ", nrow(anm), " ANM lead loci to HapMap3 rsIDs")

anm <- merge(anm, rsid_map, by.x = c("CHR", "POS"), by.y = c("CHR", "BP"))
if (nrow(anm) == 0) stop("No ANM lead loci mapped to rsIDs; cannot query OpenGWAS.")

exp_anm <- format_data(
  as.data.frame(anm),
  type = "exposure",
  snp_col = "rsid",
  beta_col = "Effect",
  se_col = "SE",
  eaf_col = "EAF",
  effect_allele_col = "Effect_Allele",
  other_allele_col = "Other_Allele",
  pval_col = "Pval",
  samplesize_col = "N"
)
exp_anm$exposure <- "Age at Natural Menopause (ANM)"
exp_anm$id.exposure <- "ANM_total_effect"
fwrite(as.data.table(exp_anm), file.path(out_dir, "anm_total_effect_exposure_instruments.tsv"), sep = "\t")

message("Extracting CAD (", cad_id, ") outcome data for ", length(unique(exp_anm$SNP)), " ANM lead-loci rsIDs")
cad_cache <- file.path(out_dir, paste0("cad_outcome_", cad_id, "_anm_lead_loci.rds"))
if (file.exists(cad_cache)) {
  out_cad <- readRDS(cad_cache)
} else {
  out_cad <- extract_outcome_data(snps = unique(exp_anm$SNP), outcomes = cad_id)
  saveRDS(out_cad, cad_cache)
}
if (is.null(out_cad) || nrow(out_cad) == 0) stop("No CAD outcome data returned from OpenGWAS.")
fwrite(as.data.table(out_cad), file.path(out_dir, "cad_outcome_at_anm_lead_loci.tsv.gz"), sep = "\t")

message("Harmonising ANM instruments against CAD outcome")
harm <- as.data.table(harmonise_data(exp_anm, out_cad))
fwrite(harm, file.path(out_dir, "anm_to_cad_harmonised.tsv.gz"), sep = "\t")

message("Running total-effect ANM -> CAD MR (", sum(harm$mr_keep), " SNPs after harmonisation)")
mr_res <- as.data.table(mr(harm, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median")))
pleio <- as.data.table(mr_pleiotropy_test(harm))
het <- as.data.table(mr_heterogeneity(harm))

fwrite(mr_res, file.path(out_dir, "anm_to_cad_total_effect_all_methods.tsv"), sep = "\t")
fwrite(pleio, file.path(out_dir, "anm_to_cad_total_effect_egger_intercept.tsv"), sep = "\t")
fwrite(het, file.path(out_dir, "anm_to_cad_total_effect_heterogeneity.tsv"), sep = "\t")

message("Done.")
print(mr_res)
print(pleio)
print(het)
