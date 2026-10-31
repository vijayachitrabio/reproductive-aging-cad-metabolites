#!/usr/bin/env Rscript
# Sensitivity analysis for the primary metabolite instruments (Supplementary Table S1).
# The primary screen used 500-kb distance pruning (01_instruments/extract_metabolite_instruments.R).
# Here the instruments of the 18 prioritized metabolites are re-selected by LD-based clumping
# (r2 < 0.001, 10 Mb, PLINK v1.9 with the 1000 Genomes Phase 3 EUR reference panel) and IVW is
# re-estimated against CAD and ANM. Variants absent from the reference panel are dropped by the clumping step.

source(file.path(Sys.getenv("RAC_ROOT", "."), "config", "paths.R"))

suppressPackageStartupMessages({
  library(TwoSampleMR); library(ieugwasr); library(data.table)
})
atlas <- REPRO_GWAS_DIR

root <- RAC_ROOT
out_dir <- file.path(root, "results/ld_clump_sensitivity"); dir.create(out_dir, showWarnings = FALSE)
mqtl <- fread(file.path(INSTRUMENT_DIR, "data/mqtl/mqtl_gwas_instruments.csv"))
ev <- fread(file.path(root, "results/cad_all_metabolites/metabolic_paper_evidence_table.tsv"))
prior <- ev[evidence_tier != "not_prioritized"]

# ---- LD clumping per metabolite ----
clumped <- rbindlist(lapply(prior$metabolite, function(m) {
  x <- mqtl[metabolite == m, .(rsid = SNP, pval, id = m)]
  k <- ld_clump(x, clump_r2 = 0.001, clump_kb = 10000,
                plink_bin = file.path(atlas, "tools/plink/plink"),
                bfile = file.path(atlas, "data/lava_ref/g1000_eur"))
  message(sprintf("%-20s distance-pruned %3d -> LD-clumped %3d", m, nrow(x), nrow(k)))
  data.table(metabolite = m, SNP = k$rsid)
}))
fwrite(clumped, file.path(out_dir, "ld_clumped_instruments.tsv"), sep = "\t")

ivw <- function(h) {
  h <- h[h$mr_keep, ]
  if (nrow(h) < 2) return(data.table(b = NA_real_, se = NA_real_, p = NA_real_, nsnp = nrow(h)))
  r <- mr(h, method_list = "mr_ivw")
  data.table(b = r$b, se = r$se, p = r$pval, nsnp = r$nsnp)
}

# ---- CAD: reuse the primary harmonised CAD data (ieu-a-7), restricted to LD-clumped SNPs ----
cadh <- fread(file.path(root, "results/cad_all_metabolites/metabolite_cad_harmonised.tsv.gz"))
cad <- rbindlist(lapply(prior$metabolite, function(m) {
  h <- as.data.frame(cadh[exposure == m & SNP %in% clumped[metabolite == m, SNP]])
  cbind(metabolite = m, outcome = "CAD", ivw(h))
}))

# ---- ANM: same harmonisation as 02_forward_mr/metabolites_to_anm_mr.R, restricted to LD-clumped SNPs ----
sub <- merge(mqtl, clumped, by = c("metabolite", "SNP"))
sub[, SNP_chrpos := paste0("chr", chr, ":", pos)]
exp_dat <- format_data(as.data.frame(sub), type = "exposure", phenotype_col = "metabolite",
  snp_col = "SNP_chrpos", beta_col = "beta", se_col = "se", eaf_col = "eaf",
  effect_allele_col = "effect_allele", other_allele_col = "other_allele", pval_col = "pval")
anm_full <- fread(file.path(REPRO_GWAS_DIR, "data/raw/ANM.txt.gz"))
anm_full <- anm_full[SNP %in% unique(exp_dat$SNP)]
out_dat <- format_data(as.data.frame(anm_full), type = "outcome", snp_col = "SNP", beta_col = "Effect",
  se_col = "SE", eaf_col = "EAF", effect_allele_col = "Effect_Allele", other_allele_col = "Other_Allele", pval_col = "Pval")
anmh <- harmonise_data(exp_dat, out_dat)
anm <- rbindlist(lapply(prior$metabolite, function(m) cbind(metabolite = m, outcome = "ANM", ivw(anmh[anmh$exposure == m, ]))))

# ---- compare with primary (distance-pruned) estimates ----
res <- rbind(cad, anm)
pri <- rbind(prior[, .(metabolite, outcome = "CAD", primary_b = cad_beta, primary_p = cad_p, primary_nsnp = cad_nsnp)],
             prior[, .(metabolite, outcome = "ANM", primary_b = anm_beta, primary_p = anm_p, primary_nsnp = anm_nsnp)])
res <- merge(pri, res, by = c("metabolite", "outcome"))
setnames(res, c("b", "se", "p", "nsnp"), c("ldclump_b", "ldclump_se", "ldclump_p", "ldclump_nsnp"))
res[, same_direction := sign(primary_b) == sign(ldclump_b)]
res[, ldclump_fdr := p.adjust(ldclump_p, "fdr"), by = outcome]
fwrite(res, file.path(out_dir, "ld_clump_vs_distance_pruning_ivw.tsv"), sep = "\t")
print(res[, .(metabolite, outcome, primary_nsnp, ldclump_nsnp, primary_b = round(primary_b, 3),
              ldclump_b = round(ldclump_b, 3), ldclump_p = signif(ldclump_p, 2), same_direction)])
cat("\nDirection concordant:", sum(res$same_direction, na.rm = TRUE), "/", nrow(res),
    "| LD-clumped P<0.05:", sum(res$ldclump_p < 0.05, na.rm = TRUE), "| FDR<0.05:", sum(res$ldclump_fdr < 0.05, na.rm = TRUE), "\n")
