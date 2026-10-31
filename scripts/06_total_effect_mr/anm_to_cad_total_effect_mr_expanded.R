#!/usr/bin/env Rscript
# Sensitivity analysis: expanded instrument set for the ANM -> CAD total-effect MR.
# The primary analysis (anm_to_cad_total_effect_mr.R) mapped the 276 ANM lead loci
# (chr:pos, GRCh37) to rsIDs only through HapMap3 LD-score positions, so 52 of 276 loci were usable.
# Here the same 276 lead-locus instruments are mapped to rsIDs through the position-keyed metabolite GWAS
# extracts produced for the reverse ANM -> metabolite MR, and looked up in a local copy of the
# same CAD GWAS (OpenGWAS ieu-a-7). Estimators, harmonisation and outcome are unchanged.

source(file.path(Sys.getenv("RAC_ROOT", "."), "config", "paths.R"))

suppressPackageStartupMessages({ library(data.table); library(TwoSampleMR) })
root <- RAC_ROOT
out  <- file.path(root, "results/mediation_anm_cad_expanded_R2"); dir.create(out, showWarnings = FALSE)

expo <- fread(file.path(root, "results/reverse_mr_anm_to_metabolites/anm_reverse_exposure_instruments.tsv"))
stopifnot(nrow(expo) == 276)
ext <- rbindlist(lapply(list.files(file.path(root, "results/reverse_mr_anm_to_metabolites/extracted_metabolite_outcomes"), full.names = TRUE),
                        function(f) fread(f, select = c("SNP", "rsid"))))
map <- unique(ext[!is.na(rsid) & grepl("^rs", rsid)])
map <- map[, .SD[1], by = SNP]                       # one rsID per chr:pos
cat("ANM lead loci with an rsID through the metabolite extracts:", sum(expo$SNP %in% map$SNP), "of", nrow(expo), "\n")
expo <- merge(expo, map, by = "SNP")
expo[, chrpos := SNP][, SNP := rsid][, rsid := NULL]

cad <- fread(file.path(root, "data/dense_sumstats/CAD_ieu-a-7.tsv.gz"))
cad <- cad[rsid %in% expo$SNP]
cat("Found in the CAD GWAS (ieu-a-7):", nrow(cad), "\n")
outc <- format_data(as.data.frame(cad), type = "outcome", snp_col = "rsid", beta_col = "beta", se_col = "se",
                    eaf_col = "eaf", effect_allele_col = "effect_allele", other_allele_col = "other_allele",
                    pval_col = "pval", samplesize_col = "n")
outc$outcome <- "CAD (ieu-a-7)"; outc$id.outcome <- "ieu-a-7"
h <- harmonise_data(as.data.frame(expo), outc)
cat("Harmonised (mr_keep):", sum(h$mr_keep), "\n")
res <- mr(h, method_list = c("mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_weighted_mode"))
het <- mr_heterogeneity(h); ple <- mr_pleiotropy_test(h)
h$F <- (h$beta.exposure / h$se.exposure)^2
res$lci <- res$b - 1.96 * res$se; res$uci <- res$b + 1.96 * res$se
res$OR <- exp(res$b); res$OR_lci <- exp(res$lci); res$OR_uci <- exp(res$uci)
fwrite(res, file.path(out, "anm_to_cad_expanded_all_methods.tsv"), sep = "\t")
fwrite(het, file.path(out, "anm_to_cad_expanded_heterogeneity.tsv"), sep = "\t")
fwrite(ple, file.path(out, "anm_to_cad_expanded_egger_intercept.tsv"), sep = "\t")
fwrite(h, file.path(out, "anm_to_cad_expanded_harmonised.tsv.gz"), sep = "\t")
s <- data.table(step = c("ANM lead loci", "rsID via metabolite extracts", "present in CAD GWAS", "retained after harmonisation"),
                n = c(276, nrow(expo), nrow(cad), sum(h$mr_keep)))
fwrite(s, file.path(out, "instrument_attrition.tsv"), sep = "\t")
ivw <- res[res$method == "Inverse variance weighted", ]
mde <- (1.96 + 0.8416) * ivw$se
cat(sprintf("IVW beta %.5f (95%% CI %.4f to %.4f), SE %.5f, P %.3g; OR per year %.3f (%.3f-%.3f)\n",
            ivw$b, ivw$lci, ivw$uci, ivw$se, ivw$pval, ivw$OR, ivw$OR_lci, ivw$OR_uci))
cat(sprintf("Approximate minimum detectable OR per year at 80%% power (alpha 0.05): %.3f or %.3f\n", exp(mde), exp(-mde)))
cat(sprintf("Mean F %.1f, min F %.1f\n", mean(h$F[h$mr_keep]), min(h$F[h$mr_keep])))
print(res[, c("method", "nsnp", "b", "se", "pval", "OR", "OR_lci", "OR_uci")]); print(het); print(ple)
