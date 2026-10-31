#!/usr/bin/env Rscript
# Parse the LDSC genetic-correlation logs into one table (Supplementary Table S5, Table 3).

source(file.path(Sys.getenv("RAC_ROOT", "."), "config", "paths.R"))

suppressPackageStartupMessages(library(data.table))

root <- RAC_ROOT
evidence <- fread(file.path(root, "results/final_integrated_evidence/metabolite_integrated_evidence.tsv"))

# read.table/fread choke on this log because the file paths themselves contain
# a space ("Mobile Documents") -- naive whitespace splitting misparses p1/p2.
# Parse by regex instead: pull the 10 trailing numeric fields directly, and
# identify trait2 by matching known trait basenames against the line text.
KNOWN_TRAITS <- c("CAD", "ApoA1", "Cholines", "DHA", "GlycA", "LA", "MUFA",
                   "Omega_3", "Omega_6_by_Omega_3", "Omega_6", "PUFA", "SFA",
                   "TG_by_PG", "Total_C", "Total_CE", "Total_FA", "Total_L",
                   "Total_PL", "Total_TG")

parse_log <- function(path, trait1_name) {
  lines <- readLines(path)
  start <- grep("^\\s*p1\\s+p2\\s+rg", lines)
  if (length(start) == 0) return(data.table())
  end <- grep("^Analysis finished", lines)[1] - 1
  body <- lines[(start + 1):end]
  body <- body[nzchar(trimws(body))]

  num_pat <- "(-?[0-9.]+(?:e-?[0-9]+)?|NA)\\s+"
  full_pat <- paste0("(", strrep(num_pat, 9), "(?:-?[0-9.]+(?:e-?[0-9]+)?|NA))\\s*$")

  rows <- lapply(body, function(ln) {
    m <- regmatches(ln, regexpr(full_pat, ln, perl = TRUE))
    if (length(m) == 0 || m == "") return(NULL)
    nums <- as.numeric(strsplit(trimws(m), "\\s+")[[1]])
    if (length(nums) != 10) return(NULL)
    prefix <- sub(full_pat, "", ln, perl = TRUE)
    # p1's own name always spuriously matches too (it's textually present in
    # the line); exclude it so a short trait2 name (DHA, LA, SFA) sharing a
    # tied nchar() with trait1_name isn't wrongly picked by which.max.
    candidates <- setdiff(KNOWN_TRAITS, trait1_name)
    hit <- candidates[sapply(candidates, function(t) grepl(paste0("/", t, "\\.sumstats\\.gz"), prefix))]
    trait2 <- if (length(hit) >= 1) hit[which.max(nchar(hit))] else NA_character_
    data.table(trait1 = trait1_name, trait2 = trait2,
               rg = nums[1], se = nums[2], z = nums[3], p = nums[4],
               h2_obs = nums[5], h2_obs_se = nums[6], h2_int = nums[7],
               h2_int_se = nums[8], gcov_int = nums[9], gcov_int_se = nums[10])
  })
  rbindlist(rows, fill = TRUE)
}

anm <- parse_log(file.path(root, "results/ldsc/rg_ANM_vs_all.log"), "ANM")
cad <- parse_log(file.path(root, "results/ldsc/rg_CAD_vs_metabolites.log"), "CAD")

anm_cad <- anm[trait2 == "CAD"]
anm_met <- anm[trait2 != "CAD"]
cad_met <- cad[trait2 != "CAD"]

anm_met[, pair_type := "ANM_vs_metabolite"]
cad_met[, pair_type := "CAD_vs_metabolite"]
anm_cad[, pair_type := "ANM_vs_CAD"]

all_pairs <- rbindlist(list(anm_cad, anm_met, cad_met), fill = TRUE)
all_pairs[, fdr_across_all_37 := p.adjust(p, method = "BH")]

all_pairs <- merge(all_pairs, evidence[, .(metabolite, evidence_tier, manuscript_priority)],
                    by.x = "trait2", by.y = "metabolite", all.x = TRUE)
all_pairs[trait2 == "CAD", `:=`(evidence_tier = NA, manuscript_priority = NA)]

setorder(all_pairs, p)
out <- all_pairs[, .(pair_type, trait1, trait2, rg, se, z, p, fdr_across_all_37,
                      h2_obs, h2_obs_se, gcov_int, gcov_int_se,
                      evidence_tier, manuscript_priority)]
fwrite(out, file.path(root, "results/ldsc/ldsc_rg_summary_all_pairs.tsv"), sep = "\t")

cat("=== ANM vs CAD ===\n")
print(out[pair_type == "ANM_vs_CAD"])
cat("\n=== Significant at FDR<0.05 (across all 37 tests) ===\n")
print(out[fdr_across_all_37 < 0.05, .(pair_type, trait1, trait2, rg, se, p, fdr_across_all_37, manuscript_priority)])
cat("\nTotal pairs tested:", nrow(out), "\n")
cat("Significant at nominal P<0.05:", sum(out$p < 0.05, na.rm=TRUE), "\n")
cat("Significant at FDR<0.05:", sum(out$fdr_across_all_37 < 0.05, na.rm=TRUE), "\n")
