#!/usr/bin/env python3
"""Compile the CAD/metabolite S-PrediXcan runs together with the ANM runs (same
8 GTEx tissues), applying a single pooled Bonferroni threshold across all
gene-tissue tests, then the QC filter n_snps_used >= 5 and pred_perf_r2 > 0
(or NA, which MASHR models commonly report) (Supplementary Table S8)."""
import glob
import os
import pandas as pd

ROOT = os.environ.get("RAC_ROOT", ".")
REPRO_GWAS_DIR = os.environ.get("REPRO_GWAS_DIR", os.path.join(ROOT, "external", "reproductive_gwas"))
NEW_DIR = f"{ROOT}/results/twas"
ANM_DIR = f"{REPRO_GWAS_DIR}/results/ovarian_twas_signal_atlas/per_tissue_results"  # S-PrediXcan results for ANM
OUT_DIR = f"{ROOT}/results/twas_compiled"
os.makedirs(OUT_DIR, exist_ok=True)
MIN_SNPS_USED = 5

TISSUES = {"Artery_Aorta", "Artery_Coronary", "Artery_Tibial", "Kidney_Cortex",
           "Pancreas", "Liver", "Adipose_Subcutaneous", "Adipose_Visceral_Omentum"}

all_results = []
for f in glob.glob(f"{NEW_DIR}/*.csv"):
    trait, tissue = os.path.basename(f).replace(".csv", "").split("__")
    df = pd.read_csv(f)
    df["Trait"], df["Tissue"] = trait, tissue
    all_results.append(df)

# Only pull ANM's 8 matching tissues (that dir also has Ovary/Uterus/Vagina
# runs from a different sub-analysis; keep this study's test universe
# restricted to the 8 tissues actually run for CAD/metabolites here).
for f in glob.glob(f"{ANM_DIR}/ANM__*.csv"):
    trait, tissue = os.path.basename(f).replace(".csv", "").split("__")
    if tissue not in TISSUES:
        continue
    df = pd.read_csv(f)
    df["Trait"], df["Tissue"] = trait, tissue
    all_results.append(df)

combined = pd.concat(all_results, ignore_index=True)
total_tests = len(combined)
bonferroni_p = 0.05 / total_tests
print(f"Total gene-tissue tests (ANM + CAD + 18 metabolites, 8 tissues, pooled): {total_tests}")
print(f"Pooled Bonferroni threshold: {bonferroni_p:.3e}")
print(f"Traits included: {sorted(combined['Trait'].unique())}")

sig = combined[combined["pvalue"] < bonferroni_p].sort_values("pvalue")
print(f"\nSignificant before QC: {len(sig)}")

qc_sig = sig.copy()
qc_sig = qc_sig[pd.to_numeric(qc_sig["n_snps_used"], errors="coerce") >= MIN_SNPS_USED]
pred_r2 = pd.to_numeric(qc_sig["pred_perf_r2"], errors="coerce")
qc_sig = qc_sig[pred_r2.isna() | (pred_r2 > 0)]
qc_sig = qc_sig.sort_values("pvalue")
print(f"Significant after QC (n_snps_used>=5, pred_perf_r2>0 or NA): {len(qc_sig)}")

combined.to_csv(f"{OUT_DIR}/twas_all_tests_pooled.csv", index=False)
sig.to_csv(f"{OUT_DIR}/twas_significant_before_qc.csv", index=False)
qc_sig.to_csv(f"{OUT_DIR}/twas_significant_after_qc.csv", index=False)

if len(qc_sig) > 0:
    print("\n--- Significant after QC ---")
    print(qc_sig[["Trait", "Tissue", "gene_name", "pvalue", "zscore", "n_snps_used"]].to_string(index=False))

# Cross-check: any gene significant (pre-QC) in BOTH ANM and CAD, or BOTH
# a metabolite and CAD, in the same tissue -- the gene-level analogue of
# this study's whole "shared genetic architecture" question.
metabolites = {"ApoA1", "Cholines", "DHA", "GlycA", "LA", "MUFA", "Omega_3",
               "Omega_6", "Omega_6_by_Omega_3", "PUFA", "SFA", "TG_by_PG",
               "Total_C", "Total_CE", "Total_FA", "Total_L", "Total_PL", "Total_TG"}
anm_genes = sig[sig["Trait"] == "ANM"][["gene_name", "Tissue"]].drop_duplicates()
cad_genes = sig[sig["Trait"] == "CAD"][["gene_name", "Tissue"]].drop_duplicates()
met_genes = sig[sig["Trait"].isin(metabolites)][["gene_name", "Tissue", "Trait"]].drop_duplicates()

anm_cad_overlap = anm_genes.merge(cad_genes, on=["gene_name", "Tissue"])
met_cad_overlap = met_genes.merge(cad_genes, on=["gene_name", "Tissue"])

print(f"\nGenes significant (pre-QC) in BOTH ANM and CAD, same tissue: {len(anm_cad_overlap)}")
if len(anm_cad_overlap) > 0:
    print(anm_cad_overlap.to_string(index=False))
print(f"\nGenes significant (pre-QC) in BOTH a prioritized metabolite and CAD, same tissue: {len(met_cad_overlap)}")
if len(met_cad_overlap) > 0:
    print(met_cad_overlap.to_string(index=False))
