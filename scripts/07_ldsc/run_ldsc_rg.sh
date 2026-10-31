#!/usr/bin/env bash
# LDSC genetic correlations (Figure 4, Table 3, Supplementary Table S5).
# The three ldsc.py calls below are those recorded in the headers of the LDSC logs
# (results/ldsc/rg_*.log). Run from the repository root after
# format_and_munge_sumstats.py; LDSC_DIR is a checkout of the original LDSC software
# (https://github.com/bulik/ldsc), run under its Python environment.
set -euo pipefail
LDSC_DIR="${LDSC_DIR:-external/ldsc}"
REPRO_GWAS_DIR="${REPRO_GWAS_DIR:-external/reproductive_gwas}"
REF="$REPRO_GWAS_DIR/data/ldsc_ref/eur_w_ld_chr/"
METS="data/munge/ApoA1.sumstats.gz,data/munge/Cholines.sumstats.gz,data/munge/DHA.sumstats.gz,data/munge/GlycA.sumstats.gz,data/munge/LA.sumstats.gz,data/munge/MUFA.sumstats.gz,data/munge/Omega_3.sumstats.gz,data/munge/Omega_6.sumstats.gz,data/munge/Omega_6_by_Omega_3.sumstats.gz,data/munge/PUFA.sumstats.gz,data/munge/SFA.sumstats.gz,data/munge/TG_by_PG.sumstats.gz,data/munge/Total_C.sumstats.gz,data/munge/Total_CE.sumstats.gz,data/munge/Total_FA.sumstats.gz,data/munge/Total_L.sumstats.gz,data/munge/Total_PL.sumstats.gz,data/munge/Total_TG.sumstats.gz"
mkdir -p results/ldsc

"$LDSC_DIR/ldsc.py" \
  --out results/ldsc/rg_ANM_vs_all \
  --rg "$REPRO_GWAS_DIR/data/munge/ANM.sumstats.gz,data/munge/CAD.sumstats.gz,$METS" \
  --ref-ld-chr "$REF" \
  --w-ld-chr "$REF"

"$LDSC_DIR/ldsc.py" \
  --out results/ldsc/rg_CAD_vs_metabolites \
  --rg "data/munge/CAD.sumstats.gz,$METS" \
  --ref-ld-chr "$REF" \
  --w-ld-chr "$REF"

"$LDSC_DIR/ldsc.py" \
  --out results/ldsc/rg_ParentalLifespan_vs_all \
  --rg "data/munge/ParentalLifespan.sumstats.gz,$REPRO_GWAS_DIR/data/munge/ANM.sumstats.gz,data/munge/CAD.sumstats.gz,$METS" \
  --ref-ld-chr "$REF" \
  --w-ld-chr "$REF"
