#!/usr/bin/env bash
# S-PrediXcan TWAS for CAD and the 18 prioritized metabolites across 8 GTEx v8
# MASHR tissues (4 vascular, 4 metabolic), the same tissues used for ANM.
# TWAS_TOOL_DIR holds MetaXcan, the MASHR weight databases and the SNP covariance file.
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/../.."

TWAS_TOOL_DIR="${TWAS_TOOL_DIR:-external/metaxcan}"
GWAS_DIR="data/twas_gwas_varid"
RESULTS_DIR="results/twas"
LOG_DIR="logs/twas"
mkdir -p "$RESULTS_DIR" "$LOG_DIR"

TISSUES=(Artery_Aorta Artery_Coronary Artery_Tibial Kidney_Cortex Pancreas Liver Adipose_Subcutaneous Adipose_Visceral_Omentum)
TRAITS=(CAD ApoA1 Cholines DHA GlycA LA MUFA Omega_3 Omega_6 Omega_6_by_Omega_3 PUFA SFA TG_by_PG Total_C Total_CE Total_FA Total_L Total_PL Total_TG)

for TRAIT in "${TRAITS[@]}"; do
  for TISSUE in "${TISSUES[@]}"; do
    OUT="${RESULTS_DIR}/${TRAIT}__${TISSUE}.csv"
    echo "=== ${TRAIT} x ${TISSUE} ==="
    python3 "${TWAS_TOOL_DIR}/MetaXcan/software/SPrediXcan.py" \
      --gwas_file "${GWAS_DIR}/${TRAIT}.tsv.gz" \
      --snp_column varID --effect_allele_column A1 --non_effect_allele_column A2 \
      --beta_column BETA --se_column SE --keep_non_rsid \
      --model_db_path "${TWAS_TOOL_DIR}/data/eqtl/mashr/mashr_${TISSUE}.db" \
      --model_db_snp_key varID \
      --covariance "${TWAS_TOOL_DIR}/data/gtex_v8_mashr_snp_covariance.txt.gz" \
      --output_file "$OUT" \
      --throw --overwrite > "${LOG_DIR}/${TRAIT}__${TISSUE}.log" 2>&1 || echo "  FAILED (see ${LOG_DIR}/${TRAIT}__${TISSUE}.log)"
  done
done
echo "All $((${#TRAITS[@]} * ${#TISSUES[@]})) trait x tissue TWAS runs complete."
