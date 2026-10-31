#!/usr/bin/env python3
"""Format CAD and the 18 prioritized metabolites' dense GWAS into GTEx varID-keyed
input for S-PrediXcan (rsid_to_varid.tsv lookup; output columns varID/A1/A2/BETA/SE)."""
import os
import csv
import gzip
from pathlib import Path

TWAS_DIR = Path(os.environ.get("TWAS_TOOL_DIR", Path(os.environ.get("RAC_ROOT", ".")) / "external" / "metaxcan"))
DENSE_DIR = Path(os.environ.get("RAC_ROOT", ".")) / "data/dense_sumstats"
OUT = Path(os.environ.get("RAC_ROOT", ".")) / "data/twas_gwas_varid"
OUT.mkdir(parents=True, exist_ok=True)

TRAITS = {
    "CAD": "CAD_ieu-a-7.tsv.gz",
    "ApoA1": "ApoA1_met-d-ApoA1.tsv.gz",
    "Cholines": "Cholines_met-d-Cholines.tsv.gz",
    "DHA": "DHA_met-d-DHA.tsv.gz",
    "GlycA": "GlycA_met-d-GlycA.tsv.gz",
    "LA": "LA_met-d-LA.tsv.gz",
    "MUFA": "MUFA_met-d-MUFA.tsv.gz",
    "Omega_3": "Omega_3_met-d-Omega_3.tsv.gz",
    "Omega_6": "Omega_6_met-d-Omega_6.tsv.gz",
    "Omega_6_by_Omega_3": "Omega_6_by_Omega_3_met-d-Omega_6_by_Omega_3.tsv.gz",
    "PUFA": "PUFA_met-d-PUFA.tsv.gz",
    "SFA": "SFA_met-d-SFA.tsv.gz",
    "TG_by_PG": "TG_by_PG_met-d-TG_by_PG.tsv.gz",
    "Total_C": "Total_C_met-d-Total_C.tsv.gz",
    "Total_CE": "Total_CE_met-d-Total_CE.tsv.gz",
    "Total_FA": "Total_FA_met-d-Total_FA.tsv.gz",
    "Total_L": "Total_L_met-d-Total_L.tsv.gz",
    "Total_PL": "Total_PL_met-d-Total_PL.tsv.gz",
    "Total_TG": "Total_TG_met-d-Total_TG.tsv.gz",
}


def load_mapping():
    m = {}
    with open(TWAS_DIR / "data" / "rsid_to_varid.tsv") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            m[row["rsid"]] = row["varID"]
    return m


def main():
    mapping = load_mapping()
    print(f"{len(mapping)} rsID->varID mappings loaded")

    for name, fname in TRAITS.items():
        path = DENSE_DIR / fname
        n_in = n_out = 0
        out_path = OUT / f"{name}.tsv.gz"
        with gzip.open(path, "rt") as fh, gzip.open(out_path, "wt", newline="") as out:
            reader = csv.DictReader(fh, delimiter="\t")
            w = csv.writer(out, delimiter="\t")
            w.writerow(["varID", "A1", "A2", "BETA", "SE"])
            for row in reader:
                n_in += 1
                rsid = row.get("rsid", "")
                varid = mapping.get(rsid)
                if not varid:
                    continue
                try:
                    beta = float(row["beta"])
                    se = float(row["se"])
                except (ValueError, KeyError):
                    continue
                if se <= 0:
                    continue
                a1 = row["effect_allele"].upper()
                a2 = row["other_allele"].upper()
                w.writerow([varid, a1, a2, beta, se])
                n_out += 1
        print(f"{name}: {n_out} / {n_in} SNPs mapped to varID -> {out_path.name}")


if __name__ == "__main__":
    main()
