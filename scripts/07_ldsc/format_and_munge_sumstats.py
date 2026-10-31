#!/usr/bin/env python3
"""
Format the dense GWAS files (CAD + 18 prioritized metabolites, in
data/dense_sumstats/) for LDSC, then munge them with munge_sumstats.py.

The metabolite files have "." (missing) in their per-row N column; the
manifest carries each trait's total N, which is used as a constant N
column before munging.
"""
import os
import gzip
import csv
import subprocess
import sys
from pathlib import Path

ROOT = Path(os.environ.get("RAC_ROOT", "."))
DENSE_DIR = ROOT / "data" / "dense_sumstats"
FORMATTED_DIR = ROOT / "data" / "formatted"
MUNGE_DIR = ROOT / "data" / "munge"
FORMATTED_DIR.mkdir(parents=True, exist_ok=True)
MUNGE_DIR.mkdir(parents=True, exist_ok=True)

REPRO_GWAS_DIR = Path(os.environ.get("REPRO_GWAS_DIR", "."))
MUNGE_PY = REPRO_GWAS_DIR / "tools" / "ldsc" / "munge_sumstats.py"
REF_SNPS = REPRO_GWAS_DIR / "data" / "ldsc_ref" / "w_hm3.snplist"

manifest = {}
with open(DENSE_DIR / "manifest.tsv") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        manifest[row["trait_name"]] = row

def format_file(trait_name, in_path, out_path, fallback_n=None):
    with gzip.open(in_path, "rt") as fin, gzip.open(out_path, "wt") as fout:
        header = fin.readline().strip().split("\t")
        idx = {c: i for i, c in enumerate(header)}
        fout.write("SNP\tA1\tA2\tEAF\tBETA\tSE\tP\tN\n")
        n_written = 0
        for line in fin:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < len(header):
                continue
            rsid = parts[idx["rsid"]]
            if not rsid.startswith("rs"):
                continue
            a1 = parts[idx["effect_allele"]]
            a2 = parts[idx["other_allele"]]
            eaf = parts[idx["eaf"]]
            beta = parts[idx["beta"]]
            se = parts[idx["se"]]
            pval = parts[idx["pval"]]
            n = parts[idx["n"]]
            if n == "." or n == "":
                if fallback_n is None:
                    continue
                n = str(fallback_n)
            if eaf == "." or beta == "." or se == "." or pval == ".":
                continue
            fout.write(f"{rsid}\t{a1}\t{a2}\t{eaf}\t{beta}\t{se}\t{pval}\t{n}\n")
            n_written += 1
    print(f"{trait_name}: wrote {n_written} rows to {out_path.name}")

def munge(trait_name, formatted_path):
    out_prefix = MUNGE_DIR / trait_name
    cmd = [
        sys.executable, str(MUNGE_PY),
        "--sumstats", str(formatted_path),
        "--out", str(out_prefix),
        "--merge-alleles", str(REF_SNPS),
        "--snp", "SNP", "--a1", "A1", "--a2", "A2", "--p", "P",
        "--signed-sumstats", "BETA,0",
        "--N-col", "N",
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    log_tail = "\n".join(result.stdout.strip().split("\n")[-15:])
    print(f"--- munge {trait_name} ---\n{log_tail}\n")
    return (MUNGE_DIR / f"{trait_name}.sumstats.gz").exists()

if __name__ == "__main__":
    results = {}
    for trait_name, row in manifest.items():
        file_path = DENSE_DIR / Path(row["file_path"]).name
        if not file_path.exists():
            print(f"MISSING: {trait_name} -> {file_path}")
            continue
        fallback_n = None if trait_name == "CAD" else row["n"]
        formatted_path = FORMATTED_DIR / f"{trait_name}.ldsc.tsv.gz"
        format_file(trait_name, file_path, formatted_path, fallback_n)
        ok = munge(trait_name, formatted_path)
        results[trait_name] = ok

    print("\n=== Summary ===")
    for t, ok in results.items():
        print(f"{t}: {'OK' if ok else 'FAILED'}")
