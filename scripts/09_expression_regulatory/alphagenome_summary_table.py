#!/usr/bin/env python3
"""Assemble Supplementary Table S7 (AlphaGenome predictions, FADS region) from the locus outputs.

Produces a TSV of the per-variant regulatory prediction and gene attribution
for the chromosome 11 FADS region, plus a per-gene summary.

Input : results/alphagenome_locus/FADS_chr11_omega3_gene_attribution.tsv
        results/alphagenome_locus/FADS_chr11_omega3_variants_grch38.tsv
        results/alphagenome_locus/FADS_chr11_omega3_variant_summary.tsv
Output: results/tables/Supplementary_Table_S21_alphagenome_fads_locus.tsv
        results/tables/Supplementary_Table_S21_alphagenome_fads_locus_per_gene.tsv
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import pandas as pd

ROOT = Path(os.environ.get("RAC_ROOT", Path(__file__).resolve().parents[2]))
LOCUS = "FADS_chr11_omega3"
SRC = ROOT / "results" / "alphagenome_locus"
OUT = ROOT / "results" / "tables"
THRESHOLD = 0.10


def main() -> None:
    attrib = pd.read_csv(SRC / f"{LOCUS}_gene_attribution.tsv", sep="\t")
    variants = pd.read_csv(SRC / f"{LOCUS}_variants_grch38.tsv", sep="\t")
    agnostic = pd.read_csv(SRC / f"{LOCUS}_variant_summary.tsv", sep="\t")

    # gene columns are everything except the annotation columns the scorer appended
    meta = {"snpId", "best_gene", "best_gene_max_abs_lfc", "margin_over_runner_up",
            "s8_genes", "pos_grch38", "ag_best_gene_in_s8"}
    genes = [c for c in attrib.columns if c not in meta]
    in_window = json.loads((SRC / f"{LOCUS}_genes_in_window.json").read_text())
    assert not (set(genes) - set(in_window)), "off-window gene leaked into the attribution table"
    assert not (set(in_window) - set(genes)) - {"LINC03144"}, "unexpected missing in-window gene"

    v = variants.drop_duplicates("snpId").set_index("snpId")
    out = pd.DataFrame(
        {
            "variant_id": attrib["snpId"],
            "chr": 11,
            "pos_grch38": attrib["snpId"].map(v["pos_grch38"]),
            "ref": attrib["snpId"].map(v["ref"]),
            "alts": attrib["snpId"].map(v["alts"]),
            "variant_class": attrib["snpId"].map(v["var_class"]),
            "most_severe_consequence": attrib["snpId"].map(v["most_severe_consequence"]),
            "pos_offset_vs_S8": attrib["snpId"].map(v["pos_offset_vs_s8"]),
            "S8_tagged_genes": attrib["s8_genes"],
            "S8_tissues": attrib["snpId"].map(v["s8_tissues"]),
            "predicted_top_gene": attrib["best_gene"],
            "predicted_top_gene_max_abs_lfc": attrib["best_gene_max_abs_lfc"].round(4),
            "margin_over_runner_up": attrib["margin_over_runner_up"].round(4),
            "top_gene_is_S8_tagged": attrib["ag_best_gene_in_s8"],
            "exceeds_0.10_threshold": attrib["best_gene_max_abs_lfc"] > THRESHOLD,
        }
    )
    for gene in genes:
        out[f"pred_abs_lfc_{gene}"] = attrib[gene].round(4)

    agnostic_max = (
        agnostic.groupby("snpId")["mean_quantile_score"].max().rename("max_mean_quantile_score")
    )
    out = out.merge(agnostic_max, left_on="variant_id", right_index=True, how="left")
    out["max_mean_quantile_score"] = out["max_mean_quantile_score"].round(4)
    out = out.sort_values("predicted_top_gene_max_abs_lfc", ascending=False)

    OUT.mkdir(parents=True, exist_ok=True)
    main_path = OUT / f"Supplementary_Table_S21_alphagenome_fads_locus.tsv"
    out.to_csv(main_path, sep="\t", index=False)

    per_gene = pd.DataFrame(
        {
            "gene": genes,
            "in_S8_tagged_set": [g in set(" ".join(attrib.s8_genes).split(",")) for g in genes],
            "n_variants_as_top_gene": [int((attrib["best_gene"] == g).sum()) for g in genes],
            "median_pred_abs_lfc": [round(float(attrib[g].median()), 4) for g in genes],
            "max_pred_abs_lfc": [round(float(attrib[g].max()), 4) for g in genes],
        }
    ).sort_values("max_pred_abs_lfc", ascending=False)
    per_gene_path = OUT / f"Supplementary_Table_S21_alphagenome_fads_locus_per_gene.tsv"
    per_gene.to_csv(per_gene_path, sep="\t", index=False)

    print(f"wrote {len(out)} rows x {out.shape[1]} cols -> {main_path}")
    print(f"wrote {len(per_gene)} rows -> {per_gene_path}")
    print(f"\nvariants exceeding {THRESHOLD}: {int(out['exceeds_0.10_threshold'].sum())} of {len(out)}")
    print("\ntop 5 by predicted effect:")
    print(out.head(5)[["variant_id", "predicted_top_gene", "predicted_top_gene_max_abs_lfc",
                       "margin_over_runner_up", "S8_tagged_genes"]].to_string(index=False))


if __name__ == "__main__":
    main()
