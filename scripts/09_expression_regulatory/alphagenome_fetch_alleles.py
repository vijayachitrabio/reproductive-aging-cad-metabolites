#!/usr/bin/env python3
"""Resolve GRCh38 REF/ALT alleles for every variant in a locus of the locus variant table.

The table tags each variant to the nearest gene, and the causal gene in a lipid-TWAS
locus is not necessarily the tagged one (FADS_chr11_omega3 spans 6 genes).

Positions are GRCh38 (rs174528 -> chr11:61,776,027 matches Ensembl GRCh38;
hg19 would be 61,514,675). Every position is cross-checked against Ensembl and any offset is
reported per variant, so a coordinate/build mismatch cannot pass silently.

Output: results/alphagenome_locus/<locus>_variants_grch38.tsv
"""

from __future__ import annotations

import argparse
import os
import csv
import sys
import time
from collections import defaultdict
from pathlib import Path

import requests

ROOT = Path(os.environ.get("RAC_ROOT", Path(__file__).resolve().parents[2]))
# Locus variant table (one row per variant: locus, rsid, GRCh38 position, tagged gene).
S8 = Path(os.environ.get("ALPHAGENOME_LOCUS_TABLE", ROOT / "data" / "alphagenome" / "locus_variants.tsv"))
OUTDIR = ROOT / "results" / "alphagenome_locus"

ENSEMBL = "https://rest.ensembl.org/variation/human"
BATCH = 200
ASSEMBLY = "GRCh38"


def read_locus(locus: str) -> tuple[dict[str, int], dict[str, set[str]], dict[str, set[str]]]:
    pos: dict[str, int] = {}
    genes: dict[str, set[str]] = defaultdict(set)
    tissues: dict[str, set[str]] = defaultdict(set)
    with S8.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["locus"] != locus:
                continue
            rsid = row["snpId"]
            pos[rsid] = int(row["pos"])
            genes[rsid].add(row["geneSymbol"])
            tissues[rsid].add(row["tissue"])
    return pos, genes, tissues


def fetch_batch(rsids: list[str]) -> dict:
    for attempt in range(4):
        try:
            response = requests.post(
                ENSEMBL,
                json={"ids": rsids},
                headers={"Content-Type": "application/json", "Accept": "application/json"},
                timeout=60,
            )
            response.raise_for_status()
            return response.json()
        except Exception as exc:
            print(f"  batch of {len(rsids)} failed ({exc}); retry {attempt + 1}/4")
            if attempt == 3:
                return {}
            time.sleep(3 * (attempt + 1))
    return {}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--locus", default="FADS_chr11_omega3")
    args = parser.parse_args()

    s8_pos, s8_genes, s8_tissues = read_locus(args.locus)
    if not s8_pos:
        raise SystemExit(f"error: no rows for locus {args.locus} in {S8}")
    OUTDIR.mkdir(parents=True, exist_ok=True)
    out = OUTDIR / f"{args.locus}_variants_grch38.tsv"

    all_genes = sorted({g for v in s8_genes.values() for g in v})
    print(f"{len(s8_pos)} unique variants in {args.locus} across genes: {', '.join(all_genes)}")

    rsids = sorted(s8_pos, key=lambda r: s8_pos[r])
    records: dict[str, dict] = {}
    for start in range(0, len(rsids), BATCH):
        chunk = rsids[start : start + BATCH]
        print(f"fetching {start + 1}-{start + len(chunk)} of {len(rsids)}")
        records.update(fetch_batch(chunk))

    rows, unresolved, mismatch, var_classes = [], [], [], defaultdict(int)

    for rsid in rsids:
        record = records.get(rsid)
        if record is None:
            unresolved.append(rsid)
            continue
        var_classes[record.get("var_class", "?")] += 1
        mappings = [m for m in record.get("mappings", []) if m.get("assembly_name") == ASSEMBLY]
        if not mappings:
            unresolved.append(rsid)
            continue
        mapping = mappings[0]
        ens_pos = int(mapping["start"])
        offset = ens_pos - s8_pos[rsid]
        if offset != 0:
            mismatch.append((rsid, s8_pos[rsid], ens_pos))
        alleles = mapping["allele_string"].split("/")
        if len(alleles) < 2:
            unresolved.append(rsid)
            continue
        rows.append(
            {
                "locus": args.locus,
                "chrom": mapping["seq_region_name"],
                "pos_grch38": mapping["start"],
                "pos_in_s8": s8_pos[rsid],
                "pos_offset_vs_s8": offset,
                "snpId": rsid,
                "s8_genes": ",".join(sorted(s8_genes[rsid])),
                "s8_tissues": ",".join(sorted(s8_tissues[rsid])),
                "ref": alleles[0],
                "alts": ",".join(alleles[1:]),
                "n_alts": len(alleles) - 1,
                "var_class": record.get("var_class", ""),
                "most_severe_consequence": record.get("most_severe_consequence", ""),
                "ambiguity": record.get("ambiguity", ""),
                "assembly": ASSEMBLY,
                "source": "Ensembl REST /variation/human",
            }
        )

    # Ensembl REST returns transient 503s. fetch_batch() swallows them after its retries and
    # returns {}, so a partial outage silently looks like "these variants are unresolvable"
    # and would otherwise overwrite a complete variants file with a fraction of it. Observed
    # in practice: a re-fetch cut the FADS locus from 230 variants to 33 with no error. Refuse
    # to replace a good file unless the resolution rate is consistent with a few genuinely
    # unresolved variants, and never leave the file half-written.
    resolved_rate = len(rows) / len(rsids)
    if resolved_rate < 0.90:
        backup = out.with_suffix(".tsv.incomplete")
        print(
            f"error: only {len(rows)}/{len(rsids)} variants resolved "
            f"({resolved_rate:.1%}); refusing to overwrite {out}.",
            file=sys.stderr,
        )
        print(
            "This almost always means Ensembl REST was failing, not that the variants "
            "vanished. Re-run when the service is healthy.",
            file=sys.stderr,
        )
        if out.exists():
            backup.write_bytes(out.read_bytes())
            print(f"  existing file preserved at {backup}", file=sys.stderr)
        sys.exit(1)

    tmp = out.with_suffix(".tsv.tmp")
    with tmp.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    tmp.replace(out)

    print(f"wrote {len(rows)}/{len(rsids)} variants -> {out}")
    print(f"var classes: {dict(var_classes)}")
    print(f"position offsets vs S8: {len(mismatch)} (Ensembl GRCh38 used as authoritative)")
    for item in mismatch[:10]:
        print(f"  OFFSET {item[0]}: S8={item[1]} Ensembl38={item[2]} delta={item[2] - item[1]:+d}")
    print(f"unresolved: {len(unresolved)} {unresolved[:10]}")


if __name__ == "__main__":
    main()
