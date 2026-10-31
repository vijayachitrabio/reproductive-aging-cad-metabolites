#!/usr/bin/env python3
"""Score all variants in a locus with AlphaGenome, in chunks.

Input : results/alphagenome_locus/<locus>_variants_grch38.tsv
Output: results/alphagenome_locus/<locus>_variant_summary.tsv        gene-agnostic, variant x output_type
        results/alphagenome_locus/<locus>_gene_attribution.tsv       per-variant best gene
        results/alphagenome_locus/<locus>_<gene>_gene_summary.tsv    per-gene variant x output_type
        results/alphagenome_locus/<locus>_<gene>_{rnaseq,splice}.tsv only with --keep-tracks
        results/alphagenome_locus/<locus>_run_meta.json
        results/alphagenome_locus/cache/<locus>_chunk<i>.pkl         per-chunk AnnData cache

Chunking: tidy_scores() materialises every track for every variant, so each chunk is
scored, reduced to per-variant summaries, then discarded. Per-track tables are opt-in
via --keep-tracks because they are the bulk of the size.

Requires an AlphaGenome API key: export ALPHA_GENOME_API_KEY=<key>
Re-tidy from cache without re-querying: --reuse-cache
"""

from __future__ import annotations

import argparse
import gc
import hashlib
import json
import os
import pickle
import sys
from collections import defaultdict
from pathlib import Path

import pandas as pd

ROOT = Path(os.environ.get("RAC_ROOT", Path(__file__).resolve().parents[2]))
OUTDIR = ROOT / "results" / "alphagenome_locus"
CACHEDIR = OUTDIR / "cache"

SCORER_KEYS = (
    "ATAC",
    "DNASE",
    "CAGE",
    "PROCAP",
    "RNA_SEQ",
    "CHIP_HISTONE",
    "CHIP_TF",
    "SPLICE_SITES",
    "SPLICE_SITE_USAGE",
    "SPLICE_JUNCTIONS",
    "POLYADENYLATION",
    "CONTACT_MAPS",
)

AGG = {
    "n_tracks": ("raw_score", "size"),
    "n_genes": ("gene_name", "nunique"),
    "mean_quantile_score": ("quantile_score", "mean"),
    "max_abs_raw_score": ("raw_score", lambda s: s.abs().max()),
    "mean_abs_raw_score": ("raw_score", lambda s: s.abs().mean()),
    "frac_tracks_abs_gt_0p1": ("raw_score", lambda s: (s.abs() > 0.1).mean()),
}


def normalize_alleles(ref: str, alt: str) -> tuple[str, str]:
    return ("" if ref == "-" else ref), ("" if alt == "-" else alt)


def genes_in_window(chrom: str, start: int, end: int, cache_path: Path) -> list[str]:
    """Gene symbols physically overlapping the scored window, from Ensembl GRCh38.

    tidy_scores() emits a gene_name for every RNA_SEQ track, including a few that do not
    overlap the window. Attributing to an off-window gene would be meaningless, so the
    candidate set is pinned to the interval rather than trusted from the model output.
    """
    if cache_path.exists():
        return json.loads(cache_path.read_text())
    import requests

    session = requests.Session()
    session.headers.update({"Content-Type": "application/json"})
    response = session.get(
        f"https://rest.ensembl.org/overlap/region/human/{chrom}:{start}-{end}",
        params={"feature": "gene"},
        timeout=120,
    )
    response.raise_for_status()
    names = sorted({g["external_name"] for g in response.json() if g.get("external_name")})
    cache_path.write_text(json.dumps(names, indent=2) + "\n")
    print(f"  {len(names)} genes overlap the window (Ensembl) -> {cache_path}")
    return names


def load_variants(path: Path) -> list[dict]:
    import csv

    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def build_variants(rows: list[dict]):
    from alphagenome.data import genome

    variants, provenance = [], []
    for row in rows:
        for alt in row["alts"].split(","):
            ref_bases, alt_bases = normalize_alleles(row["ref"], alt)
            if ref_bases == "" and alt_bases == "":
                continue
            variants.append(
                genome.Variant(
                    chromosome=f"chr{row['chrom']}",
                    position=int(row["pos_grch38"]),
                    reference_bases=ref_bases,
                    alternate_bases=alt_bases,
                )
            )
            provenance.append(
                {
                    "snpId": row["snpId"],
                    "pos_grch38": int(row["pos_grch38"]),
                    "ref": row["ref"],
                    "alt": alt,
                    "var_class": row["var_class"],
                    "most_severe_consequence": row["most_severe_consequence"],
                    "pos_offset_vs_s8": int(row["pos_offset_vs_s8"]),
                    "s8_genes": row["s8_genes"],
                    "s8_tissues": row["s8_tissues"],
                }
            )
    return variants, provenance


def reduce_chunk(tidy: pd.DataFrame, chunk_index: int, keep_tracks: set[str] | None):
    """Collapse one chunk's tidy to per-variant summaries; return them, never the tidy."""
    agnostic = (
        tidy.groupby(["snpId", "alt_allele", "var_class", "output_type"], dropna=False)
        .agg(**AGG)
        .reset_index()
    )
    per_gene = {}
    for gene, frame in tidy.dropna(subset=["gene_name"]).groupby("gene_name"):
        per_gene[gene] = (
            frame.groupby(["snpId", "output_type"], dropna=False).agg(**AGG).reset_index()
        )
    gene_wide = {}
    rnaseq = tidy[tidy["output_type"] == "RNA_SEQ"]
    for gene, frame in rnaseq.dropna(subset=["gene_name"]).groupby("gene_name"):
        gene_wide[gene] = frame.groupby("snpId")["raw_score"].apply(lambda s: s.abs().max())

    tracks = {}
    if keep_tracks:
        wanted = keep_tracks
        if keep_tracks == {"all"}:
            wanted = set(tidy["gene_name"].dropna().unique())
        for gene in wanted:
            frame = tidy[tidy["gene_name"] == gene]
            if not frame.empty:
                tracks[gene] = frame
    return agnostic, per_gene, gene_wide, tracks


def chunk_fingerprint(interval, chunk: list, scorer_keys) -> str:
    """Stable id for the exact scoring request behind one cache chunk.

    Caches were previously named by chunk index alone, so a cache written at one
    --chunk-size could be read back at another: the variant list per chunk depends on
    chunk size, so a mismatch produced either a missing chunk001 or duplicated
    variants. Naming by a fingerprint of the request itself removes the ambiguity.

    The fingerprint covers everything that changes the answer -- the scored interval,
    the scorer set, and the exact variants -- so a cache is reused only when it was
    computed from the identical request.
    """
    digest = hashlib.sha256()
    digest.update(f"interval={interval.chromosome}:{interval.start}-{interval.end}\n".encode())
    digest.update(("scorers=" + ",".join(sorted(scorer_keys)) + "\n").encode())
    for variant in sorted(
        chunk,
        key=lambda v: (v.chromosome, v.position, v.reference_bases, v.alternate_bases),
    ):
        digest.update(
            f"{variant.chromosome}:{variant.position}:"
            f"{variant.reference_bases}>{variant.alternate_bases}\n".encode()
        )
    return digest.hexdigest()[:16]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--locus", default="FADS_chr11_omega3")
    parser.add_argument("--api-key-env", default="ALPHA_GENOME_API_KEY")
    parser.add_argument("--reuse-cache", action="store_true")
    parser.add_argument("--chunk-size", type=int, default=250)
    parser.add_argument(
        "--keep-tracks",
        default="",
        help="comma-separated genes to keep full per-track tables for, or 'all'. Default: none.",
    )
    args = parser.parse_args()

    OUTDIR.mkdir(parents=True, exist_ok=True)
    CACHEDIR.mkdir(parents=True, exist_ok=True)
    variant_list = OUTDIR / f"{args.locus}_variants_grch38.tsv"
    if not variant_list.exists():
        sys.exit(f"error: {variant_list} not found — run fetch_locus_alleles.py first")

    keep_tracks = {g for g in args.keep_tracks.split(",") if g} or None
    api_key = os.environ.get(args.api_key_env)

    from alphagenome.data import genome
    from alphagenome.models import dna_client, variant_scorers

    rows = load_variants(variant_list)
    variants, provenance = build_variants(rows)
    chrom = rows[0]["chrom"]
    key = {}
    for p in provenance:
        ref_bases, alt_bases = normalize_alleles(p["ref"], p["alt"])
        key[f"chr{chrom}:{p['pos_grch38']}:{ref_bases}>{alt_bases}"] = p

    lo = min(v.position for v in variants)
    hi = max(v.position + len(v.reference_bases) for v in variants)
    # The API requires the "chr" prefix on the scored interval as well as on individual
    # variants. build_variants() already applies it to variants, but the interval was
    # being built from the bare chromosome name, which the API rejects with
    # INVALID_ARGUMENT: Unrecognized chromosome '11' for organism. The Ensembl overlap
    # call further down still uses the bare name, which is correct for REST.
    iv = genome.Interval(chromosome=f"chr{chrom}", start=lo, end=hi).resize(
        dna_client.SEQUENCE_LENGTH_1MB
    )
    if hi - lo > iv.width:
        sys.exit(
            f"error: locus span {(hi - lo) / 1e6:.2f} Mb exceeds the 1 Mb model window; "
            "split the locus and run separately"
        )
    print(
        f"{len(provenance)} observations from {len({p['snpId'] for p in provenance})} unique SNPs; "
        f"span {chrom}:{lo}-{hi} ({(hi - lo) / 1000:.1f} kb)"
    )
    print(f"interval: {iv}")

    scorers = [
        variant_scorers.RECOMMENDED_VARIANT_SCORERS[k]
        for k in SCORER_KEYS
        if k in variant_scorers.RECOMMENDED_VARIANT_SCORERS
    ]
    chunks = [variants[i : i + args.chunk_size] for i in range(0, len(variants), args.chunk_size)]
    print(f"scorers: {len(scorers)} | chunks: {len(chunks)} of <= {args.chunk_size}")

    model = None
    if not args.reuse_cache:
        if not api_key:
            sys.exit(f"error: no API key in ${args.api_key_env}. Export it or pass --reuse-cache.")
        model = dna_client.create(api_key)

    agnostic_all, per_gene_all, gene_wide_frames = [], defaultdict(list), []
    track_frames: dict[str, list] = defaultdict(list)
    scored_genes: set[str] = set()
    model_versions: set[str] = set()

    for index, chunk in enumerate(chunks):
        fingerprint = chunk_fingerprint(iv, chunk, SCORER_KEYS)
        cache = CACHEDIR / f"{args.locus}_chunk{index:03d}_{fingerprint}.pkl"
        legacy = CACHEDIR / f"{args.locus}_chunk{index:03d}.pkl"
        if args.reuse_cache and cache.exists():
            raw = pickle.loads(cache.read_bytes())
            print(f"chunk {index + 1}/{len(chunks)}: from cache {cache.name}")
        elif args.reuse_cache and legacy.exists():
            # Pre-fingerprint cache. Its name carries only the chunk index, so it cannot be
            # checked against the interval or scorers it was actually computed under; fall
            # back to a variant-count check and say so, rather than silently trusting it.
            raw = pickle.loads(legacy.read_bytes())
            n_cached = len(raw) if isinstance(raw, list) else 1
            if n_cached != len(chunk):
                sys.exit(
                    f"error: {legacy.name} holds {n_cached} variants but chunk "
                    f"{index + 1}/{len(chunks)} at --chunk-size {args.chunk_size} "
                    f"contains {len(chunk)}. This cache was written with a different "
                    f"--chunk-size; rerun with the original size, or without "
                    f"--reuse-cache to rescore."
                )
            print(
                f"chunk {index + 1}/{len(chunks)}: from legacy cache {legacy.name} "
                f"(variant count matches; interval and scorers NOT verified -- this "
                f"file predates cache fingerprinting)"
            )
        else:
            if model is None:
                model = dna_client.create(api_key)
            print(f"chunk {index + 1}/{len(chunks)}: querying {len(chunk)} variants")
            raw = model.score_variants(
                intervals=iv,
                variants=chunk,
                variant_scorers=scorers,
                organism=dna_client.Organism.HOMO_SAPIENS,
                progress_bar=True,
            )
            cache.write_bytes(pickle.dumps(raw))

        tidy = variant_scorers.tidy_scores(raw)
        if tidy is None:
            sys.exit("error: tidy_scores returned None")
        tidy["vid"] = tidy["variant_id"].map(str)
        tidy["snpId"] = tidy["vid"].map(lambda v: key.get(v, {}).get("snpId", ""))
        tidy["alt_allele"] = tidy["vid"].map(lambda v: key.get(v, {}).get("alt", ""))
        tidy["var_class"] = tidy["vid"].map(lambda v: key.get(v, {}).get("var_class", ""))
        tidy["pos_offset_vs_s8"] = tidy["vid"].map(
            lambda v: key.get(v, {}).get("pos_offset_vs_s8", "")
        )
        if (tidy["snpId"] == "").any():
            bad = tidy.loc[tidy["snpId"] == "", "vid"].unique()[:3]
            sys.exit(f"error: unmapped variant_id values: {list(bad)}")

        agnostic, per_gene, gene_wide, tracks = reduce_chunk(tidy, index, keep_tracks)
        agnostic_all.append(agnostic)
        for gene, frame in per_gene.items():
            per_gene_all[gene].append(frame)
            scored_genes.add(gene)
        for gene, series in gene_wide.items():
            gene_wide_frames.append(series.rename(gene))
        for gene, frame in tracks.items():
            track_frames[gene].append(frame)
        print(f"  reduced {len(tidy)} rows -> {len(agnostic)} agnostic, {len(per_gene)} genes")

        del tidy, raw, per_gene, gene_wide, tracks
        gc.collect()

    summary = pd.concat(agnostic_all, ignore_index=True)
    meta_by_snp = {}
    for p in provenance:
        meta_by_snp.setdefault(p["snpId"], p)
    summary["s8_genes"] = summary["snpId"].map(lambda s: meta_by_snp[s]["s8_genes"])
    summary["s8_tissues"] = summary["snpId"].map(lambda s: meta_by_snp[s]["s8_tissues"])
    summary["pos_grch38"] = summary["snpId"].map(lambda s: meta_by_snp[s]["pos_grch38"])
    summary = summary.sort_values(["output_type", "mean_quantile_score"], ascending=[True, False])
    summary_path = OUTDIR / f"{args.locus}_variant_summary.tsv"
    summary.to_csv(summary_path, sep="\t", index=False)
    print(f"wrote {len(summary)} rows -> {summary_path}")

    gene_wide = pd.concat(gene_wide_frames, axis=1).fillna(0.0)
    # A gene scored across several chunks contributes one column per chunk. Collapse
    # those to a single max per (variant, gene) before anything selects columns by name;
    # otherwise the in-window filter below keeps every copy and the written table ends up
    # with one column per gene per chunk. Symptom seen with a multi-chunk locus: the APOE
    # attribution table had 3463 columns for 61 distinct gene names.
    gene_wide = gene_wide.T.groupby(level=0).max().T
    gene_wide = gene_wide.reindex(sorted({p["snpId"] for p in provenance})).fillna(0.0)
    # Attribute across every gene physically in the scored window, not just the S8-tagged set, so
    # genes S8 never tags (FADS2, SCGB1D1, ...) can still carry an effect. S8 tags are flagged,
    # not filtered. Off-window gene_name values from tidy_scores are dropped.
    in_window = set(genes_in_window(chrom, iv.start, iv.end, OUTDIR / f"{args.locus}_genes_in_window.json"))
    gene_cols = [c for c in gene_wide.columns if c in in_window]
    dropped = sorted(set(gene_wide.columns) - in_window)
    if dropped:
        print(f"  dropped {len(dropped)} off-window gene_name values: {dropped}")
    if not gene_cols:
        sys.exit("error: no RNA_SEQ gene overlaps the scored window; check the interval")
    # Restrict the written table to in-window genes too, so the written file contains
    # only genes that are actually in the locus (off-window ENSG* labels otherwise persist).
    gene_wide = gene_wide[sorted(gene_cols)]
    unscored = sorted(in_window - set(gene_cols))
    if unscored:
        print(f"  {len(unscored)} in-window genes have no RNA_SEQ rows: {unscored}")
    gene_wide["best_gene"] = gene_wide[gene_cols].idxmax(axis=1)
    gene_wide["best_gene_max_abs_lfc"] = gene_wide[gene_cols].max(axis=1)
    gene_wide["margin_over_runner_up"] = (
        gene_wide[gene_cols].apply(lambda r: r.nlargest(2).iloc[0] - r.nlargest(2).iloc[1], axis=1)
        if len(gene_cols) > 1
        else 0.0
    )
    gene_wide.index.name = "snpId"
    gene_wide = gene_wide.reset_index()
    gene_wide["s8_genes"] = gene_wide["snpId"].map(lambda s: meta_by_snp[s]["s8_genes"])
    gene_wide["pos_grch38"] = gene_wide["snpId"].map(lambda s: meta_by_snp[s]["pos_grch38"])
    gene_wide["ag_best_gene_in_s8"] = [
        g in set(s.split(",")) for g, s in zip(gene_wide["best_gene"], gene_wide["s8_genes"])
    ]
    gene_wide = gene_wide.sort_values("best_gene_max_abs_lfc", ascending=False)
    attrib_path = OUTDIR / f"{args.locus}_gene_attribution.tsv"
    gene_wide.to_csv(attrib_path, sep="\t", index=False)
    print(
        f"wrote {len(gene_wide)} rows -> {attrib_path} "
        f"(best gene counts: {gene_wide['best_gene'].value_counts().to_dict()})"
    )

    for gene, frames in per_gene_all.items():
        gsum = pd.concat(frames, ignore_index=True)
        gsum["s8_genes"] = gsum["snpId"].map(lambda s: meta_by_snp[s]["s8_genes"])
        gsum = gsum.sort_values("max_abs_raw_score", ascending=False)
        gsum.to_csv(
            OUTDIR / f"{args.locus}_{gene.lower()}_gene_summary.tsv", sep="\t", index=False
        )
    print(f"wrote per-gene summaries for {len(per_gene_all)} genes")

    if track_frames:
        for gene, frames in track_frames.items():
            frame = pd.concat(frames, ignore_index=True)
            for label, subset in (
                ("rnaseq", frame[frame["output_type"] == "RNA_SEQ"]),
                ("splice", frame[frame["output_type"].str.startswith("SPLICE")]),
            ):
                if not subset.empty:
                    subset.to_csv(
                        OUTDIR / f"{args.locus}_{gene.lower()}_{label}.tsv", sep="\t", index=False
                    )
        print(f"wrote per-track tables for {len(track_frames)} genes: {sorted(track_frames)}")

    (OUTDIR / f"{args.locus}_run_meta.json").write_text(
        json.dumps(
            {
                "locus": args.locus,
                "interval": {"chromosome": iv.chromosome, "start": iv.start, "end": iv.end},
                "n_unique_snps": len({p["snpId"] for p in provenance}),
                "n_observations": len(provenance),
                "chunk_size": args.chunk_size,
                "n_chunks": len(chunks),
                "s8_genes": sorted({g for p in provenance for g in p["s8_genes"].split(",") if g}),
                "scored_genes": sorted(scored_genes),
                "scorers": [s.base_variant_scorer.name for s in scorers],
                "model_used": str(getattr(model, "model_version", "cached")) if model else "cached",
                "keep_tracks": sorted(keep_tracks) if keep_tracks else None,
            },
            indent=2,
        )
        + "\n"
    )
    print(f"meta -> {OUTDIR / f'{args.locus}_run_meta.json'}")


if __name__ == "__main__":
    main()
