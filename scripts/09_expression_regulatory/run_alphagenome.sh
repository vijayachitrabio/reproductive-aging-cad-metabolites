#!/usr/bin/env bash
# AlphaGenome scoring for a full Table S8 locus (default FADS_chr11_omega3).
#
# AlphaGenome is API-only (no public weights) and needs a free non-commercial key.
# Once ALPHA_GENOME_API_KEY is exported, this runs end to end.
#
# Usage: run_alphagenome_locus_run.sh [locus] [chunk_size]   (default FADS_chr11_omega3, 250)
#
# NOTE on chunk_size and the FADS cache: cache files are fingerprinted against the
# exact scoring request (interval, scorers, variant list), so a cache is only reused
# when it matches. The pre-existing FADS cache predates fingerprinting and holds all
# 366 observations in one chunk (written at --chunk-size 1000), so to re-tidy FADS
# offline you must pass 1000:
#   run_alphagenome_locus_run.sh FADS_chr11_omega3 1000
# That cache is honoured with a variant-count check and a warning that its interval
# and scorers are unverified; any other chunk size is a clean cache miss.
set -euo pipefail

cd "$(dirname "$0")/../.."

LOCUS="${1:-FADS_chr11_omega3}"
ENV_PY="${ALPHAGENOME_PYTHON:-$HOME/miniconda3/envs/alphagenome/bin/python}"

if [ ! -x "$ENV_PY" ]; then
  echo "error: AlphaGenome interpreter not found at $ENV_PY" >&2
  echo "create it with:" >&2
  echo "  conda create -n alphagenome python=3.11 -y" >&2
  echo "  $ENV_PY -m pip install -U alphagenome requests" >&2
  exit 1
fi

echo "step 1/2: resolving GRCh38 REF/ALT for $LOCUS (Ensembl REST)"
"$ENV_PY" scripts/09_expression_regulatory/alphagenome_fetch_alleles.py --locus "$LOCUS"

echo "step 2/2: scoring $LOCUS with AlphaGenome"
# Only reuse the cache when it already exists, so a first run still queries the API.
# The scorer caches per chunk under cache/, so an interrupted run resumes rather than restarting.
# Chunk size matters: a single 1 Mb call over 1,794 APOE observations would be ~40M tidy rows.
CHUNK="${2:-250}"
if compgen -G "results/alphagenome_locus/cache/${LOCUS}_chunk*.pkl" > /dev/null; then
  "$ENV_PY" scripts/09_expression_regulatory/alphagenome_score_locus.py --locus "$LOCUS" --chunk-size "$CHUNK" --reuse-cache
else
  "$ENV_PY" scripts/09_expression_regulatory/alphagenome_score_locus.py --locus "$LOCUS" --chunk-size "$CHUNK"
fi

echo
echo "outputs in results/alphagenome_locus/:"
echo "  ${LOCUS}_variants_grch38.tsv        variants with GRCh38 REF/ALT"
echo "  ${LOCUS}_variant_summary.tsv        gene-agnostic variant x output_type"
echo "  ${LOCUS}_gene_attribution.tsv       per-variant, which gene carries the effect"
echo "  ${LOCUS}_<gene>_gene_summary.tsv    per-gene variant x output_type"
echo "  ${LOCUS}_genes_in_window.json       Ensembl gene set used for attribution"
echo "  ${LOCUS}_run_meta.json              interval, scorers, genes"
echo "  cache/${LOCUS}_chunk*.pkl           per-chunk AnnData cache (--reuse-cache re-tidies free)"
echo
echo "per-gene per-track tables are opt-in: add --keep-tracks FADS1,FADS2 (large)"
echo "figure: python3 alphagenome/plot_alphagenome_fads_locus.py   (FADS locus only)"
