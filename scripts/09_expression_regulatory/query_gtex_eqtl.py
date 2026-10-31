#!/usr/bin/env python3
"""
GTEx v8 cis-eQTL lookups at the chr11 FADS and chr19 APOE loci, across tissues
relevant to lipid metabolism and CAD biology (liver, adipose, artery, adrenal,
whole blood) (Supplementary Table S6).
"""
import os
import requests
import csv
import time

GTEX_API = "https://gtexportal.org/api/v2/association/singleTissueEqtlByLocation"

# GRCh38 liftover of the GRCh37 colocalization lead positions (approximate windows,
# +/- 250kb around the lead SNP, GRCh37->GRCh38 offsets for these loci are ~1kb
# scale on chr11/chr19 so a wide window absorbs the uncertainty safely)
LOCI = {
    "FADS_chr11_omega3": {"chrom": "chr11", "start": 61_550_000, "end": 61_850_000},
    "APOE_chr19_lipids": {"chrom": "chr19", "start": 44_800_000, "end": 45_600_000},
}

TISSUES = [
    "Liver", "Adipose_Subcutaneous", "Adipose_Visceral_Omentum",
    "Artery_Coronary", "Artery_Aorta", "Artery_Tibial",
    "Adrenal_Gland", "Whole_Blood",
]

rows = []
for locus_name, region in LOCI.items():
    for tissue in TISSUES:
        params = {
            "chromosome": region["chrom"],
            "start": region["start"],
            "end": region["end"],
            "tissueSiteDetailId": tissue,
            "datasetId": "gtex_v8",
            "itemsPerPage": 5000,
        }
        try:
            resp = requests.get(GTEX_API, params=params, timeout=60)
            resp.raise_for_status()
            payload = resp.json()
            data = payload.get("singleTissueEqtl", [])
            paging = payload.get("paging_info", {})
        except Exception as e:
            print(f"FAILED {locus_name} / {tissue}: {e}")
            continue
        if paging.get("totalNumberOfItems", 0) > len(data):
            print(f"  NOTE: {locus_name}/{tissue} has {paging.get('totalNumberOfItems')} total items, only {len(data)} fetched (raise itemsPerPage)")

        # GTEx cis-eQTL scans test a much smaller local window than a GWAS, so
        # the conventional GWAS 5e-8 threshold is not the right bar here; use
        # a standard nominal cis-eQTL threshold of 1e-5 instead (note: this is
        # nominal significance, not GTEx's own permutation-based per-gene qval,
        # which this endpoint does not return -- flagged as a limitation).
        sig = [d for d in data if d.get("pValue") is not None and d["pValue"] < 1e-5]
        print(f"{locus_name} / {tissue}: {len(data)} total, {len(sig)} genome-wide-significant eQTLs")
        for d in sig:
            rows.append({
                "locus": locus_name, "tissue": tissue,
                "gencodeId": d.get("gencodeId"), "geneSymbol": d.get("geneSymbol"),
                "snpId": d.get("snpId"), "pos": d.get("pos"),
                "pValue": d.get("pValue"), "nes": d.get("nes"),
            })
        time.sleep(0.3)

out_path = os.path.join(os.environ.get("RAC_ROOT", "."), "results/tissue_enrichment_gtex/gtex_eqtl_hits_fads_apoe.tsv")
os.makedirs(os.path.dirname(out_path), exist_ok=True)
with open(out_path, "w", newline="") as f:
    if rows:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    else:
        f.write("locus\ttissue\tgencodeId\tgeneSymbol\tsnpId\tpos\tpValue\tnes\n")

print(f"\nWrote {len(rows)} significant eQTL rows to {out_path}")

# Summary: which genes recur across tissues at each locus
from collections import Counter
for locus_name in LOCI:
    genes = Counter(r["geneSymbol"] for r in rows if r["locus"] == locus_name)
    print(f"\n{locus_name} -- gene x tissue-count (top 15):")
    for gene, count in genes.most_common(15):
        print(f"  {gene}: {count} tissues")
