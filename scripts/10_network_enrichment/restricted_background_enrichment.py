#!/usr/bin/env python3
"""Run g:Profiler enrichment with the reconstructed network as custom background.

The query is each Louvain module; the background is all genes represented by that
same network setting. g:Profiler's FDR option is used, and only one term per
redundancy group and source is retained for the compact report.
"""
import os
import csv
import json
import ssl
import urllib.request
from pathlib import Path

ROOT = Path(os.environ.get("RAC_ROOT", "."))
BASE = ROOT / "results" / "ppi_network"
MEMBERSHIP = BASE / "robustness" / "robustness_module_membership.tsv"
OUT = BASE / "enrichment_restricted_background"
OUT.mkdir(parents=True, exist_ok=True)

rows = list(csv.DictReader(MEMBERSHIP.open(), delimiter="\t"))
settings = {}
for row in rows:
    settings.setdefault(row["setting"], []).append(row)

ctx = ssl._create_unverified_context()
url = "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
all_results = []

def module_label(genes):
    g = set(genes)
    if "FADS1" in g:
        return "FADS1_lipid_module"
    if g & {"AMH", "AMHR2"}:
        return "AMH_AMHR2_module"
    if g & {"FSHR", "LHCGR", "FSHB", "CGA"}:
        return "follicular_development_module"
    if g & {"PTGS1", "PTGS2"}:
        return "cyclooxygenase_module"
    if g & {"ALOX5", "ALOX5AP", "ALOX15", "ALOX15B"}:
        return "lipoxygenase_module"
    return "other_module"

for setting, setting_rows in settings.items():
    background = sorted({r["gene"] for r in setting_rows})
    modules = {}
    for r in setting_rows:
        modules.setdefault(r["module"], []).append(r["gene"])
    for module, genes in modules.items():
        payload = {
            "organism": "hsapiens",
            "query": sorted(set(genes)),
            "background": background,
            "sources": ["GO:BP", "KEGG", "REAC"],
            "user_threshold": 0.05,
            "significance_threshold_method": "fdr",
            "domain_scope": "custom_annotated",
            "no_evidences": False,
        }
        req = urllib.request.Request(url, data=json.dumps(payload).encode(), headers={"Content-Type": "application/json"})
        try:
            with urllib.request.urlopen(req, timeout=45, context=ctx) as response:
                data = json.loads(response.read().decode())
        except Exception as exc:
            (OUT / f"{setting}_module_{module}.error.txt").write_text(str(exc))
            continue
        (OUT / f"{setting}_module_{module}.json").write_text(json.dumps(data, indent=2))
        for item in data.get("result", []):
            if item.get("significant") is not True:
                continue
            all_results.append({
                "setting": setting,
                "module": module,
                "module_label": module_label(genes),
                "query_size": len(set(genes)),
                "background_size": len(background),
                "source": item.get("source", ""),
                "term_id": item.get("native", ""),
                "term": item.get("name", ""),
                "description": item.get("description", ""),
                "adjusted_p_value": item.get("p_value", ""),
                "intersection_size": item.get("intersection_size", ""),
                "term_size_in_background": item.get("term_size", ""),
                "precision": item.get("precision", ""),
                "recall": item.get("recall", ""),
                "redundancy_group": item.get("group_id", ""),
                "genes": ",".join(sorted(set(genes))),
            })

fields = ["setting", "module", "module_label", "query_size", "background_size", "source", "term_id", "term", "description", "adjusted_p_value", "intersection_size", "term_size_in_background", "precision", "recall", "redundancy_group", "genes"]
with (OUT / "all_significant_terms.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(all_results)

# Compact nonredundant view: one best term per source/redundancy group.
best = {}
for row in all_results:
    try:
        p = float(row["adjusted_p_value"])
    except (TypeError, ValueError):
        continue
    key = (row["setting"], row["module"], row["source"], row["redundancy_group"])
    if key not in best or p < float(best[key]["adjusted_p_value"]):
        best[key] = row
with (OUT / "nonredundant_terms.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(sorted(best.values(), key=lambda r: (r["setting"], int(r["module"]), r["source"], float(r["adjusted_p_value"]))))

# Focused stability table for PUFA/fatty-acid biology in the FADS1 module.
needles = ("fatty acid", "unsaturated", "polyunsaturated", "PUFA", "elongation", "linoleic", "alpha-linolenic", "lipid metabolism")
focused = [r for r in best.values() if r["module_label"] == "FADS1_lipid_module" and any(n.lower() in (r["term"] + " " + r["description"]).lower() for n in needles)]
with (OUT / "fads1_module_fatty_acid_terms.tsv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
    writer.writeheader()
    writer.writerows(sorted(focused, key=lambda r: (r["setting"], r["source"], float(r["adjusted_p_value"]))))

print(f"Wrote {len(all_results)} significant terms and {len(best)} nonredundant terms to {OUT}")
