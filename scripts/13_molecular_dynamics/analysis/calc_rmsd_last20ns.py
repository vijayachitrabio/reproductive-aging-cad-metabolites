"""Mean and SD of protein-backbone and ligand RMSD over the last 20 ns (80-100 ns) of the FADS1 MD runs.

The 80-100 ns window is the MM-GBSA window. The two 10-ns halves (80-90, 90-100 ns) are given as a
simple drift check. SD is the sample SD (ddof=1) of the frames inside the window, i.e. fluctuation
within one trajectory; neighbouring frames are correlated, so it is not an uncertainty across runs.
Ligand RMSD frames above 2.0 nm are dropped, as in the figure script (md_figures_7_S6_S7.py).

Usage:
  python calc_rmsd_last20ns.py --xvg-dir PATH_WITH_apo_c1_c2_XVG --t-dir PATH_WITH_T_XVG --out rmsd_last20ns_80-100ns.csv
Reads {key}_protein_rmsd.xvg and {key}_lig_rmsd_vs_start.xvg (time in ns, RMSD in nm).
Keys: apo, c1 (Luminespib), c2 (NU7441), t (D5D-IN-326).
"""
import argparse
import csv
from pathlib import Path

import numpy as np

SYSTEMS = [("apo", "Apo"), ("c1", "Luminespib"), ("c2", "NU7441"), ("t", "D5D-IN-326")]
LIG_SPIKE = 2.0      # nm, as in the figure script
WINDOW = (80.0, 100.0)
HALVES = ((80.0, 90.0), (90.0, 100.0))


def read_xvg(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            s = line.strip()
            if s and s[0] not in "#@":
                rows.append([float(v) for v in s.split()])
    return np.array(rows)


def window(t, y, lo, hi):
    return y[(t >= lo) & (t <= hi)]


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--xvg-dir", required=True, type=Path, help="folder with apo_, c1_, c2_ xvg files")
    ap.add_argument("--t-dir", required=True, type=Path, help="folder with t_ (D5D-IN-326) xvg files")
    ap.add_argument("--out", default="rmsd_last20ns_80-100ns.csv", type=Path)
    a = ap.parse_args()

    header = ["metric", "system", "full_run_mean_nm", "mean_80_100ns_nm", "sd_80_100ns_nm",
              "n_frames_80_100ns", "mean_80_90ns_nm", "mean_90_100ns_nm"]
    out = []
    for metric, suffix in (("Protein backbone RMSD", "protein_rmsd"), ("Ligand RMSD vs start", "lig_rmsd_vs_start")):
        for key, label in SYSTEMS:
            if metric.startswith("Ligand") and key == "apo":
                continue
            d = read_xvg((a.t_dir if key == "t" else a.xvg_dir) / f"{key}_{suffix}.xvg")
            if metric.startswith("Ligand"):
                d = d[d[:, 1] <= LIG_SPIKE]
            t, y = d[:, 0], d[:, 1]
            w = window(t, y, *WINDOW)
            out.append([metric, label, round(y.mean(), 3), round(w.mean(), 3), round(w.std(ddof=1), 3), len(w),
                        round(window(t, y, *HALVES[0]).mean(), 3), round(window(t, y, *HALVES[1]).mean(), 3)])
    with open(a.out, "w", newline="") as fh:
        csv.writer(fh).writerows([header] + out)
    for r in out:
        print(*r, sep="\t")


if __name__ == "__main__":
    main()
