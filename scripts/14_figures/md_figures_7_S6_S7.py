# Script by Naisarg Patel (contribution to this study).
# Plots the molecular-dynamics figures: Figure 7 (output Fig1_RMSD_RMSF_Rg_SASA.png),
# Supplementary Figure S6 (Fig2_LigandRMSD_HBonds.png) and Supplementary Figure S7 (Fig3_PCA_FEL.png).
# Inputs are the GROMACS .xvg files written by 13_molecular_dynamics/analysis: put the apo, c1 and c2
# files in xvg_files/ and the t (D5D-IN-326) files in T/, next to this script. The last block
# (Fig4_Summary_Tables.png) also needs fig_v2/Average_Metrics_v2.csv and MMGBSA_Results_v2.csv.
"""
plot_mds_figures_v2.py  –  FADS1 MD figures, v2: adds D5D-IN-326 (T/)

Writes to fig_v2/. Same plotting as plot_mds_figures.py; NU7441 moves to
orange so D5D-IN-326 can take red.

Figure 1 : 2x2  – Protein RMSD, RMSF, Rg, SASA      (single shared legend)
Figure 2 : 1x2  – Ligand RMSD, Protein-ligand H-bonds (single shared legend)
Figure 3 : 2x4  – PCA projections (top row) + Free Energy Landscapes (bottom)
Figure 4 : 2x1  – fig_v2/Average_Metrics_v2.csv and MMGBSA_Results_v2.csv rendered as tables
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.lines import Line2D
from pathlib import Path

# ── Configuration ─────────────────────────────────────────────────────────────
ROOT     = Path(__file__).parent
OUT_DIR  = ROOT / "fig_v2"
OUT_DIR.mkdir(exist_ok=True)

SYSTEMS = [
    ("apo", "Apo",         "black"),
    ("c1",  "Luminespib",  "#638ccc"),
    ("c2",  "NU7441",      "#e69f00"),
    ("t",   "D5D-IN-326",  "#c9556f"),
]
DATA_DIR = {"apo": ROOT / "xvg_files", "c1": ROOT / "xvg_files",
            "c2": ROOT / "xvg_files", "t": ROOT / "T"}
LIG_SYSTEMS = [s for s in SYSTEMS if s[0] != "apo"]   # apo has no ligand data

LETTERS  = "ABCDEFGH"
kT       = 2.479    # kJ/mol at ~298 K
FEL_VMAX = 15       # kJ/mol colour ceiling (shared across panels)
FEL_BINS = 60       # matches BINS in PCA.R. PC projections span ~30-35 nm here,
                    # so finer binning (80) leaves most bins with 0-1 frames and
                    # the landscape breaks up into single-frame specks.
RMSF_YMAX = 1.0     # nm – flexible N-terminus runs off-scale (see annotation)
LIG_SPIKE = 2.0     # nm – ligand-RMSD frames above this are the 3-frame PBC
                    # reimaging glitch at ~51.6 ns (C2 restart seam); dropped.
                    # D5D-IN-326 never exceeds 1.1 nm, so nothing dropped there
LIG_YMAX  = 3.0     # nm

# ── Helpers ───────────────────────────────────────────────────────────────────
def read_xvg(filepath):
    """Parse a GROMACS .xvg file, skipping comment/header lines."""
    rows = []
    with open(filepath) as fh:
        for line in fh:
            s = line.strip()
            if not s or s[0] in ("#", "@"):
                continue
            rows.append([float(v) for v in s.split()])
    return np.array(rows)


def smooth(x, y, frac=0.05, minw=50):
    """Running average via convolution; returns trimmed (x, y_smooth)."""
    w = max(minw, int(len(y) * frac))
    w = min(w, len(y) // 2)
    kernel = np.ones(w) / w
    ys = np.convolve(y, kernel, mode="valid")
    half = (w - 1) // 2
    xs = x[half: half + len(ys)]
    return xs, ys


def compute_fel(x, y, bins=FEL_BINS):
    """2-D free energy landscape from PC1/PC2 projections."""
    hist, xe, ye = np.histogram2d(x, y, bins=bins)
    hist = hist.astype(float)
    hist[hist == 0] = np.nan
    prob = hist / np.nansum(hist)
    fel  = -kT * np.log(prob)
    fel -= np.nanmin(fel)
    xc = 0.5 * (xe[:-1] + xe[1:])
    yc = 0.5 * (ye[:-1] + ye[1:])
    return xc, yc, fel.T          # transpose: rows=y, cols=x for contourf


def trace(ax, t, y, color, label, raw_alpha=0.15):
    """Raw trace (faint) + running average (bold)."""
    ax.plot(t, y, color=color, alpha=raw_alpha, linewidth=0.6)
    xs, ys = smooth(t, y)
    ax.plot(xs, ys, color=color, linewidth=1.8, label=label)


def panel_title(ax, letter, text, **kw):
    ax.set_title(f"{letter}. {text}", fontweight="bold", loc="left", **kw)


def legend_handles(systems):
    return [Line2D([0], [0], color=c, linewidth=2.2, label=l)
            for _, l, c in systems]


# ── Style ─────────────────────────────────────────────────────────────────────
plt.rcParams.update({
    "font.family":       "sans-serif",
    "font.size":         10,
    "axes.linewidth":    0.8,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "legend.frameon":    False,
})

# ═════════════════════════════════════════════════════════════════════════════
# Figure 1 – RMSD / RMSF / Rg / SASA
# ═════════════════════════════════════════════════════════════════════════════
fig1, axes1 = plt.subplots(2, 2, figsize=(11, 8))
ax_rmsd, ax_rmsf, ax_rg, ax_sasa = axes1.flat

for key, label, color in SYSTEMS:
    # A – Protein backbone RMSD (x already in ns)
    d = read_xvg(DATA_DIR[key] / f"{key}_protein_rmsd.xvg")
    trace(ax_rmsd, d[:, 0], d[:, 1], color, label)

    # B – per-residue RMSF (no smoothing)
    d = read_xvg(DATA_DIR[key] / f"{key}_rmsf.xvg")
    ax_rmsf.plot(d[:, 0], d[:, 1], color=color, linewidth=1.3, label=label)

    # C – radius of gyration (x in ps -> ns; col 1 = total Rg)
    d = read_xvg(DATA_DIR[key] / f"{key}_gyrate.xvg")
    trace(ax_rg, d[:, 0] / 1000.0, d[:, 1], color, label)

    # D – solvent accessible surface area (x already in ns)
    d = read_xvg(DATA_DIR[key] / f"{key}_area.xvg")
    trace(ax_sasa, d[:, 0], d[:, 1], color, label)

ax_rmsd.set_xlabel("Time (ns)")
ax_rmsd.set_ylabel("RMSD (nm)")
panel_title(ax_rmsd, "A", "Protein backbone RMSD")

ax_rmsf.set_xlabel("Residue")
ax_rmsf.set_ylabel("RMSF (nm)")
ax_rmsf.set_ylim(0, RMSF_YMAX)
ax_rmsf.annotate("N-terminal residues 1–10\noff-scale (up to 2.8 nm)",
                 xy=(0.30, 0.97), xycoords="axes fraction",
                 fontsize=7.5, color="0.35", va="top")
panel_title(ax_rmsf, "B", "Per-residue RMSF")

ax_rg.set_xlabel("Time (ns)")
ax_rg.set_ylabel("R$_g$ (nm)")
panel_title(ax_rg, "C", "Radius of gyration")

ax_sasa.set_xlabel("Time (ns)")
ax_sasa.set_ylabel("SASA (nm$^2$)")
panel_title(ax_sasa, "D", "Solvent accessible surface area")

fig1.tight_layout(pad=2.0)
fig1.subplots_adjust(bottom=0.12)
fig1.legend(handles=legend_handles(SYSTEMS), loc="lower center",
            ncol=len(SYSTEMS), fontsize=10, frameon=False,
            bbox_to_anchor=(0.5, 0.015))

out1 = OUT_DIR / "Fig1_RMSD_RMSF_Rg_SASA.png"
fig1.savefig(out1, dpi=600, bbox_inches="tight")
print(f"Saved: {out1}")

# ═════════════════════════════════════════════════════════════════════════════
# Figure 2 – Ligand RMSD + protein–ligand hydrogen bonds
# ═════════════════════════════════════════════════════════════════════════════
fig2, (ax_lig, ax_hb) = plt.subplots(1, 2, figsize=(11, 4.2))

for key, label, color in LIG_SYSTEMS:
    # A – ligand heavy-atom RMSD after fit on protein backbone
    d = read_xvg(DATA_DIR[key] / f"{key}_lig_rmsd_vs_start.xvg")
    spike = d[:, 1] > LIG_SPIKE
    if spike.any():
        print(f"  {key}: dropped {spike.sum()} ligand-RMSD frame(s) > "
              f"{LIG_SPIKE} nm at t = {np.round(d[spike, 0], 2).tolist()} ns")
    d = d[~spike]
    trace(ax_lig, d[:, 0], d[:, 1], color, label)

    # B – H-bond count (col 1; col 2 is pairs within 0.35 nm)
    d = read_xvg(DATA_DIR[key] / f"{key}_hbnum.xvg")
    trace(ax_hb, d[:, 0], d[:, 1], color, label, raw_alpha=0.07)

ax_lig.set_xlabel("Time (ns)")
ax_lig.set_ylabel("RMSD (nm)")
ax_lig.set_ylim(0, LIG_YMAX)
panel_title(ax_lig, "A", "Ligand RMSD")

ax_hb.set_xlabel("Time (ns)")
ax_hb.set_ylabel("Number of H-bonds")
ax_hb.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
panel_title(ax_hb, "B", "Protein–ligand hydrogen bonds")

fig2.tight_layout(pad=2.0)
fig2.subplots_adjust(bottom=0.22)
fig2.legend(handles=legend_handles(LIG_SYSTEMS), loc="lower center",
            ncol=len(LIG_SYSTEMS), fontsize=10, frameon=False,
            bbox_to_anchor=(0.5, 0.01))

out2 = OUT_DIR / "Fig2_LigandRMSD_HBonds.png"
fig2.savefig(out2, dpi=600, bbox_inches="tight")
print(f"Saved: {out2}")

# ═════════════════════════════════════════════════════════════════════════════
# Figure 3 – PCA projections (top) and free energy landscapes (bottom)
# ═════════════════════════════════════════════════════════════════════════════
fig3, axes3 = plt.subplots(2, 4, figsize=(20, 8.6))

proj   = {k: read_xvg(DATA_DIR[k] / f"{k}_2dproj.xvg") for k, _, _ in SYSTEMS}
tmax   = 100.0                                   # ns (10001 frames, 0–100 ns)
levels = np.linspace(0, FEL_VMAX, 21)

for i, (key, label, color) in enumerate(SYSTEMS):
    x, y = proj[key][:, 0], proj[key][:, 1]
    t    = np.linspace(0, tmax, len(x))

    # ── Top row: PC1 vs PC2 scatter, coloured by simulation time ─────────────
    ax = axes3[0, i]
    sc = ax.scatter(x, y, c=t, cmap="viridis", s=2, alpha=0.6,
                    vmin=0, vmax=tmax, rasterized=True)
    ax.set_xlabel("PC1 (nm)")
    ax.set_ylabel("PC2 (nm)")
    panel_title(ax, LETTERS[i], f"{label} – PCA", color=color)

    # ── Bottom row: Gibbs free energy landscape ──────────────────────────────
    ax = axes3[1, i]
    xc, yc, fel = compute_fel(x, y)
    cf = ax.contourf(xc, yc, fel, levels=levels, cmap="jet", extend="max")
    ax.contour(xc, yc, fel, levels=np.linspace(0, FEL_VMAX, 11),
               colors="white", linewidths=0.3, alpha=0.35)
    ax.set_xlabel("PC1 (nm)")
    ax.set_ylabel("PC2 (nm)")
    panel_title(ax, LETTERS[i + len(SYSTEMS)], f"{label} – FEL", color=color)

fig3.tight_layout(pad=2.0)
fig3.subplots_adjust(right=0.92)

# colour bars aligned to the right-hand column of each row
def cbar_axes(ax):
    p = ax.get_position()
    return fig3.add_axes([0.935, p.y0, 0.010, p.height])

cax_top = cbar_axes(axes3[0, -1])
cb_top  = fig3.colorbar(sc, cax=cax_top)
cb_top.set_label("Time (ns)", fontsize=9)
cb_top.ax.tick_params(labelsize=8)

cax_bot = cbar_axes(axes3[1, -1])
cb_bot  = fig3.colorbar(cf, cax=cax_bot,
                        ticks=np.arange(0, FEL_VMAX + 1, 3))
cb_bot.set_label("ΔG (kJ/mol)", fontsize=9)
cb_bot.ax.tick_params(labelsize=8)

out3 = OUT_DIR / "Fig3_PCA_FEL.png"
fig3.savefig(out3, dpi=600, bbox_inches="tight")
print(f"Saved: {out3}")


# ═════════════════════════════════════════════════════════════════════════════
# Figure 4 – the two result CSVs rendered as tables
# ═════════════════════════════════════════════════════════════════════════════
import csv

HEADER_BG = "#3c4657"
STRIPE_BG = "#f2f4f7"
TOTAL_BG  = "#dde5f0"
NA_DASH   = "\u2014"          # em dash for missing values
SYS_COLOR = {l: c for _, l, c in SYSTEMS}


def read_csv(path):
    with open(path, encoding="utf-8") as fh:
        rows = [r for r in csv.reader(fh) if r and any(f.strip() for f in r)]
    return rows[0], rows[1:]


def fmt_num(v, nd=3):
    """Trim the long floats in Average_Metrics_v2.csv to a readable precision."""
    v = v.strip()
    if v in ("", "N/A", "NA", "nan"):
        return NA_DASH
    try:
        return f"{float(v):.{nd}f}"
    except ValueError:
        return v


def col_widths(header, body, first_col_pad=3):
    """Column widths proportional to the widest string in each column."""
    w = []
    for c, h in enumerate(header):
        longest = max([len(h)] + [len(r[c]) for r in body])
        w.append(longest + (first_col_pad if c == 0 else 2))
    total = sum(w)
    return [x / total for x in w]


def draw_table(ax, header, body, title, letter, bold_rows=(),
               color_first_col=False):
    ax.axis("off")
    panel_title(ax, letter, title)

    tbl = ax.table(cellText=body, colLabels=header, cellLoc="center",
                   colWidths=col_widths(header, body),
                   bbox=[0, 0, 1, 1])          # fill the axes exactly
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)

    for (r, c), cell in tbl.get_celld().items():
        cell.set_edgecolor("white")
        cell.set_linewidth(1.2)
        if r == 0:                                   # header row
            cell.set_facecolor(HEADER_BG)
            cell.set_text_props(color="white", fontweight="bold")
            continue
        cell.set_facecolor(STRIPE_BG if r % 2 else "white")
        if c == 0:                                   # row-label column
            cell.set_text_props(ha="left", fontweight="bold")
            cell.PAD = 0.03
            if color_first_col:
                cell.set_text_props(color=SYS_COLOR.get(body[r - 1][0], "black"),
                                    ha="left", fontweight="bold")
        if r - 1 in bold_rows:                       # highlighted row (TOTAL)
            cell.set_facecolor(TOTAL_BG)
            cell.set_text_props(fontweight="bold")
    return tbl


# ── Panel A – average structural metrics ─────────────────────────────────────
hdr_a, rows_a = read_csv(OUT_DIR / "Average_Metrics_v2.csv")
body_a = [[r[0]] + [fmt_num(v, 2 if "SASA" in hdr_a[i] else 3)
                    for i, v in enumerate(r[1:], start=1)] for r in rows_a]

# ── Panel B – MM-GBSA binding free energy ────────────────────────────────────
hdr_b, rows_b = read_csv(OUT_DIR / "MMGBSA_Results_v2.csv")
body_b = [[r[0]] + [v.strip().replace("+/-", "\u00b1") for v in r[1:]]
          for r in rows_b]
total_idx = tuple(i for i, r in enumerate(rows_b) if r[0].strip().upper() == "TOTAL")

# height of each panel tracks its row count so the cells stay evenly sized
rows_A, rows_B = len(body_a) + 1, len(body_b) + 1
ROW_IN = 0.42
fig4 = plt.figure(figsize=(12, (rows_A + rows_B) * ROW_IN + 1.3))
gs = fig4.add_gridspec(2, 1, height_ratios=[rows_A, rows_B], hspace=0.30)

draw_table(fig4.add_subplot(gs[0]), hdr_a, body_a,
           "Average structural metrics (100 ns)", "A", color_first_col=True)
draw_table(fig4.add_subplot(gs[1]), hdr_b, body_b,
           "MM-GBSA binding free energy", "B", bold_rows=total_idx)

out4 = OUT_DIR / "Fig4_Summary_Tables.png"
fig4.savefig(out4, dpi=600, bbox_inches="tight", facecolor="white")
print(f"Saved: {out4}")
