"""Graphical Abstract for the IJC submission (author decision 2026-09-12, checklist item A10).

IJC 3.5 asks for a roughly square, colour image that summarises the key findings, carries only
labels and short phrases, and contains neither the manuscript title nor author names; it is read
alongside the "What's New" text.

The abstract is built from the derived outputs, so its numbers cannot drift away from the
manuscript, and it shows the data rather than describing it: a forest plot of the external
replication, the within-donor paired comparison, and the two same-sample scatters that carry the
central contrast. Rank axes are used wherever a Spearman correlation is quoted, as in Figures 2
and 3. Every connector begins and ends on a box border, so nothing overlaps.

Output: figures/Graphical_Abstract_v3.{pdf,png} — 180 x 180 mm, 600 dpi.
Run with the project interpreter (~/miniforge3/bin/python3; Python 3.9.2, SciPy 1.13.1).
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from scipy import stats

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
SRC = os.path.join(BASE, "14_Revision_v3_20260806/outputs")
FIG = os.path.join(BASE, "14_Revision_v3_20260806/figures")

plt.rcParams.update({
    "font.size": 8, "font.family": ["Arial", "DejaVu Sans"],  # IJC: Helvetica/Arial in figures
    "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold", "savefig.bbox": "tight", "pdf.fonttype": 42, "ps.fonttype": 42,
    "axes.linewidth": 0.7, "xtick.major.width": 0.7, "ytick.major.width": 0.7,
    "xtick.major.size": 2.4, "ytick.major.size": 2.4,
})
BLUE, RED, GREY, GREEN, ORANGE = "#2C5F8A", "#A81E22", "#666666", "#2E7D53", "#C4762A"
INK, PALE_BLUE, PALE_GREEN, PALE_GREY = "#15304A", "#EDF3F8", "#EDF5F0", "#F4F4F4"


def minus(s):
    return str(s).replace("-", "−")


# ---- numbers, read from the derived outputs (never transcribed) ----------------------------
meta = pd.read_csv(os.path.join(SRC, "Meta_TwoExternalCohorts_v3.csv")).iloc[0]
s9 = pd.read_csv(os.path.join(SRC, "SuppTable_External_Replication_and_Leukocyte_Sensitivity_v3.csv"))
spear = s9[s9.model == "Spearman DLL3~APM"].drop_duplicates("cohort").set_index("cohort")
ptprc = s9[s9.model == "std OLS APM~DLL3 + PTPRC"].set_index("cohort")
ccle = pd.read_csv(os.path.join(SRC, "SuppTable_CCLE_LeukocyteFree_Falsification_v3.csv"))
ccle_apm = ccle[ccle.comparison.str.startswith("DLL3 ~ APM (5-gene")].iloc[0]
depth = pd.read_csv(os.path.join(SRC, "SuppTable_DetectionDepth_Sensitivity_v3.csv"))
adj = depth[(depth.analysis == "cell_level_logistic")
            & (depth.model == "+ depth + subtype + donor")].set_index("quantity")
paired = pd.read_csv(os.path.join(SRC, "Plotdata_Chan_Paired_Observed_v3.csv"))
per = pd.read_csv(os.path.join(SRC, "Chan_PerDonor_with_epithelialAPM_v3.csv"))
ax12 = pd.read_csv(os.path.join(SRC, "SuppTable_Axis1_vs_Axis2_and_Equivalence_v3.csv"))
prin = ax12[(ax12.analysis == "Axis1_vs_Axis2_direct")
            & (ax12.y == "TNFRSF9pos_pct_of_CD8")].iloc[0]
null = ax12[ax12.analysis == "DLL3_to_Axis2_equivalence"]
null_prin = null[null.y == "TNFRSF9pos_pct_of_CD8"].iloc[0]

fig = plt.figure(figsize=(7.087, 7.087))           # 180 x 180 mm
bg = fig.add_axes([0, 0, 1, 1])
bg.set_xlim(0, 1)
bg.set_ylim(0, 1)
bg.axis("off")


def box(x0, y0, x1, y1, face, edge, lw=1.2):
    bg.add_patch(FancyBboxPatch((x0, y0), x1 - x0, y1 - y0,
                                boxstyle="round,pad=0,rounding_size=0.018",
                                facecolor=face, edgecolor=edge, linewidth=lw, zorder=1))


def say(x, y, s, size=7.6, weight="normal", color="#111111", ha="center", va="center"):
    bg.text(x, y, s, fontsize=size, fontweight=weight, color=color, ha=ha, va=va,
            zorder=4, linespacing=1.4)


def tidy(ax, size=6.2):
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.tick_params(labelsize=size, pad=1.6)
    ax.set_facecolor("none")


# ============================================================ exposure node
box(0.300, 0.900, 0.700, 0.990, PALE_GREY, "#3A3A3A", lw=1.3)
say(0.500, 0.966, "DLL3", size=15, weight="bold")
say(0.500, 0.9455, "baseline tumour expression, modelled continuously", size=7.0, color="#333333")
cb = fig.add_axes([0.385, 0.9125, 0.230, 0.015])
cb.imshow(np.linspace(0, 1, 256).reshape(1, -1), aspect="auto",
          cmap=matplotlib.colors.LinearSegmentedColormap.from_list("d", ["#DCE6EF", INK]))
cb.set_xticks([])
cb.set_yticks([])
for sp in cb.spines.values():
    sp.set_linewidth(0.6)
    sp.set_color("#3A3A3A")
say(0.378, 0.920, "low", size=6.4, color="#444444", ha="right")
say(0.622, 0.920, "high", size=6.4, color="#444444", ha="left")

# ============================================================ arrows out of the node
# each starts on the node's bottom border and ends on the panel's top border
bg.add_patch(FancyArrowPatch((0.420, 0.900), (0.262, 0.847), arrowstyle="-|>", mutation_scale=13,
                             linewidth=1.7, color=BLUE, shrinkA=0, shrinkB=0, zorder=3))
say(0.105, 0.884, "replicated in two\nexternal cohorts", size=7.2, weight="bold", color=BLUE)

bg.add_patch(FancyArrowPatch((0.580, 0.900), (0.738, 0.847), arrowstyle="-|>", mutation_scale=13,
                             linewidth=1.7, color=GREY, linestyle=(0, (3.5, 2.5)), shrinkA=0,
                             shrinkB=0, zorder=3))
say(0.897, 0.884, "association\nnot demonstrated", size=7.2, weight="bold", color=RED)

# ============================================================ feature 1
box(0.030, 0.560, 0.480, 0.845, PALE_BLUE, BLUE)
say(0.255, 0.821, "Antigen-presentation signal", size=9.4, weight="bold", color=INK)
say(0.255, 0.795, "lower in DLL3-high tumours", size=7.8)

fa = fig.add_axes([0.132, 0.658, 0.308, 0.112])
rows = [("Jiang\n$n$ = 79", spear.loc["Jiang (GSE60052)"]),
        ("Cologne\n$n$ = 81", spear.loc["George/Cologne (cBioPortal)"])]
for i, (lab, r) in enumerate(rows):
    y = 2 - i
    fa.plot([r.ci_low_fisher, r.ci_high_fisher], [y, y], color=BLUE, lw=1.3, zorder=2)
    fa.plot([r.estimate], [y], "o", color=BLUE, ms=5, mec="#0E2537", mew=0.5, zorder=3)
fa.plot([meta.ci_low, meta.ci_high], [0, 0], color=RED, lw=1.5, zorder=2)
fa.plot([meta.pooled_r], [0], "D", color=RED, ms=6, mec="#5A1214", mew=0.5, zorder=3)
fa.axvline(0, color="#888888", ls=(0, (2, 2)), lw=0.8, zorder=1)
fa.set_yticks([2, 1, 0])
fa.set_yticklabels([rows[0][0], rows[1][0], "Pooled\n$N$ = 160"], fontsize=6.0)
fa.set_ylim(-0.75, 2.75)
fa.set_xlim(-0.72, 0.14)
fa.set_xticks([-0.6, -0.4, -0.2, 0])
fa.set_xticklabels([minus("-0.6"), minus("-0.4"), minus("-0.2"), "0"])
fa.set_xlabel("Spearman $r$ (DLL3 vs antigen presentation)", fontsize=6.4, labelpad=1.5)
tidy(fa)
fa.text(0.02, 0.05, f"pooled $r$ = {minus(f'{meta.pooled_r:+.3f}')}", transform=fa.transAxes,
        fontsize=6.6, color=RED, weight="bold", ha="left")

say(0.255, 0.589, f"but $PTPRC$-adjusted $P$ = {ptprc.loc['Jiang (GSE60052)', 'p']:.2f} and "
                  f"{ptprc.loc['George/Cologne (cBioPortal)', 'p']:.2f}, and absent in 59\n"
                  f"leukocyte-free cell lines ($r$ = {minus(f'{ccle_apm.spearman_r:+.3f}')}): "
                  f"a compositional signal", size=6.9, color=GREY)

# ============================================================ feature 2
box(0.520, 0.560, 0.970, 0.845, PALE_BLUE, BLUE)
say(0.745, 0.821, "$TNFRSF9^+$ CD8 T cells", size=9.4, weight="bold", color=INK)
say(0.745, 0.795, "PD-1/TIM-3-marker-high phenotype", size=7.8)

pa = fig.add_axes([0.622, 0.658, 0.308, 0.112])
for mk, colour, x0 in [("PD-1", BLUE, 0.0), ("TIM-3", ORANGE, 1.6)]:
    lo, hi = paired[f"{mk}_neg"].values, paired[f"{mk}_pos"].values
    for a, b in zip(lo, hi):
        pa.plot([x0, x0 + 0.75], [a, b], color=colour, lw=0.7, alpha=0.5, zorder=2)
    pa.plot([x0] * len(lo), lo, "o", color=colour, ms=2.8, alpha=0.8, zorder=3)
    pa.plot([x0 + 0.75] * len(hi), hi, "o", color=colour, ms=2.8, alpha=0.8, zorder=3)
    pa.text(x0 + 0.375, 108, mk, ha="center", fontsize=6.8, color=colour, weight="bold")
pa.set_xticks([0, 0.75, 1.6, 2.35])
pa.set_xticklabels(["$TNFRSF9^-$", "$TNFRSF9^+$"] * 2, fontsize=5.5)
pa.set_xlim(-0.45, 2.8)
pa.set_ylim(-5, 120)
pa.set_yticks([0, 50, 100])
pa.set_ylabel("marker$^+$ (% of CD8)", fontsize=6.4, labelpad=1.5)
tidy(pa)

say(0.745, 0.589, f"higher in 14/15 (PD-1) and 15/15 (TIM-3) donors;\n"
                  f"depth-, subtype- and donor-adjusted odds ratio "
                  f"{adj.loc['PD-1', 'odds_ratio']:.2f} and {adj.loc['TIM-3', 'odds_ratio']:.2f}",
    size=6.9, color=GREY)

# ============================================================ the central contrast
box(0.030, 0.180, 0.970, 0.520, "#FFFFFF", "#3A3A3A")
say(0.500, 0.494, "In the same 19 donors, what tracks the abundance of these cells?",
    size=9.2, weight="bold", color="#111111")


def scatter(rect, xv, yv, colour, fit_colour, xlabel, solid_fit=True):
    ax = fig.add_axes(rect)
    rx, ry = stats.rankdata(xv), stats.rankdata(yv)
    xs = np.linspace(rx.min(), rx.max(), 10)
    sl, ic = np.polyfit(rx, ry, 1)
    ax.plot(xs, sl * xs + ic, color=fit_colour, lw=1.6 if solid_fit else 1.3,
            ls="-" if solid_fit else (0, (4, 3)), zorder=2)
    ax.plot(rx, ry, "o", color=colour, ms=5, mec="#FFFFFF", mew=0.7, alpha=0.95, zorder=3)
    ax.set_xlabel(xlabel, fontsize=6.6, labelpad=2.0)
    ax.set_ylabel("$TNFRSF9^+$ % of CD8 (rank)", fontsize=6.6, labelpad=2.0)
    ax.set_xticks([1, 5, 10, 15, 19])
    ax.set_yticks([1, 5, 10, 15, 19])
    ax.set_xlim(-0.5, 20.5)
    ax.set_ylim(-0.5, 20.5)
    tidy(ax, size=5.8)
    return ax


scatter([0.112, 0.288, 0.310, 0.166], per.epithelial_APM_z_5gene_primary.values,
        per.TNFRSF9pos_pct_of_CD8.values, GREEN, GREEN,
        "epithelial antigen-presentation score (rank)")
say(0.267, 0.234, f"$r$ = {prin.spearman_r:+.3f}, $q$ = {prin.bh_q:.3f}",
    size=8.0, weight="bold", color=GREEN)
say(0.267, 0.208, "antigen-presentation transcription tracks it", size=6.9, color=GREEN)

scatter([0.612, 0.288, 0.310, 0.166], per.DLL3_mean.values,
        per.TNFRSF9pos_pct_of_CD8.values, GREY, "#AAAAAA",
        "tumour DLL3 (rank)", solid_fit=False)
say(0.767, 0.234, f"$r$ = {null_prin.spearman_r:+.3f}, $P$ = {null_prin.p:.2f}",
    size=8.0, weight="bold", color=RED)
say(0.767, 0.206, f"DLL3 does not: equivalence only within\n"
                  f"|ρ| = {null.smallest_equivalence_margin.min():.3f} to "
                  f"{null.smallest_equivalence_margin.max():.3f}", size=6.9, color=RED)

# ============================================================ take-home
box(0.030, 0.052, 0.970, 0.152, PALE_GREEN, GREEN)
say(0.500, 0.120, "Baseline DLL3 showed no detectable association with this CD8 phenotype.",
    size=9.2, weight="bold", color="#14402A")
say(0.500, 0.082, "Antigen-presentation transcription did, and the two immune features are not "
                  "independent axes.", size=8.4, color="#14402A")

say(0.500, 0.018, "bulk RNA-seq 160 external + 15 institutional tumours  ·  single-cell 19 donors, 23 biospecimens  "
                  "·  multiplex immunofluorescence 10 patients, 0 of 51 phenotypes significant",
    size=6.5, color=GREY)

for ext in ("pdf", "png"):
    fig.savefig(os.path.join(FIG, f"Graphical_Abstract_v3.{ext}"), dpi=600)
plt.close(fig)
print("  Graphical_Abstract_v3.{pdf,png}")
