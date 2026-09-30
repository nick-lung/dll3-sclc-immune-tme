"""Graphical Abstract for the IJC submission (author decision 2026-09-12, checklist item A10).

IJC 3.5 asks for a roughly square, colour image that summarises the key findings, carries only
labels and short phrases, and contains neither the manuscript title nor author names; it is read
alongside the "What's New" text. This builds it from the same derived outputs as Figure 4, so the
numbers cannot drift away from the manuscript.

Output: figures/Graphical_Abstract_v3.{pdf,png} — 180 x 180 mm, 600 dpi.
Run with the project interpreter (~/miniforge3/bin/python3; Python 3.9.2, SciPy 1.13.1).
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
SRC = os.path.join(BASE, "14_Revision_v3_20260806/outputs")
FIG = os.path.join(BASE, "14_Revision_v3_20260806/figures")

plt.rcParams.update({
    "font.size": 9, "font.family": ["Arial", "DejaVu Sans"],  # IJC: Helvetica/Arial in figures
    "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold", "savefig.bbox": "tight", "pdf.fonttype": 42, "ps.fonttype": 42,
})
BLUE, RED, GREY, GREEN = "#2C5F8A", "#A81E22", "#666666", "#2E7D53"
PALE_BLUE, PALE_GREEN, PALE_GREY = "#EAF1F7", "#EAF4EE", "#F2F2F2"

# ---- numbers, read from the derived outputs (never transcribed) ---------------------------
meta = pd.read_csv(os.path.join(SRC, "Meta_TwoExternalCohorts_v3.csv")).iloc[0]
ax12 = pd.read_csv(os.path.join(SRC, "SuppTable_Axis1_vs_Axis2_and_Equivalence_v3.csv"))
prin = ax12[(ax12.analysis == "Axis1_vs_Axis2_direct")
            & (ax12.y == "TNFRSF9pos_pct_of_CD8")].iloc[0]
null = ax12[ax12.analysis == "DLL3_to_Axis2_equivalence"]
ccle = pd.read_csv(os.path.join(SRC, "SuppTable_CCLE_LeukocyteFree_Falsification_v3.csv"))
ccle_apm = ccle[ccle.comparison.str.startswith("DLL3 ~ APM (5-gene")].iloc[0]
s9 = pd.read_csv(os.path.join(SRC, "SuppTable_External_Replication_and_Leukocyte_Sensitivity_v3.csv"))
ptprc = s9[s9.model == "std OLS APM~DLL3 + PTPRC"].set_index("cohort")

pooled = f"pooled $r$ = {meta.pooled_r:+.3f} (95% CI {meta.ci_low:+.3f} to {meta.ci_high:+.3f})".replace("-", "−")
null_lo, null_hi = null.spearman_r.min(), null.spearman_r.max()
bound_lo, null_bound = null.smallest_equivalence_margin.min(), null.smallest_equivalence_margin.max()

fig = plt.figure(figsize=(7.087, 7.087))          # 180 x 180 mm
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")


def box(xc, yc, w, h, face, edge, lw=1.4):
    ax.add_patch(FancyBboxPatch((xc - w / 2, yc - h / 2), w, h,
                                boxstyle="round,pad=0.012,rounding_size=0.02",
                                facecolor=face, edgecolor=edge, linewidth=lw, zorder=2))


def text(x, y, s, size=8.4, weight="normal", color="#111111", ha="center", va="center", style="normal"):
    ax.text(x, y, s, fontsize=size, fontweight=weight, color=color, ha=ha, va=va,
            style=style, zorder=3, linespacing=1.45)


# ---- exposure -----------------------------------------------------------------------------
box(0.5, 0.905, 0.40, 0.105, PALE_GREY, "#444444")
text(0.5, 0.928, "DLL3", size=14, weight="bold")
text(0.5, 0.884, "baseline tumour expression, small cell lung cancer", size=8.2, color="#333333")

# ---- connectors ---------------------------------------------------------------------------
ax.add_patch(FancyArrowPatch((0.40, 0.853), (0.255, 0.755), arrowstyle="-|>", mutation_scale=15,
                             linewidth=1.8, color=BLUE, shrinkA=0, shrinkB=0, zorder=2))
text(0.115, 0.826, "replicated in two\nexternal cohorts", size=7.8, color=BLUE, weight="bold")
text(0.115, 0.787, "160 tumours", size=7.4, color=BLUE)

ax.add_patch(FancyArrowPatch((0.60, 0.853), (0.745, 0.755), arrowstyle="-|>", mutation_scale=15,
                             linewidth=1.8, color=GREY, linestyle=(0, (4, 3)), shrinkA=0, shrinkB=0, zorder=2))
text(0.885, 0.826, "association\nNOT demonstrated", size=7.8, color=RED, weight="bold")
text(0.885, 0.780, f"$r$ = {null_lo:+.3f} to {null_hi:+.3f}\nall $P$ > 0.6", size=7.4, color=GREY)

# ---- the two features ---------------------------------------------------------------------
box(0.255, 0.585, 0.45, 0.30, PALE_BLUE, BLUE)
text(0.255, 0.700, "Antigen-presentation signal", size=9.6, weight="bold", color="#0F2F4A")
text(0.255, 0.648, "lower in DLL3-high tumours", size=8.4)
text(0.255, 0.606, pooled, size=7.8, color="#0F2F4A")
text(0.255, 0.556, "but the association is markedly\nattenuated by leukocyte content", size=7.8, color=GREY)
text(0.255, 0.505, "$PTPRC$-adjusted $P$ = "
     f"{ptprc.loc['Jiang (GSE60052)', 'p']:.2f} and "
     f"{ptprc.loc['George/Cologne (cBioPortal)', 'p']:.2f}", size=7.4, color=GREY)
_ccle_r = f"{ccle_apm.spearman_r:+.3f}".replace("-", "−")
text(0.255, 0.468, f"and absent in 59 leukocyte-free\nSCLC cell lines ($r$ = {_ccle_r})",
     size=7.4, color=GREY)

box(0.745, 0.585, 0.45, 0.30, PALE_BLUE, BLUE)
text(0.745, 0.700, "$TNFRSF9^+$ CD8 T cells", size=9.6, weight="bold", color="#0F2F4A")
text(0.745, 0.648, "PD-1/TIM-3-marker-high phenotype", size=8.4)
text(0.745, 0.606, "present in SCLC tumours", size=7.8, color="#0F2F4A")
text(0.745, 0.556, "survives adjustment for sequencing\ndepth, CD8 subtype and donor", size=7.8, color=GREY)
text(0.745, 0.495, "odds ratio 1.96 (PD-1), 3.04 (TIM-3)", size=7.4, color=GREY)

# ---- the two features covary ---------------------------------------------------------------
for x in (0.255, 0.745):
    ax.plot([x, x], [0.435, 0.375], color=GREEN, linewidth=1.4, zorder=2)
ax.add_patch(FancyArrowPatch((0.255, 0.375), (0.745, 0.375), arrowstyle="<|-|>", mutation_scale=15,
                             linewidth=1.8, color=GREEN, shrinkA=0, shrinkB=0, zorder=2))
box(0.5, 0.300, 0.72, 0.095, PALE_GREEN, GREEN, lw=1.2)
text(0.5, 0.325, "the two features covary with each other", size=8.6, weight="bold", color="#1B5638")
text(0.5, 0.281, f"epithelial antigen-presentation score vs $TNFRSF9^+$ % of CD8: "
                 f"$r$ = {prin.spearman_r:+.3f}, $q$ = {prin.bh_q:.3f}", size=7.8, color="#1B5638")

# ---- what the null does and does not exclude -----------------------------------------------
text(0.5, 0.205, f"equivalence is supported only within |ρ| ≈ {bound_lo:.2f} to {null_bound:.2f}: "
                 f"a moderate association is not excluded", size=7.6, color=RED)

# ---- take-home ------------------------------------------------------------------------------
box(0.5, 0.105, 0.94, 0.115, "#FFFFFF", "#444444", lw=1.2)
text(0.5, 0.131, "Baseline DLL3 does not identify tumours enriched for this CD8 phenotype.",
     size=9.0, weight="bold")
text(0.5, 0.084, "Antigen-presentation transcription tracks it; DLL3 does not.", size=9.0, weight="bold")

text(0.5, 0.020, "bulk RNA-seq 160 external + 15 institutional tumours  ·  single-cell 19 samples  ·  "
                 "multiplex immunofluorescence 10 patients", size=6.8, color=GREY)

for ext in ("pdf", "png"):
    fig.savefig(os.path.join(FIG, f"Graphical_Abstract_v3.{ext}"), dpi=600)
plt.close(fig)
print("  Graphical_Abstract_v3.{pdf,png}")
