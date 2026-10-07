"""v3: Regenerate main Figures 1-4 from v3 derived data.

Changes required by peer review and implemented here:
  Figure 1 - cohort flow at N = 15, with exclusion by diagnosis shown separately
             from exclusion by missing exposure measurement.
  Figure 2 - both external cohorts on one panel; panel F places every estimate on
             the SAME standardized-coefficient scale (unadjusted, +NE, +PTPRC,
             +NE+PTPRC) for both cohorts, instead of mixing Spearman r with OLS
             beta; a single interval method (Fisher z) throughout.
  Figure 3 - panel C plots observed paired values per donor rather than pinning
             the TNFRSF9-negative side to zero; panel D shows the per-donor
             distribution rather than a median pie chart.
  Figure 4 - the red cross is removed. DLL3-to-feature-2 is a dashed connector
             labelled "association not demonstrated" with the equivalence bound
             stated; the newly tested feature-1-to-feature-2 relationship is a
             solid positive connector.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from scipy import stats

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
SRC = os.path.join(BASE, "14_Revision_v3_20260806/outputs")
FIG = os.path.join(BASE, "14_Revision_v3_20260806/figures")
DER = os.path.join(BASE, "05_Source_Data/06_Derived_Analysis_Tables")
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({
    "font.size": 8, "font.family": ["Arial", "DejaVu Sans"],  # IJC: Helvetica/Arial in figures
    "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold", "savefig.bbox": "tight",
    "pdf.fonttype": 42, "ps.fonttype": 42, "axes.linewidth": 0.8,
    "xtick.major.width": 0.8, "ytick.major.width": 0.8,
})
BLUE, RED, GREY, GREEN = "#2C5F8A", "#A81E22", "#666666", "#2E7D53"


def save(fig, name):
    for ext in ("png", "pdf", "tiff"):
        kw = {"pil_kwargs": {"compression": "tiff_lzw"}} if ext == "tiff" else {}
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), dpi=600, **kw)
    plt.close(fig)
    print(f"  {name}.{{png,pdf,tiff}}")


def panel(ax, letter, title):
    ax.set_title(f"({letter}) {title}", fontsize=8.5, loc="left", weight="bold", pad=6)


def fisher_ci(r, n, a=0.05):
    z, se = np.arctanh(r), 1 / np.sqrt(n - 3)
    k = stats.norm.ppf(1 - a / 2)
    return float(np.tanh(z - k * se)), float(np.tanh(z + k * se))


# ================================================================ FIGURE 1
def figure1():
    fig = plt.figure(figsize=(7.2, 6.4))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.15, 1], hspace=0.32, wspace=0.22)

    def box(ax, x, y, w, h, text, fc="#EEF3F8", ec=BLUE, fs=6.6, weight="normal"):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.012",
                                    fc=fc, ec=ec, lw=0.9))
        ax.text(x + w / 2, y + h / 2, text, ha="center", va="center",
                fontsize=fs, weight=weight, wrap=True)

    def arrow(ax, x1, y1, x2, y2):
        ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                     mutation_scale=8, lw=0.9, color=GREY))

    # (A) RNA-seq flow
    ax = fig.add_subplot(gs[0, 0]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    panel(ax, "A", "Institutional bulk RNA-seq")
    box(ax, .08, .84, .84, .12, "17 tumour specimens sequenced\n(CellCarta P1947, all DV200 ≥ 40)")
    arrow(ax, .5, .83, .5, .77)
    box(ax, .08, .60, .84, .16,
        "Excluded on pathology adjudication: 2\n21029B2866, 21029T2952 → LC-NEC",
        fc="#FBEDED", ec=RED)
    arrow(ax, .5, .59, .5, .50)
    box(ax, .08, .36, .84, .13, "PRIMARY ANALYSIS SET\n15 pathologically confirmed SCLC",
        fc="#E6F0E9", ec=GREEN, weight="bold")
    arrow(ax, .5, .35, .5, .26)
    box(ax, .08, .12, .84, .13, "Subset with DLL3 IHC H-score: 10\n(= the mIF analysis cohort)")
    ax.text(.5, .025, "The 17-specimen set is retained\nas a sensitivity analysis",
            ha="center", va="center", fontsize=6.0, style="italic", color=GREY)

    # (B) mIF flow
    ax = fig.add_subplot(gs[0, 1]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    panel(ax, "B", "Institutional multiplex immunofluorescence")
    box(ax, .04, .84, .42, .12, "Archival batch\n15 patients / 372 ROIs", fs=6.3)
    box(ax, .54, .84, .42, .12, "Prospective batch\n16 patients / 260 ROIs", fs=6.3)
    arrow(ax, .25, .83, .25, .76)
    box(ax, .04, .48, .42, .27, "EXCLUDED in full\nno DLL3 H-score;\n6 of 15 are LC-NEC",
        fc="#FBEDED", ec=RED, fs=6.0)
    arrow(ax, .75, .83, .75, .76)
    box(ax, .54, .48, .42, .27,
        "Excluded by diagnosis: 1\n(21029B2866, LC-NEC)\n\nExcluded, no H-score: 5\n(67 ROIs)",
        fc="#FBEDED", ec=RED, fs=6.0)
    arrow(ax, .75, .47, .75, .36)
    box(ax, .12, .17, .82, .18, "mIF ANALYSIS SET\n10 patients / 193 ROIs\n(ROIs per patient 1–38, median 17)",
        fc="#E6F0E9", ec=GREEN, fs=6.3, weight="bold")
    ax.text(.5, .045, "All 10 are a subset of the 15 RNA-sequenced SCLC patients:\n"
                     "the two institutional analyses are not independent",
            ha="center", va="center", fontsize=6.0, style="italic", color=GREY)

    # (C) external datasets
    ax = fig.add_subplot(gs[1, :]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    panel(ax, "C", "External datasets and the analysis each supports")
    cards = [
        (.015, "Jiang cohort\nGSE60052", "79 tumours\n(7 normal lung controls\nexcluded from all\nprimary analyses)",
         "Feature 1:\nDLL3 vs bulk APM"),
        (.255, "George/Cologne\ncBioPortal", "81 tumours\n(sclc_ucologne_2015)",
         "Feature 1:\nindependent replication"),
        (.495, "Chan SCLC atlas\nCZ CELLxGENE", "19 samples\n54,313 epithelial +\n16,475 immune cells",
         "Feature 2, the DLL3 link,\nand the axis-vs-axis test"),
        (.735, "GSE319155\n(supportive)", "7 patients\nsingle-cell", "Supplementary\nconcordance only"),
    ]
    for x, title, detail, role in cards:
        ax.add_patch(FancyBboxPatch((x, .18), .225, .68, boxstyle="round,pad=0.012",
                                    fc="#F7F9FB", ec=BLUE, lw=0.9))
        ax.text(x + .1125, .78, title, ha="center", va="center", fontsize=6.8, weight="bold")
        ax.text(x + .1125, .55, detail, ha="center", va="center", fontsize=6.2)
        ax.text(x + .1125, .29, role, ha="center", va="center", fontsize=6.2,
                style="italic", color=BLUE)
    save(fig, "Figure1_v3_Study_Design_and_Cohort_Flow")


# ================================================================ FIGURE 2
def figure2():
    disc = pd.read_csv(os.path.join(SRC, "Plotdata_Discovery_SCLC15_v3.csv"))
    d17 = pd.read_csv(os.path.join(SRC, "Plotdata_Discovery_AllSequenced17_v3.csv"))
    jia = pd.read_csv(os.path.join(SRC, "Plotdata_Jiang_v3.csv"))
    col = pd.read_csv(os.path.join(SRC, "Plotdata_Cologne_v3.csv"))
    s9 = pd.read_csv(os.path.join(SRC, "SuppTable_External_Replication_and_Leukocyte_Sensitivity_v3.csv"))
    mt = pd.read_csv(os.path.join(SRC, "Meta_TwoExternalCohorts_v3.csv")).iloc[0]

    fig = plt.figure(figsize=(7.2, 7.4))
    # external cohorts first (A-D, cited first in Results), institutional cohort last (E-G)
    gs = fig.add_gridspec(3, 3, hspace=0.52, wspace=0.42,
                          height_ratios=[1, 1.05, 1])

    # (E) discovery scatter, rank axes
    ax = fig.add_subplot(gs[2, 0])
    rx, ry = stats.rankdata(disc.DLL3_log2), stats.rankdata(disc.CD274_APM_ratio)
    r, p = stats.spearmanr(disc.DLL3_log2, disc.CD274_APM_ratio)
    ax.scatter(rx, ry, s=26, c=BLUE, edgecolor="#12314A", lw=.5, zorder=3)
    m, b = np.polyfit(rx, ry, 1)
    ax.plot([1, len(rx)], [m + b, m * len(rx) + b], color=RED, lw=1.2)
    ax.set_xlabel("DLL3 mRNA (rank)", fontsize=7)
    ax.set_ylabel("CD274/APM ratio (rank)", fontsize=7)
    ax.text(.04, .93, f"$r$ = {r:+.3f}\n$P$ = {p:.3f}", transform=ax.transAxes,
            fontsize=6.6, va="top")
    panel(ax, "E", "Discovery, $N$ = 15 SCLC")
    ax.tick_params(labelsize=6.3)

    # (F) LOO
    ax = fig.add_subplot(gs[2, 1])
    loo = [stats.spearmanr(np.delete(disc.DLL3_log2.values, i),
                           np.delete(disc.CD274_APM_ratio.values, i)).statistic
           for i in range(len(disc))]
    ax.bar(range(1, len(loo) + 1), loo, color=BLUE, width=.72)
    ax.axhline(r, color=RED, lw=1.0, ls="--")
    ax.axhline(0, color="#222", lw=.8)
    ax.set_xlabel("Case deleted", fontsize=7)
    ax.set_ylabel("Spearman $r$", fontsize=7)
    ax.set_ylim(0, .8)
    ax.text(.04, .93, f"all 15 iterations positive\n({min(loo):+.3f} to {max(loo):+.3f})",
            transform=ax.transAxes, fontsize=6.2, va="top")
    panel(ax, "F", "Leave-one-out")
    ax.tick_params(labelsize=6.3)

    # (G) adjudication sensitivity
    ax = fig.add_subplot(gs[2, 2])
    r17, p17 = stats.spearmanr(d17.DLL3_log2, d17.CD274_APM_ratio)
    ests, los, his = [], [], []
    for rr, nn in [(r, 15), (r17, 17)]:
        lo, hi = fisher_ci(rr, nn); ests.append(rr); los.append(lo); his.append(hi)
    y = [1, 0]
    ax.errorbar(ests, y, xerr=[np.array(ests) - los, np.array(his) - ests], fmt="none",
                ecolor="#333", lw=1.1, capsize=2.5)
    ax.scatter(ests, y, s=[70, 50], c=[GREEN, GREY], zorder=3, edgecolor="#111", lw=.6)
    ax.axvline(0, color="#888", lw=.8, ls="--")
    ax.set_yticks(y)
    ax.set_yticklabels(["$N$ = 15\n(primary)", "$N$ = 17\n(sensitivity)"],
                       fontsize=6.0)
    ax.set_xlabel("Spearman $r$", fontsize=7)
    ax.set_ylim(-.6, 1.6)
    panel(ax, "G", "Pathology adjudication")
    ax.tick_params(labelsize=6.3)

    # (A) Jiang scatter  (B) Cologne scatter
    for i, (dat, lab, key) in enumerate([(jia, "Jiang (GSE60052), $n$ = 79", "APM_z_5gene"),
                                         (col, "George/Cologne, $n$ = 81", "APM_z_5gene")]):
        ax = fig.add_subplot(gs[0, i])
        x, yv = stats.rankdata(dat.DLL3_log2), stats.rankdata(dat[key])
        rr, pp = stats.spearmanr(dat.DLL3_log2, dat[key])
        ax.scatter(x, yv, s=13, c=BLUE, edgecolor="#12314A", lw=.35, alpha=.85, zorder=3)
        m, b = np.polyfit(x, yv, 1)
        ax.plot([1, len(x)], [m + b, m * len(x) + b], color=RED, lw=1.2)
        ax.set_xlabel("DLL3 (rank)", fontsize=7)
        ax.set_ylabel("APM module, 5-gene (rank)", fontsize=7)
        ax.text(.04, .12, f"$r$ = {rr:+.3f}\n$P$ = {pp:.1e}", transform=ax.transAxes, fontsize=6.6)
        panel(ax, "AB"[i], lab)
        ax.tick_params(labelsize=6.3)

    # (C) meta forest, one interval method
    ax = fig.add_subplot(gs[0, 2])
    rJ = stats.spearmanr(jia.DLL3_log2, jia.APM_z_5gene).statistic
    rC = stats.spearmanr(col.DLL3_log2, col.APM_z_5gene).statistic
    est = [rJ, rC, float(mt.pooled_r)]
    lo = [fisher_ci(rJ, 79)[0], fisher_ci(rC, 81)[0], float(mt.ci_low)]
    hi = [fisher_ci(rJ, 79)[1], fisher_ci(rC, 81)[1], float(mt.ci_high)]
    yy = [2, 1, 0]
    ax.errorbar(est, yy, xerr=[np.array(est) - lo, np.array(hi) - est], fmt="none",
                ecolor="#333", lw=1.1, capsize=2.5)
    ax.scatter(est[:2], yy[:2], s=48, c=BLUE, zorder=3, edgecolor="#111", lw=.6)
    ax.scatter(est[2:], yy[2:], s=80, c=RED, marker="D", zorder=3, edgecolor="#111", lw=.6)
    ax.axvline(0, color="#888", lw=.8, ls="--")
    ax.set_yticks(yy)
    ax.set_yticklabels(["Jiang\n($n$ = 79)", "Cologne\n($n$ = 81)", "Pooled\n($N$ = 160)"], fontsize=6.2)
    ax.set_xlabel("Spearman $r$ (DLL3 vs APM)", fontsize=7)
    ax.set_xlim(-.75, .05)
    ax.set_ylim(-.6, 2.6)
    ax.text(.03, -.30, "Fisher $z$ intervals throughout; heterogeneity not detectable ($k$ = 2)",
            transform=ax.transAxes, fontsize=5.8, color=GREY)
    panel(ax, "C", "External replication")
    ax.tick_params(labelsize=6.3)

    # (D) stepwise adjustment, ONE estimand scale
    ax = fig.add_subplot(gs[1, :])
    models = ["unadjusted", "+ NE lineage", "+ PTPRC", "+ NE lineage + PTPRC"]
    coh = {"Jiang (GSE60052)": BLUE, "George/Cologne (cBioPortal)": "#C4762A"}
    h = 0.34
    for ci, (cname, colr) in enumerate(coh.items()):
        sub = s9[(s9.cohort == cname) & s9.model.str.startswith("std OLS")]
        vals, los_, his_ = [], [], []
        for m in models:
            row = sub[sub.model == f"std OLS APM~DLL3 {m}"].iloc[0]
            vals.append(row.estimate); los_.append(row.ci_low_fisher); his_.append(row.ci_high_fisher)
        ypos = np.arange(len(models))[::-1] + (h / 2 if ci == 0 else -h / 2)
        ax.errorbar(vals, ypos, xerr=[np.array(vals) - los_, np.array(his_) - vals],
                    fmt="none", ecolor=colr, lw=1.1, capsize=2.5)
        ax.scatter(vals, ypos, s=46, c=colr, zorder=3, edgecolor="#111", lw=.5,
                   label=cname.split(" (")[0])
        for v, yp, pv in zip(vals, ypos, [sub[sub.model == f"std OLS APM~DLL3 {m}"].iloc[0].p
                                          for m in models]):
            ax.text(v, yp + 0.13, f"β = {v:+.3f}, $P$ = {pv:.3f}" if pv >= .001
                    else f"β = {v:+.3f}, $P$ < 0.001", ha="center", fontsize=5.6, color=colr)
    ax.axvline(0, color="#222", lw=.9, ls="--")
    ax.set_yticks(np.arange(len(models))[::-1])
    ax.set_yticklabels(models, fontsize=7)
    ax.set_xlabel("Standardized OLS coefficient for DLL3 (same scale for every row)", fontsize=7.2)
    ax.set_xlim(-0.85, 0.45)
    ax.legend(fontsize=6.4, loc="lower left", frameon=False)
    ax.text(.985, .06, "PTPRC models are post hoc", transform=ax.transAxes,
            fontsize=6.0, style="italic", color=GREY, ha="right")
    panel(ax, "D", "Stepwise multivariable adjustment in both cohorts, one estimand")
    ax.tick_params(labelsize=6.3)

    save(fig, "Figure2_v3_Antigen_Presentation_and_Leukocyte_Content")


# ================================================================ FIGURE 5
def figure5():
    ax_s = pd.read_csv(os.path.join(SRC, "SuppTable_Axis1_vs_Axis2_and_Equivalence_v3.csv"))
    a12 = ax_s[ax_s.analysis == "Axis1_vs_Axis2_direct"]
    a12_main = a12[a12.endpoint_role == "principal exploratory endpoint"].iloc[0]
    a12_support = a12[a12.endpoint_role == "supportive nested endpoint"]
    eqv = ax_s[ax_s.analysis == "DLL3_to_Axis2_equivalence"]
    apm_ctx = ax_s[(ax_s.analysis == "DLL3_to_epithelial_APM_context") &
                   (ax_s.endpoint_role == "primary 5-gene APM context")].iloc[0]

    fig, ax = plt.subplots(figsize=(7.2, 4.3))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")

    def node(x, y, w, h, title, body, fc, ec):
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.014",
                                    fc=fc, ec=ec, lw=1.1))
        ax.text(x + w / 2, y + h * 0.74, title, ha="center", va="center",
                fontsize=7.6, weight="bold")
        ax.text(x + w / 2, y + h * 0.34, body, ha="center", va="center", fontsize=6.3)

    node(.015, .56, .30, .28, "DLL3 expression",
         "continuous, all platforms", "#F2F2F2", "#444444")
    node(.40, .70, .56, .24, "Feature 1\nAntigen-presentation signal",
         "Bulk: lower with DLL3; replicated but $PTPRC$-sensitive\n"
         f"Chan epithelial five-gene APM: DLL3 $r$ = {apm_ctx.spearman_r:+.3f}, "
         f"$P$ = {apm_ctx.p:.2f}\nused for the within-Chan green connector", "#EEF3F8", BLUE)
    node(.40, .12, .56, .24, "Feature 2\n$TNFRSF9$-defined, PD-1/TIM-3-marker-high CD8",
         "present in SCLC; robust to sequencing depth\n"
         "(adjusted OR 1.96 for PD-1, 3.04 for TIM-3)", "#EEF3F8", BLUE)

    # DLL3 -> feature 1 : solid
    ax.add_patch(FancyArrowPatch((.325, .76), (.395, .82), arrowstyle="-|>",
                                 mutation_scale=10, lw=1.6, color=BLUE))
    ax.text(.352, .875, "bulk signal replicated\nin 2 cohorts", fontsize=6.0, color=BLUE,
            ha="center", va="center")

    # DLL3 -> feature 2 : dashed, association not demonstrated
    ax.add_patch(FancyArrowPatch((.325, .64), (.395, .28), arrowstyle="-|>",
                                 mutation_scale=10, lw=1.4, color=GREY,
                                 linestyle=(0, (4, 3))))
    ax.text(.045, .34,
            "association NOT demonstrated\n"
            f"$r$ = {eqv.spearman_r.min():+.3f} to {eqv.spearman_r.max():+.3f}, all $P$ > 0.6\n"
            f"equivalence supported only within\n"
            f"|ρ| = {eqv.smallest_equivalence_margin.min():.3f} to "
            f"{eqv.smallest_equivalence_margin.max():.3f}\n"
            "a moderate association remains possible",
            fontsize=6.2, color=GREY, ha="left", va="center")

    # feature 1 <-> feature 2 : solid positive, POST HOC
    ax.add_patch(FancyArrowPatch((.60, .70), (.60, .36), arrowstyle="<|-|>",
                                 mutation_scale=10, lw=1.7, color=GREEN))
    ax.text(.625, .53,
            "POST HOC, within Chan\n"
            "epithelial APM vs $TNFRSF9$+ %CD8\n"
            f"$r$ = {a12_main.spearman_r:+.3f}, BH $q$ = {a12_main.bh_q:.3f}\n"
            f"nested PD-1/TIM-3 endpoints: $q$ = "
            f"{a12_support.bh_q.min():.3f} to {a12_support.bh_q.max():.3f}\n"
            "the two features covary in this dataset",
            fontsize=5.9, color=GREEN, ha="left", va="center", weight="bold")

    ax.text(.5, .035,
            "These baseline cross-sectional data do not test the mechanism or merit of DLL3-directed plus checkpoint-directed combination therapy.",
            ha="center", fontsize=6.0, style="italic", color="#333333")
    save(fig, "Figure5_v3_What_Is_And_Is_Not_Established")


if __name__ == "__main__":
    print("Generating v3 main figures ->", FIG)
    figure1()
    figure2()
    figure5()
    print("done")
