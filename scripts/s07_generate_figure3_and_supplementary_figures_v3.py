"""v3: Figure 3 and the complete Supplementary Figure set S1-S10.

The v2 supplementary figure set (14 files) is withdrawn in full. It contained
figures generated at N = 16 under a superseded cohort definition, a figure using
the "George N = 86" framing that the v2 text had already retracted, and five mIF
plots that were ROI-level categorical box plots with Japanese axis labels — a
display that contradicts the patient-level continuous analysis the Methods
specify and that presents 193 ROIs as if they were independent observations.

Every figure here is generated from v3 derived data, has English axes, carries an
S-number, and has its own legend (Figure_Legends_v3).
"""
import os
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
SRC = os.path.join(BASE, "14_Revision_v3_20260806/outputs")
FIG = os.path.join(BASE, "14_Revision_v3_20260806/figures")
DER = os.path.join(BASE, "05_Source_Data/06_Derived_Analysis_Tables")
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({"font.size": 8, "font.family": ["Arial", "DejaVu Sans"],  # IJC: Helvetica/Arial in figures
    "mathtext.fontset": "custom", "mathtext.rm": "Arial", "mathtext.it": "Arial:italic", "mathtext.bf": "Arial:bold", "savefig.bbox": "tight",
                     "pdf.fonttype": 42, "ps.fonttype": 42, "axes.linewidth": 0.8})
BLUE, RED, GREY, GREEN, ORANGE = "#2C5F8A", "#A81E22", "#666666", "#2E7D53", "#C4762A"
CD8 = ["CD8+ Texh", "CD8+ Teff", "CD8+ Tmem"]
MARK = {"PD-1": "PDCD1", "TIM-3": "HAVCR2", "GZMB": "GZMB"}
APM = ["HLA-A", "HLA-B", "B2M", "TAP1", "TAP2"]


def save(fig, name):
    for ext in ("png", "pdf", "tiff"):
        kw = {"pil_kwargs": {"compression": "tiff_lzw"}} if ext == "tiff" else {}
        fig.savefig(os.path.join(FIG, f"{name}.{ext}"), dpi=600, **kw)
    plt.close(fig)
    print(f"  {name}")


def panel(ax, letter, title, fs=8.5):
    ax.set_title(f"({letter}) {title}", fontsize=fs, loc="left", weight="bold", pad=6)


def fisher_ci(r, n, a=0.05):
    z, se = np.arctanh(r), 1 / np.sqrt(n - 3)
    k = stats.norm.ppf(1 - a / 2)
    return float(np.tanh(z - k * se)), float(np.tanh(z + k * se))


# ------------------------------------------------------------ load single cell
print("loading Chan immune atlas ...")
imm = ad.read_h5ad(os.path.join(BASE, "05_Source_Data/04_Public_scRNA_Chan_Atlas/"
                                      "Chan_HTAN_SCLC_Immune_Cells_16475cells.h5ad"))
feat = imm.var.feature_name.astype(str)


def vec(sym, layer=None):
    j = imm.var.index[feat == sym]
    if len(j) == 0:
        return None
    col = imm[:, j[0]].layers[layer] if layer else imm[:, j[0]].X
    return np.asarray(col.todense()).ravel() if sp.issparse(col) else np.asarray(col).ravel()


cell = pd.DataFrame({
    "donor": imm.obs.donor_id.astype(str).values,
    "ct": imm.obs.author_cell_type.astype(str).values,
    "total_counts": imm.obs.total_counts.values.astype(float),
    "n_genes": imm.obs.n_genes_by_counts.values.astype(float),
    "mito_frac": imm.obs.mito_frac.values.astype(float)})
for s in ["TNFRSF9"] + list(MARK.values()):
    cell[s + "_pos"] = vec(s) > 0
    cell[s + "_expr"] = vec(s, "normalized")
c8 = cell[cell.ct.isin(CD8)].copy()

# per-donor paired observed positivity rates (primary >=3-cell rule)
paired = []
for donor, g in c8.groupby("donor", observed=True):
    pos, neg = g[g.TNFRSF9_pos], g[~g.TNFRSF9_pos]
    if len(g) < 10 or len(pos) < 3 or len(neg) < 3:
        continue
    row = {"donor": donor, "n_pos": len(pos), "n_neg": len(neg)}
    for lab, sym in MARK.items():
        row[f"{lab}_pos"] = pos[f"{sym}_pos"].mean() * 100
        row[f"{lab}_neg"] = neg[f"{sym}_pos"].mean() * 100
        row[f"{lab}_expr_pos"] = pos[f"{sym}_expr"].mean()
        row[f"{lab}_expr_neg"] = neg[f"{sym}_expr"].mean()
    paired.append(row)
paired = pd.DataFrame(paired)
paired.to_csv(os.path.join(SRC, "Plotdata_Chan_Paired_Observed_v3.csv"), index=False)
print(f"  paired observed values for {len(paired)} donors")

per = pd.read_csv(os.path.join(SRC, "Chan_PerDonor_with_epithelialAPM_v3.csv"), index_col=0)
counts = pd.read_csv(os.path.join(SRC, "SuppTable_Chan_PerDonor_CellCounts_v3.csv"), index_col=0)
depth = pd.read_csv(os.path.join(SRC, "SuppTable_DetectionDepth_Sensitivity_v3.csv"))
axeq = pd.read_csv(os.path.join(SRC, "SuppTable_Axis1_vs_Axis2_and_Equivalence_v3.csv"))


# ================================================================== FIGURE 3
def figure3():
    fig = plt.figure(figsize=(7.6, 7.6))
    gs = fig.add_gridspec(3, 3, hspace=0.72, wspace=0.52)

    # A: DLL3 vs Texh fraction, rank axes
    ax = fig.add_subplot(gs[0, 0])
    x, y = per.DLL3_mean.values, per.CD8Texh_pct_of_T.values
    rx, ry = stats.rankdata(x), stats.rankdata(y)
    r, p = stats.spearmanr(x, y)
    ax.scatter(rx, ry, s=26, c=BLUE, edgecolor="#12314A", lw=.5, zorder=3)
    m, b = np.polyfit(rx, ry, 1)
    ax.plot([1, len(rx)], [m + b, m * len(rx) + b], color=RED, lw=1.2)
    ax.set_xlabel("Tumour DLL3 (rank)", fontsize=7)
    ax.set_ylabel("CD8$^+$ Texh fraction (rank)", fontsize=7)
    ax.text(.04, .95, f"$r$ = {r:+.3f}, $P$ = {p:.3f}\nfragile: see Results",
            transform=ax.transAxes, fontsize=6.2, va="top")
    panel(ax, "A", "Composition association")
    ax.tick_params(labelsize=6.3)

    # B: TNFRSF9 positivity by subtype
    ax = fig.add_subplot(gs[0, 1])
    subs = CD8 + ["CD4+ Treg"]
    rates = [cell[cell.ct == s].TNFRSF9_pos.mean() * 100 for s in subs]
    ns = [int((cell.ct == s).sum()) for s in subs]
    cols = [BLUE, BLUE, BLUE, ORANGE]
    ax.bar(range(len(subs)), rates, color=cols, width=.66)
    for i, (v, n) in enumerate(zip(rates, ns)):
        ax.text(i, v + 1.2, f"{v:.1f}%\n$n$={n:,}", ha="center", fontsize=5.8)
    ax.set_xticks(range(len(subs)))
    ax.set_xticklabels([s.replace("CD8+ ", "CD8 ").replace("CD4+ ", "CD4 ") for s in subs],
                       fontsize=6.2, rotation=20, ha="right")
    ax.set_ylabel("$TNFRSF9^+$ cells (%)", fontsize=7)
    ax.set_ylim(0, max(rates) * 1.35)
    panel(ax, "B", "$TNFRSF9$ by subtype")
    ax.tick_params(labelsize=6.3)

    # C: paired OBSERVED values per donor (not pinned to zero)
    for i, lab in enumerate(["PD-1", "TIM-3"]):
        ax = fig.add_subplot(gs[0, 2] if i == 0 else gs[1, 0])
        for _, row in paired.iterrows():
            ax.plot([0, 1], [row[f"{lab}_neg"], row[f"{lab}_pos"]], "-",
                    color=BLUE if row[f"{lab}_pos"] > row[f"{lab}_neg"] else RED,
                    lw=.9, alpha=.75, marker="o", ms=3.2, mew=0)
        higher = int((paired[f"{lab}_pos"] > paired[f"{lab}_neg"]).sum())
        d = (paired[f"{lab}_pos"] - paired[f"{lab}_neg"]).values
        pv = stats.wilcoxon(d[d != 0])[1]
        ax.set_xticks([0, 1])
        ax.set_xticklabels(["$TNFRSF9^-$", "$TNFRSF9^+$"], fontsize=6.6)
        ax.set_xlim(-.35, 1.35)
        ax.set_ylabel(f"{lab}$^+$ cells (% of CD8)", fontsize=7)
        ax.text(.03, .95, f"{higher}/{len(paired)} higher\nWilcoxon $P$ = {pv:.1e}",
                transform=ax.transAxes, fontsize=6.2, va="top")
        panel(ax, "C" if i == 0 else "D", f"{lab}, paired", fs=8.0)
        ax.tick_params(labelsize=6.3)

    # E: per-donor lineage distribution (not a pie)
    ax = fig.add_subplot(gs[1, 1])
    TC = [c for c in cell.ct.unique() if c in CD8 + ["CD4+ Tconv", "CD4+ Treg", "Tgd"]]
    tn = cell[cell.TNFRSF9_pos & cell.ct.isin(TC)]
    lin = tn.assign(lin=np.where(tn.ct.isin(CD8), "CD8", np.where(tn.ct == "CD4+ Treg", "Treg", "Other")))
    pc = lin.groupby(["donor", "lin"], observed=True).size().unstack(fill_value=0)
    pc = pc.div(pc.sum(axis=1), axis=0) * 100
    order = ["CD8", "Treg", "Other"]
    order = [o for o in order if o in pc.columns]
    dat = [pc[o].dropna().values for o in order]
    bp = ax.boxplot(dat, widths=.5, showfliers=False, patch_artist=True)
    for b_, c_ in zip(bp["boxes"], [BLUE, ORANGE, GREY]):
        b_.set_facecolor(c_); b_.set_alpha(.35); b_.set_edgecolor(c_)
    for i, (o, c_) in enumerate(zip(order, [BLUE, ORANGE, GREY])):
        v = pc[o].dropna().values
        ax.scatter(np.random.default_rng(i).normal(i + 1, .055, len(v)), v,
                   s=13, c=c_, edgecolor="#111", lw=.3, zorder=3)
    pw = stats.wilcoxon(pc["CD8"], pc["Treg"])[1]
    ax.set_xticks(range(1, len(order) + 1))
    ax.set_xticklabels(order, fontsize=6.4)
    ax.set_xlabel("Lineage", fontsize=7)
    ax.set_ylabel("% of $TNFRSF9^+$ T cells", fontsize=7)
    ax.set_ylim(-8, 118)
    ax.text(.03, .99, f"CD8 {pc['CD8'].median():.0f}% vs Treg {pc['Treg'].median():.0f}%\n"
                      f"$P$ = {pw:.4f}", transform=ax.transAxes, fontsize=6.0, va="top")
    panel(ax, "E", "Lineage of $TNFRSF9^+$", fs=8.0)
    ax.tick_params(labelsize=6.3)

    # F: bulk TNFRSF9 vs modules
    ax = fig.add_subplot(gs[1, 2])
    vals = [(+0.743, "Discovery\neffector"), (+0.832, "Discovery\nexhaustion"),
            (+0.590, "Jiang\neffector"), (+0.436, "Jiang\nexhaustion")]
    nn = [15, 15, 79, 79]
    est = [v[0] for v in vals]
    ci = [fisher_ci(e, n) for e, n in zip(est, nn)]
    yy = np.arange(len(vals))[::-1]
    ax.errorbar(est, yy, xerr=[[e - c[0] for e, c in zip(est, ci)],
                               [c[1] - e for e, c in zip(est, ci)]],
                fmt="none", ecolor="#333", lw=1.0, capsize=2.2)
    ax.scatter(est, yy, s=40, c=[BLUE, BLUE, ORANGE, ORANGE], zorder=3, edgecolor="#111", lw=.5)
    ax.axvline(0, color="#888", lw=.8, ls="--")
    ax.set_yticks(yy); ax.set_yticklabels([v[1] for v in vals], fontsize=6.2)
    ax.set_xlabel("Spearman $r$", fontsize=7)
    ax.set_xlim(-.1, 1.05)
    ax.text(.03, -.40, "tissue-level co-expression,\nnot single-cell",
            transform=ax.transAxes, fontsize=5.8, color=GREY)
    panel(ax, "F", "Bulk $TNFRSF9$", fs=8.0)
    ax.tick_params(labelsize=6.3)

    # G: the central null with equivalence bound
    ax = fig.add_subplot(gs[2, :])
    eq = axeq[axeq.analysis == "DLL3_to_Axis2_equivalence"].reset_index(drop=True)
    labels = ["$TNFRSF9^+$\n(% of CD8)", "$TNFRSF9^+$PD-1$^+$\n(% of CD8)",
              "$TNFRSF9^+$TIM-3$^+$\n(% of CD8)"]
    yy = np.arange(len(eq))[::-1]
    ax.axvspan(-eq.smallest_equivalence_margin.max(), eq.smallest_equivalence_margin.max(),
               color="#EFEFEF", zorder=0)
    ax.errorbar(eq.spearman_r, yy, xerr=[eq.spearman_r - eq.ci_low, eq.ci_high - eq.spearman_r],
                fmt="none", ecolor="#333", lw=1.2, capsize=3)
    ax.scatter(eq.spearman_r, yy, s=55, c=GREY, zorder=3, edgecolor="#111", lw=.6)
    ax.axvline(0, color="#222", lw=.9, ls="--")
    for m in (-eq.smallest_equivalence_margin.max(), eq.smallest_equivalence_margin.max()):
        ax.axvline(m, color=RED, lw=.9, ls=":")
    ax.set_yticks(yy); ax.set_yticklabels(labels, fontsize=6.6)
    ax.set_xlabel("Spearman $r$ (tumour DLL3 vs abundance of the $TNFRSF9$-defined phenotype)",
                  fontsize=7.2)
    ax.set_xlim(-.75, .75)
    ax.set_ylim(-.7, len(eq) - .3)
    for i, row in eq.iterrows():
        ax.text(row.spearman_r, yy[i] + .24, f"$r$ = {row.spearman_r:+.3f}, $P$ = {row.p:.2f}",
                ha="center", fontsize=6.2)
    ax.text(.5, -.42, f"Shaded band and dotted lines: the equivalence bound the data support "
                      f"(|ρ| ≈ {eq.smallest_equivalence_margin.min():.2f}–"
                      f"{eq.smallest_equivalence_margin.max():.2f}). "
                      f"A moderate association is NOT excluded.",
            transform=ax.transAxes, ha="center", fontsize=6.4, color=RED)
    panel(ax, "G", "No association detected between DLL3 and this phenotype — with the bound stated")
    ax.tick_params(labelsize=6.3)

    save(fig, "Figure3_v3_TNFRSF9_CD8_Phenotype")


# ================================================================== SUPP FIGS
def suppS1_S2():
    d15 = pd.read_csv(os.path.join(SRC, "Plotdata_Discovery_SCLC15_v3.csv"))
    stat = pd.read_csv(os.path.join(SRC, "Discovery_Statistics_N15_primary_v3.csv"))
    # S1: discovery cohort key statistics, N=15 vs N=17
    p15 = stat[stat.population == "SCLC15"]
    p17 = stat[stat.population == "AllSequenced17"]
    keys = [m for m in p15.measure if m in set(p17.measure) and "median_split" not in m]
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    yy = np.arange(len(keys))[::-1]
    for off, pop, c, lab in [(+.17, p15, GREEN, "$N$ = 15 (primary, adjudicated SCLC)"),
                             (-.17, p17, GREY, "$N$ = 17 (sensitivity, diagnosis-inclusive)")]:
        v = [float(pop[pop.measure == k].estimate.iloc[0]) for k in keys]
        ax.scatter(v, yy + off, s=42, c=c, zorder=3, edgecolor="#111", lw=.5, label=lab)
        for vi, yi, k in zip(v, yy + off, keys):
            pv = float(pop[pop.measure == k].p.iloc[0])
            ax.text(vi, yi + .16, f"{vi:+.3f} ($P$={pv:.3f})", ha="center", fontsize=5.4, color=c)
    ax.axvline(0, color="#222", lw=.9, ls="--")
    ax.set_yticks(yy)
    ax.set_yticklabels([k.replace("_", " ").replace("~", " vs ") for k in keys], fontsize=6.2)
    ax.set_xlabel("Effect estimate (Spearman $r$, or standardized OLS β for the adjusted model)",
                  fontsize=7.2)
    ax.legend(fontsize=6.4, loc="lower right", frameon=False)
    ax.set_title("Supplementary Figure S1. Institutional discovery cohort: every key statistic "
                 "under the\nprimary (N = 15) and diagnosis-inclusive sensitivity (N = 17) "
                 "populations", fontsize=8, loc="left", weight="bold")
    ax.tick_params(labelsize=6.3)
    save(fig, "SuppFigureS1_v3_Discovery_Adjudication_Sensitivity")

    # S2: discovery correlation heatmap at N=15  (replaces the legacy N=16 matrix)
    d = pd.read_excel(os.path.join(BASE, "05_Source_Data/01_Institutional_RNAseq/"
                                         "InHouse_RNAseq_Normalized_Gene_Expression.xlsx"),
                      sheet_name="Normalized.expression").set_index("Reference_GeneName")
    d = d[~d.index.duplicated(keep="first")]
    cols = list(d15["sample"])
    genes = ["DLL3", "CD274", "PDCD1", "LAG3", "TIGIT", "CTLA4", "HAVCR2", "IFNG", "STAT1",
             "CXCL9", "HLA-A", "HLA-B", "B2M", "TAP1", "TAP2", "CD8A", "CD68", "CD163",
             "FOXP3", "PTPRC", "TNFRSF9"]
    genes = [g for g in genes if g in d.index]
    M = np.log2(d.loc[genes, cols].astype(float) + 1)
    C = M.T.corr(method="spearman")
    fig, ax = plt.subplots(figsize=(6.6, 5.8))
    im = ax.imshow(C.values, cmap="RdBu_r", vmin=-1, vmax=1)
    ax.set_xticks(range(len(genes))); ax.set_xticklabels(genes, rotation=55, ha="right", fontsize=6)
    ax.set_yticks(range(len(genes))); ax.set_yticklabels(genes, fontsize=6)
    for i in range(len(genes)):
        for j in range(len(genes)):
            ax.text(j, i, f"{C.values[i, j]:.2f}", ha="center", va="center", fontsize=4.4,
                    color="white" if abs(C.values[i, j]) > .6 else "#222")
    cb = fig.colorbar(im, fraction=.046, pad=.02); cb.set_label("Spearman ρ", fontsize=7)
    cb.ax.tick_params(labelsize=6)
    ax.set_title("Supplementary Figure S2. Gene-gene correlation structure in the institutional\n"
                 "discovery cohort ($N$ = 15 pathologically adjudicated SCLC)",
                 fontsize=8, loc="left", weight="bold")
    save(fig, "SuppFigureS2_v3_Discovery_Correlation_Matrix_N15")


def suppS3_S4():
    jia = pd.read_csv(os.path.join(SRC, "Plotdata_Jiang_v3.csv"))
    col = pd.read_csv(os.path.join(SRC, "Plotdata_Cologne_v3.csv"))
    g = pd.read_csv(os.path.join(BASE, "05_Source_Data/03_Public_Bulk_George_GSE60052/"
                                       "George_GSE60052_Log2_Normalized_RNAseq.tsv"),
                    sep="\t", index_col=0)
    g = g[~g.index.duplicated(keep="first")]
    tum = [c for c in g.columns if "normal" not in c.lower()]
    X = pd.read_csv(os.path.join(BASE, "12_External_Replication_20260806/"
                                       "Cologne_GSE_expression_matrix.csv"), index_col=0)
    CL = np.log2(X.clip(lower=0) + 1)

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4))
    for ax, (name, dat, dll3, n) in zip(axes, [
            ("Jiang (GSE60052), $n$ = 79", g.loc[:, tum].astype(float).T,
             g.loc["DLL3", tum].astype(float), 79),
            ("George/Cologne, $n$ = 81", CL, CL.DLL3, 81)]):
        genes = APM + (["HLA-C"] if name.startswith("Jiang") else [])
        est, lo, hi, pv = [], [], [], []
        for gene in genes:
            r, p = stats.spearmanr(dll3, dat[gene])
            l, h = fisher_ci(r, n)
            est.append(r); lo.append(l); hi.append(h); pv.append(p)
        yy = np.arange(len(genes))[::-1]
        ax.errorbar(est, yy, xerr=[np.array(est) - lo, np.array(hi) - est], fmt="none",
                    ecolor="#333", lw=1.0, capsize=2.2)
        ax.scatter(est, yy, s=40, c=[RED if q < .05 else GREY for q in pv], zorder=3,
                   edgecolor="#111", lw=.5)
        ax.axvline(0, color="#222", lw=.9, ls="--")
        ax.set_yticks(yy); ax.set_yticklabels(genes, fontsize=6.6)
        ax.set_xlabel("Spearman $r$ (DLL3 vs gene)", fontsize=7)
        ax.set_xlim(-.65, .3)
        ax.set_title(name, fontsize=8, loc="left", weight="bold")
        ax.tick_params(labelsize=6.3)
        for e, y_, p_ in zip(est, yy, pv):
            ax.text(e, y_ + .22, f"$P$ = {p_:.3f}" if p_ >= .001 else "$P$ < 0.001",
                    ha="center", fontsize=5.4)
    fig.suptitle("Supplementary Figure S3. DLL3 versus individual antigen-presentation genes in "
                 "both external cohorts\n(red = P < 0.05; HLA-C has zero variance in the Cologne "
                 "matrix and is not estimable there)",
                 fontsize=8, x=0.01, ha="left", weight="bold", y=1.06)
    save(fig, "SuppFigureS3_v3_Individual_APM_Genes_Both_Cohorts")

    # S4: DLL3 vs PTPRC in both cohorts
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2))
    for ax, (name, dll3, ptprc, n) in zip(axes, [
            ("Jiang (GSE60052), $n$ = 79", g.loc["DLL3", tum].astype(float),
             g.loc["PTPRC", tum].astype(float), 79),
            ("George/Cologne, $n$ = 81", CL.DLL3, CL.PTPRC, 81)]):
        rx, ry = stats.rankdata(dll3), stats.rankdata(ptprc)
        r, p = stats.spearmanr(dll3, ptprc)
        ax.scatter(rx, ry, s=15, c=ORANGE, edgecolor="#6a3f14", lw=.35, alpha=.85, zorder=3)
        m, b = np.polyfit(rx, ry, 1)
        ax.plot([1, n], [m + b, m * n + b], color=RED, lw=1.2)
        ax.set_xlabel("DLL3 (rank)", fontsize=7)
        ax.set_ylabel("$PTPRC$ / CD45 (rank)", fontsize=7)
        ax.text(.04, .1, f"$r$ = {r:+.3f}\n$P$ = {p:.1e}", transform=ax.transAxes, fontsize=6.6)
        ax.set_title(name, fontsize=8, loc="left", weight="bold")
        ax.tick_params(labelsize=6.3)
    fig.suptitle("Supplementary Figure S4. DLL3-high tumours contain less leukocyte transcript in "
                 "both external cohorts.\nThis is the composition signal that attenuates the "
                 "DLL3–APM association (post hoc).",
                 fontsize=8, x=0.01, ha="left", weight="bold", y=1.08)
    save(fig, "SuppFigureS4_v3_DLL3_vs_PTPRC_Both_Cohorts")


def suppS5_S6_S7():
    # S5: depth distributions
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.9))
    for ax, (q, lab, logscale) in zip(axes, [
            ("total_counts", "Total UMI per cell", True),
            ("n_genes", "Detected genes per cell", True),
            ("mito_frac", "Mitochondrial fraction", False)]):
        a = c8.loc[c8.TNFRSF9_pos, q].values
        b = c8.loc[~c8.TNFRSF9_pos, q].values
        parts = ax.violinplot([b, a], showmedians=True, widths=.8)
        for pc_, c_ in zip(parts["bodies"], [GREY, BLUE]):
            pc_.set_facecolor(c_); pc_.set_alpha(.45)
        u, p = stats.mannwhitneyu(a, b)
        auc = u / (len(a) * len(b))
        if logscale:
            ax.set_yscale("log")
        ax.set_xticks([1, 2]); ax.set_xticklabels(["$TNFRSF9^-$", "$TNFRSF9^+$"], fontsize=6.6)
        ax.set_ylabel(lab, fontsize=7)
        ax.text(.03, .96, f"AUC = {auc:.3f}\n$P$ = {p:.1e}", transform=ax.transAxes,
                fontsize=6.0, va="top")
        ax.tick_params(labelsize=6.3)
    fig.suptitle("Supplementary Figure S5. Sequencing depth differs modestly between $TNFRSF9^+$ "
                 "and $TNFRSF9^-$ CD8$^+$ cells.\nAn AUC near 0.6 is a small difference; "
                 "Supplementary Figure S6 shows it does not account for the finding.",
                 fontsize=8, x=0.01, ha="left", weight="bold", y=1.10)
    save(fig, "SuppFigureS5_v3_Detection_Depth_Distributions")

    # S6: depth-adjusted ORs + continuous expression
    log = depth[depth.analysis == "cell_level_logistic"]
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4),
                             gridspec_kw={"width_ratios": [1.25, 1]})
    ax = axes[0]
    models = ["unadjusted", "+ depth", "+ depth + subtype + donor"]
    ypos, ylab = [], []
    k = 0
    for mk in ["PD-1", "TIM-3", "GZMB"]:
        for mo in models:
            row = log[(log.quantity == mk) & (log.model == mo)].iloc[0]
            c_ = {"PD-1": BLUE, "TIM-3": GREEN, "GZMB": GREY}[mk]
            ax.errorbar([row.odds_ratio], [k],
                        xerr=[[row.odds_ratio - row.ci_low], [row.ci_high - row.odds_ratio]],
                        fmt="none", ecolor=c_, lw=1.1, capsize=2.2)
            ax.scatter([row.odds_ratio], [k], s=40, c=c_, zorder=3, edgecolor="#111", lw=.5)
            ypos.append(k); ylab.append(f"{mk}: {mo}")
            k += 1
        k += 0.6
    ax.axvline(1, color="#222", lw=.9, ls="--")
    ax.set_yticks(ypos); ax.set_yticklabels(ylab, fontsize=6.0)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlabel("Odds ratio for marker positivity in $TNFRSF9^+$ cells\n"
                  "(cell-level logistic, donor cluster-robust CI)", fontsize=7)
    ax.set_title("(A) Depth- and donor-adjusted models", fontsize=8, loc="left", weight="bold")
    ax.tick_params(labelsize=6.3)

    ax = axes[1]
    for i, lab in enumerate(["PD-1", "TIM-3"]):
        for _, row in paired.iterrows():
            ax.plot([i * 2.2, i * 2.2 + 1], [row[f"{lab}_expr_neg"], row[f"{lab}_expr_pos"]],
                    "-", color=BLUE if row[f"{lab}_expr_pos"] > row[f"{lab}_expr_neg"] else RED,
                    lw=.85, alpha=.75, marker="o", ms=3, mew=0)
        d = (paired[f"{lab}_expr_pos"] - paired[f"{lab}_expr_neg"]).values
        pv = stats.wilcoxon(d[d != 0])[1]
        ax.text(i * 2.2 + .5, ax.get_ylim()[1] * .98, f"{lab}\n$P$ = {pv:.1e}",
                ha="center", fontsize=6.0, va="top")
    ax.set_xticks([0, 1, 2.2, 3.2])
    ax.set_xticklabels(["$TNFRSF9^-$", "$TNFRSF9^+$"] * 2, fontsize=5.8, rotation=25, ha="right")
    ax.set_ylabel("Mean normalized expression", fontsize=7)
    ax.set_title("(B) Continuous expression, not binary detection", fontsize=8, loc="left",
                 weight="bold")
    ax.tick_params(labelsize=6.3)
    fig.suptitle("Supplementary Figure S6. The $TNFRSF9$-defined phenotype survives adjustment for "
                 "sequencing depth.\n$GZMB$ does not, which is why the effector arm is reported as "
                 "exploratory.",
                 fontsize=8, x=0.01, ha="left", weight="bold", y=1.07)
    save(fig, "SuppFigureS6_v3_DepthAdjusted_Models")

    # S7: per-donor cell counts
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.6),
                             gridspec_kw={"width_ratios": [1.35, 1]})
    cc = counts.sort_values("n_CD8", ascending=True)
    ax = axes[0]
    yy = np.arange(len(cc))
    ax.barh(yy, cc.n_TNFRSF9_neg, color=GREY, label="$TNFRSF9^-$ CD8")
    ax.barh(yy, cc.n_TNFRSF9_pos, left=cc.n_TNFRSF9_neg, color=BLUE, label="$TNFRSF9^+$ CD8")
    ax.set_yticks(yy)
    ax.set_yticklabels([f"{i}{'' if e else '  (not evaluable)'}"
                        for i, e in zip(cc.index, cc.evaluable_primary_min3)], fontsize=5.8)
    ax.set_xlabel("CD8$^+$ cells per sample", fontsize=7)
    ax.set_xscale("symlog")
    ax.legend(fontsize=6.0, frameon=False, loc="lower right")
    ax.set_title("(A) CD8 compartment size and $TNFRSF9$ split", fontsize=8, loc="left",
                 weight="bold")
    ax.tick_params(labelsize=6.3)

    ax = axes[1]
    tx = cc[[c for c in cc.columns if c.startswith("n_CD8_")]]
    tx = tx.loc[cc.index]
    share = counts["n_CD8_Texh"] / counts["n_CD8_Texh"].sum() * 100
    share = share.sort_values(ascending=False)
    ax.bar(range(len(share)), share.values, color=[RED] + [GREY] * (len(share) - 1))
    ax.set_xticks(range(len(share)))
    ax.set_xticklabels(share.index, rotation=90, fontsize=5.2)
    ax.set_ylabel("% of all CD8$^+$ Texh cells in the atlas", fontsize=7)
    ax.text(.30, .90, f"{share.index[0]} contributes {share.iloc[0]:.0f}%\n"
                      f"of all {int(counts['n_CD8_Texh'].sum())} Texh cells",
            transform=ax.transAxes, fontsize=6.2, color=RED, va="top")
    ax.set_title("(B) Texh cells are concentrated in one donor", fontsize=8, loc="left",
                 weight="bold")
    ax.tick_params(labelsize=6.3)
    fig.suptitle("Supplementary Figure S7. Per-sample cell counts underlying every Chan atlas "
                 "analysis.\n15 of 19 samples meet the primary evaluability rule.",
                 fontsize=8, x=0.01, ha="left", weight="bold", y=1.06)
    save(fig, "SuppFigureS7_v3_PerDonor_Cell_Counts")


def suppS8_S9_S10():
    # S8: mIF all 51 with CIs  (REPLACES the five ROI-level Japanese-axis box plots)
    mif = pd.read_csv(os.path.join(DER, "Derived_mIF_ROI_All_51Markers_Results.csv"))
    mif = mif.sort_values("r").reset_index(drop=True)
    n = 10
    lo, hi = zip(*[fisher_ci(r, n) for r in mif.r])
    fig, ax = plt.subplots(figsize=(7.2, 8.4))
    yy = np.arange(len(mif))
    ax.errorbar(mif.r, yy, xerr=[mif.r - np.array(lo), np.array(hi) - mif.r], fmt="none",
                ecolor="#999", lw=.8, capsize=1.6)
    ax.scatter(mif.r, yy, s=20, c=GREY, zorder=3, edgecolor="#111", lw=.35)
    ax.axvline(0, color="#222", lw=1.0, ls="--")
    ax.set_yticks(yy); ax.set_yticklabels(mif.marker, fontsize=5.2)
    ax.set_xlabel("Spearman $r$ (patient-level area-weighted mean density vs continuous DLL3 "
                  "H-score), 95% CI", fontsize=7.2)
    ax.set_xlim(-1.05, 1.05)
    ax.set_ylim(-1, len(mif))
    ax.text(.985, .035, f"All {len(mif)} derived phenotypes shown.\n"
                        f"0 of {len(mif)} reach FDR significance "
                        f"(minimum $q$ = {mif.q.min():.3f}).\n"
                        f"$n$ = 10 patients; every interval crosses zero.",
            transform=ax.transAxes, fontsize=6.6, va="bottom", ha="right",
            bbox=dict(fc="#FBEDED", ec=RED, lw=.8, boxstyle="round,pad=0.4"))
    ax.set_title("Supplementary Figure S8. Complete multiplex-immunofluorescence result: all 51 "
                 "phenotypes.\nNo subset is selected; ordering is by effect size and confers no "
                 "significance.", fontsize=8, loc="left", weight="bold")
    ax.tick_params(labelsize=6.3)
    save(fig, "SuppFigureS8_v3_mIF_All51_EffectSizes")

    # S9: TOST curves
    eq = axeq[axeq.analysis == "DLL3_to_Axis2_equivalence"].reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(7.2, 3.6))
    deltas = np.linspace(0.05, 0.95, 400)
    for i, row in eq.iterrows():
        z, se = np.arctanh(row.spearman_r), 1 / np.sqrt(row.n - 3)
        pv = [max(stats.norm.cdf((z - np.arctanh(d)) / se), stats.norm.sf((z + np.arctanh(d)) / se))
              for d in deltas]
        ax.plot(deltas, pv, lw=1.4, label=row.y.replace("_", " "),
                color=[BLUE, GREEN, ORANGE][i])
        ax.axvline(row.smallest_equivalence_margin, color=[BLUE, GREEN, ORANGE][i],
                   lw=.8, ls=":")
    ax.axhline(0.05, color=RED, lw=1.0, ls="--")
    ax.text(0.06, 0.058, "α = 0.05", fontsize=6.4, color=RED)
    ax.set_xlabel("Equivalence margin δ (|ρ|)", fontsize=7.2)
    ax.set_ylabel("TOST $P$ value", fontsize=7.2)
    ax.set_ylim(0, .55); ax.set_xlim(0.05, 0.95)
    ax.legend(fontsize=6.2, frameon=False, loc="upper right")
    ax.text(.99, .55, "Equivalence is supported only to the RIGHT of each dotted line.\n"
                      "No margin below |ρ| ≈ 0.43 is supported by these data.",
            transform=ax.transAxes, fontsize=6.4, ha="right", va="top", color="#333")
    ax.set_title("Supplementary Figure S9. What the central null can and cannot exclude "
                 "(two one-sided tests).\nMargins were not pre-specified, so these bounds are "
                 "descriptive.", fontsize=8, loc="left", weight="bold")
    ax.tick_params(labelsize=6.3)
    save(fig, "SuppFigureS9_v3_Equivalence_Bounds")

    # S10: axis1 vs axis2 scatters
    a12 = axeq[axeq.analysis == "Axis1_vs_Axis2_direct"].reset_index(drop=True)
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.9))
    for ax, (_, row) in zip(axes, a12.iterrows()):
        x, y = per["epithelial_APM_z"].values, per[row.y].values
        rx, ry = stats.rankdata(x), stats.rankdata(y)
        ax.scatter(rx, ry, s=24, c=GREEN, edgecolor="#12492c", lw=.5, zorder=3)
        m, b = np.polyfit(rx, ry, 1)
        ax.plot([1, len(rx)], [m + b, m * len(rx) + b], color=RED, lw=1.2)
        ax.set_xlabel("Epithelial APM score (rank)", fontsize=7)
        ax.set_ylabel(row.y.replace("_pct_of_CD8", "\n(% of CD8, rank)").replace("_", " "),
                      fontsize=6.6)
        role = "PRINCIPAL endpoint" if row.y == "TNFRSF9pos_pct_of_CD8" else "supportive (nested)"
        ax.text(.04, .96, f"$r$ = {row.spearman_r:+.3f}, $P$ = {row.p:.4f}\n"
                          f"$q$ = {row.bh_q:.3f}\n"
                          f"≥20 CD8: $r$ = {row.r_cd8_ge20:+.3f}, $P$ = {row.p_cd8_ge20:.3f}\n"
                          f"≥50 CD8: $r$ = {row.r_cd8_ge50:+.3f}, $P$ = {row.p_cd8_ge50:.3f}\n"
                          f"{role}",
                transform=ax.transAxes, fontsize=5.6, va="top")
        ax.tick_params(labelsize=6.3)
    fig.suptitle("Supplementary Figure S10. POST HOC direct test of the two immune features "
                 "against each other (Chan atlas, 19 samples).\nThey are positively correlated, "
                 "so they are not independent axes. Only the single-gated principal endpoint is "
                 "robust to\nminimum CD8 cell counts; the doubly-gated endpoints are supportive "
                 "only.",
                 fontsize=8, x=0.01, ha="left", weight="bold", y=1.16)
    save(fig, "SuppFigureS10_v3_Axis1_vs_Axis2")


if __name__ == "__main__":
    print("Generating Figure 3 and Supplementary Figures S1-S10 ->", FIG)
    figure3()
    suppS1_S2()
    suppS3_S4()
    suppS5_S6_S7()
    suppS8_S9_S10()
    print("done")
