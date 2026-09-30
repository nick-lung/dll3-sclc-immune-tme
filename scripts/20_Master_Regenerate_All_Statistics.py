#!/usr/bin/env python3
"""
Master analysis script — regenerates EVERY number reported in the manuscript,
figures and tables from the raw source data, into one versioned CSV.

Written 2026-07-25 for the IJC resubmission (v2). Addresses the root cause
identified in 11_Draft_Check_20260720/Root_Cause_Analysis_JP.md: main tables and
figures had no generator and drifted out of sync with the analysis.

ANALYSIS SPECIFICATION (fixed here, mirrored verbatim in the Methods):
  * Jiang GSE60052 PRIMARY population  = 79 tumours (7 adjacent-normal samples
    excluded; the GEO matrix is literally 79tumor.7normal). Normal-inclusive
    N=86 is reported only as a sensitivity analysis.
    ATTRIBUTION (corrected 2026-07-25): GSE60052 is the cohort of Jiang et al.,
    PLoS Genet 2016;12(4):e1005895 (99 Chinese SCLC patients; GEO series title
    "The RNAseq of 79 small cell lung cancer (sclc) and 7 normal control").
    It is NOT the George et al. Nature 2015 cohort, whose data are deposited in
    EGA (EGAS00001000925), not GEO. Earlier versions of this project called it
    "the George cohort", which was a misattribution.
  * Module score = mean of gene-wise z-scored log2 expression, z-scored WITHIN
    the analysis population. (Previously the Methods said z-score but the
    implementation used raw log2 means.)
  * Institutional discovery cohort PRIMARY = all 17 sequenced samples.
  * Chan paired tests PRIMARY = >=10 CD8 cells per donor AND >=3 cells in each
    TNFRSF9 subgroup (n=15). >=1 cell (n=18) is the sensitivity analysis.
    (Previously the main text quoted the >=1 row that Supplementary Table S6
    itself labelled "Sensitivity".)
  * mIF FDR family = all 51 derived phenotypes, one Benjamini-Hochberg pass.

Usage:  python3 20_Master_Regenerate_All_Statistics.py
Output: 05_Source_Data/06_Derived_Analysis_Tables/Derived_MASTER_Statistics_v2.csv
"""
import os
import sys
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SEED = 20260725
N_BOOT = 2000
rng = np.random.default_rng(SEED)
ROWS = []


def emit(dataset, n, measure, effect, value, p=None, ci_low=None, ci_high=None,
         extra="", role=""):
    ROWS.append(dict(dataset=dataset, n=n, measure=measure, effect=effect,
                     value=value, ci_low=ci_low, ci_high=ci_high, p=p,
                     sig=("" if p is None else ("SIG" if p < 0.05 else "NS")),
                     extra=extra, role=role))


def hr(title):
    print("\n" + "=" * 86)
    print(title)
    print("=" * 86)


def boot_spearman(x, y, n_boot=N_BOOT):
    x, y = np.asarray(x, float), np.asarray(y, float)
    idx = rng.integers(0, len(x), (n_boot, len(x)))
    vals = np.array([stats.spearmanr(x[i], y[i])[0] for i in idx])
    vals = vals[~np.isnan(vals)]
    return np.percentile(vals, [2.5, 97.5])


def loo_spearman(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    v = [stats.spearmanr(np.delete(x, i), np.delete(y, i))[0] for i in range(len(x))]
    return float(np.min(v)), float(np.max(v))


def fisher_ci(r, n, alpha=0.05):
    if n < 4:
        return np.nan, np.nan
    z, se = np.arctanh(r), 1 / np.sqrt(n - 3)
    k = stats.norm.ppf(1 - alpha / 2)
    return float(np.tanh(z - k * se)), float(np.tanh(z + k * se))


def zmodule(df, genes, cols):
    """Mean of gene-wise z-scored expression, z-scored within `cols`."""
    g = [x for x in genes if x in df.index]
    s = df.loc[g, cols].astype(float)
    return s.sub(s.mean(axis=1), axis=0).div(s.std(axis=1, ddof=1), axis=0).mean(axis=0)


def std_ols(y, X):
    """Standardised OLS; returns (beta, ci_lo, ci_hi, p, r2) for the first column."""
    Xs = (X - X.mean()) / X.std(ddof=1)
    ys = (y - y.mean()) / y.std(ddof=1)
    m = sm.OLS(ys, sm.add_constant(Xs)).fit()
    k = X.columns[0]
    lo, hi = m.conf_int().loc[k]
    return float(m.params[k]), float(lo), float(hi), float(m.pvalues[k]), float(m.rsquared)


APM = ["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2"]
M2 = ["MRC1", "MARCO", "MSR1", "CD163"]
NE4 = ["ASCL1", "NEUROD1", "POU2F3", "YAP1"]
EFF = ["GZMB", "PRF1", "NKG7", "IFNG"]
EXH = ["PDCD1", "LAG3", "TIGIT", "HAVCR2", "CTLA4", "TOX"]

# ===========================================================================
# 1. Jiang GSE60052
# ===========================================================================
hr("1. JIANG GSE60052  (primary = 79 tumours; 86 incl. 7 normals = sensitivity)")

# NB: the archived directory/file still carry the legacy "George_" prefix from
# before the attribution was corrected; the contents are GSE60052 (Jiang et al.).
g = pd.read_csv(os.path.join(BASE, "05_Source_Data/03_Public_Bulk_George_GSE60052/"
                                   "George_GSE60052_Log2_Normalized_RNAseq.tsv"),
                sep="\t", index_col=0)
g.columns = [c.strip() for c in g.columns]
g = g[~g.index.duplicated(keep="first")]
NORMALS = [c for c in g.columns if c.lower().endswith(".normal")]
TUMORS = [c for c in g.columns if c not in NORMALS]
print(f"total {g.shape[1]} samples = {len(TUMORS)} tumours + {len(NORMALS)} adjacent normals")
emit("Jiang", 86, "cohort_composition", "n_tumor/n_normal", f"{len(TUMORS)}/{len(NORMALS)}",
     extra="GEO file GSE60052_79tumor.7normal.normalized.log2.data", role="audit")

for pop_label, cols, role in [("Jiang_tumor79", TUMORS, "PRIMARY"),
                              ("Jiang_all86", list(g.columns), "sensitivity")]:
    dll3 = g.loc["DLL3", cols].astype(float)
    apm = zmodule(g, APM, cols)
    m2 = zmodule(g, M2, cols)
    n = len(cols)

    r, p = stats.spearmanr(dll3, apm)
    lo, hi = boot_spearman(dll3, apm)
    l0, l1 = loo_spearman(dll3, apm)
    print(f"\n[{role}] {pop_label} (n={n})")
    print(f"  DLL3~APM (z-score module)  r={r:+.4f} P={p:.5g}  boot95%[{lo:+.3f},{hi:+.3f}]  LOO[{l0:+.3f},{l1:+.3f}]")
    emit(pop_label, n, "DLL3~APM_module_z", "Spearman_r", round(r, 4), p, round(lo, 4), round(hi, 4),
         extra=f"LOO_min={l0:.4f};LOO_max={l1:.4f}", role=role)

    r2_, p2_ = stats.spearmanr(dll3, g.loc[APM, cols].astype(float).mean(axis=0))
    emit(pop_label, n, "DLL3~APM_raw_log2_mean", "Spearman_r", round(r2_, 4), p2_,
         extra="legacy definition used in v1 manuscript", role="audit")
    print(f"  DLL3~APM (raw log2 mean, legacy)  r={r2_:+.4f} P={p2_:.5g}")

    for gene in APM + ["CD274", "TNFRSF9"] + M2:
        rr, pp = stats.spearmanr(dll3, g.loc[gene, cols].astype(float))
        emit(pop_label, n, f"DLL3~{gene}", "Spearman_r", round(rr, 4), pp, role=role)
        if pop_label == "Jiang_tumor79":
            print(f"    DLL3~{gene:<8} r={rr:+.4f} P={pp:.4g} {'' if pp < 0.05 else '(NS)'}")

    rm, pm = stats.spearmanr(dll3, m2)
    emit(pop_label, n, "DLL3~M2_module_z", "Spearman_r", round(rm, 4), pm, role=role)
    print(f"  DLL3~M2 module  r={rm:+.4f} P={pm:.5g}")

    # NE-lineage-adjusted models — CI and P are mandatory columns from now on
    for yname, yv in [("APM", apm), ("M2", m2)]:
        X = g.loc[["DLL3"] + NE4, cols].T.astype(float)
        b, blo, bhi, bp, r2 = std_ols(yv.astype(float), X)
        print(f"  {yname} ~ DLL3 + NE4 :  beta={b:+.4f} [{blo:+.4f},{bhi:+.4f}] P={bp:.4f} R2={r2:.4f}"
              f"  {'*** NOT SIGNIFICANT' if bp >= .05 else ''}")
        emit(pop_label, n, f"{yname}~DLL3+NE4", "std_OLS_beta_DLL3", round(b, 4), bp,
             round(blo, 4), round(bhi, 4), extra=f"R2={r2:.4f}", role=role)

    for gene in ["TNFRSF9"]:
        t = g.loc[gene, cols].astype(float)
        for mn, mg in [("effector", EFF), ("exhaustion", EXH)]:
            rr, pp = stats.spearmanr(t, zmodule(g, mg, cols))
            emit(pop_label, n, f"{gene}~{mn}_module_z", "Spearman_r", round(rr, 4), pp, role=role)

# tumour vs normal contrast — documents WHY the normals inflate the correlation
hr("1b. Why the 7 normals matter (they are DLL3-null and APM-high)")
for gene in ["DLL3"] + APM:
    t = g.loc[gene, TUMORS].astype(float)
    nn = g.loc[gene, NORMALS].astype(float)
    u, pu = stats.mannwhitneyu(t, nn)
    print(f"  {gene:<7} tumour median {t.median():7.3f} | normal median {nn.median():7.3f} | MWU P={pu:.4f}")
    emit("Jiang_tumor_vs_normal", 86, f"{gene}_tumor_vs_normal", "MannWhitneyU_P",
         f"{t.median():.3f} vs {nn.median():.3f}", pu, role="audit")

# ===========================================================================
# 2. Institutional discovery cohort
# ===========================================================================
hr("2. INSTITUTIONAL DISCOVERY COHORT  (primary = all 17 sequenced samples)")

d = pd.read_excel(os.path.join(BASE, "05_Source_Data/01_Institutional_RNAseq/"
                                     "InHouse_RNAseq_Normalized_Gene_Expression.xlsx"),
                  sheet_name="Normalized.expression")
SAMP = [c for c in d.columns if str(c).startswith("21029")]
d = d.set_index("Reference_GeneName")
d = d[~d.index.duplicated(keep="first")]
print(f"samples sequenced: {len(SAMP)}")
emit("Discovery", len(SAMP), "cohort_composition", "n_sequenced_n_analysed",
     f"{len(SAMP)}/{len(SAMP)}", extra="all 17 SCLC, all DV200>=40; no exclusions", role="audit")

dll3_log = np.log2(d.loc["DLL3", SAMP].astype(float))
ratio = np.log2(d.loc["CD274", SAMP].astype(float)) - np.log2(d.loc[APM, SAMP].astype(float)).mean(axis=0)

r, p = stats.spearmanr(dll3_log, ratio)
lo, hi = boot_spearman(dll3_log, ratio)
l0, l1 = loo_spearman(dll3_log, ratio)
print(f"  CD274/APM ratio ~ DLL3   r={r:+.4f} P={p:.4f}  boot95%[{lo:+.3f},{hi:+.3f}]  LOO[{l0:+.3f},{l1:+.3f}]")
emit("Discovery", len(SAMP), "CD274_APM_ratio~DLL3", "Spearman_r", round(r, 4), p,
     round(lo, 4), round(hi, 4), extra=f"LOO_min={l0:.4f};LOO_max={l1:.4f}", role="EXPLORATORY")

# cut-off sensitivity
med = dll3_log.median()
hi_g, lo_g = ratio[dll3_log >= med], ratio[dll3_log < med]
u, pmed = stats.mannwhitneyu(hi_g, lo_g)
print(f"  median split (DLL3 mRNA)  P={pmed:.4f}")
emit("Discovery", len(SAMP), "CD274_APM_ratio_median_split", "MannWhitneyU_P", round(u, 1), pmed,
     role="EXPLORATORY")

for lab, genes in [("APM_module_z", APM)]:
    rr, pp = stats.spearmanr(dll3_log, zmodule(d, genes, SAMP))
    print(f"  DLL3~{lab} (direct)  r={rr:+.4f} P={pp:.4f}")
    emit("Discovery", len(SAMP), f"DLL3~{lab}", "Spearman_r", round(rr, 4), pp, role="context")
rr, pp = stats.spearmanr(dll3_log, d.loc["CD274", SAMP].astype(float))
emit("Discovery", len(SAMP), "DLL3~CD274", "Spearman_r", round(rr, 4), pp, role="context")
print(f"  DLL3~CD274 (direct)  r={rr:+.4f} P={pp:.4f}")

t = np.log2(d.loc["TNFRSF9", SAMP].astype(float) + 1)
for mn, mg in [("effector", EFF), ("exhaustion", EXH)]:
    rr, pp = stats.spearmanr(t, zmodule(d, mg, SAMP))
    print(f"  TNFRSF9~{mn:<10} r={rr:+.4f} P={pp:.4g}")
    emit("Discovery", len(SAMP), f"TNFRSF9~{mn}_module_z", "Spearman_r", round(rr, 4), pp, role="SECONDARY")

# ridge-adjusted DLL3 coefficient
try:
    from sklearn.linear_model import Ridge
    covs = ["IFNG"] + NE4
    X = pd.concat([dll3_log.rename("DLL3"), np.log2(d.loc[covs, SAMP].astype(float) + 1).T], axis=1)
    Xs = (X - X.mean()) / X.std(ddof=1)
    ys = (ratio - ratio.mean()) / ratio.std(ddof=1)
    rg = Ridge(alpha=1.0).fit(Xs, ys)
    b = float(rg.coef_[0])
    print(f"  ridge (alpha=1) DLL3 coefficient on CD274/APM ratio = {b:+.4f}")
    emit("Discovery", len(SAMP), "CD274_APM_ratio~DLL3+IFNG+NE4_ridge", "ridge_beta_DLL3",
         round(b, 4), extra="alpha=1.0, predictors standardized", role="EXPLORATORY")
except ImportError:
    print("  (sklearn unavailable — ridge skipped)")

# ===========================================================================
# 3. Chan single-cell atlas
# ===========================================================================
hr("3. CHAN SCLC ATLAS  (19 donor_ids = 18 named donors + 1 pleural-effusion sample)")

import anndata as ad

epi = ad.read_h5ad(os.path.join(BASE, "05_Source_Data/04_Public_scRNA_Chan_Atlas/"
                                      "Chan_HTAN_SCLC_Epithelial_Cells_54313cells.h5ad"), backed="r")
epi_ids = {}
for s in ["DLL3"] + APM:
    j = epi.var.index[epi.var.feature_name.astype(str) == s]
    if len(j):
        epi_ids[s] = j[0]
sub = epi[:, list(epi_ids.values())].to_memory()
Xe = np.asarray(sub.X.todense()) if hasattr(sub.X, "todense") else np.asarray(sub.X)
epi_df = pd.DataFrame(Xe, columns=list(epi_ids.keys()))
epi_df["donor"] = epi.obs.donor_id.values
epi_mean = epi_df.groupby("donor", observed=True).mean()
dll3_donor = epi_mean["DLL3"]
dll3_pct = (pd.DataFrame({"donor": epi.obs.donor_id.values, "v": Xe[:, 0] > 0})
            .groupby("donor", observed=True).v.mean() * 100)
print(f"  per-donor mean epithelial DLL3 computed for {len(dll3_donor)} donors")

# tumour-intrinsic APM at the epithelial single-cell level (different estimand from bulk)
apm_epi = epi_mean[[c for c in APM if c in epi_mean.columns]]
apm_epi_z = apm_epi.apply(lambda c: (c - c.mean()) / c.std(ddof=1)).mean(axis=1)
r_e, p_e = stats.spearmanr(dll3_donor, apm_epi_z)
print(f"  Chan EPITHELIAL DLL3~APM (z-module, per donor) r={r_e:+.4f} P={p_e:.4f}"
      f"   [v1 manuscript claimed +0.254 / P=0.29 — not reproduced]")
emit("Chan_epithelial", len(apm_epi_z), "DLL3~APM_module_z", "Spearman_r", round(r_e, 4), p_e,
     extra="tumour-cell estimand; contrast with bulk tissue estimand", role="context")
del epi, sub, Xe, epi_df

imm = ad.read_h5ad(os.path.join(BASE, "05_Source_Data/04_Public_scRNA_Chan_Atlas/"
                                      "Chan_HTAN_SCLC_Immune_Cells_16475cells.h5ad"))
obs = imm.obs
CD8 = ["CD8+ Texh", "CD8+ Teff", "CD8+ Tmem"]
TCELL = [c for c in obs.author_cell_type.cat.categories
         if str(c) in CD8 + ["CD4+ Tconv", "CD4+ Treg", "Tgd"]]


def gene_vec(sym):
    j = imm.var.index[imm.var.feature_name.astype(str) == sym]
    if len(j) == 0:
        return None
    col = imm[:, j[0]].X
    return (np.asarray(col.todense()).ravel() if hasattr(col, "todense") else np.asarray(col).ravel())


mk = {s: gene_vec(s) for s in ["TNFRSF9", "PDCD1", "HAVCR2", "GZMB"]}
cell = pd.DataFrame({"donor": obs.donor_id.values, "ct": obs.author_cell_type.astype(str).values,
                     **{k: (v > 0) for k, v in mk.items() if v is not None}})

# --- composition correlations
tot_T = cell[cell.ct.isin([str(x) for x in TCELL])].groupby("donor", observed=True).size()
texh = cell[cell.ct == "CD8+ Texh"].groupby("donor", observed=True).size().reindex(tot_T.index, fill_value=0)
pct_texh = (texh / tot_T * 100)
cd8 = cell[cell.ct.isin(CD8)]
cd8n = cd8.groupby("donor", observed=True).size()


def pct_of_cd8(mask):
    s = cd8[mask].groupby("donor", observed=True).size().reindex(cd8n.index, fill_value=0)
    return s / cd8n * 100


series = {
    "CD8Texh_pct_of_T": pct_texh,
    "TNFRSF9pos_pct_of_CD8": pct_of_cd8(cd8.TNFRSF9),
    "TNFRSF9+PD1+_pct_of_CD8": pct_of_cd8(cd8.TNFRSF9 & cd8.PDCD1),
    "TNFRSF9+TIM3+_pct_of_CD8": pct_of_cd8(cd8.TNFRSF9 & cd8.HAVCR2),
}
print()
for name, s in series.items():
    common = [x for x in s.index if x in dll3_donor.index]
    r, p = stats.spearmanr(dll3_donor[common], s[common])
    flo, fhi = fisher_ci(r, len(common))
    l0, l1 = loo_spearman(dll3_donor[common].values, s[common].values)
    role = "PRIMARY" if name == "CD8Texh_pct_of_T" else "KEY_NEGATIVE"
    print(f"  DLL3~{name:<26} r={r:+.4f} P={p:.4f}  95%CI[{flo:+.3f},{fhi:+.3f}]  LOO[{l0:+.3f},{l1:+.3f}]")
    emit("Chan_scRNA", len(common), f"DLL3~{name}", "Spearman_r", round(r, 4), p,
         round(flo, 4), round(fhi, 4), extra=f"LOO_min={l0:.4f};LOO_max={l1:.4f}", role=role)
    # metric-dependence check: %DLL3+ cells instead of mean expression
    r2_, p2_ = stats.spearmanr(dll3_pct[common], s[common])
    emit("Chan_scRNA", len(common), f"DLL3pctpos~{name}", "Spearman_r", round(r2_, 4), p2_,
         extra="sensitivity: DLL3 as % positive cells", role="sensitivity")
    print(f"      [sensitivity: DLL3 as %positive cells] r={r2_:+.4f} P={p2_:.4f}")

# --- TNFRSF9 positivity by subtype
print()
for ct in CD8 + ["CD4+ Treg"]:
    sl = cell[cell.ct == ct]
    if len(sl):
        pctp = sl.TNFRSF9.mean() * 100
        print(f"  TNFRSF9+ rate in {ct:<12} = {pctp:5.1f}%  (n={len(sl)} cells)")
        emit("Chan_scRNA", len(sl), f"TNFRSF9pos_rate_{ct}", "percent", round(pctp, 1), role="descriptive")

print(f"\n  CD8 compartment sizes: " + ", ".join(
    f"{c}={int((cell.ct == c).sum())}" for c in CD8))
tx = cell[cell.ct == "CD8+ Texh"].groupby("donor", observed=True).size().sort_values(ascending=False)
print(f"  Texh concentration: top donor {tx.index[0]} contributes {tx.iloc[0]}/{tx.sum()} "
      f"({100*tx.iloc[0]/tx.sum():.0f}%) of all Texh cells")
emit("Chan_scRNA", int(tx.sum()), "Texh_top_donor_share", "percent",
     round(100 * tx.iloc[0] / tx.sum(), 1), extra=f"donor={tx.index[0]}", role="audit")

# --- paired within-donor tests
def paired_test(pool_mask, min_sub, min_cells, label, role):
    sl = cell[pool_mask]
    out = []
    for donor, gdf in sl.groupby("donor", observed=True):
        pos, neg = gdf[gdf.TNFRSF9], gdf[~gdf.TNFRSF9]
        if len(gdf) < min_cells or len(pos) < min_sub or len(neg) < min_sub:
            continue
        row = {"donor": donor}
        for m, col in [("PD-1", "PDCD1"), ("TIM-3", "HAVCR2"), ("GZMB", "GZMB")]:
            row[m] = pos[col].mean() * 100 - neg[col].mean() * 100
        out.append(row)
    res = pd.DataFrame(out)
    for m in ["PD-1", "TIM-3", "GZMB"]:
        v = res[m].dropna().values
        if len(v) < 3:
            continue
        nz = v[v != 0]
        p = stats.wilcoxon(nz)[1] if len(nz) else 1.0
        higher = int((v > 0).sum())
        # matched-pairs rank-biserial
        rb = np.nan
        if len(nz):
            ranks = stats.rankdata(np.abs(nz))
            rb = (ranks[nz > 0].sum() - ranks[nz < 0].sum()) / ranks.sum()
        print(f"  [{label}] {m:<6} n={len(v):2d} higher={higher:2d}/{len(v):<2d} "
              f"median_delta={np.median(v):+6.2f}pp  rb_r={rb:+.3f}  P={p:.4g}")
        emit("Chan_scRNA", len(v), f"{label}_{m}", "paired_Wilcoxon_median_delta_pp",
             round(float(np.median(v)), 2), p, extra=f"higher={higher}/{len(v)};rank_biserial={rb:.3f}",
             role=role)
    return res


print()
cd8_mask = cell.ct.isin(CD8)
res_primary = paired_test(cd8_mask, 3, 10, "pooledCD8_PRIMARY_min3", "PRIMARY")
res_sens = paired_test(cd8_mask, 1, 1, "pooledCD8_sensitivity_min1", "sensitivity")
print()
strat = {}
for ct in CD8:
    strat[ct] = paired_test(cell.ct == ct, 3, 3, f"{ct.replace('+','').replace(' ','')}_only", "stratified")

# --- the decisive interaction test (Tmem vs Texh within the same donor)
print()
a_, b_ = strat["CD8+ Tmem"], strat["CD8+ Texh"]
if len(a_) and len(b_):
    j = a_.merge(b_, on="donor", suffixes=("_tmem", "_texh"))
    print(f"  INTERACTION TEST — donors evaluable in BOTH Tmem and Texh: n={len(j)}")
    for m in ["PD-1", "TIM-3", "GZMB"]:
        dvals = (j[f"{m}_tmem"] - j[f"{m}_texh"]).dropna().values
        nz = dvals[dvals != 0]
        p = stats.wilcoxon(nz)[1] if len(nz) else 1.0
        print(f"    {m:<6} median(Tmem-Texh)={np.median(dvals):+6.2f}pp  "
              f"{int((dvals>0).sum())}/{len(dvals)}  P={p:.4g}   "
              f"{'-> NO localisation difference' if p >= .05 else ''}")
        emit("Chan_scRNA", len(dvals), f"INTERACTION_Tmem_minus_Texh_{m}",
             "paired_Wilcoxon_median_delta_pp", round(float(np.median(dvals)), 2), p,
             extra=f"higher_in_Tmem={int((dvals>0).sum())}/{len(dvals)}", role="KEY_INTERACTION")

# lineage origin of TNFRSF9+ T cells
tn = cell[cell.TNFRSF9 & cell.ct.isin([str(x) for x in TCELL])]
lin = tn.assign(lin=np.where(tn.ct.isin(CD8), "CD8", np.where(tn.ct == "CD4+ Treg", "Treg", "Other")))
per = lin.groupby(["donor", "lin"], observed=True).size().unstack(fill_value=0)
per = per.div(per.sum(axis=1), axis=0) * 100
if {"CD8", "Treg"} <= set(per.columns):
    p = stats.wilcoxon(per["CD8"], per["Treg"])[1]
    print(f"\n  TNFRSF9+ lineage: CD8 median {per['CD8'].median():.0f}% vs Treg {per['Treg'].median():.0f}%  P={p:.4f}")
    emit("Chan_scRNA", len(per), "TNFRSF9pos_lineage_CD8_vs_Treg", "Wilcoxon_P",
         f"CD8 {per['CD8'].median():.0f}% vs Treg {per['Treg'].median():.0f}%", p, role="descriptive")

# ===========================================================================
# 4. mIF
# ===========================================================================
hr("4. INSTITUTIONAL mIF  (10 patients, 193 ROIs; FDR family = all 51 phenotypes)")

mif = pd.read_csv(os.path.join(BASE, "05_Source_Data/06_Derived_Analysis_Tables/"
                                     "Derived_mIF_Patient_AreaWeighted_Means.csv"))
print(f"  patient-level table: {mif.shape[0]} patients x {mif.shape[1]} columns")
idc = [c for c in mif.columns if mif[c].dtype == object or "atient" in c or "ROI" in c or "score" in c.lower()]
print(f"  non-marker columns: {idc[:8]}")
emit("mIF", mif.shape[0], "cohort_composition", "n_patients",
     mif.shape[0], extra="193 ROIs; selected from 31 patients/632 ROIs", role="audit")

res51 = pd.read_csv(os.path.join(BASE, "05_Source_Data/06_Derived_Analysis_Tables/"
                                       "Derived_mIF_ROI_All_51Markers_Results.csv"))
print(f"  51-marker result table: {res51.shape}; columns {list(res51.columns)}")
if "q" in res51.columns:
    print(f"  markers reaching q<0.05: {int((res51.q < 0.05).sum())} of {len(res51)}")
    emit("mIF", 10, "n_markers_FDR_significant", "count", int((res51.q < 0.05).sum()),
         extra=f"family size = {len(res51)} phenotypes, BH", role="PRIMARY")
    for _, row in res51.sort_values("p").head(6).iterrows():
        emit("mIF", 10, f"DLL3~{row['marker']}", "Spearman_r", round(row["r"], 4), row["p"],
             extra=f"q={row['q']:.4f};FC={row.get('fc', np.nan):.2f}", role="supportive")
        print(f"    {row['marker']:<22} r={row['r']:+.3f} P={row['p']:.3f} q={row['q']:.3f}")

# ===========================================================================
# WRITE
# ===========================================================================
out = pd.DataFrame(ROWS)
DD = os.path.join(BASE, "05_Source_Data/06_Derived_Analysis_Tables")
dest = os.path.join(DD, "Derived_MASTER_Statistics_v2.csv")
out.to_csv(dest, index=False)

# --- plotting caches, so figure scripts never re-read the 8.9 GB h5ad --------
plot_chan = pd.DataFrame({"DLL3_mean": dll3_donor, "DLL3_pct_pos": dll3_pct})
for nm, s in series.items():
    plot_chan[nm] = s
plot_chan.to_csv(os.path.join(DD, "Plotdata_Chan_PerDonor_v2.csv"))
res_primary.to_csv(os.path.join(DD, "Plotdata_Chan_Paired_PRIMARY_min3_v2.csv"), index=False)

pd.DataFrame({
    "sample": list(g.columns),
    "is_tumor": [c not in NORMALS for c in g.columns],
    "DLL3_log2": g.loc["DLL3"].astype(float).values,
    "APM_z_tumoronly": zmodule(g, APM, TUMORS).reindex(g.columns).values,
    "APM_raw_log2_mean": g.loc[APM].astype(float).mean(axis=0).values,
}).to_csv(os.path.join(DD, "Plotdata_Jiang_v2.csv"), index=False)

pd.DataFrame({
    "sample": SAMP,
    "DLL3_log2": dll3_log.values,
    "CD274_APM_ratio": ratio.values,
}).to_csv(os.path.join(DD, "Plotdata_Discovery_v2.csv"), index=False)
print("  plotting caches written (Plotdata_*_v2.csv)")
hr("DONE")
print(f"  {len(out)} rows written to\n  {dest}")
print(f"  seed={SEED}, bootstrap iterations={N_BOOT}")
print(f"  python={sys.version.split()[0]}, pandas={pd.__version__}, "
      f"numpy={np.__version__}, scipy={stats.__name__.split('.')[0]}")
