"""v3 P1-1: Detection-depth sensitivity for the TNFRSF9-defined CD8 phenotype.

The peer review raised a specific technical objection: TNFRSF9 positivity is
defined as detected transcript (normalized expression > 0), and PD-1/TIM-3
positivity is defined the same way in the same cells. Cells with more total UMIs
and more detected genes will co-detect several low-abundance transcripts more
often, so part of the very large rank-biserial effect could be global detection
bias rather than a biological co-expression programme.

This script quantifies that. It reports:
  1. total UMI, detected genes and mitochondrial fraction in TNFRSF9+ vs
     TNFRSF9- CD8 cells (the exposure check).
  2. A cell-level logistic model of marker positivity on TNFRSF9 positivity,
     adjusted for log10 total counts, CD8 subtype and donor, with donor
     cluster-robust standard errors.
  3. A continuous-expression sensitivity analysis using the normalized layer
     rather than binary detection, as a paired within-donor test.
  4. A depth-matched sensitivity analysis in which TNFRSF9- cells are matched to
     TNFRSF9+ cells on total-count decile within donor.
  5. Per-donor cell counts for every quantity used in the paired analyses.
"""
import os
import numpy as np
import pandas as pd
import anndata as ad
import scipy.sparse as sp
import statsmodels.api as sm
import statsmodels.formula.api as smf
from scipy import stats

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(BASE, "14_Revision_v3_20260806/outputs")
os.makedirs(OUT, exist_ok=True)
CD8 = ["CD8+ Texh", "CD8+ Teff", "CD8+ Tmem"]
MARKERS = {"PD-1": "PDCD1", "TIM-3": "HAVCR2", "GZMB": "GZMB"}
RNG = np.random.default_rng(20260806)

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
    "mito_frac": imm.obs.mito_frac.values.astype(float),
    "assay": imm.obs.assay.astype(str).values,
})
for sym in ["TNFRSF9"] + list(MARKERS.values()):
    cell[sym + "_pos"] = vec(sym) > 0
    cell[sym + "_expr"] = vec(sym, "normalized")

c8 = cell[cell.ct.isin(CD8)].copy()
c8["log10_counts"] = np.log10(c8.total_counts)
print(f"CD8 pool: {len(c8)} cells, {c8.donor.nunique()} donors")

rows = []

# --- 1. exposure check -----------------------------------------------------------
print("\n=== 1. Sequencing depth in TNFRSF9+ vs TNFRSF9- CD8 cells ===")
for q, lab in [("total_counts", "total UMI"), ("n_genes", "detected genes"),
               ("mito_frac", "mitochondrial fraction")]:
    a = c8.loc[c8.TNFRSF9_pos, q]
    b = c8.loc[~c8.TNFRSF9_pos, q]
    u, p = stats.mannwhitneyu(a, b)
    cles = u / (len(a) * len(b))          # common-language effect size / AUC
    print(f"  {lab:<24} TNFRSF9+ median {a.median():10.4f} | TNFRSF9- median {b.median():10.4f}"
          f" | AUC={cles:.3f} | P={p:.3g}")
    rows.append(dict(analysis="depth_exposure_check", quantity=lab,
                     tnfrsf9_pos_median=round(float(a.median()), 4),
                     tnfrsf9_neg_median=round(float(b.median()), 4),
                     auc=round(float(cles), 4), p=p, n_pos=len(a), n_neg=len(b)))

# --- 2. cell-level logistic model, depth- and donor-adjusted ---------------------
print("\n=== 2. Cell-level logistic model (donor cluster-robust SE) ===")
print("    marker_pos ~ TNFRSF9_pos + log10_counts + CD8 subtype + donor")
for lab, sym in MARKERS.items():
    d = c8[[f"{sym}_pos", "TNFRSF9_pos", "log10_counts", "ct", "donor"]].copy()
    d.columns = ["y", "tn", "depth", "ct", "donor"]
    d["y"] = d.y.astype(int)
    d["tn"] = d.tn.astype(int)
    for formula, tag in [("y ~ tn", "unadjusted"),
                         ("y ~ tn + depth", "+ depth"),
                         ("y ~ tn + depth + C(ct) + C(donor)", "+ depth + subtype + donor")]:
        m = smf.glm(formula, data=d, family=sm.families.Binomial()).fit(
            cov_type="cluster", cov_kwds={"groups": d.donor})
        b = float(m.params["tn"])
        lo, hi = [float(x) for x in m.conf_int().loc["tn"]]
        p = float(m.pvalues["tn"])
        print(f"  {lab:<6} {tag:<28} OR={np.exp(b):6.2f}  95% CI [{np.exp(lo):.2f}, {np.exp(hi):.2f}]"
              f"  P={p:.3g}")
        rows.append(dict(analysis="cell_level_logistic", quantity=lab, model=tag,
                         odds_ratio=round(float(np.exp(b)), 3),
                         ci_low=round(float(np.exp(lo)), 3), ci_high=round(float(np.exp(hi)), 3),
                         p=p, n_cells=len(d)))

# --- 3. continuous expression, paired within donor -------------------------------
print("\n=== 3. Continuous normalized expression (not binary detection), paired by donor ===")
for lab, sym in MARKERS.items():
    recs = []
    for donor, gdf in c8.groupby("donor", observed=True):
        pos, neg = gdf[gdf.TNFRSF9_pos], gdf[~gdf.TNFRSF9_pos]
        if len(gdf) < 10 or len(pos) < 3 or len(neg) < 3:
            continue
        recs.append(float(pos[f"{sym}_expr"].mean() - neg[f"{sym}_expr"].mean()))
    v = np.array(recs)
    nz = v[v != 0]
    p = float(stats.wilcoxon(nz)[1]) if len(nz) else 1.0
    ranks = stats.rankdata(np.abs(nz))
    rb = float((ranks[nz > 0].sum() - ranks[nz < 0].sum()) / ranks.sum()) if len(nz) else np.nan
    print(f"  {lab:<6} n={len(v):2d}  higher in TNFRSF9+ {int((v>0).sum())}/{len(v)}  "
          f"median delta={np.median(v):+.4f}  rank-biserial={rb:+.3f}  P={p:.4g}")
    rows.append(dict(analysis="continuous_expression_paired", quantity=lab, n_donors=len(v),
                     higher_in_pos=int((v > 0).sum()), median_delta=round(float(np.median(v)), 4),
                     rank_biserial=round(rb, 3), p=p))

# --- 4. depth-matched sensitivity -------------------------------------------------
print("\n=== 4. Depth-matched sensitivity (TNFRSF9- matched on total-count decile within donor) ===")
for lab, sym in MARKERS.items():
    recs = []
    for donor, gdf in c8.groupby("donor", observed=True):
        gdf = gdf.copy()
        try:
            gdf["dec"] = pd.qcut(gdf.total_counts, 10, labels=False, duplicates="drop")
        except ValueError:
            continue
        pos = gdf[gdf.TNFRSF9_pos]
        neg = gdf[~gdf.TNFRSF9_pos]
        if len(pos) < 3 or len(neg) < 3:
            continue
        keep = []
        for dec, npos in pos.dec.value_counts().items():
            pool = neg[neg.dec == dec]
            if len(pool) == 0:
                continue
            k = min(len(pool), int(npos) * 3)
            keep.append(pool.sample(k, random_state=int(RNG.integers(1e6))))
        if not keep:
            continue
        negm = pd.concat(keep)
        recs.append(float(pos[f"{sym}_pos"].mean() * 100 - negm[f"{sym}_pos"].mean() * 100))
    v = np.array(recs)
    nz = v[v != 0]
    p = float(stats.wilcoxon(nz)[1]) if len(nz) else 1.0
    print(f"  {lab:<6} n={len(v):2d}  higher in TNFRSF9+ {int((v>0).sum())}/{len(v)}  "
          f"median delta={np.median(v):+.2f} pp  P={p:.4g}")
    rows.append(dict(analysis="depth_matched_paired", quantity=lab, n_donors=len(v),
                     higher_in_pos=int((v > 0).sum()), median_delta=round(float(np.median(v)), 2), p=p))

# --- 5. per-donor cell counts -----------------------------------------------------
cnt = c8.groupby("donor", observed=True).agg(
    n_CD8=("ct", "size"),
    n_TNFRSF9_pos=("TNFRSF9_pos", "sum"),
    n_TNFRSF9_neg=("TNFRSF9_pos", lambda s: int((~s).sum())),
    n_TNFRSF9_PD1_pos=("PDCD1_pos", "sum"),
    median_total_counts=("total_counts", "median"),
    median_genes=("n_genes", "median"))
for ct in CD8:
    cnt[f"n_{ct.replace('+','').replace(' ','_')}"] = (
        c8[c8.ct == ct].groupby("donor", observed=True).size().reindex(cnt.index, fill_value=0))
cnt["evaluable_primary_min3"] = (cnt.n_CD8 >= 10) & (cnt.n_TNFRSF9_pos >= 3) & (cnt.n_TNFRSF9_neg >= 3)
epi_n = pd.read_csv(os.path.join(OUT, "Chan_PerDonor_with_epithelialAPM_v3.csv"), index_col=0)
if "n_epithelial_cells" in epi_n.columns:
    cnt["n_epithelial"] = epi_n["n_epithelial_cells"].reindex(cnt.index)
cnt.to_csv(os.path.join(OUT, "SuppTable_Chan_PerDonor_CellCounts_v3.csv"))
print(f"\n=== 5. Per-donor cell counts ({len(cnt)} donors; "
      f"{int(cnt.evaluable_primary_min3.sum())} evaluable under the primary rule) ===")
print(cnt.to_string())

pd.DataFrame(rows).to_csv(os.path.join(OUT, "SuppTable_DetectionDepth_Sensitivity_v3.csv"), index=False)
print(f"\nwritten -> {OUT}/SuppTable_DetectionDepth_Sensitivity_v3.csv")
print(f"written -> {OUT}/SuppTable_Chan_PerDonor_CellCounts_v3.csv")
