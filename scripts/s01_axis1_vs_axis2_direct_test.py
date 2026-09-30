"""v3 P0-3: Direct test of Axis 1 (epithelial APM) vs Axis 2 (TNFRSF9-defined CD8
abundance) within the Chan atlas, plus formal equivalence (TOST) bounds for the
central DLL3 -> Axis 2 null.

Rationale: the v2 Figure 4 legend stated "The two axes were tested for association
and none was detected (r = +0.042 to +0.105)". Those r values are DLL3-versus-Axis-2
correlations; the association BETWEEN Axis 1 and Axis 2 had never been tested.
This script performs the missing test, so that the claim is either made true or
withdrawn on the basis of evidence.
"""
import os
import numpy as np
import pandas as pd
import anndata as ad
from scipy import stats

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(BASE, "14_Revision_v3_20260806/outputs")
os.makedirs(OUT, exist_ok=True)
APM5 = ["HLA-A", "HLA-B", "B2M", "TAP1", "TAP2"]
APM6 = ["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2"]
AX2 = ["TNFRSF9pos_pct_of_CD8", "TNFRSF9+PD1+_pct_of_CD8", "TNFRSF9+TIM3+_pct_of_CD8"]
RNG = np.random.default_rng(20260806)
N_BOOT = 2000


def fisher_ci(r, n, a=0.05):
    z, se = np.arctanh(r), 1 / np.sqrt(n - 3)
    k = stats.norm.ppf(1 - a / 2)
    return float(np.tanh(z - k * se)), float(np.tanh(z + k * se))


def tost_min_delta(r, n, a=0.05):
    """Smallest equivalence margin delta for which two one-sided tests both reject."""
    z, se = np.arctanh(r), 1 / np.sqrt(n - 3)
    k = stats.norm.ppf(1 - a)
    return float(np.tanh(max(abs(z + k * se), abs(z - k * se))))


def tost_p(r, n, d):
    z, se = np.arctanh(r), 1 / np.sqrt(n - 3)
    zd = np.arctanh(d)
    return float(max(stats.norm.cdf((z - zd) / se), stats.norm.sf((z + zd) / se)))


def bootstrap_ci(x, y, n_boot=N_BOOT):
    """Percentile bootstrap interval for a donor-level Spearman correlation."""
    x, y = np.asarray(x, float), np.asarray(y, float)
    values = []
    for idx in RNG.integers(0, len(x), size=(n_boot, len(x))):
        r = stats.spearmanr(x[idx], y[idx]).statistic
        if np.isfinite(r):
            values.append(r)
    return tuple(float(v) for v in np.percentile(values, [2.5, 97.5]))


def loo_spearman(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    values = [stats.spearmanr(np.delete(x, i), np.delete(y, i)).statistic
              for i in range(len(x))]
    return float(np.min(values)), float(np.max(values))


def bh_adjust(p_values):
    """Benjamini-Hochberg adjusted P values, returned in the original order."""
    p = np.asarray(p_values, float)
    order = np.argsort(p)
    ranked = p[order]
    adjusted = ranked * len(p) / np.arange(1, len(p) + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    out = np.empty_like(adjusted)
    out[order] = np.minimum(adjusted, 1.0)
    return out


# --- per-donor epithelial APM from the Chan atlas -------------------------------
epi = ad.read_h5ad(os.path.join(BASE, "05_Source_Data/04_Public_scRNA_Chan_Atlas/"
                                      "Chan_HTAN_SCLC_Epithelial_Cells_54313cells.h5ad"),
                   backed="r")
feat = epi.var.feature_name.astype(str)
ids = {s: epi.var.index[feat == s][0] for s in ["DLL3"] + APM6 if (feat == s).any()}
sub = epi[:, list(ids.values())].to_memory()
X = np.asarray(sub.X.todense()) if hasattr(sub.X, "todense") else np.asarray(sub.X)
df = pd.DataFrame(X, columns=list(ids.keys()))
df["donor"] = epi.obs.donor_id.values
m = df.groupby("donor", observed=True).mean()
apm5_cols = [c for c in APM5 if c in m.columns]
apm6_cols = [c for c in APM6 if c in m.columns]
apm5_z = m[apm5_cols].apply(lambda c: (c - c.mean()) / c.std(ddof=1)).mean(axis=1)
apm6_z = m[apm6_cols].apply(lambda c: (c - c.mean()) / c.std(ddof=1)).mean(axis=1)

per = pd.read_csv(os.path.join(BASE, "05_Source_Data/06_Derived_Analysis_Tables/"
                                     "Plotdata_Chan_PerDonor_v2.csv"), index_col=0)
per["epithelial_APM_z"] = apm5_z
per["epithelial_APM_z_5gene_primary"] = apm5_z
per["epithelial_APM_z_6gene_sensitivity"] = apm6_z
per["n_epithelial_cells"] = df.groupby("donor", observed=True).size()

# CD8 evaluability sensitivity is derived independently from the source immune object.
imm = ad.read_h5ad(os.path.join(BASE, "05_Source_Data/04_Public_scRNA_Chan_Atlas/"
                                      "Chan_HTAN_SCLC_Immune_Cells_16475cells.h5ad"),
                   backed="r")
cd8_types = {"CD8+ Texh", "CD8+ Teff", "CD8+ Tmem"}
imm_obs = imm.obs[["donor_id", "author_cell_type"]].copy()
imm_obs["author_cell_type"] = imm_obs.author_cell_type.astype(str)
n_cd8 = imm_obs[imm_obs.author_cell_type.isin(cd8_types)].groupby(
    "donor_id", observed=True).size()
per["n_CD8"] = n_cd8
per.to_csv(os.path.join(OUT, "Chan_PerDonor_with_epithelialAPM_v3.csv"))
print(f"Primary epithelial APM genes: {apm5_cols}")
print(f"Six-gene sensitivity APM genes: {apm6_cols}")
print(f"per-donor table: {per.shape[0]} donors")

rows = []
direct_row_indices = []
direct_p_values = []

print("\n=== DIRECT TEST  Axis 1 (epithelial APM) vs Axis 2 (TNFRSF9-defined CD8 abundance) ===")
for a2 in AX2:
    d = per[["epithelial_APM_z", a2]].dropna()
    n = len(d)
    r, p = stats.spearmanr(d.iloc[:, 0], d.iloc[:, 1])
    lo, hi = fisher_ci(r, n)
    blo, bhi = bootstrap_ci(d.iloc[:, 0], d.iloc[:, 1])
    llo, lhi = loo_spearman(d.iloc[:, 0], d.iloc[:, 1])
    threshold_results = {}
    for threshold in (20, 50):
        ds = per.loc[per.n_CD8 >= threshold, ["epithelial_APM_z", a2]].dropna()
        sr, sp = stats.spearmanr(ds.iloc[:, 0], ds.iloc[:, 1])
        threshold_results[f"n_cd8_ge{threshold}"] = len(ds)
        threshold_results[f"r_cd8_ge{threshold}"] = round(sr, 4)
        threshold_results[f"p_cd8_ge{threshold}"] = round(sp, 6)
    print(f"  APM_epi ~ {a2:<26} n={n:2d}  r={r:+.4f}  P={p:.4f}  "
          f"Fisher 95% CI [{lo:+.3f}, {hi:+.3f}]  "
          f"bootstrap 95% CI [{blo:+.3f}, {bhi:+.3f}]  LOO [{llo:+.3f}, {lhi:+.3f}]")
    direct_row_indices.append(len(rows))
    direct_p_values.append(p)
    rows.append(dict(analysis="Axis1_vs_Axis2_direct", x="epithelial_APM_z", y=a2, n=n,
                     spearman_r=round(r, 4), ci_low=round(lo, 4), ci_high=round(hi, 4),
                     bootstrap_ci_low=round(blo, 4), bootstrap_ci_high=round(bhi, 4),
                     loo_min_r=round(llo, 4), loo_max_r=round(lhi, 4),
                     p=round(p, 6), bh_q=None,
                     endpoint_role=("principal exploratory endpoint" if a2 == AX2[0]
                                    else "supportive nested endpoint"),
                     tost_p_delta_0p5=None, tost_p_delta_0p3=None,
                     smallest_equivalence_margin=None, **threshold_results))

direct_q_values = bh_adjust(direct_p_values)
for row_index, q in zip(direct_row_indices, direct_q_values):
    rows[row_index]["bh_q"] = round(float(q), 6)
print("  BH-adjusted q values:", ", ".join(f"{q:.4f}" for q in direct_q_values))

print("\n=== SENSITIVITY  six-gene epithelial APM including HLA-C ===")
for a2 in AX2:
    d = per[["epithelial_APM_z_6gene_sensitivity", a2]].dropna()
    r, p = stats.spearmanr(d.iloc[:, 0], d.iloc[:, 1])
    lo, hi = fisher_ci(r, len(d))
    print(f"  APM6_epi ~ {a2:<26} n={len(d):2d}  r={r:+.4f}  P={p:.4f}")
    rows.append(dict(analysis="Axis1_vs_Axis2_6gene_sensitivity",
                     x="epithelial_APM_z_6gene_sensitivity", y=a2, n=len(d),
                     spearman_r=round(r, 4), ci_low=round(lo, 4), ci_high=round(hi, 4),
                     bootstrap_ci_low=None, bootstrap_ci_high=None,
                     loo_min_r=None, loo_max_r=None, p=round(p, 6), bh_q=None,
                     endpoint_role="APM-definition sensitivity",
                     tost_p_delta_0p5=None, tost_p_delta_0p3=None,
                     smallest_equivalence_margin=None,
                     n_cd8_ge20=None, r_cd8_ge20=None, p_cd8_ge20=None,
                     n_cd8_ge50=None, r_cd8_ge50=None, p_cd8_ge50=None))

print("\n=== CONTEXT  DLL3 vs epithelial APM in the Chan atlas ===")
for x, role in [("epithelial_APM_z_5gene_primary", "primary 5-gene APM context"),
                ("epithelial_APM_z_6gene_sensitivity", "six-gene APM sensitivity context")]:
    d = per[["DLL3_mean", x]].dropna()
    r, p = stats.spearmanr(d.iloc[:, 0], d.iloc[:, 1])
    lo, hi = fisher_ci(r, len(d))
    print(f"  DLL3 ~ {x:<38} n={len(d):2d}  r={r:+.4f}  P={p:.4f}")
    rows.append(dict(analysis="DLL3_to_epithelial_APM_context", x="DLL3_mean", y=x,
                     n=len(d), spearman_r=round(r, 4), ci_low=round(lo, 4),
                     ci_high=round(hi, 4), bootstrap_ci_low=None,
                     bootstrap_ci_high=None, loo_min_r=None, loo_max_r=None,
                     p=round(p, 6), bh_q=None, endpoint_role=role,
                     tost_p_delta_0p5=None, tost_p_delta_0p3=None,
                     smallest_equivalence_margin=None,
                     n_cd8_ge20=None, r_cd8_ge20=None, p_cd8_ge20=None,
                     n_cd8_ge50=None, r_cd8_ge50=None, p_cd8_ge50=None))

print("\n=== EQUIVALENCE BOUNDS  central null: DLL3 -> Axis 2 abundance ===")
for a2 in AX2:
    d = per[["DLL3_mean", a2]].dropna()
    n = len(d)
    r, p = stats.spearmanr(d.iloc[:, 0], d.iloc[:, 1])
    lo, hi = fisher_ci(r, n)
    md = tost_min_delta(r, n)
    print(f"  DLL3 ~ {a2:<26} n={n:2d}  r={r:+.4f}  P={p:.4f}  95% CI [{lo:+.3f}, {hi:+.3f}]  "
          f"TOST P(d=0.5)={tost_p(r, n, 0.5):.4f}  TOST P(d=0.3)={tost_p(r, n, 0.3):.4f}  "
          f"smallest margin={md:.3f}")
    rows.append(dict(analysis="DLL3_to_Axis2_equivalence", x="DLL3_mean", y=a2, n=n,
                     spearman_r=round(r, 4), ci_low=round(lo, 4), ci_high=round(hi, 4),
                     bootstrap_ci_low=None, bootstrap_ci_high=None,
                     loo_min_r=None, loo_max_r=None,
                     p=round(p, 6), bh_q=None, endpoint_role="bounded null analysis",
                     tost_p_delta_0p5=round(tost_p(r, n, 0.5), 4),
                     tost_p_delta_0p3=round(tost_p(r, n, 0.3), 4),
                     smallest_equivalence_margin=round(md, 3),
                     n_cd8_ge20=None, r_cd8_ge20=None, p_cd8_ge20=None,
                     n_cd8_ge50=None, r_cd8_ge50=None, p_cd8_ge50=None))

res = pd.DataFrame(rows)
res.to_csv(os.path.join(OUT, "SuppTable_Axis1_vs_Axis2_and_Equivalence_v3.csv"), index=False)
print(f"\nwritten -> {OUT}/SuppTable_Axis1_vs_Axis2_and_Equivalence_v3.csv")
