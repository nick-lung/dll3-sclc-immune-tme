"""v3 P0-4: External replication (Jiang + George/Cologne) and leukocyte-content
sensitivity, with the harmonisation the peer review required.

Corrections applied relative to 12_External_Replication_20260806/s14_paper1_addendum.py:
  1. APM module is computed from the SAME five genes in BOTH cohorts
     (HLA-A, HLA-B, B2M, TAP1, TAP2). HLA-C has zero variance across all 81
     Cologne samples, so the earlier version compared a 6-gene Jiang module with
     a 5-gene Cologne module. The 6-gene Jiang module is retained as a
     sensitivity analysis.
  2. Confidence intervals in the forest plot use ONE method (Fisher z) for both
     cohorts and the pooled estimate. Percentile bootstrap CIs are still
     reported in the table as a robustness check, but methods are no longer
     mixed within a single figure.
  3. The analysis is labelled POST HOC throughout; it was performed after the
     v2 results were assembled.
  4. PTPRC is described as attenuating the association, not as accounting for or
     abolishing it. PTPRC may be a confounder, a mediator, or a composition
     proxy, and this design cannot distinguish them.
  5. With k = 2 cohorts, Cochran's Q has minimal power; heterogeneity is
     reported as "not detectable", never as "absent".
"""
import os
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(BASE, "14_Revision_v3_20260806/outputs")
os.makedirs(OUT, exist_ok=True)
RNG = np.random.default_rng(20260806)
N_BOOT = 2000

APM6 = ["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2"]
APM5 = ["HLA-A", "HLA-B", "B2M", "TAP1", "TAP2"]          # harmonised primary
NE4 = ["ASCL1", "NEUROD1", "POU2F3", "YAP1"]


def zm(df, genes):
    s = df[[g for g in genes if g in df.columns]]
    return ((s - s.mean()) / s.std(ddof=1)).mean(axis=1)


def boot_ci(x, y, n_boot=N_BOOT):
    x, y = np.asarray(x, float), np.asarray(y, float)
    idx = RNG.integers(0, len(x), (n_boot, len(x)))
    v = np.array([stats.spearmanr(x[i], y[i]).statistic for i in idx])
    v = v[~np.isnan(v)]
    return float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))


def fisher_ci(r, n, a=0.05):
    z, se = np.arctanh(r), 1 / np.sqrt(n - 3)
    k = stats.norm.ppf(1 - a / 2)
    return float(np.tanh(z - k * se)), float(np.tanh(z + k * se))


def loo(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    v = [stats.spearmanr(np.delete(x, i), np.delete(y, i)).statistic for i in range(len(x))]
    return float(np.min(v)), float(np.max(v))


def fit(y, dat, cov):
    Xd = pd.concat([dat["DLL3"].rename("DLL3")] + [dat[c] for c in cov], axis=1)
    Xs = (Xd - Xd.mean()) / Xd.std(ddof=1)
    ys = (y - y.mean()) / y.std(ddof=1)
    f = sm.OLS(ys, sm.add_constant(Xs)).fit()
    ci = f.conf_int().loc["DLL3"]
    return float(f.params["DLL3"]), float(ci[0]), float(ci[1]), float(f.pvalues["DLL3"]), float(f.rsquared)


# --- Cologne (cBioPortal sclc_ucologne_2015), cached matrix ---------------------
X = pd.read_csv(os.path.join(BASE, "12_External_Replication_20260806/"
                                   "Cologne_GSE_expression_matrix.csv"), index_col=0)
zero_var = [g for g in APM6 if X[g].nunique() <= 1]
COL = np.log2(X.clip(lower=0) + 1)
print(f"Cologne: n={len(COL)} tumours; zero-variance APM genes in the deposited matrix: {zero_var}")

# --- Jiang (GSE60052), tumours only ---------------------------------------------
g = pd.read_csv(os.path.join(BASE, "05_Source_Data/03_Public_Bulk_George_GSE60052/"
                                   "George_GSE60052_Log2_Normalized_RNAseq.tsv"),
                sep="\t", index_col=0)
g = g[~g.index.duplicated(keep="first")]
tum = [c for c in g.columns if "normal" not in c.lower()]
JIA = g.loc[:, tum].astype(float).T
print(f"Jiang:   n={len(JIA)} tumours (7 adjacent-normal samples excluded)")

COHORTS = {"Jiang (GSE60052)": JIA, "George/Cologne (cBioPortal)": COL}
rows, meta_in = [], {}

print("\n=== PRIMARY (harmonised 5-gene APM: HLA-A, HLA-B, B2M, TAP1, TAP2) ===")
for name, dat in COHORTS.items():
    apm = zm(dat, APM5)
    dll3 = dat["DLL3"]
    n = len(dat)
    r, p = stats.spearmanr(dll3, apm)
    fl, fh = fisher_ci(r, n)
    bl, bh = boot_ci(dll3, apm)
    l0, l1 = loo(dll3, apm)
    meta_in[name] = (r, n)
    print(f"  {name:<28} n={n:3d}  r={r:+.4f}  P={p:.3g}")
    print(f"      Fisher 95% CI [{fl:+.3f}, {fh:+.3f}]   bootstrap 95% CI [{bl:+.3f}, {bh:+.3f}]"
          f"   LOO [{l0:+.3f}, {l1:+.3f}]")
    rows.append(dict(cohort=name, n=n, apm_definition="5-gene (harmonised)",
                     model="Spearman DLL3~APM", estimate=round(r, 4),
                     ci_low_fisher=round(fl, 4), ci_high_fisher=round(fh, 4),
                     ci_low_boot=round(bl, 4), ci_high_boot=round(bh, 4),
                     p=p, loo_min=round(l0, 4), loo_max=round(l1, 4),
                     role="EXTERNAL REPLICATION (post hoc)"))

print("\n=== SENSITIVITY: 6-gene APM where available ===")
for name, dat, genes in [("Jiang (GSE60052)", JIA, APM6),
                         ("George/Cologne (cBioPortal)", COL, APM6)]:
    apm = zm(dat, genes)
    r, p = stats.spearmanr(dat["DLL3"], apm)
    used = [x for x in genes if x in dat.columns and dat[x].nunique() > 1]
    print(f"  {name:<28} 6-gene request -> genes with variance: {len(used)}  r={r:+.4f}  P={p:.3g}")
    rows.append(dict(cohort=name, n=len(dat), apm_definition=f"6-gene requested ({len(used)} informative)",
                     model="Spearman DLL3~APM", estimate=round(r, 4),
                     ci_low_fisher=None, ci_high_fisher=None, ci_low_boot=None, ci_high_boot=None,
                     p=p, loo_min=None, loo_max=None, role="sensitivity (APM definition)"))

print("\n=== Stepwise adjustment (harmonised 5-gene APM, identical models) ===")
MODELS = [("unadjusted", []), ("+ NE lineage", NE4), ("+ PTPRC", ["PTPRC"]),
          ("+ NE lineage + PTPRC", NE4 + ["PTPRC"])]
for name, dat in COHORTS.items():
    apm = zm(dat, APM5)
    print(f"  -- {name}")
    for label, cov in MODELS:
        b, lo, hi, p, r2 = fit(apm, dat, cov)
        print(f"     {label:<22} beta(DLL3)={b:+.4f} [{lo:+.4f}, {hi:+.4f}]  P={p:.4f}  R2={r2:.3f}"
              f"  {'' if p < .05 else '(NS)'}")
        rows.append(dict(cohort=name, n=len(dat), apm_definition="5-gene (harmonised)",
                         model=f"std OLS APM~DLL3 {label}", estimate=round(b, 4),
                         ci_low_fisher=round(lo, 4), ci_high_fisher=round(hi, 4),
                         ci_low_boot=None, ci_high_boot=None, p=round(p, 4),
                         loo_min=None, loo_max=None,
                         role="PRIMARY model" if label == "unadjusted" else "sensitivity (post hoc)"))
    rp, pp = stats.spearmanr(dat["DLL3"], dat["PTPRC"])
    print(f"     DLL3 ~ PTPRC (CD45)    r={rp:+.4f}  P={pp:.3g}")
    rows.append(dict(cohort=name, n=len(dat), apm_definition="-", model="Spearman DLL3~PTPRC",
                     estimate=round(rp, 4), ci_low_fisher=None, ci_high_fisher=None,
                     ci_low_boot=None, ci_high_boot=None, p=pp, loo_min=None, loo_max=None,
                     role="composition audit (post hoc)"))

# --- fixed-effect meta-analysis, Fisher z ---------------------------------------
rs = [meta_in[k][0] for k in COHORTS]
ns = [meta_in[k][1] for k in COHORTS]
z = np.arctanh(rs)
w = np.array([n - 3 for n in ns], float)
mu = float((z * w).sum() / w.sum())
se = float(np.sqrt(1 / w.sum()))
Q = float((w * (z - mu) ** 2).sum())
meta = dict(pooled_r=round(float(np.tanh(mu)), 4),
            ci_low=round(float(np.tanh(mu - 1.96 * se)), 4),
            ci_high=round(float(np.tanh(mu + 1.96 * se)), 4),
            p=float(2 * stats.norm.sf(abs(mu / se))),
            Q=round(Q, 4), df=1, p_heterogeneity=round(float(stats.chi2.sf(Q, 1)), 4),
            note="k=2; Cochran's Q has minimal power, so heterogeneity is reported as "
                 "not detectable rather than absent")
print(f"\n=== Fixed-effect meta (Fisher z, k=2, N={sum(ns)}) ===")
print(f"  pooled r = {meta['pooled_r']:+.4f}  95% CI [{meta['ci_low']:+.4f}, {meta['ci_high']:+.4f}]"
      f"  P = {meta['p']:.3g}")
print(f"  Q = {meta['Q']} (df=1), P = {meta['p_heterogeneity']}  -> heterogeneity NOT DETECTABLE (k=2, low power)")
rows.append(dict(cohort="Fixed-effect pooled (Jiang + Cologne)", n=sum(ns),
                 apm_definition="5-gene (harmonised)", model="pooled Spearman (Fisher z)",
                 estimate=meta["pooled_r"], ci_low_fisher=meta["ci_low"], ci_high_fisher=meta["ci_high"],
                 ci_low_boot=None, ci_high_boot=None, p=meta["p"], loo_min=None, loo_max=None,
                 role="EXTERNAL REPLICATION (post hoc)"))

S = pd.DataFrame(rows)
S.to_csv(os.path.join(OUT, "SuppTable_External_Replication_and_Leukocyte_Sensitivity_v3.csv"), index=False)
pd.DataFrame([meta]).to_csv(os.path.join(OUT, "Meta_TwoExternalCohorts_v3.csv"), index=False)
pd.DataFrame({"sampleId": COL.index, "DLL3_log2": COL["DLL3"].values,
              "APM_z_5gene": zm(COL, APM5).values}).to_csv(
    os.path.join(OUT, "Plotdata_Cologne_v3.csv"), index=False)
pd.DataFrame({"sampleId": JIA.index, "DLL3_log2": JIA["DLL3"].values,
              "APM_z_5gene": zm(JIA, APM5).values,
              "APM_z_6gene": zm(JIA, APM6).values}).to_csv(
    os.path.join(OUT, "Plotdata_Jiang_v3.csv"), index=False)
print(f"\nwritten -> {OUT}/SuppTable_External_Replication_and_Leukocyte_Sensitivity_v3.csv")
