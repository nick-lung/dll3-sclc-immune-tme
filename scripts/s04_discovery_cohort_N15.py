"""v3 P0-1: Institutional discovery cohort re-analysis after pathology adjudication.

Pathology adjudication (2026-08-06) confirmed that two of the 17 sequenced
specimens are large-cell neuroendocrine carcinoma (LC-NEC), not SCLC:
    21029B2866, 21029T2952
Both are therefore excluded from the primary analysis set, which becomes N = 15
pathologically confirmed SCLC.

This aligns the RNA-seq inclusion criterion with the mIF inclusion criterion
(the archival mIF batch was already excluded in part because 6 of its 15 patients
were LC-NEC). v2 stated "all 17 histologically confirmed SCLC; no exclusions",
which was incorrect.

The N=17 (diagnosis-inclusive) analysis is retained as a sensitivity analysis so
that the effect of the adjudication is fully visible.
"""
import os
import numpy as np
import pandas as pd
from scipy import stats
import statsmodels.api as sm

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(BASE, "14_Revision_v3_20260806/outputs")
os.makedirs(OUT, exist_ok=True)
SEED = 20260806
rng = np.random.default_rng(SEED)
N_BOOT = 2000

LC_NEC = ["21029B2866", "21029T2952"]          # pathology-adjudicated LC-NEC
APM = ["HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2"]
NE4 = ["ASCL1", "NEUROD1", "POU2F3", "YAP1"]
EFF = ["GZMB", "PRF1", "NKG7", "IFNG"]
EXH = ["PDCD1", "LAG3", "TIGIT", "HAVCR2", "CTLA4", "TOX"]
M2 = ["MRC1", "MARCO", "MSR1", "CD163"]


def boot_ci(x, y, n_boot=N_BOOT):
    x, y = np.asarray(x, float), np.asarray(y, float)
    idx = rng.integers(0, len(x), (n_boot, len(x)))
    v = np.array([stats.spearmanr(x[i], y[i]).statistic for i in idx])
    v = v[~np.isnan(v)]
    return float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5))


def loo(x, y):
    x, y = np.asarray(x, float), np.asarray(y, float)
    v = [stats.spearmanr(np.delete(x, i), np.delete(y, i)).statistic for i in range(len(x))]
    return float(np.min(v)), float(np.max(v))


def zmodule(df, genes, cols):
    gl = [x for x in genes if x in df.index]
    s = df.loc[gl, cols].astype(float)
    return s.sub(s.mean(axis=1), axis=0).div(s.std(axis=1, ddof=1), axis=0).mean(axis=0)


d = pd.read_excel(os.path.join(BASE, "05_Source_Data/01_Institutional_RNAseq/"
                                     "InHouse_RNAseq_Normalized_Gene_Expression.xlsx"),
                  sheet_name="Normalized.expression")
ALL17 = [c for c in d.columns if str(c).startswith("21029")]
d = d.set_index("Reference_GeneName")
d = d[~d.index.duplicated(keep="first")]
SCLC15 = [s for s in ALL17 if s not in LC_NEC]

assert len(ALL17) == 17 and len(SCLC15) == 15, (len(ALL17), len(SCLC15))
print(f"sequenced: {len(ALL17)}")
print(f"LC-NEC excluded by pathology adjudication: {LC_NEC}")
print(f"PRIMARY analysis set (pathologically confirmed SCLC): N={len(SCLC15)}")

rows = []


def analyse(cols, label, role):
    n = len(cols)
    dll3 = np.log2(d.loc["DLL3", cols].astype(float))
    ratio = (np.log2(d.loc["CD274", cols].astype(float))
             - np.log2(d.loc[APM, cols].astype(float)).mean(axis=0))
    print(f"\n--- {label} (n={n}) [{role}] ---")

    r, p = stats.spearmanr(dll3, ratio)
    lo, hi = boot_ci(dll3, ratio)
    l0, l1 = loo(dll3, ratio)
    print(f"  CD274/APM ratio ~ DLL3        r={r:+.4f}  P={p:.4f}  boot95%[{lo:+.3f},{hi:+.3f}]"
          f"  LOO[{l0:+.3f},{l1:+.3f}]")
    rows.append(dict(population=label, n=n, measure="CD274_APM_ratio~DLL3", estimate=round(r, 4),
                     ci_low=round(lo, 4), ci_high=round(hi, 4), p=round(p, 4),
                     extra=f"LOO {l0:+.4f} to {l1:+.4f}", role=role))

    med = dll3.median()
    u, pmed = stats.mannwhitneyu(ratio[dll3 >= med], ratio[dll3 < med])
    print(f"  median split (DLL3 mRNA)      U={u:.1f}  P={pmed:.4f}")
    rows.append(dict(population=label, n=n, measure="CD274_APM_ratio_median_split",
                     estimate=round(float(u), 1), ci_low=None, ci_high=None, p=round(pmed, 4),
                     extra="Mann-Whitney U", role=role))

    for lab, genes in [("APM_module_z", APM), ("M2_module_z", M2)]:
        rr, pp = stats.spearmanr(dll3, zmodule(d, genes, cols))
        print(f"  DLL3 ~ {lab:<16}      r={rr:+.4f}  P={pp:.4f}")
        rows.append(dict(population=label, n=n, measure=f"DLL3~{lab}", estimate=round(rr, 4),
                         ci_low=None, ci_high=None, p=round(pp, 4), extra="", role=role))

    for gene in ["CD274", "TNFRSF9", "PTPRC"]:
        if gene in d.index:
            rr, pp = stats.spearmanr(dll3, d.loc[gene, cols].astype(float))
            print(f"  DLL3 ~ {gene:<16}      r={rr:+.4f}  P={pp:.4f}")
            rows.append(dict(population=label, n=n, measure=f"DLL3~{gene}", estimate=round(rr, 4),
                             ci_low=None, ci_high=None, p=round(pp, 4), extra="", role=role))

    t = np.log2(d.loc["TNFRSF9", cols].astype(float) + 1)
    for mn, mg in [("effector", EFF), ("exhaustion", EXH)]:
        rr, pp = stats.spearmanr(t, zmodule(d, mg, cols))
        print(f"  TNFRSF9 ~ {mn:<13}      r={rr:+.4f}  P={pp:.4g}")
        rows.append(dict(population=label, n=n, measure=f"TNFRSF9~{mn}_module_z",
                         estimate=round(rr, 4), ci_low=None, ci_high=None, p=pp,
                         extra="tissue-level co-expression", role=role))

    # NE-adjusted OLS with full reporting (replaces the unreported ridge)
    y = (ratio - ratio.mean()) / ratio.std(ddof=1)
    X = pd.concat([dll3.rename("DLL3"),
                   np.log2(d.loc[["IFNG"] + NE4, cols].astype(float) + 1).T], axis=1)
    Xs = (X - X.mean()) / X.std(ddof=1)
    m = sm.OLS(y, sm.add_constant(Xs)).fit()
    b = float(m.params["DLL3"])
    clo, chi = [float(v) for v in m.conf_int().loc["DLL3"]]
    print(f"  ratio ~ DLL3 + IFNG + NE4     beta={b:+.4f} [{clo:+.4f}, {chi:+.4f}]"
          f"  P={float(m.pvalues['DLL3']):.4f}  R2={m.rsquared:.3f}  (OLS, fully reported)")
    rows.append(dict(population=label, n=n, measure="CD274_APM_ratio~DLL3+IFNG+NE4",
                     estimate=round(b, 4), ci_low=round(clo, 4), ci_high=round(chi, 4),
                     p=round(float(m.pvalues["DLL3"]), 4),
                     extra=f"std OLS; R2={m.rsquared:.4f}; replaces unreported ridge", role=role))

    pd.DataFrame({"sample": cols, "DLL3_log2": dll3.values,
                  "CD274_APM_ratio": ratio.values}).to_csv(
        os.path.join(OUT, f"Plotdata_Discovery_{label}_v3.csv"), index=False)
    return r, p


print("\n" + "=" * 88)
print("PRIMARY: pathologically confirmed SCLC only")
print("=" * 88)
r15, p15 = analyse(SCLC15, "SCLC15", "PRIMARY")

print("\n" + "=" * 88)
print("SENSITIVITY: diagnosis-inclusive set as reported in v2")
print("=" * 88)
r17, p17 = analyse(ALL17, "AllSequenced17", "sensitivity (diagnosis-inclusive)")

print("\n" + "=" * 88)
print(f"Effect of the pathology adjudication on the headline exploratory statistic:")
print(f"  N=17 (v2, incorrect)  r={r17:+.4f}  P={p17:.4f}")
print(f"  N=15 (v3, primary)    r={r15:+.4f}  P={p15:.4f}")
print("=" * 88)

# per-sample diagnosis table for the cohort-flow supplement
diag = pd.DataFrame({"sample_id": ALL17})
diag["adjudicated_diagnosis"] = np.where(diag.sample_id.isin(LC_NEC),
                                         "LC-NEC (large-cell neuroendocrine carcinoma)", "SCLC")
diag["in_primary_analysis_set"] = ~diag.sample_id.isin(LC_NEC)
diag["exclusion_reason"] = np.where(diag.sample_id.isin(LC_NEC),
                                    "Diagnosis: LC-NEC on pathology adjudication (2026-08-06)", "")
diag = diag.sort_values("sample_id")
diag.to_csv(os.path.join(OUT, "SuppTable_RNAseq_Diagnosis_Adjudication_v3.csv"), index=False)

pd.DataFrame(rows).to_csv(os.path.join(OUT, "Discovery_Statistics_N15_primary_v3.csv"), index=False)
print(f"\nwritten -> {OUT}/Discovery_Statistics_N15_primary_v3.csv")
print(f"written -> {OUT}/SuppTable_RNAseq_Diagnosis_Adjudication_v3.csv")
