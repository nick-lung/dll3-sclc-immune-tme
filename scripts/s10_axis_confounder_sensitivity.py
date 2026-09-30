"""v3: Technical-confounder sensitivity for the post hoc epithelial APM vs
TNFRSF9-defined CD8 abundance association.

Both quantities in that correlation are derived from the same single-cell object,
so a reviewer will reasonably ask whether a shared technical factor drives both.
The manuscript already adjusts the within-cell CD8 phenotype analysis for
sequencing depth, but the per-donor correlation was not adjusted for anything.
This script closes that gap.

Candidate shared drivers tested:
  * per-donor epithelial sequencing depth (median total UMI, median detected genes)
  * per-donor epithelial cell number
  * per-donor immune cell number
  * immune fraction of all profiled cells (i.e. how immune-rich the dissociation
    was), which is the single-cell analogue of the PTPRC confounder the paper
    identifies for the bulk cohorts
  * per-donor immune-compartment sequencing depth

Adjustment is by rank-based partial correlation: both variables and the
covariates are rank-transformed, each variable is regressed on the covariates by
OLS, and the Spearman correlation of the residuals is reported. This keeps the
estimator on the same rank scale as the unadjusted Spearman statistic.
"""
import os
import numpy as np
import pandas as pd
import anndata as ad
import statsmodels.api as sm
from scipy import stats

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(BASE, "14_Revision_v3_20260806/outputs")
ATLAS = os.path.join(BASE, "05_Source_Data/04_Public_scRNA_Chan_Atlas")
QC = ["total_counts", "n_genes_by_counts", "mito_frac"]
AX2 = ["TNFRSF9pos_pct_of_CD8", "TNFRSF9+PD1+_pct_of_CD8", "TNFRSF9+TIM3+_pct_of_CD8"]
X = "epithelial_APM_z"


def per_donor_qc(path, prefix):
    a = ad.read_h5ad(path, backed="r")
    obs = a.obs
    cols = [c for c in QC if c in obs.columns]
    d = obs.groupby("donor_id", observed=True)[cols].median()
    d.columns = [f"{prefix}_{c}" for c in cols]
    d[f"n_{prefix}"] = obs.groupby("donor_id", observed=True).size()
    return d


per = pd.read_csv(os.path.join(OUT, "Chan_PerDonor_with_epithelialAPM_v3.csv"), index_col=0)
epi = per_donor_qc(os.path.join(ATLAS, "Chan_HTAN_SCLC_Epithelial_Cells_54313cells.h5ad"), "epi")
imm = per_donor_qc(os.path.join(ATLAS, "Chan_HTAN_SCLC_Immune_Cells_16475cells.h5ad"), "imm")
d = per.join(epi, how="left").join(imm, how="left")
d["immune_fraction"] = d.n_imm / (d.n_imm + d.n_epi)


def partial_spearman(df, x, y, covs):
    sub = df[[x, y] + covs].dropna()
    R = sub.rank()
    M = sm.add_constant(R[covs])
    rx = sm.OLS(R[x], M).fit().resid
    ry = sm.OLS(R[y], M).fit().resid
    res = stats.spearmanr(rx, ry)
    return float(res.statistic), float(res.pvalue), len(sub)


MODELS = [
    ([], "unadjusted"),
    (["epi_total_counts"], "+ epithelial sequencing depth"),
    (["epi_n_genes_by_counts"], "+ epithelial detected genes"),
    (["n_epi"], "+ epithelial cell number"),
    (["n_imm"], "+ immune cell number"),
    (["immune_fraction"], "+ immune fraction of profiled cells"),
    (["imm_total_counts"], "+ immune-compartment depth"),
    (["epi_total_counts", "immune_fraction"], "+ depth and immune fraction"),
    (["epi_total_counts", "immune_fraction", "n_epi"], "+ depth, immune fraction, cell number"),
]

rows = []
print("Shared-technical-driver check for the post hoc axis-versus-axis association")
print("=" * 92)
for y in AX2:
    role = "PRINCIPAL" if y == AX2[0] else "supportive"
    print(f"\n--- {X} ~ {y}   [{role}] ---")
    for covs, label in MODELS:
        if covs and not all(c in d.columns for c in covs):
            continue
        if covs:
            r, p, n = partial_spearman(d, X, y, covs)
        else:
            sub = d[[X, y]].dropna()
            res = stats.spearmanr(sub[X], sub[y])
            r, p, n = float(res.statistic), float(res.pvalue), len(sub)
        flag = "" if p < 0.05 else "   <- not significant"
        print(f"    {label:<42} n={n:2d}  r={r:+.3f}  P={p:.4f}{flag}")
        rows.append(dict(y=y, endpoint_role=role, model=label,
                         n=n, partial_spearman_r=round(r, 4), p=round(p, 4),
                         covariates=";".join(covs) if covs else "none"))

print("\n" + "=" * 92)
print("Association of each candidate driver with the two correlated quantities")
print("=" * 92)
drv = []
for c in ["epi_total_counts", "epi_n_genes_by_counts", "n_epi", "n_imm",
          "immune_fraction", "imm_total_counts"]:
    if c not in d.columns:
        continue
    sub = d[[X, AX2[0], c]].dropna()
    r1 = stats.spearmanr(sub[c], sub[X])
    r2 = stats.spearmanr(sub[c], sub[AX2[0]])
    print(f"  {c:<26} vs epithelial APM r={r1.statistic:+.3f} (P={r1.pvalue:.3f}) | "
          f"vs TNFRSF9+ %CD8 r={r2.statistic:+.3f} (P={r2.pvalue:.3f})")
    drv.append(dict(driver=c, vs_epithelial_APM_r=round(float(r1.statistic), 4),
                    vs_epithelial_APM_p=round(float(r1.pvalue), 4),
                    vs_TNFRSF9pos_r=round(float(r2.statistic), 4),
                    vs_TNFRSF9pos_p=round(float(r2.pvalue), 4)))

pd.DataFrame(rows).to_csv(os.path.join(OUT, "SuppTable_Axis_ConfounderSensitivity_v3.csv"), index=False)
pd.DataFrame(drv).to_csv(os.path.join(OUT, "SuppTable_Axis_CandidateDrivers_v3.csv"), index=False)
d.to_csv(os.path.join(OUT, "Chan_PerDonor_with_QC_v3.csv"))
print(f"\nwritten -> {OUT}/SuppTable_Axis_ConfounderSensitivity_v3.csv")
print(f"written -> {OUT}/SuppTable_Axis_CandidateDrivers_v3.csv")
