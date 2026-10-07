"""v3 addendum (2026-09-12): does the unit of analysis change the post hoc result?

The Chan atlas holds 23 biospecimens from 19 donors: three donors contributed more than one
biospecimen (RU1108 three lung resection pieces; RU1144 and RU1181 a lung and a lymph-node
specimen each). The published analysis pools cells to the donor, which is the unit that avoids
pseudoreplication, but the manuscript described that unit as a "sample". This script recomputes
the principal post hoc endpoint and the central null under three alternative units so the
aggregation can be reported with evidence rather than asserted:

  * per biospecimen (n = 23; donors repeat, so the rows are not independent)
  * per donor, restricted to the 16 donors with a single biospecimen
  * one biospecimen per donor (n = 19), taking the biospecimen with the most epithelial cells

Both quantities are recomputed from the atlas at each unit, with the same definitions as s01:
the epithelial APM score is the mean of gene-wise z-scored per-unit mean expression of the five
harmonised APM genes, and the outcome is the percentage of CD8+ cells with detected TNFRSF9.

Output: outputs/SuppTable_Chan_UnitOfAnalysis_Sensitivity_v3.csv
Run with the project interpreter (~/miniforge3/bin/python3; Python 3.9.2, SciPy 1.13.1).
"""
import os

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(BASE, "14_Revision_v3_20260806/outputs")
ATLAS = os.path.join(BASE, "05_Source_Data/04_Public_scRNA_Chan_Atlas")
APM5 = ["HLA-A", "HLA-B", "B2M", "TAP1", "TAP2"]
CD8 = {"CD8+ Texh", "CD8+ Teff", "CD8+ Tmem"}

epi = ad.read_h5ad(os.path.join(ATLAS, "Chan_HTAN_SCLC_Epithelial_Cells_54313cells.h5ad"),
                   backed="r")
feat = epi.var.feature_name.astype(str)
ids = {s: epi.var.index[feat == s][0] for s in ["DLL3"] + APM5}
sub = epi[:, list(ids.values())].to_memory()
X = np.asarray(sub.X.todense()) if hasattr(sub.X, "todense") else np.asarray(sub.X)
e = pd.DataFrame(X, columns=list(ids.keys()))
e["donor"] = epi.obs.donor_id.astype(str).values
e["spec"] = epi.obs.HTAN_Biospecimen_ID.astype(str).values
epi.file.close()

imm = ad.read_h5ad(os.path.join(ATLAS, "Chan_HTAN_SCLC_Immune_Cells_16475cells.h5ad"),
                   backed="r")
featI = imm.var.feature_name.astype(str)
subI = imm[:, [imm.var.index[featI == "TNFRSF9"][0]]].to_memory()
XI = np.asarray(subI.X.todense()) if hasattr(subI.X, "todense") else np.asarray(subI.X)
i = pd.DataFrame({"TNFRSF9": XI[:, 0],
                  "donor": imm.obs.donor_id.astype(str).values,
                  "spec": imm.obs.HTAN_Biospecimen_ID.astype(str).values,
                  "ct": imm.obs.author_cell_type.astype(str).values})
imm.file.close()
i = i[i.ct.isin(CD8)]

n_spec_per_donor = e.groupby("donor", observed=True).spec.nunique()
print(f"atlas structure: {e.spec.nunique()} biospecimens from {e.donor.nunique()} donors; "
      f"{int((n_spec_per_donor > 1).sum())} donors with more than one biospecimen")


def build(key):
    """Per-unit epithelial APM score, mean epithelial DLL3 and TNFRSF9+ % of CD8."""
    m = e.groupby(key, observed=True)[APM5].mean()
    apm = m.apply(lambda c: (c - c.mean()) / c.std(ddof=1)).mean(axis=1)
    g = i.groupby(key, observed=True)
    out = pd.DataFrame({
        "epithelial_APM_z": apm,
        "DLL3_mean": e.groupby(key, observed=True)["DLL3"].mean(),
        "TNFRSF9pos_pct_of_CD8": g.TNFRSF9.apply(lambda v: 100 * (v > 0).mean()),
        "n_epithelial_cells": e.groupby(key, observed=True).size(),
        "n_CD8": g.size()})
    return out.dropna(subset=["epithelial_APM_z", "TNFRSF9pos_pct_of_CD8"])


by_spec = build("spec")
by_donor = build("donor")
by_spec["donor"] = (pd.concat([e[["donor", "spec"]], i[["donor", "spec"]]])
                    .drop_duplicates().set_index("spec").donor)
singles = by_donor.loc[n_spec_per_donor[n_spec_per_donor == 1].index]
one_each = by_spec.sort_values("n_epithelial_cells").groupby("donor", observed=True).tail(1)

rows = []
for unit, d in [("per donor (published primary unit)", by_donor),
                ("per biospecimen (donors repeat; not independent)", by_spec),
                ("per donor, single-biospecimen donors only", singles),
                ("one biospecimen per donor (most epithelial cells)", one_each)]:
    for x, label in [("epithelial_APM_z", "epithelial APM vs TNFRSF9+ %CD8"),
                     ("DLL3_mean", "tumour DLL3 vs TNFRSF9+ %CD8")]:
        r, p = stats.spearmanr(d[x], d.TNFRSF9pos_pct_of_CD8)
        rows.append({"comparison": label, "unit_of_analysis": unit, "n": len(d),
                     "spearman_r": round(float(r), 4), "p": float(p)})
        print(f"  {label:<34s} {unit:<50s} n = {len(d):2d}  r = {r:+.3f}  P = {p:.4f}")

res = pd.DataFrame(rows)
res["note"] = np.where(res.unit_of_analysis.str.startswith("per donor (published"),
                       "primary analysis reported in the manuscript", "sensitivity analysis")
res.to_csv(os.path.join(OUT, "SuppTable_Chan_UnitOfAnalysis_Sensitivity_v3.csv"), index=False)

spec_detail = by_spec.loc[by_spec.donor.isin(n_spec_per_donor[n_spec_per_donor > 1].index)]
spec_detail = spec_detail.sort_values(["donor", "n_epithelial_cells"])
spec_detail.to_csv(os.path.join(OUT, "SuppTable_Chan_MultiBiospecimen_Donors_v3.csv"))
print(f"\nCD8 cells per biospecimen in the three multi-biospecimen donors: "
      f"{int(spec_detail.n_CD8.min())} to {int(spec_detail.n_CD8.max())}")
print("written ->", OUT)
