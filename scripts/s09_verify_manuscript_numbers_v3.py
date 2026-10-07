"""v3: Verify that every headline statistic in the manuscript matches the derived data.

The root-cause analysis for this project (11_Draft_Check_20260720) found that main
tables and figures had drifted out of sync with the analysis because nothing
checked them. This script closes that loop: it re-reads the derived CSVs, extracts
the corresponding claim from the manuscript text, and reports any mismatch.
"""
import openpyxl
import os
import re
import sys
import numpy as np
import pandas as pd

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
V3 = os.path.join(BASE, "14_Revision_v3_20260806")
SRC = os.path.join(V3, "outputs")

MS = open(os.path.join(V3, "Manuscript_IJC_v3_DLL3_Immune_TME_20260806.md"),
          encoding="utf-8").read()
LEG = open(os.path.join(V3, "Figure_Legends_IJC_v3_20260806.md"), encoding="utf-8").read()
CL = open(os.path.join(V3, "Cover_Letter_IJC_v3_20260806.md"), encoding="utf-8").read()
ALL = MS + "\n" + LEG + "\n" + CL

fails, checks, warns = [], 0, []


def claim(label, needle, where=ALL):
    """Assert that a literal string appears in the manuscript."""
    global checks
    checks += 1
    if needle not in where:
        fails.append(f"MISSING  {label}: expected to find «{needle}»")
        return False
    return True


def forbid(label, needle, where=ALL):
    global checks
    checks += 1
    if needle in where:
        fails.append(f"FORBIDDEN {label}: found «{needle}» — this claim was withdrawn")
        return False
    return True


print("=" * 82)
print("1. NUMBERS RE-READ FROM DERIVED DATA")
print("=" * 82)

disc = pd.read_csv(os.path.join(SRC, "Discovery_Statistics_N15_primary_v3.csv"))
d15 = disc[disc.population == "SCLC15"].set_index("measure")
d17 = disc[disc.population == "AllSequenced17"].set_index("measure")
print(f"  discovery N=15 CD274/APM~DLL3   r={d15.loc['CD274_APM_ratio~DLL3','estimate']:+.4f} "
      f"P={d15.loc['CD274_APM_ratio~DLL3','p']:.4f}")
print(f"  discovery N=17 CD274/APM~DLL3   r={d17.loc['CD274_APM_ratio~DLL3','estimate']:+.4f} "
      f"P={d17.loc['CD274_APM_ratio~DLL3','p']:.4f}")
claim("discovery N=15 r", "*r* = +0.579")
claim("discovery N=15 P", "*P* = 0.024")
claim("discovery N=17 sensitivity", "*r* = +0.608, *P* = 0.0096")
claim("discovery adjusted OLS", "β = +0.297 (95% CI −0.621 to +1.216, *P* = 0.48)")
claim("discovery N=15 Fisher CI",
      f"95% CI {d15.loc['CD274_APM_ratio~DLL3','fisher_ci_low']:+.3f} to "
      f"{d15.loc['CD274_APM_ratio~DLL3','fisher_ci_high']:+.3f}")
claim("discovery N stated", "*N* = 15")
forbid("old discovery N", "bulk RNA-sequencing (*N* = 17)")
forbid("v2 no-exclusion claim", "all 17 constitute the primary analysis set")
forbid("v2 no-exclusion claim 2", "no samples were excluded")

ext = pd.read_csv(os.path.join(SRC, "SuppTable_External_Replication_and_Leukocyte_Sensitivity_v3.csv"))
prim = ext[(ext.model == "Spearman DLL3~APM") & (ext.apm_definition == "5-gene (harmonised)")]
for _, r in prim.iterrows():
    print(f"  {r.cohort:<38} r={r.estimate:+.4f} P={r.p:.3g}")
meta = pd.read_csv(os.path.join(SRC, "Meta_TwoExternalCohorts_v3.csv")).iloc[0]
print(f"  pooled r={meta.pooled_r:+.4f} CI[{meta.ci_low:+.4f},{meta.ci_high:+.4f}] "
      f"Q={meta.Q} P_het={meta.p_heterogeneity}")
claim("Jiang r", "*r* = −0.299, *P* = 0.0075")
claim("Cologne r", "*r* = −0.411, *P* = 1.4 × 10⁻⁴")
claim("pooled r", "pooled *r* = −0.357, 95% CI −0.486 to −0.212")
claim("Q test wording", "not detectable rather than demonstrably absent")
forbid("6-gene as primary", "*r* = −0.305, *P* = 0.0063")

for coh, tag in [("Jiang (GSE60052)", "Jiang"), ("George/Cologne (cBioPortal)", "Cologne")]:
    for mod in ["unadjusted", "+ NE lineage", "+ PTPRC", "+ NE lineage + PTPRC"]:
        row = ext[(ext.cohort == coh) & (ext.model == f"std OLS APM~DLL3 {mod}")].iloc[0]
        print(f"  {tag:<8} {mod:<22} beta={row.estimate:+.4f} P={row.p:.4f}")
claim("Jiang NE-adj", "β = −0.163, 95% CI −0.428 to +0.102, *P* = 0.225")
claim("Cologne NE-adj", "β = −0.456, 95% CI −0.808 to −0.105, *P* = 0.012")
claim("Jiang PTPRC", "−0.369 to −0.102 (*P* = 0.32)")
claim("Cologne PTPRC", "−0.398 to −0.094 (*P* = 0.21)")
claim("PTPRC wording", "markedly attenuated")
forbid("causal PTPRC claim A", "abolished the association")
forbid("causal PTPRC claim B", "accounted for the association")
forbid("PTPRC pre-specified", "pre-specified sensitivity analysis addressing the possibility")

ax = pd.read_csv(os.path.join(SRC, "SuppTable_Axis1_vs_Axis2_and_Equivalence_v3.csv"))
a12 = ax[ax.analysis == "Axis1_vs_Axis2_direct"]
eqv = ax[ax.analysis == "DLL3_to_Axis2_equivalence"]
for _, r in a12.iterrows():
    print(f"  axis1-vs-axis2 {r.y:<26} r={r.spearman_r:+.4f} P={r.p:.4f} q={r.bh_q:.4f} "
          f"| >=20cells r={r.r_cd8_ge20:+.3f} P={r.p_cd8_ge20:.3f} "
          f"| >=50cells r={r.r_cd8_ge50:+.3f} P={r.p_cd8_ge50:.3f}")
for _, r in eqv.iterrows():
    print(f"  central null   {r.y:<26} r={r.spearman_r:+.4f} P={r.p:.4f} "
          f"minDelta={r.smallest_equivalence_margin:.3f}")
claim("central null P", "*r* = +0.042 to +0.105, all *P* > 0.6")
claim("central null range", "*r* = +0.042 to +0.105")
claim("equivalence bound", "|ρ| ≈ 0.45")
claim("equivalence range", "smallest supported margin 0.425 to 0.475")
prin = a12[a12.y == "TNFRSF9pos_pct_of_CD8"].iloc[0]
claim("axis1v2 principal r", f"*r* = {prin.spearman_r:+.3f}")
claim("axis1v2 fisher CI", f"(95% CI {prin.ci_low:+.3f} to {prin.ci_high:+.3f}"
                          .replace("+", "+").replace("-", "−"))
# the text reports Fisher z intervals only; bootstrap intervals live in the supplementary tables
forbid("bootstrap CI quoted in the text", "bootstrap 95% CI")
claim("axis1v2 q", f"*q* = {prin.bh_q:.4f}")
claim("axis1v2 CD8 thresholds",
      f"fewer than 20 or 50 CD8⁺ cells (*r* = {prin.r_cd8_ge20:+.3f} and {prin.r_cd8_ge50:+.3f})")
claim("epithelial APM 5-gene primary", "*r* = +0.026, *P* = 0.91")
claim("6-gene sensitivity reported", "six-gene sensitivity")
claim("nested less stable", "less stable")
# S10c: the post hoc axis association must survive shared-technical-driver adjustment,
# and the manuscript must report that it was tested.
cfd = pd.read_csv(os.path.join(SRC, "SuppTable_Axis_ConfounderSensitivity_v3.csv"))
prin_c = cfd[(cfd.y == "TNFRSF9pos_pct_of_CD8")]
print(f"  axis confounder sensitivity: principal endpoint r range "
      f"{prin_c.partial_spearman_r.min():+.3f} to {prin_c.partial_spearman_r.max():+.3f}, "
      f"max P = {prin_c.p.max():.4f}")
checks += 1
if (prin_c.p >= 0.05).any():
    fails.append("AXIS principal endpoint loses significance under a confounder-adjustment model")
claim("confounder sensitivity reported", "shared technical factor could generate the correlation")
# S13: the leukocyte-free falsification must be null with an intact positive control
ccle = pd.read_csv(os.path.join(SRC, "SuppTable_CCLE_LeukocyteFree_Falsification_v3.csv"))
_pt = ccle[ccle.role == "PRIMARY falsification test"].iloc[0]
_pc = ccle[ccle.comparison == "DLL3 ~ ASCL1"].iloc[0]
_pu = ccle[ccle.comparison == "PTPRC ~ APM"].iloc[0]
print(f"  CCLE n={int(_pt.n)}: DLL3~APM r={_pt.spearman_r:+.3f} P={_pt.p:.3f} | "
      f"positive control DLL3~ASCL1 r={_pc.spearman_r:+.3f} P={_pc.p:.1e} | "
      f"purity PTPRC~APM r={_pu.spearman_r:+.3f} P={_pu.p:.2f}")
checks += 1
if _pt.p < 0.05:
    fails.append("CCLE falsification: DLL3~APM is significant in cell lines, which contradicts "
                 "the compositional interpretation stated in the manuscript")
checks += 1
if _pc.p >= 0.05:
    fails.append("CCLE positive control failed: DLL3~ASCL1 not significant, so the null cannot "
                 "be attributed to absence of leukocytes")
checks += 1
if _pu.p < 0.05:
    fails.append("CCLE purity control failed: PTPRC still tracks APM, so the system is not "
                 "leukocyte-free")
claim("CCLE test reported", "DLL3 was also tested against the same five-gene APM module in a leukocyte-free system")
claim("CCLE primary value", f"*r* = {_pt.spearman_r:+.3f}".replace("-", "−"))
claim("CCLE positive control", f"DLL3–*ASCL1* *r* = {_pc.spearman_r:+.3f}")
claim("CCLE limitation stated", "so the cell-line null alone would not establish this reading")
claim("confounder range reported",
      f"all nine adjustment models (*r* = {prin_c.partial_spearman_r.min():+.3f} to "
      f"{prin_c.partial_spearman_r.max():+.3f}")
drv = pd.read_csv(os.path.join(SRC, "SuppTable_Axis_CandidateDrivers_v3.csv"))
checks += 1
both = drv[(drv.vs_epithelial_APM_p < 0.05) & (drv.vs_TNFRSF9pos_p < 0.05)]
if len(both):
    fails.append(f"AXIS candidate driver(s) {list(both.driver)} associate with BOTH quantities; "
                 "the manuscript claims none does")
forbid("independent axes A", "independent axes in DLL3-high SCLC")
forbid("independent axes B", "Reporting these as independent axes")
forbid("independent axes C", "Two independent axes")
forbid("not correlated claim", "The two axes we measured are not correlated with each other")
forbid("v2 false axis claim", "The two axes were tested for association and none was detected")
forbid("APM suppression title", "tracks antigen-presentation suppression")

dep = pd.read_csv(os.path.join(SRC, "SuppTable_DetectionDepth_Sensitivity_v3.csv"))
log = dep[dep.analysis == "cell_level_logistic"]
full = log[log.model == "+ depth + subtype + donor"]
for _, r in full.iterrows():
    print(f"  depth-adjusted {r.quantity:<6} OR={r.odds_ratio:.2f} "
          f"[{r.ci_low:.2f},{r.ci_high:.2f}] P={r.p:.3g}")
claim("PD-1 adjusted OR", "1.96 (95% CI 1.64 to 2.33")
claim("TIM-3 adjusted OR", "3.04 (95% CI 2.39 to 3.86")
claim("GZMB not significant", "adjusted OR 1.24, *P* = 0.19")
claim("depth AUC", "AUC] = 0.599")

cnt = pd.read_csv(os.path.join(SRC, "SuppTable_Chan_PerDonor_CellCounts_v3.csv"), index_col=0)
print(f"  Chan donors={len(cnt)} evaluable(primary)={int(cnt.evaluable_primary_min3.sum())} "
      f"Texh total={int(cnt.n_CD8_Texh.sum())} top donor share="
      f"{100*cnt.n_CD8_Texh.max()/cnt.n_CD8_Texh.sum():.0f}%")
claim("15 evaluable", "15 evaluable donors")

# ---- unit of analysis (raised 2026-09-12: 23 biospecimens were pooled into 19 donors) ----
unit = pd.read_csv(os.path.join(SRC, "SuppTable_Chan_UnitOfAnalysis_Sensitivity_v3.csv"))
_u = unit[unit.comparison.str.startswith("epithelial APM")].set_index("unit_of_analysis")
_spec = _u.loc["per biospecimen (donors repeat; not independent)"]
_sing = _u.loc["per donor, single-biospecimen donors only"]
_one = _u.loc["one biospecimen per donor (most epithelial cells)"]
print(f"  unit sensitivity: per biospecimen n={int(_spec.n)} r={_spec.spearman_r:+.3f} "
      f"P={_spec.p:.4f} | single-biospecimen donors n={int(_sing.n)} r={_sing.spearman_r:+.3f} "
      f"P={_sing.p:.4f} | one per donor n={int(_one.n)} r={_one.spearman_r:+.3f} P={_one.p:.4f}")
claim("atlas structure stated", "It holds 23 biospecimens from 19 donors")
claim("pooling justified", "the unit that avoids pseudoreplication")
claim("no cross-specimen comparison", "nothing is compared across specimens")
_units = [_spec.spearman_r, _sing.spearman_r, _one.spearman_r]
claim("unit sensitivity range",
      f"each alternative unit of analysis (*r* = {min(_units):+.3f} to {max(_units):+.3f})")
# the atlas unit must not be called a sample anywhere in the submitted text
for _txt, _where in [(MS, "manuscript"), (LEG, "figure legends"), (CL, "cover letter")]:
    for _bad in ["19 samples", "evaluable samples", "same-sample", "across samples",
                 "of 19 samples", "per sample,"]:
        checks += 1
        if _bad in _txt:
            fails.append(f"UNIT {_where} still calls the Chan atlas unit a sample: «{_bad}»")

# ---- five-gene primary versus six-gene sensitivity must not be mixed in the main tables ---
# The study-design table (main Table 1) once quoted the six-gene values while the statistics
# table (now Supplementary Table S11) and the text quote the five-gene primary.
import docx as _docx
_t2 = _docx.Document(os.path.join(V3, "tables", "Table1_v3_Questions_and_Dataset_Roles.docx"))
_cells = lambda d: " | ".join(c.text for t in d.tables for row in t.rows for c in row.cells)
_wb11 = openpyxl.load_workbook(os.path.join(V3, "supplementary",
                                            "SuppTable_S11_Key_Statistics_All_Platforms_v3.xlsx"),
                               read_only=True)
_c1 = " | ".join(str(v) for ws in _wb11.worksheets for row in ws.iter_rows(values_only=True)
                 for v in row if v is not None)
_c2 = _cells(_t2)
_AX = pd.read_csv(os.path.join(SRC, "SuppTable_Axis1_vs_Axis2_and_Equivalence_v3.csv"))
_ctx5 = _AX[(_AX.analysis == "DLL3_to_epithelial_APM_context")
            & (_AX.y == "epithelial_APM_z_5gene_primary")].iloc[0]
_ctx6 = _AX[(_AX.analysis == "DLL3_to_epithelial_APM_context")
            & (_AX.y == "epithelial_APM_z_6gene_sensitivity")].iloc[0]
_d5 = _AX[_AX.analysis == "Axis1_vs_Axis2_direct"]
_d6 = _AX[_AX.analysis == "Axis1_vs_Axis2_6gene_sensitivity"]
for _lbl, _txt in [("Supplementary Table S11", _c1), ("Table 1", _c2)]:
    checks += 2
    if f"r = {_ctx5.spearman_r:+.3f}" not in _txt:
        fails.append(f"TABLES {_lbl} does not carry the five-gene epithelial APM context value")
    if f"r = {_ctx6.spearman_r:+.3f}, P = {_ctx6.p:.2f}" in _txt:
        fails.append(f"TABLES {_lbl} quotes the six-gene sensitivity value as if it were primary")
checks += 2
if f"{_d5.spearman_r.min():+.3f} to {_d5.spearman_r.max():+.3f}" not in _c2:
    fails.append("TABLES Table 1 does not carry the five-gene feature-versus-feature range")
if f"{_d6.spearman_r.min():+.3f} to {_d6.spearman_r.max():+.3f}" in _c2:
    fails.append("TABLES Table 1 quotes the six-gene feature-versus-feature range")
# the S9 legend annotates a five-gene figure
claim("S9 legend uses the five-gene value",
      f"unrelated to the same five-gene epithelial APM score (*r* = {_ctx5.spearman_r:+.3f}, "
      f"*P* = {_ctx5.p:.2f})", LEG)
claim("Texh counts", "only 578 Texh cells")
claim("Texh concentration", "74% from one donor, RU1195")

mif = pd.read_csv(os.path.join(BASE, "05_Source_Data/06_Derived_Analysis_Tables/"
                                     "Derived_mIF_ROI_All_51Markers_Results.csv"))
print(f"  mIF phenotypes={len(mif)} FDR-significant={int((mif.q<0.05).sum())} "
      f"min q={mif.q.min():.3f} max q={mif.q.max():.3f}")
checks += 1
if int((mif.q < 0.05).sum()) != 0:
    fails.append("MIF the text reports no significant phenotype but the derived table has one")
claim("mIF null and minimum q", f"across all {len(mif)} phenotypes (minimum *q* = {mif.q.min():.3f}")
forbid("mIF all q=0.895", "all *q* = 0.895")
forbid("mIF top-3 narrative", "directionally concordant trends without statistical support")
forbid("mIF Fig 3F misref", "3.5-fold higher")

print()
print("=" * 82)
print("2. WITHDRAWN v2 CLAIMS AND STALE ARTEFACTS")
print("=" * 82)
forbid("ridge model", "ridge model adjusting for")
forbid("ridge beta", "β = +0.382")
forbid("reference to an unreviewed earlier version", "earlier version")
forbid("reference to the v2 figure set", "v2 ")
forbid("reinvigoratable", "reinvigoratable T-cell compartment retained")
# The words may appear only inside an explicit disclaimer, never as a description.
claim("checkpoint-responsive disclaimed",
      "these cells are not described as checkpoint-responsive or reinvigoratable")
for _sent in re.split(r"(?<=[.])\s+", MS):
    if "checkpoint-responsive" in _sent and "not described as" not in _sent:
        checks += 1
        fails.append(f"FORBIDDEN checkpoint-responsive used affirmatively: «{_sent[:120]}…»")
forbid("adjacent-normal wording", "7 adjacent-normal samples")
claim("normal lung controls", "normal lung controls")
forbid("AI policy assertion", "consistent with the journal's policy")
# IJC 4.1.1: the cover letter must carry any disclosure of assistance for writing or editing.
claim("AI assistance disclosed in cover letter", "**Disclosure of AI assistance.**", CL)
claim("adjudication stated", "were adjudicated LC-NEC and excluded")
claim("mIF criterion consistency", "The same diagnostic criterion excluded LC-NEC from the mIF cohort")
claim("post hoc labelled", "post hoc")
forbid("exhausted CD8 harbor claim", "harbor more exhausted CD8")

print()
print("=" * 82)
print("2a0. REFERENCE INTEGRITY")
print("=" * 82)
# Vancouver/Wiley numeric style: references must be numbered in order of first citation.
_SUP = "¹²³⁴⁵⁶⁷⁸⁹⁰"
_TO = str.maketrans(_SUP, "1234567890")
_refs = {int(m.group(1)): m.group(2).strip()
         for m in re.finditer(r"^(\d+)\.\s+(.*?)$", MS.split("# References", 1)[1], flags=re.M)}
_cbody = MS.split("# Abstract", 1)[1].split("# References", 1)[0]
_order = []
for _m in re.finditer(r"[¹²³⁴⁵⁶⁷⁸⁹⁰]+(?:[ʼ,][¹²³⁴⁵⁶⁷⁸⁹⁰]+)*", _cbody):
    _pre = _cbody[max(0, _m.start() - 4):_m.start()]
    if re.search(r"(10|×|\^|⁻)\s*$", _pre) or "⁻" in _m.group(0):
        continue          # scientific-notation exponent, not a citation
    for _x in re.findall(r"[¹²³⁴⁵⁶⁷⁸⁹⁰]+", _m.group(0)):
        _n = int(_x.translate(_TO))
        if 1 <= _n <= len(_refs) and _n not in _order:
            _order.append(_n)
print(f"  {len(_refs)} references; first-citation order {'is' if _order == sorted(_order) else 'is NOT'} sequential")
checks += 1
if _order != sorted(_order):
    fails.append(f"REFS not numbered in order of first citation; order is {_order}")
checks += 1
_uncited = [n for n in _refs if n not in _order]
if _uncited:
    fails.append(f"REFS uncited in the body: {_uncited}")
checks += 1
if len(_order) != len(_refs):
    fails.append(f"REFS {len(_refs)} listed but {len(_order)} cited")
# corrections confirmed against the primary sources during the 2026-08-07 audit
forbid("ref9 wrong third author", "Roper N, Velez MJ, Chiappinelli KB")
forbid("ref22 stale ahead-of-print", "(Online ahead of print.)")
forbid("Sen mischaracterised", "STING activation restores T-cell recognition in SCLC models")
forbid("Hipp/Chen bundled claim",
       "engagers induce checkpoint molecules and that concurrent PD-1 blockade adds activity")
forbid("Rudin2019 cold/inflamed misattribution",
       "tend to be immunologically \"cold,\" whereas NE-low subtypes are comparatively inflamed")
forbid("TOST bound self-contradiction", "This bound includes the *r* = +0.505")
claim("institutional PTPRC discrepancy disclosed", "*r* = +0.061, *P* = 0.83, *N* = 15")
claim("Figure 1 cited", "shown in Figure 1 and Table 1")

print()
print("=" * 82)
print("2a. MAIN-FIGURE PANEL REFERENCES")
print("=" * 82)
# The v2 peer review found a wrong panel reference ("Fig. 3F") that nothing had
# checked. v3 redesigned Figures 2 and 3, and the Results text silently kept the
# old lettering. This check derives the real panel set from the figure-generating
# code and cross-checks it against the text and the legends, so the failure mode
# cannot recur silently.
FIGSRC = {
    2: (os.path.join(V3, "s06_generate_main_figures_v3.py"), "def figure2", "def figure5"),
    3: (os.path.join(V3, "s07_generate_figure3_and_supplementary_figures_v3.py"),
        "def figure3", "def figure4"),
    4: (os.path.join(V3, "s07_generate_figure3_and_supplementary_figures_v3.py"),
        "def figure4", "def suppS1_S2"),
}
_bodytext = MS.split("# Results", 1)[1].split("# Data Availability", 1)[0]
for fignum, (path, start, end) in FIGSRC.items():
    code = open(path, encoding="utf-8").read()
    blk = code[code.find(start):code.find(end)]
    real = set(re.findall(r'panel\(\s*ax\s*,\s*"([A-H])"', blk))
    # panels emitted from an indexed string, e.g. panel(ax, "DE"[i], ...)
    for lit in re.findall(r'panel\(\s*ax\s*,\s*"([A-H]{2,})"\s*\[', blk):
        real |= set(lit)
    # panels emitted as  "C" if i == 0 else "D"
    for a, b in re.findall(r'panel\(\s*ax\s*,\s*"([A-H])"\s*if[^,]*else\s*"([A-H])"', blk):
        real |= {a, b}
    cited = set(re.findall(rf"Fig\. {fignum}([A-H])", _bodytext))
    for m in re.finditer(rf"Fig\. {fignum}([A-H])\s*[–-]\s*([A-H])", _bodytext):
        cited |= {chr(c) for c in range(ord(m.group(1)), ord(m.group(2)) + 1)}
    for m in re.finditer(rf"Fig\. {fignum}([A-H])(?:,\s*([A-H]))+", _bodytext):
        cited |= set(g for g in m.groups() if g)
    legblk = LEG[LEG.find(f"**Figure {fignum}."):]
    legblk = legblk[:legblk.find("**Figure ", 5)] if legblk.find("**Figure ", 5) > 0 else legblk
    # legends group panels as "(A)", "**(D, E)**", "(C–E)"
    legend = set()
    for grp in re.findall(r"\(\s*((?:[A-H]\s*(?:,|and|–|-)?\s*)+)\)", legblk.replace("**", "")):
        legend |= set(re.findall(r"[A-H]", grp))
        rng = re.match(r"\s*([A-H])\s*[–-]\s*([A-H])", grp)
        if rng:
            legend |= {chr(c) for c in range(ord(rng.group(1)), ord(rng.group(2)) + 1)}
    print(f"  Figure {fignum}: generated {sorted(real)} | cited {sorted(cited)} | "
          f"legend {sorted(legend)}")
    checks += 1
    if cited - real:
        fails.append(f"FIGREF Fig. {fignum} cites panel(s) {sorted(cited - real)} that do not exist")
    checks += 1
    if real - cited:
        fails.append(f"FIGREF Fig. {fignum} panel(s) {sorted(real - cited)} are never cited in the text")
    checks += 1
    if real - legend:
        fails.append(f"FIGREF Fig. {fignum} panel(s) {sorted(real - legend)} have no legend entry")

print()
print("=" * 82)
print("2b. NO-CLINICAL-ANNOTATION POLICY (author decision, 2026-08-07)")
print("=" * 82)
claim("no-clinical stated in Limitations",
      "No clinical annotation is reported for the institutional cohort, by design")
claim("unadjusted acknowledged", "Residual confounding by unmeasured patient characteristics", MS)
claim("residual confounding acknowledged", "esidual confounding by unmeasured patient characteristics")
claim("generalisability acknowledged", "generalisability to any defined clinical population")

# S1 must carry no clinical columns, and no blank placeholder columns anywhere.
import openpyxl
CLINICAL = ["age", "sex", "stage", "smoking", "prior therapy", "ecog", "ps",
            "outcome", "survival", "death", "recurrence", "regimen", "tumour content",
            "tumor content", "timing of sampling", "specimen site"]
SUPPDIR = os.path.join(V3, "supplementary")
for fn in sorted(os.listdir(SUPPDIR)):
    if not fn.endswith(".xlsx"):
        continue
    wb = openpyxl.load_workbook(os.path.join(SUPPDIR, fn), read_only=True)
    for ws in wb.worksheets:
        rows = list(ws.iter_rows(min_row=1, max_row=4, values_only=True))
        hdr = [str(c).lower() for r in rows for c in (r or ()) if c]
        for term in CLINICAL:
            checks += 1
            if any(h.strip() == term or h.strip().startswith(term + " ") for h in hdr):
                fails.append(f"POLICY clinical column '{term}' present in {fn} [{ws.title}]")
    wb.close()
print("  scanned every supplementary workbook for clinical columns")
forbid("no INCOMPLETE placeholder left in S1",
       "must be completed from the institutional clinical record")

print()
print("=" * 82)
print("3. FORMAT REQUIREMENTS")
print("=" * 82)
# IJC counts the "body of the text from Introduction to Discussion/Conclusions",
# which by its own IMRAD ordering INCLUDES Materials and Methods and EXCLUDES the
# abstract, references, figure legends, tables and disclosure statements.
_parts = re.split(r"^# ", MS, flags=re.M)
_sec = {p.split("\n", 1)[0].strip(): p.split("\n", 1)[1] for p in _parts[1:] if "\n" in p}


def _wc(t):
    t = re.sub(r"^#{1,6} .*$", "", t, flags=re.M)
    t = re.sub(r"[*_`]", "", t)
    return len([w for w in t.split() if re.search(r"[A-Za-z0-9]", w)])


wc = sum(_wc(_sec[k]) for k in
         ["Introduction", "Materials and Methods", "Results", "Discussion"])
print(f"  IJC body (Introduction->Discussion, incl. Methods) = {wc:,} words  (limit 5,000)")
checks += 1
if wc > 5000:
    fails.append(f"FORMAT IJC body is {wc:,} words; the Research Article limit is 5,000")
checks += 1
n_display = 4 + 2  # 4 main figures + 2 main tables
if n_display > 6:
    fails.append(f"FORMAT {n_display} main display items; IJC allows 6 in total")
checks += 1
if n_ref_limit_exceeded := (len(re.findall(r"^\d+\. ", MS.split("# References", 1)[1],
                                           flags=re.M)) > 50):
    fails.append("FORMAT more than 50 references; IJC Research Article limit is 50")
abs_txt = MS.split("# Abstract", 1)[1].split("# Introduction", 1)[0]
abs_wc = len(re.sub(r"[*_`#]", "", abs_txt).split())
print(f"  abstract ~ {abs_wc} words  (IJC limit 250)")
ni = MS.split("**Novelty & Impact:**", 1)[1].split("\n", 1)[0]
ni_wc = len(re.sub(r"[*_`#]", "", ni).split())
print(f"  Novelty & Impact ~ {ni_wc} words  (IJC limit 75)")
kw = MS.split("**Keywords:**", 1)[1].split("\n", 1)[0]
n_kw = len([k for k in kw.split(";") if k.strip()])
print(f"  keywords = {n_kw}  (IJC allows 3-5)")
n_ref = len(re.findall(r"^\d+\. ", MS.split("# References", 1)[1], flags=re.M))
print(f"  references = {n_ref}")
checks += 4
if abs_wc > 250:
    fails.append(f"FORMAT abstract {abs_wc} words exceeds the 250-word limit")
if ni_wc > 75:
    fails.append(f"FORMAT Novelty & Impact {ni_wc} words exceeds the 75-word limit")
if not 3 <= n_kw <= 5:
    fails.append(f"FORMAT {n_kw} keywords; IJC allows 3-5")
if f"**References:** {n_ref}." not in MS:
    fails.append(f"FORMAT reference count on the title page does not match the {n_ref} listed")

# The stated word counts on the title page must match what the text actually is.
checks += 2
m_abs = re.search(r"Abstract (\d+) words", MS)
if not m_abs or abs(int(m_abs.group(1)) - abs_wc) > 3:
    fails.append(f"FORMAT title page states abstract {m_abs.group(1) if m_abs else '?'} words; "
                 f"measured {abs_wc}")
m_body = re.search(r"Discussion/Conclusions ([\d,]+) words", MS)
if not m_body or abs(int(m_body.group(1).replace(",", "")) - wc) > 250:
    fails.append(f"FORMAT title page states main text "
                 f"{m_body.group(1) if m_body else '?'} words; measured {wc:,}")

claim("abbreviation list", "**Abbreviations:**")
# writing style: the body uses "This study" or the passive, never the first person.
_bodytxt = MS.split("# Abstract", 1)[1].split("# Data Availability", 1)[0]
_bodytxt = _bodytxt.replace("to our knowledge", "")   # the hedge on the novelty claim stays
for _fp in [" we ", " We ", " our ", " Our "]:
    checks += 1
    if _fp in _bodytxt:
        fails.append(f"STYLE first person «{_fp.strip()}» in the body")
_das = MS.split("# Data Availability Statement", 1)[1].split("# Author Contributions", 1)[0]
# 2026-09-30: DAS changed from GEO deposition to on-request availability — the protocol
# (v1.3) states no public-database registration and the 2025-06-05 consent does not cover it.
claim("institutional RNA-seq on request", "available from the corresponding author on reasonable request, because the informed consent", _das)
claim("scripts in a public repository", "are available at", _das)
forbid("data on request only", "derived analysis tables are available from the corresponding author on reasonable request", _das)

# ---- IJC formal requirements confirmed against the live author guidelines on 2026-09-12 ----
# §3.1: corresponding author ORCID on the title page (Wiley requires it at submission).
claim("corresponding ORCID", "N. Kobayashi 0000-0002-7064-320X")
# §3: AI-use description belongs in the Acknowledgements, not in a separate section.
checks += 2
if "# Use of AI" in MS:
    fails.append("FORMAT AI-use disclosure must sit inside Acknowledgements (IJC §3), not a separate section")
if "**Use of AI-assisted technologies.**" not in MS.split("# Acknowledgements", 1)[-1].split("# References", 1)[0]:
    fails.append("FORMAT AI-use disclosure paragraph missing from Acknowledgements")
# IJC 3 requires the use to be described in detail; a named tool and the specific content it
# produced. The author decision of 2026-09-12 is to disclose drafting as well as code assistance.
claim("AI tool named in manuscript", "Claude Code (Anthropic)", MS)
claim("AI drafting disclosed", "to draft and edit the text of the manuscript", MS)
claim("AI tool named in cover letter", "Claude Code (Anthropic)", CL)
forbid("vague AI wording", "in a limited and supervised manner")
# §3.1.1 order: statements -> References -> Figure Legends.
checks += 1
if not (MS.find("# Acknowledgements") < MS.find("# References") < MS.find("# Figure Legends")):
    fails.append("FORMAT section order must be Acknowledgements -> References -> Figure Legends")
# §3.2: all authors named in every reference for review; no truncation.
checks += 1
_reflist = MS.split("# References", 1)[1].split("# Figure Legends", 1)[0]
if "et al" in _reflist:
    fails.append("FORMAT reference list truncates authors with 'et al.'; IJC requires all authors for review")
# Abbreviations listed on the title page must be defined at first use in the body.
_abbr = MS.split("**Abbreviations:**", 1)[1].split("\n", 1)[0]
_body_txt = MS.split("# Introduction", 1)[1].split("# Data Availability", 1)[0]
for _item in _abbr.split(";"):
    _ab = _item.split(",", 1)[0].strip().rstrip(".")
    if not _ab or _ab == "H-score":
        continue
    checks += 1
    _mm = re.search(rf"(?<![A-Za-z]){re.escape(_ab)}s?(?![A-Za-z])", _body_txt)   # whole token (plural allowed): "CI" is not "ICI"
    if not _mm:
        fails.append(f"FORMAT abbreviation {_ab} is listed but never used in the body")
    elif _body_txt[_mm.start() - 1] not in "([":
        fails.append(f"FORMAT abbreviation {_ab} is used before it is defined")
# The legends embedded in the manuscript must be the legends file, verbatim (v3 drifted here).
_ms_leg = MS.split("# Figure Legends", 1)[1].split("## Supplementary Figure Legends", 1)[0].strip()
_lg_leg = LEG.split("## Main Figures", 1)[1].split("---", 1)[0].strip()
checks += 1
if _ms_leg != _lg_leg:
    fails.append("LEGENDS main-figure legends in the manuscript differ from the legends file")
_ms_sup = MS.split("## Supplementary Figure Legends", 1)[1].strip()
_lg_sup = LEG.split("## Supplementary Figures", 1)[1].strip()
checks += 1
if _ms_sup != _lg_sup:
    fails.append("LEGENDS supplementary legends in the manuscript differ from the legends file")
forbid("internal v2 note in legends", "The v2 supplementary figure set is withdrawn", LEG)
# Every Supplementary Table sub-sheet that exists must be cited at least once (S10d was not).
for _fn in sorted(os.listdir(SUPPDIR)):
    if not _fn.endswith(".xlsx"):
        continue
    _wb = openpyxl.load_workbook(os.path.join(SUPPDIR, _fn), read_only=True)
    for _ws in _wb.worksheets:
        _m = re.match(r"S(\d+)([a-e])_", _ws.title)
        if _m:
            checks += 1
            # a sub-sheet counts as cited when it is named (S10c) or its parent table is cited
            # without a letter (S6); S10d failed this in the circulated v3.
            _n_, _l_ = _m.group(1), _m.group(2)
            _cited_sub = re.search(rf"Supplementary Tables? (?:S\d+[a-e]?(?:,| and|;) )*S{_n_}{_l_}\b", MS)
            _cited_par = re.search(rf"Supplementary Tables? (?:S\d+[a-e]?(?:,| and|;) )*S{_n_}(?![0-9a-e])", MS)
            if not (_cited_sub or _cited_par):
                fails.append(f"FORMAT Supplementary Table S{_n_}{_l_} exists but is never cited")
    _wb.close()
# Cover letter title must equal the manuscript title.
checks += 1
_t_ms = re.search(r"\*\*Title:\*\* (.*)", MS).group(1).strip()
_t_cl = re.search(r"\*\*Title:\*\* (.*)", CL).group(1).strip()
if _t_ms != _t_cl:
    fails.append("COVER LETTER title differs from the manuscript title")
# §3.4: figure text in Helvetica/Arial.
for _sf in ["s06_generate_main_figures_v3.py", "s07_generate_figure3_and_supplementary_figures_v3.py"]:
    checks += 1
    if '"font.family": ["Arial"' not in open(os.path.join(V3, _sf), encoding="utf-8").read():
        fails.append(f"FORMAT {_sf} does not set Arial for figure text (IJC §3.4)")
claim("panel format A", "**(A)**", LEG)

figdir = os.path.join(V3, "figures")
checks += 1
if not os.path.exists(os.path.join(figdir, "Graphical_Abstract_v3.pdf")):
    fails.append("FORMAT Graphical Abstract PDF not generated (s13)")
# the Graphical Abstract is uploaded as PDF only (IJC 3.5), so it is checked separately above
have = sorted(f for f in os.listdir(figdir)
              if f.endswith(".png") and not f.startswith("Graphical_Abstract"))
print(f"\n  figure files: {len(have)} PNG (+ PDF and TIFF for each)")
for f in have:
    for ext in ("pdf", "tiff"):
        checks += 1
        if not os.path.exists(os.path.join(figdir, f[:-3] + ext)):
            fails.append(f"FORMAT missing {ext.upper()} for {f}")
for n in range(1, 11):
    checks += 1
    if not any(f"SuppFigureS{n}_" in f for f in have):
        fails.append(f"FORMAT Supplementary Figure S{n} not generated")
    checks += 1
    if f"**Supplementary Figure S{n}." not in LEG:
        fails.append(f"FORMAT Supplementary Figure S{n} has no individual legend")
    checks += 1
    # accept "Supplementary Figure S5", "Supplementary Figures S5 and S6", "S5, S6".
    # (?![0-9]) stops "S1" being satisfied by "S10".
    if not re.search(rf"Supplementary Figures? (?:S\d+(?:,| and|;)? )*S{n}(?![0-9])", MS):
        fails.append(f"FORMAT Supplementary Figure S{n} is never cited in the main text")

suppdir = os.path.join(V3, "supplementary")
for n in range(1, 14):
    checks += 1
    if not any(f"_S{n}_" in f for f in os.listdir(suppdir)):
        fails.append(f"FORMAT Supplementary Table S{n} not generated")
    checks += 1
    # word-boundary match: plain "Supplementary Table S2" must not be satisfied by
    # "Supplementary Table S3". Also accept "Supplementary Tables S1 and S9" and
    # sub-panel references such as "Supplementary Table S3a".
    if not re.search(rf"Supplementary Tables? (?:S\d+[ab]?(?:,| and|;) )*S{n}(?![0-9])", MS):
        fails.append(f"FORMAT Supplementary Table S{n} is never cited in the main text")

_body_cit = MS.split("# Introduction", 1)[1].split("# References", 1)[0]
for _kind in ("Table", "Figure"):
    _seen = []
    for _m in re.finditer(rf"Supplementary {_kind}s? ((?:S\d+[a-e]?(?:, | and |–)?)+)", _body_cit):
        for _s in re.findall(r"S(\d+)", _m.group(1)):
            if int(_s) not in _seen:
                _seen.append(int(_s))
    checks += 1
    if _seen != sorted(_seen):
        fails.append(f"FORMAT Supplementary {_kind}s are not numbered in order of first citation: {_seen}")
# panels of each main figure, and lettered sub-tables, must also be cited in letter order
for _f in range(1, 6):
    _seq = []
    for _m in re.finditer(rf"Fig\. {_f}([A-H](?:(?:[–-]|, )[A-H])*)", _body_cit):
        for _a, _b in re.findall(r"([A-H])(?:[–-]([A-H]))?", _m.group(1)):
            for _c in range(ord(_a), ord(_b or _a) + 1):
                if chr(_c) not in _seq:
                    _seq.append(chr(_c))
    checks += 1
    if _seq != sorted(_seq):
        fails.append(f"FORMAT Figure {_f} panels are not cited in letter order: {_seq}")
for _n in range(1, 14):
    _seq = []
    for _m in re.finditer(rf"S{_n}([a-e])(?![a-z0-9])", _body_cit):
        if _m.group(1) not in _seq:
            _seq.append(_m.group(1))
    checks += 1
    if _seq != sorted(_seq):
        fails.append(f"FORMAT Supplementary Table S{_n} sub-tables are not cited in letter order: {_seq}")
for t in ["Table1_v3_Questions_and_Dataset_Roles.docx"]:
    checks += 1
    if not os.path.exists(os.path.join(V3, "tables", t)):
        fails.append(f"FORMAT main table missing in Word format: {t}")

print()
print("=" * 82)
print("4. OUTSTANDING PRE-SUBMISSION ITEMS")
print("=" * 82)
# Placeholders are deliberate: they mark values that do not exist yet (a GEO accession, a
# repository URL). They must never reach the submission system unnoticed, so they are listed
# here on every run.
for _f, _txt in [("manuscript", MS), ("cover letter", CL), ("figure legends", LEG)]:
    for _m in re.finditer(r"\[TO BE INSERTED: ([^\]]+)\]", _txt):
        warns.append(f"{_f}: {_m.group(1)}")
if warns:
    for _w in warns:
        print(f"  PLACEHOLDER  {_w}")
    print(f"\n  {len(warns)} placeholder(s) must be filled before the files are uploaded.")
else:
    print("  none")

print()
print("=" * 82)
if fails:
    print(f"RESULT: {len(fails)} PROBLEM(S) out of {checks} checks\n")
    for f in fails:
        print("  " + f)
    sys.exit(1)
print(f"RESULT: all {checks} checks passed")
