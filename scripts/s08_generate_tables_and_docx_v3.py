"""v3: Main Tables 1-2 as Word documents, plus the manuscript and figure-legend DOCX.

IJC requires main tables in doc/docx/rtf; xlsx is for supplementary material only.
Tables 1 and 2 are therefore emitted as Word tables, not spreadsheets.

Table 1 changes required by peer review and implemented here:
  * the mIF row reads "Tested; none significant (0 of 51)", not "Not tested";
  * the DLL3-Texh fraction row is no longer the KEY TEST for the link between the
    two features - the central link is the abundance of the TNFRSF9-defined
    phenotype, and the axis-versus-axis test is listed separately;
  * bulk TNFRSF9 co-expression is labelled supportive tissue-level evidence, not
    validation of a single-cell state;
  * every N reflects the pathology adjudication (discovery N = 15).
"""
import os
import pandas as pd
from docx import Document
from docx.shared import Pt, Inches, RGBColor
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.enum.table import WD_TABLE_ALIGNMENT
import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
V3 = os.path.join(BASE, "14_Revision_v3_20260806")
OUT = os.path.join(V3, "tables")
SRC = os.path.join(V3, "outputs")
os.makedirs(OUT, exist_ok=True)

# Read the axis-versus-axis and equivalence results rather than transcribing them.
# Transcribed values are how the v1/v2 tables drifted away from the analysis; the
# root-cause report for this project (11_Draft_Check_20260720) traced every main
# table discrepancy to a number that had been typed rather than generated.
_AX = pd.read_csv(os.path.join(SRC, "SuppTable_Axis1_vs_Axis2_and_Equivalence_v3.csv"))
_A12 = _AX[_AX.analysis == "Axis1_vs_Axis2_direct"].set_index("y")
_CTX = _AX[_AX.analysis == "DLL3_to_epithelial_APM_context"].set_index("y")


def _fmt_p(p):
    return f"{p:.4f}" if p >= 0.001 else f"{p:.1e}"


def ax_row(key, label, role, status_fn):
    r = _A12.loc[key]
    return (f"Feature 1 vs feature 2 (post hoc){'' if role.startswith('Post') else ''}",
            f"Chan atlas ({int(r.n)} donors)",
            f"Epithelial APM ~ {label}, Spearman",
            f"r = {r.spearman_r:+.3f}",
            f"{r.ci_low:+.3f} to {r.ci_high:+.3f}",
            f"{_fmt_p(r.p)} (q = {r.bh_q:.4f})",
            role,
            status_fn(r))


_PRIN = _A12.loc["TNFRSF9pos_pct_of_CD8"]
_A12_R_MIN, _A12_R_MAX = _A12.spearman_r.min(), _A12.spearman_r.max()
_EPI_CTX = _CTX.loc["epithelial_APM_z_5gene_primary"]
_DISC = pd.read_csv(os.path.join(SRC, "Discovery_Statistics_N15_primary_v3.csv"))
_DISC = _DISC[(_DISC.population == "SCLC15") & (_DISC.measure == "CD274_APM_ratio~DLL3")].iloc[0]


def styled_doc(landscape=False):
    d = Document()
    st = d.styles["Normal"]
    st.font.name = "Times New Roman"
    st.font.size = Pt(10)
    for s in d.sections:
        if landscape:
            s.page_width, s.page_height = s.page_height, s.page_width
        s.left_margin = s.right_margin = Inches(0.7)
    return d


def add_table(doc, df, title, notes, colwidths=None):
    p = doc.add_paragraph()
    r = p.add_run(title)
    r.bold = True
    r.font.size = Pt(10)
    t = doc.add_table(rows=1, cols=len(df.columns))
    t.style = "Table Grid"
    t.alignment = WD_TABLE_ALIGNMENT.CENTER
    for i, c in enumerate(df.columns):
        cell = t.rows[0].cells[i]
        cell.text = ""
        run = cell.paragraphs[0].add_run(str(c))
        run.bold = True
        run.font.size = Pt(8)
    for _, row in df.iterrows():
        cells = t.add_row().cells
        for i, v in enumerate(row):
            cells[i].text = ""
            run = cells[i].paragraphs[0].add_run("" if pd.isna(v) else str(v))
            run.font.size = Pt(8)
    if colwidths:
        for row in t.rows:
            for i, w in enumerate(colwidths):
                row.cells[i].width = Inches(w)
    doc.add_paragraph()
    for n in notes:
        p = doc.add_paragraph()
        run = p.add_run(n)
        run.font.size = Pt(8)
        p.paragraph_format.space_after = Pt(2)
    return doc


# ------------------------------------------------------------------ SUPPLEMENTARY TABLE S11
T1 = pd.DataFrame([
    # Feature 1
    ("Feature 1: lower bulk APM signal", "Jiang tumours (n = 79)", "DLL3 ~ APM module, Spearman",
     "r = −0.299", "−0.488 to −0.083", "0.0075", "External replication (primary)", "Significant"),
    ("Feature 1: lower bulk APM signal", "George/Cologne (n = 81)", "DLL3 ~ APM module, Spearman",
     "r = −0.411", "−0.578 to −0.212", "1.4 × 10⁻⁴", "External replication (primary)", "Significant"),
    ("Feature 1: pooled", "Jiang + Cologne (N = 160)", "Fixed-effect pooled Spearman",
     "r = −0.357", "−0.486 to −0.212", "3.6 × 10⁻⁶", "External replication", "Significant"),
    ("Feature 1: NE-adjusted", "Jiang tumours (n = 79)", "Standardized OLS β(DLL3), + NE lineage",
     "β = −0.163", "−0.428 to +0.102", "0.225", "Subtype sensitivity", "Not significant"),
    ("Feature 1: NE-adjusted", "George/Cologne (n = 81)", "Standardized OLS β(DLL3), + NE lineage",
     "β = −0.456", "−0.808 to −0.105", "0.012", "Subtype sensitivity", "Significant"),
    ("Feature 1: leukocyte-adjusted", "Jiang tumours (n = 79)", "Standardized OLS β(DLL3), + PTPRC",
     "β = −0.102", "−0.304 to +0.101", "0.321", "Composition sensitivity (post hoc)",
     "Not significant"),
    ("Feature 1: leukocyte-adjusted", "George/Cologne (n = 81)",
     "Standardized OLS β(DLL3), + PTPRC", "β = −0.094", "−0.240 to +0.053", "0.207",
     "Composition sensitivity (post hoc)", "Not significant"),
    ("Feature 1: tumour-cell level", f"Chan atlas ({int(_EPI_CTX.n)} donors)",
     "DLL3 ~ epithelial APM (5-gene), Spearman", f"r = {_EPI_CTX.spearman_r:+.3f}",
     "not estimated", f"{_EPI_CTX.p:.2f}", "Estimand contrast", "Not significant"),
    ("Feature 1: discovery (exploratory)", "Institutional (N = 15 SCLC)",
     "DLL3 ~ CD274/APM ratio, Spearman", "r = +0.579",
     f"{_DISC.fisher_ci_low:+.3f} to {_DISC.fisher_ci_high:+.3f}", "0.024",
     "Exploratory; not replicated", "Significant"),
    ("Feature 1: discovery, adjusted", "Institutional (N = 15 SCLC)",
     "Standardized OLS β(DLL3), + IFNG + NE lineage", "β = +0.297", "−0.621 to +1.216", "0.48",
     "Exploratory", "Not significant"),
    # Feature 2
    ("Feature 2: phenotype exists", "Chan atlas (15 evaluable)",
     "PD-1⁺ in TNFRSF9⁺ vs TNFRSF9⁻ CD8, paired", "median Δ +18.5 pp", "rank-biserial +0.97",
     "1.8 × 10⁻⁴", "Primary", "Significant"),
    ("Feature 2: phenotype exists", "Chan atlas (15 evaluable)",
     "TIM-3⁺ in TNFRSF9⁺ vs TNFRSF9⁻ CD8, paired", "median Δ +26.0 pp", "rank-biserial +1.00",
     "6.1 × 10⁻⁵", "Primary", "Significant"),
    ("Feature 2: depth-adjusted", "Chan atlas (5,285 CD8 cells)",
     "PD-1⁺ odds ratio, + depth + subtype + donor", "OR = 1.96", "1.64 to 2.33", "9.8 × 10⁻¹⁴",
     "Technical sensitivity", "Significant"),
    ("Feature 2: depth-adjusted", "Chan atlas (5,285 CD8 cells)",
     "TIM-3⁺ odds ratio, + depth + subtype + donor", "OR = 3.04", "2.39 to 3.86", "1.0 × 10⁻¹⁹",
     "Technical sensitivity", "Significant"),
    ("Feature 2: effector arm", "Chan atlas (15 evaluable)",
     "GZMB⁺ in TNFRSF9⁺ vs TNFRSF9⁻ CD8, paired", "median Δ +12.9 pp", "rank-biserial +0.58",
     "0.048", "Exploratory; threshold-dependent", "Not robust"),
    ("Feature 2: compartment", "Chan atlas (8 evaluable in both)",
     "Tmem minus Texh effect, TIM-3, within donor", "median +0.7 pp", "not estimated", "1.00",
     "Interaction test", "Not significant"),
    ("Feature 2: supportive tissue-level", "Institutional (N = 15) / Jiang (n = 79)",
     "Bulk TNFRSF9 ~ exhaustion module, Spearman", "r = +0.832 / +0.436", "not estimated",
     "1.2 × 10⁻⁴ / 5.9 × 10⁻⁵",
     "Supportive tissue-level co-expression; NOT validation of a single-cell state", "Significant"),
    # The central link
    ("CENTRAL LINK: DLL3 vs feature 2", "Chan atlas (19 donors)",
     "DLL3 ~ TNFRSF9⁺PD-1⁺ (% of CD8), Spearman", "r = +0.072", "−0.395 to +0.510", "0.77",
     "KEY TEST", "No association detected"),
    ("CENTRAL LINK: DLL3 vs feature 2", "Chan atlas (19 donors)",
     "DLL3 ~ TNFRSF9⁺ (% of CD8), Spearman", "r = +0.105", "−0.366 to +0.534", "0.67",
     "KEY TEST", "No association detected"),
    ("CENTRAL LINK: DLL3 vs feature 2", "Chan atlas (19 donors)",
     "DLL3 ~ TNFRSF9⁺TIM-3⁺ (% of CD8), Spearman", "r = +0.042", "−0.420 to +0.487", "0.86",
     "KEY TEST", "No association detected"),
    ("CENTRAL LINK: equivalence bound", "Chan atlas (19 donors)",
     "Smallest supported equivalence margin (TOST)", "|ρ| = 0.425 to 0.475", "descriptive",
     "—", "Bounds the null", "Moderate association NOT excluded"),
    # Axis vs axis  (values read from the derived table, never transcribed)
    ax_row("TNFRSF9pos_pct_of_CD8", "TNFRSF9⁺ (% of CD8)",
           "Post hoc direct test — principal endpoint", lambda r: "Significant"),
    ("Feature 1 vs feature 2, cell-count sensitivity",
     f"Chan atlas ({int(_PRIN.n_cd8_ge20)} and {int(_PRIN.n_cd8_ge50)} donors)",
     "Same, restricted to ≥20 and ≥50 CD8⁺ cells",
     f"r = {_PRIN.r_cd8_ge20:+.3f} / {_PRIN.r_cd8_ge50:+.3f}", "not estimated",
     f"{_PRIN.p_cd8_ge20:.4f} / {_PRIN.p_cd8_ge50:.3f}",
     "Robustness of the principal endpoint", "Significant at both thresholds"),
    ax_row("TNFRSF9+TIM3+_pct_of_CD8", "TNFRSF9⁺TIM-3⁺ (% of CD8)",
           "Supportive nested endpoint",
           lambda r: f"Not robust (P = {r.p_cd8_ge50:.3f} at ≥50 cells)"),
    ax_row("TNFRSF9+PD1+_pct_of_CD8", "TNFRSF9⁺PD-1⁺ (% of CD8)",
           "Supportive nested endpoint",
           lambda r: f"Not robust (P = {r.p_cd8_ge50:.3f} at ≥50 cells)"),
    # mIF
    ("mIF phenotypes", "Institutional (10 patients, 193 ROIs)",
     "All 51 derived phenotypes vs continuous DLL3 H-score, patient-level",
     "0 of 51 significant", "minimum q = 0.895", "—",
     "Tested; none significant", "Not significant"),
], columns=["Finding", "Dataset", "Measure", "Effect", "95% CI", "P", "Role", "Status"])

# The cross-platform statistics table is a supplementary table (S11); the main text keeps
# only the study-design table. Same layout as s05: title row, table from row 3, footnotes below.
_S11 = os.path.join(V3, "supplementary", "SuppTable_S11_Key_Statistics_All_Platforms_v3.xlsx")
_S11_NOTES = ["APM, antigen-presentation machinery; CI, confidence interval; NE, neuroendocrine; "
           "OLS, ordinary least squares; OR, odds ratio; pp, percentage points; ROI, region of "
           "interest; TOST, two one-sided tests.",
           "The institutional discovery cohort is N = 15 pathologically adjudicated SCLC. Two of "
           "the 17 sequenced specimens were adjudicated as large-cell neuroendocrine carcinoma "
           "and excluded (Supplementary Tables S4 and S3a).",
           "The APM module uses the same five genes (HLA-A, HLA-B, B2M, TAP1, TAP2) in both "
           "external cohorts.",
           "Analyses labelled post hoc were performed after the primary results were assembled.",
           "The mIF row reports the complete 51-phenotype family. No phenotype is singled out, "
           "because none is significant; the full result is in Supplementary Table S12 and "
           "Supplementary Figure S10.",
           "Bulk TNFRSF9 co-expression is supportive tissue-level evidence. It does not "
           "demonstrate co-expression within individual cells and is not validation of the "
           "single-cell phenotype.",
           "For the feature-1-versus-feature-2 test, only the single-gated TNFRSF9⁺ abundance "
           "endpoint is robust to minimum CD8 cell counts. The doubly-gated endpoints are "
           "directionally concordant but are estimated from progressively smaller cell pools and "
           "lose significance at a ≥50-cell threshold; they are reported as supportive only."]
with pd.ExcelWriter(_S11, engine="openpyxl") as w:
    pd.DataFrame([["Supplementary Table S11. All key statistics across platforms, with the role and "
                   "inferential status of each."]]).to_excel(w, sheet_name="S11", index=False,
                                                            header=False, startrow=0)
    T1.to_excel(w, sheet_name="S11", index=False, startrow=2)
    pd.DataFrame([[n] for n in _S11_NOTES]).to_excel(w, sheet_name="S11", index=False, header=False,
                                                     startrow=2 + len(T1) + 3)
print("  SuppTable_S11_Key_Statistics_All_Platforms_v3.xlsx")

# ------------------------------------------------------------------ TABLE 1
T2 = pd.DataFrame([
    ("Feature 1. DLL3-high SCLC shows a lower bulk antigen-presentation signal.",
     "Institutional RNA-seq (N = 15)", "Discovery / exploratory",
     "CD274/APM ratio rises with DLL3 (r = +0.579); direct DLL3–APM correlation not significant"),
    ("", "Jiang GSE60052 (n = 79 tumours)", "External replication (primary)",
     "r = −0.299, P = 0.0075"),
    ("", "George/Cologne (n = 81 tumours)", "External replication (independent)",
     "r = −0.411, P = 1.4 × 10⁻⁴; pooled r = −0.357"),
    ("", "Both external cohorts, + PTPRC", "Composition sensitivity (post hoc)",
     "Association markedly attenuated in both (P = 0.32 and 0.21)"),
    ("", "Chan atlas, epithelial cells", "Estimand contrast",
     f"No tumour-cell-level association (r = {_EPI_CTX.spearman_r:+.3f}, "
     f"P = {_EPI_CTX.p:.2f})"),
    ("Feature 2. A TNFRSF9-defined, PD-1/TIM-3-marker-high CD8 phenotype is present in SCLC.",
     "Chan atlas, within-donor paired (15 evaluable donors)", "Primary",
     "PD-1 +18.5 pp (P = 1.8 × 10⁻⁴); TIM-3 +26.0 pp (P = 6.1 × 10⁻⁵)"),
    ("", "Chan atlas, cell-level depth-adjusted models", "Technical sensitivity",
     "Survives adjustment for library size, subtype and donor (OR 1.96 and 3.04)"),
    ("", "Institutional and Jiang bulk RNA-seq", "Supportive tissue-level",
     "Bulk TNFRSF9 tracks effector and exhaustion modules; not single-cell co-expression"),
    ("Central question. Does baseline DLL3 predict the baseline abundance of feature 2?",
     "Chan atlas (19 donors)", "KEY TEST",
     "No association detected (r = +0.042 to +0.105, all P > 0.6)"),
    ("", "Chan atlas, TOST", "Bounds the null",
     "Equivalence supported only within |ρ| ≈ 0.45; a moderate association is not excluded"),
    ("", "GSE319155 (n = 7)", "Supportive concordance",
     "Directionally concordant null (r = 0.000)"),
    ("Are the two features independent of each other?",
     "Chan atlas (19 donors)", "Post hoc direct test",
     f"No. They are positively correlated (r = {_A12_R_MIN:+.3f} to {_A12_R_MAX:+.3f}, "
     f"all P < 0.05)"),
    ("Orthogonal phenotype context. mIF.", "Institutional (10 patients, 193 ROIs)",
     "Descriptive only; supports no inference",
     "0 of 51 phenotypes significant after correction across the full family"),
], columns=["Question or feature", "Dataset", "Role in the study", "Result"])

doc = styled_doc(landscape=True)
add_table(doc, T2,
          "Table 1. Study questions, the role of each dataset, and what each contributes.",
          ["TOST, two one-sided tests.",
           "Blank cells in the first column continue the feature or question named in the row "
           "above.",
           "The discovery and mIF cohorts are not independent: all 10 mIF patients are among the "
           "15 RNA-sequenced SCLC patients."],
          colwidths=[2.2, 1.9, 1.9, 3.6])
doc.save(os.path.join(OUT, "Table1_v3_Questions_and_Dataset_Roles.docx"))
print("  Table1_v3_Questions_and_Dataset_Roles.docx")


# ------------------------------------------------------------------ MANUSCRIPT DOCX
from md2docx_v3 import md_to_docx


md_to_docx(os.path.join(V3, "Manuscript_IJC_v3_DLL3_Immune_TME_20260806.md"),
           os.path.join(V3, "Manuscript_IJC_v3_DLL3_Immune_TME_20260806.docx"))
print("\ndone ->", OUT)
