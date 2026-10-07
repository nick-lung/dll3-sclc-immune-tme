"""v3 P0-5: Regenerate the complete supplementary table set from v3 source data.

The v2 supplementary set is superseded in full. It contained statements
incompatible with the v2 main text (S4), a stale primary-analysis N (S5), and a
cohort definition that the 2026-08-06 pathology adjudication has now corrected.
Tables are numbered S1-S13 (S11 is written by s08); each carries a title row and explicit footnotes.
"""
import os
import numpy as np
import pandas as pd

BASE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
OUT = os.path.join(BASE, "14_Revision_v3_20260806/supplementary")
SRC = os.path.join(BASE, "14_Revision_v3_20260806/outputs")
DER = os.path.join(BASE, "05_Source_Data/06_Derived_Analysis_Tables")
os.makedirs(OUT, exist_ok=True)

LC_NEC = ["21029B2866", "21029T2952"]


def write(name, sheets):
    """sheets: {sheet_name: (title, dataframe, [footnotes])}"""
    path = os.path.join(OUT, name)
    with pd.ExcelWriter(path, engine="openpyxl") as w:
        for sheet, (title, df, notes) in sorted(sheets.items()):   # sub-tables in letter order
            pd.DataFrame([[title]]).to_excel(w, sheet_name=sheet, index=False,
                                             header=False, startrow=0)
            df.to_excel(w, sheet_name=sheet, index=False, startrow=2)
            pd.DataFrame([[n] for n in notes]).to_excel(
                w, sheet_name=sheet, index=False, header=False, startrow=2 + len(df) + 3)
    print(f"  wrote {name}")


meta = pd.read_csv(os.path.join(
    BASE, "05_Source_Data/01_Institutional_RNAseq/CellCarta_Reference/"
          "SCLC_Prospective_17samples_Metadata_20260622.csv"))

# ------------------------------------------------------------------ S2
hs = pd.read_excel(os.path.join(
    BASE, "13_Review_Packet_20260806/04_補足資料/"
          "SuppTable1_Institutional_Patient_Characteristics.xlsx"))
allids = set(meta.SampleID)
hmap = {}
for pid, h, sp, nr in zip(hs.Patient_ID, hs.DLL3_Hscore, hs.Specimen_type, hs.N_ROIs):
    for pref in ("21029B", "21029T"):
        if f"{pref}{pid}" in allids:
            hmap[f"{pref}{pid}"] = (h, sp, nr)

# S1 reports assay-level specimen metadata only. Patient demographics, stage,
# smoking history, treatment and outcome are deliberately NOT reported: this
# study is an assay-level analysis of archived specimens and draws no clinical
# inference, so no clinical annotation is presented for any patient. The
# limitation this creates is stated explicitly in the manuscript Limitations.
s1 = pd.DataFrame({"Sample ID": meta.SampleID})
s1["Sequencing-provider histology label"] = meta.SampleOrigin.values
s1["Adjudicated diagnosis"] = np.where(s1["Sample ID"].isin(LC_NEC), "LC-NEC", "SCLC")
s1["In RNA-seq primary analysis set"] = np.where(s1["Sample ID"].isin(LC_NEC), "No", "Yes")
s1["DV200"] = meta.DV200.values
s1["DLL3 IHC H-score"] = s1["Sample ID"].map(lambda x: hmap[x][0] if x in hmap else None)
s1["DLL3 IHC group"] = s1["DLL3 IHC H-score"].map(
    lambda v: "" if pd.isna(v) else ("High" if v >= 100 else "Low"))
s1["Specimen procedure (mIF)"] = s1["Sample ID"].map(
    lambda x: hmap[x][1] if x in hmap else None)
s1["mIF ROIs analysed"] = s1["Sample ID"].map(lambda x: hmap[x][2] if x in hmap else None)
s1["In mIF analysis set"] = np.where(s1["mIF ROIs analysed"].notna(), "Yes", "No")
s1 = s1.sort_values("Sample ID").reset_index(drop=True)

write("SuppTable_S2_Specimen_and_Assay_Metadata_v3.xlsx", {"S2": (
    "Supplementary Table S2. Specimen and assay metadata for the institutional cohort "
    "(17 specimens sequenced; 15 in the RNA-seq primary analysis set after pathology "
    "adjudication; 10 of these also in the mIF analysis set).",
    s1,
    ["Adjudicated diagnosis: independent review of the original pathology reports and slides by "
     "the institutional pathology service, performed without reference to any analysis result.",
     "21029B2866 and 21029T2952 were adjudicated as large-cell neuroendocrine carcinoma (LC-NEC) "
     "and are excluded from the RNA-seq primary analysis set.",
     "DLL3 IHC H-score, IHC group, specimen procedure and mIF ROI count are populated only for "
     "the 10 patients in the mIF analysis set; the DLL3 IHC group threshold is H-score >= 100.",
     "NO CLINICAL ANNOTATION IS REPORTED. Age, sex, stage, smoking history, treatment and outcome "
     "are deliberately not presented, because this study is an assay-level analysis of archived "
     "specimens that draws no clinical inference and performs no outcome, subgroup or "
     "clinically-adjusted analysis. Every reported estimate is unadjusted for patient "
     "characteristics. The resulting limitation - that residual confounding by unmeasured "
     "clinical variables cannot be assessed - is stated in the manuscript Limitations.",
     "The 'SCLC subtype' column present in the v2 version of this table is removed. All 10 "
     "patients had been labelled SCLC-I by an automated classifier, including specimens with "
     "ASCL1 above 7,000, which indicates classifier failure rather than a biological finding."])})

# ------------------------------------------------------------------ S8
s2 = pd.DataFrame([
    ("APM (antigen-presentation machinery)", "HLA-A, HLA-B, B2M, TAP1, TAP2",
     "PRIMARY. Harmonized 5-gene definition, used identically in both external cohorts and in "
     "the discovery cohort."),
    ("APM, 6-gene sensitivity", "HLA-A, HLA-B, HLA-C, B2M, TAP1, TAP2",
     "Sensitivity only. HLA-C has zero variance across all 81 George/Cologne samples, so a 6-gene "
     "module is not computable there; using 6 genes in one cohort and 5 in the other would make "
     "the cohorts non-comparable."),
    ("Effector", "GZMB, PRF1, NKG7, IFNG", "Secondary"),
    ("Exhaustion", "PDCD1, LAG3, TIGIT, HAVCR2, CTLA4, TOX", "Secondary"),
    ("NE lineage (adjustment set)", "ASCL1, NEUROD1, POU2F3, YAP1", "Covariate set"),
    ("Leukocyte content (adjustment)", "PTPRC", "Post hoc covariate"),
    ("CD274/APM ratio", "log2(CD274) - mean(log2(APM genes))",
     "Exploratory, discovery cohort only"),
], columns=["Module", "Genes", "Role"])
write("SuppTable_S8_Gene_Module_Definitions_v3.xlsx", {"S8": (
    "Supplementary Table S8. Gene module definitions.", s2,
    ["Module score = mean of gene-wise z-scored log2 expression, z-scored within the analysis "
     "population.",
     "CHANGED FROM v2: the APM module is now the harmonized 5-gene definition in all cohorts."])})

# ------------------------------------------------------------------ S12
mif = pd.read_csv(os.path.join(DER, "Derived_mIF_ROI_All_51Markers_Results.csv"))
mif = mif.sort_values("p").reset_index(drop=True)
mif.insert(0, "Rank by P", np.arange(1, len(mif) + 1))
mif.columns = ["Rank by P", "Phenotype", "Spearman r", "P", "Mean density, DLL3-high",
               "Mean density, DLL3-low", "Fold change", "q (Benjamini-Hochberg)",
               "FDR significant"]
write("SuppTable_S12_mIF_All51_Full_Results_v3.xlsx", {"S12": (
    "Supplementary Table S12. Complete multiplex-immunofluorescence results: all 51 derived "
    "phenotypes, patient-level area-weighted mean densities, 10 patients.", mif,
    [f"0 of {len(mif)} phenotypes reached FDR significance. Minimum q = "
     f"{mif['q (Benjamini-Hochberg)'].min():.3f}; maximum q = "
     f"{mif['q (Benjamini-Hochberg)'].max():.3f}.",
     "The Benjamini-Hochberg family is all 51 phenotypes, corrected in a single pass.",
     "The complete table is reported so that no subset is selected after inspecting the ranking. "
     "Ordering by P is presentational only and confers no significance on the top rows.",
     "DLL3 H-score is modelled continuously in the Spearman analysis. The DLL3-high and DLL3-low "
     "density columns are descriptive summaries at the median split and are not the basis of any "
     "test."])})

# ------------------------------------------------------------------ S5
qc = meta.rename(columns={
    "SampleID": "Sample ID", "CellCarta_ID": "CellCarta ID", "TissueType": "Tissue type",
    "Raw Number of Reads": "Raw reads (M)",
    "Percentage Uniquely Mapped Reads": "% uniquely mapped",
    "Percentage mRNA Bases": "% mRNA bases", "Percentage Duplicates": "% duplicates",
    "BAM_Filename": "BAM filename"})
qc["Adjudicated diagnosis"] = np.where(qc["Sample ID"].isin(LC_NEC), "LC-NEC", "SCLC")
qc["In primary analysis set"] = np.where(qc["Sample ID"].isin(LC_NEC), "No", "Yes")
qc = qc[["Sample ID", "CellCarta ID", "Tissue type", "Adjudicated diagnosis",
         "In primary analysis set", "DV200", "Raw reads (M)", "% uniquely mapped",
         "% mRNA bases", "% duplicates", "BAM filename"]]
write("SuppTable_S5_RNAseq_QC_Metrics_v3.xlsx", {"S5": (
    "Supplementary Table S5. Per-sample RNA-sequencing quality metrics and adjudicated diagnosis "
    "(CellCarta P1947).", qc,
    ["17 specimens were sequenced; all passed the pre-specified DV200 >= 40 threshold.",
     "THE PRIMARY ANALYSIS SET IS N = 15. 21029B2866 and 21029T2952 were adjudicated as LC-NEC "
     "and excluded on diagnosis.",
     "CORRECTION: the v2 version of this table carried a footnote stating that N = 16 specimens "
     "'with complete paired clinical and DLL3 quantification metadata were included in the primary "
     "analysis'. That footnote was incorrect and is withdrawn. Exclusion is by adjudicated "
     "diagnosis, not by metadata completeness, and the resulting N is 15.",
     "RNA extraction: RNeasy FFPE Kit (Qiagen); library preparation: TruSeq RNA Exome Kit "
     "(Illumina). Sequencing: NovaSeq 6000 (Illumina), target 25-40 M reads per sample. "
     "Pipeline: CellCarta in-house RNA-seq pipeline; reference GENCODE GRCh38.",
     "Source: CellCarta P1947 study; metadata delivered 2026-06-22."])})

# ------------------------------------------------------------------ S7
counts = pd.read_csv(os.path.join(SRC, "SuppTable_Chan_PerDonor_CellCounts_v3.csv"))
eff = pd.DataFrame([
    ("Primary (>=3 cells/subgroup)", "PD-1", 15, "14/15", 18.5, 0.967, 1.8e-04),
    ("Primary (>=3 cells/subgroup)", "TIM-3", 15, "15/15", 26.0, 1.000, 6.1e-05),
    ("Primary (>=3 cells/subgroup)", "GZMB", 15, "11/15", 12.9, 0.583, 4.8e-02),
    ("Sensitivity (>=1 cell/subgroup)", "PD-1", 18, "16/18", 22.2, 0.974, 4.6e-05),
    ("Sensitivity (>=1 cell/subgroup)", "TIM-3", 18, "16/18", 24.9, 0.953, 5.3e-05),
    ("Sensitivity (>=1 cell/subgroup)", "GZMB", 18, "11/18", 8.3, 0.294, 3.1e-01),
], columns=["Analysis", "Marker", "N evaluable", "Higher in TNFRSF9+", "Median delta (pp)",
            "Rank-biserial r", "Wilcoxon P"])
write("SuppTable_S7_Chan_PerDonor_Counts_and_EffectSizes_v3.xlsx", {
    "S7a_cell_counts": (
        "Supplementary Table S7a. Per-donor cell counts for every quantity entering the Chan "
        "atlas analyses.", counts,
        ["Evaluability under the primary rule: >= 10 CD8+ cells per donor AND >= 3 cells in each "
         "TNFRSF9 subgroup. 15 of 19 donors qualify.",
         "The atlas contains 578 CD8+ Texh cells in total; donor RU1195 contributes 430 of them "
         "(74%). This is why only 8 donors are evaluable for the Texh-stratified analysis.",
         "The PleuralEffusion specimen has no TNFRSF9+ CD8 cells and is not evaluable under any "
         "rule."]),
    "S7c_biospecimens": (
        "Supplementary Table S7c. Biospecimen structure of the Chan atlas and the three donors "
        "who contributed more than one biospecimen.",
        pd.read_csv(os.path.join(SRC, "SuppTable_Chan_MultiBiospecimen_Donors_v3.csv")),
        ["The atlas holds 23 biospecimens from 19 donors (HTAN_Biospecimen_ID and donor_id in the "
         "CELLxGENE object). Sixteen donors contributed one biospecimen each.",
         "RU1108 contributed three lung resection pieces; RU1144 and RU1181 each contributed a "
         "lung and a lymph-node specimen. This table lists those seven biospecimens.",
         "Cells were pooled to the donor for every analysis. Donor-level pseudo-bulk is the unit "
         "that avoids pseudoreplication, and four of these seven biospecimens hold fewer than 20 "
         "CD8+ cells, so a per-biospecimen analysis would lose most of them to the evaluability "
         "rule.",
         "Every biospecimen contributes both epithelial and immune cells, so no quantity is "
         "compared across specimens. Supplementary Table S3d repeats the principal post hoc test "
         "under three alternative units of analysis."]),
    "S7b_effect_sizes": (
        "Supplementary Table S7b. Within-donor paired comparison of marker positivity in TNFRSF9+ "
        "versus TNFRSF9- CD8+ T cells.", eff,
        ["Positivity = detected transcript (normalized expression > 0). CD8+ pool = CD8+ Texh, "
         "Teff and Tmem (author_cell_type, Chan atlas immune object).",
         "Effect size = matched-pairs rank-biserial correlation.",
         "GZMB is threshold-dependent and is reported as exploratory. See Supplementary Table S9 "
         "for depth-adjusted analyses in which GZMB is not significant."])})

# ------------------------------------------------------------------ S1
s7 = pd.DataFrame([
    ("Institutional cohort composition", "15 pathologically adjudicated SCLC",
     "All 17 sequenced specimens (diagnosis-inclusive)",
     "CD274/APM ratio r = +0.579 (P = 0.024)", "r = +0.608 (P = 0.0096)",
     "Pathology adjudication (2026-08-06) reclassified 21029B2866 and 21029T2952 as LC-NEC. v2 "
     "stated 'all 17 histologically confirmed SCLC; no exclusions', which was incorrect. "
     "Excluding LC-NEC also makes the RNA-seq criterion identical to the mIF criterion."),
    ("Jiang population", "79 tumours only", "All 86 including 7 normal lung controls",
     "r = -0.299 (P = 0.0075)", "r = -0.410 (P = 8.7e-05)",
     "Controls are DLL3-null and APM-high, so retaining them inflates the correlation without "
     "tumour-biological basis."),
    ("APM module definition", "5 genes (HLA-A, HLA-B, B2M, TAP1, TAP2), identical in both cohorts",
     "6 genes including HLA-C (Jiang only)", "Jiang r = -0.299", "Jiang r = -0.305",
     "HLA-C has zero variance in all 81 Cologne samples. v2 compared a 6-gene Jiang module with a "
     "5-gene Cologne module, which was not comparable."),
    ("APM score computation", "Mean of gene-wise z-scored log2 expression",
     "Mean of raw log2 expression", "Jiang r = -0.299", "Jiang r = -0.301",
     "Under the harmonized 5-gene definition this choice is immaterial."),
    ("Chan paired-test evaluability",
     ">= 10 CD8 cells/donor AND >= 3 per TNFRSF9 subgroup (15 donors)",
     ">= 1 cell per subgroup (18 donors)", "PD-1 median delta +18.5 pp (P = 1.8e-04)",
     "PD-1 median delta +22.2 pp (P = 4.6e-05)", "Unchanged from v2."),
    ("Chan GZMB (effector arm)", "Reported as EXPLORATORY", "-", "P = 0.048 (>= 3 cells)",
     "P = 0.31 (>= 1 cell); depth-adjusted OR 1.24, P = 0.19",
     "Not significant in any depth-adjusted analysis. The phenotype is NOT described as an "
     "effector-exhaustion dual program."),
    ("Tumour DLL3 metric (single cell)", "Per-donor mean normalized epithelial DLL3",
     "Per-donor % DLL3-positive epithelial cells", "DLL3~Texh% r = +0.505 (P = 0.028)",
     "r = +0.346 (P = 0.15)",
     "The composition association is metric-dependent and is reported as fragile."),
    ("Multivariable estimator (discovery)", "OLS with beta, 95% CI and P reported", "-",
     "beta = +0.297 [-0.621, +1.216], P = 0.48", "-",
     "v2 reported a ridge coefficient (beta = +0.382) with no penalty justification, no interval "
     "and no P value. The ridge model is withdrawn and replaced by a fully reported OLS fit, "
     "which does not support an independent DLL3 contribution."),
    ("Leukocyte-content adjustment", "POST HOC sensitivity analysis", "-",
     "Jiang beta = -0.102 (P = 0.32); Cologne beta = -0.094 (P = 0.21)", "-",
     "Performed after the primary results were assembled. Not pre-specified."),
    ("Equivalence bounds on the central null", "TOST on Fisher z, descriptive", "-",
     "Smallest supported margin |rho| = 0.425 to 0.475", "-",
     "Margins were NOT pre-specified, so these bounds are descriptive, not confirmatory."),
    ("Axis-versus-axis test", "POST HOC", "-",
     "five-gene epithelial APM vs TNFRSF9+ %CD8 r = +0.614 "
     "(P = 0.0052, BH q = 0.0155)",
     "six-gene sensitivity r = +0.612 (P = 0.0053)",
     "Performed in response to peer review, which identified that v2 claimed to have tested the "
     "two axes against one another when in fact only DLL3-versus-axis-2 had been tested."),
    ("mIF FDR family", "All 51 derived phenotypes, one BH pass", "No subset correction",
     "0 of 51 reach q < 0.05 (min q = 0.895)", "-",
     "No individual phenotype is singled out, because all 51 are non-significant."),
], columns=["Analytic choice", "PRIMARY", "Sensitivity", "Value under primary",
            "Value under sensitivity", "Rationale / what changed from v2"])
write("SuppTable_S1_Analysis_Specification_v3.xlsx", {"S1": (
    "Supplementary Table S1. Analysis specification: primary and sensitivity analyses, and every "
    "change from v2.", s7,
    ["Analyses marked POST HOC were performed after the primary results were assembled and are "
     "labelled as post hoc wherever they appear in the manuscript."])})

# ------------------------------------------------------------------ S4
s8 = pd.DataFrame([
    ("Bulk RNA-seq", "Tumours sequenced (CellCarta P1947)", 17, "-",
     "All 17 passed the pre-specified DV200 >= 40 threshold."),
    ("Bulk RNA-seq", "EXCLUDED - diagnosis", 2, "-",
     "21029B2866 and 21029T2952, adjudicated as large-cell neuroendocrine carcinoma (LC-NEC) "
     "rather than SCLC on independent review of the original pathology reports and slides, "
     "performed without reference to any analysis result (2026-08-06)."),
    ("Bulk RNA-seq", "ANALYSED (primary)", 15, "-",
     "Pathologically confirmed SCLC. Continuous-DLL3 analyses use all 15."),
    ("Bulk RNA-seq", "Sensitivity population", 17, "-",
     "The diagnosis-inclusive 17-specimen set is retained as a sensitivity analysis throughout "
     "(Supplementary Table S1)."),
    ("Bulk RNA-seq", "Subset with DLL3 IHC H-score", 10, "-",
     "H-scores 0, 6, 10, 10, 50, 100, 180, 270, 270, 300. These 10 are the mIF analysis cohort."),
    ("mIF", "Archival batch (slides 2023-04)", 15, 372,
     "EXCLUDED IN FULL. No linked DLL3 H-score; early-stage surgical specimens; 6 of the 15 are "
     "LC-NEC rather than SCLC."),
    ("mIF", "Prospective batch (slides 2025-01)", 16, 260, "-"),
    ("mIF", "Prospective, EXCLUDED - diagnosis", 1, "-",
     "21029B2866 (LC-NEC). Reported separately from the H-score exclusions below, because "
     "exclusion by diagnosis and exclusion by missing exposure measurement carry different "
     "selection-bias implications."),
    ("mIF", "Prospective, EXCLUDED - no DLL3 H-score", 5, 67,
     "21029B2861, 21029B2863, 21029T2948, 21029T2949, 21029T2950."),
    ("mIF", "ANALYSED (primary)", 10, 193,
     "21029B2857, B2858, B2859, B2860, B2864, B2865, T2945, T2951, T2953, T2954. ROI count per "
     "patient ranges 1-38 (median 17), which is why patient-level area-weighted means are the "
     "unit of analysis."),
    ("Overlap", "mIF patients also in the RNA-seq primary set", 10, "-",
     "All 10 mIF patients are a strict subset of the 15 analysed RNA-seq patients; the two "
     "institutional analyses are therefore not independent."),
    ("Jiang (external)", "Samples in the deposited GSE60052 matrix", 86, "-", "-"),
    ("Jiang (external)", "Normal lung controls, EXCLUDED", 7, "-",
     "B08-3483NA, B08-3758NA, B08-4386NA, B08-4579NA, N08-3503A, N08-3758A, N08-4497A. Median "
     "log2 DLL3 0.00 vs 4.50 in tumours (P = 0.0007); median APM 8.43 vs 6.87 (P = 1.1e-05)."),
    ("Jiang (external)", "ANALYSED (primary)", 79, "-", "SCLC tumours only."),
    ("George/Cologne (external)", "ANALYSED (primary)", 81, "-",
     "cBioPortal study sclc_ucologne_2015, molecular profile sclc_ucologne_2015_rna_seq_mrna. "
     "No samples shared with the Jiang cohort."),
    ("Chan (external)", "Biospecimens deposited", 23, "-",
     "23 HTAN_Biospecimen_IDs. Three donors contributed more than one (Supplementary Table S7c)."),
    ("Chan (external)", "Distinct donors, the unit of analysis", 19, "-",
     "18 named donors plus one donor labelled only 'PleuralEffusion'. Cells were pooled to the "
     "donor; the manuscript reports donors, not biospecimens."),
    ("Chan (external)", "Evaluable for the primary paired test", 15, "-",
     ">= 10 CD8 cells and >= 3 cells in each TNFRSF9 subgroup."),
    ("Chan (external)", "Evaluable for the sensitivity paired test", 18, "-",
     ">= 1 cell in each subgroup; one donor has no TNFRSF9+ CD8 cells."),
], columns=["Cohort", "Step", "Patients / samples", "ROIs", "Reason / detail"])
write("SuppTable_S4_Cohort_Flow_v3.xlsx", {"S4": (
    "Supplementary Table S4. Cohort flow.", s8,
    ["CHANGED FROM v2: v2 recorded 'Excluded: 0 - none. All 17 are histologically confirmed "
     "SCLC'. Pathology adjudication on 2026-08-06 established that two specimens are LC-NEC. The "
     "primary RNA-seq analysis set is therefore 15, not 17.",
     "This also resolves an internal inconsistency in v2, in which LC-NEC was grounds for "
     "exclusion in the mIF cohort but not in the RNA-seq cohort."])})

# ------------------------------------------------------------------ S6
s9 = pd.read_csv(os.path.join(SRC, "SuppTable_External_Replication_and_Leukocyte_Sensitivity_v3.csv"))
meta9 = pd.read_csv(os.path.join(SRC, "Meta_TwoExternalCohorts_v3.csv"))
write("SuppTable_S6_External_Replication_and_Leukocyte_v3.xlsx", {
    "S6a_models": (
        "Supplementary Table S6a. External replication in two independent bulk cohorts, and "
        "leukocyte-content sensitivity analysis.", s9,
        ["The APM module uses the same five genes (HLA-A, HLA-B, B2M, TAP1, TAP2) in both "
         "cohorts. HLA-C has zero variance across all 81 George/Cologne samples.",
         "Fisher z and percentile bootstrap intervals are both reported here. Figures use Fisher "
         "z intervals only, so that a single interval method appears within any one display.",
         "PTPRC models are POST HOC. PTPRC adjustment markedly attenuates the DLL3 coefficient. "
         "PTPRC may be a confounder, a mediator, or a proxy for tissue composition, and this "
         "design cannot distinguish among them."]),
    "S6b_meta": (
        "Supplementary Table S6b. Fixed-effect meta-analysis of the two external cohorts.", meta9,
        ["With k = 2 cohorts, Cochran's Q has minimal power. Heterogeneity is reported as not "
         "detectable, never as absent."])})

# ------------------------------------------------------------------ S3
s10 = pd.read_csv(os.path.join(SRC, "SuppTable_Axis1_vs_Axis2_and_Equivalence_v3.csv"))
adj = pd.read_csv(os.path.join(SRC, "SuppTable_RNAseq_Diagnosis_Adjudication_v3.csv"))
write("SuppTable_S3_Equivalence_and_AxisVsAxis_v3.xlsx", {
    "S3b_axis_and_equivalence": (
        "Supplementary Table S3b. Equivalence bounds on the central null, and the post hoc "
        "direct test of the two immune features against each other (Chan atlas, 19 donors).",
        s10,
        ["'smallest_equivalence_margin' is the smallest margin delta for which two one-sided "
         "tests both reject at alpha = 0.05, that is, the tightest bound the data support.",
         "Equivalence margins were NOT pre-specified. These bounds are descriptive.",
         "The axis-versus-axis analysis is POST HOC. It was performed because peer review "
         "identified that the v2 Figure 4 legend claimed the two axes had been tested against "
         "each other when in fact only DLL3-versus-axis-2 had been tested.",
         "The primary epithelial APM score uses the same five genes as the harmonized bulk score "
         "(HLA-A, HLA-B, B2M, TAP1, TAP2). Adding HLA-C as a six-gene sensitivity analysis "
         "produces materially unchanged estimates.",
         "The Benjamini-Hochberg family contains the three prespecified abundance metrics in this "
         "post hoc direct-test family. TNFRSF9+ as %CD8 is the principal exploratory endpoint; "
         "the nested PD-1+ and TIM-3+ fractions are supportive endpoints.",
         "The table reports Fisher z and percentile-bootstrap intervals, leave-one-out ranges, "
         "and sensitivity analyses restricted to donors with at least 20 or 50 CD8 cells.",
         "Epithelial APM is a different estimand from the bulk APM score that carries the "
         "external replication. Within the Chan atlas DLL3 itself is unrelated to the primary "
         "five-gene epithelial APM score (r = +0.026, P = 0.91)."]),
    "S3d_unit_of_analysis": (
        "Supplementary Table S3d. The principal post hoc test and the central null under three "
        "alternative units of analysis.",
        pd.read_csv(os.path.join(SRC, "SuppTable_Chan_UnitOfAnalysis_Sensitivity_v3.csv")),
        ["The published unit is the donor (19). The atlas holds 23 biospecimens, because three "
         "donors contributed more than one (Supplementary Table S7c).",
         "The per-biospecimen rows are descriptive only: the three multi-biospecimen donors "
         "appear more than once, so the 23 rows are not independent.",
         "Both quantities are recomputed from the atlas at each unit with the definitions used "
         "throughout: epithelial APM is the mean of gene-wise z-scored per-unit mean expression "
         "of the five harmonised genes, and the outcome is the percentage of CD8+ cells with "
         "detected TNFRSF9.",
         "The direction and significance of the principal association, and the central null, "
         "hold under every unit."]),
    "S3c_confounder_sensitivity": (
        "Supplementary Table S3c. Shared-technical-driver sensitivity for the post hoc "
        "epithelial APM versus TNFRSF9-defined CD8 abundance association (Chan atlas, 19 donors).",
        pd.read_csv(os.path.join(SRC, "SuppTable_Axis_ConfounderSensitivity_v3.csv")),
        ["Both correlated quantities are derived from the same single-cell object, so a shared "
         "technical factor could in principle generate the association. Each row is a rank-based "
         "partial correlation: both variables and the covariates are rank-transformed, each "
         "variable is regressed on the covariates by OLS, and the Spearman correlation of the "
         "residuals is reported.",
         "The principal single-gated endpoint remains significant under every adjustment model "
         "(r = +0.493 to +0.714). The doubly gated endpoints do not, which is consistent with "
         "their designation as supportive.",
         "Sheet S10d lists the association of each candidate driver with the two correlated "
         "quantities; none is significantly associated with both."]),
    "S3e_candidate_drivers": (
        "Supplementary Table S3e. Association of each candidate shared technical driver with "
        "epithelial APM score and with TNFRSF9-positive abundance.",
        pd.read_csv(os.path.join(SRC, "SuppTable_Axis_CandidateDrivers_v3.csv")),
        ["No candidate driver is significantly associated with both quantities "
         "(all |r| <= 0.38, P >= 0.11).",
         "Sequencing depth is the strongest single candidate and is the covariate that attenuates "
         "the association most (from r = +0.614 to r = +0.509)."]),
    "S3a_adjudication": (
        "Supplementary Table S3a. Per-specimen pathology adjudication of the 17 sequenced "
        "specimens.", adj,
        ["Adjudication was performed by the institutional pathology service against the original "
         "pathology reports and slides, without reference to any analysis result."])})

# ------------------------------------------------------------------ S9
s11 = pd.read_csv(os.path.join(SRC, "SuppTable_DetectionDepth_Sensitivity_v3.csv"))
write("SuppTable_S9_DetectionDepth_Sensitivity_v3.xlsx", {"S9": (
    "Supplementary Table S9. Sequencing-depth sensitivity analyses for the TNFRSF9-defined CD8 "
    "phenotype (Chan atlas).", s11,
    ["TNFRSF9 positivity and marker positivity are both defined by transcript detection in the "
     "same cells, so global detection depth could in principle generate co-detection. These "
     "analyses test that possibility directly.",
     "Depth does contribute: TNFRSF9+ CD8 cells have modestly higher total UMI (AUC = 0.599) and "
     "more detected genes (AUC = 0.606).",
     "Depth does not account for the finding: PD-1 and TIM-3 remain significant after adjustment "
     "for log10 total counts, CD8 subtype and donor; on continuous expression; and under "
     "within-donor total-count-decile matching.",
     "GZMB is not significant in any of these analyses, which reinforces its exploratory status."])})

# ------------------------------------------------------------------ S10
s13 = pd.read_csv(os.path.join(SRC, "SuppTable_CCLE_LeukocyteFree_Falsification_v3.csv"))
write("SuppTable_S10_CCLE_LeukocyteFree_Falsification_v3.xlsx", {"S10": (
    "Supplementary Table S10. Leukocyte-free falsification test of the bulk DLL3-APM "
    "association in 59 DepMap/CCLE small cell lung cancer cell lines.",
    s13,
    ["Rationale: if the bulk DLL3-APM inverse association is compositional, it must disappear in "
     "a system containing tumour cells and no leukocytes; if it is tumour-cell-intrinsic, it must "
     "persist. Cell lines provide that system.",
     "PRIMARY TEST: DLL3 versus the same harmonized five-gene APM module is null "
     "(r = -0.053, P = 0.69), against r = -0.299 (Jiang) and r = -0.411 (George/Cologne) in "
     "tissue.",
     "POSITIVE CONTROL: tumour-intrinsic DLL3 biology is preserved (DLL3-ASCL1 r = +0.667, "
     "P = 7.9e-09), so the null is not an artefact of the cell-line system failing to express "
     "DLL3 biology.",
     "PURITY CONTROLS: PTPRC median expression is 0.098 and PTPRC no longer tracks APM "
     "(r = +0.048, P = 0.72), against r = +0.599 and +0.814 in the two tissue cohorts. The "
     "compositional mechanism is therefore genuinely absent from this system.",
     "LIMITATION: cell lines lack the in vivo interferon and Notch milieu, so a reviewer may "
     "argue the null reflects loss of microenvironmental signalling rather than absence of a "
     "tumour-intrinsic link. This test is presented as one of three convergent estimands "
     "(PTPRC adjustment in tissue; epithelial single-cell; leukocyte-free cell lines), not as "
     "proof on its own."])})

# ------------------------------------------------------------------ S13
s12 = pd.DataFrame([
    ("s01_axis1_vs_axis2_direct_test.py",
     "Chan per-donor epithelial APM; axis-versus-axis direct test; TOST equivalence bounds", "S3"),
    ("s02_external_replication_harmonised.py",
     "Jiang + George/Cologne harmonized 5-gene APM; stepwise NE and PTPRC adjustment; "
     "fixed-effect meta-analysis", "S6, Figure 2A-C"),
    ("s03_detection_depth_sensitivity.py",
     "Sequencing-depth exposure check; cell-level logistic models; continuous-expression and "
     "depth-matched paired tests; per-donor cell counts", "S7, S9"),
    ("s04_discovery_cohort_N15.py",
     "Institutional discovery cohort at N = 15 adjudicated SCLC, with the N = 17 "
     "diagnosis-inclusive sensitivity analysis", "S1, S5, Figure 2E-G"),
    ("s05_generate_supplementary_tables_v3.py",
     "This script; generates Supplementary Tables S1-S13 from the derived outputs", "S1-S13"),
    ("s06_generate_main_figures_v3.py", "Main Figures 1, 2 and 5 from the derived outputs",
     "Figures 1, 2, 5"),
    ("s07_generate_figure3_and_supplementary_figures_v3.py",
     "Main Figures 3 and 4 and Supplementary Figures S1-S10; per-donor paired plot data",
     "Figures 3, 4, Supplementary Figures S1-S10"),
    ("s08_generate_tables_and_docx_v3.py", "Main Table 1 (Word), Supplementary Table S11 and the manuscript DOCX",
     "Table 1, S11"),
    ("s09_verify_manuscript_numbers_v3.py",
     "Verification: every quoted statistic, panel reference, citation and format rule is "
     "re-checked against the derived outputs; fails on any mismatch", "none (verification)"),
    ("s10_axis_confounder_sensitivity.py",
     "Shared-technical-factor sensitivity for the epithelial APM versus TNFRSF9+ CD8 test "
     "(rank-based partial correlations; candidate-driver associations)", "S3c, S3e"),
    ("s14_specimen_unit_sensitivity.py",
     "Unit-of-analysis sensitivity: the principal post hoc test and the central null computed per "
     "biospecimen, in single-biospecimen donors, and with one biospecimen per donor", "S7c, S3d"),
    ("s13_graphical_abstract.py", "Graphical Abstract", "none (submission item)"),
    ("s11_leukocyte_free_falsification.py",
     "DepMap/CCLE SCLC cell lines: DLL3 versus the five-gene APM module in a leukocyte-free "
     "system, with DLL3-ASCL1 positive control and PTPRC purity checks", "S10"),
    ("20_Master_Regenerate_All_Statistics.py",
     "v2 master statistics (Jiang individual genes, Chan paired tests, mIF), retained unchanged "
     "where the pathology adjudication does not apply", "S7b, S11, S12"),
], columns=["Script", "What it computes", "Supplementary items produced"])
write("SuppTable_S13_Script_Manifest_v3.xlsx", {"S13": (
    "Supplementary Table S13. Analysis script manifest and reproducibility.", s12,
    ["Environment: Python 3.9.2; numpy 1.26.4; pandas 2.3.3; scipy 1.13.1; statsmodels 0.14.6; "
     "anndata 0.10.9; matplotlib 3.9.4; scanpy; openpyxl. Wilcoxon signed-rank P values depend on "
     "the SciPy version's default method (SciPy 1.17 gives P = 0.0094 for the Figure 3D lineage "
     "test); every reported value is from SciPy 1.13.1.",
     "Random seed 20260806 for all bootstrap and matching procedures in v3 scripts; seed "
     "20260725 for retained v2 outputs.",
     "Scripts and derived tables will be deposited in a public versioned archive on acceptance."])})

print("\nAll supplementary tables written to:", OUT)
