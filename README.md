# Analysis code for the DLL3 / SCLC immune-microenvironment study

Code and derived data for:

> Kobayashi N, Somekawa K, Kajita Y, Muraoka S, Kaneko A, Kubo S, Hirose T, Teranishi S, Irie H,
> Fujii E, Yamazaki M, Kamata-Sakurai M, Kashima K. DLL3 expression and its associations with
> antigen-presentation genes and *TNFRSF9*-expressing CD8 T cells in small cell lung cancer.
> Submitted.

Every statistic, figure and table in the manuscript is produced by the scripts in `scripts/` from
the files in `derived_outputs/`. Nothing in the manuscript is typed in by hand: `s09` re-reads the
derived files, extracts the corresponding claim from the manuscript text, and fails on any
mismatch.

## Layout

| Path | Contents |
|---|---|
| `scripts/` | Analysis, figure, table, verification and packaging scripts (`s01`–`s14`) |
| `derived_outputs/` | The derived tables every reported number is computed from |
| `supplementary_tables/` | Supplementary Tables S1–S13 as published |

## Running the scripts

The reported values were produced with **Python 3.9.2**: numpy 1.26.4, pandas 2.3.3, scipy 1.13.1,
statsmodels 0.14.6, anndata 0.10.9, matplotlib 3.9.4, scanpy, openpyxl, python-docx, reportlab.

Use that environment. Wilcoxon signed-rank *P* values depend on the SciPy version's default method:
SciPy 1.17 returns *P* = 0.0094 for the Figure 3D lineage test where SciPy 1.13.1 returns the
reported *P* = 0.0082. All bootstrap and matching procedures use seed 20260806 (20260725 for the
retained v2 outputs), so runs are reproducible within a version.

```bash
python scripts/s09_verify_manuscript_numbers_v3.py   # verifies the manuscript against the data
python scripts/s06_generate_main_figures_v3.py       # Figures 1, 2, 5
python scripts/s13_graphical_abstract.py             # Graphical Abstract
```

`s01`–`s04`, `s07`, `s10`, `s11` and `s14` recompute the derived tables from the primary data rather than
reading `derived_outputs/`, so they additionally need the source files listed below. The scripts
resolve paths relative to their own location and expect the original project layout
(`../05_Source_Data/...`); adjust `BASE` if you place them elsewhere.

## Data sources

**Public, not redistributed here.**

| Dataset | Location |
|---|---|
| Jiang cohort (86 samples; 79 tumours used) | Gene Expression Omnibus, GSE60052 |
| George/Cologne cohort (81 tumours) | cBioPortal, study `sclc_ucologne_2015`, profile `sclc_ucologne_2015_rna_seq_mrna` |
| Chan SCLC atlas (54,313 epithelial, 16,475 immune cells) | CZ CELLxGENE Discover |
| SCLC cell lines (59) | DepMap / CCLE, OncotreeCode SCLC |
| Supportive concordance check (7 patients) | Gene Expression Omnibus, GSE319155 |

**Institutional.** Bulk RNA-sequencing of 15 pathologically adjudicated SCLC specimens (17
sequenced) is available from the corresponding author on reasonable request, as stated in the
manuscript's Data Availability Statement; the informed consent under which the specimens were
collected does not cover deposition in public repositories. The multiplex-immunofluorescence measurements are included here at both the
region-of-interest and the patient level (`Derived_mIF_ROI_All_51Markers_Results.csv`,
`Derived_mIF_Patient_AreaWeighted_Means.csv`); the primary images are held by the service provider.

No clinical annotation is included anywhere in this repository, by design: no analysis is
stratified by, adjusted for, or interpreted against a clinical variable. Specimen identifiers are
the assay identifiers already published in Supplementary Table S2.

## What each script does

| Script | Produces |
|---|---|
| `s01_axis1_vs_axis2_direct_test.py` | Per-donor epithelial APM score; the feature-1-versus-feature-2 test; TOST equivalence bounds (S3) |
| `s02_external_replication_harmonised.py` | Harmonised five-gene APM in both external cohorts; stepwise NE and *PTPRC* adjustment; fixed-effect meta-analysis (S6, Figure 2A–C) |
| `s03_detection_depth_sensitivity.py` | Sequencing-depth exposure check; cell-level logistic models; continuous-expression and depth-matched paired tests (S7, S9) |
| `s04_discovery_cohort_N15.py` | Institutional discovery cohort at *N* = 15, with the *N* = 17 sensitivity analysis (S1, S5, Figure 2E–G) |
| `s05_generate_supplementary_tables_v3.py` | Supplementary Tables S1–S13 except S11 |
| `s06_generate_main_figures_v3.py` | Figures 1, 2 and 5 |
| `s07_generate_figure3_and_supplementary_figures_v3.py` | Figures 3 and 4 and Supplementary Figures S1–S10 |
| `s08_generate_tables_and_docx_v3.py` | Main Table 1, Supplementary Table S11 and the manuscript DOCX |
| `s09_verify_manuscript_numbers_v3.py` | Verification: every quoted statistic, panel reference, citation and journal format rule |
| `s10_axis_confounder_sensitivity.py` | Shared-technical-factor sensitivity for the feature-versus-feature test (S3c, S3e) |
| `s11_leukocyte_free_falsification.py` | DepMap/CCLE leukocyte-free falsification test with its positive and purity controls (S10) |
| `s12_build_submission_package.py` | Assembles the submission package |
| `s14_specimen_unit_sensitivity.py` | Unit-of-analysis sensitivity: per biospecimen, single-biospecimen donors, one biospecimen per donor (S3d, S7c) |
| `s13_graphical_abstract.py` | Graphical Abstract |
| `20_Master_Regenerate_All_Statistics.py` | Retained v2 master statistics, unchanged where the pathology adjudication does not apply (S7b, S11, S12) |
| `md2docx_v3.py` | Markdown-to-DOCX converter used by `s08` and `s12` |

## Use of AI-assisted technologies

Claude Code (Anthropic) was used to write and refine these scripts and to draft and edit the
manuscript text, as disclosed in the Acknowledgements of the manuscript. The authors read all code,
re-executed it against the source data and verified every reported statistic.

## Contact

Nobuaki Kobayashi, Department of Pulmonology, Yokohama City University Graduate School of Medicine —
nkobayas@yokohama-cu.ac.jp
