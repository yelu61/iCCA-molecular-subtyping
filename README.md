# Multi-Cohort Discovery and Validation of Intrahepatic Cholangiocarcinoma Molecular Subtypes

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![R: >=4.3](https://img.shields.io/badge/R-%3E%3D4.3-brightgreen.svg)](https://www.r-project.org/)
[![Python: 3.10](https://img.shields.io/badge/Python-3.10-blue.svg)](https://www.python.org/)

![Project Overview](assets/project_overview.png)

This repository contains the open-source code, analysis pipelines, staged data, and figure generation workflows for the multi-omics discovery, biological characterization, and clinical validation of **intrahepatic cholangiocarcinoma (iCCA) molecular subtypes**.

---

## Scientific Summary & Locked Evidence

- **Consensus Discovery (k = 3)**: Unsupervised multi-cohort consensus clustering across 586 RNA-seq samples identified **497 principal-state tumors** and 89 intermediate profiles, establishing three robust molecular subtypes:
  - **C1 (Metabolic)**: Characterized by preserved hepatic/metabolic programs, xenobiotic metabolism, organic acid transport, high *IDH1/BAP1* mutation rates, and favorable prognosis.
  - **C2 (Immune/Ductal)**: Characterized by stromal remodeling, inflammatory signaling, complement activation, and small/large duct intermediate features.
  - **C3 (Proliferative)**: Characterized by cell-cycle acceleration, DNA replication, *TP53/KRAS* alterations, high early mortality, and predicted sensitivity to MEK/ERK inhibitors.
- **Locked Primary Classifier**: **114-Gene Consensus Linear SVM** (`candidate_models.rds`), selected via three-way consensus (limma-voom, Random Forest, multinomial LASSO) strictly within training folds:
  - **Repeated Nested CV Balanced Accuracy**: **98.79%** (95% CI: 97.78% – 99.64%)
  - **Macro-F1**: **0.9881** (95% CI: 0.9785 – 0.9964)
  - **Multiclass Brier Score**: **0.0272** (95% CI: 0.0186 – 0.0377)
  - **Class AUROC**: C1 = 0.9996, C2 = 0.9996, C3 = 0.9997
  - **Leave-One-Study-Out (LOSO) Transfer**: Cross-study balanced accuracy up to 100% across independent discovery cohorts.
- **Exploratory Minimal Panel**: **17-gene panel** (`compact_candidate.rds`) achieving 97.63% Balanced Accuracy (1.156% drop vs. 114-gene primary).
- **Single-Cell Validation**: Evaluated across **34,227 cells (17 patients)** from GSE151530 and GSE138709 using native UCell signature scores and matched tumor vs. adjacent tissue comparisons without premature malignancy assumptions.
- **Prognostic Impact**: Study-stratified survival analysis demonstrates that C3 is associated with significantly elevated early hazard (time-varying Cox HR > 3.0), providing substantial prognostic discrimination beyond TNM stage and pathology duct type.

---

## Repository Structure

```
├── README.md                      # Comprehensive user & reproducibility guide
├── LICENSE                        # MIT License
├── CITATION.cff.template          # Software citation metadata template
├── environment/                   # Reproducibility configurations
│   ├── environment.yml            # Conda environment definition
│   ├── requirements.txt           # Python package requirements
│   ├── install_packages.R         # R package installer
│   └── sessionInfo_R45.txt        # Exact R 4.5 session details
├── data/                          # Governed metadata and dictionaries
│   ├── subtype_master_table.csv   # Canonical 497-sample patient metadata & subtype calls
│   ├── study_color_dictionary.csv # Cohort palette standards
│   └── subtype_label_color_dictionary.csv # Subtype palette standards
├── pipelines/                     # Interactive Jupyter notebooks for upstream analyses (Steps 0–8)
│   ├── 0-data-qc.ipynb            # Multi-cohort data ingestion, QC & batch evaluation
│   ├── 1-clustering.ipynb         # Consensus clustering & subtype discovery
│   ├── 2-transcriptomic-characterization.ipynb # GSEA, GSVA, WGCNA & regulons
│   ├── 3-duct-type-analysis.ipynb # Observed pathology duct type associations
│   ├── 4-multiomic.ipynb          # Mutations, CPTAC proteomics & microenvironment
│   ├── 5-classifier-validation.ipynb # 114-gene consensus SVM, nested CV & LOSO
│   ├── 6-clinical-survival.ipynb  # Study-stratified Cox, non-proportional hazards
│   ├── 7-single-cell-validation.ipynb # 34,227-cell atlas & native UCell scoring
│   └── 8-visualization.ipynb      # Interactive generation & inspection of all 6 figures
├── figures/                       # Standalone figure generation scripts & staged data
│   ├── run_all_figures.R          # Master script to reproduce all 6 figures in sequence
│   ├── plot_figure1.R             # Figure 1: Consensus clustering & discovery
│   ├── plot_figure2.R             # Figure 2: Biological programs & signatures
│   ├── plot_figure3.R             # Figure 3: Single-cell transcriptomic validation
│   ├── plot_figure4.R             # Figure 4: Multi-omics landscape & proteomics
│   ├── plot_figure5.R             # Figure 5: Machine learning classifier validation
│   ├── plot_figure6.R             # Figure 6: Clinical prognosis & therapy hypotheses
│   ├── data/                      # Self-contained staged inputs for all figures (~19MB)
│   │   ├── figure1/
│   │   ├── figure2/
│   │   ├── figure3/
│   │   ├── figure4/
│   │   ├── figure5/
│   │   └── figure6/
│   └── outputs/                   # Generated individual panels and full composite figures
└── docs/                          # Detailed methodology notes & data dictionary
    ├── DATA_DICTIONARY.md
    └── final_panel_manifest.csv
```

---

## Quick Start — Reproducing Publication Figures

All 6 primary publication figures can be generated in under 1 minute from the staged data included in this repository.

### 1. Environment Setup

#### Option A: Using Conda (Recommended)
```bash
git clone https://github.com/your-username/iCCA-molecular-subtyping.git
cd iCCA-molecular-subtyping

conda env create -f environment/environment.yml
conda activate icca-subtyping
```

#### Option B: Using Existing R Installation (R >= 4.3)
```bash
Rscript environment/install_packages.R
```

### 2. Generate All Figures

Run the master figure pipeline:
```bash
Rscript figures/run_all_figures.R
```

This sequentially generates:
- **`figures/outputs/figure1/Figure_1_composite.pdf` & `.png`** (Multi-omics discovery)
- **`figures/outputs/figure2/Figure_2_composite.pdf` & `.png`** (Biological characterization)
- **`figures/outputs/figure3/Figure_3_composite.pdf` & `.png`** (Single-cell validation)
- **`figures/outputs/figure4/Figure_4_composite.pdf` & `.png`** (Multi-omics landscape)
- **`figures/outputs/figure5/Figure_5_composite.pdf` & `.png`** (114-gene classifier validation)
- **`figures/outputs/figure6/Figure_6_composite.pdf` & `.png`** (Survival & therapeutic hypotheses)

All composite figures are strictly rendered to publication standards (180 mm × 270 mm, 300/600 dpi).

### 3. Generate Specific Figures Individually

Each script can also be executed independently:
```bash
Rscript figures/plot_figure1.R
Rscript figures/plot_figure2.R
Rscript figures/plot_figure3.R
Rscript figures/plot_figure4.R
Rscript figures/plot_figure5.R
Rscript figures/plot_figure6.R
```

---

## Upstream Pipeline Walkthrough

The `pipelines/` directory contains complete interactive Jupyter notebooks reflecting the step-by-step analytical sequence:

| Notebook | Focus | Primary Outputs |
|---|---|---|
| [`0-data-qc.ipynb`](pipelines/0-data-qc.ipynb) | Multi-cohort ingestion, log2-transformation, MAD filtering, batch QC | Cross-study expression matrices |
| [`1-clustering.ipynb`](pipelines/1-clustering.ipynb) | ConsensusClusterPlus, NMF, k-means, silhouette widths (k=2..6) | Subtype assignments (C1, C2, C3) |
| [`2-transcriptomic-characterization.ipynb`](pipelines/2-transcriptomic-characterization.ipynb) | Corrected `subtype - others` DEGs, GSEA, GSVA, WGCNA, DoRothEA TFs | Subtype functional signatures |
| [`3-duct-type-analysis.ipynb`](pipelines/3-duct-type-analysis.ipynb) | Association with observed pathology duct types | Small-duct vs. large-duct distribution |
| [`4-multiomic.ipynb`](pipelines/4-multiomic.ipynb) | Somatic mutation oncoplots, CPTAC proteomic validation, xCell immune deconvolution | Multi-omics convergence |
| [`5-classifier-validation.ipynb`](pipelines/5-classifier-validation.ipynb) | 114-gene consensus SVM, 5x5 nested CV, LOSO cross-study transfer, 17-gene minimal panel | Model objects & validation tables |
| [`6-clinical-survival.ipynb`](pipelines/6-clinical-survival.ipynb) | Study-stratified multivariable Cox models, non-proportional hazards, time-varying effects | Clinical hazard ratios & KM curves |
| [`7-single-cell-validation.ipynb`](pipelines/7-single-cell-validation.ipynb) | 34,227-cell atlas (GSE151530 + GSE138709), native UCell signature scoring, source holds | Single-cell program localization |
| [`8-visualization.ipynb`](pipelines/8-visualization.ipynb) | Interactive runner and inline viewer for Figures 1 through 6 | Composite figure validation |

---

## Data Governance & Master Metadata

The canonical sample classifications and patient covariates are stored in [`data/subtype_master_table.csv`](data/subtype_master_table.csv).
All downstream analyses reference this table directly. Key variables:
- `Sample`: Unique tumor identifier.
- `Study`: Cohort source (2019Ahn, 2022Dong, 2023Deng, 2023Jeon, 2023Lee, 2024Fan, TCGA-CHOL).
- `FinalSubtype`: Subtype classification (`C1`, `C2`, `C3`, or intermediate `C1_C2`, `C1_C3`, `C2_C3`).
- `SubtypeClass`: `Main` (principal state, n=497) or `Intermediate` (n=89).
- `OS_months` & `OS_status`: Overall survival follow-up and event indicator.
- `DuctType_observed`: Pathologically observed duct type (small duct / large duct).

---

## Citation & License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

If you find this work or code useful, please cite our manuscript:
```bibtex
@article{iCCA_molecular_subtypes_2026,
  title={Multi-Cohort Discovery and Validation of Intrahepatic Cholangiocarcinoma Molecular Subtypes},
  author={Research Consortium},
  journal={Pending Publication},
  year={2026}
}
```
