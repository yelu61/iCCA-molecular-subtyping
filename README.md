# Molecular Subtyping of Intrahepatic Cholangiocarcinoma (iCCA)

Multi-omics integrative analysis identifying and validating three molecular subtypes (C1, C2, C3) of intrahepatic cholangiocarcinoma through consensus clustering.

## Overview

This project integrates transcriptomic, proteomic, and genomic data from multiple independent studies (2012–2025) together with TCGA-CHOL to establish a robust molecular classification system for iCCA. Three subtypes were identified:

| Subtype | Characteristics |
|---------|-----------------|
| **C1** | Metabolic subtype |
| **C2** | Immune subtype |
| **C3** | Proliferative subtype |

## Analysis Pipeline

The analysis follows a strict sequential workflow. Notebooks must be run in order because downstream steps load `.Rdata`/`.rds` files produced by upstream steps.

```
0-Data/ (raw & processed data)
    ↓
0-data.ipynb — Data preparation and integration
    ↓
1-clustering.ipynb — Consensus clustering (C1/C2/C3 subtypes)
    ↓
2-classfiler.ipynb — Classifier building (SVM / Random Forest)
    ↓
3-multiomic.ipynb — Proteomic validation & multi-omics integration
    ↓
4-Validation.ipynb — External validation (microarray, CCLE, DepMap, prognostic model)
```

### Step 0: Data Preparation
- Filters genes expressed in ≥20% of samples
- Handles batch effects (FPKM→TPM conversion, log2 transformation)
- Produces: `0-Data-rnaseq.Rdata`, microarray validation datasets

### Step 1: Clustering
- Batch correction (ComBat)
- ConsensusClusterPlus + NMF clustering
- Differential expression per subtype
- Tumor microenvironment analysis (ESTIMATE/EPIC)
- Drug sensitivity prediction (oncoPredict)
- **Key output**: `1-Clustering/data/sample_info_with_clusters.xlsx`

### Step 2: Classifier
- SVM / Random Forest feature selection
- LASSO for large-duct vs small-duct classification
- **Key output**: `2-Classfiler/data/pred.csv`

### Step 3: Multi-omics Integration
- Proteomic validation
- Pathway enrichment comparison (proteomics vs transcriptomics)
- mRNA + protein + CNV integration

### Step 4: Validation
- Microarray data validation (Illumina & Affymetrix)
- CCLE drug sensitivity with predicted subtypes
- DepMap expression data validation
- LASSO-Cox prognostic risk score model

## Repository Structure

```
.
├── 0-data.ipynb                  # Data preparation
├── 1-clustering.ipynb            # Consensus clustering
├── 2-classfiler.ipynb            # Classifier building
├── 3-multiomic.ipynb             # Proteomic & multi-omics
├── 4-Validation.ipynb            # External validation
├── CLAUDE.md                     # Project documentation
├── iCCA_Subtype_Analysis_Workflow.md  # Detailed workflow (Chinese)
├── sample_info_with_predictions.csv   # Final sample subtype predictions
├── sample_info_with_subtype.xlsx      # Final sample subtype assignments
├── cluster_compare_multiomic.csv      # Multi-omics comparison results
├── 0-Data/                       # Raw & processed data (not tracked)
├── 1-Clustering/                 # Clustering outputs (not tracked)
├── 2-Classfiler/                 # Classifier outputs (not tracked)
├── 3-Proteomic/                  # Proteomic results (not tracked)
├── 4-Validation/                 # Validation outputs (not tracked)
├── Meeting/                      # Presentation materials (not tracked)
└── sources/                      # Reference materials
```

## Key Output Files

| File | Description |
|------|-------------|
| `sample_info_with_subtype.xlsx` | Final subtype assignments for all samples |
| `sample_info_with_predictions.csv` | Classification predictions |
| `cluster_compare_multiomic.csv` | Multi-omics subtype comparison |
| `single_gene_cox_results.csv` | Single-gene Cox regression results |

## Dependencies

The analysis is implemented in **R** via Jupyter notebooks. Key packages:

```r
# Core analysis
ConsensusClusterPlus    # Consensus clustering
sva                     # Batch correction (ComBat)
limma                   # Differential expression
NMF                     # NMF clustering
WGCNA                   # Weighted gene co-expression
GSVA                    # Gene set variation analysis
estimate                # TME scores

# Visualization
pheatmap
ComplexHeatmap
ggplot2
EnhancedVolcano
survminer

# Machine learning
e1071                   # SVM
randomForest
glmnet                  # LASSO / LASSO-Cox
caret

# Data manipulation
tidyverse
data.table
plyr
```

## Data Sources

The project integrates data from the following studies:

| Study | Year | Platform |
|-------|------|----------|
| Andersen et al. | 2012 | Microarray |
| Sia et al. | 2013 | Microarray |
| Zou et al. | 2014 | Microarray |
| Chaisaingmongkol et al. | 2017 | RNA-seq |
| Jusakual et al. | 2017 | RNA-seq |
| Maeve et al. | 2018 | RNA-seq |
| Wardell et al. | 2018 | RNA-seq |
| Ahn et al. | 2019 | RNA-seq |
| Job et al. | 2020 | RNA-seq |
| Zhang et al. | 2020 | RNA-seq |
| Boerner et al. | 2021 | RNA-seq |
| Bao, Dong, Gao, Lin, Ma, Song, Wang et al. | 2022 | RNA-seq |
| Cho, Deng, Jeon, Lee et al. | 2023 | RNA-seq |
| Chen, Fan, Liu, Zhou et al. | 2024 | RNA-seq |
| Ruan et al. | 2025 | RNA-seq |
| TCGA-CHOL | — | RNA-seq |

## Reproducibility Notes

- **Random seed**: `123` is the dominant seed for clustering; `42` for classifier training; `61` appears once in a single cell.
- **Expression format**: Matrices are **genes × samples** (rows = genes, columns = samples).
- **Working directory**: Notebooks assume the project root as working directory.
- **Execution time**: `1-clustering.ipynb` contains long-running cells (ConsensusClusterPlus with 100 reps, NMF with 20 runs) that may take tens of minutes to hours.

## Citation

If you use this code or data, please cite the corresponding publication (TBD).

## License

TBD — to be determined upon publication.
