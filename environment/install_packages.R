# ==============================================================================
# R Package Installation Script for iCCA Molecular Subtyping
# Compatible with R >= 4.3 (tested on R 4.5)
# ==============================================================================

options(repos = c(CRAN = "https://cloud.r-project.org/"))

cran_packages <- c(
  "ragg", "circlize", "png", "data.table", "ggplot2", "cowplot",
  "survival", "survminer", "readr", "dplyr", "pROC", "patchwork",
  "scales", "jsonlite", "e1071", "randomForest", "glmnet"
)

new_cran <- cran_packages[!(cran_packages %in% installed.packages()[, "Package"])]
if (length(new_cran) > 0) {
  cat("Installing CRAN packages:", paste(new_cran, collapse = ", "), "\n")
  install.packages(new_cran)
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_packages <- c("ComplexHeatmap")
new_bioc <- bioc_packages[!(bioc_packages %in% installed.packages()[, "Package"])]
if (length(new_bioc) > 0) {
  cat("Installing Bioconductor packages:", paste(new_bioc, collapse = ", "), "\n")
  BiocManager::install(new_bioc)
}

cat("All required R packages are installed.\n")
