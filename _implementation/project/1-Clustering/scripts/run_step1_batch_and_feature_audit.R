#!/usr/bin/env Rscript

# Reproduce and diagnose the historical study correction without overwriting
# either locked expression matrix. Publication-facing tables and figures only.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(limma)
  library(readxl)
})

set.seed(123)
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
data_dir <- file.path(root, "1-Clustering", "data")
table_dir <- file.path(root, "1-Clustering", "publication", "tables")
figure_dir <- file.path(root, "1-Clustering", "publication", "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

load_object <- function(path, object_name) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  stopifnot(exists(object_name, envir = env, inherits = FALSE))
  get(object_name, envir = env, inherits = FALSE)
}

adjusted_rand_index <- function(x, y) {
  tab <- table(x, y)
  choose2 <- function(z) z * (z - 1) / 2
  n <- sum(tab)
  index <- sum(choose2(tab))
  a <- sum(choose2(rowSums(tab)))
  b <- sum(choose2(colSums(tab)))
  expected <- a * b / choose2(n)
  maximum <- (a + b) / 2
  if (maximum == expected) return(1)
  (index - expected) / (maximum - expected)
}

cramers_v <- function(x, y) {
  tab <- table(x, y)
  test <- suppressWarnings(chisq.test(tab, correct = FALSE))
  sqrt(unname(test$statistic) / (sum(tab) * min(nrow(tab) - 1, ncol(tab) - 1)))
}

before <- as.matrix(load_object(file.path(data_dir, "combined_expr_before_correction.rds"), "combined_expr"))
after <- as.matrix(load_object(file.path(data_dir, "combined_expr_after_correction.rds"), "corrected_expr_limma"))
sample_info <- as.data.table(read_excel(file.path(data_dir, "sample_info_with_clusters.xlsx")))
sample_info[, Sample := as.character(Sample)]
sample_info <- sample_info[match(colnames(before), Sample)]
stopifnot(
  identical(dim(before), c(10531L, 586L)),
  identical(dim(after), dim(before)),
  identical(rownames(before), rownames(after)),
  identical(colnames(before), colnames(after)),
  !anyNA(sample_info$Study)
)

# Exact reconstruction is held in memory only and compared with the locked file.
rederived_after <- removeBatchEffect(before, batch = factor(sample_info$Study))
reproduction_difference <- rederived_after - after

within_study_residuals <- function(matrix, study) {
  result <- matrix
  for (current_study in unique(study)) {
    idx <- which(study == current_study)
    result[, idx] <- sweep(matrix[, idx, drop = FALSE], 1L, rowMeans(matrix[, idx, drop = FALSE]), "-")
  }
  result
}
before_residual <- within_study_residuals(before, sample_info$Study)
after_residual <- within_study_residuals(after, sample_info$Study)

correction_contract <- data.table(
  Field = c(
    "Input scale", "Input dimensions", "Batch variable", "Correction function",
    "Biological design", "Locked output policy", "Maximum absolute reconstruction difference",
    "Reconstruction RMSE", "Reconstruction correlation", "Within-study residual correlation"
  ),
  Value = c(
    "Per-study processed log2-scale expression", "10,531 genes x 586 samples", "Study (7 discovery cohorts)",
    "limma::removeBatchEffect", "No biological design supplied in the historical discovery fit",
    "Audit reconstruction is memory-only; locked matrix is never overwritten",
    format(max(abs(reproduction_difference)), scientific = TRUE, digits = 5),
    format(sqrt(mean(reproduction_difference^2)), scientific = TRUE, digits = 5),
    sprintf("%.12f", cor(as.vector(rederived_after), as.vector(after))),
    sprintf("%.12f", cor(as.vector(before_residual), as.vector(after_residual)))
  )
)
fwrite(correction_contract, file.path(table_dir, "step1_batch_correction_contract.csv"))

pca_one <- function(matrix, state) {
  fit <- prcomp(t(matrix), center = TRUE, scale. = TRUE, rank. = 10L)
  scores <- as.data.table(fit$x[, seq_len(10L), drop = FALSE], keep.rownames = "Sample")
  scores[, `:=`(CorrectionState = state, Study = sample_info$Study[match(Sample, sample_info$Sample)])]
  variance <- 100 * fit$sdev^2 / sum(fit$sdev^2)
  association <- rbindlist(lapply(seq_len(10L), function(index) {
    score <- fit$x[, index]
    fit_aov <- aov(score ~ factor(sample_info$Study))
    table_aov <- summary(fit_aov)[[1L]]
    data.table(
      CorrectionState = state,
      PC = paste0("PC", index),
      PCIndex = index,
      VarianceExplainedPercent = variance[index],
      StudyEtaSquared = table_aov$`Sum Sq`[1L] / sum(table_aov$`Sum Sq`),
      PValue = table_aov$`Pr(>F)`[1L]
    )
  }))
  list(scores = scores, association = association)
}

pca_before <- pca_one(before, "Before correction")
pca_after <- pca_one(after, "Locked limma correction")
pca_scores <- rbindlist(list(pca_before$scores, pca_after$scores), use.names = TRUE)
pca_association <- rbindlist(list(pca_before$association, pca_after$association), use.names = TRUE)
fwrite(pca_scores, file.path(table_dir, "step1_batch_pca_scores.csv"))
fwrite(pca_association, file.path(table_dir, "step1_batch_pc_study_association.csv"))

study_palette <- fread(file.path(root, "0-Data", "qc", "study_color_dictionary.csv"))
study_colors <- setNames(study_palette$ColorHex, study_palette$Study)
pca_plot_data <- pca_scores[, .(Sample, Study, CorrectionState, PC1, PC2)]
pca_centroids <- pca_plot_data[, .(PC1 = mean(PC1), PC2 = mean(PC2)), by = .(CorrectionState, Study)]
variance_labels <- pca_association[, .(
  PC1Variance = VarianceExplainedPercent[PCIndex == 1L],
  PC2Variance = VarianceExplainedPercent[PCIndex == 2L]
), by = CorrectionState]
pca_plot_data <- merge(pca_plot_data, variance_labels, by = "CorrectionState", all.x = TRUE)

p_pca <- ggplot(pca_plot_data, aes(PC1, PC2, color = Study)) +
  geom_point(size = 1.4, alpha = 0.68) +
  geom_point(data = pca_centroids, size = 3.2, shape = 4, stroke = 1.1, show.legend = FALSE) +
  facet_wrap(~CorrectionState, scales = "free", nrow = 1L) +
  scale_color_manual(values = study_colors) +
  labs(
    title = "Study-associated structure before and after the locked correction",
    subtitle = "Crosses denote study centroids; panels use independent PCA bases",
    x = "PC1", y = "PC2", color = "Study"
  ) +
  theme_classic(base_size = 9) +
  theme(legend.position = "right", strip.text = element_text(face = "bold"))

p_association <- ggplot(pca_association, aes(factor(PCIndex), StudyEtaSquared, fill = CorrectionState)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  geom_point(aes(color = CorrectionState), position = position_dodge(width = 0.75), size = 1.4, show.legend = FALSE) +
  scale_fill_manual(values = c("Before correction" = "#7A7A7A", "Locked limma correction" = "#0072B2")) +
  scale_color_manual(values = c("Before correction" = "#7A7A7A", "Locked limma correction" = "#0072B2")) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.03))) +
  labs(
    title = "Study association across the first ten principal components",
    subtitle = "Eta-squared quantifies the fraction of each PC score explained by Study",
    x = "Principal component", y = "Study eta-squared", fill = NULL
  ) +
  theme_classic(base_size = 9) +
  theme(legend.position = "top")

for (extension in c("pdf", "png")) {
  ggsave(file.path(figure_dir, paste0("step1_batch_pca_before_after.", extension)),
         p_pca, width = 10.2, height = 4.9, dpi = 300, bg = "white")
  ggsave(file.path(figure_dir, paste0("step1_batch_pc_study_association.", extension)),
         p_association, width = 7.6, height = 4.3, dpi = 300, bg = "white")
}

# Feature selection and a bounded sensitivity comparison for batch handling.
locked_feature_table <- fread(file.path(data_dir, "top_variable_genes.csv"))
locked_genes <- as.character(locked_feature_table[[1L]])
after_variance <- apply(after, 1L, var)
before_variance <- apply(before, 1L, var)
before_top_genes <- names(sort(before_variance, decreasing = TRUE))[seq_len(2000L)]
feature_summary <- data.table(
  Item = c(
    "Rank rule", "Selected genes", "Locked variance cutoff", "Locked genes present in corrected matrix",
    "Overlap with uncorrected top 2,000"
  ),
  Value = c(
    "Variance across the 586 corrected discovery samples", "2000",
    format(min(after_variance[locked_genes]), digits = 8),
    as.character(sum(locked_genes %in% rownames(after))),
    as.character(length(intersect(locked_genes, before_top_genes)))
  )
)
fwrite(feature_summary, file.path(table_dir, "step1_feature_selection_summary.csv"))

variance_plot_data <- data.table(Gene = names(after_variance), Variance = after_variance)
p_variance <- ggplot(variance_plot_data, aes(Variance)) +
  geom_histogram(bins = 80, fill = "#8FBBD9", color = "white", linewidth = 0.15) +
  geom_vline(xintercept = min(after_variance[locked_genes]), color = "#D55E00", linetype = 2) +
  scale_x_log10() +
  labs(
    title = "Variance-based feature selection for unsupervised discovery",
    subtitle = "Dashed line marks the locked top-2,000-gene cutoff",
    x = "Gene variance across corrected samples (log10 scale)", y = "Genes"
  ) +
  theme_classic(base_size = 9)
for (extension in c("pdf", "png")) {
  ggsave(file.path(figure_dir, paste0("step1_feature_selection_variance.", extension)),
         p_variance, width = 6.8, height = 4.3, dpi = 300, bg = "white")
}

kmeans_fit <- function(matrix, genes, label) {
  x <- t(matrix[genes, sample_info$Sample, drop = FALSE])
  set.seed(123)
  fit <- kmeans(x, centers = 3L, nstart = 25L, iter.max = 100L)
  data.table(
    Analysis = label,
    Features = length(genes),
    ARIvsLockedKMeans = adjusted_rand_index(fit$cluster, sample_info$KMeans),
    ARIvsLockedNMF = adjusted_rand_index(fit$cluster, sample_info$NMF),
    StudyCramersV = cramers_v(fit$cluster, sample_info$Study),
    SmallestCluster = min(table(fit$cluster)),
    LargestCluster = max(table(fit$cluster))
  )
}

batch_cluster_sensitivity <- rbindlist(list(
  kmeans_fit(before, locked_genes, "Uncorrected matrix; locked 2,000 genes"),
  kmeans_fit(before, before_top_genes, "Uncorrected matrix; own top 2,000 genes"),
  kmeans_fit(after, locked_genes, "Locked limma correction; locked 2,000 genes")
))
fwrite(batch_cluster_sensitivity, file.path(table_dir, "step1_batch_handling_cluster_sensitivity.csv"))

cat("Step1 batch/feature audit completed without writing either locked matrix.\n")
