#!/usr/bin/env Rscript

# Sensitivity analysis for the operational NMF/KMeans three-state structure.
# Outputs are diagnostic only and never overwrite locked subtype assignments.

suppressPackageStartupMessages({
  library(cluster)
  library(data.table)
  library(readxl)
})

options(width = 100)
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
table_dir <- file.path(root, "1-Clustering", "publication", "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

load_rdata_object <- function(path, object) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  get(object, envir = env, inherits = FALSE)
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

expression <- load_rdata_object(
  file.path(root, "1-Clustering", "data", "combined_expr_after_correction.rds"),
  "corrected_expr_limma"
)
genes <- fread(file.path(root, "1-Clustering", "data", "top_variable_genes.csv"))[[1L]]
sample_info <- as.data.table(read_excel(
  file.path(root, "1-Clustering", "data", "sample_info_with_clusters.xlsx")
))
sample_info[, Sample := as.character(Sample)]
x <- t(expression[genes, sample_info$Sample, drop = FALSE])
stopifnot(nrow(x) == 586L, ncol(x) == 2000L)
euclidean <- dist(x)

# KMeans k diagnostics under the explicit historical seed.
kmeans_rank <- rbindlist(lapply(2:6, function(k) {
  set.seed(123)
  fit <- kmeans(x, centers = k, nstart = 25, iter.max = 100)
  sil <- silhouette(fit$cluster, euclidean)
  data.table(
    K = k, Seed = 123L, NStart = 25L,
    TotalWithinSS = fit$tot.withinss,
    BetweenSSFraction = fit$betweenss / fit$totss,
    MeanSilhouette = mean(sil[, "sil_width"]),
    SmallestCluster = min(table(fit$cluster)),
    LargestCluster = max(table(fit$cluster)),
    ARIvsHistoricalKMeans = adjusted_rand_index(fit$cluster, sample_info$KMeans),
    ARIvsHistoricalNMF = adjusted_rand_index(fit$cluster, sample_info$NMF)
  )
}))
fwrite(kmeans_rank, file.path(table_dir, "step1_kmeans_rank_sensitivity.csv"))

kmeans_seed <- rbindlist(lapply(1:100, function(seed) {
  set.seed(seed)
  fit <- kmeans(x, centers = 3L, nstart = 25L, iter.max = 100L)
  data.table(
    Seed = seed,
    TotalWithinSS = fit$tot.withinss,
    ARIvsHistoricalKMeans = adjusted_rand_index(fit$cluster, sample_info$KMeans),
    ARIvsHistoricalNMF = adjusted_rand_index(fit$cluster, sample_info$NMF),
    SmallestCluster = min(table(fit$cluster)),
    LargestCluster = max(table(fit$cluster))
  )
}))
fwrite(kmeans_seed, file.path(table_dir, "step1_kmeans_k3_seed_sensitivity.csv"))

# Leave-one-study-out KMeans refitting. Held-out samples are assigned to the
# nearest training centroid; ARI is label-permutation invariant.
kmeans_loso <- rbindlist(lapply(sort(unique(sample_info$Study)), function(held_out) {
  train <- sample_info$Study != held_out
  set.seed(123)
  fit <- kmeans(x[train, , drop = FALSE], centers = 3L, nstart = 25L, iter.max = 100L)
  squared_distance <- sapply(seq_len(3L), function(k) {
    rowSums((x[!train, , drop = FALSE] - matrix(
      fit$centers[k, ], nrow = sum(!train), ncol = ncol(x), byrow = TRUE
    ))^2)
  })
  predicted <- integer(nrow(x))
  predicted[train] <- fit$cluster
  predicted[!train] <- max.col(-squared_distance, ties.method = "first")
  data.table(
    HeldOutStudy = held_out,
    TrainingN = sum(train),
    HeldOutN = sum(!train),
    ARIAllVsHistoricalKMeans = adjusted_rand_index(predicted, sample_info$KMeans),
    ARITrainingVsHistoricalKMeans = adjusted_rand_index(fit$cluster, sample_info$KMeans[train]),
    ARIHeldOutVsHistoricalKMeans = adjusted_rand_index(predicted[!train], sample_info$KMeans[!train]),
    ARIAllVsHistoricalNMF = adjusted_rand_index(predicted, sample_info$NMF),
    SmallestClusterAll = min(table(predicted)),
    LargestClusterAll = max(table(predicted))
  )
}))
fwrite(kmeans_loso, file.path(table_dir, "step1_kmeans_leave_one_study_out.csv"))

# NMF rank sensitivity. The historical non-negativity rule is reproduced:
# negative corrected values are truncated to zero and 1e-6 is added.
nmf_data <- t(x)
nmf_data[nmf_data < 0] <- 0
nmf_data <- nmf_data + 1e-6
nmf_rank <- rbindlist(lapply(2:6, function(k) {
  fits <- parallel::mclapply(1:10, function(run_seed) {
    fit <- NMF::nmf(
      nmf_data, rank = k, method = "brunet", nrun = 1L, seed = run_seed,
      .options = ""
    )
    list(labels = as.integer(NMF::predict(fit)), rss = as.numeric(NMF::rss(fit, target = nmf_data)))
  }, mc.cores = 4L, mc.preschedule = FALSE)
  labels_matrix <- do.call(cbind, lapply(fits, `[[`, "labels"))
  consensus <- Reduce(`+`, lapply(seq_len(ncol(labels_matrix)), function(i) {
    outer(labels_matrix[, i], labels_matrix[, i], FUN = "==")
  })) / ncol(labels_matrix)
  best_run <- which.min(vapply(fits, `[[`, numeric(1), "rss"))
  labels <- labels_matrix[, best_run]
  values <- consensus[upper.tri(consensus)]
  sil <- silhouette(labels, euclidean)
  consensus_distance <- as.dist(1 - consensus)
  consensus_tree <- hclust(consensus_distance, method = "average")
  data.table(
    K = k,
    Runs = 10L,
    Seed = 123L,
    Cophenetic = cor(consensus_distance, cophenetic(consensus_tree)),
    Dispersion = mean(4 * (consensus - 0.5)^2),
    StandardPAC_0.1_0.9 = mean(values > 0.1 & values < 0.9),
    MeanSilhouetteEuclidean = mean(sil[, "sil_width"]),
    SmallestCluster = min(table(labels)),
    LargestCluster = max(table(labels)),
    ARIvsHistoricalNMF = adjusted_rand_index(labels, sample_info$NMF),
    ARIvsHistoricalKMeans = adjusted_rand_index(labels, sample_info$KMeans)
  )
}))
fwrite(nmf_rank, file.path(table_dir, "step1_nmf_rank_sensitivity_10run.csv"))

summary <- data.table(
  Item = c(
    "KMeans k3 seed reproducibility", "KMeans leave-one-study-out",
    "KMeans best silhouette k", "NMF best cophenetic k", "NMF k3 reproducibility"
  ),
  Value = c(
    sprintf("median ARI vs locked KMeans %.3f; minimum %.3f",
            median(kmeans_seed$ARIvsHistoricalKMeans), min(kmeans_seed$ARIvsHistoricalKMeans)),
    sprintf("minimum all-sample ARI vs locked KMeans %.3f",
            min(kmeans_loso$ARIAllVsHistoricalKMeans)),
    as.character(kmeans_rank[which.max(MeanSilhouette), K]),
    as.character(nmf_rank[which.max(Cophenetic), K]),
    sprintf("ARI vs locked NMF %.3f; cophenetic %.3f; PAC %.3f",
            nmf_rank[K == 3, ARIvsHistoricalNMF], nmf_rank[K == 3, Cophenetic],
            nmf_rank[K == 3, StandardPAC_0.1_0.9])
  )
)
fwrite(summary, file.path(table_dir, "step1_operational_structure_sensitivity_summary.csv"))

cat("Step 1 operational structure sensitivity completed\n")
