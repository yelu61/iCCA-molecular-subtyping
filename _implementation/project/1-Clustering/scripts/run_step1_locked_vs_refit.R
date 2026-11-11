#!/usr/bin/env Rscript

# Independent refit audit for the locked Step 1 subtype interface. This script
# never overwrites the locked expression matrix or sample assignments.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(NMF)
  library(readxl)
})

set.seed(123)
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
publication_dir <- file.path(root, "1-Clustering", "publication")
table_dir <- file.path(publication_dir, "tables")
figure_dir <- file.path(publication_dir, "figures")
model_dir <- file.path(publication_dir, "models")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

load_rdata_object <- function(path, object) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  if (!exists(object, envir = env, inherits = FALSE)) stop("Missing object ", object)
  get(object, envir = env, inherits = FALSE)
}

choose2 <- function(x) x * (x - 1) / 2

adjusted_rand_index <- function(x, y) {
  tab <- table(x, y)
  n <- sum(tab)
  index <- sum(choose2(tab))
  a <- sum(choose2(rowSums(tab)))
  b <- sum(choose2(colSums(tab)))
  expected <- a * b / choose2(n)
  maximum <- (a + b) / 2
  if (maximum == expected) return(1)
  (index - expected) / (maximum - expected)
}

normalized_mutual_information <- function(x, y) {
  tab <- table(x, y)
  pxy <- tab / sum(tab)
  px <- rowSums(pxy)
  py <- colSums(pxy)
  nz <- which(pxy > 0, arr.ind = TRUE)
  mi <- sum(vapply(seq_len(nrow(nz)), function(i) {
    r <- nz[i, 1L]
    c <- nz[i, 2L]
    pxy[r, c] * log(pxy[r, c] / (px[r] * py[c]))
  }, numeric(1)))
  hx <- -sum(px[px > 0] * log(px[px > 0]))
  hy <- -sum(py[py > 0] * log(py[py > 0]))
  if (hx == 0 || hy == 0) return(1)
  mi / sqrt(hx * hy)
}

permutations <- function(x) {
  if (length(x) == 1L) return(matrix(x, nrow = 1L))
  do.call(rbind, lapply(seq_along(x), function(i) {
    cbind(x[i], permutations(x[-i]))
  }))
}

align_to_locked <- function(refit, locked) {
  refit_levels <- sort(unique(refit))
  locked_levels <- sort(unique(locked))
  if (length(refit_levels) != length(locked_levels)) stop("Cluster number mismatch")
  candidate <- permutations(locked_levels)
  scores <- apply(candidate, 1L, function(target) {
    mapped <- target[match(refit, refit_levels)]
    sum(mapped == locked)
  })
  best <- candidate[which.max(scores), ]
  mapping <- setNames(best, refit_levels)
  list(labels = unname(mapping[as.character(refit)]), mapping = mapping, matches = max(scores))
}

pair_state <- function(nmf_cluster, kmeans_cluster) {
  nmf_identity <- c(`1` = "C1", `2` = "C2", `3` = "C3")[as.character(nmf_cluster)]
  kmeans_identity <- c(`1` = "C2", `2` = "C3", `3` = "C1")[as.character(kmeans_cluster)]
  vapply(seq_along(nmf_identity), function(i) {
    pair <- unique(c(nmf_identity[i], kmeans_identity[i]))
    if (length(pair) == 1L) return(pair)
    paste(sort(pair), collapse = "_")
  }, character(1))
}

expression <- load_rdata_object(
  file.path(root, "1-Clustering", "data", "combined_expr_after_correction.rds"),
  "corrected_expr_limma"
)
top_variable <- fread(file.path(root, "1-Clustering", "data", "top_variable_genes.csv"))
genes <- intersect(as.character(top_variable[[1L]]), rownames(expression))
sample_info <- as.data.table(read_excel(
  file.path(root, "1-Clustering", "data", "sample_info_with_clusters.xlsx")
))
sample_info[, Sample := as.character(Sample)]
sample_info <- sample_info[match(colnames(expression), Sample)]
stopifnot(length(genes) == 2000L, nrow(sample_info) == 586L, !anyNA(sample_info$Sample))

x <- t(expression[genes, sample_info$Sample, drop = FALSE])

# Exact KMeans refit using the locked feature set and historical parameters.
set.seed(123)
kmeans_fit <- kmeans(x, centers = 3L, nstart = 25L, iter.max = 100L)
saveRDS(kmeans_fit, file.path(model_dir, "step1_kmeans_k3_exact_refit.rds"))
kmeans_alignment <- align_to_locked(kmeans_fit$cluster, sample_info$KMeans)
kmeans_aligned <- as.integer(kmeans_alignment$labels)
kmeans_distances <- sapply(seq_len(3L), function(k) {
  rowSums((x - matrix(kmeans_fit$centers[k, ], nrow(x), ncol(x), byrow = TRUE))^2)
})
kmeans_sorted <- t(apply(kmeans_distances, 1L, sort))
kmeans_margin <- (kmeans_sorted[, 2L] - kmeans_sorted[, 1L]) / pmax(kmeans_sorted[, 2L], .Machine$double.eps)

# Exact k=3 NMF refit. The run count defaults to the historical final 50 runs.
nmf_runs <- as.integer(Sys.getenv("SCLUCID_STEP1_NMF_RUNS", "50"))
if (!is.finite(nmf_runs) || nmf_runs < 1L) stop("Invalid SCLUCID_STEP1_NMF_RUNS")
nmf_data <- t(x)
nmf_data[nmf_data < 0] <- 0
nmf_data <- nmf_data + 1e-6
nmf_fit <- NMF::nmf(
  nmf_data, rank = 3L, method = "brunet", nrun = nmf_runs, seed = 123,
  .options = ""
)
saveRDS(nmf_fit, file.path(model_dir, paste0("step1_nmf_k3_refit_", nmf_runs, "run.rds")))
nmf_refit <- as.integer(NMF::predict(nmf_fit))
nmf_alignment <- align_to_locked(nmf_refit, sample_info$NMF)
nmf_aligned <- as.integer(nmf_alignment$labels)
nmf_coefficients <- NMF::coef(nmf_fit)
nmf_membership <- sweep(nmf_coefficients, 2L, pmax(colSums(nmf_coefficients), .Machine$double.eps), "/")
nmf_sorted <- apply(nmf_membership, 2L, sort, decreasing = TRUE)
nmf_top_membership <- nmf_sorted[1L, ]
nmf_margin <- nmf_sorted[1L, ] - nmf_sorted[2L, ]

refit_final <- pair_state(nmf_aligned, kmeans_aligned)
locked_final_from_methods <- pair_state(sample_info$NMF, sample_info$KMeans)
stopifnot(identical(locked_final_from_methods, as.character(sample_info$Ensemble)))

comparison <- data.table(
  Sample = sample_info$Sample,
  Study = sample_info$Study,
  LockedNMF = as.integer(sample_info$NMF),
  RefitNMFAligned = nmf_aligned,
  NMFAgreement = nmf_aligned == sample_info$NMF,
  NMFTopMembership = nmf_top_membership,
  NMFMembershipMargin = nmf_margin,
  LockedKMeans = as.integer(sample_info$KMeans),
  RefitKMeansAligned = kmeans_aligned,
  KMeansAgreement = kmeans_aligned == sample_info$KMeans,
  KMeansDistanceMargin = kmeans_margin,
  LockedFinalSubtype = as.character(sample_info$Ensemble),
  RefitFinalSubtype = refit_final,
  FinalSubtypeAgreement = refit_final == as.character(sample_info$Ensemble)
)
comparison[, ReviewFlag := !NMFAgreement | !KMeansAgreement | !FinalSubtypeAgreement]
fwrite(comparison, file.path(table_dir, "step1_locked_vs_refit_sample_comparison.csv"))

summarize_comparison <- function(method, locked, refit) {
  data.table(
    Method = method,
    N = length(locked),
    ARI = adjusted_rand_index(locked, refit),
    NMI = normalized_mutual_information(locked, refit),
    ExactAgreementAfterAlignment = mean(locked == refit),
    ChangedSamples = sum(locked != refit)
  )
}

summary <- rbindlist(list(
  summarize_comparison(paste0("NMF k=3 exact refit (", nmf_runs, " runs)"), sample_info$NMF, nmf_aligned),
  summarize_comparison("KMeans k=3 exact refit", sample_info$KMeans, kmeans_aligned),
  summarize_comparison("Governed six-state assignment", sample_info$Ensemble, refit_final)
))
summary[, RefitRuns := c(nmf_runs, NA_integer_, NA_integer_)]
fwrite(summary, file.path(table_dir, "step1_locked_vs_refit_summary.csv"))

make_confusion <- function(method, locked, refit) {
  data.table(Method = method, Locked = as.character(locked), Refit = as.character(refit))[
    , .(Samples = .N), by = .(Method, Locked, Refit)
  ]
}
confusion <- rbindlist(list(
  make_confusion("NMF", sample_info$NMF, nmf_aligned),
  make_confusion("KMeans", sample_info$KMeans, kmeans_aligned),
  make_confusion("SixState", sample_info$Ensemble, refit_final)
))
fwrite(confusion, file.path(table_dir, "step1_locked_vs_refit_confusion_long.csv"))

study_switches <- comparison[, .(
  N = .N,
  NMFChanges = sum(!NMFAgreement),
  KMeansChanges = sum(!KMeansAgreement),
  FinalSubtypeChanges = sum(!FinalSubtypeAgreement),
  FinalSubtypeAgreement = mean(FinalSubtypeAgreement)
), by = Study][order(Study)]
fwrite(study_switches, file.path(table_dir, "step1_locked_vs_refit_by_study.csv"))

# KMeans feature-number sensitivity; feature ranking is recomputed from the
# locked corrected matrix without changing the locked 2,000-gene interface.
variance_rank <- names(sort(apply(expression, 1L, var), decreasing = TRUE))
feature_sensitivity <- rbindlist(lapply(c(1000L, 2000L, 5000L), function(n_features) {
  feature_genes <- head(variance_rank, n_features)
  feature_x <- t(expression[feature_genes, sample_info$Sample, drop = FALSE])
  set.seed(123)
  fit <- kmeans(feature_x, centers = 3L, nstart = 25L, iter.max = 100L)
  aligned <- align_to_locked(fit$cluster, sample_info$KMeans)$labels
  data.table(
    Features = n_features,
    ARIvsLockedKMeans = adjusted_rand_index(sample_info$KMeans, aligned),
    NMIvsLockedKMeans = normalized_mutual_information(sample_info$KMeans, aligned),
    ExactAgreementAfterAlignment = mean(sample_info$KMeans == aligned),
    ChangedSamples = sum(sample_info$KMeans != aligned),
    SmallestCluster = min(table(aligned)),
    LargestCluster = max(table(aligned))
  )
}))
fwrite(feature_sensitivity, file.path(table_dir, "step1_kmeans_feature_number_sensitivity.csv"))

state_levels <- c("C1", "C2", "C3", "C1_C2", "C1_C3", "C2_C3")
six_confusion <- confusion[Method == "SixState"]
six_confusion[, `:=`(
  Locked = factor(Locked, levels = state_levels),
  Refit = factor(Refit, levels = state_levels)
)]
p_confusion <- ggplot(six_confusion, aes(Refit, Locked, fill = Samples)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = ifelse(Samples > 0, Samples, "")), size = 3.1) +
  scale_fill_gradient(low = "#F7FBFF", high = "#08519C") +
  coord_equal() +
  labs(
    title = "Locked versus independently refitted six-state assignments",
    subtitle = paste0("NMF refit: ", nmf_runs, " runs; KMeans refit: seed 123, 25 starts"),
    x = "Refitted assignment", y = "Locked assignment", fill = "Samples"
  ) +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))

uncertainty <- rbindlist(list(
  comparison[, .(Sample, Study, FinalSubtype = LockedFinalSubtype, Method = "NMF membership margin", Margin = NMFMembershipMargin, Changed = !NMFAgreement)],
  comparison[, .(Sample, Study, FinalSubtype = LockedFinalSubtype, Method = "KMeans distance margin", Margin = KMeansDistanceMargin, Changed = !KMeansAgreement)]
))
p_uncertainty <- ggplot(uncertainty, aes(Changed, Margin, fill = Changed)) +
  geom_boxplot(width = 0.55, outlier.size = 0.7) +
  facet_wrap(~Method, scales = "free_y", nrow = 1L) +
  scale_x_discrete(labels = c(`FALSE` = "Stable", `TRUE` = "Changed")) +
  scale_fill_manual(values = c(`FALSE` = "#4D4D4D", `TRUE` = "#D55E00"), guide = "none") +
  labs(
    title = "Refit changes concentrate in lower-certainty samples",
    x = NULL, y = "Assignment margin"
  ) +
  theme_classic(base_size = 9)

for (extension in c("pdf", "png")) {
  ggsave(file.path(figure_dir, paste0("step1_locked_vs_refit_confusion.", extension)),
         p_confusion, width = 6.4, height = 5.2, dpi = 300)
  ggsave(file.path(figure_dir, paste0("step1_refit_assignment_uncertainty.", extension)),
         p_uncertainty, width = 7.2, height = 4.2, dpi = 300)
}

# Redraw the public-facing panels from the saved source tables. This guarantees
# complete zero cells in the confusion matrix and avoids a changed/stable split
# when the exact refit contains no changed samples.
sys.source(
  file.path(root, "1-Clustering", "scripts", "plot_step1_clustering_refit.R"),
  envir = new.env(parent = globalenv())
)

cat("Step 1 locked-vs-refit audit completed with ", nmf_runs, " NMF runs\n", sep = "")
print(summary)
