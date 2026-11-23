#!/usr/bin/env Rscript

# Reader-facing method-selection and state-definition figures from saved Step1
# source tables. This script performs no clustering and modifies no inputs.

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(readxl)
})

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
table_dir <- file.path(root, "1-Clustering", "publication", "tables")
figure_dir <- file.path(root, "1-Clustering", "publication", "figures")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

consensus <- fread(file.path(table_dir, "step1_standard_consensus_pac_cdf.csv"))
consensus_sizes <- fread(file.path(table_dir, "step1_standard_consensus_cluster_sizes.csv"))
nmf <- fread(file.path(table_dir, "step1_nmf_rank_sensitivity_10run.csv"))
kmeans <- fread(file.path(table_dir, "step1_kmeans_rank_sensitivity.csv"))

method_metrics <- rbindlist(list(
  consensus[, .(Method = "Consensus", K, Metric = "PAC (lower is better)", Value = StandardPAC_0.1_0.9)],
  consensus_sizes[, .(Value = max(Fraction)), by = K][, .(Method = "Consensus", K, Metric = "Largest cluster fraction", Value)],
  nmf[, .(Method = "NMF", K, Metric = "Cophenetic (higher is better)", Value = Cophenetic)],
  nmf[, .(Method = "NMF", K, Metric = "PAC (lower is better)", Value = StandardPAC_0.1_0.9)],
  kmeans[, .(Method = "KMeans", K, Metric = "Mean silhouette (higher is better)", Value = MeanSilhouette)]
), use.names = TRUE)
method_metrics[, Facet := paste(Method, Metric, sep = " | ")]
method_metrics[, ScaleFloor := fifelse(grepl("Largest cluster fraction", Metric), 0, NA_real_)]
fwrite(method_metrics, file.path(table_dir, "step1_method_selection_metrics_long.csv"))

p_selection <- ggplot(method_metrics, aes(K, Value, group = 1L)) +
  geom_blank(aes(y = ScaleFloor)) +
  geom_vline(xintercept = 3, color = "#B8B8B8", linetype = 2) +
  geom_line(color = "#3B5F8A", linewidth = 0.65) +
  geom_point(aes(fill = K == 3), shape = 21, color = "#3B5F8A", size = 2.2, stroke = 0.45) +
  scale_fill_manual(values = c(`FALSE` = "white", `TRUE` = "#E69F00"), guide = "none") +
  scale_x_continuous(breaks = 2:6) +
  facet_wrap(~Facet, scales = "free_y", ncol = 3) +
  labs(
    title = "Candidate cluster structure is evaluated separately by method",
    subtitle = "The vertical line marks k=3; degenerate consensus structure is disclosed rather than counted as supporting evidence",
    x = "Candidate number of clusters (k)", y = NULL
  ) +
  theme_classic(base_size = 9) +
  theme(strip.text = element_text(size = 8.2, face = "bold"))

pairwise <- fread(file.path(table_dir, "step1_stored_method_pairwise_agreement.csv"))
method_levels <- c("ConsensusCluster", "NMF", "KMeans", "HC")
agreement_grid <- CJ(MethodA = method_levels, MethodB = method_levels)
agreement_grid[, AdjustedRandIndex := fifelse(MethodA == MethodB, 1, NA_real_)]
for (index in seq_len(nrow(pairwise))) {
  agreement_grid[MethodA == pairwise$MethodA[index] & MethodB == pairwise$MethodB[index], AdjustedRandIndex := pairwise$AdjustedRandIndex[index]]
  agreement_grid[MethodA == pairwise$MethodB[index] & MethodB == pairwise$MethodA[index], AdjustedRandIndex := pairwise$AdjustedRandIndex[index]]
}
fwrite(agreement_grid, file.path(table_dir, "step1_method_agreement_matrix_long.csv"))
agreement_grid[, `:=`(
  MethodA = factor(MethodA, levels = method_levels),
  MethodB = factor(MethodB, levels = rev(method_levels))
)]
p_agreement <- ggplot(agreement_grid, aes(MethodA, MethodB, fill = AdjustedRandIndex)) +
  geom_tile(color = "white", linewidth = 0.55) +
  geom_text(aes(label = sprintf("%.2f", AdjustedRandIndex)), size = 3) +
  scale_fill_gradient2(low = "#F2F2F2", mid = "#F6C65B", high = "#0072B2", midpoint = 0.5, limits = c(0, 1)) +
  coord_equal() +
  labs(
    title = "Stored method agreement identifies NMF and KMeans as the operational pair",
    subtitle = "Adjusted Rand index; Consensus/HC are degenerate two-cluster results",
    x = NULL, y = NULL, fill = "ARI"
  ) +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

method_disposition <- data.table(
  Method = c("ConsensusClusterPlus", "NMF", "KMeans", "Hierarchical clustering", "Historical four-method ensemble"),
  CandidateK = c("2-6", "3", "3", "2", "6"),
  Decision = c("EXCLUDED_DEGENERATE", "RETAINED_OPERATIONAL", "RETAINED_OPERATIONAL", "EXCLUDED_DEGENERATE", "HISTORICAL_ONLY"),
  Evidence = c(
    sprintf("Matched standard k=3 largest cluster %.1f%%", 100 * consensus_sizes[K == 3, max(Fraction)]),
    sprintf("k=3 cophenetic %.3f; PAC %.3f", nmf[K == 3, Cophenetic], nmf[K == 3, StandardPAC_0.1_0.9]),
    sprintf("k=3 mean silhouette %.3f; ARI vs locked KMeans %.3f", kmeans[K == 3, MeanSilhouette], kmeans[K == 3, ARIvsHistoricalKMeans]),
    "Stored solution has k=2 and matches the degenerate consensus partition",
    "Not used to define C1/C2/C3; retained only as historical provenance"
  )
)
fwrite(method_disposition, file.path(table_dir, "step1_method_disposition.csv"))

mapping <- fread(file.path(table_dir, "step1_nmf_kmeans_label_mapping.csv"))
definition <- merge(CJ(NMF = 1:3, KMeans = 1:3), mapping, by = c("NMF", "KMeans"), all.x = TRUE)
definition[is.na(N), `:=`(N = 0L, Ensemble = "None", StateClass = "Empty")]
palette <- fread(file.path(root, "docs", "subtype_color_dictionary.csv"))
state_colors <- setNames(palette$ColorHex, palette$FinalSubtype)
state_colors <- c(state_colors, None = "#F2F2F2")
definition[, `:=`(
  NMFLabel = factor(NMF, levels = 3:1, labels = c("NMF 3 (C3-aligned)", "NMF 2 (C2-aligned)", "NMF 1 (C1-aligned)")),
  KMeansLabel = factor(KMeans, levels = 1:3, labels = c("KMeans 1 (C2-aligned)", "KMeans 2 (C3-aligned)", "KMeans 3 (C1-aligned)")),
  Annotation = fifelse(N > 0, paste0(Ensemble, "\nn=", N), "n=0")
)]
fwrite(definition[, .(NMF, KMeans, Ensemble, N, StateClass)], file.path(table_dir, "step1_nmf_kmeans_state_definition_full.csv"))
p_definition <- ggplot(definition, aes(KMeansLabel, NMFLabel, fill = Ensemble)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = Annotation), size = 3.2, lineheight = 0.95) +
  scale_fill_manual(values = state_colors, drop = FALSE) +
  coord_equal() +
  labs(
    title = "Concordant NMF-KMeans assignments define principal states",
    subtitle = "Off-diagonal disagreements are retained as explicit pairwise intermediate states",
    x = NULL, y = NULL, fill = "Final state"
  ) +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 25, hjust = 1))

for (extension in c("pdf", "png")) {
  ggsave(file.path(figure_dir, paste0("step1_method_selection_summary.", extension)),
         p_selection, width = 10.2, height = 6.0, dpi = 300, bg = "white")
  ggsave(file.path(figure_dir, paste0("step1_cross_method_agreement.", extension)),
         p_agreement, width = 6.6, height = 5.4, dpi = 300, bg = "white")
  ggsave(file.path(figure_dir, paste0("step1_nmf_kmeans_state_definition.", extension)),
         p_definition, width = 8.2, height = 6.0, dpi = 300, bg = "white")
}

cat("Step1 reader-facing discovery figures generated from saved source tables.\n")
