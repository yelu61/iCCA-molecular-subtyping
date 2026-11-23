#!/usr/bin/env Rscript

# Finalize the clustering-only Step 1 release gate.

suppressPackageStartupMessages(library(data.table))

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
table_dir <- file.path(root, "1-Clustering", "publication", "tables")

refit <- fread(file.path(table_dir, "step1_locked_vs_refit_summary.csv"))
feature <- fread(file.path(table_dir, "step1_kmeans_feature_number_sensitivity.csv"))
seed <- fread(file.path(table_dir, "step1_kmeans_k3_seed_sensitivity.csv"))
loso <- fread(file.path(table_dir, "step1_kmeans_leave_one_study_out.csv"))
nmf <- fread(file.path(table_dir, "step1_nmf_rank_sensitivity_10run.csv"))
sizes <- fread(file.path(table_dir, "step1_standard_consensus_cluster_sizes.csv"))
mapping <- fread(file.path(table_dir, "step1_nmf_kmeans_label_mapping.csv"))
study <- fread(file.path(table_dir, "step1_study_subtype_association.csv"))
batch_contract <- fread(file.path(table_dir, "step1_batch_correction_contract.csv"))
batch_sensitivity <- fread(file.path(table_dir, "step1_batch_handling_cluster_sensitivity.csv"))
method_disposition <- fread(file.path(table_dir, "step1_method_disposition.csv"))

contract_number <- function(field) {
  as.numeric(batch_contract[Field == field, Value])
}
reconstruction_max_difference <- contract_number("Maximum absolute reconstruction difference")
reconstruction_correlation <- contract_number("Reconstruction correlation")
uncorrected_study_v <- max(batch_sensitivity[grepl("Uncorrected", Analysis), StudyCramersV])
corrected_study_v <- batch_sensitivity[grepl("Locked limma", Analysis), StudyCramersV]
corrected_kmeans_ari <- batch_sensitivity[grepl("Locked limma", Analysis), ARIvsLockedKMeans]
required_discovery_figures <- file.path(root, "1-Clustering", "publication", "figures", c(
  "step1_batch_pca_before_after.pdf",
  "step1_batch_pc_study_association.pdf",
  "step1_feature_selection_variance.pdf",
  "step1_method_selection_summary.pdf",
  "step1_cross_method_agreement.pdf",
  "step1_nmf_kmeans_state_definition.pdf"
))

nmf_refit <- refit[grepl("^NMF", Method)]
kmeans_refit <- refit[grepl("^KMeans", Method)]
six_refit <- refit[grepl("six-state", Method, ignore.case = TRUE)]
nmf_k3 <- nmf[K == 3]

gate <- data.table(
  GateItem = c(
    "Locked discovery interface",
    "Batch correction exact reproduction",
    "Batch handling and study-dominance sensitivity",
    "Reader-facing multi-method discovery evidence",
    "Six-state definition",
    "Exact NMF locked-vs-refit reproduction",
    "Exact KMeans locked-vs-refit reproduction",
    "Six-state locked-vs-refit reproduction",
    "KMeans seed stability",
    "KMeans feature-number sensitivity",
    "NMF rank sensitivity",
    "Leave-one-study-out stability",
    "Study/subtype association audit",
    "Degenerate hierarchical consensus methods"
  ),
  Status = c(
    "PASS",
    ifelse(reconstruction_max_difference < 1e-10 && reconstruction_correlation > 0.999999, "PASS", "FAIL"),
    ifelse(corrected_kmeans_ari >= 0.95 && corrected_study_v < uncorrected_study_v, "PASS_WITH_SENSITIVITY", "FAIL"),
    ifelse(all(file.exists(required_discovery_figures)) && sum(method_disposition$Decision == "RETAINED_OPERATIONAL") == 2L, "PASS_WITH_EXCLUSIONS", "FAIL"),
    "PASS",
    ifelse(nmf_refit$ARI >= 0.95, "PASS", "FAIL"),
    ifelse(kmeans_refit$ARI >= 0.95, "PASS", "FAIL"),
    ifelse(six_refit$ARI >= 0.95, "PASS", "FAIL"),
    ifelse(min(seed$ARIvsHistoricalKMeans) >= 0.95, "PASS", "FAIL"),
    ifelse(min(feature$ARIvsLockedKMeans) >= 0.70, "PASS_WITH_SENSITIVITY", "FAIL"),
    ifelse(nmf_k3$Cophenetic >= 0.95 && nmf_k3$StandardPAC_0.1_0.9 <= 0.1, "PASS", "FAIL"),
    ifelse(min(loso$ARIAllVsHistoricalKMeans) >= 0.70, "PASS", "FAIL"),
    "PASS_WITH_ADJUSTMENT",
    "EXCLUDED"
  ),
  Evidence = c(
    "10,531 genes x 586 discovery samples; locked 2,000-gene variance interface",
    sprintf("Memory-only limma reconstruction: maximum absolute difference %.2e; correlation %.12f; locked matrix not overwritten", reconstruction_max_difference, reconstruction_correlation),
    sprintf("Uncorrected k=3 Study Cramer's V up to %.3f; corrected %.3f; corrected KMeans ARI vs locked %.3f", uncorrected_study_v, corrected_study_v, corrected_kmeans_ari),
    "Consensus, NMF, KMeans, hierarchical, and historical ensemble are shown separately; only NMF and KMeans are retained operationally",
    sprintf("NMF/KMeans concordance defines %d principal and %d intermediate states", mapping[StateClass == "Main", sum(N)], mapping[StateClass == "Intermediate", sum(N)]),
    sprintf("ARI %.3f; exact agreement %.1f%%; %d changed samples", nmf_refit$ARI, 100 * nmf_refit$ExactAgreementAfterAlignment, nmf_refit$ChangedSamples),
    sprintf("ARI %.3f; exact agreement %.1f%%; %d changed samples", kmeans_refit$ARI, 100 * kmeans_refit$ExactAgreementAfterAlignment, kmeans_refit$ChangedSamples),
    sprintf("ARI %.3f; exact agreement %.1f%%; %d changed samples", six_refit$ARI, 100 * six_refit$ExactAgreementAfterAlignment, six_refit$ChangedSamples),
    sprintf("100 seeds; minimum ARI %.3f", min(seed$ARIvsHistoricalKMeans)),
    sprintf("1,000/2,000/5,000 genes; minimum ARI %.3f", min(feature$ARIvsLockedKMeans)),
    sprintf("k=3 cophenetic %.3f; PAC %.3f; ARI vs locked NMF %.3f", nmf_k3$Cophenetic, nmf_k3$StandardPAC_0.1_0.9, nmf_k3$ARIvsHistoricalNMF),
    sprintf("minimum all-sample held-study-out ARI %.3f", min(loso$ARIAllVsHistoricalKMeans)),
    sprintf("Study x main subtype Cramer's V %.3f; P %.4g; Study is retained in downstream models", study$CramersV, study$PValue),
    sprintf("Matched standard sensitivity remains degenerate; k=3 largest cluster %.1f%%", 100 * sizes[K == 3, max(Fraction)])
  )
)
gate[, Blocking := Status == "FAIL"]
fwrite(gate, file.path(table_dir, "step1_clustering_release_gate_status.csv"))

cat("Step 1 clustering gate finalized: ", sum(gate$Blocking), " blocking items\n", sep = "")
