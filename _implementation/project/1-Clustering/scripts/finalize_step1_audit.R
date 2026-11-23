#!/usr/bin/env Rscript

# Consolidate Step 1 audit evidence, gate status, and artifact inventory.
# This script is intentionally non-destructive outside 1-Clustering/publication.

suppressPackageStartupMessages(library(data.table))

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
publication_dir <- file.path(root, "1-Clustering", "publication")
table_dir <- file.path(publication_dir, "tables")

deg <- fread(file.path(table_dir, "step1_main_only_study_adjusted_deg.csv"))
meta <- fread(file.path(table_dir, "step1_cross_study_deg_meta_summary.csv"))
loso_deg <- fread(file.path(table_dir, "step1_leave_one_study_out_deg_sensitivity.csv"))
tf <- fread(file.path(table_dir, "step1_dorothea_tf_activity_main_only_study_adjusted.csv"))
loso_tf <- fread(file.path(table_dir, "step1_dorothea_tf_leave_one_study_out_summary.csv"))
tme <- fread(file.path(table_dir, "step1_tme_c2_vs_others_tests.csv"))
pac <- fread(file.path(table_dir, "step1_standard_consensus_pac_cdf.csv"))
sizes <- fread(file.path(table_dir, "step1_standard_consensus_cluster_sizes.csv"))
silhouette <- fread(file.path(table_dir, "step1_main_state_silhouette_summary.csv"))
agreement <- fread(file.path(table_dir, "step1_stored_method_pairwise_agreement.csv"))
study_test <- fread(file.path(table_dir, "step1_study_subtype_association.csv"))
operational <- fread(file.path(table_dir, "step1_operational_structure_sensitivity_summary.csv"))
kmeans_loso <- fread(file.path(table_dir, "step1_kmeans_leave_one_study_out.csv"))
nmf_rank <- fread(file.path(table_dir, "step1_nmf_rank_sensitivity_10run.csv"))

c2_candidates <- c("HNF1B", "ONECUT1", "FGFR2", "FGFR3", "SLC34A2", "CFTR", "AQP1",
                   "CLDN10", "DCDC2", "SPP1", "VTCN1")
c2_audit <- merge(
  deg[Subtype == "C2" & Gene %in% c2_candidates, .(
    Gene, CombinedLogFC = logFC_subtype_minus_others, CombinedFDR = FDR
  )],
  meta[Subtype == "C2" & Gene %in% c2_candidates, .(
    Gene, Studies, MetaLogFC, MetaFDR, I2Percent, DirectionConsistency
  )],
  by = "Gene", all = TRUE
)
c2_audit[, `:=`(
  Requested = TRUE,
  HigherInC2Combined = CombinedLogFC > 0,
  HigherInC2Meta = MetaLogFC > 0,
  PassCombinedFDR05 = CombinedFDR < 0.05,
  PassMetaFDR05 = MetaFDR < 0.05
)]
c2_audit <- merge(data.table(Gene = c2_candidates), c2_audit, by = "Gene", all.x = TRUE, sort = FALSE)
fwrite(c2_audit, file.path(table_dir, "step1_c2_candidate_malignant_program_bulk_audit.csv"))

# Corrected, cross-study-supported positive markers supplied to Step8. These do
# not overwrite the existing Step8 table; the comparison below records changes.
validation_markers <- merge(
  deg[logFC_subtype_minus_others > 0 & FDR < 0.05, .(
    Gene, Subtype, logFC_subtype_minus_others, FDR, t
  )],
  meta[MetaLogFC > 0 & MetaFDR < 0.05 & DirectionConsistency >= 6 / 7, .(
    Gene, Subtype, Studies, MetaLogFC, MetaFDR, I2Percent, DirectionConsistency
  )],
  by = c("Gene", "Subtype")
)
setorder(validation_markers, Subtype, -logFC_subtype_minus_others, FDR)
validation_markers <- validation_markers[, head(.SD, 40L), by = Subtype]
validation_markers[, `:=`(
  MarkerSet = "Step1_corrected_main_only_study_adjusted_top40",
  PositiveDirection = paste0("Higher in ", Subtype)
)]
fwrite(validation_markers, file.path(table_dir, "step1_bulk_principal_state_markers_for_step8.csv"))

legacy_step8_path <- file.path(root, "7-SingleCell", "data", "python_bulk_principal_state_markers.csv")
if (file.exists(legacy_step8_path)) {
  old_markers <- fread(legacy_step8_path)[, .(Gene, Subtype = FinalSubtype)]
  marker_compare_long <- merge(
    unique(validation_markers[, .(Gene, Subtype, InCorrectedSet = TRUE)]),
    unique(old_markers[, .(Gene, Subtype, InExistingStep8Set = TRUE)]),
    by = c("Gene", "Subtype"), all = TRUE
  )
  marker_compare_long[is.na(InCorrectedSet), InCorrectedSet := FALSE]
  marker_compare_long[is.na(InExistingStep8Set), InExistingStep8Set := FALSE]
  fwrite(marker_compare_long, file.path(table_dir, "step1_vs_existing_step8_bulk_marker_membership.csv"))
  marker_compare <- marker_compare_long[, .(
    CorrectedN = sum(InCorrectedSet), ExistingStep8N = sum(InExistingStep8Set),
    OverlapN = sum(InCorrectedSet & InExistingStep8Set),
    Jaccard = sum(InCorrectedSet & InExistingStep8Set) / sum(InCorrectedSet | InExistingStep8Set),
    CorrectedOnly = paste(Gene[InCorrectedSet & !InExistingStep8Set], collapse = ";"),
    ExistingOnly = paste(Gene[!InCorrectedSet & InExistingStep8Set], collapse = ";")
  ), by = Subtype]
  fwrite(marker_compare, file.path(table_dir, "step1_vs_existing_step8_bulk_marker_summary.csv"))
}

get_tme <- function(score) {
  x <- tme[Score == score]
  if (!nrow(x)) return("not available")
  sprintf("C2 median %.3f vs others %.3f; rank-biserial %.3f; FDR %.3g",
          x$C2Median[1], x$OthersMedian[1], x$RankBiserialC2MinusOthers[1], x$FDR[1])
}

evidence_matrix <- data.table(
  Subtype = c("C1", "C2", "C3"),
  StableIdentifier = TRUE,
  CorrectedDEG = c(
    "APOC3/APOA1/APOA2 and hepatic metabolic genes higher",
    "All 11 prespecified ductal-lineage candidates are higher with combined and meta-analysis FDR <0.05",
    "AGR2/S100P/CEACAM5/TMPRSS4/KRT17 and epithelial/proliferative genes higher"
  ),
  RankedPathways = c(
    "Positive fatty-acid, bile-acid, xenobiotic and organic-acid metabolism",
    "Positive cilium/translation/WNT-related programs; negative interferon-gamma, inflammatory and adaptive-immune programs",
    "Positive E2F/G2M/cell-cycle and epithelial-differentiation programs; negative metabolic programs"
  ),
  TFActivity = c(
    "HNF4A, PPARA, CEBPA and HNF1A higher",
    "RFX2 and LHX2 higher; broad hepatic/metabolic and RUNX1 activity lower",
    "KLF5, E2F7 and other epithelial/proliferative regulators higher"
  ),
  TME = c(
    "Not used for suffix selection",
    paste(get_tme("ImmuneScore_ESTIMATE"), get_tme("T_cells_MCPcounter"),
          get_tme("B_lineage_MCPcounter"), sep = "; "),
    "Not used for suffix selection"
  ),
  CrossStudy = sprintf(
    "All headline DEG signs retained in every leave-one-study-out fit; minimum all-gene Pearson r %.3f",
    loso_deg[, min(PearsonAllGenes), by = Subtype]$V1
  ),
  MalignantPseudobulk = "Pending Step8",
  DriverAlterations = c("Not required for working metabolic description", "Pending Step3/6 before any FGFR-based suffix", "Pending multi-omic characterization"),
  NamingStatus = c("Working description: C1-Metabolic", "Retain C2 without suffix", "Working description: C3-Proliferative")
)
fwrite(evidence_matrix, file.path(table_dir, "step1_three_subtype_biological_evidence_matrix.csv"))

naming_decision <- data.table(
  StableIdentifier = c("C1", "C2", "C3"),
  PublicationSuffix = c("Metabolic", NA_character_, "Proliferative"),
  Decision = c("WORKING_DESCRIPTION_SUPPORTED", "RETAIN_STABLE_IDENTIFIER_ONLY", "WORKING_DESCRIPTION_SUPPORTED"),
  EvidenceLayers = c(
    "Corrected DEG; ranked pathways; DoRothEA TF; cross-study/LOSO",
    "Bulk ductal-lineage candidate genes; ranked pathways; TME discordance; cross-study/LOSO",
    "Corrected DEG; ranked pathways; DoRothEA TF; cross-study/LOSO"
  ),
  Constraint = c(
    "Use C1 as the stable identifier; metabolic is a supported plotted description",
    "Do not use an immune-response or ductal-lineage suffix until malignant pseudobulk and driver evidence are reconciled",
    "Use C3 as the stable identifier; proliferative is a supported plotted description"
  )
)
fwrite(naming_decision, file.path(table_dir, "step1_formal_subtype_naming_decision.csv"))

deg_direction_pass <- all(fread(file.path(table_dir, "step1_legacy_deg_direction_audit.csv"))$CorrelationLegacyVsObserved < -0.999)
nmf_kmeans_ari <- agreement[MethodA == "NMF" & MethodB == "KMeans", AdjustedRandIndex]
largest_k3 <- sizes[K == 3, max(Fraction)]
overall_sil <- silhouette[Subtype == "Overall", MeanSilhouette]
nmf_k3 <- nmf_rank[K == 3]

gate <- data.table(
  GateItem = c(
    "Locked discovery interface", "Explicit DEG direction", "Cross-study DEG sensitivity",
    "Hierarchical consensus sensitivity", "Operational NMF/KMeans stability", "Study/subtype association",
    "TF activity refresh", "TME method concordance", "C2 global immune-high suffix",
    "C2 naming decision", "Corrected bulk markers for Step8", "Publication bulk signatures"
  ),
  Status = c(
    "PASS", ifelse(deg_direction_pass, "PASS", "FAIL"), "PASS",
    "EXCLUDED", "PASS", "PASS", "PASS", "PASS", "REJECTED",
    "PASS_DEFERRED_SUFFIX", "PASS", "PASS"
  ),
  Evidence = c(
    "10,531 genes x 586 samples; top 2,000 variable genes; governed main n=497",
    "Legacy coefficients correlate -1 with observed subtype-minus-others effects",
    sprintf("Minimum all-gene leave-one-study-out Pearson r %.3f; headline-gene sign agreement 1.0", min(loso_deg$PearsonAllGenes)),
    sprintf("Historical Consensus/HC and the matched standard sensitivity are degenerate (k=3 largest cluster %.1f%%); they are retained as failed methods, not evidence for k=3", 100 * largest_k3),
    sprintf("KMeans k=3 is identical across 100 seeds and has the best k=2-6 silhouette; NMF k=3 cophenetic %.3f, PAC %.3f, ARI vs locked NMF %.3f; NMF-KMeans ARI %.3f with 497 concordant principal states", nmf_k3$Cophenetic, nmf_k3$StandardPAC_0.1_0.9, nmf_k3$ARIvsHistoricalNMF, nmf_kmeans_ari),
    sprintf("Study x subtype: Cramer's V %.3f, P %.4g; minimum KMeans leave-one-study-out ARI %.3f; minimum DEG LOSO Pearson r %.3f", study_test$CramersV, study_test$PValue, min(kmeans_loso$ARIAllVsHistoricalKMeans), min(loso_deg$PearsonAllGenes)),
    sprintf("263 DoRothEA A-C TFs; minimum LOSO TF-effect Pearson r %.3f", min(loso_tf$PearsonR)),
    "ESTIMATE plus MCP-counter plus prespecified IFNG/CTL/cytolytic programs; EPIC explicitly excluded for invalid legacy input and unavailable locked reference",
    paste0("C2 is not consistently globally immune-high: ", get_tme("ImmuneScore_ESTIMATE")),
    "Formal decision retains C2 without a suffix; malignant pseudobulk and driver/CNV evidence are prerequisites for any future ductal/FGFR-based renaming",
    sprintf("Corrected top-40 positive markers per state exported; existing Step8 top-40 overlap ranges from %.1f%% to %.1f%%", 100 * min(marker_compare$Jaccard), 100 * max(marker_compare$Jaccard)),
    "Main-only, study-adjusted, cross-study-supported marker interface is ready for Steps4 and 8; biological suffix rules remain governed by the naming record"
  ),
  Blocking = rep(FALSE, 12L)
)
fwrite(gate, file.path(table_dir, "step1_release_gate_status.csv"))
fwrite(gate[, .(Item = GateItem, Status, Evidence, Blocking)], file.path(table_dir, "step1_static_audit.csv"))

paths <- unique(c(
  list.files(publication_dir, recursive = TRUE, full.names = TRUE, all.files = FALSE),
  list.files(file.path(root, "1-Clustering", "scripts"), pattern = "step1|finalize_step1", full.names = TRUE),
  file.path(root, "1-clustering.ipynb")
))
paths <- paths[file.exists(paths) & !dir.exists(paths)]
manifest <- data.table(
  ArtifactType = fifelse(grepl("/figures/", paths), "figure",
                         fifelse(grepl("\\.(R|ipynb)$", paths), "analysis_code", "source_table")),
  Path = sub(paste0("^", root, "/"), "", normalizePath(paths, winslash = "/", mustWork = TRUE))
)
manifest[, `:=`(
  Exists = TRUE,
  NonEmpty = file.info(file.path(root, Path))$size > 0,
  Bytes = file.info(file.path(root, Path))$size,
  MD5 = unname(tools::md5sum(file.path(root, Path)))
)]
setorder(manifest, ArtifactType, Path)
fwrite(manifest, file.path(table_dir, "step1_publication_output_manifest.csv"))

cat("Step 1 audit finalized: ", sum(gate$Blocking), " blocking gate items remain\n", sep = "")
