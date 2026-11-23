#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(openxlsx)
})

options(stringsAsFactors = FALSE, warn = 1)
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
out_dir <- file.path(root, "governance")
data_dir <- file.path(out_dir, "data")
report_dir <- file.path(out_dir, "reports")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

assignment_path <- file.path(root, "1-Clustering/data/sample_info_with_clusters.xlsx")
dictionary_path <- file.path(root, "2-TranscriptomicCharacterization/data/subtype_label_color_dictionary.csv")
study_dictionary_path <- file.path(root, "0-Data/qc/study_color_dictionary.csv")
stopifnot(file.exists(assignment_path), file.exists(dictionary_path), file.exists(study_dictionary_path))

main_states <- c("C1", "C2", "C3")
intermediate_states <- c("C1_C2", "C1_C3", "C2_C3")
expected_states <- c(main_states, intermediate_states)
expected_counts <- c(C1 = 112L, C2 = 223L, C3 = 162L,
                     C1_C2 = 41L, C1_C3 = 39L, C2_C3 = 9L)

assignments <- as.data.table(read_excel(assignment_path))
if (!all(c("Sample", "Ensemble") %in% names(assignments))) {
  stop("Step1 assignment file lacks Sample or Ensemble")
}
assignments[, Sample := as.character(Sample)]
if (anyDuplicated(assignments$Sample)) stop("Step1 assignment file contains duplicate Sample values")
if (!setequal(unique(assignments$Ensemble), expected_states)) stop("Unexpected Step1 subtype state")

dictionary <- fread(dictionary_path, showProgress = FALSE)
required_dict <- c("FinalSubtype", "DisplayLabel", "SubtypeClass", "Color")
if (!all(required_dict %in% names(dictionary))) stop("Subtype dictionary is incomplete")
dictionary <- dictionary[match(expected_states, FinalSubtype)]
if (anyNA(dictionary$FinalSubtype)) stop("Subtype dictionary is missing a locked state")
if (any(grepl("Immune", dictionary$DisplayLabel, ignore.case = TRUE))) {
  stop("Canonical dictionary contains a prohibited C2 immune suffix")
}

master <- copy(assignments)
master[, OriginalEnsemble := Ensemble]
master[, FinalSubtype := Ensemble]
master[, SubtypeClass := fifelse(FinalSubtype %chin% main_states, "Main", "Intermediate")]
master[, MainSubtypeForPlot := fifelse(SubtypeClass == "Main", FinalSubtype, NA_character_)]
master[, IntermediatePairA := fifelse(SubtypeClass == "Intermediate",
                                      tstrsplit(FinalSubtype, "_", fixed = TRUE)[[1]], NA_character_)]
master[, IntermediatePairB := fifelse(SubtypeClass == "Intermediate",
                                      tstrsplit(FinalSubtype, "_", fixed = TRUE)[[2]], NA_character_)]
master <- merge(master, dictionary[, .(FinalSubtype, SubtypeDisplayLabel = DisplayLabel)],
                by = "FinalSubtype", all.x = TRUE, sort = FALSE)
setcolorder(master, c("Sample", "OriginalEnsemble", "FinalSubtype", "SubtypeClass",
                      "MainSubtypeForPlot", "IntermediatePairA", "IntermediatePairB",
                      "SubtypeDisplayLabel",
                      setdiff(names(master), c("Sample", "OriginalEnsemble", "FinalSubtype",
                                               "SubtypeClass", "MainSubtypeForPlot",
                                               "IntermediatePairA", "IntermediatePairB",
                                               "SubtypeDisplayLabel"))))
master <- master[match(assignments$Sample, Sample)]

counts <- master[, .(N = .N), by = .(FinalSubtype, SubtypeClass)]
counts[, ExpectedN := expected_counts[FinalSubtype]]
if (!all(counts$N == counts$ExpectedN)) stop("Locked subtype counts do not match Step1")
counts <- counts[match(expected_states, FinalSubtype)]

fwrite(master, file.path(data_dir, "subtype_master_table.csv"))
write.xlsx(master, file.path(data_dir, "subtype_master_table.xlsx"), overwrite = TRUE)
fwrite(master[SubtypeClass == "Main"], file.path(data_dir, "main_subtype_samples.csv"))
fwrite(master[SubtypeClass == "Intermediate"], file.path(data_dir, "intermediate_subtype_samples.csv"))
fwrite(counts, file.path(data_dir, "subtype_counts.csv"))
fwrite(dictionary, file.path(data_dir, "subtype_label_color_dictionary.csv"))

# Synchronize only lightweight dictionaries/tables consumed by the current
# visualization layer. No downstream statistics are regenerated here.
viz_table_dir <- file.path(root, "8-Visualization/tables")
if (dir.exists(viz_table_dir)) {
  fwrite(dictionary, file.path(viz_table_dir, "subtype_label_color_dictionary.csv"))
  plot_counts <- merge(counts[, .(FinalSubtype, N)], dictionary,
                       by = "FinalSubtype", all.x = TRUE, sort = FALSE)
  plot_counts <- plot_counts[match(expected_states, FinalSubtype)]
  setcolorder(plot_counts, c("FinalSubtype", "DisplayLabel", "SubtypeClass", "N", "Color"))
  fwrite(plot_counts, file.path(viz_table_dir, "subtype_counts_for_figures.csv"))
}

study_dictionary <- fread(study_dictionary_path, showProgress = FALSE)
study_col <- intersect(c("Study", "study"), names(study_dictionary))
if (!length(study_col)) stop("Study color dictionary lacks Study column")
missing_studies <- setdiff(unique(master$Study), study_dictionary[[study_col[[1]]]])
if (length(missing_studies)) stop("Studies missing from study color dictionary: ", paste(missing_studies, collapse = ", "))

contract <- data.table(
  Item = c("Canonical assignment source", "Canonical subtype dictionary", "Study dictionary",
           "Master-table rows", "Primary key", "Main-state policy", "Intermediate-state policy",
           "Observed duct field", "Predicted duct field", "Mutation burden policy"),
  Specification = c(
    "1-Clustering/data/sample_info_with_clusters.xlsx",
    "2-TranscriptomicCharacterization/data/subtype_label_color_dictionary.csv",
    "0-Data/qc/study_color_dictionary.csv",
    as.character(nrow(master)), "Sample is unique",
    "Main-text inference uses SubtypeClass == Main",
    "Intermediate states are retained only as continuity/sensitivity evidence",
    "iCCA classification (observed pathology where available)",
    "lasso_pred or Step3 OOF/full projection fields; never merged with observed pathology",
    "Report nonsynonymous variant counts; no harmonized cross-study TMB"
  )
)
fwrite(contract, file.path(data_dir, "subtype_governance_contract.csv"))

report <- c(
  "# Subtype governance build report", "",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")), "",
  "The governance layer does not re-cluster samples or recompute biological statistics. It freezes the Step1 sample assignments, canonical labels/colors, and main/intermediate plotting policy for downstream joins.", "",
  "## Locked counts", "",
  paste(sprintf("- %s: %d", counts$FinalSubtype, counts$N), collapse = "\n"), "",
  "## Canonical naming", "",
  paste(sprintf("- %s: %s (%s)", dictionary$FinalSubtype, dictionary$DisplayLabel, dictionary$Color), collapse = "\n"), "",
  "C2 intentionally has no biological suffix. Observed and predicted duct variables remain separate."
)
writeLines(report, file.path(report_dir, "subtype_governance_build_report.md"))
cat("Subtype governance master table rebuilt:", nrow(master), "samples\n")
