#!/usr/bin/env Rscript

# Publication-facing transcription-factor activity analysis for Step 1.
# Reads the locked corrected expression interface and writes only to
# 1-Clustering/publication; no Step 0 or historical Step 1 object is modified.

suppressPackageStartupMessages({
  library(data.table)
  library(decoupleR)
  library(dorothea)
  library(limma)
  library(readxl)
})

set.seed(123)
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
table_dir <- file.path(root, "1-Clustering", "publication", "tables")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

load_rdata_object <- function(path, object) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  if (!exists(object, envir = env, inherits = FALSE)) {
    stop("Missing object ", object, " in ", path)
  }
  get(object, envir = env, inherits = FALSE)
}

sample_info <- as.data.table(read_excel(
  file.path(root, "1-Clustering", "data", "sample_info_with_clusters.xlsx")
))
expression <- load_rdata_object(
  file.path(root, "1-Clustering", "data", "combined_expr_after_correction.rds"),
  "corrected_expr_limma"
)
main_subtypes <- c("C1", "C2", "C3")
main_info <- sample_info[Ensemble %in% main_subtypes, .(
  Sample = as.character(Sample), Subtype = Ensemble, Study = as.character(Study)
)]
expression <- expression[, main_info$Sample, drop = FALSE]
stopifnot(nrow(main_info) == 497L, ncol(expression) == 497L, !anyNA(expression))

# High- and intermediate-confidence regulons only. Duplicate TF-target edges are
# collapsed deterministically, retaining the mean mode of regulation.
network <- as.data.table(dorothea_hs)[confidence %in% c("A", "B", "C"), .(
  mor = mean(mor)
), by = .(source = tf, target)]
network <- network[target %in% rownames(expression)]
regulon_coverage <- network[, .(TargetsAvailable = uniqueN(target)), by = source][order(-TargetsAvailable)]
fwrite(regulon_coverage, file.path(table_dir, "step1_dorothea_regulon_coverage.csv"))

ulm <- as.data.table(run_ulm(
  mat = expression,
  network = network,
  .source = source,
  .target = target,
  .mor = mor,
  minsize = 5L
))
stopifnot(all(c("source", "condition", "score", "p_value") %in% names(ulm)))
setnames(ulm, c("source", "condition", "score", "p_value"), c("TF", "Sample", "Activity", "ActivityPValue"))
ulm <- merge(ulm, main_info, by = "Sample", all.x = TRUE, sort = FALSE)
stopifnot(nrow(ulm) == uniqueN(ulm$TF) * 497L, !anyNA(ulm$Subtype), all(is.finite(ulm$Activity)))
fwrite(ulm, file.path(table_dir, "step1_dorothea_ulm_activity_long.csv"))

activity_matrix <- dcast(ulm, TF ~ Sample, value.var = "Activity")
tf_names <- activity_matrix$TF
activity <- as.matrix(activity_matrix[, -"TF"])
rownames(activity) <- tf_names

fit_tf_contrast <- function(subtype, metadata = main_info, include_study = TRUE) {
  metadata <- copy(metadata)
  metadata[, Group := factor(ifelse(Subtype == subtype, subtype, "Others"), levels = c("Others", subtype))]
  metadata[, Study := factor(Study)]
  x <- activity[, metadata$Sample, drop = FALSE]
  design <- if (include_study && uniqueN(metadata$Study) > 1L) {
    model.matrix(~ 0 + Group + Study, data = metadata)
  } else {
    model.matrix(~ 0 + Group, data = metadata)
  }
  subtype_column <- grep(paste0("^Group", subtype, "$"), colnames(design), value = TRUE)
  others_column <- grep("^GroupOthers$", colnames(design), value = TRUE)
  stopifnot(length(subtype_column) == 1L, length(others_column) == 1L)
  contrast <- matrix(0, nrow = ncol(design), ncol = 1L,
                     dimnames = list(colnames(design), paste0(subtype, "_minus_others")))
  contrast[subtype_column, 1L] <- 1
  contrast[others_column, 1L] <- -1
  fit <- eBayes(contrasts.fit(lmFit(x, design), contrast), robust = TRUE)
  result <- as.data.table(topTable(fit, number = Inf, sort.by = "P"), keep.rownames = "TF")
  setnames(result, c("logFC", "P.Value", "adj.P.Val"), c(
    "ActivityDifferenceSubtypeMinusOthers", "PValue", "FDR"
  ))
  result[, `:=`(
    Subtype = subtype,
    Contrast = paste0(subtype, "_minus_others"),
    PositiveDirection = paste0("Higher inferred TF activity in ", subtype),
    NSubtype = metadata[Subtype == subtype, .N],
    NOthers = metadata[Subtype != subtype, .N],
    NTotal = nrow(metadata),
    StudyAdjusted = include_study && uniqueN(metadata$Study) > 1L
  )]
  result
}

tf_differential <- rbindlist(lapply(main_subtypes, fit_tf_contrast))
setcolorder(tf_differential, c(
  "TF", "Subtype", "Contrast", "PositiveDirection", "ActivityDifferenceSubtypeMinusOthers",
  "AveExpr", "t", "PValue", "FDR", "B", "NSubtype", "NOthers", "NTotal", "StudyAdjusted"
))
fwrite(tf_differential, file.path(table_dir, "step1_dorothea_tf_activity_main_only_study_adjusted.csv"))

top_tf <- tf_differential[order(Subtype, FDR, -abs(ActivityDifferenceSubtypeMinusOthers)), head(.SD, 20L), by = Subtype]
fwrite(top_tf, file.path(table_dir, "step1_dorothea_top_tf_by_subtype.csv"))

# Leave-one-study-out sensitivity for TFs significant in the complete fit.
headline_tf <- tf_differential[FDR < 0.05, unique(TF)]
loso <- rbindlist(lapply(sort(unique(main_info$Study)), function(held_out) {
  reduced_info <- main_info[Study != held_out]
  rbindlist(lapply(main_subtypes, function(subtype) {
    reduced <- fit_tf_contrast(subtype, reduced_info, include_study = TRUE)
    full <- tf_differential[Subtype == subtype & TF %in% headline_tf]
    joined <- merge(
      full[, .(TF, FullDifference = ActivityDifferenceSubtypeMinusOthers, FullFDR = FDR)],
      reduced[, .(TF, ReducedDifference = ActivityDifferenceSubtypeMinusOthers, ReducedFDR = FDR)],
      by = "TF"
    )
    joined[, `:=`(
      Subtype = subtype,
      HeldOutStudy = held_out,
      SignAgreement = sign(FullDifference) == sign(ReducedDifference)
    )]
    joined
  }))
}))
fwrite(loso, file.path(table_dir, "step1_dorothea_tf_leave_one_study_out.csv"))

loso_summary <- loso[, .(
  TFs = .N,
  SignAgreementRate = mean(SignAgreement),
  PearsonR = cor(FullDifference, ReducedDifference),
  SignificantAgainAtFDR05 = mean(ReducedFDR < 0.05)
), by = .(Subtype, HeldOutStudy)]
fwrite(loso_summary, file.path(table_dir, "step1_dorothea_tf_leave_one_study_out_summary.csv"))

method <- data.table(
  Field = c("Expression input", "Sample scope", "Regulon", "Activity method", "Contrast", "Multiplicity", "Sensitivity"),
  Value = c(
    "Locked limma-corrected 10,531-gene expression interface",
    "497 main-state samples only (C1/C2/C3)",
    paste0("DoRothEA human confidence A-C; ", uniqueN(network$source), " TFs with >=5 available targets"),
    "decoupleR univariate linear model (ULM)",
    "Explicit subtype minus others with Study covariate",
    "Benjamini-Hochberg FDR within each subtype-versus-others contrast",
    "Leave-one-study-out refitting for full-fit significant TFs"
  )
)
fwrite(method, file.path(table_dir, "step1_dorothea_tf_method_contract.csv"))

cat("Step 1 TF activity completed: ", uniqueN(ulm$TF), " TFs x 497 samples\n", sep = "")
