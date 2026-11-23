#!/usr/bin/env Rscript

# Reproducible, study-aware oncoPredict refresh for Step5.
# Predicted response scores are hypothesis-generating; lower values indicate
# greater predicted sensitivity. They are not observed patient responses.

options(stringsAsFactors = FALSE, warn = 1)
root <- normalizePath(getwd(), mustWork = TRUE)
.libPaths(c(file.path(root, ".r-library"), .libPaths()))
suppressPackageStartupMessages({
  library(oncoPredict)
  library(sandwich)
  library(ggplot2)
})

module_dir <- file.path(root, "5-ClassifierValidation")
table_dir <- file.path(module_dir, "publication", "tables")
figure_dir <- file.path(module_dir, "publication", "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

load_single_object <- function(path, expected = NULL) {
  env <- new.env(parent = emptyenv())
  names <- load(path, envir = env)
  if (!is.null(expected) && expected %in% names) return(env[[expected]])
  if (length(names) != 1) stop("Expected one object in ", path, "; found: ", paste(names, collapse = ", "))
  env[[names[[1]]]]
}

bh <- function(p) p.adjust(p, method = "BH")
md5 <- function(path) unname(tools::md5sum(path))
truthy <- function(x) tolower(as.character(x)) %in% c("1", "true", "yes")

expr_path <- file.path(root, "1-Clustering", "data", "combined_expr_before_correction.rds")
master_path <- file.path(root, "governance", "data", "subtype_master_table.csv")
train_expr_path <- file.path(root, "oncoPredict_data", "Training Data", "GDSC2_Expr (RMA Normalized and Log Transformed).rds")
train_response_path <- file.path(root, "oncoPredict_data", "Training Data", "GDSC2_Res.rds")
stopifnot(file.exists(expr_path), file.exists(master_path), file.exists(train_expr_path), file.exists(train_response_path))

expr <- as.matrix(load_single_object(expr_path, "combined_expr"))
storage.mode(expr) <- "double"
master <- read.csv(master_path, check.names = FALSE)
master$Sample <- as.character(master$Sample)
master <- master[master$SubtypeClass == "Main", c("Sample", "FinalSubtype", "Study")]
master$FinalSubtype <- factor(master$FinalSubtype, levels = c("C1", "C2", "C3"))
common <- master$Sample[master$Sample %in% colnames(expr)]
if (length(common) != 497L || length(unique(common)) != 497L) stop("Expected 497 unique main-subtype samples; found ", length(unique(common)))
master <- master[match(common, master$Sample), ]
test_expr <- expr[, common, drop = FALSE]
if (any(!is.finite(test_expr))) stop("Test expression contains non-finite values")

training_expr <- as.matrix(readRDS(train_expr_path))
training_response <- as.matrix(readRDS(train_response_path))
shared_training <- intersect(colnames(training_expr), rownames(training_response))
training_expr <- training_expr[, shared_training, drop = FALSE]
training_response <- training_response[shared_training, , drop = FALSE]
shared_genes <- intersect(rownames(training_expr), rownames(test_expr))
if (length(shared_genes) < 5000L) stop("Too few genes shared by GDSC2 and iCCA: ", length(shared_genes))
# Training CV must use the same transferable gene universe as the tumor
# projection; otherwise GDSC2-only genes can overstate cross-platform quality.
training_expr <- training_expr[shared_genes, , drop = FALSE]
test_expr <- test_expr[shared_genes, , drop = FALSE]
if (ncol(training_expr) < 100 || ncol(training_response) < 50) stop("Unexpected GDSC2 training dimensions")

cores <- suppressWarnings(as.integer(Sys.getenv("SCLUCID_ONCOPREDICT_CORES", "4")))
if (!is.finite(cores) || cores < 1) cores <- 1L
set.seed(20260721)

quality_path <- file.path(table_dir, "step5_oncopredict_gdsc2_cv_quality.csv")
cv_manifest_path <- file.path(table_dir, "step5_oncopredict_gdsc2_cv_manifest.csv")
expected_cv_manifest <- data.frame(
  Field = c("training_expression_md5", "training_response_md5", "test_expression_md5",
            "n_shared_genes", "seed", "folds"),
  Value = c(md5(train_expr_path), md5(train_response_path), md5(expr_path),
            length(shared_genes), 20260721, 5), stringsAsFactors = FALSE
)
reuse_cv <- truthy(Sys.getenv("SCLUCID_ONCOPREDICT_REUSE_CV", "0")) && file.exists(quality_path)
if (reuse_cv) {
  if (!file.exists(cv_manifest_path)) stop("CV reuse requested but the CV provenance manifest is missing")
  observed_cv_manifest <- read.csv(cv_manifest_path, check.names = FALSE, colClasses = "character")
  expected_as_character <- expected_cv_manifest
  expected_as_character$Value <- as.character(expected_as_character$Value)
  if (!identical(observed_cv_manifest, expected_as_character)) stop("CV reuse requested but provenance does not match current inputs")
  message("Reusing completed shared-gene GDSC2 five-fold CV quality table")
  drug_quality <- read.csv(quality_path, check.names = FALSE)
} else {
  message("Running 5-fold GDSC2 cross-validation for drug-level reliability...")
  cv <- predictionAccuracybyCV(
    trainingExprData = training_expr,
    trainingPtype = training_response,
    testExprData = NULL,
    batchCorrect = "none",
    powerTransformPhenotype = TRUE,
    removeLowVaryingGenes = 0.2,
    minNumSamples = 10,
    selection = 1,
    printOutput = FALSE,
    pcr = FALSE,
    removeLowVaringGenesFrom = "homogenizeData",
    cvFold = 5,
    parallel = cores > 1,
    cores = cores
  )
  cv_pred <- as.matrix(cv$cvPtype)
  cv_real <- as.matrix(cv$realPtype)
  if (!identical(dim(cv_pred), dim(cv_real))) stop("CV prediction and truth dimensions differ")
  drug_quality <- do.call(rbind, lapply(seq_len(ncol(cv_real)), function(j) {
    ok <- is.finite(cv_real[, j]) & is.finite(cv_pred[, j])
    data.frame(
      Drug = colnames(cv_real)[j], N_training = sum(ok),
      Pearson_r = if (sum(ok) >= 10) cor(cv_real[ok, j], cv_pred[ok, j], method = "pearson") else NA_real_,
      Spearman_rho = if (sum(ok) >= 10) cor(cv_real[ok, j], cv_pred[ok, j], method = "spearman") else NA_real_,
      RMSE = if (sum(ok) >= 10) sqrt(mean((cv_real[ok, j] - cv_pred[ok, j])^2)) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  drug_quality$ReliabilityFlag <- with(drug_quality, N_training >= 50 & is.finite(Pearson_r) & Pearson_r >= 0.30)
  write.csv(drug_quality, quality_path, row.names = FALSE)
  write.csv(expected_cv_manifest, cv_manifest_path, row.names = FALSE)
}

# `standardize` transforms train and test independently gene by gene. The
# standardized GDSC2 training matrix is identical for every study, so fit the
# shared drug models once after retaining separate per-study test transforms.
study_hom <- lapply(unique(master$Study), function(study) {
  ids <- master$Sample[master$Study == study]
  message("Standardizing study ", study, " (n=", length(ids), ")")
  homogenizeData(test_expr[, ids, drop = FALSE], training_expr,
                 batchCorrect = "standardize", selection = 1, printOutput = FALSE)
})
training_standardized <- study_hom[[1]]$train
if (!all(vapply(study_hom[-1], function(x) isTRUE(all.equal(x$train, training_standardized, tolerance = 1e-12)), logical(1)))) {
  stop("Study-wise standardization unexpectedly changed the shared training matrix")
}
test_standardized <- do.call(cbind, lapply(study_hom, `[[`, "test"))
test_standardized <- test_standardized[, master$Sample, drop = FALSE]
message("Fitting one shared GDSC2 model set to the separately standardized studies...")
primary <- as.matrix(calcPhenotype(
  trainingExprData = training_standardized,
  trainingPtype = training_response,
  testExprData = test_standardized,
  batchCorrect = "none",
  powerTransformPhenotype = TRUE,
  removeLowVaryingGenes = 0.2,
  minNumSamples = 10,
  selection = 1,
  printOutput = FALSE,
  pcr = FALSE,
  removeLowVaringGenesFrom = "homogenizeData",
  report_pc = FALSE,
  cc = FALSE,
  rsq = FALSE,
  folder = FALSE,
  parallel = cores > 1,
  cores = cores
))
primary <- primary[master$Sample, , drop = FALSE]
if (nrow(primary) != 497L || any(!is.finite(primary))) stop("Primary oncoPredict matrix is incomplete")
write.csv(data.frame(Sample = rownames(primary), primary, check.names = FALSE),
          file.path(table_dir, "step5_oncopredict_gdsc2_studywise_predictions.csv"), row.names = FALSE)

# Remove each sample's global prediction tendency so that drug-specific signals
# can be checked alongside raw predictions.
sample_median <- apply(primary, 1, median, na.rm = TRUE)
centered <- primary - sample_median
write.csv(data.frame(Sample = rownames(centered), GlobalMedian = sample_median, centered, check.names = FALSE),
          file.path(table_dir, "step5_oncopredict_gdsc2_sample_centered_predictions.csv"), row.names = FALSE)

robust_fit_rows <- function(y, drug, scale_type) {
  dat <- data.frame(y = as.numeric(y), Subtype = master$FinalSubtype, Study = factor(master$Study))
  fit <- lm(y ~ Subtype + Study, data = dat)
  v <- vcovHC(fit, type = "HC3")
  b <- coef(fit)
  contrasts <- list("C2-C1" = c(SubtypeC2 = 1), "C3-C1" = c(SubtypeC3 = 1),
                    "C3-C2" = c(SubtypeC3 = 1, SubtypeC2 = -1))
  do.call(rbind, lapply(names(contrasts), function(label) {
    L <- rep(0, length(b)); names(L) <- names(b)
    terms <- contrasts[[label]]
    L[names(terms)] <- terms
    est <- sum(L * b); se <- sqrt(as.numeric(t(L) %*% v %*% L)); z <- est / se
    data.frame(Drug = drug, ScoreType = scale_type, Contrast = label, N = nrow(dat),
               Estimate = est, SE_HC3 = se, CI_low = est - 1.96 * se, CI_high = est + 1.96 * se,
               P = 2 * pnorm(abs(z), lower.tail = FALSE), stringsAsFactors = FALSE)
  }))
}

message("Fitting Study-adjusted subtype associations...")
effect_rows <- vector("list", ncol(primary) * 2L)
k <- 0L
for (j in seq_len(ncol(primary))) {
  k <- k + 1L; effect_rows[[k]] <- robust_fit_rows(primary[, j], colnames(primary)[j], "Raw")
  k <- k + 1L; effect_rows[[k]] <- robust_fit_rows(centered[, j], colnames(centered)[j], "SampleCentered")
}
effects <- do.call(rbind, effect_rows)
effects$FDR_within_score_contrast <- ave(effects$P, interaction(effects$ScoreType, effects$Contrast), FUN = bh)
effects$FDR_global <- bh(effects$P)
effects <- merge(effects, drug_quality, by = "Drug", all.x = TRUE, sort = FALSE)
effects$DirectionInterpretation <- ifelse(effects$Estimate < 0, "first subtype predicted more sensitive", "first subtype predicted less sensitive")
write.csv(effects, file.path(table_dir, "step5_oncopredict_study_adjusted_effects.csv"), row.names = FALSE)

study_direction <- do.call(rbind, lapply(colnames(primary), function(drug) {
  do.call(rbind, lapply(split(seq_len(nrow(master)), master$Study), function(ix) {
    c3 <- ix[master$FinalSubtype[ix] == "C3"]
    other <- ix[master$FinalSubtype[ix] != "C3"]
    if (min(length(c3), length(other)) < 2) return(NULL)
    est <- mean(centered[c3, drug]) - mean(centered[other, drug])
    data.frame(Drug = drug, Study = master$Study[ix[1]], N_C3 = length(c3), N_nonC3 = length(other),
               C3_minus_nonC3 = est, C3_more_sensitive = est < 0, stringsAsFactors = FALSE)
  }))
}))
write.csv(study_direction, file.path(table_dir, "step5_oncopredict_c3_study_direction.csv"), row.names = FALSE)
direction_summary <- aggregate(C3_more_sensitive ~ Drug, study_direction, function(x) c(Studies = length(x), Consistent = sum(x), Fraction = mean(x)))
direction_summary <- data.frame(Drug = direction_summary$Drug,
  Studies = direction_summary$C3_more_sensitive[, "Studies"],
  Studies_C3_more_sensitive = direction_summary$C3_more_sensitive[, "Consistent"],
  Fraction_C3_more_sensitive = direction_summary$C3_more_sensitive[, "Fraction"], check.names = FALSE)

candidate_pattern <- paste(c("trametinib", "selumetinib", "ulixertinib", "binimetinib", "cobimetinib", "dabrafenib",
                             "gemcitabine", "cisplatin", "oxaliplatin", "fluorouracil", "irinotecan", "paclitaxel"), collapse = "|")
candidates <- effects[effects$ScoreType == "SampleCentered" & effects$Contrast == "C3-C1" &
                        grepl(candidate_pattern, effects$Drug, ignore.case = TRUE), ]
candidates <- merge(candidates, direction_summary, by = "Drug", all.x = TRUE, sort = FALSE)
candidates <- candidates[order(candidates$FDR_within_score_contrast, candidates$Estimate), ]
write.csv(candidates, file.path(table_dir, "step5_oncopredict_prespecified_candidate_summary.csv"), row.names = FALSE)

run_sensitivity <- truthy(Sys.getenv("SCLUCID_ONCOPREDICT_RUN_SENSITIVITY", "1"))
sensitivity_summary <- data.frame()
if (run_sensitivity) {
  message("Running pooled standardize/ComBat sensitivity predictions...")
  pooled_standardize <- as.matrix(calcPhenotype(training_expr, training_response, test_expr,
    batchCorrect = "standardize", minNumSamples = 10, printOutput = FALSE,
    removeLowVaringGenesFrom = "homogenizeData", parallel = cores > 1, cores = cores))
  pooled_eb <- as.matrix(calcPhenotype(training_expr, training_response, test_expr,
    batchCorrect = "eb", minNumSamples = 10, printOutput = FALSE,
    removeLowVaringGenesFrom = "homogenizeData", parallel = cores > 1, cores = cores))
  pooled_standardize <- pooled_standardize[master$Sample, , drop = FALSE]
  pooled_eb <- pooled_eb[master$Sample, , drop = FALSE]
  sensitivity_summary <- do.call(rbind, lapply(colnames(primary), function(drug) data.frame(
    Drug = drug,
    Primary_vs_pooled_standardize_Spearman = cor(primary[, drug], pooled_standardize[, drug], method = "spearman"),
    Primary_vs_pooled_eb_Spearman = cor(primary[, drug], pooled_eb[, drug], method = "spearman"),
    stringsAsFactors = FALSE)))
  write.csv(sensitivity_summary, file.path(table_dir, "step5_oncopredict_batch_method_sensitivity.csv"), row.names = FALSE)
}

candidate_evidence <- merge(candidates, sensitivity_summary, by = "Drug", all.x = TRUE, sort = FALSE)
candidate_evidence <- candidate_evidence[order(candidate_evidence$FDR_within_score_contrast, candidate_evidence$Estimate), ]
write.csv(candidate_evidence, file.path(table_dir, "step5_oncopredict_candidate_evidence_matrix.csv"), row.names = FALSE)

plot_candidates <- candidates[is.finite(candidates$Estimate), ]
if (nrow(plot_candidates)) {
  plot_candidates$Drug <- factor(plot_candidates$Drug, levels = rev(unique(plot_candidates$Drug)))
  p <- ggplot(plot_candidates, aes(Estimate, Drug, color = ReliabilityFlag)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey60") +
    geom_errorbar(aes(xmin = CI_low, xmax = CI_high), orientation = "y", width = 0) +
    geom_point(size = 2.4) +
    scale_color_manual(values = c(`TRUE` = "#0072B2", `FALSE` = "#999999"), na.value = "#999999") +
    labs(x = "Adjusted C3 - C1 predicted response score\n(lower = greater predicted sensitivity)", y = NULL,
         color = "GDSC2 CV\nPearson r >= 0.30") + theme_bw(base_size = 10)
  ggsave(file.path(figure_dir, "step5_oncopredict_c3_candidates.pdf"), p, width = 7, height = max(4, .28 * nrow(plot_candidates) + 2))
  ggsave(file.path(figure_dir, "step5_oncopredict_c3_candidates.png"), p, width = 7, height = max(4, .28 * nrow(plot_candidates) + 2), dpi = 300)
}

manifest <- data.frame(
  Field = c("analysis", "oncoPredict_version", "R_version", "training_expression", "training_response",
            "training_expression_md5", "training_response_md5", "test_expression_md5", "master_table_md5",
            "n_samples", "n_studies", "n_shared_genes", "n_drugs", "primary_batch_strategy", "training_quality",
            "sensitivity_analysis", "interpretation"),
  Value = c("Step5 study-aware oncoPredict refresh", as.character(packageVersion("oncoPredict")), R.version.string,
            sub(paste0("^", root, "/"), "", train_expr_path), sub(paste0("^", root, "/"), "", train_response_path),
            md5(train_expr_path), md5(train_response_path), md5(expr_path), md5(master_path),
            nrow(primary), length(unique(master$Study)), length(shared_genes), ncol(primary),
            "GDSC2 microarray training; each iCCA Study standardized separately, then one identical shared GDSC2 drug-model fit",
            "seeded 5-fold CV in GDSC2; per-drug Pearson/Spearman/RMSE; reliability flag Pearson r>=0.30 and n>=50",
            ifelse(run_sensitivity, "pooled standardize and pooled empirical-Bayes/ComBat", "not requested"),
            "hypothesis-generating predicted response; lower score means predicted sensitivity; not patient treatment response"),
  stringsAsFactors = FALSE
)
write.csv(manifest, file.path(table_dir, "step5_oncopredict_run_manifest.csv"), row.names = FALSE)

gates <- data.frame(
  Check = c("main_sample_count_497", "exact_unique_sample_ids", "all_three_main_subtypes", "all_seven_studies",
            "shared_gene_universe_ge_5000", "all_gdsc2_drugs_predicted", "finite_predictions", "gdsc2_5fold_cv_complete",
            "study_adjusted_effects_written", "cv_cache_provenance_matches", "candidate_evidence_matrix_written", "score_direction_explicit", "hypothesis_generating_language_explicit"),
  Pass = c(nrow(primary) == 497L, length(unique(rownames(primary))) == 497L,
           setequal(as.character(unique(master$FinalSubtype)), c("C1", "C2", "C3")), length(unique(master$Study)) == 7L,
           length(shared_genes) >= 5000L, ncol(primary) == ncol(training_response), all(is.finite(primary)), nrow(drug_quality) == ncol(primary),
           nrow(effects) == ncol(primary) * 2L * 3L, file.exists(cv_manifest_path), nrow(candidate_evidence) == nrow(candidates), TRUE, TRUE),
  Detail = c(nrow(primary), length(unique(rownames(primary))), paste(levels(droplevels(master$FinalSubtype)), collapse = ","),
             length(unique(master$Study)), length(shared_genes), ncol(primary), sum(!is.finite(primary)), nrow(drug_quality), nrow(effects), cv_manifest_path, nrow(candidate_evidence),
             "lower=response score means greater predicted sensitivity", "explicit in manifest and notebook"),
  stringsAsFactors = FALSE
)
write.csv(gates, file.path(table_dir, "step5_oncopredict_release_gate_status.csv"), row.names = FALSE)
if (!all(gates$Pass)) stop("One or more Step5 oncoPredict release checks failed")
cat(sprintf('{"status":"PASS","samples":%d,"studies":%d,"drugs":%d,"tables":"%s"}\n',
            nrow(primary), length(unique(master$Study)), ncol(primary), table_dir))
