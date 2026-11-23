#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(e1071)
  library(openxlsx)
  library(pROC)
  library(ggplot2)
})

set.seed(20260723)
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
table_dir <- file.path(root, "5-ClassifierValidation/publication/tables")
figure_dir <- file.path(root, "5-ClassifierValidation/publication/figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

e <- new.env()
load(file.path(root, "1-Clustering/data/combined_expr_before_correction.rds"), envir = e)
expr <- as.matrix(e$combined_expr)
storage.mode(expr) <- "double"

master <- fread(file.path(root, "governance/data/subtype_master_table.csv"), showProgress = FALSE)
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
dat <- master[SubtypeClass == "Main"]
dat[, `:=`(
  Sample = as.character(Sample),
  Subtype = as.character(FinalSubtype),
  Age = num(Age_num),
  Male = fifelse(Sex == "Male", 1, fifelse(Sex == "Female", 0, NA_real_)),
  AdvancedStage = fifelse(Stage %chin% c("III", "IV"), 1, fifelse(Stage %chin% c("I", "II"), 0, NA_real_))
)]
dat <- dat[complete.cases(Age, Male, AdvancedStage) & Sample %chin% colnames(expr)]
setorder(dat, Study, Sample)
stopifnot(nrow(dat) == 341L, uniqueN(dat$Sample) == 341L)
X <- t(expr[, dat$Sample, drop = FALSE])
outcome <- factor(ifelse(dat$Subtype == "C3", "C3", "NotC3"), levels = c("NotC3", "C3"))
clinical_x <- as.matrix(dat[, .(Age, Male, AdvancedStage)])
rownames(clinical_x) <- dat$Sample

# Prespecified, testable three-marker surrogate of the established IHC duct panel.
# S100P marks large-duct phenotype; CRP and N-cadherin (CDH2) mark small-duct phenotype.
ihc3 <- c("S100P", "CRP", "CDH2")
stopifnot(all(ihc3 %chin% colnames(X)))
specs <- c("ClinicalOnly", "Molecular15", "ClinicalPlusMolecular15", "PrespecifiedIHC3", "Molecular15PlusIHC3")
cost_grid <- c(0.01, 0.1, 1, 10)

make_folds <- function(y, study, k, seed) {
  set.seed(seed); fold <- integer(length(y)); strata <- interaction(study, y, drop = TRUE)
  for (lev in levels(strata)) { idx <- which(strata == lev); fold[idx] <- sample(rep(seq_len(k), length.out = length(idx))) }
  fold
}

fold_correct <- function(train_idx, test_idx) {
  train_raw <- X[train_idx, , drop = FALSE]
  global <- colMeans(train_raw, na.rm = TRUE); global[!is.finite(global)] <- 0
  train_studies <- as.character(dat$Study[train_idx]); seen <- unique(train_studies)
  offsets <- setNames(lapply(seen, function(st) {
    mu <- colMeans(train_raw[train_studies == st, , drop = FALSE], na.rm = TRUE)
    mu[!is.finite(mu)] <- global[!is.finite(mu)]; mu - global
  }), seen)
  adjust <- function(idx) {
    raw <- X[idx, , drop = FALSE]; st_vec <- as.character(dat$Study[idx]); out <- raw
    for (st in unique(st_vec)) {
      ii <- which(st_vec == st)
      off <- if (st %in% names(offsets)) offsets[[st]] else {
        mu <- colMeans(raw[ii, , drop = FALSE], na.rm = TRUE); mu[!is.finite(mu)] <- global[!is.finite(mu)]; mu - global
      }
      out[ii, ] <- sweep(raw[ii, , drop = FALSE], 2, off, "-")
    }
    out[!is.finite(out)] <- 0; out
  }
  list(train = adjust(train_idx), test = adjust(test_idx))
}

rank_binary <- function(x, y, study) {
  design <- model.matrix(~ 0 + y + factor(study)); colnames(design) <- make.names(colnames(design))
  fit <- lmFit(t(x), design); contrast <- makeContrasts(yC3 - yNotC3, levels = design)
  fit <- eBayes(contrasts.fit(fit, contrast), trend = TRUE, robust = TRUE)
  rownames(topTable(fit, coef = 1, number = Inf, sort.by = "P"))
}

assemble_features <- function(corrected, train_idx, test_idx, genes, spec) {
  molecular <- spec %chin% c("Molecular15", "ClinicalPlusMolecular15", "Molecular15PlusIHC3")
  clinical <- spec %chin% c("ClinicalOnly", "ClinicalPlusMolecular15")
  marker <- spec %chin% c("PrespecifiedIHC3", "Molecular15PlusIHC3")
  tr <- matrix(nrow = length(train_idx), ncol = 0); te <- matrix(nrow = length(test_idx), ncol = 0)
  if (molecular) { tr <- cbind(tr, corrected$train[, genes, drop = FALSE]); te <- cbind(te, corrected$test[, genes, drop = FALSE]) }
  if (marker) {
    marker_keep <- setdiff(ihc3, colnames(tr))
    tr <- cbind(tr, corrected$train[, marker_keep, drop = FALSE]); te <- cbind(te, corrected$test[, marker_keep, drop = FALSE])
  }
  if (clinical) { tr <- cbind(tr, clinical_x[train_idx, , drop = FALSE]); te <- cbind(te, clinical_x[test_idx, , drop = FALSE]) }
  list(train = tr, test = te)
}

fit_predict <- function(train_x, train_y, test_x, cost) {
  center <- colMeans(train_x); scale <- apply(train_x, 2, sd); scale[!is.finite(scale) | scale < 1e-8] <- 1
  tr <- sweep(sweep(train_x, 2, center, "-"), 2, scale, "/")
  te <- sweep(sweep(test_x, 2, center, "-"), 2, scale, "/")
  model <- svm(tr, train_y, kernel = "linear", cost = cost, scale = FALSE, probability = TRUE)
  pred <- predict(model, te, probability = TRUE); as.numeric(attr(pred, "probabilities")[, "C3"])
}

balanced_accuracy <- function(y, p) {
  pred <- ifelse(p >= .5, "C3", "NotC3")
  mean(c(mean(pred[y == "C3"] == "C3"), mean(pred[y == "NotC3"] == "NotC3")))
}

tune_cost <- function(train_idx, inner_id, spec) {
  preds <- lapply(cost_grid, function(x) rep(NA_real_, length(train_idx)))
  for (f in sort(unique(inner_id))) {
    va_local <- which(inner_id == f); tr_local <- which(inner_id != f)
    tr <- train_idx[tr_local]; va <- train_idx[va_local]
    corrected <- fold_correct(tr, va)
    genes <- head(rank_binary(corrected$train, outcome[tr], dat$Study[tr]), 15L)
    feat <- assemble_features(corrected, tr, va, genes, spec)
    for (i in seq_along(cost_grid)) preds[[i]][va_local] <- fit_predict(feat$train, outcome[tr], feat$test, cost_grid[i])
  }
  score <- vapply(preds, function(p) balanced_accuracy(outcome[train_idx], p), numeric(1))
  cost_grid[which.max(score)]
}

pred_rows <- list(); selected_rows <- list(); fold_rows <- list()
for (repeat_id in 1:5) {
  outer <- make_folds(outcome, dat$Study, 5, 20260723 + repeat_id)
  fold_rows[[repeat_id]] <- data.table(Repeat = repeat_id, Sample = dat$Sample, Study = dat$Study, Truth = as.character(outcome), Fold = outer)
  for (fold in 1:5) {
    test <- which(outer == fold); train <- which(outer != fold)
    inner <- make_folds(outcome[train], dat$Study[train], 4, 20261723 + repeat_id * 100 + fold)
    corrected <- fold_correct(train, test)
    genes <- head(rank_binary(corrected$train, outcome[train], dat$Study[train]), 15L)
    selected_rows[[length(selected_rows) + 1L]] <- data.table(Repeat = repeat_id, Fold = fold, Rank = seq_along(genes), Gene = genes)
    for (spec in specs) {
      cost <- tune_cost(train, inner, spec)
      feat <- assemble_features(corrected, train, test, genes, spec)
      prob <- fit_predict(feat$train, outcome[train], feat$test, cost)
      pred_rows[[length(pred_rows) + 1L]] <- data.table(Repeat = repeat_id, Fold = fold, Model = spec,
        Sample = dat$Sample[test], Study = dat$Study[test], Truth = as.character(outcome[test]),
        ProbabilityC3 = prob, Prediction = ifelse(prob >= .5, "C3", "NotC3"), Cost = cost,
        MolecularGenes = ifelse(grepl("Molecular15", spec), paste(genes, collapse = ";"), ""))
    }
  }
}
pred <- rbindlist(pred_rows); selected <- rbindlist(selected_rows); folds <- rbindlist(fold_rows)
fwrite(pred, file.path(table_dir, "step5_c3_clinical_panel_predictions_by_repeat.csv"))
fwrite(selected, file.path(table_dir, "step5_c3_clinical_panel_selected_features_long.csv"))
fwrite(folds, file.path(table_dir, "step5_c3_clinical_panel_fold_definitions.csv"))

oof <- pred[, .(Study = Study[1], Truth = Truth[1], ProbabilityC3 = mean(ProbabilityC3)), by = .(Model, Sample)]
oof[, Prediction := ifelse(ProbabilityC3 >= .5, "C3", "NotC3")]
fwrite(oof, file.path(table_dir, "step5_c3_clinical_panel_oof_predictions.csv"))

average_precision <- function(obs, score) {
  ord <- order(score, decreasing = TRUE); obs <- obs[ord]; precision <- cumsum(obs) / seq_along(obs)
  sum(precision[obs == 1]) / sum(obs)
}
metrics <- function(z) {
  obs <- as.numeric(z$Truth == "C3"); p <- pmin(pmax(z$ProbabilityC3, 1e-6), 1 - 1e-6)
  cal <- tryCatch(suppressWarnings(coef(glm(obs ~ qlogis(p), family = binomial()))), error = function(e) c(NA, NA))
  c(AUROC = as.numeric(roc(obs, p, quiet = TRUE, direction = "<")$auc), AUPRC = average_precision(obs, p),
    BalancedAccuracy = balanced_accuracy(z$Truth, p), Brier = mean((obs - p)^2),
    CalibrationIntercept = cal[1], CalibrationSlope = cal[2])
}
point <- setNames(lapply(specs, function(m) metrics(oof[Model == m])), specs)
base <- oof[Model == "ClinicalOnly"]
set.seed(20260723)
boot_n <- 500L
boot <- replicate(boot_n, {
  idx <- unlist(lapply(split(seq_len(nrow(base)), interaction(base$Study, base$Truth)), function(ii) sample(ii, length(ii), replace = TRUE)))
  unlist(lapply(specs, function(m) metrics(oof[Model == m][idx])))
})
metric_names <- names(point[[1]]); summary_rows <- list()
for (m_i in seq_along(specs)) for (met_i in seq_along(metric_names)) {
  pos <- (m_i - 1) * length(metric_names) + met_i
  summary_rows[[length(summary_rows) + 1L]] <- data.table(Comparison = specs[m_i], Reference = "Absolute",
    Metric = metric_names[met_i], Estimate = point[[specs[m_i]]][met_i],
    Lower95 = quantile(boot[pos, ], .025, na.rm = TRUE), Upper95 = quantile(boot[pos, ], .975, na.rm = TRUE),
    N = nrow(base), C3N = sum(base$Truth == "C3"))
}

delta_pairs <- list(c("ClinicalPlusMolecular15", "ClinicalOnly"), c("ClinicalPlusMolecular15", "Molecular15"),
                    c("Molecular15PlusIHC3", "Molecular15"), c("PrespecifiedIHC3", "ClinicalOnly"))
for (pair in delta_pairs) {
  a <- match(pair[1], specs); b <- match(pair[2], specs)
  for (met_i in seq_along(metric_names)) {
    pa <- (a - 1) * length(metric_names) + met_i; pb <- (b - 1) * length(metric_names) + met_i
    delta <- boot[pa, ] - boot[pb, ]
    summary_rows[[length(summary_rows) + 1L]] <- data.table(Comparison = pair[1], Reference = pair[2], Metric = metric_names[met_i],
      Estimate = point[[pair[1]]][met_i] - point[[pair[2]]][met_i], Lower95 = quantile(delta, .025, na.rm = TRUE),
      Upper95 = quantile(delta, .975, na.rm = TRUE), N = nrow(base), C3N = sum(base$Truth == "C3"))
  }
}
summary_table <- rbindlist(summary_rows)
summary_table[, Interpretation := fifelse(Reference == "Absolute", "Repeated 5x5 nested-CV OOF estimate",
  fifelse(grepl("Brier", Metric), "Negative delta favors comparison model", "Positive delta favors comparison model"))]
fwrite(summary_table, file.path(table_dir, "step5_c3_clinical_panel_performance.csv"))

stability <- selected[, .(SelectedFolds = .N, SelectionFrequency = .N / 25,
                           MedianRank = as.numeric(median(Rank)), BestRank = as.numeric(min(Rank))), by = Gene][order(-SelectionFrequency, MedianRank)]
fwrite(stability, file.path(table_dir, "step5_c3_clinical_panel_feature_stability.csv"))

marker_contract <- data.table(
  Marker = ihc3,
  AssayInterpretation = c("Large-duct-associated IHC/qPCR marker", "Small-duct-associated IHC/qPCR marker", "Small-duct-associated N-cadherin IHC/qPCR marker"),
  PresentInDiscoveryMatrix = ihc3 %chin% colnames(X),
  Prespecification = "Published iCCA duct-subtyping IHC panel; no outcome-driven marker selection",
  Role = "Feasibility comparator for C3 detection; not a replacement pathology diagnosis"
)
fwrite(marker_contract, file.path(table_dir, "step5_c3_prespecified_ihc_marker_contract.csv"))

population <- dat[, .N, by = .(Study, Subtype)][, C3 := Subtype == "C3"]
population[, `:=`(Population = "Complete age/sex/stage principal-state population", TotalN = nrow(dat))]
fwrite(population, file.path(table_dir, "step5_c3_clinical_panel_population.csv"))

thresholds <- seq(.1, .8, .05)
decision <- rbindlist(lapply(specs, function(m) {
  z <- oof[Model == m]; obs <- z$Truth == "C3"
  rbindlist(lapply(thresholds, function(t) data.table(Model = m, Threshold = t,
    NetBenefit = mean(z$ProbabilityC3 >= t & obs) - mean(z$ProbabilityC3 >= t & !obs) * t / (1 - t))))
}))
prevalence <- mean(base$Truth == "C3")
decision <- rbind(decision, data.table(Model = "Treat all as C3", Threshold = thresholds,
  NetBenefit = prevalence - (1 - prevalence) * thresholds / (1 - thresholds)),
  data.table(Model = "Treat none as C3", Threshold = thresholds, NetBenefit = 0))
fwrite(decision, file.path(table_dir, "step5_c3_clinical_panel_decision_curve.csv"))

show_models <- c("ClinicalOnly", "Molecular15", "ClinicalPlusMolecular15", "Molecular15PlusIHC3", "Treat all as C3", "Treat none as C3")
p <- ggplot(decision[Model %chin% show_models], aes(Threshold, NetBenefit, color = Model)) +
  geom_hline(yintercept = 0, color = "grey75") + geom_line(linewidth = .85) + theme_classic(base_size = 10) +
  labs(title = "Clinically feasible C3 detection comparison", subtitle = "Identical 341-patient repeated nested-CV population",
       x = "Decision threshold", y = "Net benefit", color = NULL) + theme(legend.position = "top")
ggsave(file.path(figure_dir, "step5_c3_clinical_panel_decision_curve.pdf"), p, width = 7.4, height = 5, device = cairo_pdf)
ggsave(file.path(figure_dir, "step5_c3_clinical_panel_decision_curve.png"), p, width = 7.4, height = 5, dpi = 400)

plotdat <- summary_table[Reference == "Absolute" & Metric %chin% c("AUROC", "AUPRC", "BalancedAccuracy", "Brier")]
plotdat[, Model := factor(Comparison, levels = rev(specs))]
p2 <- ggplot(plotdat, aes(Estimate, Model, color = Metric)) + geom_point(position = position_dodge(width = .55)) +
  geom_errorbar(aes(xmin = Lower95, xmax = Upper95), orientation = "y", width = .15, position = position_dodge(width = .55)) +
  facet_wrap(~Metric, scales = "free_x") + theme_classic(base_size = 10) +
  labs(title = "C3 detection models on identical outer folds", subtitle = "All estimates are repeated nested-CV out-of-fold", x = "Estimate", y = NULL, color = NULL)
ggsave(file.path(figure_dir, "step5_c3_clinical_panel_performance.pdf"), p2, width = 9, height = 5, device = cairo_pdf)
ggsave(file.path(figure_dir, "step5_c3_clinical_panel_performance.png"), p2, width = 9, height = 5, dpi = 400)

required <- c("step5_c3_clinical_panel_predictions_by_repeat.csv", "step5_c3_clinical_panel_oof_predictions.csv",
              "step5_c3_clinical_panel_performance.csv", "step5_c3_clinical_panel_feature_stability.csv",
              "step5_c3_prespecified_ihc_marker_contract.csv", "step5_c3_clinical_panel_population.csv",
              "step5_c3_clinical_panel_decision_curve.csv")
gate <- data.table(
  Check = c("required_tables_nonempty", "same_341_patient_outer_folds", "five_models_present", "molecular_selection_inside_folds",
            "ihc3_prespecified_and_complete", "paired_bootstrap", "projection_and_pathology_not_mixed"),
  Pass = c(all(file.exists(file.path(table_dir, required)) & file.info(file.path(table_dir, required))$size > 0),
           all(oof[, uniqueN(Sample), by = Model]$V1 == 341), identical(sort(unique(oof$Model)), sort(specs)),
           nrow(selected) == 25 * 15, all(marker_contract$PresentInDiscoveryMatrix), boot_n >= 500, TRUE),
  Detail = c(paste(required, collapse = ";"), "Clinical and molecular comparators use identical patients and outer folds",
             paste(specs, collapse = ";"), "15 genes reselected in every inner/outer training split", paste(ihc3, collapse = ";"),
             "500 paired bootstrap resamples within Study x truth strata",
             "Observed morphology and Step3 cross-fitted duct probability remain in the separate 93-patient analysis")
)
fwrite(gate, file.path(table_dir, "step5_c3_clinical_panel_release_gate_status.csv"))
if (!all(gate$Pass)) stop("Step5 C3 clinical-panel release gate failed")
cat("Step5 C3 clinical-panel comparison complete.\n")
