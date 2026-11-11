#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(limma)
  library(e1071)
  library(pROC)
  library(openxlsx)
  library(ggplot2)
})

set.seed(20260715)
classes <- c("C1", "C2", "C3")
panel_sizes <- c(Operational123 = 123L, Nested31 = 31L, Compact15 = 15L)
cost_grid <- c(0.01, 0.1, 1, 10)
outer_repeats <- 5L
outer_folds <- 5L
inner_folds <- 4L
bootstrap_reps <- 1000L

root <- "5-ClassifierValidation/publication"
table_dir <- file.path(root, "tables")
figure_dir <- file.path(root, "figures")
model_dir <- file.path(root, "models")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)

subtype_colors <- read.csv("docs/subtype_color_dictionary.csv", check.names = FALSE)
subtype_colors <- setNames(subtype_colors$ColorHex[subtype_colors$FinalSubtype %in% classes],
                           subtype_colors$FinalSubtype[subtype_colors$FinalSubtype %in% classes])

load_saved_object <- function(path, object_name = NULL) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  nms <- ls(e)
  if (!is.null(object_name) && object_name %in% nms) return(get(object_name, e))
  if (length(nms) != 1L) stop("Expected one object in ", path, "; found: ", paste(nms, collapse = ", "))
  get(nms, e)
}

expr <- load_saved_object("1-Clustering/data/combined_expr_before_correction.rds", "combined_expr")
meta <- read.xlsx("1-Clustering/data/sample_info_with_clusters.xlsx")
meta <- meta[meta$Ensemble %in% classes, c("Sample", "Study", "Ensemble")]
colnames(meta) <- c("Sample", "Study", "Subtype")
meta$Sample <- as.character(meta$Sample)
meta$Study <- as.character(meta$Study)
meta$Subtype <- factor(meta$Subtype, levels = classes)

common_samples <- intersect(meta$Sample, colnames(expr))
meta <- meta[match(common_samples, meta$Sample), ]
expr <- expr[, common_samples, drop = FALSE]
stopifnot(identical(meta$Sample, colnames(expr)), nrow(meta) == 497L)
X <- t(as.matrix(expr))
storage.mode(X) <- "double"
rm(expr)

# Fold-safe study correction.  Study offsets are estimated using training
# samples only.  An unseen held-out study is centered with its own unlabeled
# cohort mean (transductive transfer sensitivity), never with subtype labels.
fold_correct <- function(train_idx, test_idx) {
  train_raw <- X[train_idx, , drop = FALSE]
  global <- colMeans(train_raw, na.rm = TRUE)
  global[!is.finite(global)] <- 0
  train_studies <- as.character(meta$Study[train_idx])
  seen <- unique(train_studies)
  offsets <- setNames(lapply(seen, function(st) {
    mu <- colMeans(train_raw[train_studies == st, , drop = FALSE], na.rm = TRUE)
    mu[!is.finite(mu)] <- global[!is.finite(mu)]
    mu - global
  }), seen)
  adjust <- function(idx) {
    raw <- X[idx, , drop = FALSE]
    st_vec <- as.character(meta$Study[idx])
    out <- raw
    for (st in unique(st_vec)) {
      ii <- which(st_vec == st)
      off <- if (st %in% names(offsets)) offsets[[st]] else {
        # Held-out study: unlabeled cohort centering only.
        mu <- colMeans(raw[ii, , drop = FALSE], na.rm = TRUE)
        mu[!is.finite(mu)] <- global[!is.finite(mu)]
        mu - global
      }
      out[ii, ] <- sweep(raw[ii, , drop = FALSE], 2L, off, "-")
    }
    out[!is.finite(out)] <- 0
    out
  }
  list(train = adjust(train_idx), test = adjust(test_idx),
       policy = "training-study offsets; held-out study unlabeled cohort centering")
}

legacy_env <- new.env()
load("archive/notebook_backups/5-ClassifierValidation/selected_signature_genes.rds", envir = legacy_env)
legacy_123 <- get("genes_in_at_least_two", legacy_env)
legacy_31 <- get("intersect_all_three", legacy_env)

write.csv(data.frame(Gene = legacy_31, Role = "Historical triple-overlap interpretive core",
                     ModelUse = "Not used as the fixed feature set for headline validation"),
          file.path(table_dir, "step5_historical_31_gene_interpretive_core.csv"), row.names = FALSE)

make_folds <- function(y, study, k, seed) {
  set.seed(seed)
  strata <- interaction(study, y, drop = TRUE)
  fold <- integer(length(y))
  for (lev in levels(strata)) {
    idx <- which(strata == lev)
    fold[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }
  fold
}

rank_features <- function(x, y, study) {
  y <- droplevels(factor(y, levels = classes))
  study <- droplevels(factor(study))
  design <- model.matrix(~ 0 + y + study)
  colnames(design) <- make.names(colnames(design))
  fit <- lmFit(t(x), design)
  contrast_matrix <- makeContrasts(
    C1_vs_C2 = yC1 - yC2,
    C1_vs_C3 = yC1 - yC3,
    C2_vs_C3 = yC2 - yC3,
    levels = design
  )
  fit <- eBayes(contrasts.fit(fit, contrast_matrix), trend = TRUE, robust = TRUE)
  tab <- topTable(fit, coef = seq_len(ncol(contrast_matrix)), number = Inf, sort.by = "F")
  rownames(tab)
}

scale_train_test <- function(train_x, test_x) {
  center <- colMeans(train_x)
  scale <- apply(train_x, 2, sd)
  scale[!is.finite(scale) | scale < 1e-8] <- 1
  train_scaled <- sweep(sweep(train_x, 2, center, "-"), 2, scale, "/")
  test_scaled <- sweep(sweep(test_x, 2, center, "-"), 2, scale, "/")
  list(train = train_scaled, test = test_scaled, center = center, scale = scale)
}

fit_predict_svm <- function(train_x, train_y, test_x, cost) {
  scaled <- scale_train_test(train_x, test_x)
  model <- svm(x = scaled$train, y = factor(train_y, levels = classes), kernel = "linear",
               cost = cost, scale = FALSE, probability = TRUE)
  pred <- predict(model, scaled$test, probability = TRUE)
  prob <- attr(pred, "probabilities")
  aligned <- matrix(0, nrow = nrow(test_x), ncol = length(classes),
                    dimnames = list(rownames(test_x), classes))
  aligned[, intersect(classes, colnames(prob))] <- prob[, intersect(classes, colnames(prob)), drop = FALSE]
  aligned <- aligned / rowSums(aligned)
  list(prob = aligned, pred = classes[max.col(aligned, ties.method = "first")],
       model = model, center = scaled$center, scale = scaled$scale)
}

balanced_accuracy <- function(truth, pred) {
  mean(vapply(classes, function(cl) mean(pred[truth == cl] == cl), numeric(1)), na.rm = TRUE)
}

tune_cost <- function(train_idx, panel_n, fold_id, seed_offset = 0L) {
  truth_all <- as.character(meta$Subtype[train_idx])
  pred_by_cost <- lapply(cost_grid, function(z) rep(NA_character_, length(train_idx)))
  for (f in sort(unique(fold_id))) {
    inner_test_local <- which(fold_id == f)
    inner_train_local <- which(fold_id != f)
    inner_train <- train_idx[inner_train_local]
    inner_test <- train_idx[inner_test_local]
    corrected <- fold_correct(inner_train, inner_test)
    ranked <- rank_features(corrected$train, meta$Subtype[inner_train], meta$Study[inner_train])
    genes <- head(ranked, panel_n)
    for (ci in seq_along(cost_grid)) {
      fit <- fit_predict_svm(corrected$train[, genes, drop = FALSE], meta$Subtype[inner_train],
                             corrected$test[, genes, drop = FALSE], cost_grid[[ci]])
      pred_by_cost[[ci]][inner_test_local] <- fit$pred
    }
  }
  score <- vapply(pred_by_cost, function(p) balanced_accuracy(truth_all, p), numeric(1))
  cost_grid[which.max(score)]
}

raw_predictions <- list()
selected_features <- list()
fold_definitions <- list()
counter <- 0L

for (repeat_id in seq_len(outer_repeats)) {
  outer_id <- make_folds(meta$Subtype, meta$Study, outer_folds, 20260715 + repeat_id)
  for (fold in seq_len(outer_folds)) {
    test_idx <- which(outer_id == fold)
    train_idx <- which(outer_id != fold)
    inner_id <- make_folds(meta$Subtype[train_idx], meta$Study[train_idx], inner_folds,
                           20261715 + repeat_id * 100 + fold)
    outer_corrected <- fold_correct(train_idx, test_idx)
    fold_definitions[[length(fold_definitions) + 1L]] <- data.frame(
      Analysis = "Repeated nested CV", Repeat = repeat_id, Fold = fold,
      Sample = meta$Sample, Study = meta$Study, Subtype = as.character(meta$Subtype),
      Role = ifelse(seq_len(nrow(meta)) %in% test_idx, "Outer test", "Outer train")
    )
    outer_ranking <- rank_features(outer_corrected$train, meta$Subtype[train_idx], meta$Study[train_idx])
    for (panel in names(panel_sizes)) {
      n_gene <- panel_sizes[[panel]]
      best_cost <- tune_cost(train_idx, n_gene, inner_id)
      genes <- head(outer_ranking, n_gene)
      fit <- fit_predict_svm(outer_corrected$train[, genes, drop = FALSE], meta$Subtype[train_idx],
                             outer_corrected$test[, genes, drop = FALSE], best_cost)
      counter <- counter + 1L
      raw_predictions[[counter]] <- data.frame(
        Analysis = "Repeated nested CV", Panel = panel, Repeat = repeat_id, Fold = fold,
        Sample = meta$Sample[test_idx], Study = meta$Study[test_idx],
        Truth = as.character(meta$Subtype[test_idx]), Prediction = fit$pred,
        ProbabilityC1 = fit$prob[, "C1"], ProbabilityC2 = fit$prob[, "C2"],
        ProbabilityC3 = fit$prob[, "C3"], Cost = best_cost, stringsAsFactors = FALSE
      )
      selected_features[[length(selected_features) + 1L]] <- data.frame(
        Analysis = "Repeated nested CV", Panel = panel, Repeat = repeat_id, Fold = fold,
        Cost = best_cost, Rank = seq_along(genes), Gene = genes, stringsAsFactors = FALSE
      )
    }
    cat(sprintf("Nested CV repeat %d/%d fold %d/%d complete\n", repeat_id, outer_repeats, fold, outer_folds))
  }
}

raw_predictions <- do.call(rbind, raw_predictions)
selected_features <- do.call(rbind, selected_features)
fold_definitions <- do.call(rbind, fold_definitions)
write.csv(raw_predictions, file.path(table_dir, "step5_nested_cv_predictions_by_repeat.csv"), row.names = FALSE)
write.csv(selected_features, file.path(table_dir, "step5_nested_cv_selected_features_long.csv"), row.names = FALSE)
write.csv(fold_definitions, file.path(table_dir, "step5_nested_cv_fold_definitions.csv"), row.names = FALSE)

aggregate_oof <- do.call(rbind, lapply(split(raw_predictions, interaction(raw_predictions$Panel, raw_predictions$Sample)), function(z) {
  p <- colMeans(z[, c("ProbabilityC1", "ProbabilityC2", "ProbabilityC3")])
  data.frame(Analysis = "Repeated nested CV (mean of five outer-test probabilities)", Panel = z$Panel[[1]],
             Sample = z$Sample[[1]], Study = z$Study[[1]], Truth = z$Truth[[1]],
             Prediction = classes[which.max(p)], ProbabilityC1 = p[[1]], ProbabilityC2 = p[[2]],
             ProbabilityC3 = p[[3]], MaxProbability = max(p), ProbabilityMargin = sort(p, decreasing = TRUE)[1] - sort(p, decreasing = TRUE)[2])
}))
row.names(aggregate_oof) <- NULL
write.csv(aggregate_oof, file.path(table_dir, "step5_nested_cv_oof_predictions.csv"), row.names = FALSE)

average_precision <- function(obs, score) {
  ord <- order(score, decreasing = TRUE)
  obs <- obs[ord]
  if (sum(obs) == 0) return(NA_real_)
  precision <- cumsum(obs) / seq_along(obs)
  sum(precision[obs == 1]) / sum(obs)
}

metric_vector <- function(dat) {
  truth <- factor(dat$Truth, levels = classes)
  pred <- factor(dat$Prediction, levels = classes)
  recalls <- vapply(classes, function(cl) mean(pred[truth == cl] == cl), numeric(1))
  precision <- vapply(classes, function(cl) {
    denom <- sum(pred == cl)
    if (denom == 0) return(NA_real_)
    sum(pred == cl & truth == cl) / denom
  }, numeric(1))
  f1 <- 2 * precision * recalls / (precision + recalls)
  prob <- as.matrix(dat[, paste0("Probability", classes)])
  onehot <- sapply(classes, function(cl) as.numeric(truth == cl))
  out <- c(BalancedAccuracy = mean(recalls), MacroF1 = mean(f1, na.rm = TRUE),
           MulticlassBrier = mean(rowSums((prob - onehot)^2)))
  for (cl in classes) {
    obs <- as.numeric(truth == cl)
    score <- dat[[paste0("Probability", cl)]]
    auc <- as.numeric(roc(obs, score, quiet = TRUE, direction = "<")$auc)
    out[paste0(cl, "_Sensitivity")] <- recalls[[cl]]
    out[paste0(cl, "_Precision")] <- precision[[cl]]
    out[paste0(cl, "_AUROC")] <- auc
    out[paste0(cl, "_AUPRC")] <- average_precision(obs, score)
  }
  out
}

metric_table <- function(dat, analysis, panel, boot = bootstrap_reps, seed = 1L) {
  estimate <- metric_vector(dat)
  set.seed(seed)
  boots <- replicate(boot, {
    idx <- unlist(lapply(classes, function(cl) {
      w <- which(dat$Truth == cl)
      sample(w, length(w), replace = TRUE)
    }))
    metric_vector(dat[idx, , drop = FALSE])
  })
  data.frame(Analysis = analysis, Panel = panel, Metric = names(estimate), Estimate = as.numeric(estimate),
             Lower95 = apply(boots, 1, quantile, 0.025, na.rm = TRUE),
             Upper95 = apply(boots, 1, quantile, 0.975, na.rm = TRUE), row.names = NULL)
}

performance <- do.call(rbind, lapply(names(panel_sizes), function(panel) {
  metric_table(subset(aggregate_oof, Panel == panel), "Repeated nested CV", panel,
               seed = 20262715 + panel_sizes[[panel]])
}))

confusion <- do.call(rbind, lapply(names(panel_sizes), function(panel) {
  z <- subset(aggregate_oof, Panel == panel)
  tab <- as.data.frame(table(Truth = factor(z$Truth, levels = classes),
                             Prediction = factor(z$Prediction, levels = classes)))
  tab$Panel <- panel
  tab$RowFraction <- ave(tab$Freq, tab$Truth, FUN = function(x) x / sum(x))
  tab
}))
write.csv(confusion, file.path(table_dir, "step5_nested_cv_confusion_matrix.csv"), row.names = FALSE)

calibration <- do.call(rbind, lapply(names(panel_sizes), function(panel) {
  z <- subset(aggregate_oof, Panel == panel)
  do.call(rbind, lapply(classes, function(cl) {
    obs <- as.numeric(z$Truth == cl)
    prob <- pmin(pmax(z[[paste0("Probability", cl)]], 1e-6), 1 - 1e-6)
    logit <- qlogis(prob)
    fit <- glm(obs ~ logit, family = binomial())
    data.frame(Analysis = "Repeated nested CV", Panel = panel, Class = cl,
               N = nrow(z), Events = sum(obs), Brier = mean((obs - prob)^2),
               CalibrationIntercept = unname(coef(fit)[1]), CalibrationSlope = unname(coef(fit)[2]))
  }))
}))
write.csv(calibration, file.path(table_dir, "step5_nested_cv_calibration_metrics.csv"), row.names = FALSE)

calibration_bins <- do.call(rbind, lapply(names(panel_sizes), function(panel) {
  z <- subset(aggregate_oof, Panel == panel)
  do.call(rbind, lapply(classes, function(cl) {
    prob <- z[[paste0("Probability", cl)]]
    bin <- cut(prob, breaks = seq(0, 1, 0.1), include.lowest = TRUE)
    a <- aggregate(cbind(Predicted = prob, Observed = as.numeric(z$Truth == cl)), list(Bin = bin), mean)
    a$N <- as.integer(table(bin)[as.character(a$Bin)])
    a$Panel <- panel
    a$Class <- cl
    a
  }))
}))
write.csv(calibration_bins, file.path(table_dir, "step5_nested_cv_calibration_bins.csv"), row.names = FALSE)

# Cohort-standardized transfer sensitivity. Study offsets for training cohorts
# are estimated without the held-out study. The held-out cohort contributes
# only its unlabeled gene means for transductive centering, so this is not a
# single-sample or independently labelled validation claim.
loso_predictions <- list()
loso_features <- list()
for (held_study in unique(meta$Study)) {
  test_idx <- which(meta$Study == held_study)
  train_idx <- which(meta$Study != held_study)
  inner_id <- as.integer(factor(meta$Study[train_idx], levels = unique(meta$Study[train_idx])))
  best_cost <- tune_cost(train_idx, 123L, inner_id)
  loso_corrected <- fold_correct(train_idx, test_idx)
  ranked <- rank_features(loso_corrected$train, meta$Subtype[train_idx], meta$Study[train_idx])
  genes <- head(ranked, 123L)
  fit <- fit_predict_svm(loso_corrected$train[, genes, drop = FALSE], meta$Subtype[train_idx],
                         loso_corrected$test[, genes, drop = FALSE], best_cost)
  loso_predictions[[length(loso_predictions) + 1L]] <- data.frame(
    Analysis = "Cohort-standardized transfer sensitivity", HeldOutStudy = held_study, Sample = meta$Sample[test_idx],
    Truth = as.character(meta$Subtype[test_idx]), Prediction = fit$pred,
    ProbabilityC1 = fit$prob[, "C1"], ProbabilityC2 = fit$prob[, "C2"], ProbabilityC3 = fit$prob[, "C3"],
    Cost = best_cost, stringsAsFactors = FALSE
  )
  loso_features[[length(loso_features) + 1L]] <- data.frame(HeldOutStudy = held_study, Cost = best_cost,
                                                            Rank = seq_along(genes), Gene = genes)
  cat("LOSO complete: ", held_study, "\n", sep = "")
}
loso_predictions <- do.call(rbind, loso_predictions)
loso_features <- do.call(rbind, loso_features)
write.csv(loso_predictions, file.path(table_dir, "step5_loso_predictions.csv"), row.names = FALSE)
write.csv(loso_features, file.path(table_dir, "step5_loso_selected_features_long.csv"), row.names = FALSE)
performance <- rbind(performance,
                     metric_table(transform(loso_predictions, Panel = "Operational123"),
                                  "Pooled cohort-standardized transfer sensitivity", "Operational123", seed = 20263715))

loso_by_study <- do.call(rbind, lapply(split(loso_predictions, loso_predictions$HeldOutStudy), function(z) {
  m <- metric_vector(z)
  data.frame(Study = z$HeldOutStudy[[1]], N = nrow(z),
             C1N = sum(z$Truth == "C1"), C2N = sum(z$Truth == "C2"), C3N = sum(z$Truth == "C3"),
             BalancedAccuracy = m[["BalancedAccuracy"]], MacroF1 = m[["MacroF1"]],
             MulticlassBrier = m[["MulticlassBrier"]])
}))
write.csv(loso_by_study, file.path(table_dir, "step5_loso_per_study_performance.csv"), row.names = FALSE)
write.csv(performance, file.path(table_dir, "step5_classifier_performance_with_bootstrap_ci.csv"), row.names = FALSE)

# Fit deployment models only after unbiased performance estimation is complete.
tune_grouped_all <- function(panel_n) {
  fold_id <- as.integer(factor(meta$Study, levels = unique(meta$Study)))
  tune_cost(seq_len(nrow(meta)), panel_n, fold_id)
}

full_corrected <- fold_correct(seq_len(nrow(meta)), seq_len(nrow(meta)))
full_ranking <- rank_features(full_corrected$train, meta$Subtype, meta$Study)
final_models <- list()
final_features <- list()
for (panel in names(panel_sizes)) {
  n_gene <- panel_sizes[[panel]]
  genes <- head(full_ranking, n_gene)
  best_cost <- tune_grouped_all(n_gene)
  scaled <- scale_train_test(full_corrected$train[, genes, drop = FALSE], full_corrected$train[, genes, drop = FALSE])
  model <- svm(x = scaled$train, y = meta$Subtype, kernel = "linear", cost = best_cost,
               scale = FALSE, probability = TRUE)
  final_models[[panel]] <- list(
    model = model, genes = genes, center = scaled$center, scale = scaled$scale,
    classes = classes, cost = best_cost, training_samples = meta$Sample,
    training_studies = unique(meta$Study), feature_selection = "Study-adjusted limma moderated F statistic",
    preprocessing_policy = "study mean-offset correction from before-correction matrix; new-study cohort-transductive centering",
    performance_source = "step5_nested_cv_oof_predictions.csv and step5_loso_predictions.csv"
  )
  final_features[[length(final_features) + 1L]] <- data.frame(Panel = panel, Rank = seq_along(genes), Gene = genes,
                                                              Cost = best_cost)
}
final_features <- do.call(rbind, final_features)
write.csv(final_features, file.path(table_dir, "step5_final_model_features.csv"), row.names = FALSE)
write.csv(subset(final_features, Panel == "Operational123"),
          file.path(table_dir, "step5_final_123_gene_signature.csv"), row.names = FALSE)
saveRDS(final_models, file.path(model_dir, "step5_final_linear_svm_models.rds"), compress = "xz")
saveRDS(final_models$Operational123, "5-ClassifierValidation/iCCA_Subtype_SVM_Classifier.rds", compress = "xz")
saveRDS(final_models, "5-ClassifierValidation/final_trained_models.rds", compress = "xz")
saveRDS(list(operational_123 = final_models$Operational123$genes,
             nested_31 = final_models$Nested31$genes,
             compact_15 = final_models$Compact15$genes,
             historical_interpretive_31 = legacy_31,
             historical_123_for_comparison_only = legacy_123),
        "5-ClassifierValidation/selected_signature_genes.rds", compress = "xz")
saveRDS(list(expression_source = "1-Clustering/data/combined_expr_before_correction.rds",
             preprocessing_policy = "training-study offsets; held-out/new study unlabeled cohort centering",
             sample_metadata = meta, panel_features = final_features),
        "5-ClassifierValidation/input_validation_data.rds", compress = "xz")

feature_stability <- aggregate(list(OuterFoldsSelected = selected_features$Gene),
                               list(Panel = selected_features$Panel, Gene = selected_features$Gene), length)
feature_stability$SelectionFraction <- feature_stability$OuterFoldsSelected / (outer_repeats * outer_folds)
feature_stability <- feature_stability[order(feature_stability$Panel, -feature_stability$SelectionFraction,
                                             feature_stability$Gene), ]
write.csv(feature_stability, file.path(table_dir, "step5_nested_cv_feature_stability.csv"), row.names = FALSE)

save_pair <- function(plot, stem, width, height) {
  ggsave(file.path(figure_dir, paste0(stem, ".pdf")), plot, width = width, height = height, device = cairo_pdf)
  ggsave(file.path(figure_dir, paste0(stem, ".png")), plot, width = width, height = height, dpi = 400)
}

operational_conf <- subset(confusion, Panel == "Operational123")
p_conf <- ggplot(operational_conf, aes(Prediction, Truth, fill = RowFraction)) +
  geom_tile(color = "white") + geom_text(aes(label = sprintf("%d\n%.1f%%", Freq, 100 * RowFraction)), size = 3.3) +
  scale_fill_gradient(low = "white", high = "#3B6FB6", limits = c(0, 1)) +
  coord_equal() + labs(title = "Operational 123-gene classifier", subtitle = "Repeated nested-CV out-of-fold predictions",
                       x = "Predicted subtype", y = "Locked subtype", fill = "Row fraction") +
  theme_classic(base_size = 10) + theme(legend.position = "right")
save_pair(p_conf, "step5_nested_cv_confusion_matrix", 5.5, 4.8)

perf_plot <- subset(performance, Analysis == "Repeated nested CV" & Metric %in% c("BalancedAccuracy", "MacroF1"))
perf_plot$Panel <- factor(perf_plot$Panel, levels = c("Operational123", "Nested31", "Compact15"))
perf_plot$MetricLabel <- factor(ifelse(perf_plot$Metric == "BalancedAccuracy", "Balanced accuracy", "Macro-F1"),
                                levels = c("Balanced accuracy", "Macro-F1"))
p_perf <- ggplot(perf_plot, aes(Panel, Estimate, ymin = Lower95, ymax = Upper95, color = MetricLabel)) +
  geom_pointrange(position = position_dodge(width = 0.35)) +
  scale_y_continuous(limits = c(0.90, 1.005), breaks = seq(0.90, 1.00, 0.025)) +
  scale_color_manual(values = c("Balanced accuracy" = "#0072B2", "Macro-F1" = "#D55E00")) +
  labs(title = "Leakage-controlled classifier comparison", subtitle = "Focused 0.90–1.00 scale; identical outer splits and within-fold feature selection",
       x = NULL, y = "Out-of-fold estimate (95% bootstrap CI)", color = NULL) +
  theme_classic(base_size = 10) + theme(legend.position = "top")
save_pair(p_perf, "step5_nested_cv_panel_comparison", 6.4, 4.5)

p_cal <- ggplot(subset(calibration_bins, Panel == "Operational123"), aes(Predicted, Observed, color = Class)) +
  geom_abline(slope = 1, intercept = 0, linetype = 2, color = "grey55") +
  geom_line() + geom_point(aes(size = N), alpha = 0.9) +
  scale_color_manual(values = subtype_colors) + coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = "Operational classifier calibration", subtitle = "Repeated nested-CV out-of-fold probabilities",
       x = "Mean predicted probability", y = "Observed frequency", color = NULL, size = "Bin N") +
  theme_classic(base_size = 10) + theme(legend.position = "top")
save_pair(p_cal, "step5_nested_cv_calibration", 5.5, 5.0)

p_loso <- ggplot(loso_by_study, aes(reorder(Study, BalancedAccuracy), BalancedAccuracy)) +
  geom_col(fill = "#4D4D4D", width = 0.72) + geom_text(aes(label = paste0("N=", N)), hjust = -0.1, size = 3) +
  coord_flip(ylim = c(0, 1.05)) +
  labs(title = "Cohort-standardized transfer sensitivity", subtitle = "123 genes reselected without held-out labels; unlabeled held-out cohort centering",
       x = NULL, y = "Balanced accuracy") + theme_classic(base_size = 10)
save_pair(p_loso, "step5_loso_per_study_performance", 6.5, 4.6)

contract <- data.frame(
  Item = c("Training population", "Locked labels", "Operational algorithm", "Feature selection",
           "Preprocessing", "Hyperparameter tuning", "Headline internal validation", "Cross-study transfer sensitivity",
           "Batch-correction provenance",
           "Historical 31-gene set", "Prognostic model"),
  Specification = c("497 principal-state tumors", "C1/C2/C3 from Step1", "123-gene linear SVM",
                    "Study-adjusted limma moderated F statistic inside every training fold",
                    "Fold-safe study-offset correction from before-correction matrix; training-fold centering/scaling; held-out study uses unlabeled cohort centering",
                    "Inner CV over C = 0.01, 0.1, 1, 10", "Five-times repeated five-fold nested CV",
                    "Held-out-study transfer sensitivity; tuning and feature selection use training studies only",
                    "Before-correction expression is used; seen-study offsets are training-only and held-out-study centering is unlabeled/transductive",
                    "Interpretive core only; not the fixed headline feature set",
                    "Deferred to the later clinical gate; not mixed with subtype-classifier validation"),
  stringsAsFactors = FALSE
)
write.csv(contract, file.path(table_dir, "step5_classifier_analysis_contract.csv"), row.names = FALSE)

cat("Step5 core classifier analysis complete.\n")
