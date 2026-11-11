#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(limma)
  library(e1071)
  library(openxlsx)
  library(pROC)
  library(ggplot2)
})

set.seed(20260715)
table_dir <- "5-ClassifierValidation/publication/tables"
figure_dir <- "5-ClassifierValidation/publication/figures"
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

e <- new.env()
load("1-Clustering/data/combined_expr_before_correction.rds", envir = e)
expr <- e$combined_expr
meta <- read.xlsx("1-Clustering/data/sample_info_with_clusters.xlsx")
meta <- meta[meta$Ensemble %in% c("C1", "C2", "C3"), c("Sample", "Study", "Ensemble")]
colnames(meta) <- c("Sample", "Study", "Subtype")
meta$Sample <- as.character(meta$Sample)

duct <- read.csv("3-DuctType/publication/tables/step3_nested_cv_oof_predictions.csv", check.names = FALSE)
duct$Sample <- as.character(duct$Sample)
dat <- merge(meta, duct[, c("Sample", "ProbabilityLargeDuct")], by = "Sample", all = FALSE)
dat <- dat[match(intersect(dat$Sample, colnames(expr)), dat$Sample), ]
X <- t(as.matrix(e$combined_expr[, dat$Sample, drop = FALSE]))
storage.mode(X) <- "double"
outcome <- factor(ifelse(dat$Subtype == "C3", "C3", "NotC3"), levels = c("NotC3", "C3"))
stopifnot(nrow(dat) == 93L, identical(rownames(X), dat$Sample))

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
        mu <- colMeans(raw[ii, , drop = FALSE], na.rm = TRUE)
        mu[!is.finite(mu)] <- global[!is.finite(mu)]; mu - global
      }
      out[ii, ] <- sweep(raw[ii, , drop = FALSE], 2L, off, "-")
    }
    out[!is.finite(out)] <- 0; out
  }
  list(train = adjust(train_idx), test = adjust(test_idx))
}

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

rank_binary <- function(x, y, study) {
  design <- model.matrix(~ 0 + y + factor(study))
  colnames(design) <- make.names(colnames(design))
  fit <- lmFit(t(x), design)
  contrast <- makeContrasts(yC3 - yNotC3, levels = design)
  fit <- eBayes(contrasts.fit(fit, contrast), trend = TRUE, robust = TRUE)
  rownames(topTable(fit, coef = 1, number = Inf, sort.by = "P"))
}

fit_predict <- function(train_x, train_y, test_x, cost) {
  center <- colMeans(train_x)
  scale <- apply(train_x, 2, sd)
  scale[!is.finite(scale) | scale < 1e-8] <- 1
  tr <- sweep(sweep(train_x, 2, center, "-"), 2, scale, "/")
  te <- sweep(sweep(test_x, 2, center, "-"), 2, scale, "/")
  model <- svm(tr, train_y, kernel = "linear", cost = cost, scale = FALSE, probability = TRUE)
  pred <- predict(model, te, probability = TRUE)
  prob <- attr(pred, "probabilities")
  as.numeric(prob[, "C3"])
}

balanced_accuracy <- function(y, p) {
  pred <- ifelse(p >= 0.5, "C3", "NotC3")
  mean(c(mean(pred[y == "C3"] == "C3"), mean(pred[y == "NotC3"] == "NotC3")))
}

tune_cost <- function(train_idx, inner_id, include_duct) {
  grid <- c(0.01, 0.1, 1, 10)
  predictions <- lapply(grid, function(x) rep(NA_real_, length(train_idx)))
  for (f in sort(unique(inner_id))) {
    val_local <- which(inner_id == f)
    tr_local <- which(inner_id != f)
    tr <- train_idx[tr_local]
    va <- train_idx[val_local]
    corrected <- fold_correct(tr, va)
    genes <- head(rank_binary(corrected$train, outcome[tr], dat$Study[tr]), 15L)
    xtr <- corrected$train[, genes, drop = FALSE]
    xva <- corrected$test[, genes, drop = FALSE]
    if (include_duct) {
      xtr <- cbind(xtr, CrossFittedDuctProbability = dat$ProbabilityLargeDuct[tr])
      xva <- cbind(xva, CrossFittedDuctProbability = dat$ProbabilityLargeDuct[va])
    }
    for (i in seq_along(grid)) predictions[[i]][val_local] <- fit_predict(xtr, outcome[tr], xva, grid[[i]])
  }
  scores <- vapply(predictions, function(p) balanced_accuracy(outcome[train_idx], p), numeric(1))
  grid[which.max(scores)]
}

predictions <- list()
features <- list()
for (repeat_id in 1:5) {
  outer <- make_folds(outcome, dat$Study, 5, 20260715 + repeat_id)
  for (fold in 1:5) {
    test <- which(outer == fold)
    train <- which(outer != fold)
    inner <- make_folds(outcome[train], dat$Study[train], 4, 20261715 + repeat_id * 100 + fold)
    corrected <- fold_correct(train, test)
    genes <- head(rank_binary(corrected$train, outcome[train], dat$Study[train]), 15L)
    for (model_name in c("Molecular15", "Molecular15PlusCrossFittedDuct")) {
      include_duct <- model_name == "Molecular15PlusCrossFittedDuct"
      cost <- tune_cost(train, inner, include_duct)
    xtr <- corrected$train[, genes, drop = FALSE]
    xte <- corrected$test[, genes, drop = FALSE]
      if (include_duct) {
        xtr <- cbind(xtr, CrossFittedDuctProbability = dat$ProbabilityLargeDuct[train])
        xte <- cbind(xte, CrossFittedDuctProbability = dat$ProbabilityLargeDuct[test])
      }
      p <- fit_predict(xtr, outcome[train], xte, cost)
      predictions[[length(predictions) + 1L]] <- data.frame(
        Repeat = repeat_id, Fold = fold, Model = model_name, Sample = dat$Sample[test], Study = dat$Study[test],
        Truth = as.character(outcome[test]), ProbabilityC3 = p, Prediction = ifelse(p >= 0.5, "C3", "NotC3"),
        CrossFittedDuctProbability = dat$ProbabilityLargeDuct[test], Cost = cost
      )
    }
    features[[length(features) + 1L]] <- data.frame(Repeat = repeat_id, Fold = fold,
                                                    Rank = seq_along(genes), Gene = genes)
  }
}
predictions <- do.call(rbind, predictions)
features <- do.call(rbind, features)
write.csv(predictions, file.path(table_dir, "step5_c3_incremental_nested_cv_predictions_by_repeat.csv"), row.names = FALSE)
write.csv(features, file.path(table_dir, "step5_c3_incremental_selected_features_long.csv"), row.names = FALSE)

oof <- do.call(rbind, lapply(split(predictions, interaction(predictions$Model, predictions$Sample)), function(z) {
  p <- mean(z$ProbabilityC3)
  data.frame(Model = z$Model[[1]], Sample = z$Sample[[1]], Study = z$Study[[1]], Truth = z$Truth[[1]],
             ProbabilityC3 = p, Prediction = ifelse(p >= 0.5, "C3", "NotC3"),
             CrossFittedDuctProbability = z$CrossFittedDuctProbability[[1]])
}))
row.names(oof) <- NULL
write.csv(oof, file.path(table_dir, "step5_c3_incremental_nested_cv_oof_predictions.csv"), row.names = FALSE)

average_precision <- function(obs, score) {
  ord <- order(score, decreasing = TRUE)
  obs <- obs[ord]
  precision <- cumsum(obs) / seq_along(obs)
  sum(precision[obs == 1]) / sum(obs)
}

metrics <- function(z) {
  obs <- as.numeric(z$Truth == "C3")
  c(AUROC = as.numeric(roc(obs, z$ProbabilityC3, quiet = TRUE, direction = "<")$auc),
    AUPRC = average_precision(obs, z$ProbabilityC3),
    BalancedAccuracy = balanced_accuracy(z$Truth, z$ProbabilityC3),
    Brier = mean((obs - z$ProbabilityC3)^2))
}

models <- c("Molecular15", "Molecular15PlusCrossFittedDuct")
est <- lapply(models, function(m) metrics(subset(oof, Model == m)))
names(est) <- models
set.seed(20260715)
boot <- replicate(2000, {
  idx <- unlist(lapply(c("NotC3", "C3"), function(cl) {
    w <- which(subset(oof, Model == models[[1]])$Truth == cl)
    sample(w, length(w), replace = TRUE)
  }))
  a <- metrics(subset(oof, Model == models[[1]])[idx, ])
  b <- metrics(subset(oof, Model == models[[2]])[idx, ])
  c(a, b, Delta = b - a)
})

summary_rows <- list()
for (metric in names(est[[1]])) {
  for (model in models) {
    prefix <- if (model == models[[1]]) metric else paste0(metric, "1")
    # replicate() makes duplicate names unique only when converted; use positional indices explicitly.
    pos <- match(metric, names(est[[1]])) + ifelse(model == models[[1]], 0, length(est[[1]]))
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      Comparison = model, Metric = metric, Estimate = est[[model]][[metric]],
      Lower95 = quantile(boot[pos, ], 0.025, na.rm = TRUE), Upper95 = quantile(boot[pos, ], 0.975, na.rm = TRUE)
    )
  }
  pos_delta <- 2 * length(est[[1]]) + match(metric, names(est[[1]]))
  summary_rows[[length(summary_rows) + 1L]] <- data.frame(
    Comparison = "Delta: plus duct minus molecular", Metric = metric,
    Estimate = est[[models[[2]]]][[metric]] - est[[models[[1]]]][[metric]],
    Lower95 = quantile(boot[pos_delta, ], 0.025, na.rm = TRUE),
    Upper95 = quantile(boot[pos_delta, ], 0.975, na.rm = TRUE)
  )
}
summary_table <- do.call(rbind, summary_rows)
summary_table$N <- nrow(subset(oof, Model == models[[1]]))
summary_table$C3N <- sum(subset(oof, Model == models[[1]])$Truth == "C3")
summary_table$Interpretation <- ifelse(grepl("Delta", summary_table$Comparison),
                                      "Paired sample-bootstrap delta; exploratory observed-duct subset; fold-safe before-correction study offsets",
                                      "Repeated nested-CV out-of-fold estimate; fold-safe before-correction study offsets")
write.csv(summary_table, file.path(table_dir, "step5_c3_incremental_value_summary.csv"), row.names = FALSE)

thresholds <- seq(0.1, 0.8, 0.05)
decision <- do.call(rbind, lapply(models, function(m) {
  z <- subset(oof, Model == m)
  obs <- z$Truth == "C3"
  do.call(rbind, lapply(thresholds, function(t) {
    pred <- z$ProbabilityC3 >= t
    data.frame(Model = m, Threshold = t,
               NetBenefit = sum(pred & obs) / length(obs) - sum(pred & !obs) / length(obs) * t / (1 - t))
  }))
}))
prevalence <- mean(subset(oof, Model == models[[1]])$Truth == "C3")
decision <- rbind(decision,
                  data.frame(Model = "Treat all as C3", Threshold = thresholds,
                             NetBenefit = prevalence - (1 - prevalence) * thresholds / (1 - thresholds)),
                  data.frame(Model = "Treat none as C3", Threshold = thresholds, NetBenefit = 0))
write.csv(decision, file.path(table_dir, "step5_c3_incremental_decision_curve.csv"), row.names = FALSE)

p <- ggplot(decision, aes(Threshold, NetBenefit, color = Model)) +
  geom_hline(yintercept = 0, color = "grey75") + geom_line(linewidth = 0.85) +
  scale_color_manual(values = c("Molecular15" = "#0072B2", "Molecular15PlusCrossFittedDuct" = "#D55E00",
                                "Treat all as C3" = "grey45", "Treat none as C3" = "grey75"),
                     labels = c("Molecular15" = "15-gene molecular",
                                "Molecular15PlusCrossFittedDuct" = "15-gene + cross-fitted duct",
                                "Treat all as C3" = "Treat all as C3",
                                "Treat none as C3" = "Treat none as C3")) +
  labs(title = "Incremental value of cross-fitted duct probability", subtitle = "C3 detection in 93 tumors with observed duct annotation",
       x = "Decision threshold", y = "Net benefit", color = NULL) +
  theme_classic(base_size = 10) + theme(legend.position = "top")
ggsave(file.path(figure_dir, "step5_c3_incremental_decision_curve.pdf"), p, width = 6.8, height = 4.8, device = cairo_pdf)
ggsave(file.path(figure_dir, "step5_c3_incremental_decision_curve.png"), p, width = 6.8, height = 4.8, dpi = 400)

cat("Step5 C3 incremental-value analysis complete.\n")
