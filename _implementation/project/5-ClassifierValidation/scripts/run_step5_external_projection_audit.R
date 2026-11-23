#!/usr/bin/env Rscript

# External projection audit for Step 5.  This is deliberately separate from
# accuracy validation: the external cohorts have no independent locked labels.
# External platforms are standardized gene-wise within each cohort before being
# passed to the fixed SVM.  This is a transductive cohort projection (not a
# single-sample deployment rule and not accuracy validation).  Missing genes
# are imputed to zero in that cohort-standardized space.
options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({library(openxlsx); library(e1071); library(ggplot2)})

out_dir <- "5-ClassifierValidation/publication/tables"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
model <- readRDS("5-ClassifierValidation/iCCA_Subtype_SVM_Classifier.rds")
stopifnot(length(model$genes) == 123L, identical(model$classes, c("C1", "C2", "C3")))

align_predict <- function(z, source, study = NA_character_) {
  z <- as.matrix(z)
  storage.mode(z) <- "double"
  genes <- model$genes
  available <- intersect(genes, rownames(z))
  if (length(available) == 0L) stop("No deployment genes found in ", source)
  # The deployment SVM was fit on gene-wise standardized training features.
  # For a platform with a different assay scale, standardize each available gene
  # within this cohort.  This intentionally uses the whole external cohort and
  # must be labelled transductive; it cannot be used for single-sample calls.
  sample_names <- colnames(z)
  z <- z[available, , drop = FALSE]
  z <- t(apply(z, 1L, function(v) {
    v <- as.numeric(v)
    v[!is.finite(v)] <- NA_real_
    mu <- mean(v, na.rm = TRUE)
    sig <- sd(v, na.rm = TRUE)
    if (!is.finite(mu)) mu <- 0
    if (!is.finite(sig) || sig < 1e-8) sig <- 1
    (v - mu) / sig
  }))
  dimnames(z) <- list(available, sample_names)
  z[!is.finite(z)] <- 0
  x <- matrix(0, nrow = ncol(z), ncol = length(genes),
              dimnames = list(colnames(z), genes))
  x[, available] <- t(z)
  pred <- predict(model$model, x, probability = TRUE)
  pp <- attr(pred, "probabilities")
  prob <- matrix(0, nrow = nrow(x), ncol = 3L,
                 dimnames = list(rownames(x), model$classes))
  prob[, intersect(colnames(pp), model$classes)] <- pp[, intersect(colnames(pp), model$classes), drop = FALSE]
  prob <- prob / rowSums(prob)
  ps <- apply(prob, 1L, function(v) sort(v, decreasing = TRUE)[1L] - sort(v, decreasing = TRUE)[2L])
  out <- data.frame(Platform = source, Sample = rownames(x),
                    MatchedFeatures = length(available), TotalFeatures = length(genes),
                    MissingFeatures = length(genes) - length(available),
                    MissingRate = 1 - length(available) / length(genes),
                    ScalePolicy = "cohort-wise gene z-score (transductive)",
                    ImputationRule = "0 in cohort-standardized space",
                    Prediction = model$classes[max.col(prob, ties.method = "first")],
                    ProbabilityC1 = prob[, "C1"], ProbabilityC2 = prob[, "C2"],
                    ProbabilityC3 = prob[, "C3"], MaxProbability = apply(prob, 1L, max),
                    ProbabilityMargin = ps, ProjectionOnly = TRUE,
                    IndependentTruthAvailable = FALSE, stringsAsFactors = FALSE)
  out
}

tables <- list()
add_env <- function(path, keep = NULL) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  nms <- ls(e)
  # Use the explicitly processed log2 expression objects; raw platform objects
  # are retained in the RData only for provenance and must not be projected.
  nms <- nms[grepl("^log2dat_", nms)]
  nms <- nms[vapply(nms, function(n) {
    x <- get(n, e)
    is.matrix(x) || is.data.frame(x)
  }, logical(1))]
  if (!is.null(keep)) nms <- intersect(nms, keep)
  for (n in nms) {
    x <- as.matrix(get(n, e))
    if (is.null(rownames(x)) || is.null(colnames(x))) next
    # Collapse duplicate probe/gene rows deterministically.
    if (anyDuplicated(rownames(x))) {
      x <- rowsum(x, group = rownames(x), reorder = FALSE) /
        as.vector(table(rownames(x))[unique(rownames(x))])
    }
    tables[[length(tables) + 1L]] <<- align_predict(x, n)
  }
}

add_env("0-Data/qc/0-Data-rnaseq-microarray1.Rdata")
add_env("0-Data/qc/0-Data-rnaseq-microarray2.Rdata")

# DepMap/CCLE is an optional, large input.  Read only the 123 deployment genes
# when requested; this keeps the default notebook audit lightweight.
if (identical(Sys.getenv("SCLUCID_STEP4_INCLUDE_CCLE", "0"), "1")) {
  if (!requireNamespace("data.table", quietly = TRUE)) stop("data.table is required for CCLE projection")
  hdr <- names(data.table::fread("5-ClassifierValidation/OmicsExpressionTPMLogp1HumanAllGenes.csv", nrows = 0L, showProgress = FALSE))
  hit <- vapply(model$genes, function(g) {
    z <- hdr[startsWith(hdr, paste0(g, " ("))]
    if (length(z)) z[[1L]] else NA_character_
  }, character(1))
  hit <- hit[is.finite(match(hit, hdr))]
  dat <- data.table::fread("5-ClassifierValidation/OmicsExpressionTPMLogp1HumanAllGenes.csv",
                           select = c("ModelID", "is_default_entry", unname(hit)),
                           showProgress = FALSE)
  allowed <- read.csv("5-ClassifierValidation/Model.csv", check.names = FALSE)
  allowed <- allowed$ModelID[allowed$OncotreeSubtype == "Intrahepatic Cholangiocarcinoma"]
  dat <- dat[dat$ModelID %in% allowed & dat$is_default_entry == TRUE, , drop = FALSE]
  x <- as.data.frame(dat)
  rownames(x) <- x$ModelID
  x$ModelID <- NULL; x$is_default_entry <- NULL
  colnames(x) <- sub(" \\(.*$", "", colnames(x))
  x <- t(as.matrix(x))
  tables[[length(tables) + 1L]] <- align_predict(x, "DepMap/CCLE")
}

sample_out <- do.call(rbind, tables)
sample_out <- sample_out[order(sample_out$Platform, sample_out$Sample), ]
write.csv(sample_out, file.path(out_dir, "step5_external_projection_sample_predictions.csv"), row.names = FALSE)

cohort_out <- do.call(rbind, lapply(split(sample_out, sample_out$Platform), function(z) {
  data.frame(Platform = z$Platform[[1L]], N = nrow(z), MatchedFeatures = z$MatchedFeatures[[1L]],
             TotalFeatures = z$TotalFeatures[[1L]], MissingFeatures = z$MissingFeatures[[1L]],
             MissingRate = z$MissingRate[[1L]], PredictedC1 = sum(z$Prediction == "C1"),
             PredictedC2 = sum(z$Prediction == "C2"), PredictedC3 = sum(z$Prediction == "C3"),
             MeanMaxProbability = mean(z$MaxProbability), MedianMargin = median(z$ProbabilityMargin),
             ScalePolicy = z$ScalePolicy[[1L]], ImputationRule = z$ImputationRule[[1L]], ProjectionOnly = TRUE,
             IndependentTruthAvailable = FALSE, stringsAsFactors = FALSE)
}))
write.csv(cohort_out, file.path(out_dir, "step5_external_projection_preservation.csv"), row.names = FALSE)

# Projection coverage/uncertainty only: no external cohort has independent
# subtype truth, so the panel deliberately avoids an accuracy interpretation.
plot_dat <- rbind(
  data.frame(Platform = cohort_out$Platform, N = cohort_out$N,
             Metric = "Matched deployment features", Value = 1 - cohort_out$MissingRate),
  data.frame(Platform = cohort_out$Platform, N = cohort_out$N,
             Metric = "Median assignment margin", Value = cohort_out$MedianMargin)
)
plot_dat$Platform <- sub("^log2dat_", "", plot_dat$Platform)
plot_dat$Platform[plot_dat$Platform == "2017Jusakual1"] <- "2017Jusakual (array 1)"
plot_dat$Platform[plot_dat$Platform == "2017Jusakual2"] <- "2017Jusakual (array 2)"
plot_dat$Metric <- factor(plot_dat$Metric,
                          levels = c("Matched deployment features",
                                     "Median assignment margin"))
p <- ggplot(plot_dat, aes(Value, reorder(Platform, Value), color = Metric)) +
  geom_segment(aes(x = 0, xend = Value, yend = Platform), color = "grey85", linewidth = 0.55) +
  geom_point(size = 2.4) +
  geom_text(aes(label = sprintf("%.2f  (N=%d)", Value, N)), hjust = -0.08,
            color = "#333333", size = 3) +
  facet_grid(Metric ~ ., scales = "free_y", space = "free_y") +
  scale_x_continuous(limits = c(0, 1.18), breaks = seq(0, 1, 0.25)) +
  scale_color_manual(values = c("Matched deployment features" = "#0072B2",
                                "Median assignment margin" = "#E69F00"), guide = "none") +
  labs(title = "External projection coverage and assignment margins",
       subtitle = "Cohort-standardized projections; margins are not external accuracy or calibration",
       x = "Fraction", y = NULL) +
  theme_classic(base_size = 10) +
  theme(strip.background = element_rect(fill = "grey95", color = "grey80"),
        strip.text = element_text(face = "bold"), plot.title.position = "plot")
figure_dir <- "5-ClassifierValidation/publication/figures"
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(figure_dir, "step5_external_projection_coverage_uncertainty.pdf"), p,
       width = 7.2, height = 6.8, device = cairo_pdf)
ggsave(file.path(figure_dir, "step5_external_projection_coverage_uncertainty.png"), p,
       width = 7.2, height = 6.8, dpi = 400)
cat("External projection audit complete:", nrow(sample_out), "samples across", nrow(cohort_out), "platform cohorts.\n")
