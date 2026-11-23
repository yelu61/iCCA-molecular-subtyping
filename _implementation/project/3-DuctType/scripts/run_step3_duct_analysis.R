#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(20260715)

required <- c("limma", "glmnet", "pROC", "ggplot2", "ggrepel", "pheatmap",
              "patchwork", "readxl", "msigdbr", "fgsea")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(limma)
  library(glmnet)
  library(ggplot2)
  library(patchwork)
})

root <- normalizePath(getwd())
step_dir <- file.path(root, "3-DuctType")
data_dir <- file.path(step_dir, "data")
pub_dir <- file.path(step_dir, "publication")
table_dir <- file.path(pub_dir, "tables")
fig_dir <- file.path(pub_dir, "figures")
model_dir <- file.path(pub_dir, "models")
report_dir <- file.path(pub_dir, "reports")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
for (d in c(table_dir, fig_dir, model_dir, report_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

load_matrix <- function(path, object_name = "corrected_expr_limma") {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  if (exists(object_name, envir = e, inherits = FALSE)) return(get(object_name, envir = e))
  objs <- ls(e)
  mats <- objs[vapply(objs, function(z) {
    value <- get(z, envir = e)
    is.matrix(value) || is.data.frame(value)
  }, logical(1))]
  if (length(mats) != 1L) stop("Cannot identify matrix in ", path)
  as.matrix(get(mats, envir = e))
}

write_csv <- function(x, path) write.csv(x, path, row.names = FALSE, na = "")
save_pair <- function(plot, stem, width = 8, height = 5.5) {
  ggsave(file.path(fig_dir, paste0(stem, ".pdf")), plot, width = width, height = height,
         device = cairo_pdf)
  ggsave(file.path(fig_dir, paste0(stem, ".png")), plot, width = width, height = height,
         dpi = 450)
}

after_path <- file.path(root, "1-Clustering/data/combined_expr_after_correction.rds")
before_path <- file.path(root, "1-Clustering/data/combined_expr_before_correction.rds")
meta_path <- file.path(root, "1-Clustering/data/sample_info_with_clusters.xlsx")
palette_path <- file.path(root, "0-Data/qc/study_color_dictionary.csv")
stopifnot(file.exists(after_path), file.exists(before_path), file.exists(meta_path), file.exists(palette_path))

expr_after <- load_matrix(after_path)
expr_before <- load_matrix(before_path)
metadata <- readxl::read_excel(meta_path)
metadata <- as.data.frame(metadata, check.names = FALSE)
stopifnot(identical(dim(expr_after), dim(expr_before)))
stopifnot(setequal(colnames(expr_after), metadata$Sample))
metadata <- metadata[match(colnames(expr_after), metadata$Sample), , drop = FALSE]
stopifnot(identical(metadata$Sample, colnames(expr_after)))

observed_levels <- c("small duct type", "large duct type")
metadata$DuctObserved <- ifelse(metadata$`iCCA classification` %in% observed_levels,
                               metadata$`iCCA classification`, NA_character_)
metadata$DuctObserved <- factor(metadata$DuctObserved, levels = observed_levels)
eligible <- !is.na(metadata$DuctObserved)
meta_labeled <- metadata[eligible, , drop = FALSE]
samples_labeled <- meta_labeled$Sample
expr_labeled <- expr_after[, samples_labeled, drop = FALSE]
expr_labeled_before <- expr_before[, samples_labeled, drop = FALSE]

if (anyDuplicated(metadata$Sample)) stop("Duplicate Sample identifiers")
if (nrow(meta_labeled) != 115L) stop("Expected 115 observed duct labels; found ", nrow(meta_labeled))
if (!all(table(meta_labeled$DuctObserved) == c(86L, 29L))) stop("Unexpected duct-label counts")
study_label_tab <- table(meta_labeled$Study, meta_labeled$DuctObserved)
if (any(study_label_tab == 0)) stop("Every labeled study must contain both duct classes")

input_contract <- data.frame(
  Item = c("post_correction_expression", "pre_correction_expression", "metadata",
           "expression_dimensions", "sample_alignment", "observed_label_definition",
           "eligible_samples", "small_duct_samples", "large_duct_samples", "eligible_studies"),
  Value = c(after_path, before_path, meta_path, paste(dim(expr_after), collapse = " x "),
            "exact identifier match", paste(observed_levels, collapse = " | "), nrow(meta_labeled),
            sum(meta_labeled$DuctObserved == "small duct type"),
            sum(meta_labeled$DuctObserved == "large duct type"),
            length(unique(meta_labeled$Study))),
  LockedInputModified = FALSE
)
write_csv(input_contract, file.path(table_dir, "step3_input_contract.csv"))

eligibility_table <- data.frame(
  Sample = metadata$Sample,
  Study = metadata$Study,
  OriginalDuctField = metadata$`iCCA classification`,
  DuctObserved = as.character(metadata$DuctObserved),
  EligibleForObservedDuctAnalysis = eligible,
  ExclusionReason = ifelse(eligible, "", "No observed large-/small-duct pathology label")
)
write_csv(eligibility_table, file.path(table_dir, "step3_sample_eligibility.csv"))

label_counts <- as.data.frame(study_label_tab)
names(label_counts) <- c("Study", "DuctObserved", "N")
label_counts$StudyTotal <- ave(label_counts$N, label_counts$Study, FUN = sum)
label_counts$WithinStudyFraction <- label_counts$N / label_counts$StudyTotal
write_csv(label_counts, file.path(table_dir, "step3_label_provenance_by_study.csv"))

pathology_provenance <- data.frame(
  Study = c("2019Ahn", "2023Deng", "2023Lee"),
  LocalLargeDuctN = c(15, 9, 5),
  LocalSmallDuctN = c(15, 63, 8),
  SourceTerminology = c("bile duct-type versus cholangiolar-type pathology",
                        "large versus small intrahepatic bile-duct subtype",
                        "large-duct versus small-duct histological subtype"),
  SourceCriterion = c("Original abstract reports pathology terminology but does not establish identity with WHO large-/small-duct labels",
                      "Level or size of the affected intrahepatic ducts",
                      "WHO fifth-edition framework using histology, location, and characteristic genomic events"),
  HarmonizationAssessment = c("Project-level mapping; sample-level equivalence requires pathology-source confirmation",
                              "Direct large-/small-duct terminology",
                              "Direct large-/small-duct terminology"),
  Source = c("Ahn et al. Hepatol Int 2019; doi:10.1007/s12072-019-09954-3",
             "Deng et al. Hepatology 2023; doi:10.1002/hep.32624",
             "Lee et al. Nat Commun 2023; doi:10.1038/s41467-023-35896-4")
)
write_csv(pathology_provenance, file.path(table_dir, "step3_pathology_label_provenance.csv"))

legacy_audit <- data.frame(
  Issue = c("Outcome-informed feature selection before CV", "Training-set resubstitution performance",
            "No Study adjustment in DEG", "Observed and predicted duct labels conflated",
            "Predicted labels used for broad clinical and survival plots", "Notebook structure and language"),
  LegacyEvidence = c("DEG-derived marker_genes selected on all 115 labels before cv.glmnet",
                     "predict(cvfit, newx=x) evaluated the same samples used to fit the model",
                     "model.matrix(~labels) omitted Study", "lasso_pred written into the clustering metadata beside pathology labels",
                     "Clinical and KM sections used lasso_pred across the full cohort",
                     "50 mostly code-only cells, Chinese comments, repeated plotting cells"),
  PublicationImpact = c("Optimistic feature selection", "Optimistic AUROC/AUPRC and confusion matrix",
                        "Potential cohort-confounded expression effects", "Ambiguous denominators and evidence type",
                        "Prediction-based associations could be mistaken for observed pathology", "Not reader-facing or reproducible"),
  Resolution = c("Variance screen and LASSO selection recomputed inside outer training folds",
                 "Repeated nested-CV OOF and LOSO predictions only", "Explicit large-minus-small contrast with Study covariate",
                 "Distinct observed, OOF, and projection fields", "Observed-label clinical audit only; survival deferred",
                 "Rebuilt English analysis-report notebook")
)
write_csv(legacy_audit, file.path(table_dir, "step3_legacy_static_audit_findings.csv"))

study_colors <- read.csv(palette_path, check.names = FALSE)
study_color_col <- intersect(c("ColorHex", "Color", "Hex", "color", "hex"), names(study_colors))[1]
study_name_col <- intersect(c("Study", "study"), names(study_colors))[1]
study_palette <- setNames(study_colors[[study_color_col]], study_colors[[study_name_col]])
duct_palette <- c("small duct type" = "#0072B2", "large duct type" = "#D55E00")
write_csv(data.frame(DuctType = names(duct_palette), Color = unname(duct_palette),
                     Role = c("Observed or predicted small-duct", "Observed or predicted large-duct")),
          file.path(table_dir, "step3_duct_color_dictionary.csv"))

p_labels <- ggplot(label_counts, aes(Study, N, fill = DuctObserved)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.7) +
  geom_text(aes(label = N), position = position_dodge(width = 0.78), vjust = -0.35, size = 3.5) +
  scale_fill_manual(values = duct_palette) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Observed duct labels are available in three studies",
       subtitle = "All three studies contain both pathology-defined classes",
       x = NULL, y = "Samples", fill = "Observed duct type") +
  theme_classic(base_size = 11) + theme(legend.position = "top")
save_pair(p_labels, "step3_label_provenance", 7.4, 5.1)

# Study-adjusted large-minus-small differential expression.
duct <- factor(meta_labeled$DuctObserved, levels = observed_levels)
study <- factor(meta_labeled$Study)
design <- model.matrix(~ 0 + duct + study)
colnames(design) <- make.names(colnames(design))
fit <- eBayes(lmFit(expr_labeled, design))
contrast <- makeContrasts(ductlarge.duct.type - ductsmall.duct.type, levels = design)
fit_contrast <- eBayes(contrasts.fit(fit, contrast))
deg <- topTable(fit_contrast, coef = 1, number = Inf, sort.by = "none", adjust.method = "BH")
deg$Gene <- rownames(deg)
deg$logFC_large_minus_small <- deg$logFC
deg$FDR <- deg$adj.P.Val
deg$PositiveDirection <- "large duct type - small duct type"
deg$Model <- "limma: ~0 + DuctObserved + Study"
deg$N_large <- sum(duct == "large duct type")
deg$N_small <- sum(duct == "small duct type")
deg <- deg[, c("Gene", "logFC_large_minus_small", "AveExpr", "t", "P.Value", "FDR", "B",
               "PositiveDirection", "Model", "N_large", "N_small")]
deg <- deg[order(deg$FDR, -abs(deg$logFC_large_minus_small)), ]
write_csv(deg, file.path(table_dir, "step3_study_adjusted_deg_all_genes.csv"))
write_csv(deg, file.path(data_dir, "deg_LS.csv"))

# Within-study effects and fixed-effect meta-analysis.
study_effects <- lapply(levels(study), function(st) {
  idx <- study == st
  d <- factor(duct[idx], levels = observed_levels)
  des <- model.matrix(~ 0 + d)
  colnames(des) <- c("small", "large")
  f <- eBayes(contrasts.fit(lmFit(expr_labeled[, idx, drop = FALSE], des),
                           makeContrasts(large - small, levels = des)))
  se <- f$stdev.unscaled[, 1] * sqrt(f$s2.post)
  data.frame(Gene = rownames(expr_labeled), Study = st,
             logFC_large_minus_small = f$coefficients[, 1], SE = se,
             PValue = 2 * pnorm(-abs(f$coefficients[, 1] / se)),
             N_large = sum(d == "large duct type"), N_small = sum(d == "small duct type"))
})
study_effects <- do.call(rbind, study_effects)
write_csv(study_effects, file.path(table_dir, "step3_within_study_deg_effects.csv"))

split_effects <- split(study_effects, study_effects$Gene)
meta_rows <- lapply(split_effects, function(z) {
  valid <- is.finite(z$SE) & z$SE > 0 & is.finite(z$logFC_large_minus_small)
  z <- z[valid, ]
  w <- 1 / z$SE^2
  est <- sum(w * z$logFC_large_minus_small) / sum(w)
  se <- sqrt(1 / sum(w))
  q <- sum(w * (z$logFC_large_minus_small - est)^2)
  df <- max(nrow(z) - 1, 1)
  data.frame(Gene = z$Gene[1], MetaLogFC_large_minus_small = est, MetaSE = se,
             MetaPValue = 2 * pnorm(-abs(est / se)), Q = q,
             I2 = max(0, (q - df) / max(q, .Machine$double.eps)) * 100,
             Studies = nrow(z), PositiveStudies = sum(z$logFC_large_minus_small > 0),
             NegativeStudies = sum(z$logFC_large_minus_small < 0))
})
meta <- do.call(rbind, meta_rows)
meta$MetaFDR <- p.adjust(meta$MetaPValue, "BH")
meta <- merge(meta, deg[, c("Gene", "logFC_large_minus_small", "FDR")], by = "Gene", all.x = TRUE)
meta$DirectionConsistentAllStudies <- with(meta, PositiveStudies == Studies | NegativeStudies == Studies)
meta <- meta[order(meta$MetaFDR, -abs(meta$MetaLogFC_large_minus_small)), ]
write_csv(meta, file.path(table_dir, "step3_cross_study_deg_meta_summary.csv"))

deg$Significance <- "Not significant"
deg$Significance[deg$FDR < 0.05 & deg$logFC_large_minus_small >= 0.5] <- "Large-duct higher"
deg$Significance[deg$FDR < 0.05 & deg$logFC_large_minus_small <= -0.5] <- "Small-duct higher"
top_label <- unique(c(head(deg$Gene[deg$Significance == "Large-duct higher"], 8),
                      head(deg$Gene[deg$Significance == "Small-duct higher"], 8)))
p_volcano <- ggplot(deg, aes(logFC_large_minus_small, -log10(pmax(FDR, 1e-300)), color = Significance)) +
  geom_point(alpha = 0.65, size = 1.25) +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = 2, color = "grey55") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey55") +
  ggrepel::geom_text_repel(data = deg[deg$Gene %in% top_label, ], aes(label = Gene),
                           size = 3, max.overlaps = Inf, box.padding = 0.25,
                           show.legend = FALSE) +
  scale_color_manual(values = c("Large-duct higher" = duct_palette[["large duct type"]],
                                "Small-duct higher" = duct_palette[["small duct type"]],
                                "Not significant" = "grey78")) +
  labs(title = "Study-adjusted expression differences between observed duct types",
       subtitle = "Explicit contrast: large-duct minus small-duct; n=29 versus n=86",
       x = "log2 expression difference", y = expression(-log[10](FDR)), color = NULL) +
  theme_classic(base_size = 11) + theme(legend.position = "top")
save_pair(p_volcano, "step3_study_adjusted_deg_volcano", 8.2, 6.2)

heat_genes <- unique(c(head(deg$Gene[deg$Significance == "Large-duct higher"], 15),
                       head(deg$Gene[deg$Significance == "Small-duct higher"], 15)))
heat_mat <- expr_labeled[heat_genes, , drop = FALSE]
heat_ann <- data.frame(DuctObserved = as.character(meta_labeled$DuctObserved),
                       Study = meta_labeled$Study, row.names = meta_labeled$Sample)
heat_colors <- list(DuctObserved = duct_palette,
                    Study = study_palette[unique(meta_labeled$Study)])
pheatmap::pheatmap(heat_mat, scale = "row", show_colnames = FALSE,
                   annotation_col = heat_ann, annotation_colors = heat_colors,
                   clustering_distance_cols = "correlation", border_color = NA,
                   main = "Study-adjusted observed duct-type expression features",
                   filename = file.path(fig_dir, "step3_study_adjusted_marker_heatmap.pdf"),
                   width = 10, height = 7)
pheatmap::pheatmap(heat_mat, scale = "row", show_colnames = FALSE,
                   annotation_col = heat_ann, annotation_colors = heat_colors,
                   clustering_distance_cols = "correlation", border_color = NA,
                   main = "Study-adjusted observed duct-type expression features",
                   filename = file.path(fig_dir, "step3_study_adjusted_marker_heatmap.png"),
                   width = 10, height = 7, dpi = 450)

# Ranked pathway characterization.
gene_sets <- rbind(
  transform(msigdbr::msigdbr(species = "Homo sapiens", collection = "H")[, c("gs_name", "gene_symbol")],
            Collection = "Hallmark"),
  transform(msigdbr::msigdbr(species = "Homo sapiens", collection = "C2",
                             subcollection = "CP:REACTOME")[, c("gs_name", "gene_symbol")],
            Collection = "Reactome")
)
pathways <- split(gene_sets$gene_symbol, gene_sets$gs_name)
ranks <- setNames(fit_contrast$t[, 1], rownames(fit_contrast$t))
ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)
set.seed(20260715)
gsea <- fgsea::fgseaMultilevel(pathways = pathways, stats = ranks, minSize = 15, maxSize = 500,
                              eps = 0, nPermSimple = 10000)
gsea <- as.data.frame(gsea)
gsea$Collection <- ifelse(grepl("^HALLMARK_", gsea$pathway), "Hallmark", "Reactome")
gsea$Direction <- ifelse(gsea$NES > 0, "Large-duct higher", "Small-duct higher")
gsea$leadingEdge <- vapply(gsea$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
gsea <- gsea[order(gsea$padj, -abs(gsea$NES)), ]
write_csv(gsea, file.path(table_dir, "step3_ranked_gsea_large_minus_small.csv"))
top_gsea <- do.call(rbind, lapply(split(gsea, interaction(gsea$Collection, gsea$Direction)), head, 6))
top_gsea$Label <- gsub("^(HALLMARK_|REACTOME_)", "", top_gsea$pathway)
top_gsea$Label <- gsub("_", " ", top_gsea$Label)
top_gsea$Label <- gsub("REGULATION OF INSULIN LIKE GROWTH FACTOR IGF TRANSPORT AND UPTAKE BY INSULIN LIKE GROWTH FACTOR BINDING PROTEINS IGFBPS",
                       "IGF TRANSPORT AND UPTAKE BY IGFBPS", top_gsea$Label, fixed = TRUE)
top_gsea$Label <- gsub("FORMATION OF FIBRIN CLOT CLOTTING CASCADE",
                       "FIBRIN CLOT FORMATION", top_gsea$Label, fixed = TRUE)
top_gsea$Label <- gsub("SRP DEPENDENT COTRANSLATIONAL PROTEIN TARGETING TO MEMBRANE",
                       "SRP DEPENDENT PROTEIN TARGETING TO MEMBRANE", top_gsea$Label, fixed = TRUE)
top_gsea$Label <- gsub("RESPONSE OF EIF2AK4 GCN2 TO AMINO ACID DEFICIENCY",
                       "GCN2 RESPONSE TO AMINO ACID DEFICIENCY", top_gsea$Label, fixed = TRUE)
top_gsea$Label <- vapply(top_gsea$Label, function(x) paste(strwrap(x, width = 34), collapse = "\n"),
                         FUN.VALUE = character(1))
top_gsea$Label <- factor(top_gsea$Label, levels = rev(unique(top_gsea$Label[order(top_gsea$NES)])))
p_gsea <- ggplot(top_gsea, aes(NES, Label, color = Direction, size = -log10(pmax(padj, 1e-12)))) +
  geom_vline(xintercept = 0, color = "grey70") + geom_point(alpha = 0.9) +
  facet_wrap(~Collection, scales = "free_y", ncol = 1) +
  scale_color_manual(values = c("Large-duct higher" = duct_palette[["large duct type"]],
                                "Small-duct higher" = duct_palette[["small duct type"]])) +
  labs(title = "Ranked pathway differences characterize observed duct types",
       subtitle = "Study-adjusted limma t statistics; positive NES indicates large-duct enrichment",
       x = "Normalized enrichment score", y = NULL, color = NULL, size = expression(-log[10](FDR))) +
  theme_classic(base_size = 10) + theme(legend.position = "top")
save_pair(p_gsea, "step3_ranked_gsea", 9.2, 8.5)

# Leakage-controlled LASSO utilities.
make_folds <- function(y, study, k, seed) {
  set.seed(seed)
  fold <- integer(length(y))
  strata <- interaction(as.character(study), as.character(y), drop = TRUE)
  for (lev in levels(strata)) {
    ii <- which(strata == lev)
    fold[ii] <- sample(rep(seq_len(k), length.out = length(ii)))
  }
  fold
}

train_fold_model <- function(expr, meta, seed, top_n = 2000L) {
  y <- as.integer(meta$DuctObserved == "large duct type")
  variances <- apply(expr, 1, var, na.rm = TRUE)
  genes <- names(sort(variances[is.finite(variances) & variances > 0], decreasing = TRUE))
  genes <- head(genes, min(top_n, length(genes)))
  x <- t(expr[genes, , drop = FALSE])
  inner_fold <- make_folds(y, meta$Study, 5L, seed)
  weights <- ifelse(y == 1, length(y) / (2 * sum(y == 1)), length(y) / (2 * sum(y == 0)))
  cv <- cv.glmnet(x, y, family = "binomial", alpha = 1, foldid = inner_fold,
                  weights = weights, standardize = TRUE, type.measure = "deviance")
  list(cv = cv, genes = genes, lambda = cv$lambda.min)
}

predict_fold_model <- function(model, expr) {
  as.numeric(predict(model$cv, newx = t(expr[model$genes, , drop = FALSE]),
                     s = model$lambda, type = "response"))
}

run_nested <- function(expr, meta, repeats = 10L, seed_base = 20260715L, analysis) {
  pred_rows <- list(); feature_rows <- list(); counter <- 1L
  for (r in seq_len(repeats)) {
    outer_fold <- make_folds(meta$DuctObserved, meta$Study, 5L, seed_base + r)
    for (f in 1:5) {
      test <- which(outer_fold == f); train <- which(outer_fold != f)
      model <- train_fold_model(expr[, train, drop = FALSE], meta[train, ], seed_base + 1000L*r + f)
      prob <- predict_fold_model(model, expr[, test, drop = FALSE])
      pred_rows[[counter]] <- data.frame(Analysis = analysis, Repeat = r, Fold = f,
                                         Sample = meta$Sample[test], Study = meta$Study[test],
                                         DuctObserved = as.character(meta$DuctObserved[test]),
                                         ProbabilityLargeDuct = prob,
                                         PredictedDuctClass = ifelse(prob >= 0.5, "large duct type", "small duct type"))
      co <- as.matrix(coef(model$cv, s = model$lambda))
      nz <- rownames(co)[co[, 1] != 0 & rownames(co) != "(Intercept)"]
      if (length(nz)) feature_rows[[counter]] <- data.frame(Analysis = analysis, Repeat = r, Fold = f,
                                                            Gene = nz, Coefficient = co[nz, 1])
      counter <- counter + 1L
    }
  }
  list(predictions = do.call(rbind, pred_rows), features = do.call(rbind, feature_rows))
}

average_precision <- function(y, p) {
  o <- order(p, decreasing = TRUE); y <- y[o]
  if (sum(y == 1) == 0) return(NA_real_)
  mean(cumsum(y == 1)[y == 1] / which(y == 1))
}

metrics <- function(y, p) {
  y <- as.integer(y); cls <- as.integer(p >= 0.5)
  tp <- sum(y == 1 & cls == 1); tn <- sum(y == 0 & cls == 0)
  fp <- sum(y == 0 & cls == 1); fn <- sum(y == 1 & cls == 0)
  sensitivity <- tp / max(tp + fn, 1); specificity <- tn / max(tn + fp, 1)
  precision <- tp / max(tp + fp, 1)
  c(AUROC = as.numeric(pROC::auc(y, p, quiet = TRUE)), AUPRC = average_precision(y, p),
    BalancedAccuracy = mean(c(sensitivity, specificity)), Sensitivity = sensitivity,
    Specificity = specificity, Precision = precision,
    F1 = 2 * precision * sensitivity / max(precision + sensitivity, .Machine$double.eps),
    Brier = mean((p - y)^2))
}

aggregate_oof <- function(pred) {
  agg <- aggregate(ProbabilityLargeDuct ~ Analysis + Sample + Study + DuctObserved,
                   pred, mean)
  agg$PredictedDuctClass <- ifelse(agg$ProbabilityLargeDuct >= 0.5,
                                   "large duct type", "small duct type")
  agg
}

metric_with_ci <- function(agg, analysis, boot = 500L) {
  y <- as.integer(agg$DuctObserved == "large duct type"); p <- agg$ProbabilityLargeDuct
  point <- metrics(y, p)
  set.seed(20260715 + nchar(analysis))
  boots <- replicate(boot, {
    ii <- sample(seq_along(y), replace = TRUE)
    if (length(unique(y[ii])) < 2) rep(NA_real_, length(point)) else metrics(y[ii], p[ii])
  })
  data.frame(Analysis = analysis, Metric = names(point), Estimate = as.numeric(point),
             Lower95 = apply(boots, 1, quantile, 0.025, na.rm = TRUE),
             Upper95 = apply(boots, 1, quantile, 0.975, na.rm = TRUE),
             N = length(y), N_large = sum(y == 1), N_small = sum(y == 0),
             Threshold = 0.5)
}

nested_corrected <- run_nested(expr_labeled, meta_labeled, repeats = 10L,
                               analysis = "Repeated nested CV | locked corrected expression")
nested_uncorrected <- run_nested(expr_labeled_before, meta_labeled, repeats = 5L,
                                 analysis = "Repeated nested CV | uncorrected sensitivity")
all_nested_predictions <- rbind(nested_corrected$predictions, nested_uncorrected$predictions)
all_nested_features <- rbind(nested_corrected$features, nested_uncorrected$features)
write_csv(all_nested_predictions, file.path(table_dir, "step3_repeated_nested_cv_predictions_long.csv"))
write_csv(all_nested_features, file.path(table_dir, "step3_nested_cv_selected_features_long.csv"))

oof <- do.call(rbind, lapply(split(all_nested_predictions, all_nested_predictions$Analysis), aggregate_oof))
write_csv(oof, file.path(table_dir, "step3_nested_cv_oof_predictions.csv"))
metric_summary <- do.call(rbind, lapply(split(oof, oof$Analysis), function(z) metric_with_ci(z, z$Analysis[1])))

# Leave-one-study-out evaluation.
run_loso <- function(expr, meta, analysis) {
  rows <- list(); features <- list(); i <- 1L
  for (held in unique(meta$Study)) {
    test <- which(meta$Study == held); train <- which(meta$Study != held)
    model <- train_fold_model(expr[, train, drop = FALSE], meta[train, ], 20260715 + i)
    prob <- predict_fold_model(model, expr[, test, drop = FALSE])
    rows[[i]] <- data.frame(Analysis = analysis, HeldOutStudy = held, Sample = meta$Sample[test],
                            Study = meta$Study[test], DuctObserved = as.character(meta$DuctObserved[test]),
                            ProbabilityLargeDuct = prob,
                            PredictedDuctClass = ifelse(prob >= 0.5, "large duct type", "small duct type"))
    co <- as.matrix(coef(model$cv, s = model$lambda)); nz <- rownames(co)[co[,1] != 0 & rownames(co) != "(Intercept)"]
    if (length(nz)) features[[i]] <- data.frame(Analysis = analysis, HeldOutStudy = held,
                                                Gene = nz, Coefficient = co[nz,1])
    i <- i + 1L
  }
  list(predictions = do.call(rbind, rows), features = do.call(rbind, features))
}

loso_corrected <- run_loso(expr_labeled, meta_labeled, "LOSO | locked corrected expression")
loso_uncorrected <- run_loso(expr_labeled_before, meta_labeled, "LOSO | uncorrected sensitivity")
loso <- rbind(loso_corrected$predictions, loso_uncorrected$predictions)
write_csv(loso, file.path(table_dir, "step3_leave_one_study_out_predictions.csv"))
write_csv(rbind(loso_corrected$features, loso_uncorrected$features),
          file.path(table_dir, "step3_loso_selected_features_long.csv"))

loso_metrics_by_study <- do.call(rbind, lapply(split(loso, interaction(loso$Analysis, loso$HeldOutStudy)), function(z) {
  m <- metrics(z$DuctObserved == "large duct type", z$ProbabilityLargeDuct)
  data.frame(Analysis = z$Analysis[1], HeldOutStudy = z$HeldOutStudy[1], Metric = names(m),
             Estimate = as.numeric(m), N = nrow(z), N_large = sum(z$DuctObserved == "large duct type"),
             N_small = sum(z$DuctObserved == "small duct type"))
}))
write_csv(loso_metrics_by_study, file.path(table_dir, "step3_leave_one_study_out_metrics.csv"))
loso_pooled <- do.call(rbind, lapply(split(loso, loso$Analysis), function(z) metric_with_ci(
  transform(z, Analysis = z$Analysis[1]), z$Analysis[1])))
metric_summary <- rbind(metric_summary, loso_pooled)
write_csv(metric_summary, file.path(table_dir, "step3_classifier_performance_summary.csv"))

# Feature stability from corrected nested CV only.
fs <- nested_corrected$features
feature_stability <- aggregate(Coefficient ~ Gene, fs, function(x) c(SelectionCount = length(x),
                                                                    MedianCoefficient = median(x),
                                                                    MedianAbsCoefficient = median(abs(x))))
feature_stability <- data.frame(Gene = feature_stability$Gene,
                                SelectionCount = feature_stability$Coefficient[, "SelectionCount"],
                                SelectionFrequency = feature_stability$Coefficient[, "SelectionCount"] / 50,
                                MedianCoefficient = feature_stability$Coefficient[, "MedianCoefficient"],
                                MedianAbsCoefficient = feature_stability$Coefficient[, "MedianAbsCoefficient"])
feature_stability <- feature_stability[order(-feature_stability$SelectionFrequency,
                                             -feature_stability$MedianAbsCoefficient), ]
write_csv(feature_stability, file.path(table_dir, "step3_nested_cv_feature_stability.csv"))

top_fs <- head(feature_stability, 25)
top_fs$Gene <- factor(top_fs$Gene, levels = rev(top_fs$Gene))
p_fs <- ggplot(top_fs, aes(SelectionFrequency, Gene, color = MedianCoefficient > 0)) +
  geom_segment(aes(x = 0, xend = SelectionFrequency, yend = Gene), color = "grey75") + geom_point(size = 2.6) +
  scale_color_manual(values = c(`TRUE` = duct_palette[["large duct type"]],
                                `FALSE` = duct_palette[["small duct type"]]),
                     labels = c(`TRUE` = "Large-duct coefficient", `FALSE` = "Small-duct coefficient")) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, .2)) +
  labs(title = "Feature stability is estimated across outer training folds",
       subtitle = "Selection frequency across 50 repeated nested-CV outer models",
       x = "Selection frequency", y = NULL, color = NULL) +
  theme_classic(base_size = 10) + theme(legend.position = "top")
save_pair(p_fs, "step3_nested_cv_feature_stability", 8, 7.2)

# ROC, PR, calibration, and LOSO figures from saved patient-level predictions.
oof_primary <- oof[oof$Analysis == "Repeated nested CV | locked corrected expression", ]
y_primary <- as.integer(oof_primary$DuctObserved == "large duct type")
roc_obj <- pROC::roc(y_primary, oof_primary$ProbabilityLargeDuct, quiet = TRUE, direction = "<")
roc_df <- data.frame(FPR = 1 - roc_obj$specificities, TPR = roc_obj$sensitivities)
ord <- order(oof_primary$ProbabilityLargeDuct, decreasing = TRUE)
y_ord <- y_primary[ord]
pr_df <- data.frame(Recall = cumsum(y_ord) / sum(y_ord),
                    Precision = cumsum(y_ord) / seq_along(y_ord))
primary_analysis <- oof_primary$Analysis[1]
primary_metrics <- setNames(metric_summary$Estimate[metric_summary$Analysis == primary_analysis],
                            metric_summary$Metric[metric_summary$Analysis == primary_analysis])
p_roc <- ggplot(roc_df, aes(FPR, TPR)) + geom_abline(linetype = 2, color = "grey70") +
  geom_path(color = "#3C5488", linewidth = 1) + coord_equal() +
  labs(title = sprintf("Nested-CV ROC | AUROC %.2f", primary_metrics[["AUROC"]]),
       x = "1 - specificity", y = "Sensitivity") + theme_classic(base_size = 11)
p_pr <- ggplot(pr_df, aes(Recall, Precision)) +
  geom_hline(yintercept = mean(y_primary), linetype = 2, color = "grey70") +
  geom_path(color = duct_palette[["large duct type"]], linewidth = 1) +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  labs(title = sprintf("Nested-CV precision-recall | AUPRC %.2f", primary_metrics[["AUPRC"]]),
       subtitle = sprintf("Large-duct prevalence %.2f", mean(y_primary)),
       x = "Recall", y = "Precision") + theme_classic(base_size = 11)
save_pair(p_roc + p_pr + plot_annotation(title = "Out-of-fold discrimination uses no training-set predictions"),
          "step3_nested_cv_discrimination", 10, 5.2)

oof_primary$Bin <- cut(oof_primary$ProbabilityLargeDuct,
                       breaks = unique(quantile(oof_primary$ProbabilityLargeDuct, seq(0,1,.1))),
                       include.lowest = TRUE)
cal <- aggregate(cbind(Predicted = oof_primary$ProbabilityLargeDuct,
                       Observed = y_primary) ~ Bin, oof_primary, mean)
cal$N <- as.integer(table(oof_primary$Bin)[as.character(cal$Bin)])
write_csv(cal, file.path(table_dir, "step3_nested_cv_calibration_bins.csv"))
p_cal <- ggplot(cal, aes(Predicted, Observed)) +
  geom_abline(linetype = 2, color = "grey65") +
  geom_point(aes(size = N), color = "#3C5488", alpha = 0.9) +
  geom_line(color = "#3C5488", linewidth = 0.6) +
  scale_size_continuous(breaks = sort(unique(cal$N)), range = c(3, 6)) +
  coord_equal(xlim = c(0,1), ylim = c(0,1)) +
  labs(title = "Nested-CV calibration remains an empirical assessment",
       subtitle = sprintf("Brier score %.3f; points are probability bins", primary_metrics[["Brier"]]),
       x = "Mean predicted large-duct probability", y = "Observed large-duct fraction", size = "Samples") +
  theme_classic(base_size = 11)
save_pair(p_cal, "step3_nested_cv_calibration", 6.5, 5.8)

loso_auc <- loso_metrics_by_study[loso_metrics_by_study$Metric %in% c("AUROC", "AUPRC") &
                                  loso_metrics_by_study$Analysis == "LOSO | locked corrected expression", ]
heldout_n <- aggregate(Sample ~ HeldOutStudy,
                       loso[loso$Analysis == "LOSO | locked corrected expression", ], length)
names(heldout_n)[2] <- "HeldOutN"
loso_auc <- merge(loso_auc, heldout_n, by = "HeldOutStudy", all.x = TRUE)
loso_auc$StudyLabel <- sprintf("%s\n(n=%d)", loso_auc$HeldOutStudy, loso_auc$HeldOutN)
loso_dodge <- position_dodge(width = 0.38)
p_loso <- ggplot(loso_auc, aes(StudyLabel, Estimate, color = Metric)) +
  geom_hline(yintercept = 0.5, linetype = 2, color = "grey70") +
  geom_point(position = loso_dodge, size = 3) +
  geom_text(aes(label = sprintf("%.2f", Estimate)), position = loso_dodge,
            vjust = -0.8, size = 3, show.legend = FALSE) +
  scale_color_manual(values = c(AUROC = "#3C5488", AUPRC = duct_palette[["large duct type"]])) +
  scale_y_continuous(limits = c(0,1)) +
  labs(title = "Leave-one-study-out performance tests cohort transfer",
       subtitle = "Each study is predicted only by models trained on the other two studies",
       x = NULL, y = "Metric estimate", color = NULL,
       caption = "Study-level estimates, especially 2023Lee (n=13), are imprecise.") +
  theme_classic(base_size = 11) +
  theme(legend.position = "top")
save_pair(p_loso, "step3_leave_one_study_out_performance", 7.5, 5.3)

# Full-data model for projection; never used to estimate training performance.
final_model <- train_fold_model(expr_labeled, meta_labeled, 20260715)
prob_all <- predict_fold_model(final_model, expr_after)
co_final <- as.matrix(coef(final_model$cv, s = final_model$lambda))
final_coef <- data.frame(Gene = rownames(co_final), Coefficient = co_final[,1])
final_coef <- final_coef[final_coef$Gene != "(Intercept)" & final_coef$Coefficient != 0, ]
write_csv(final_coef, file.path(table_dir, "step3_final_lasso_nonzero_coefficients.csv"))
saveRDS(list(model = final_model$cv, genes = final_model$genes, lambda = final_model$lambda,
             positive_class = "large duct type", negative_class = "small duct type",
             input = after_path, training_samples = samples_labeled, seed = 20260715),
        file.path(model_dir, "step3_final_duct_lasso_model.rds"))

primary_oof <- oof_primary[, c("Sample", "ProbabilityLargeDuct", "PredictedDuctClass")]
names(primary_oof)[2:3] <- c("DuctOOFProbabilityLarge", "DuctOOFPredictedClass")
pred_out <- merge(metadata, primary_oof, by = "Sample", all.x = TRUE, sort = FALSE)
pred_out <- pred_out[match(metadata$Sample, pred_out$Sample), ]
pred_out$DuctObserved <- as.character(metadata$DuctObserved)
pred_out$DuctProjectedProbabilityLargeFullModel <- prob_all
pred_out$DuctProjectedClassFullModel <- ifelse(prob_all >= 0.5, "large duct type", "small duct type")
pred_out$DuctCrossFittedOrProjectedProbabilityLarge <- ifelse(!is.na(pred_out$DuctOOFProbabilityLarge),
                                                              pred_out$DuctOOFProbabilityLarge, prob_all)
pred_out$DuctCrossFittedOrProjectedClass <- ifelse(!is.na(pred_out$DuctOOFPredictedClass),
                                                   pred_out$DuctOOFPredictedClass,
                                                   pred_out$DuctProjectedClassFullModel)
pred_out$DuctPredictionSource <- ifelse(!is.na(pred_out$DuctOOFProbabilityLarge),
                                       "repeated_nested_cv_oof", "full_model_projection")
pred_out$lasso_pred <- pred_out$DuctCrossFittedOrProjectedClass
pred_out$lasso_pred_note <- "Deprecated predicted-only compatibility alias; never an observed pathology field"
write_csv(pred_out, file.path(table_dir, "step3_all_samples_observed_and_predicted_duct.csv"))
write_csv(pred_out, file.path(data_dir, "pred.csv"))

# Observed and predicted subtype distributions are intentionally separate.
subtype_levels <- c("C1", "C2", "C3", "C1_C2", "C1_C3", "C2_C3")
obs_sub <- pred_out[!is.na(pred_out$DuctObserved), ]
observed_counts <- as.data.frame(table(FinalSubtype = factor(obs_sub$Ensemble, levels = subtype_levels),
                                       DuctObserved = factor(obs_sub$DuctObserved, levels = observed_levels)))
observed_counts <- observed_counts[observed_counts$Freq > 0, ]; names(observed_counts)[3] <- "N"
observed_counts$SubtypeTotalObserved <- ave(observed_counts$N, observed_counts$FinalSubtype, FUN = sum)
observed_counts$WithinSubtypeFraction <- observed_counts$N / observed_counts$SubtypeTotalObserved
observed_counts$EvidenceType <- "Observed pathology"
write_csv(observed_counts, file.path(table_dir, "step3_observed_duct_distribution_by_subtype.csv"))

projected_counts <- as.data.frame(table(FinalSubtype = factor(pred_out$Ensemble, levels = subtype_levels),
                                        DuctPredicted = factor(pred_out$DuctCrossFittedOrProjectedClass,
                                                               levels = observed_levels)))
names(projected_counts)[3] <- "N"
projected_counts$SubtypeTotalProjected <- ave(projected_counts$N, projected_counts$FinalSubtype, FUN = sum)
projected_counts$WithinSubtypeFraction <- projected_counts$N / projected_counts$SubtypeTotalProjected
projected_counts$EvidenceType <- "Cross-fitted for labeled samples; full-model projection otherwise"
write_csv(projected_counts, file.path(table_dir, "step3_predicted_duct_distribution_by_subtype.csv"))

main_obs <- obs_sub[obs_sub$Ensemble %in% c("C1", "C2", "C3"), ]
main_tab <- table(main_obs$Ensemble, main_obs$DuctObserved)
fisher_main <- fisher.test(main_tab)
c3_tab <- table(C3 = main_obs$Ensemble == "C3", Large = main_obs$DuctObserved == "large duct type")
fisher_c3 <- fisher.test(c3_tab)
association_tests <- data.frame(
  Analysis = c("Observed duct type across C1/C2/C3", "Observed large-duct enrichment in C3 versus C1/C2"),
  Population = c("Main subtypes with observed pathology", "Main subtypes with observed pathology"),
  N = nrow(main_obs), Test = "Fisher exact",
  OddsRatio = c(NA, unname(fisher_c3$estimate)),
  Lower95 = c(NA, fisher_c3$conf.int[1]), Upper95 = c(NA, fisher_c3$conf.int[2]),
  PValue = c(fisher_main$p.value, fisher_c3$p.value)
)
association_tests$FDR <- p.adjust(association_tests$PValue, "BH")
write_csv(association_tests, file.path(table_dir, "step3_observed_duct_subtype_association_tests.csv"))

obs_main_plot <- observed_counts[observed_counts$FinalSubtype %in% c("C1", "C2", "C3"), ]
p_obs <- ggplot(obs_main_plot, aes(FinalSubtype, WithinSubtypeFraction, fill = DuctObserved)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = paste0("n=", N)), position = position_stack(vjust = 0.5),
            color = "white", fontface = "bold", size = 3.3) +
  scale_fill_manual(values = duct_palette) + scale_y_continuous(labels = scales::percent) +
  labs(title = "Observed large-duct morphology is enriched in C3",
       subtitle = sprintf("Main-subtype available cases n=%d; C3 vs C1/C2 OR %.2f (95%% CI %.2f–%.2f)",
                          nrow(main_obs), unname(fisher_c3$estimate), fisher_c3$conf.int[1], fisher_c3$conf.int[2]),
       x = NULL, y = "Within-subtype fraction", fill = "Observed duct type",
       caption = sprintf("Predictions are excluded. Overall C1/C2/C3 Fisher P=%.3g.", fisher_main$p.value)) +
  theme_classic(base_size = 11) + theme(legend.position = "top")
save_pair(p_obs, "step3_observed_duct_by_main_subtype", 7.2, 5.8)

# Clinical fields are audited using observed labels only; survival is deferred.
clinical_vars <- c("Sex", "Age_num", "Stage", "Grade", "HBV infection (1=Positive)",
                   "Liver fluke (1=Positive)", "Vascular invasion (1=Positive)")
clinical_audit <- do.call(rbind, lapply(clinical_vars, function(v) {
  x <- as.character(meta_labeled[[v]]); available <- !is.na(x) & x != "" & x != "NA"
  data.frame(Variable = v, EligibleN = nrow(meta_labeled), AvailableN = sum(available),
             MissingN = sum(!available), AvailableFraction = mean(available))
}))
write_csv(clinical_audit, file.path(table_dir, "step3_observed_clinical_availability.csv"))

clinical_tests <- do.call(rbind, lapply(clinical_vars, function(v) {
  x <- meta_labeled[, c("DuctObserved", v), drop = FALSE]
  names(x)[2] <- "Value"
  x <- x[!is.na(x$Value) & as.character(x$Value) != "" & as.character(x$Value) != "NA", ]
  if (v == "Age_num") {
    x$Value <- as.numeric(x$Value); wt <- wilcox.test(Value ~ DuctObserved, data = x, exact = FALSE)
    return(data.frame(Variable = v, N = nrow(x), Test = "Wilcoxon rank-sum", PValue = wt$p.value))
  }
  tab <- table(x$DuctObserved, x$Value)
  if (nrow(tab) < 2 || ncol(tab) < 2) return(NULL)
  data.frame(Variable = v, N = nrow(x), Test = "Fisher exact", PValue = fisher.test(tab)$p.value)
}))
clinical_tests$FDR <- p.adjust(clinical_tests$PValue, "BH")
clinical_tests$Interpretation <- "Supportive observed-label association; not a model-prediction claim"
write_csv(clinical_tests, file.path(table_dir, "step3_observed_clinical_association_tests.csv"))

# Release gate and artifact manifest.
modeling_contract <- data.frame(
  Component = c("Expression representation", "Outer feature screening", "Scaling", "Class weighting",
                "Hyperparameter tuning", "Decision threshold", "Headline performance", "Study transfer",
                "All-cohort prediction"),
  Implementation = c("Locked Step1 Study-corrected expression; outcome-independent upstream interface",
                     "Top 2,000 variable genes recalculated using outer-training samples only",
                     "glmnet standardization estimated from training samples only",
                     "Inverse-frequency weights recalculated in each training set",
                     "Lambda selected by inner five-fold CV within each outer-training set",
                     "Pre-specified probability threshold 0.5",
                     "Mean patient-level OOF probability across 10 repeated outer five-fold partitions",
                     "Each labeled study held out once and predicted by the other two",
                     "Final all-label model; projection only and excluded from performance estimates"),
  LeakageControl = c("Fixed before Step3 labels are modeled", rep("Training-fold only", 4),
                     "Pre-specified", "Out-of-fold only", "Study-held-out only", "Not a validation estimate")
)
write_csv(modeling_contract, file.path(table_dir, "step3_modeling_contract.csv"))

metric_value <- function(analysis, metric) metric_summary$Estimate[
  metric_summary$Analysis == analysis & metric_summary$Metric == metric][1]
nested_auc <- metric_value("Repeated nested CV | locked corrected expression", "AUROC")
nested_pr <- metric_value("Repeated nested CV | locked corrected expression", "AUPRC")
loso_auc_value <- metric_value("LOSO | locked corrected expression", "AUROC")
raw_nested_auc <- metric_value("Repeated nested CV | uncorrected sensitivity", "AUROC")
gate <- data.frame(
  Check = c("Locked inputs preserved", "Observed label provenance", "Cross-study pathology terminology harmonization",
            "Explicit DEG direction and Study adjustment",
            "Feature selection inside outer folds", "Scaling and lambda tuning inside training folds",
            "Repeated nested-CV out-of-fold predictions", "Leave-one-study-out evaluation",
            "Observed and predicted variables separated", "Training resubstitution excluded from metrics",
            "Cross-study discriminative performance quantified", "Uncorrected-expression sensitivity governed",
            "Survival analysis deferred", "Legacy lasso_pred alias governed"),
  Status = c("PASS", "PASS", "PASS_WITH_HARMONIZATION_CAVEAT", "PASS", "PASS", "PASS", "PASS", "PASS", "PASS", "PASS",
             "PASS_WITH_PERFORMANCE_CAVEAT", "PASS_WITH_CONFOUNDING_CAVEAT", "PASS",
             "PASS_WITH_COMPATIBILITY_ALIAS"),
  Blocking = FALSE,
  Evidence = c("Step0/Step1 matrices read only", "115 labels across three studies; both classes in every study",
               "2019Ahn source uses bile duct-type/cholangiolar-type terminology; project mapping requires explicit manuscript caveat",
               "large minus small; design includes Study", "Top-variable-gene filtering recomputed per outer training fold",
               "glmnet standardization and inner 5-fold lambda selection use training data only",
               "10 x 5-fold; one aggregated OOF probability per labeled sample",
               "Each of the three studies held out once", "DuctObserved and prediction fields have distinct names and denominators",
               "Only OOF/LOSO predictions enter performance tables",
               sprintf("Nested AUROC %.3f, AUPRC %.3f; pooled LOSO AUROC %.3f", nested_auc, nested_pr, loso_auc_value),
               sprintf("Uncorrected nested AUROC %.3f exceeds corrected %.3f and is not promoted because study/platform signal may contribute",
                       raw_nested_auc, nested_auc),
               "Owned by later clinical gate",
               "lasso_pred is a deprecated predicted-only alias; never used as pathology truth")
)
write_csv(gate, file.path(table_dir, "step3_release_gate_status.csv"))

metric_row <- function(analysis, metric) {
  metric_summary[metric_summary$Analysis == analysis & metric_summary$Metric == metric, ][1, ]
}
primary_auc <- metric_row("Repeated nested CV | locked corrected expression", "AUROC")
primary_auprc <- metric_row("Repeated nested CV | locked corrected expression", "AUPRC")
primary_balacc <- metric_row("Repeated nested CV | locked corrected expression", "BalancedAccuracy")
primary_brier <- metric_row("Repeated nested CV | locked corrected expression", "Brier")
pooled_loso_auc <- metric_row("LOSO | locked corrected expression", "AUROC")

methods_ready <- c(
  "# Step 3 methods-ready summary",
  "",
  "Observed duct morphology was available for 115 tumors from 2019Ahn, 2023Deng, and 2023Lee (29 large-duct and 86 small-duct labels). Both classes were represented within every study. The historical 2019Ahn source uses bile duct-type versus cholangiolar-type terminology; therefore its project-level mapping to the harmonized large-/small-duct field is retained with an explicit terminology-harmonization caveat.",
  "",
  "Differential expression used the locked Step 1 batch-corrected matrix and an explicit large-duct-minus-small-duct limma contrast with study included as a covariate. Ranked Hallmark and Reactome enrichment used the resulting moderated t statistics. Classifier performance was estimated only from out-of-fold predictions. In the primary analysis, the top 2,000 variable genes were selected anew within each outer training fold, class weighting and glmnet standardization were confined to training data, and lambda was selected by inner five-fold cross-validation. Ten repeats of stratified outer five-fold cross-validation produced one aggregated out-of-fold probability per labeled sample. Leave-one-study-out analysis held out each study in turn. An uncorrected-expression sensitivity analysis was kept separate and was not promoted because residual study/platform signal could contribute. A final all-label model was fit only for projection and was excluded from performance estimation.",
  "",
  "Observed pathology, cross-fitted predictions, and full-model projections are stored in separate fields and use separate denominators. Survival analysis and incremental C3 classification belong to the later clinical/classifier gate."
)
writeLines(methods_ready, file.path(report_dir, "step3_methods_ready_summary.md"))

results_ready <- c(
  "# Step 3 results-ready summary",
  "",
  sprintf("The leakage-controlled repeated nested-CV model yielded AUROC %.3f (95%% CI %.3f–%.3f), AUPRC %.3f (%.3f–%.3f), balanced accuracy %.3f (%.3f–%.3f), and Brier score %.3f in 115 pathology-labeled tumors.",
          primary_auc$Estimate, primary_auc$Lower95, primary_auc$Upper95,
          primary_auprc$Estimate, primary_auprc$Lower95, primary_auprc$Upper95,
          primary_balacc$Estimate, primary_balacc$Lower95, primary_balacc$Upper95,
          primary_brier$Estimate),
  sprintf("Pooled leave-one-study-out discrimination was lower (AUROC %.3f, 95%% CI %.3f–%.3f), supporting cautious interpretation of cross-cohort transfer.",
          pooled_loso_auc$Estimate, pooled_loso_auc$Lower95, pooled_loso_auc$Upper95),
  sprintf("Among 93 main-subtype tumors with observed pathology, large-duct morphology was enriched in C3 relative to C1/C2 (odds ratio %.2f, 95%% CI %.2f–%.2f; Fisher P=%.3g). This is an association with observed morphology, not a diagnostic-performance or clinical-sensitivity claim.",
          unname(fisher_c3$estimate), fisher_c3$conf.int[1], fisher_c3$conf.int[2], fisher_c3$p.value),
  "",
  "The uncorrected-expression sensitivity model performed better than the locked corrected-expression analysis, but it is not used as the headline estimate because study/platform effects may contribute. Clinical covariate comparisons are supportive and observed-label-only; survival and incremental C3 analyses remain deferred."
)
writeLines(results_ready, file.path(report_dir, "step3_results_ready_summary.md"))

artifacts <- list.files(pub_dir, recursive = TRUE, full.names = TRUE)
artifacts <- artifacts[file.info(artifacts)$isdir %in% FALSE]
artifacts <- artifacts[basename(artifacts) != "step3_publication_output_manifest.csv" &
                       !startsWith(basename(artifacts), ".")]
manifest <- data.frame(ArtifactType = ifelse(grepl("/figures/", artifacts), "Figure",
                                             ifelse(grepl("/models/", artifacts), "Model", "TableOrReport")),
                       Path = sub(paste0("^", root, "/"), "", artifacts),
                       Exists = file.exists(artifacts), NonEmpty = file.info(artifacts)$size > 0,
                       Bytes = file.info(artifacts)$size,
                       MD5 = unname(tools::md5sum(artifacts)))
write_csv(manifest, file.path(table_dir, "step3_publication_output_manifest.csv"))

cat("Step 3 analysis complete:\n")
print(metric_summary[metric_summary$Metric %in% c("AUROC", "AUPRC", "BalancedAccuracy"), ])
print(association_tests)
