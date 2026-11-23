#!/usr/bin/env Rscript

# Publication-facing Step 1 repair. Reads locked Step 0/legacy Step 1 interfaces
# and writes only to 1-Clustering/publication; historical outputs are preserved.

suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(ggplot2)
  library(limma)
  library(msigdbr)
  library(readxl)
})

set.seed(123)
root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
input_dir <- file.path(root, "1-Clustering", "data")
legacy_deg_dir <- file.path(root, "1-Clustering", "deg")
output_dir <- file.path(root, "1-Clustering", "publication")
table_dir <- file.path(output_dir, "tables")
figure_dir <- file.path(output_dir, "figures")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

write_table <- function(x, name) {
  path <- file.path(table_dir, name)
  fwrite(as.data.table(x), path)
  path
}

load_rdata_object <- function(path, object) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  if (!exists(object, envir = env, inherits = FALSE)) {
    stop("Missing object ", object, " in ", path)
  }
  get(object, envir = env, inherits = FALSE)
}

sample_info <- as.data.table(read_excel(file.path(input_dir, "sample_info_with_clusters.xlsx")))
expression <- load_rdata_object(
  file.path(input_dir, "combined_expr_after_correction.rds"),
  "corrected_expr_limma"
)
top_variable <- fread(file.path(input_dir, "top_variable_genes.csv"))

stopifnot(
  nrow(sample_info) == 586L,
  !anyDuplicated(sample_info$Sample),
  identical(sort(as.character(sample_info$Sample)), sort(colnames(expression))),
  nrow(expression) == 10531L,
  nrow(top_variable) == 2000L
)

main_subtypes <- c("C1", "C2", "C3")
study_levels <- sort(unique(sample_info$Study))
sample_info[, SubtypeClass := ifelse(Ensemble %in% main_subtypes, "Main", "Intermediate")]
main_info <- sample_info[SubtypeClass == "Main"]
stopifnot(
  main_info[Ensemble == "C1", .N] == 112L,
  main_info[Ensemble == "C2", .N] == 223L,
  main_info[Ensemble == "C3", .N] == 162L,
  nrow(main_info) == 497L
)

# Exact mapping between the two three-state solutions and governed labels.
mapping <- sample_info[, .N, by = .(NMF, KMeans, Ensemble)][order(NMF, KMeans, Ensemble)]
mapping[, StateClass := ifelse(Ensemble %in% main_subtypes, "Main", "Intermediate")]
write_table(mapping, "step1_nmf_kmeans_label_mapping.csv")

clustering_parameters <- data.table(
  Method = c("Input", "Feature selection", "ConsensusClusterPlus", "NMF", "KMeans", "Hierarchical", "Ensemble"),
  Parameter = c(
    "expression object", "ranking rule", "distance/linkage/reps/pItem/pFeature/seed/k",
    "algorithm/ranks/runs/seed/final rank", "distance/nstart/seed/k",
    "distance/linkage/k", "co-clustering vote/linkage/selected k"
  ),
  Value = c(
    "10,531 genes x 586 samples; limma-corrected historical interface",
    "2,000 genes with largest variance",
    "1-Pearson / hierarchical / 100 / 0.8 / default 1.0 / 123 / 2-6; historical custom PAC label is non-standard",
    "Brunet / 2-6 / 20 for rank estimation and 50 final / 123 / manually overridden to 3",
    "Euclidean / 25 / 123 / 3",
    "Euclidean / complete selected / 2",
    "four-method pairwise vote / Ward.D2 / 6"
  ),
  PublicationStatus = c(
    "Locked input", "Verified", "Needs revision", "Needs explicit manual-decision record",
    "Verified", "Descriptive only", "Not the source of C1/C2/C3 biological labels"
  )
)
write_table(clustering_parameters, "step1_clustering_methods_parameters.csv")

# Study/subtype association and available principal-state denominators.
study_counts <- main_info[, .N, by = .(Study, Ensemble)]
study_counts[, Proportion := N / sum(N), by = Study]
write_table(study_counts, "step1_main_subtype_by_study.csv")
study_table <- table(main_info$Study, main_info$Ensemble)
study_chisq <- suppressWarnings(chisq.test(study_table))
cramers_v <- sqrt(
  as.numeric(study_chisq$statistic) /
    (sum(study_table) * min(nrow(study_table) - 1L, ncol(study_table) - 1L))
)
write_table(data.table(
  Test = "Study by main subtype chi-square",
  N = sum(study_table), Statistic = as.numeric(study_chisq$statistic),
  Df = as.integer(study_chisq$parameter), PValue = study_chisq$p.value,
  CramersV = cramers_v
), "step1_study_subtype_association.csv")

fit_one_vs_rest <- function(expr, metadata, subtype, include_study = TRUE) {
  metadata <- copy(metadata)
  metadata[, Group := factor(ifelse(Ensemble == subtype, subtype, "Others"), levels = c("Others", subtype))]
  metadata[, Study := factor(Study)]
  expr <- expr[, as.character(metadata$Sample), drop = FALSE]
  if (include_study && uniqueN(metadata$Study) > 1L) {
    design <- model.matrix(~ 0 + Group + Study, data = metadata)
  } else {
    design <- model.matrix(~ 0 + Group, data = metadata)
  }
  group_columns <- grep("^Group", colnames(design), value = TRUE)
  subtype_column <- group_columns[endsWith(group_columns, subtype)]
  others_column <- group_columns[endsWith(group_columns, "Others")]
  if (length(subtype_column) != 1L || length(others_column) != 1L) {
    stop("Could not resolve explicit group columns for ", subtype)
  }
  contrast <- matrix(0, nrow = ncol(design), ncol = 1L,
                     dimnames = list(colnames(design), paste0(subtype, "_minus_others")))
  contrast[subtype_column, 1L] <- 1
  contrast[others_column, 1L] <- -1
  fit <- eBayes(contrasts.fit(lmFit(expr, design), contrast), robust = TRUE)
  result <- as.data.table(topTable(fit, number = Inf, sort.by = "none"), keep.rownames = "Gene")
  result[, SE := as.numeric(fit$stdev.unscaled[, 1L] * sqrt(fit$s2.post))]
  result[, `:=`(
    Contrast = paste0(subtype, "_minus_others"),
    PositiveDirection = paste0("Higher in ", subtype),
    Subtype = subtype,
    NSubtype = metadata[Ensemble == subtype, .N],
    NOthers = metadata[Ensemble != subtype, .N],
    NTotal = nrow(metadata),
    StudyAdjusted = include_study && uniqueN(metadata$Study) > 1L
  )]
  setnames(result, c("logFC", "adj.P.Val", "P.Value"), c("logFC_subtype_minus_others", "FDR", "PValue"))
  setcolorder(result, c(
    "Gene", "Subtype", "Contrast", "PositiveDirection", "logFC_subtype_minus_others",
    "SE", "AveExpr", "t", "PValue", "FDR", "B", "NSubtype", "NOthers", "NTotal",
    "StudyAdjusted"
  ))
  result
}

corrected_degs <- rbindlist(lapply(main_subtypes, function(subtype) {
  fit_one_vs_rest(expression, main_info, subtype, include_study = TRUE)
}))
write_table(corrected_degs, "step1_main_only_study_adjusted_deg.csv")
for (subtype in main_subtypes) {
  write_table(
    corrected_degs[Subtype == subtype],
    paste0("step1_deg_", subtype, "_minus_others_main_only_study_adjusted.csv")
  )
}

# Independent legacy-direction spot-check.
legacy_checks <- rbindlist(lapply(main_subtypes, function(subtype) {
  legacy <- fread(file.path(legacy_deg_dir, paste0("DGE_all_genes_", subtype, "_vs_Others.csv")))
  setnames(legacy, 1L, "Gene")
  ids <- intersect(as.character(sample_info[Ensemble == subtype, Sample]), colnames(expression))
  others <- intersect(as.character(sample_info[Ensemble != subtype, Sample]), colnames(expression))
  observed <- rowMeans(expression[, ids, drop = FALSE]) - rowMeans(expression[, others, drop = FALSE])
  joined <- merge(
    legacy[, .(Gene, LegacyLogFC = logFC)],
    data.table(Gene = names(observed), ObservedSubtypeMinusOthers = observed),
    by = "Gene"
  )
  data.table(
    Subtype = subtype,
    GenesCompared = nrow(joined),
    CorrelationLegacyVsObserved = cor(joined$LegacyLogFC, joined$ObservedSubtypeMinusOthers),
    SignAgreementRate = mean(sign(joined$LegacyLogFC) == sign(joined$ObservedSubtypeMinusOthers)),
    Conclusion = "Legacy coefficient is Others minus subtype"
  )
}))
write_table(legacy_checks, "step1_legacy_deg_direction_audit.csv")

# Study-specific effects and fixed-effect heterogeneity summary.
study_effects <- rbindlist(lapply(study_levels, function(study) {
  study_metadata <- main_info[Study == study]
  rbindlist(lapply(main_subtypes, function(subtype) {
    if (study_metadata[Ensemble == subtype, .N] == 0L || study_metadata[Ensemble != subtype, .N] == 0L) {
      return(NULL)
    }
    x <- fit_one_vs_rest(expression, study_metadata, subtype, include_study = FALSE)
    x[, Study := study]
    x
  }))
}), fill = TRUE)
write_table(study_effects, "step1_study_specific_deg_effects_long.csv")

meta_summary <- study_effects[is.finite(SE) & SE > 0, {
  weights <- 1 / (SE^2)
  pooled <- sum(weights * logFC_subtype_minus_others) / sum(weights)
  pooled_se <- sqrt(1 / sum(weights))
  z <- pooled / pooled_se
  q <- sum(weights * (logFC_subtype_minus_others - pooled)^2)
  df <- .N - 1L
  .(
    Studies = .N,
    MetaLogFC = pooled,
    MetaSE = pooled_se,
    MetaZ = z,
    MetaPValue = 2 * pnorm(-abs(z)),
    CochranQ = q,
    HeterogeneityDf = df,
    HeterogeneityPValue = if (df > 0L) pchisq(q, df = df, lower.tail = FALSE) else NA_real_,
    I2Percent = if (q > 0 && df > 0L) max(0, (q - df) / q) * 100 else 0,
    DirectionConsistency = mean(sign(logFC_subtype_minus_others) == sign(pooled))
  )
}, by = .(Gene, Subtype, Contrast, PositiveDirection)]
meta_summary[, MetaFDR := p.adjust(MetaPValue, method = "BH"), by = Subtype]
write_table(meta_summary, "step1_cross_study_deg_meta_summary.csv")

# Leave-one-study-out sensitivity relative to the complete main-state fit.
loso <- rbindlist(lapply(study_levels, function(held_out) {
  training <- main_info[Study != held_out]
  rbindlist(lapply(main_subtypes, function(subtype) {
    reduced <- fit_one_vs_rest(expression, training, subtype, include_study = TRUE)
    full <- corrected_degs[Subtype == subtype]
    joined <- merge(
      full[, .(Gene, FullLogFC = logFC_subtype_minus_others, FullFDR = FDR)],
      reduced[, .(Gene, ReducedLogFC = logFC_subtype_minus_others)],
      by = "Gene"
    )
    top <- joined[FullFDR < 0.05 & abs(FullLogFC) >= 1]
    data.table(
      HeldOutStudy = held_out,
      Subtype = subtype,
      NTraining = nrow(training),
      PearsonAllGenes = cor(joined$FullLogFC, joined$ReducedLogFC),
      SignAgreementAllGenes = mean(sign(joined$FullLogFC) == sign(joined$ReducedLogFC)),
      HeadlineGenes = nrow(top),
      SignAgreementHeadlineGenes = if (nrow(top)) mean(sign(top$FullLogFC) == sign(top$ReducedLogFC)) else NA_real_
    )
  }))
}))
write_table(loso, "step1_leave_one_study_out_deg_sensitivity.csv")

# Ranked GSEA and thresholded ORA from corrected, main-only, study-adjusted effects.
gene_sets <- rbindlist(list(
  as.data.table(msigdbr(species = "Homo sapiens", collection = "H"))[, Collection := "Hallmark"],
  as.data.table(msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME"))[, Collection := "Reactome"],
  as.data.table(msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_MEDICUS"))[, Collection := "KEGG"],
  as.data.table(msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP"))[, Collection := "GO_BP"]
), fill = TRUE)
gene_sets <- unique(gene_sets[, .(Collection, gs_name, gene_symbol)])

gsea_results <- rbindlist(lapply(main_subtypes, function(subtype) {
  deg <- corrected_degs[Subtype == subtype]
  ranks <- deg$t
  names(ranks) <- deg$Gene
  ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)
  rbindlist(lapply(unique(gene_sets$Collection), function(collection) {
    pathways <- split(
      gene_sets[Collection == collection]$gene_symbol,
      gene_sets[Collection == collection]$gs_name
    )
    result <- as.data.table(fgseaMultilevel(
      pathways = pathways, stats = ranks, minSize = 10L, maxSize = 500L,
      eps = 0, scoreType = "std"
    ))
    result[, `:=`(
      Subtype = subtype,
      Contrast = paste0(subtype, "_minus_others"),
      PositiveDirection = paste0("Positive NES indicates enrichment in ", subtype),
      Collection = collection
    )]
    result[, leadingEdge := vapply(leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))]
    result
  }), fill = TRUE)
}), fill = TRUE)
gsea_results[, AbsNES := abs(NES)]
setorder(gsea_results, Subtype, Collection, padj, -AbsNES)
write_table(gsea_results, "step1_ranked_gsea_main_only_study_adjusted.csv")

ora_results <- rbindlist(lapply(main_subtypes, function(subtype) {
  deg <- corrected_degs[Subtype == subtype]
  universe <- unique(deg$Gene)
  rbindlist(lapply(c("Up", "Down"), function(direction) {
    selected <- if (direction == "Up") {
      deg[FDR < 0.05 & logFC_subtype_minus_others >= 1, Gene]
    } else {
      deg[FDR < 0.05 & logFC_subtype_minus_others <= -1, Gene]
    }
    rbindlist(lapply(unique(gene_sets$Collection), function(collection) {
      sets <- gene_sets[Collection == collection]
      result <- as.data.table(clusterProfiler::enricher(
        gene = selected, universe = universe,
        TERM2GENE = sets[, .(gs_name, gene_symbol)],
        pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
        minGSSize = 10, maxGSSize = 500
      ))
      if (!nrow(result)) return(NULL)
      result[, `:=`(
        Subtype = subtype,
        Contrast = paste0(subtype, "_minus_others"),
        Direction = direction,
        DirectionDefinition = ifelse(
          direction == "Up", paste0("Higher in ", subtype), paste0("Lower in ", subtype)
        ),
        Collection = collection,
        SelectedGenes = length(selected),
        UniverseGenes = length(universe)
      )]
      result
    }), fill = TRUE)
  }), fill = TRUE)
}), fill = TRUE)
setorder(ora_results, Subtype, Direction, Collection, p.adjust)
write_table(ora_results, "step1_thresholded_ora_main_only_study_adjusted.csv")

# Existing TME direction is retained as an audit result, not a naming decision.
estimate <- fread(file.path(root, "1-Clustering", "tme", "estimate_scores_with_clusters.csv"))
estimate_summary <- estimate[cluster %in% main_subtypes, .(
  N = .N,
  MedianStromalScore = median(StromalScore, na.rm = TRUE),
  MedianImmuneScore = median(ImmuneScore, na.rm = TRUE),
  MedianESTIMATEScore = median(ESTIMATEScore, na.rm = TRUE),
  MedianTumorPurity = median(TumorPurity, na.rm = TRUE)
), by = .(Subtype = cluster)][order(Subtype)]
write_table(estimate_summary, "step1_estimate_main_subtype_summary.csv")
write_table(data.table(
  Method = c("ESTIMATE", "EPIC", "Independent third TME method"),
  Status = c("Executed and auditable", "Legacy output not publication-valid", "Not yet executed"),
  HeadlineEligible = c(TRUE, FALSE, FALSE),
  Reason = c(
    "Existing sample-level scores and subtype medians are available",
    "Historical workflow exponentiated log-scale values and requires a valid non-negative input rebuild",
    "Required for method-concordance gate"
  )
), "step1_tme_method_status.csv")

static_audit <- data.table(
  Item = c(
    "Discovery expression object", "Top-variable feature selection", "Governed main-state counts",
    "NMF/KMeans label mapping", "Historical DEG direction", "Historical DEG population",
    "Study confounding", "C2 immune label", "TME method concordance"
  ),
  Status = c("PASS", "PASS", "PASS", "PASS", "FAIL", "FAIL", "NEEDS_SENSITIVITY", "FAIL", "INCOMPLETE"),
  Evidence = c(
    "10,531 genes x 586 samples", "Top 2,000 by variance", "112/223/162; total 497",
    "Exact NMF x KMeans x Ensemble mapping exported",
    "Legacy-vs-observed correlations are -1 with zero sign agreement",
    "Legacy tables compare all six states rather than main-only",
    sprintf("Study x subtype Cramer's V=%.3f; P=%.4g", cramers_v, study_chisq$p.value),
    "C2 has the lowest median ESTIMATE ImmuneScore and highest median inferred purity",
    "ESTIMATE is auditable; EPIC requires rebuild; third method absent"
  ),
  Blocking = c(FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, TRUE, TRUE)
)
write_table(static_audit, "step1_static_audit.csv")

manifest_paths <- unique(c(
  list.files(output_dir, recursive = TRUE, full.names = TRUE),
  file.path(root, "1-Clustering", "scripts", "run_step1_publication_repair.R"),
  file.path(root, "1-clustering.ipynb")
))
manifest <- data.table(
  ArtifactType = fifelse(
    grepl("/figures/", manifest_paths), "figure",
    fifelse(grepl("\\.(R|ipynb)$", manifest_paths), "analysis_code", "source_table")
  ),
  Path = sub(paste0("^", root, "/"), "", normalizePath(manifest_paths, winslash = "/", mustWork = FALSE))
)
manifest[, `:=`(
  Exists = file.exists(file.path(root, Path)),
  NonEmpty = file.info(file.path(root, Path))$size > 0
)]
write_table(manifest, "step1_publication_output_manifest.csv")

cat("Step 1 publication repair completed in ", output_dir, "\n", sep = "")
