#!/usr/bin/env Rscript

# Reproducible, study-aware and subtype-symmetric WGCNA refresh for Step2.
#
# The largest single discovery study (2022Dong) defines signed co-expression
# modules without using subtype labels. Fixed reference PC1 loadings are then
# projected into every study. C1, C2 and C3 use identical association,
# pairwise-specificity, outlier, meta-analysis and leave-one-study-out rules.
# Results support co-expression programs; they do not independently validate a
# lineage, mechanism, treatment response or subtype suffix.

options(stringsAsFactors = FALSE, warn = 1)

root <- normalizePath(getwd(), mustWork = TRUE)
.libPaths(c(file.path(root, ".r-library"), .libPaths()))

suppressPackageStartupMessages({
  library(WGCNA)
  library(limma)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(msigdbr)
  library(metafor)
  library(ggplot2)
  library(pheatmap)
})

set.seed(123)

preservation_permutations_target <- 200L
test_permutations <- Sys.getenv("SCLUCID_WGCNA_TEST_PERMUTATIONS", unset = "")
preservation_permutations <- if (nzchar(test_permutations)) as.integer(test_permutations) else preservation_permutations_target
if (!is.finite(preservation_permutations) || preservation_permutations < 2L) {
  stop("SCLUCID_WGCNA_TEST_PERMUTATIONS must be an integer >=2")
}
test_mode <- preservation_permutations != preservation_permutations_target

module_dir <- file.path(root, "2-TranscriptomicCharacterization")
output_root <- if (test_mode) file.path(tempdir(), "step2-wgcna-test") else module_dir
table_dir <- if (test_mode) file.path(output_root, "tables") else file.path(module_dir, "publication", "tables")
figure_dir <- if (test_mode) file.path(output_root, "figures") else file.path(module_dir, "publication", "figures")
model_dir <- if (test_mode) file.path(output_root, "models") else file.path(module_dir, "publication", "models")
report_dir <- if (test_mode) file.path(output_root, "reports") else file.path(module_dir, "reports")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(model_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)

expr_path <- file.path(root, "1-Clustering", "data", "combined_expr_after_correction.rds")
master_path <- file.path(root, "governance", "data", "subtype_master_table.csv")
deg_path <- file.path(root, "1-Clustering", "publication", "tables", "step1_main_only_study_adjusted_deg.csv")
meta_deg_path <- file.path(root, "1-Clustering", "publication", "tables", "step1_cross_study_deg_meta_summary.csv")
protein_path <- file.path(root, "4-Multiomic", "publication", "tables", "step4_mass_spec_study_adjusted_deg_main_subtypes.csv")
signature_path <- file.path(root, "5-ClassifierValidation", "selected_signature_genes.rds")
palette_path <- file.path(root, "governance", "data", "subtype_label_color_dictionary.csv")

required_inputs <- c(expr_path, master_path, deg_path, meta_deg_path, protein_path,
                     signature_path, palette_path)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) stop("Missing required inputs: ", paste(missing_inputs, collapse = ", "))

reference_study <- "2022Dong"
expected_main_n <- 497L
expected_reference_n <- 223L
max_genes <- 8000L
mad_quantile <- 0.25
power_candidates <- 1:20
scale_free_target <- 0.80
scale_free_minimum <- 0.70
min_module_size <- 30L
merge_cut_height <- 0.20
deep_split <- 2L
max_block_size <- 5000L
preservation_max_module_size <- max_block_size
preservation_max_gold_module_size <- max_block_size
threads <- 2L
effect_threshold <- 0.30
shared_margin <- 0.30
reference_outlier_z_threshold <- -3
reference_sensitivity_ari_min <- 0.75
reference_sensitivity_coverage_min <- 0.75
reference_sensitivity_non_grey_min <- 0.75
reference_specific_module_coverage_min <- 0.50

load_saved_object <- function(path, expected = NULL) {
  env <- new.env(parent = emptyenv())
  object_names <- load(path, envir = env)
  if (!is.null(expected) && expected %in% object_names) return(env[[expected]])
  if (length(object_names) != 1L) {
    stop("Expected one saved object in ", path, "; found: ", paste(object_names, collapse = ", "))
  }
  env[[object_names[[1]]]]
}

write_table <- function(x, filename) {
  path <- file.path(table_dir, filename)
  temporary <- tempfile(pattern = paste0(".", filename, "."), tmpdir = table_dir)
  on.exit(unlink(temporary), add = TRUE)
  write.csv(x, temporary, row.names = FALSE, na = "")
  if (!file.rename(temporary, path)) stop("Atomic replacement failed: ", path)
  path
}

atomic_ggsave <- function(path, plot, ...) {
  extension <- tools::file_ext(path)
  temporary <- tempfile(pattern = paste0(".", basename(path), "."),
                        tmpdir = dirname(path), fileext = paste0(".", extension))
  on.exit(unlink(temporary), add = TRUE)
  ggsave(temporary, plot = plot, ...)
  if (!file.rename(temporary, path)) stop("Atomic replacement failed: ", path)
  invisible(path)
}

atomic_save_rds <- function(object, path, ...) {
  temporary <- tempfile(pattern = paste0(".", basename(path), "."), tmpdir = dirname(path))
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary, ...)
  if (!file.rename(temporary, path)) stop("Atomic replacement failed: ", path)
  invisible(path)
}

object_md5 <- function(object) {
  temporary <- tempfile(fileext = ".rds")
  on.exit(unlink(temporary), add = TRUE)
  saveRDS(object, temporary, version = 3)
  md5(temporary)
}

md5 <- function(path) unname(tools::md5sum(path))
bh <- function(p) p.adjust(p, method = "BH")

expr <- as.matrix(load_saved_object(expr_path, "corrected_expr_limma"))
storage.mode(expr) <- "double"
master <- read.csv(master_path, check.names = FALSE, colClasses = "character")
master <- master[master$SubtypeClass == "Main" & master$FinalSubtype %in% c("C1", "C2", "C3"),
                 c("Sample", "FinalSubtype", "Study")]
master$Sample <- as.character(master$Sample)
master$FinalSubtype <- factor(master$FinalSubtype, levels = c("C1", "C2", "C3"))

common <- master$Sample[master$Sample %in% colnames(expr)]
if (length(common) != expected_main_n || length(unique(common)) != expected_main_n) {
  stop("Expected ", expected_main_n, " unique main-state samples; found ", length(unique(common)))
}
master <- master[match(common, master$Sample), ]
expr <- expr[, common, drop = FALSE]
if (!identical(colnames(expr), master$Sample)) stop("Expression and governance sample order do not match")
if (any(!is.finite(expr))) stop("Expression matrix contains non-finite values")

study_counts <- as.data.frame.matrix(table(master$Study, master$FinalSubtype))
study_counts$Study <- rownames(study_counts)
rownames(study_counts) <- NULL
study_counts$Total <- rowSums(study_counts[c("C1", "C2", "C3")])
study_counts <- study_counts[, c("Study", "C1", "C2", "C3", "Total")]
write_table(study_counts, "step2_wgcna_study_subtype_counts.csv")

largest_study <- study_counts$Study[which.max(study_counts$Total)]
if (largest_study != reference_study) {
  stop("Prespecified reference study is no longer the largest: expected ", reference_study,
       ", observed ", largest_study)
}
reference_samples <- master$Sample[master$Study == reference_study]
validation_samples <- master$Sample[master$Study != reference_study]
if (length(reference_samples) != expected_reference_n) {
  stop("Expected ", expected_reference_n, " reference samples; found ", length(reference_samples))
}

gene_qc <- data.frame(
  Gene = rownames(expr),
  MissingFraction = rowMeans(!is.finite(expr)),
  SD_reference = apply(expr[, reference_samples, drop = FALSE], 1, sd),
  MAD_reference = apply(expr[, reference_samples, drop = FALSE], 1, mad),
  stringsAsFactors = FALSE
)
gene_qc$Eligible <- gene_qc$MissingFraction == 0 & is.finite(gene_qc$SD_reference) &
  gene_qc$SD_reference > 0 & is.finite(gene_qc$MAD_reference) & gene_qc$MAD_reference > 0
mad_cutoff <- unname(quantile(gene_qc$MAD_reference[gene_qc$Eligible], mad_quantile, na.rm = TRUE))
gene_qc$PassMAD <- gene_qc$Eligible & gene_qc$MAD_reference >= mad_cutoff
candidate_genes <- gene_qc$Gene[gene_qc$PassMAD]
if (length(candidate_genes) > max_genes) {
  candidate_genes <- gene_qc$Gene[order(gene_qc$MAD_reference, decreasing = TRUE)][seq_len(max_genes)]
}
gene_qc$IncludedInNetwork <- gene_qc$Gene %in% candidate_genes
gene_qc <- gene_qc[order(gene_qc$IncludedInNetwork, gene_qc$MAD_reference, decreasing = TRUE), ]
if (length(candidate_genes) < 5000L) stop("Fewer than 5,000 genes passed phenotype-independent filtering")

dat_expr_all <- t(expr[candidate_genes, master$Sample, drop = FALSE])
dat_expr_reference <- dat_expr_all[reference_samples, , drop = FALSE]
dat_expr_validation <- dat_expr_all[validation_samples, , drop = FALSE]

gsg <- goodSamplesGenes(dat_expr_reference, verbose = 0)
if (!gsg$allOK) {
  if (any(!gsg$goodSamples)) stop("Reference discovery contains samples failing goodSamplesGenes")
  candidate_genes <- candidate_genes[gsg$goodGenes]
  dat_expr_all <- dat_expr_all[, candidate_genes, drop = FALSE]
  dat_expr_reference <- dat_expr_reference[, candidate_genes, drop = FALSE]
  dat_expr_validation <- dat_expr_validation[, candidate_genes, drop = FALSE]
}
gsg_all <- goodSamplesGenes(dat_expr_all, verbose = 0)
if (!gsg_all$allOK) stop("All-study matrix contains samples or genes failing goodSamplesGenes")
gene_qc$IncludedInNetwork <- gene_qc$Gene %in% candidate_genes
write_table(gene_qc, "step2_wgcna_gene_eligibility.csv")

# Connectivity is standardized within study so a study shift cannot itself
# create a sample-outlier flag. No subtype label is used in this calculation.
sample_audit <- do.call(rbind, lapply(unique(master$Study), function(study) {
  samples <- master$Sample[master$Study == study]
  similarity <- bicor(t(dat_expr_all[samples, , drop = FALSE]),
                      use = "pairwise.complete.obs", maxPOutliers = 0.05)
  adjacency <- (1 + similarity) / 2
  diag(adjacency) <- 0
  connectivity <- rowSums(adjacency)
  connectivity_z <- as.numeric(scale(connectivity))
  data.frame(
    Sample = samples,
    FinalSubtype = as.character(master$FinalSubtype[match(samples, master$Sample)]),
    Study = study,
    NetworkRole = ifelse(study == reference_study, "Reference discovery",
                         "Held-out-study replication"),
    SampleConnectivityWithinStudy = connectivity,
    ConnectivityZWithinStudy = connectivity_z,
    PrespecifiedOutlierFlag = connectivity_z < reference_outlier_z_threshold,
    Excluded = FALSE,
    stringsAsFactors = FALSE
  )
}))
sample_audit$Excluded <- sample_audit$PrespecifiedOutlierFlag &
  sample_audit$NetworkRole == "Reference discovery"
write_table(sample_audit, "step2_wgcna_sample_audit.csv")
flagged_reference_samples <- sample_audit$Sample[
  sample_audit$PrespecifiedOutlierFlag & sample_audit$NetworkRole == "Reference discovery"
]
reference_network_samples <- setdiff(reference_samples, flagged_reference_samples)
dat_expr_reference_network <- dat_expr_reference[reference_network_samples, , drop = FALSE]
gsg_reference_network <- goodSamplesGenes(dat_expr_reference_network, verbose = 0)
if (!gsg_reference_network$allOK) {
  stop("QC-pass reference discovery network fails goodSamplesGenes")
}
gc()

allowWGCNAThreads(nThreads = threads)
sft <- pickSoftThreshold(
  dat_expr_reference_network,
  powerVector = power_candidates,
  RsquaredCut = scale_free_target,
  networkType = "signed",
  corFnc = "bicor",
  corOptions = list(use = "pairwise.complete.obs", maxPOutliers = 0.05),
  verbose = 3
)
soft_table <- as.data.frame(sft$fitIndices)
soft_table$SignedR2 <- -sign(soft_table$slope) * soft_table$SFT.R.sq
target_rows <- which(soft_table$SignedR2 >= scale_free_target & soft_table$mean.k. >= 5)
if (length(target_rows)) {
  selected_row <- target_rows[which.min(soft_table$Power[target_rows])]
  power_rule <- "lowest power with signed R2 >= 0.80 and mean connectivity >= 5"
} else {
  viable <- which(soft_table$mean.k. >= 1)
  if (!length(viable)) viable <- seq_len(nrow(soft_table))
  selected_row <- viable[which.max(soft_table$SignedR2[viable])]
  power_rule <- "fallback: maximum signed R2 among powers with mean connectivity >= 1"
}
soft_power <- soft_table$Power[selected_row]
soft_table$Selected <- soft_table$Power == soft_power
soft_table$SelectionRule <- power_rule
write_table(soft_table, "step2_wgcna_soft_threshold_diagnostics.csv")

soft_plot <- ggplot(soft_table, aes(Power, SignedR2)) +
  geom_hline(yintercept = scale_free_target, linetype = 2, color = "grey45") +
  geom_line(color = "#3B6FB6", linewidth = 0.7) +
  geom_point(aes(color = Selected), size = 2.3) +
  geom_text(aes(label = Power), nudge_y = 0.025, size = 3) +
  scale_color_manual(values = c(`FALSE` = "#3B6FB6", `TRUE` = "#D55E00"), guide = "none") +
  coord_cartesian(ylim = c(min(-0.05, min(soft_table$SignedR2)), 1)) +
  labs(x = "Soft-threshold power", y = "Signed scale-free topology fit", title = "WGCNA soft-threshold selection",
       subtitle = paste0("Reference: ", reference_study, "; selected power = ", soft_power)) +
  theme_bw(base_size = 11)
atomic_ggsave(file.path(figure_dir, "step2_wgcna_soft_threshold.pdf"), soft_plot, width = 7, height = 5)
atomic_ggsave(file.path(figure_dir, "step2_wgcna_soft_threshold.png"), soft_plot, width = 7, height = 5, dpi = 300)

net <- blockwiseModules(
  dat_expr_reference_network,
  maxBlockSize = max_block_size,
  randomSeed = 123,
  corType = "bicor",
  maxPOutliers = 0.05,
  power = soft_power,
  networkType = "signed",
  TOMType = "signed",
  deepSplit = deep_split,
  minModuleSize = min_module_size,
  mergeCutHeight = merge_cut_height,
  reassignThreshold = 0,
  pamRespectsDendro = FALSE,
  numericLabels = TRUE,
  saveTOMs = FALSE,
  nThreads = threads,
  verbose = 3
)

module_labels <- net$colors
names(module_labels) <- colnames(dat_expr_reference_network)
module_colors <- labels2colors(module_labels)
names(module_colors) <- names(module_labels)

# Discovery-cohort connectivity flags are excluded from primary module
# construction. Reintroduce them in an all-reference sensitivity network to
# quantify why the phenotype-blind QC exclusion matters.
adjusted_rand_index <- function(x, y) {
  tab <- table(x, y)
  choose2 <- function(z) z * (z - 1) / 2
  observed <- sum(choose2(tab))
  row_pairs <- sum(choose2(rowSums(tab)))
  column_pairs <- sum(choose2(colSums(tab)))
  total_pairs <- choose2(sum(tab))
  expected <- row_pairs * column_pairs / total_pairs
  maximum <- (row_pairs + column_pairs) / 2
  if (maximum == expected) return(1)
  (observed - expected) / (maximum - expected)
}

if (length(flagged_reference_samples)) {
  sensitivity_reference_samples <- reference_samples
  sensitivity_net <- blockwiseModules(
    dat_expr_reference[sensitivity_reference_samples, , drop = FALSE],
    maxBlockSize = max_block_size,
    randomSeed = 123,
    corType = "bicor",
    maxPOutliers = 0.05,
    power = soft_power,
    networkType = "signed",
    TOMType = "signed",
    deepSplit = deep_split,
    minModuleSize = min_module_size,
    mergeCutHeight = merge_cut_height,
    reassignThreshold = 0,
    pamRespectsDendro = FALSE,
    numericLabels = TRUE,
    saveTOMs = FALSE,
    nThreads = threads,
    verbose = 2
  )
  sensitivity_module_colors <- labels2colors(sensitivity_net$colors)
  names(sensitivity_module_colors) <- colnames(dat_expr_reference)
} else {
  sensitivity_reference_samples <- reference_network_samples
  sensitivity_module_colors <- module_colors
}

primary_module_colors <- sort(unique(module_colors[module_colors != "grey"]))
reference_outlier_module_mapping <- do.call(rbind, lapply(primary_module_colors, function(color) {
  genes <- names(module_colors)[module_colors == color]
  overlaps <- table(sensitivity_module_colors[genes])
  overlaps <- overlaps[names(overlaps) != "grey"]
  if (!length(overlaps)) {
    return(data.frame(
      Module = paste0("ME", color), SensitivityModule = NA_character_,
      PrimaryModuleSize = length(genes), SensitivityModuleSize = 0L,
      BestOverlap = 0L, PrimaryCoverage = 0, Jaccard = 0,
      HypergeometricP = 1, stringsAsFactors = FALSE
    ))
  }
  best_color <- names(overlaps)[which.max(overlaps)]
  best_overlap <- unname(max(overlaps))
  sensitivity_size <- sum(sensitivity_module_colors == best_color)
  universe_size <- length(module_colors)
  data.frame(
    Module = paste0("ME", color), SensitivityModule = paste0("ME", best_color),
    PrimaryModuleSize = length(genes), SensitivityModuleSize = sensitivity_size,
    BestOverlap = best_overlap, PrimaryCoverage = best_overlap / length(genes),
    Jaccard = best_overlap / (length(genes) + sensitivity_size - best_overlap),
    HypergeometricP = phyper(best_overlap - 1, sensitivity_size,
                            universe_size - sensitivity_size, length(genes),
                            lower.tail = FALSE),
    stringsAsFactors = FALSE
  )
}))
reference_outlier_module_mapping$HypergeometricFDR <- bh(reference_outlier_module_mapping$HypergeometricP)
write_table(reference_outlier_module_mapping,
            "step2_wgcna_reference_outlier_module_mapping.csv")

primary_non_grey <- module_colors != "grey"
reference_outlier_network_sensitivity <- data.frame(
  ReferenceSamplesPrimary = length(reference_network_samples),
  ReferenceSamplesSensitivity = length(sensitivity_reference_samples),
  FlaggedReferenceSamples = length(flagged_reference_samples),
  FlaggedSampleIDs = paste(flagged_reference_samples, collapse = ";"),
  AdjustedRandIndexPrimaryNonGrey = adjusted_rand_index(
    module_colors[primary_non_grey], sensitivity_module_colors[primary_non_grey]
  ),
  WeightedBestMatchCoverage = sum(reference_outlier_module_mapping$BestOverlap) /
    sum(reference_outlier_module_mapping$PrimaryModuleSize),
  PrimaryNonGreyRetainedFraction = mean(sensitivity_module_colors[primary_non_grey] != "grey"),
  ARIThreshold = reference_sensitivity_ari_min,
  CoverageThreshold = reference_sensitivity_coverage_min,
  NonGreyRetentionThreshold = reference_sensitivity_non_grey_min,
  stringsAsFactors = FALSE
)
reference_outlier_network_sensitivity$NetworkSensitivityPass <- with(
  reference_outlier_network_sensitivity,
  AdjustedRandIndexPrimaryNonGrey >= ARIThreshold &
    WeightedBestMatchCoverage >= CoverageThreshold &
    PrimaryNonGreyRetainedFraction >= NonGreyRetentionThreshold
)
write_table(reference_outlier_network_sensitivity,
            "step2_wgcna_reference_outlier_network_sensitivity.csv")

# Learn one PC1 loading vector per module in the reference study and reuse it
# unchanged in every cohort. This avoids comparing separately refitted PCs.
project_reference_eigengenes <- function(reference_dat, all_dat, colors) {
  module_names <- paste0("ME", sort(unique(colors[colors != "grey"])))
  scores <- matrix(NA_real_, nrow = nrow(all_dat), ncol = length(module_names),
                   dimnames = list(rownames(all_dat), module_names))
  loading_rows <- list()
  quality_rows <- list()
  for (module_name in module_names) {
    color <- sub("^ME", "", module_name)
    genes <- intersect(names(colors)[colors == color], colnames(reference_dat))
    centers <- colMeans(reference_dat[, genes, drop = FALSE])
    scales <- apply(reference_dat[, genes, drop = FALSE], 2, sd)
    if (any(!is.finite(scales) | scales <= 0)) stop("Invalid reference SD in ", module_name)
    z_reference <- sweep(reference_dat[, genes, drop = FALSE], 2, centers, "-")
    z_reference <- sweep(z_reference, 2, scales, "/")
    pc <- prcomp(z_reference, center = FALSE, scale. = FALSE, rank. = 1)
    loading <- pc$rotation[, 1]
    reference_score <- as.numeric(z_reference %*% loading)
    if (bicor(reference_score, rowMeans(z_reference), use = "pairwise.complete.obs",
              maxPOutliers = 0.05) < 0) loading <- -loading
    z_all <- sweep(all_dat[, genes, drop = FALSE], 2, centers, "-")
    z_all <- sweep(z_all, 2, scales, "/")
    scores[, module_name] <- as.numeric(z_all %*% loading)
    loading_rows[[module_name]] <- data.frame(
      Gene = genes, Module = module_name, ModuleColor = color,
      ReferenceCenter = centers[genes], ReferenceScale = scales[genes],
      PC1Loading = loading[genes], stringsAsFactors = FALSE
    )
    quality_rows[[module_name]] <- data.frame(
      Module = module_name,
      ReferencePC1VarianceExplained = pc$sdev[1]^2 / sum(pc$sdev^2),
      stringsAsFactors = FALSE
    )
  }
  list(scores = scores, loadings = do.call(rbind, loading_rows),
       quality = do.call(rbind, quality_rows))
}

standardize_scores_within_study <- function(scores, metadata) {
  result <- scores
  for (study in unique(metadata$Study)) {
    samples <- metadata$Sample[metadata$Study == study]
    z <- scale(scores[samples, , drop = FALSE])
    if (any(!is.finite(z))) stop("Non-finite projected module score in study ", study)
    result[samples, ] <- z
  }
  result
}

projection <- project_reference_eigengenes(dat_expr_reference_network, dat_expr_all, module_colors)
write_table(projection$loadings, "step2_wgcna_reference_eigengene_loadings.csv")
scores_all <- standardize_scores_within_study(projection$scores, master)
mes_reference <- scores_all[reference_network_samples, , drop = FALSE]
mes_reference_with_flags <- scores_all[reference_samples, , drop = FALSE]
mes_validation <- scores_all[validation_samples, , drop = FALSE]

association_table <- function(scores, metadata, dataset_name, adjust_study) {
  metadata <- metadata[match(rownames(scores), metadata$Sample), ]
  stopifnot(identical(rownames(scores), metadata$Sample))
  metadata$Subtype <- factor(metadata$FinalSubtype, levels = c("C1", "C2", "C3"))
  if (adjust_study) {
    metadata$StudyFactor <- factor(metadata$Study)
    design <- model.matrix(~ 0 + Subtype + StudyFactor, data = metadata)
  } else {
    design <- model.matrix(~ 0 + Subtype, data = metadata)
  }
  colnames(design) <- make.names(colnames(design))
  contrast_matrix <- makeContrasts(
    C1 = SubtypeC1 - (SubtypeC2 + SubtypeC3) / 2,
    C2 = SubtypeC2 - (SubtypeC1 + SubtypeC3) / 2,
    C3 = SubtypeC3 - (SubtypeC1 + SubtypeC2) / 2,
    levels = design
  )
  fit <- eBayes(contrasts.fit(lmFit(t(scores), design), contrast_matrix))
  rows <- do.call(rbind, lapply(seq_len(ncol(contrast_matrix)), function(j) {
    data.frame(
      Dataset = dataset_name,
      Module = rownames(fit$coefficients),
      Subtype = colnames(contrast_matrix)[j],
      Contrast = paste0(colnames(contrast_matrix)[j], "_minus_mean_other_two"),
      N = nrow(scores),
      Studies = length(unique(metadata$Study)),
      EffectSD = fit$coefficients[, j],
      SE = fit$stdev.unscaled[, j] * sqrt(fit$s2.post),
      T = fit$t[, j],
      PValue = fit$p.value[, j],
      stringsAsFactors = FALSE
    )
  }))
  rows$FDR_within_dataset_subtype <- ave(rows$PValue, rows$Dataset, rows$Subtype, FUN = bh)
  rows
}

pairwise_table <- function(scores, metadata, dataset_name, adjust_study) {
  metadata <- metadata[match(rownames(scores), metadata$Sample), ]
  stopifnot(identical(rownames(scores), metadata$Sample))
  metadata$Subtype <- factor(metadata$FinalSubtype, levels = c("C1", "C2", "C3"))
  if (adjust_study) {
    metadata$StudyFactor <- factor(metadata$Study)
    design <- model.matrix(~ 0 + Subtype + StudyFactor, data = metadata)
  } else {
    design <- model.matrix(~ 0 + Subtype, data = metadata)
  }
  colnames(design) <- make.names(colnames(design))
  contrast_matrix <- makeContrasts(
    C1_minus_C2 = SubtypeC1 - SubtypeC2,
    C1_minus_C3 = SubtypeC1 - SubtypeC3,
    C2_minus_C3 = SubtypeC2 - SubtypeC3,
    levels = design
  )
  fit <- eBayes(contrasts.fit(lmFit(t(scores), design), contrast_matrix))
  rows <- do.call(rbind, lapply(seq_len(ncol(contrast_matrix)), function(j) {
    data.frame(
      Dataset = dataset_name,
      Module = rownames(fit$coefficients),
      Contrast = colnames(contrast_matrix)[j],
      N = nrow(scores),
      Studies = length(unique(metadata$Study)),
      EffectSD = fit$coefficients[, j],
      SE = fit$stdev.unscaled[, j] * sqrt(fit$s2.post),
      T = fit$t[, j],
      PValue = fit$p.value[, j],
      stringsAsFactors = FALSE
    )
  }))
  rows$FDR_within_dataset_contrast <- ave(rows$PValue, rows$Dataset, rows$Contrast, FUN = bh)
  rows
}

assoc_reference <- association_table(
  mes_reference,
  master[master$Sample %in% reference_network_samples, ],
  paste0("Reference_", reference_study),
  adjust_study = FALSE
)
assoc_reference_with_flags <- association_table(
  mes_reference_with_flags,
  master[master$Sample %in% reference_samples, ],
  "Reference_with_connectivity_flags_reintroduced",
  adjust_study = FALSE
)
assoc_validation <- association_table(
  mes_validation,
  master[master$Study != reference_study, ],
  "Validation_non_reference_studies",
  adjust_study = TRUE
)
assoc_all <- association_table(
  scores_all,
  master,
  "All_main_study_centered_scores",
  adjust_study = TRUE
)
associations <- rbind(assoc_reference, assoc_validation, assoc_all)

pairwise_reference <- pairwise_table(
  mes_reference,
  master[master$Sample %in% reference_network_samples, ],
  paste0("Reference_", reference_study),
  adjust_study = FALSE
)
pairwise_reference_with_flags <- pairwise_table(
  mes_reference_with_flags,
  master[master$Sample %in% reference_samples, ],
  "Reference_with_connectivity_flags_reintroduced",
  adjust_study = FALSE
)
pairwise_validation <- pairwise_table(
  mes_validation,
  master[master$Study != reference_study, ],
  "Validation_non_reference_studies",
  adjust_study = TRUE
)
pairwise_all <- pairwise_table(
  scores_all,
  master,
  "All_main_study_centered_scores",
  adjust_study = TRUE
)
pairwise_associations <- rbind(pairwise_reference, pairwise_validation, pairwise_all)
write_table(pairwise_associations, "step2_wgcna_pairwise_associations.csv")

flagged_validation_samples <- sample_audit$Sample[
  sample_audit$PrespecifiedOutlierFlag & sample_audit$NetworkRole == "Held-out-study replication"
]
keep_validation <- !rownames(mes_validation) %in% flagged_validation_samples
assoc_validation_no_outlier <- association_table(
  mes_validation[keep_validation, , drop = FALSE],
  master[master$Sample %in% rownames(mes_validation)[keep_validation], ],
  "Validation_without_connectivity_outliers",
  adjust_study = TRUE
)
pairwise_validation_no_outlier <- pairwise_table(
  mes_validation[keep_validation, , drop = FALSE],
  master[master$Sample %in% rownames(mes_validation)[keep_validation], ],
  "Validation_without_connectivity_outliers",
  adjust_study = TRUE
)
outlier_sensitivity <- assoc_validation_no_outlier
write_table(outlier_sensitivity, "step2_wgcna_validation_outlier_sensitivity.csv")
write_table(assoc_reference_with_flags,
            "step2_wgcna_reference_flag_reintroduction_association_sensitivity.csv")
write_table(pairwise_reference_with_flags,
            "step2_wgcna_reference_flag_reintroduction_pairwise_sensitivity.csv")

module_sizes <- as.data.frame(table(module_colors), stringsAsFactors = FALSE)
names(module_sizes) <- c("ModuleColor", "ModuleSize")
module_sizes$Module <- paste0("ME", module_sizes$ModuleColor)
module_sizes <- module_sizes[module_sizes$ModuleColor != "grey", ]
associations <- merge(associations, module_sizes[, c("Module", "ModuleColor", "ModuleSize")],
                      by = "Module", all.x = TRUE, sort = FALSE)
write_table(associations, "step2_wgcna_module_subtype_associations.csv")

preservation_table_path <- file.path(table_dir, "step2_wgcna_module_preservation.csv")
checkpoint_path <- file.path(model_dir, "step2_wgcna_preservation_checkpoint.rds")
checkpoint_parameters <- list(
  seed = 123L, reference_study = reference_study, max_genes = max_genes,
  mad_quantile = mad_quantile, power_candidates = power_candidates,
  scale_free_target = scale_free_target, scale_free_minimum = scale_free_minimum,
  soft_power = soft_power, min_module_size = min_module_size,
  merge_cut_height = merge_cut_height, deep_split = deep_split,
  max_block_size = max_block_size, threads = threads,
  reference_outlier_z_threshold = reference_outlier_z_threshold,
  reference_qc_policy = "exclude flagged reference samples from module discovery",
  preservation_permutations = preservation_permutations,
  network_type = "signed", TOM_type = "signed", correlation = "bicor",
  maxPOutliers = 0.05,
  blockwise_reassign_threshold = 0, blockwise_pam_respects_dendro = FALSE,
  blockwise_numeric_labels = TRUE, blockwise_save_TOMs = FALSE,
  preservation_data_is_expression = TRUE,
  preservation_reference_networks = 1L, preservation_test_networks = 2L,
  preservation_include_kME_all_in_summary = FALSE,
  preservation_restrict_summary_for_general_networks = TRUE,
  preservation_calculate_qvalue = FALSE, preservation_random_seed = 123L,
  preservation_max_gold_module_size = preservation_max_gold_module_size,
  preservation_max_module_size = preservation_max_module_size,
  preservation_quick_cor = 1, preservation_cc_tuplet_size = 2,
  preservation_calculate_cor_kIM_all = FALSE,
  preservation_calculate_cluster_coefficient = FALSE,
  preservation_use_interpolation = FALSE, preservation_check_data = TRUE,
  preservation_grey_name = NULL, preservation_gold_name = NULL,
  preservation_save_permuted_statistics = FALSE,
  preservation_load_permuted_statistics = FALSE,
  preservation_discard_invalid_output = TRUE,
  preservation_parallel = FALSE
)
checkpoint_fingerprint <- list(
  RequiredInputMD5 = setNames(unname(tools::md5sum(required_inputs)),
                              sub(paste0("^", root, "/"), "", required_inputs)),
  RVersion = as.character(getRversion()),
  WGCNAVersion = as.character(packageVersion("WGCNA")),
  Parameters = checkpoint_parameters,
  ReferenceSamplesHash = object_md5(reference_network_samples),
  ValidationSamplesHash = object_md5(validation_samples),
  CandidateGenesHash = object_md5(candidate_genes),
  CompleteGeneModuleColorHash = object_md5(module_colors),
  FixedReferenceLoadingsHash = object_md5(projection$loadings)
)
checkpoint_fingerprint_hash <- object_md5(checkpoint_fingerprint)
resume_requested <- identical(Sys.getenv("SCLUCID_WGCNA_RESUME_PRESERVATION", unset = "0"), "1")
adopt_completed_requested <- identical(Sys.getenv("SCLUCID_WGCNA_ADOPT_COMPLETED_PRESERVATION", unset = "0"), "1")
preservation_table <- NULL

if (!test_mode && resume_requested && file.exists(checkpoint_path)) {
  checkpoint <- readRDS(checkpoint_path)
  checkpoint_valid <- isTRUE(checkpoint$Completed) &&
    identical(checkpoint$Permutations, preservation_permutations_target) &&
    identical(checkpoint$FingerprintHash, checkpoint_fingerprint_hash) &&
    identical(checkpoint$PreservationTableObjectHash, object_md5(checkpoint$PreservationTable)) &&
    file.exists(preservation_table_path) &&
    identical(checkpoint$PreservationTableFileMD5, md5(preservation_table_path))
  if (checkpoint_valid) {
    preservation_table <- checkpoint$PreservationTable
    message("Reusing strictly validated 200-permutation preservation checkpoint")
  } else {
    message("Preservation checkpoint mismatch; recomputing all permutations")
  }
}

# Recovery is restricted to an explicitly requested run after the operator has
# observed permutation 200/200 in the same session. The current in-memory
# fingerprint and complete module mapping must match the recovered table.
if (!test_mode && is.null(preservation_table) && adopt_completed_requested &&
    file.exists(preservation_table_path)) {
  recovered <- read.csv(preservation_table_path, check.names = FALSE)
  recovered <- recovered[order(recovered$Module), ]
  expected <- module_sizes[order(module_sizes$Module), c("Module", "ModuleSize")]
  recovered_valid <- identical(as.character(recovered$Module), as.character(expected$Module)) &&
    identical(as.integer(recovered$ModuleSize), as.integer(expected$ModuleSize)) &&
    "PreservationEvaluationSize" %in% names(recovered) &&
    identical(as.integer(recovered$PreservationEvaluationSize),
              pmin(as.integer(expected$ModuleSize), preservation_max_module_size)) &&
    all(is.finite(recovered$ZsummaryPreservation))
  if (!recovered_valid) stop("Recovered preservation table does not match current module mapping")
  preservation_table <- recovered
  checkpoint <- list(
    Completed = TRUE, Permutations = preservation_permutations_target,
    RecoverySource = "Observed 200/200 session interrupted after preservation by iCloud figure overwrite",
    Fingerprint = checkpoint_fingerprint, FingerprintHash = checkpoint_fingerprint_hash,
    PreservationTable = preservation_table,
    PreservationTableObjectHash = object_md5(preservation_table),
    PreservationTableFileMD5 = md5(preservation_table_path),
    Created = as.character(Sys.time())
  )
  atomic_save_rds(checkpoint, checkpoint_path, compress = "xz")
  message("Adopted completed 200-permutation table after strict current-state validation")
}

if (is.null(preservation_table)) {
  preservation <- modulePreservation(
    multiData = list(
      Reference = list(data = dat_expr_reference_network),
      Validation = list(data = dat_expr_validation)
    ),
    multiColor = list(Reference = module_colors),
    referenceNetworks = 1,
    testNetworks = 2,
    nPermutations = preservation_permutations,
    networkType = "signed",
    corFnc = "bicor",
    corOptions = "use = 'p', maxPOutliers = 0.05",
    randomSeed = 123,
    maxGoldModuleSize = preservation_max_gold_module_size,
    maxModuleSize = preservation_max_module_size,
    savePermutedStatistics = FALSE,
    parallelCalculation = FALSE,
    verbose = 3
  )
  z_pres <- preservation$preservation$Z[[1]][[2]]
  observed_pres <- preservation$preservation$observed[[1]][[2]]
  preservation_table <- data.frame(
    ModuleColor = rownames(z_pres), Module = paste0("ME", rownames(z_pres)),
    ModuleSize = module_sizes$ModuleSize[match(paste0("ME", rownames(z_pres)), module_sizes$Module)],
    PreservationEvaluationSize = z_pres$moduleSize,
    ZsummaryPreservation = z_pres$Zsummary.pres,
    ZdensityPreservation = z_pres$Zdensity.pres,
    ZconnectivityPreservation = z_pres$Zconnectivity.pres,
    MedianRankPreservation = observed_pres[rownames(z_pres), "medianRank.pres"],
    Interpretation = ifelse(z_pres$Zsummary.pres >= 10, "Strong preservation",
                            ifelse(z_pres$Zsummary.pres >= 2, "Moderate preservation", "Weak/no preservation")),
    stringsAsFactors = FALSE
  )
  preservation_table <- preservation_table[!preservation_table$ModuleColor %in% c("grey", "gold"), ]
  write_table(preservation_table, "step2_wgcna_module_preservation.csv")
  if (!test_mode) {
    checkpoint <- list(
      Completed = TRUE, Permutations = preservation_permutations,
      RecoverySource = "Direct modulePreservation execution",
      Fingerprint = checkpoint_fingerprint, FingerprintHash = checkpoint_fingerprint_hash,
      PreservationTable = preservation_table,
      PreservationTableObjectHash = object_md5(preservation_table),
      PreservationTableFileMD5 = md5(preservation_table_path),
      Created = as.character(Sys.time())
    )
    atomic_save_rds(checkpoint, checkpoint_path, compress = "xz")
  }
}

subtypes <- c("C1", "C2", "C3")

# Per-study effects use the same equal-weighted estimand as the pooled models.
# A study contributes to direction and meta-analysis only when all three groups
# have at least two tumors; this excludes the one-C2 2023Lee comparison.
study_effects <- do.call(rbind, lapply(colnames(scores_all), function(module) {
  do.call(rbind, lapply(unique(master$Study), function(study) {
    samples <- master$Sample[master$Study == study]
    subtype <- as.character(master$FinalSubtype[match(samples, master$Sample)])
    counts <- table(factor(subtype, levels = subtypes))
    do.call(rbind, lapply(subtypes, function(target) {
      others <- setdiff(subtypes, target)
      values <- lapply(subtypes, function(x) scores_all[samples[subtype == x], module])
      names(values) <- subtypes
      eligible <- all(counts >= 2L)
      effect <- mean(values[[target]]) - mean(vapply(values[others], mean, numeric(1)))
      se <- if (eligible) sqrt(
        var(values[[target]]) / counts[[target]] +
          sum(vapply(others, function(x) var(values[[x]]) / counts[[x]], numeric(1))) / 4
      ) else NA_real_
      data.frame(
        Module = module, Study = study, TargetSubtype = target,
        N_Target = unname(counts[[target]]), N_Other1 = unname(counts[[others[1]]]),
        N_Other2 = unname(counts[[others[2]]]), EligibleForInference = eligible,
        EffectSD = effect, SE = se,
        CI_low = ifelse(is.finite(se), effect - 1.96 * se, NA_real_),
        CI_high = ifelse(is.finite(se), effect + 1.96 * se, NA_real_),
        Direction = ifelse(effect > 0, paste0(target, " higher"),
                           ifelse(effect < 0, paste0(target, " lower"), "No difference")),
        stringsAsFactors = FALSE
      )
    }))
  }))
}))
write_table(study_effects, "step2_wgcna_subtype_study_effects.csv")

meta_results <- do.call(rbind, lapply(
  split(study_effects, interaction(study_effects$Module, study_effects$TargetSubtype, drop = TRUE)),
  function(d) {
    d_valid <- d[d$EligibleForInference & is.finite(d$SE) & d$SE > 0, ]
    if (nrow(d_valid) < 2L) return(NULL)
    model <- rma.uni(yi = d_valid$EffectSD, sei = d_valid$SE, method = "REML",
                     control = list(stepadj = 0.25, maxiter = 1000))
    data.frame(
      Module = d$Module[1], TargetSubtype = d$TargetSubtype[1],
      DirectionEligibleStudies = nrow(d_valid),
      PositiveEligibleStudies = sum(d_valid$EffectSD > 0),
      DirectionConsistency = mean(d_valid$EffectSD > 0),
      MetaMethod = "REML(stepadj=0.25,maxiter=1000)",
      RandomEffect = as.numeric(model$b), SE = model$se,
      CI_low = model$ci.lb, CI_high = model$ci.ub, PValue = model$pval,
      Tau2 = model$tau2, I2Percent = model$I2, stringsAsFactors = FALSE
    )
  }
))
meta_results$FDR <- ave(meta_results$PValue, meta_results$TargetSubtype, FUN = bh)
write_table(meta_results, "step2_wgcna_subtype_random_effects_meta.csv")

loso_results <- do.call(rbind, lapply(unique(master$Study), function(held_out) {
  keep <- master$Study != held_out
  result <- association_table(
    scores_all[master$Sample[keep], , drop = FALSE], droplevels(master[keep, ]),
    paste0("LOSO_without_", held_out), adjust_study = TRUE
  )
  result$HeldOutStudy <- held_out
  result
}))
write_table(loso_results, "step2_wgcna_subtype_loso_sensitivity.csv")

association_extract <- function(table, subtype, prefix) {
  d <- table[table$Subtype == subtype, c("Module", "EffectSD", "FDR_within_dataset_subtype")]
  names(d)[2:3] <- paste0(prefix, c("Effect", "FDR"))
  d
}

directed_pairwise_extract <- function(table, dataset, target, other, prefix) {
  pair <- paste(sort(c(target, other)), collapse = "_")
  contrast <- switch(pair, C1_C2 = "C1_minus_C2", C1_C3 = "C1_minus_C3",
                     C2_C3 = "C2_minus_C3")
  sign <- ifelse(substr(contrast, 1, 2) == target, 1, -1)
  d <- table[table$Dataset == dataset & table$Contrast == contrast,
             c("Module", "EffectSD", "FDR_within_dataset_contrast")]
  d$EffectSD <- sign * d$EffectSD
  names(d)[2:3] <- paste0(prefix, c("Effect", "FDR"))
  d
}

loso_simple <- do.call(rbind, lapply(
  split(loso_results, interaction(loso_results$Module, loso_results$Subtype, drop = TRUE)),
  function(d) data.frame(
    Module = d$Module[1], Subtype = d$Subtype[1],
    LOSOMinEffect = min(d$EffectSD), LOSOMaxEffect = max(d$EffectSD),
    LOSOPositive = sum(d$EffectSD > 0), LOSOMaxFDR = max(d$FDR_within_dataset_subtype),
    stringsAsFactors = FALSE
  )
))

build_subtype_evidence <- function(target) {
  others <- setdiff(subtypes, target)
  pieces <- list(
    association_extract(assoc_reference, target, "Reference"),
    association_extract(assoc_reference_with_flags, target, "ReferenceWithFlags"),
    association_extract(assoc_validation, target, "Validation"),
    association_extract(assoc_all, target, "AllMain"),
    association_extract(assoc_validation_no_outlier, target, "ValidationNoOutlier")
  )
  datasets <- c(Reference = paste0("Reference_", reference_study),
                ReferenceWithFlags = "Reference_with_connectivity_flags_reintroduced",
                Validation = "Validation_non_reference_studies",
                AllMain = "All_main_study_centered_scores",
                ValidationNoOutlier = "Validation_without_connectivity_outliers")
  pair_tables <- list(
    Reference = pairwise_associations,
    ReferenceWithFlags = pairwise_reference_with_flags,
    Validation = pairwise_associations,
    AllMain = pairwise_associations,
    ValidationNoOutlier = pairwise_validation_no_outlier
  )
  for (prefix in names(datasets)) {
    source <- pair_tables[[prefix]]
    for (other in others) {
      pieces[[length(pieces) + 1L]] <- directed_pairwise_extract(
        source, datasets[[prefix]], target, other, paste0(prefix, "_vs_", other, "_")
      )
    }
  }
  result <- Reduce(function(x, y) merge(x, y, by = "Module", all = TRUE), pieces)
  result$Subtype <- target
  result <- merge(result,
                  preservation_table[, c("Module", "ModuleColor", "ModuleSize",
                                         "ZsummaryPreservation", "Interpretation")],
                  by = "Module", all.x = TRUE)
  result <- merge(result, projection$quality, by = "Module", all.x = TRUE)
  result <- merge(result,
                  reference_outlier_module_mapping[, c("Module", "SensitivityModule",
                                                        "PrimaryCoverage", "Jaccard",
                                                        "HypergeometricFDR")],
                  by = "Module", all.x = TRUE)
  result <- merge(result, meta_results[meta_results$TargetSubtype == target,
                                      setdiff(names(meta_results), "TargetSubtype")],
                  by = "Module", all.x = TRUE)
  result <- merge(result, loso_simple[loso_simple$Subtype == target,
                                    setdiff(names(loso_simple), "Subtype")],
                  by = "Module", all.x = TRUE)
  result$ReplicatedAssociation <- with(result,
    ReferenceEffect >= effect_threshold & ReferenceFDR < 0.05 &
      ValidationEffect >= effect_threshold & ValidationFDR < 0.05 &
      AllMainEffect >= effect_threshold & AllMainFDR < 0.05)
  pairwise_pass <- rep(TRUE, nrow(result))
  reference_flag_reintroduction_pairwise_pass <- rep(TRUE, nrow(result))
  outlier_pairwise_pass <- rep(TRUE, nrow(result))
  for (other in others) {
    for (prefix in c("Reference", "Validation", "AllMain")) {
      pairwise_pass <- pairwise_pass &
        result[[paste0(prefix, "_vs_", other, "_Effect")]] > 0 &
        result[[paste0(prefix, "_vs_", other, "_FDR")]] < 0.05
    }
    reference_flag_reintroduction_pairwise_pass <- reference_flag_reintroduction_pairwise_pass &
      result[[paste0("ReferenceWithFlags_vs_", other, "_Effect")]] > 0 &
      result[[paste0("ReferenceWithFlags_vs_", other, "_FDR")]] < 0.05
    outlier_pairwise_pass <- outlier_pairwise_pass &
      result[[paste0("ValidationNoOutlier_vs_", other, "_Effect")]] > 0 &
      result[[paste0("ValidationNoOutlier_vs_", other, "_FDR")]] < 0.05
  }
  result$PairwiseSpecificity <- pairwise_pass
  result$ReferenceFlagReintroductionAveragePass <- with(result,
    ReferenceWithFlagsEffect >= effect_threshold & ReferenceWithFlagsFDR < 0.05)
  result$ReferenceFlagReintroductionPairwiseSpecificity <-
    reference_flag_reintroduction_pairwise_pass
  result$ReferenceFlagReintroductionSensitivityPass <-
    result$ReferenceFlagReintroductionAveragePass &
    result$ReferenceFlagReintroductionPairwiseSpecificity
  result$OutlierAveragePass <- with(result,
    ValidationNoOutlierEffect >= effect_threshold & ValidationNoOutlierFDR < 0.05)
  result$OutlierPairwiseSpecificity <- outlier_pairwise_pass
  result$OutlierSensitivityPass <- result$OutlierAveragePass & result$OutlierPairwiseSpecificity
  result$StructurePreserved <- result$ZsummaryPreservation >= 2
  result$UnsupervisedScoreEstimable <- with(result,
    is.finite(ReferencePC1VarianceExplained) & ReferencePC1VarianceExplained > 0)
  result$MetaSupported <- with(result,
    DirectionEligibleStudies >= 5 & RandomEffect > 0 & FDR < 0.05)
  result$CrossStudyStable <- with(result,
    DirectionEligibleStudies >= 5 & DirectionConsistency >= 5 / 6 & LOSOPositive == 7)
  result$BaseReplicatedProgram <- with(result,
    ReplicatedAssociation & ReferenceFlagReintroductionAveragePass & OutlierAveragePass &
      PrimaryCoverage >= reference_specific_module_coverage_min & StructurePreserved &
      UnsupervisedScoreEstimable & MetaSupported & CrossStudyStable)
  result$SupportiveSubtypeSpecific <- result$BaseReplicatedProgram &
    result$PairwiseSpecificity & result$ReferenceFlagReintroductionPairwiseSpecificity &
    result$OutlierPairwiseSpecificity
  result
}

rbind_fill <- function(tables) {
  all_names <- unique(unlist(lapply(tables, names)))
  tables <- lapply(tables, function(d) {
    missing <- setdiff(all_names, names(d))
    for (name in missing) d[[name]] <- NA
    d[, all_names, drop = FALSE]
  })
  do.call(rbind, tables)
}
subtype_evidence <- rbind_fill(lapply(subtypes, build_subtype_evidence))
subtype_evidence <- subtype_evidence[order(subtype_evidence$Subtype,
                                           -subtype_evidence$ValidationEffect), ]
write_table(subtype_evidence, "step2_wgcna_subtype_evidence_summary.csv")

pair_value <- function(module, dataset, target, other, field = c("Effect", "FDR"),
                       table = pairwise_associations) {
  field <- match.arg(field)
  d <- directed_pairwise_extract(table, dataset, target, other, "Value")
  value <- d[d$Module == module, paste0("Value", field)]
  ifelse(length(value), value[1], NA_real_)
}

shared_pair_pass <- function(module, a, b, third) {
  base_a <- subtype_evidence[subtype_evidence$Module == module & subtype_evidence$Subtype == a, ]
  base_b <- subtype_evidence[subtype_evidence$Module == module & subtype_evidence$Subtype == b, ]
  if (!nrow(base_a) || !nrow(base_b) || !base_a$BaseReplicatedProgram || !base_b$BaseReplicatedProgram) return(FALSE)
  for (dataset in c(paste0("Reference_", reference_study),
                    "Validation_non_reference_studies", "All_main_study_centered_scores")) {
    if (!(pair_value(module, dataset, a, third, "Effect") > 0 &&
          pair_value(module, dataset, a, third, "FDR") < 0.05 &&
          pair_value(module, dataset, b, third, "Effect") > 0 &&
          pair_value(module, dataset, b, third, "FDR") < 0.05 &&
           abs(pair_value(module, dataset, a, b, "Effect")) <= shared_margin)) return(FALSE)
  }
  reference_flag_dataset <- "Reference_with_connectivity_flags_reintroduced"
  if (!(pair_value(module, reference_flag_dataset, a, third, "Effect",
                   pairwise_reference_with_flags) > 0 &&
        pair_value(module, reference_flag_dataset, a, third, "FDR",
                   pairwise_reference_with_flags) < 0.05 &&
        pair_value(module, reference_flag_dataset, b, third, "Effect",
                   pairwise_reference_with_flags) > 0 &&
        pair_value(module, reference_flag_dataset, b, third, "FDR",
                   pairwise_reference_with_flags) < 0.05 &&
        abs(pair_value(module, reference_flag_dataset, a, b, "Effect",
                       pairwise_reference_with_flags)) <= shared_margin)) return(FALSE)
  outlier_table <- pairwise_validation_no_outlier
  dataset <- "Validation_without_connectivity_outliers"
  pair_value(module, dataset, a, third, "Effect", outlier_table) > 0 &&
    pair_value(module, dataset, a, third, "FDR", outlier_table) < 0.05 &&
    pair_value(module, dataset, b, third, "Effect", outlier_table) > 0 &&
    pair_value(module, dataset, b, third, "FDR", outlier_table) < 0.05 &&
    abs(pair_value(module, dataset, a, b, "Effect", outlier_table)) <= shared_margin
}

module_classification <- do.call(rbind, lapply(module_sizes$Module, function(module) {
  d <- subtype_evidence[subtype_evidence$Module == module, ]
  specific <- d$Subtype[d$SupportiveSubtypeSpecific %in% TRUE]
  shared <- character()
  if (shared_pair_pass(module, "C1", "C2", "C3")) shared <- c(shared, "C1_C2")
  if (shared_pair_pass(module, "C1", "C3", "C2")) shared <- c(shared, "C1_C3")
  if (shared_pair_pass(module, "C2", "C3", "C1")) shared <- c(shared, "C2_C3")
  preserved <- unique(d$StructurePreserved) %in% TRUE
  replicated_any <- any(d$BaseReplicatedProgram %in% TRUE)
  association_any <- any(d$ReferenceFDR < 0.05 & d$ValidationFDR < 0.05 &
                           sign(d$ReferenceEffect) == sign(d$ValidationEffect), na.rm = TRUE)
  if (length(specific) == 1L && !length(shared)) {
    class <- paste0(specific, "_specific")
    assigned <- specific
  } else if (!length(specific) && length(shared) == 1L) {
    class <- paste0("SharedHighCompatible_", shared)
    assigned <- shared
  } else if (!preserved) {
    class <- "Unstable_or_unpreserved"
    assigned <- NA_character_
  } else if (replicated_any || association_any || length(specific) > 1L || length(shared) > 1L) {
    class <- "Replicated_gradient_or_ambiguous"
    assigned <- paste(c(specific, shared), collapse = ";")
  } else {
    class <- "Non_specific"
    assigned <- NA_character_
  }
  data.frame(Module = module, SpecificityClass = class, AssignedSubtypeOrPair = assigned,
             SpecificCandidates = paste(specific, collapse = ";"),
             BivalentCandidates = paste(shared, collapse = ";"), stringsAsFactors = FALSE)
}))
module_classification <- merge(module_sizes, module_classification, by = "Module", all.x = TRUE)
module_classification <- merge(
  module_classification,
  reference_outlier_module_mapping[, c("Module", "SensitivityModule", "PrimaryCoverage",
                                        "Jaccard", "HypergeometricFDR")],
  by = "Module", all.x = TRUE
)
write_table(module_classification, "step2_wgcna_module_subtype_classification.csv")

deg <- read.csv(deg_path, check.names = FALSE)
meta_deg <- read.csv(meta_deg_path, check.names = FALSE)
protein <- read.csv(protein_path, check.names = FALSE)

membership_rows <- do.call(rbind, lapply(unique(module_colors[module_colors != "grey"]), function(color) {
  module <- paste0("ME", color)
  genes <- names(module_colors)[module_colors == color]
  kme <- as.numeric(bicor(dat_expr_reference_network[, genes, drop = FALSE],
                          mes_reference[, module, drop = FALSE],
                          use = "pairwise.complete.obs", maxPOutliers = 0.05))
  data.frame(Gene = genes, Module = module, ModuleColor = color, kME = kme,
             stringsAsFactors = FALSE)
}))
membership <- merge(membership_rows, gene_qc[, c("Gene", "SD_reference", "MAD_reference")],
                    by = "Gene", all.x = TRUE)
for (subtype in subtypes) {
  d <- deg[deg$Subtype == subtype, c("Gene", "logFC_subtype_minus_others", "FDR")]
  names(d)[2:3] <- paste0(subtype, c("_RNA_LogFC", "_RNA_FDR"))
  membership <- merge(membership, d, by = "Gene", all.x = TRUE)
  d <- meta_deg[meta_deg$Subtype == subtype,
                c("Gene", "MetaLogFC", "MetaFDR", "DirectionConsistency", "I2Percent")]
  names(d)[2:5] <- paste0(subtype, c("_MetaLogFC", "_MetaFDR",
                                    "_MetaDirectionConsistency", "_MetaI2Percent"))
  membership <- merge(membership, d, by = "Gene", all.x = TRUE)
  d <- protein[protein$Subtype == subtype, c("Protein", "StandardizedDifference", "FDR")]
  names(d) <- c("Gene", paste0(subtype, "_ProteinEffect"), paste0(subtype, "_ProteinFDR"))
  membership <- merge(membership, d, by = "Gene", all.x = TRUE)
}
membership <- membership[order(membership$Module, -membership$kME), ]
write_table(membership, "step2_wgcna_gene_module_membership.csv")

eigengene_cor <- bicor(mes_reference, use = "pairwise.complete.obs", maxPOutliers = 0.05)
diag(eigengene_cor) <- NA_real_
quality_kme <- do.call(rbind, lapply(split(membership, membership$Module), function(d) {
  data.frame(Module = d$Module[1], MedianKME = median(d$kME),
             Q25KME = unname(quantile(d$kME, 0.25)),
             CoreGeneFractionKME05 = mean(d$kME >= 0.5), MaxKME = max(d$kME),
             stringsAsFactors = FALSE)
}))
module_quality <- Reduce(function(x, y) merge(x, y, by = "Module", all = TRUE), list(
  module_sizes, projection$quality, quality_kme,
  data.frame(Module = rownames(eigengene_cor),
             MaxAbsCorrelationWithOtherModule = apply(abs(eigengene_cor), 1, max, na.rm = TRUE)),
  preservation_table[, c("Module", "ZsummaryPreservation", "MedianRankPreservation", "Interpretation")]
))
module_quality$UnsupervisedMetricsComplete <- with(module_quality,
  is.finite(ReferencePC1VarianceExplained) & is.finite(MedianKME) &
    is.finite(Q25KME) & is.finite(CoreGeneFractionKME05) & is.finite(MaxKME) &
    is.finite(MaxAbsCorrelationWithOtherModule) & is.finite(ZsummaryPreservation) &
    is.finite(MedianRankPreservation))
write_table(module_quality, "step2_wgcna_unsupervised_module_quality.csv")

hub_genes <- do.call(rbind, lapply(module_sizes$Module, function(module) {
  d <- membership[membership$Module == module, ]
  d <- d[order(-d$kME), ]
  d <- head(d, 25)
  merge(d, module_classification[, c("Module", "SpecificityClass", "AssignedSubtypeOrPair")],
        by = "Module", all.x = TRUE)
}))
write_table(hub_genes, "step2_wgcna_subtype_hub_genes.csv")

signature_lists <- readRDS(signature_path)
ductal_11 <- c("HNF1B", "ONECUT1", "FGFR2", "FGFR3", "SLC34A2", "CFTR",
               "AQP1", "CLDN10", "DCDC2", "SPP1", "VTCN1")
historical_magenta <- scan(file.path(root, "1-Clustering", "pattern", "ME4_genes.txt"),
                           what = character(), quiet = TRUE)
historical_pink <- scan(file.path(root, "1-Clustering", "pattern", "ME5_genes.txt"),
                        what = character(), quiet = TRUE)
overlap_sets <- list(
  Prespecified_ductal_11 = ductal_11,
  Operational_classifier_123 = signature_lists$operational_123,
  Historical_magenta_set = historical_magenta,
  Historical_pink_set = historical_pink
)
for (subtype in subtypes) {
  d <- deg[deg$Subtype == subtype, ]
  overlap_sets[[paste0(subtype, "_RNA_up_FDR05")]] <-
    d$Gene[is.finite(d$FDR) & d$FDR < 0.05 & d$logFC_subtype_minus_others > 0]
  overlap_sets[[paste0(subtype, "_RNA_down_FDR05")]] <-
    d$Gene[is.finite(d$FDR) & d$FDR < 0.05 & d$logFC_subtype_minus_others < 0]
  p <- protein[protein$Subtype == subtype, ]
  overlap_sets[[paste0(subtype, "_protein_up_FDR05")]] <-
    p$Protein[is.finite(p$FDR) & p$FDR < 0.05 & p$StandardizedDifference > 0]
}

overlap_test <- function(module_genes, set_genes, universe) {
  module_genes <- intersect(module_genes, universe)
  set_genes <- intersect(set_genes, universe)
  tab <- table(factor(universe %in% module_genes, levels = c(FALSE, TRUE)),
               factor(universe %in% set_genes, levels = c(FALSE, TRUE)))
  result <- fisher.test(tab)
  data.frame(
    UniverseSize = length(universe),
    ModuleGenes = length(module_genes),
    SetGenes = length(set_genes),
    Overlap = length(intersect(module_genes, set_genes)),
    ModuleFraction = ifelse(length(module_genes), length(intersect(module_genes, set_genes)) / length(module_genes), NA_real_),
    OddsRatio = unname(result$estimate),
    PValue = result$p.value,
    OverlapGenes = paste(sort(intersect(module_genes, set_genes)), collapse = ";"),
    stringsAsFactors = FALSE
  )
}

overlap_rows <- do.call(rbind, lapply(module_sizes$Module, function(module) {
  module_genes <- membership$Gene[membership$Module == module]
  do.call(rbind, lapply(names(overlap_sets), function(set_name) {
    universe <- candidate_genes
    if (grepl("_protein_up_FDR05$", set_name)) {
      universe <- intersect(candidate_genes, unique(protein$Protein))
    }
    result <- overlap_test(module_genes, overlap_sets[[set_name]], universe)
    result$Module <- module
    result$GeneSet <- set_name
    result
  }))
}))
if (is.null(overlap_rows)) {
  overlap_rows <- data.frame(Module = character(), GeneSet = character(), UniverseSize = integer(),
                             ModuleGenes = integer(), SetGenes = integer(), Overlap = integer(),
                             ModuleFraction = numeric(), OddsRatio = numeric(), PValue = numeric(),
                             OverlapGenes = character())
} else {
  overlap_rows$FDR_within_gene_set <- ave(overlap_rows$PValue, overlap_rows$GeneSet, FUN = bh)
}
write_table(overlap_rows, "step2_wgcna_module_gene_set_overlap.csv")

symbol_map <- bitr(candidate_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
universe_entrez <- unique(symbol_map$ENTREZID)
hallmark <- unique(msigdbr(species = "Homo sapiens", collection = "H")[, c("gs_name", "gene_symbol")])
reactome <- unique(msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")[,
                            c("gs_name", "gene_symbol")])

enrichment_rows <- list()
for (module in module_sizes$Module) {
  genes <- membership$Gene[membership$Module == module]
  entrez <- unique(symbol_map$ENTREZID[symbol_map$SYMBOL %in% genes])
  go <- as.data.frame(enrichGO(
    gene = entrez, universe = universe_entrez, OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
    pvalueCutoff = 0.05, qvalueCutoff = 0.20, minGSSize = 10, maxGSSize = 500,
    readable = TRUE
  ))
  if (nrow(go)) {
    enrichment_rows[[length(enrichment_rows) + 1L]] <- data.frame(
      Module = module, Collection = "GO_BP", Term = go$Description,
      GeneRatio = go$GeneRatio, Count = go$Count, PValue = go$pvalue,
      FDR = go$p.adjust, Genes = go$geneID, stringsAsFactors = FALSE
    )
  }
  for (collection in c("HALLMARK", "REACTOME")) {
    term2gene <- if (collection == "HALLMARK") hallmark else reactome
    enriched <- as.data.frame(enricher(
      gene = genes, universe = candidate_genes, TERM2GENE = term2gene,
      pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.20,
      minGSSize = 10, maxGSSize = 500
    ))
    if (nrow(enriched)) {
      enrichment_rows[[length(enrichment_rows) + 1L]] <- data.frame(
        Module = module, Collection = collection, Term = enriched$Description,
        GeneRatio = enriched$GeneRatio, Count = enriched$Count, PValue = enriched$pvalue,
        FDR = enriched$p.adjust, Genes = enriched$geneID, stringsAsFactors = FALSE
      )
    }
  }
}
enrichment <- if (length(enrichment_rows)) do.call(rbind, enrichment_rows) else data.frame(
  Module = character(), Collection = character(), Term = character(), GeneRatio = character(),
  Count = integer(), PValue = numeric(), FDR = numeric(), Genes = character()
)
enrichment$FDR_global_within_collection <- ave(enrichment$PValue, enrichment$Collection, FUN = bh)
enrichment <- enrichment[order(enrichment$Module, enrichment$Collection, enrichment$FDR, -enrichment$Count), ]
write_table(enrichment, "step2_wgcna_module_enrichment.csv")

top_annotations <- do.call(rbind, lapply(module_sizes$Module, function(module) {
  d <- enrichment[enrichment$Module == module & enrichment$FDR_global_within_collection < 0.05, ]
  if (!nrow(d)) {
    return(data.frame(Module = module, TopTerms = NA_character_, stringsAsFactors = FALSE))
  }
  top_terms <- do.call(rbind, lapply(split(d, d$Collection), head, 3))
  data.frame(
    Module = module,
    TopTerms = paste(paste0(top_terms$Collection, ": ", top_terms$Term), collapse = " | "),
    stringsAsFactors = FALSE
  )
}))
if (is.null(top_annotations)) top_annotations <- data.frame(Module = character(), TopTerms = character())
module_classification <- merge(module_classification, top_annotations, by = "Module", all.x = TRUE)
module_classification <- merge(module_classification,
                               module_quality[, c("Module", "ReferencePC1VarianceExplained", "MedianKME",
                                                  "CoreGeneFractionKME05", "ZsummaryPreservation",
                                                  "UnsupervisedMetricsComplete")],
                               by = "Module", all.x = TRUE)
write_table(module_classification, "step2_wgcna_module_subtype_classification.csv")

message("Writing symmetric and compatibility evidence tables")
subtype_evidence <- merge(subtype_evidence, top_annotations, by = "Module", all.x = TRUE)
write_table(subtype_evidence, "step2_wgcna_subtype_evidence_summary.csv")

# Backward-compatible C2 views are derived from the symmetric canonical tables.
c2_evidence <- subtype_evidence[subtype_evidence$Subtype == "C2", ]
c2_evidence$C2HighestInAllAnalyses <- c2_evidence$PairwiseSpecificity
c2_evidence$C2SupportiveModule <- c2_evidence$SupportiveSubtypeSpecific
write_table(c2_evidence, "step2_wgcna_c2_evidence_summary.csv")
c2_pairwise <- rbind(
  transform(pairwise_associations[pairwise_associations$Contrast == "C1_minus_C2", ],
            Contrast = "C2_minus_C1", EffectSD = -EffectSD, T = -T),
  pairwise_associations[pairwise_associations$Contrast == "C2_minus_C3", ]
)
write_table(c2_pairwise, "step2_wgcna_c2_pairwise_associations.csv")
c2_study <- study_effects[study_effects$TargetSubtype == "C2", ]
c2_study$C2_minus_others <- c2_study$EffectSD
write_table(c2_study, "step2_wgcna_c2_study_effects.csv")
write_table(meta_results[meta_results$TargetSubtype == "C2", ],
            "step2_wgcna_c2_random_effects_meta.csv")
write_table(loso_results[loso_results$Subtype == "C2", ],
            "step2_wgcna_c2_loso_sensitivity.csv")
c2_modules <- module_classification$Module[module_classification$SpecificityClass == "C2_specific"]
write_table(hub_genes[hub_genes$Module %in% c2_modules, ], "step2_wgcna_c2_hub_genes.csv")
write_table(overlap_rows[overlap_rows$Module %in% c2_modules, ],
            "step2_wgcna_c2_module_gene_set_overlap.csv")
write_table(enrichment[enrichment$Module %in% c2_modules, ],
            "step2_wgcna_c2_module_enrichment.csv")

message("Rendering WGCNA figures")
association_heatmap <- reshape(
  assoc_all[, c("Module", "Subtype", "EffectSD")],
  idvar = "Module", timevar = "Subtype", direction = "wide"
)
rownames(association_heatmap) <- association_heatmap$Module
association_heatmap$Module <- NULL
colnames(association_heatmap) <- sub("EffectSD\\.", "", colnames(association_heatmap))
association_heatmap <- as.matrix(association_heatmap[, c("C1", "C2", "C3"), drop = FALSE])
classification_order <- match(rownames(association_heatmap), module_classification$Module)
association_heatmap <- association_heatmap[
  order(module_classification$SpecificityClass[classification_order],
        -apply(association_heatmap, 1, max)), , drop = FALSE]
heat_colors <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(101)
heatmap_pdf <- file.path(figure_dir, "step2_wgcna_module_subtype_heatmap.pdf")
heatmap_pdf_tmp <- tempfile(pattern = ".step2_wgcna_module_subtype_heatmap.",
                            tmpdir = figure_dir, fileext = ".pdf")
pdf(heatmap_pdf_tmp, width = 6.5,
    height = max(5, 0.25 * nrow(association_heatmap) + 2))
pheatmap(association_heatmap, cluster_rows = FALSE, cluster_cols = FALSE, color = heat_colors,
         border_color = NA, main = "Study-adjusted module associations",
         display_numbers = matrix(sprintf("%.2f", association_heatmap), nrow = nrow(association_heatmap)),
         number_color = "black", fontsize_number = 7)
dev.off()
if (!file.rename(heatmap_pdf_tmp, heatmap_pdf)) stop("Atomic replacement failed: ", heatmap_pdf)
heatmap_png <- file.path(figure_dir, "step2_wgcna_module_subtype_heatmap.png")
heatmap_png_tmp <- tempfile(pattern = ".step2_wgcna_module_subtype_heatmap.",
                            tmpdir = figure_dir, fileext = ".png")
png(heatmap_png_tmp, width = 1800,
    height = max(1500, 80 * nrow(association_heatmap) + 500), res = 300)
pheatmap(association_heatmap, cluster_rows = FALSE, cluster_cols = FALSE, color = heat_colors,
         border_color = NA, main = "Study-adjusted module associations",
         display_numbers = matrix(sprintf("%.2f", association_heatmap), nrow = nrow(association_heatmap)),
         number_color = "black", fontsize_number = 7)
dev.off()
if (!file.rename(heatmap_png_tmp, heatmap_png)) stop("Atomic replacement failed: ", heatmap_png)

palette <- read.csv(palette_path, check.names = FALSE)
subtype_palette <- setNames(palette$Color[palette$SubtypeClass == "Main"],
                            palette$FinalSubtype[palette$SubtypeClass == "Main"])
top_specific <- do.call(rbind, lapply(subtypes, function(subtype) {
  d <- subtype_evidence[subtype_evidence$Subtype == subtype &
                          subtype_evidence$SupportiveSubtypeSpecific %in% TRUE, ]
  head(d[order(-d$ValidationEffect), c("Module", "Subtype")], 2)
}))
if (is.null(top_specific)) top_specific <- data.frame(Module = character(), Subtype = character())
top_shared <- head(module_classification[
  grepl("^SharedHighCompatible_", module_classification$SpecificityClass),
  c("Module", "AssignedSubtypeOrPair")], 2)
top_modules <- unique(c(top_specific$Module, top_shared$Module))
if (length(top_modules)) {
  score_long <- do.call(rbind, lapply(top_modules, function(module) {
    data.frame(Sample = rownames(scores_all), Module = module, Score = scores_all[, module], stringsAsFactors = FALSE)
  }))
  score_long <- merge(score_long, master, by = "Sample", all.x = TRUE, sort = FALSE)
  score_plot <- ggplot(score_long, aes(FinalSubtype, Score, fill = FinalSubtype)) +
    geom_boxplot(outlier.shape = NA, width = 0.68) +
    geom_jitter(width = 0.16, alpha = 0.18, size = 0.7) +
    facet_wrap(~ Module, scales = "free_y", ncol = 2) +
    scale_fill_manual(values = subtype_palette, guide = "none") +
    labs(x = NULL, y = "Within-study standardized module score",
         title = "Replicated subtype-associated co-expression modules") +
    theme_bw(base_size = 11)
  atomic_ggsave(file.path(figure_dir, "step2_wgcna_subtype_module_scores.pdf"), score_plot,
                width = 10, height = max(6.5, ceiling(length(top_modules) / 2) * 3))
  atomic_ggsave(file.path(figure_dir, "step2_wgcna_subtype_module_scores.png"), score_plot,
                width = 10, height = max(6.5, ceiling(length(top_modules) / 2) * 3), dpi = 300)

  forest_keys <- top_specific[!duplicated(top_specific$Subtype), ]
  forest_data <- merge(study_effects, forest_keys,
                       by.x = c("Module", "TargetSubtype"), by.y = c("Module", "Subtype"))
  forest_plot <- ggplot(forest_data, aes(EffectSD, Study)) +
    geom_vline(xintercept = 0, linetype = 2, color = "grey50") +
    geom_errorbar(aes(xmin = CI_low, xmax = CI_high), orientation = "y", width = 0.18, na.rm = TRUE) +
    geom_point(aes(color = TargetSubtype), size = 2.1) +
    facet_wrap(~ paste0(Module, " / ", TargetSubtype), scales = "free_x") +
    scale_color_manual(values = subtype_palette, guide = "none") +
    labs(x = "Target subtype minus equal-weighted mean of other two (SD)", y = NULL,
         title = "Cross-study direction of subtype-module effects") +
    theme_bw(base_size = 10)
  atomic_ggsave(file.path(figure_dir, "step2_wgcna_subtype_cross_study_effects.pdf"),
                forest_plot, width = 11, height = 5.5)
  atomic_ggsave(file.path(figure_dir, "step2_wgcna_subtype_cross_study_effects.png"),
                forest_plot, width = 11, height = 5.5, dpi = 300)

  c2_top <- head(top_specific$Module[top_specific$Subtype == "C2"], 4)
  if (length(c2_top)) {
    c2_score_plot <- score_plot + score_long[score_long$Module %in% c2_top, ] +
      labs(title = "C2-specific co-expression modules (compatibility view)")
    atomic_ggsave(file.path(figure_dir, "step2_wgcna_c2_module_scores.pdf"), c2_score_plot,
                  width = 8.5, height = 6.5)
    atomic_ggsave(file.path(figure_dir, "step2_wgcna_c2_module_scores.png"), c2_score_plot,
                  width = 8.5, height = 6.5, dpi = 300)
    c2_forest <- forest_plot + study_effects[
      study_effects$Module %in% head(c2_top, 3) & study_effects$TargetSubtype == "C2", ] +
      labs(title = "C2 module effects (compatibility view)")
    atomic_ggsave(file.path(figure_dir, "step2_wgcna_c2_cross_study_effects.pdf"), c2_forest,
                  width = 10, height = 5.5)
    atomic_ggsave(file.path(figure_dir, "step2_wgcna_c2_cross_study_effects.png"), c2_forest,
                  width = 10, height = 5.5, dpi = 300)
  }
}

if (nrow(enrichment)) {
  annotated_modules <- module_classification$Module[
    !module_classification$SpecificityClass %in% c("Non_specific", "Unstable_or_unpreserved")]
  annotated_modules <- head(annotated_modules, 12)
  enrichment_plot_data <- enrichment[enrichment$Module %in% annotated_modules &
                                       enrichment$FDR_global_within_collection < 0.05, ]
  enrichment_top <- do.call(rbind, lapply(
    split(enrichment_plot_data, interaction(enrichment_plot_data$Module,
                                            enrichment_plot_data$Collection)), head, 3))
  if (!is.null(enrichment_top) && nrow(enrichment_top)) {
  enrichment_top$GeneRatioNumeric <- vapply(strsplit(enrichment_top$GeneRatio, "/"),
                                             function(x) as.numeric(x[1]) / as.numeric(x[2]), numeric(1))
  enrichment_top$Term <- factor(enrichment_top$Term,
                                levels = rev(unique(enrichment_top$Term[order(enrichment_top$FDR)])))
  enrichment_plot <- ggplot(enrichment_top, aes(GeneRatioNumeric, Term, size = Count, color = -log10(FDR))) +
    geom_point() +
    facet_grid(Module ~ Collection, scales = "free_y", space = "free_y") +
    scale_color_viridis_c(option = "C", end = 0.9) +
    labs(x = "Gene ratio", y = NULL, color = expression(-log[10](FDR)),
         title = "Functional annotation of subtype-associated modules") +
    theme_bw(base_size = 9) +
    theme(axis.text.y = element_text(size = 7))
  atomic_ggsave(file.path(figure_dir, "step2_wgcna_subtype_module_enrichment.pdf"), enrichment_plot,
                width = 12, height = max(6, 1.0 * length(annotated_modules) + 4))
  atomic_ggsave(file.path(figure_dir, "step2_wgcna_subtype_module_enrichment.png"), enrichment_plot,
                width = 12, height = max(6, 1.0 * length(annotated_modules) + 4), dpi = 300)
  c2_enrichment_top <- enrichment_top[enrichment_top$Module %in% c2_modules, ]
  if (nrow(c2_enrichment_top)) {
    c2_enrichment_plot <- enrichment_plot + c2_enrichment_top +
      labs(title = "C2 module annotation (compatibility view)")
    atomic_ggsave(file.path(figure_dir, "step2_wgcna_c2_module_enrichment.pdf"),
                  c2_enrichment_plot, width = 12, height = max(6, length(c2_modules) + 4))
    atomic_ggsave(file.path(figure_dir, "step2_wgcna_c2_module_enrichment.png"),
                  c2_enrichment_plot, width = 12, height = max(6, length(c2_modules) + 4), dpi = 300)
  }
  }
}

dendrogram_pdf <- file.path(figure_dir, "step2_wgcna_module_dendrograms.pdf")
dendrogram_pdf_tmp <- tempfile(pattern = ".step2_wgcna_module_dendrograms.",
                               tmpdir = figure_dir, fileext = ".pdf")
pdf(dendrogram_pdf_tmp, width = 12, height = 6)
for (block in seq_along(net$dendrograms)) {
  plotDendroAndColors(
    net$dendrograms[[block]],
    module_colors[net$blockGenes[[block]]],
    "Module colors",
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05,
    main = paste0("WGCNA module dendrogram: block ", block)
  )
}
dev.off()
if (!file.rename(dendrogram_pdf_tmp, dendrogram_pdf)) stop("Atomic replacement failed: ", dendrogram_pdf)
for (block in seq_along(net$dendrograms)) {
  dendrogram_png <- file.path(figure_dir, paste0("step2_wgcna_module_dendrogram_block", block, ".png"))
  dendrogram_png_tmp <- tempfile(pattern = paste0(".step2_wgcna_module_dendrogram_block", block, "."),
                                 tmpdir = figure_dir, fileext = ".png")
  png(dendrogram_png_tmp,
      width = 3000, height = 1500, res = 300)
  plotDendroAndColors(
    net$dendrograms[[block]], module_colors[net$blockGenes[[block]]], "Module colors",
    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
    main = paste0("WGCNA module dendrogram: block ", block)
  )
  dev.off()
  if (!file.rename(dendrogram_png_tmp, dendrogram_png)) {
    stop("Atomic replacement failed: ", dendrogram_png)
  }
}

module_summary <- merge(module_classification, module_quality, by = c("Module", "ModuleColor", "ModuleSize"),
                        all.x = TRUE, sort = FALSE)
write_table(module_summary, "step2_wgcna_module_summary.csv")

network_model <- list(
  analysis = "Step2 study-aware subtype-symmetric WGCNA",
  reference_study = reference_study,
  reference_source_samples = reference_samples,
  reference_samples = reference_network_samples,
  validation_samples = validation_samples,
  network_genes = candidate_genes,
  module_labels = module_labels,
  module_colors = module_colors,
  reference_eigengene_loadings = projection$loadings,
  module_eigengenes_reference = mes_reference,
  module_eigengenes_validation = mes_validation,
  fixed_loading_projected_study_standardized_scores = scores_all,
  reference_outlier_sensitivity = list(
    flagged_samples = flagged_reference_samples,
    sensitivity_samples = sensitivity_reference_samples,
    sensitivity_module_colors = sensitivity_module_colors,
    network_summary = reference_outlier_network_sensitivity,
    module_mapping = reference_outlier_module_mapping
  ),
  soft_power = soft_power,
  soft_threshold_diagnostics = soft_table,
  block_genes = net$blockGenes,
  dendrograms = net$dendrograms,
  parameters = list(
    seed = 123, networkType = "signed", TOMType = "signed", corType = "bicor",
    maxPOutliers = 0.05, max_genes = max_genes, mad_quantile = mad_quantile,
    minModuleSize = min_module_size, mergeCutHeight = merge_cut_height,
    deepSplit = deep_split, maxBlockSize = max_block_size,
    preservation_permutations = preservation_permutations,
    subtype_effect_threshold = effect_threshold, shared_descriptive_margin = shared_margin
  )
)
atomic_save_rds(network_model, file.path(model_dir, "step2_wgcna_network.rds"), compress = "xz")

parameter_table <- data.frame(
  Parameter = c(
    "analysis_date", "R_version", "WGCNA_version", "seed", "reference_study",
    "reference_source_n", "reference_network_n", "validation_n", "validation_studies", "main_n", "network_genes",
    "gene_filter", "soft_power", "power_selection_rule", "selected_signed_R2",
    "selected_mean_connectivity", "network_type", "TOM_type", "correlation",
    "maxPOutliers", "min_module_size", "merge_cut_height", "deep_split",
    "max_block_size", "threads", "preservation_permutations",
    "preservation_max_module_size", "preservation_max_gold_module_size",
    "reference_outlier_z_threshold", "reference_sensitivity_ari_min",
    "reference_sensitivity_coverage_min", "reference_sensitivity_non_grey_min",
    "reference_specific_module_coverage_min", "reference_outlier_flags",
    "module_score_definition",
    "subtype_estimand", "subtype_effect_threshold", "shared_descriptive_margin",
    "experimental_unit", "claim_boundary"
  ),
  Value = c(
    as.character(Sys.Date()), as.character(getRversion()), as.character(packageVersion("WGCNA")),
    "123", reference_study, length(reference_samples), length(reference_network_samples), length(validation_samples),
    length(unique(master$Study[master$Study != reference_study])), nrow(master), length(candidate_genes),
    paste0("finite non-zero genes; reference MAD >= Q", mad_quantile * 100, "; cap ", max_genes),
    soft_power, power_rule, signif(soft_table$SignedR2[selected_row], 6),
    signif(soft_table$mean.k.[selected_row], 6), "signed", "signed", "bicor",
    "0.05", min_module_size, merge_cut_height, deep_split, max_block_size, threads,
    preservation_permutations, preservation_max_module_size, preservation_max_gold_module_size,
    reference_outlier_z_threshold, reference_sensitivity_ari_min,
    reference_sensitivity_coverage_min, reference_sensitivity_non_grey_min,
    reference_specific_module_coverage_min, length(flagged_reference_samples),
    "reference-study PC1 loadings projected unchanged; scores standardized within Study",
    "target subtype minus equal-weighted mean of the other two",
    effect_threshold, shared_margin, "patient/tumor sample",
    "E1 supportive co-expression association; not lineage, mechanism, treatment-response, or naming validation"
  ),
  stringsAsFactors = FALSE
)
write_table(parameter_table, "step2_wgcna_analysis_contract.csv")

scale_free_status <- if (soft_table$SignedR2[selected_row] >= scale_free_minimum) "PASS" else "WARN"
replicated_n <- sum(c2_evidence$ReplicatedAssociation %in% TRUE, na.rm = TRUE)
pairwise_specific_n <- sum(c2_evidence$C2HighestInAllAnalyses %in% TRUE, na.rm = TRUE)
outlier_sensitivity_n <- sum(c2_evidence$OutlierSensitivityPass %in% TRUE, na.rm = TRUE)
supportive_n <- sum(c2_evidence$C2SupportiveModule %in% TRUE, na.rm = TRUE)
specific_counts <- setNames(vapply(subtypes, function(x) {
  sum(module_classification$SpecificityClass == paste0(x, "_specific"))
}, integer(1)), subtypes)
shared_n <- sum(grepl("^SharedHighCompatible_", module_classification$SpecificityClass))
gradient_n <- sum(module_classification$SpecificityClass == "Replicated_gradient_or_ambiguous")
nonspecific_n <- sum(module_classification$SpecificityClass == "Non_specific")
unstable_n <- sum(module_classification$SpecificityClass == "Unstable_or_unpreserved")
strict_modules <- module_classification$Module[
  grepl("^C[123]_specific$", module_classification$SpecificityClass)
]
strict_mapping_pass <- !length(strict_modules) || all(
  module_classification$PrimaryCoverage[module_classification$Module %in% strict_modules] >=
    reference_specific_module_coverage_min
)
strict_rows <- module_classification[
  grepl("^C[123]_specific$", module_classification$SpecificityClass),
  c("Module", "SpecificityClass")
]
strict_rows$Subtype <- sub("_specific$", "", strict_rows$SpecificityClass)
strict_evidence <- merge(strict_rows, subtype_evidence, by = c("Module", "Subtype"), all.x = TRUE)
strict_reintroduction_pass <- !nrow(strict_evidence) || all(
  strict_evidence$ReferenceFlagReintroductionSensitivityPass %in% TRUE
)
reference_qc_exclusion_complete <- identical(
  sort(reference_network_samples), sort(setdiff(reference_samples, flagged_reference_samples))
)
reference_outlier_gate_status <- if (!length(flagged_reference_samples)) {
  "PASS"
} else if (reference_qc_exclusion_complete && strict_mapping_pass && strict_reintroduction_pass) {
  "QUALIFIED"
} else {
  "FAIL"
}
gate <- data.frame(
  Check = c(
    "main_sample_contract", "reference_discovery_contract", "held_out_study_replication_contract",
    "phenotype_independent_gene_filter", "all_matrix_good_samples_genes",
    "reference_connectivity_outlier_check", "scale_free_fit", "non_grey_modules",
    "fixed_reference_loading_projection", "unsupervised_module_metrics_complete",
    "module_preservation_completed", "symmetric_average_contrasts_complete",
    "symmetric_pairwise_contrasts_complete", "study_meta_eligibility_respected",
    "symmetric_loso_complete", "symmetric_module_classification_complete",
    "subtype_specific_modules_observed", "shared_high_compatible_modules_observed",
    "compatibility_c2_view_derived", "historical_module_mapping_not_reused",
    "claim_boundary_recorded", "subtype_suffix_naming_gate"
  ),
  Status = c(
    ifelse(nrow(master) == expected_main_n, "PASS", "FAIL"),
    ifelse(length(reference_samples) == expected_reference_n &&
             length(reference_network_samples) == expected_reference_n - length(flagged_reference_samples),
           "PASS", "FAIL"),
    ifelse(length(validation_samples) == 274L && length(unique(master$Study[master$Study != reference_study])) == 6L, "PASS", "FAIL"),
    ifelse(length(candidate_genes) >= 5000L, "PASS", "FAIL"),
    ifelse(gsg_all$allOK && gsg_reference_network$allOK, "PASS", "FAIL"),
    reference_outlier_gate_status,
    scale_free_status,
    ifelse(nrow(module_sizes) >= 3L, "PASS", "FAIL"),
    ifelse(all(is.finite(scores_all)) && nrow(projection$loadings) == sum(module_sizes$ModuleSize), "PASS", "FAIL"),
    ifelse(nrow(module_quality) == nrow(module_sizes) &&
             isTRUE(all(module_quality$UnsupervisedMetricsComplete)), "PASS", "FAIL"),
    ifelse(nrow(preservation_table) == nrow(module_sizes) && !test_mode &&
             identical(as.integer(preservation_table$ModuleSize[match(module_sizes$Module, preservation_table$Module)]),
                       as.integer(module_sizes$ModuleSize)) &&
             identical(as.integer(preservation_table$PreservationEvaluationSize[match(module_sizes$Module, preservation_table$Module)]),
                       as.integer(module_sizes$ModuleSize)), "PASS", "FAIL"),
    ifelse(nrow(associations) == nrow(module_sizes) * 3L * 3L, "PASS", "FAIL"),
    ifelse(nrow(pairwise_associations) == nrow(module_sizes) * 3L * 3L, "PASS", "FAIL"),
    ifelse(!any(!study_effects$EligibleForInference & is.finite(study_effects$SE)), "PASS", "FAIL"),
    ifelse(nrow(loso_results) == nrow(module_sizes) * 3L * length(unique(master$Study)), "PASS", "FAIL"),
    ifelse(nrow(module_classification) == nrow(module_sizes) && !anyDuplicated(module_classification$Module), "PASS", "FAIL"),
    ifelse(sum(specific_counts) > 0, "PASS", "NEGATIVE"),
    ifelse(shared_n > 0, "PASS", "NEGATIVE"),
    ifelse(all(c2_evidence$Subtype == "C2"), "PASS", "FAIL"),
    "PASS", "PASS", "NOT_MET"
  ),
  Blocking = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, TRUE, TRUE, TRUE,
               TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, TRUE, TRUE,
               TRUE, FALSE),
  Detail = c(
    paste0("N=", nrow(master), "; C1=", sum(master$FinalSubtype == "C1"),
           ", C2=", sum(master$FinalSubtype == "C2"), ", C3=", sum(master$FinalSubtype == "C3")),
    paste0(reference_study, " source N=", length(reference_samples), "; QC-pass network N=",
           length(reference_network_samples), "; modules discovered without subtype labels"),
    paste0("Held-out non-reference studies N=", length(validation_samples), " across ",
           length(unique(master$Study[master$Study != reference_study])),
           " studies; internal cross-study replication, not independent subtype truth"),
    paste0(length(candidate_genes), " genes selected by reference MAD only"),
    "goodSamplesGenes passed for the all-study matrix and QC-pass reference network",
    paste0(length(flagged_reference_samples), " reference flags (",
           paste(flagged_reference_samples, collapse = ";"),
           ") excluded from phenotype-blind module discovery and reintroduced for association sensitivity; all-sample alternative-network ARI=",
           signif(reference_outlier_network_sensitivity$AdjustedRandIndexPrimaryNonGrey, 4),
           ", weighted best-match coverage=",
           signif(reference_outlier_network_sensitivity$WeightedBestMatchCoverage, 4),
           ", non-grey retention=",
           signif(reference_outlier_network_sensitivity$PrimaryNonGreyRetainedFraction, 4),
           "; all-sample network stability threshold pass=",
           reference_outlier_network_sensitivity$NetworkSensitivityPass,
           "; strict-module mapping threshold=", reference_specific_module_coverage_min,
           "; strict flag-reintroduction pass=", strict_reintroduction_pass,
           "; ", length(flagged_validation_samples),
           " validation flags retained in primary analysis and removed in sensitivity"),
    paste0("Power=", soft_power, "; signed R2=", signif(soft_table$SignedR2[selected_row], 4),
           "; rule=", power_rule),
    paste0(nrow(module_sizes), " non-grey modules"),
    paste0(nrow(projection$loadings), " gene loadings fixed in ", reference_study,
           " and projected unchanged to all studies"),
    paste0(sum(module_quality$UnsupervisedMetricsComplete), "/", nrow(module_quality),
           " modules have complete finite PC1VE, kME, eigengene-separation and preservation diagnostics; continuous values are reported without an arbitrary quality cutoff"),
    paste0(nrow(preservation_table), " modules evaluated with ", preservation_permutations,
           " pooled held-out permutations using every module member (maxModuleSize=",
           preservation_max_module_size, ")"),
    paste0(nrow(associations), " module-by-subtype average contrasts across three analysis sets"),
    paste0(nrow(pairwise_associations), " module-by-pair contrasts across three analysis sets"),
    paste0("Only studies with all C1/C2/C3 n>=2 enter meta-analysis; eligible studies per contrast=",
           paste(sort(unique(meta_results$DirectionEligibleStudies)), collapse = ",")),
    paste0(nrow(loso_results), " module-by-subtype leave-one-study-out estimates"),
    paste0(nrow(module_classification), " mutually exclusive module classifications"),
    paste0("C1=", specific_counts[["C1"]], ", C2=", specific_counts[["C2"]],
           ", C3=", specific_counts[["C3"]], " strict subtype-specific modules"),
    paste0(shared_n, " shared-high-compatible modules; margin=", shared_margin,
           " SD is descriptive and is not an equivalence test"),
    paste0("C2 compatibility tables are exact subsets/transformations of symmetric outputs; C2 strict=",
           supportive_n),
    "Numeric labels and WGCNA colors are stored directly; legacy ME4/ME5 assumptions are not used",
    "Results are E1 supportive associations from the same expression resource",
    "WGCNA alone cannot establish a C1/C2/C3 biological suffix"
  ),
  stringsAsFactors = FALSE
)
write_table(gate, "step2_wgcna_release_gate_status.csv")
blocking_failures <- gate$Check[gate$Blocking & gate$Status == "FAIL"]
if (length(blocking_failures)) {
  message("Blocking WGCNA gate failures: ", paste(blocking_failures, collapse = "; "))
  stop("One or more blocking WGCNA release checks failed")
}

reported_modules <- module_classification[
  !module_classification$SpecificityClass %in% c("Non_specific", "Unstable_or_unpreserved"), ]
top_module_lines <- apply(reported_modules, 1, function(x) {
  paste0("- `", x[["Module"]], "` (", x[["ModuleSize"]], " genes): **",
         x[["SpecificityClass"]], "**; preservation Z=",
         sprintf("%.1f", as.numeric(x[["ZsummaryPreservation"]])), "; ",
         ifelse(is.na(x[["TopTerms"]]), "no globally FDR-significant term",
                x[["TopTerms"]]), ".")
})
if (!length(top_module_lines)) top_module_lines <- "- No module met a replicated subtype-associated class."

report_lines <- c(
  "# Step2 WGCNA results summary",
  "",
  paste0("Generated: ", Sys.Date()),
  "",
  "## Decision context",
  "",
  "This refresh asks which unsupervised co-expression modules show reproducible C1, C2 or C3 specificity, shared-high compatibility, gradients or no subtype specificity. No subtype is used as the prior selection target.",
  "",
  "## Design",
  "",
  paste0("- Experimental unit: patient/tumor sample; 497 main-state tumors (C1=", sum(master$FinalSubtype == "C1"),
         ", C2=", sum(master$FinalSubtype == "C2"), ", C3=", sum(master$FinalSubtype == "C3"), ")."),
  paste0("- Reference discovery source: ", reference_study, " (n=", length(reference_samples),
         "); phenotype-blind QC-pass module network n=", length(reference_network_samples), "."),
  paste0("- Held-out-study replication: six non-reference studies (n=", length(validation_samples), "); this is internal cross-study replication, not independent subtype-truth validation."),
  paste0("- Network: ", length(candidate_genes), " phenotype-independent MAD-filtered genes; signed bicor network; power ",
         soft_power, "; minimum module size ", min_module_size, "; merge cut height ", merge_cut_height, "."),
  paste0("- Structural preservation: ", preservation_permutations,
         " permutations in pooled held-out studies; maxModuleSize=",
         preservation_max_module_size, " so every member of every discovered module is evaluated."),
  "- Preservation boundary: pooled Zsummary depends on the composition, sample sizes and residual technical structure of the six held-out studies; it is structural support, not cohort-invariant proof.",
  "- Score contract: PC1 loadings learned once in 2022Dong, projected unchanged, then standardized within Study.",
  "- Symmetric estimand: each target subtype minus the equal-weighted mean of the other two; validation and all-main models adjust for Study.",
  "",
  "## Release result",
  "",
  paste0("The technical WGCNA gate passed with ", nrow(module_sizes), " non-grey modules: C1-specific=",
         specific_counts[["C1"]], ", C2-specific=", specific_counts[["C2"]],
         ", C3-specific=", specific_counts[["C3"]], ", shared-high-compatible=", shared_n,
         ", replicated-gradient/ambiguous=", gradient_n, ", non-specific=", nonspecific_n,
         ", unstable/unpreserved=", unstable_n, "."),
  paste0(length(flagged_reference_samples), " reference-discovery samples had within-study connectivity Z < ",
         reference_outlier_z_threshold, ". They were excluded from primary module construction and reintroduced for association sensitivity; the all-sample alternative network had adjusted Rand index ",
         signif(reference_outlier_network_sensitivity$AdjustedRandIndexPrimaryNonGrey, 3),
         ", weighted best-match coverage ",
         signif(reference_outlier_network_sensitivity$WeightedBestMatchCoverage, 3),
         " and non-grey retention ",
         signif(reference_outlier_network_sensitivity$PrimaryNonGreyRetainedFraction, 3),
         ". The reference-outlier gate is ", reference_outlier_gate_status, "."),
  paste0(length(flagged_validation_samples), " held-out-study samples had within-study connectivity Z < -3. They were retained in the primary analysis and removed only in sensitivity analyses."),
  "",
  "## Symmetric module classification",
  "",
  top_module_lines,
  "",
  "## Classification definitions",
  "",
  paste0("- **Subtype-specific:** effect >= ", effect_threshold, " SD and FDR<0.05 for target-versus-average plus target above both other subtypes in reference, held-out pooled and all-main analyses; preserved, meta-supported, LOSO-stable, robust to reference/validation connectivity flags, and mapped with primary-module coverage >=", reference_specific_module_coverage_min, " in the reference network rebuilt without flagged samples."),
  paste0("- **Shared-high-compatible:** two subtype scores both exceed the third in all analysis sets while their absolute difference is <=", shared_margin, " SD. The margin is descriptive; this category means pairwise-unresolved/shared enrichment and is not a formal equivalence claim."),
  "- **Replicated gradient or ambiguous:** replicated subtype association exists but strict unique/shared-high rules do not resolve one class.",
  "- **Non-specific:** no replicated subtype association meets the prespecified evidence rule.",
  "- **Unstable or unpreserved:** structural preservation fails even if a pooled association is present.",
  "",
  "## Interpretation boundary",
  "",
  "- **Supported (E1):** strict modules are reproducible subtype-associated co-expression programs within this integrated bulk resource.",
  "- **Exploratory:** hub genes and enriched pathways prioritize biological descriptions and follow-up targets.",
  "- **Not supported:** direct regulation, stable lineage identity, immune sensitivity, therapeutic response, or a formal C2 biological suffix.",
  "- Module membership and subtype labels ultimately originate from the same project expression resource; cross-study replication improves robustness but is not independent subtype-truth validation.",
  "",
  "## Canonical artifacts",
  "",
  "- Script: `2-TranscriptomicCharacterization/scripts/run_step2_wgcna_refresh.R`",
  "- Analysis contract: `publication/tables/step2_wgcna_analysis_contract.csv`",
  "- Module associations: `publication/tables/step2_wgcna_module_subtype_associations.csv`",
  "- Symmetric pairwise associations: `publication/tables/step2_wgcna_pairwise_associations.csv`",
  "- Symmetric evidence and classification: `publication/tables/step2_wgcna_subtype_evidence_summary.csv` and `step2_wgcna_module_subtype_classification.csv`",
  "- Module preservation: `publication/tables/step2_wgcna_module_preservation.csv`",
  "- Reference-outlier network sensitivity: `publication/tables/step2_wgcna_reference_outlier_network_sensitivity.csv` and `step2_wgcna_reference_outlier_module_mapping.csv`",
  "- Strict preservation checkpoint: `publication/models/step2_wgcna_preservation_checkpoint.rds`",
  "- Hub genes and enrichment: `publication/tables/step2_wgcna_subtype_hub_genes.csv` and `step2_wgcna_module_enrichment.csv`",
  "- C2-prefixed artifacts are compatibility views only; symmetric tables are canonical.",
  "- Release gate: `publication/tables/step2_wgcna_release_gate_status.csv`",
  "- Saved network interface: `publication/models/step2_wgcna_network.rds`"
)
report_path <- file.path(report_dir, "step2_wgcna_results_summary.md")
report_tmp <- tempfile(pattern = ".step2_wgcna_results_summary.", tmpdir = report_dir)
writeLines(report_lines, report_tmp)
if (!file.rename(report_tmp, report_path)) stop("Atomic replacement failed: ", report_path)

session_path <- file.path(report_dir, "step2_wgcna_session_info.txt")
session_tmp <- tempfile(pattern = ".step2_wgcna_session_info.", tmpdir = report_dir)
capture.output(sessionInfo(), file = session_tmp)
if (!file.rename(session_tmp, session_path)) stop("Atomic replacement failed: ", session_path)

run_manifest <- data.frame(
  Analysis = "Step2 study-aware subtype-symmetric WGCNA refresh",
  Status = "PASS",
  Date = as.character(Sys.Date()),
  ExpressionInput = sub(paste0("^", root, "/"), "", expr_path),
  ExpressionMD5 = md5(expr_path),
  MasterInput = sub(paste0("^", root, "/"), "", master_path),
  MasterMD5 = md5(master_path),
  ReferenceStudy = reference_study,
  MainSamples = nrow(master),
  ReferenceSourceSamples = length(reference_samples),
  ReferenceNetworkSamples = length(reference_network_samples),
  ReferenceFlaggedSamples = length(flagged_reference_samples),
  ValidationSamples = length(validation_samples),
  NetworkGenes = length(candidate_genes),
  NonGreyModules = nrow(module_sizes),
  C1SpecificModules = specific_counts[["C1"]],
  C2SpecificModules = specific_counts[["C2"]],
  C3SpecificModules = specific_counts[["C3"]],
  SharedHighCompatibleModules = shared_n,
  GradientAmbiguousModules = gradient_n,
  NonSpecificModules = nonspecific_n,
  UnstableModules = unstable_n,
  ReplicatedC2Modules = replicated_n,
  PairwiseSpecificC2Modules = pairwise_specific_n,
  SupportiveC2Modules = supportive_n,
  ClaimLevel = "E1 supportive association",
  stringsAsFactors = FALSE
)
write_table(run_manifest, "step2_wgcna_run_manifest.csv")

manifest_candidates <- c(
  list.files(table_dir, pattern = "^step2_wgcna_", full.names = TRUE),
  list.files(figure_dir, pattern = "^step2_wgcna_", full.names = TRUE),
  file.path(model_dir, c("step2_wgcna_network.rds", "step2_wgcna_preservation_checkpoint.rds")),
  file.path(report_dir, c("step2_wgcna_results_summary.md", "step2_wgcna_session_info.txt"))
)
manifest_candidates <- manifest_candidates[file.exists(manifest_candidates)]
manifest_candidates <- manifest_candidates[basename(manifest_candidates) != "step2_wgcna_output_manifest.csv"]
output_manifest <- data.frame(
  Path = sub(paste0("^", root, "/"), "", manifest_candidates),
  Bytes = file.info(manifest_candidates)$size,
  MD5 = unname(tools::md5sum(manifest_candidates)),
  NonEmpty = file.info(manifest_candidates)$size > 0,
  stringsAsFactors = FALSE
)
write_table(output_manifest, "step2_wgcna_output_manifest.csv")

cat(paste0(
  "{\"status\":\"PASS\",\"main_samples\":", nrow(master),
  ",\"reference_source_samples\":", length(reference_samples),
  ",\"reference_network_samples\":", length(reference_network_samples),
  ",\"validation_samples\":", length(validation_samples),
  ",\"network_genes\":", length(candidate_genes),
  ",\"modules\":", nrow(module_sizes),
  ",\"c1_specific_modules\":", specific_counts[["C1"]],
  ",\"c2_specific_modules\":", specific_counts[["C2"]],
  ",\"c3_specific_modules\":", specific_counts[["C3"]],
  ",\"shared_high_compatible_modules\":", shared_n,
  ",\"replicated_c2_modules\":", replicated_n,
  ",\"pairwise_specific_c2_modules\":", pairwise_specific_n,
  ",\"supportive_c2_modules\":", supportive_n, "}\n"
))
