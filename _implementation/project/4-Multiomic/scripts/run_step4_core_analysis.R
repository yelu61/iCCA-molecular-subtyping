options(stringsAsFactors = FALSE)
set.seed(20260715)

required <- c("readxl", "limma", "ggplot2", "ggrepel", "patchwork", "pheatmap", "fgsea", "msigdbr")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing R packages: ", paste(missing, collapse = ", "))

suppressPackageStartupMessages({
  library(limma)
  library(ggplot2)
  library(patchwork)
})

root <- normalizePath(getwd())
step_dir <- file.path(root, "4-Multiomic")
data_dir <- file.path(step_dir, "data")
pub_dir <- file.path(step_dir, "publication")
table_dir <- file.path(pub_dir, "tables")
fig_dir <- file.path(pub_dir, "figures")
model_dir <- file.path(pub_dir, "models")
report_dir <- file.path(pub_dir, "reports")
input_dir <- file.path(step_dir, "mofa_inputs")
for (d in c(data_dir, table_dir, fig_dir, model_dir, report_dir, input_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

write_csv <- function(x, path) write.csv(x, path, row.names = FALSE, na = "")
save_pair <- function(plot, stem, width = 8, height = 5.5) {
  ggsave(file.path(fig_dir, paste0(stem, ".pdf")), plot, width = width, height = height,
         device = cairo_pdf)
  ggsave(file.path(fig_dir, paste0(stem, ".png")), plot, width = width, height = height,
         dpi = 450)
}
load_object <- function(path, object) {
  e <- new.env(parent = emptyenv())
  load(path, envir = e)
  if (!exists(object, envir = e, inherits = FALSE)) stop("Missing object ", object, " in ", path)
  get(object, envir = e)
}
row_zscore <- function(x) {
  z <- t(scale(t(x)))
  z[!is.finite(z)] <- 0
  z
}
row_median_impute <- function(x) {
  med <- apply(x, 1, median, na.rm = TRUE)
  for (i in seq_len(nrow(x))) x[i, is.na(x[i, ])] <- med[i]
  x
}
clean_feature <- function(x) {
  x <- trimws(as.character(x))
  valid <- !is.na(x) & x != "" &
    !grepl("^[0-9]{1,2}-(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)$", x,
           ignore.case = TRUE)
  x[!valid] <- NA_character_
  x
}
aggregate_rows_mean <- function(x, groups) {
  keep <- !is.na(groups) & groups != ""
  x <- x[keep, , drop = FALSE]
  groups <- groups[keep]
  storage.mode(x) <- "numeric"
  sums <- rowsum(x, groups, reorder = FALSE, na.rm = TRUE)
  counts <- rowsum(1L * !is.na(x), groups, reorder = FALSE)
  sums / pmax(counts, 1)
}
bh <- function(p) p.adjust(p, method = "BH")

meta_path <- file.path(root, "1-Clustering/data/sample_info_with_clusters.xlsx")
expr_path <- file.path(root, "1-Clustering/data/combined_expr_after_correction.rds")
protein_path <- file.path(root, "0-Data/qc/0-Data-protein.Rdata")
maf_path <- file.path(root, "0-Data/qc/0-Data-combined_maf_data.csv")
cnv_path <- file.path(root, "0-Data/raw/2022Dong/WES_CNV_GISTIC2_raw.txt")
subtype_palette_path <- file.path(root, "docs/subtype_color_dictionary.csv")
study_palette_path <- file.path(root, "0-Data/qc/study_color_dictionary.csv")
stopifnot(all(file.exists(c(meta_path, expr_path, protein_path, maf_path, cnv_path,
                            subtype_palette_path, study_palette_path))))

metadata <- as.data.frame(readxl::read_excel(meta_path), check.names = FALSE)
stopifnot(nrow(metadata) == length(unique(metadata$Sample)))
metadata$FinalSubtype <- as.character(metadata$Ensemble)
metadata$SubtypeClass <- ifelse(metadata$FinalSubtype %in% c("C1", "C2", "C3"), "Main", "Intermediate")
subtype_dict <- read.csv(subtype_palette_path, check.names = FALSE)
subtype_palette <- setNames(subtype_dict$ColorHex, subtype_dict$FinalSubtype)
study_dict <- read.csv(study_palette_path, check.names = FALSE)
study_palette <- setNames(study_dict$ColorHex, study_dict$Study)

expr <- load_object(expr_path, "corrected_expr_limma")
stopifnot(setequal(colnames(expr), metadata$Sample))
metadata <- metadata[match(colnames(expr), metadata$Sample), , drop = FALSE]
stopifnot(identical(metadata$Sample, colnames(expr)))

input_contract <- data.frame(
  Input = c("Locked Step1 expression", "Locked Step1 subtype metadata", "Step0 harmonized mutation records",
            "Step0 protein objects", "2022Dong continuous GISTIC2 CNV", "Subtype color dictionary",
            "Study color dictionary"),
  Path = c(expr_path, meta_path, maf_path, protein_path, cnv_path, subtype_palette_path, study_palette_path),
  ReadOnly = TRUE
)
write_csv(input_contract, file.path(table_dir, "step4_input_contract.csv"))

# -----------------------------------------------------------------------------
# DNA: construct the assayed-sample roster before assigning mutation-negative calls.
# -----------------------------------------------------------------------------
maf <- read.csv(maf_path, check.names = FALSE)
maf$Study <- c(Dong2022 = "2022Dong", Deng2023 = "2023Deng", TCGA = "TCGA-CHOL")[maf$Dataset]
maf$JoinKey <- ifelse(maf$Study == "TCGA-CHOL", substr(maf$Tumor_Sample_Barcode, 1, 12),
                      maf$Tumor_Sample_Barcode)
metadata$JoinKey <- ifelse(metadata$Study == "TCGA-CHOL", substr(metadata$Sample, 1, 12), metadata$Sample)
key_map <- metadata[!duplicated(paste(metadata$Study, metadata$JoinKey)), c("Study", "JoinKey", "Sample")]
maf <- merge(maf, key_map, by = c("Study", "JoinKey"), all.x = FALSE, all.y = FALSE)

dna_roster <- rbind(
  metadata[metadata$Study == "2022Dong" & metadata[["Genome sequencing"]] == "YES", ],
  metadata[metadata$Study == "2023Deng" & metadata[["Genome sequencing"]] == "YES", ],
  metadata[metadata$Study == "TCGA-CHOL", ]
)
dna_roster <- dna_roster[!duplicated(dna_roster$Sample), ]
stopifnot(nrow(dna_roster) == 359L)

nonsilent <- maf[maf$Variant_Classification != "Silent" & !is.na(maf$Hugo_Symbol) & maf$Hugo_Symbol != "", ]
nonsilent$VariantKey <- paste(nonsilent$Hugo_Symbol, nonsilent$Chromosome,
                              nonsilent$Start_Position, nonsilent$Reference_Allele,
                              nonsilent$Tumor_Seq_Allele2, sep = ":")
burden_counts <- aggregate(VariantKey ~ Sample, nonsilent, function(z) length(unique(z)))
names(burden_counts)[2] <- "NonsynonymousVariantCount"
mutation_burden <- merge(dna_roster[, c("Sample", "Study", "FinalSubtype", "SubtypeClass")],
                         burden_counts, by = "Sample", all.x = TRUE)
mutation_burden$NonsynonymousVariantCount[is.na(mutation_burden$NonsynonymousVariantCount)] <- 0L
write_csv(mutation_burden, file.path(table_dir, "step4_nonsynonymous_variant_count_by_sample.csv"))

drivers <- c("TP53", "KRAS", "IDH1", "IDH2", "BAP1", "ARID1A", "PBRM1", "FGFR2",
             "EPHA2", "SMAD4", "STK11", "ATM", "ARID2", "NF1", "TERT")
driver_events <- unique(nonsilent[nonsilent$Hugo_Symbol %in% drivers, c("Sample", "Hugo_Symbol")])
driver_grid <- expand.grid(Sample = dna_roster$Sample, Gene = drivers, stringsAsFactors = FALSE)
driver_grid$Mutated <- paste(driver_grid$Sample, driver_grid$Gene) %in%
  paste(driver_events$Sample, driver_events$Hugo_Symbol)
driver_grid <- merge(driver_grid,
                     dna_roster[, c("Sample", "Study", "FinalSubtype", "SubtypeClass")],
                     by = "Sample", all.x = TRUE)
write_csv(driver_grid, file.path(table_dir, "step4_driver_snv_indel_sample_long.csv"))

main_driver <- driver_grid[driver_grid$SubtypeClass == "Main", ]
driver_freq <- aggregate(Mutated ~ Gene + FinalSubtype, main_driver,
                         function(z) c(N = length(z), Altered = sum(z), Fraction = mean(z)))
driver_freq <- data.frame(Gene = driver_freq$Gene, FinalSubtype = driver_freq$FinalSubtype,
                          N = driver_freq$Mutated[, "N"], Altered = driver_freq$Mutated[, "Altered"],
                          Fraction = driver_freq$Mutated[, "Fraction"])
write_csv(driver_freq, file.path(table_dir, "step4_driver_snv_indel_frequency_main_subtypes.csv"))

driver_tests <- do.call(rbind, lapply(c("C1", "C2", "C3"), function(st) {
  do.call(rbind, lapply(drivers, function(g) {
    z <- main_driver[main_driver$Gene == g, ]
    target <- z$FinalSubtype == st
    crude <- fisher.test(table(factor(target, levels = c(TRUE, FALSE)),
                               factor(z$Mutated, levels = c(TRUE, FALSE))))
    arr <- array(0, dim = c(2, 2, length(unique(z$Study))),
                 dimnames = list(Target = c("Target", "Others"), Altered = c("Yes", "No"),
                                 Study = unique(z$Study)))
    for (s in dimnames(arr)$Study) {
      zz <- z[z$Study == s, ]
      arr[, , s] <- table(factor(zz$FinalSubtype == st, levels = c(TRUE, FALSE)),
                           factor(zz$Mutated, levels = c(TRUE, FALSE)))
    }
    cmh <- tryCatch(mantelhaen.test(arr, correct = FALSE, conf.level = 0.95), error = function(e) NULL)
    data.frame(Subtype = st, Gene = g,
               TargetN = sum(target), OtherN = sum(!target),
               TargetAltered = sum(z$Mutated[target]), OtherAltered = sum(z$Mutated[!target]),
               CrudeOddsRatio = unname(crude$estimate), CrudePValue = crude$p.value,
               CMHOddsRatio = if (is.null(cmh)) NA_real_ else unname(cmh$estimate),
               CMHLower95 = if (is.null(cmh) || is.null(cmh$conf.int)) NA_real_ else unname(cmh$conf.int[1]),
               CMHUpper95 = if (is.null(cmh) || is.null(cmh$conf.int)) NA_real_ else unname(cmh$conf.int[2]),
               CMHPValue = if (is.null(cmh)) NA_real_ else cmh$p.value)
  }))
}))
driver_tests$CMHFDR <- ave(driver_tests$CMHPValue, driver_tests$Subtype, FUN = bh)
write_csv(driver_tests, file.path(table_dir, "step4_driver_snv_indel_study_stratified_tests.csv"))

# Keep the figure annotation synchronized with the actual multiplicity-corrected
# results.  The previous hard-coded statement that no contrast was significant
# became stale when the full assay roster (including mutation-negative samples)
# was restored.  FDR is controlled within each subtype across the prespecified
# driver panel; the global count is reported in the publication table/summary.
driver_fdr_n <- sum(driver_tests$CMHFDR < 0.05, na.rm = TRUE)
driver_fdr_caption <- if (driver_fdr_n == 0L) {
  "No prespecified contrast met FDR < 0.05; estimates are descriptive."
} else {
  sprintf("%d of %d prespecified subtype-versus-other contrasts met within-subtype FDR < 0.05; interpret with study and panel limitations.",
          driver_fdr_n, nrow(driver_tests))
}

p_driver <- ggplot(driver_freq, aes(FinalSubtype, Gene)) +
  geom_point(aes(size = Fraction, fill = FinalSubtype), shape = 21, color = "grey20", stroke = 0.35) +
  geom_text(aes(label = ifelse(Altered > 0, sprintf("%d", Altered), "")), size = 2.4) +
  scale_fill_manual(values = subtype_palette[c("C1", "C2", "C3")], guide = "none") +
  scale_size_continuous(range = c(1.5, 8), labels = scales::percent, name = "Altered fraction") +
  labs(title = "Driver SNV/indel frequencies across principal subtypes",
       subtitle = sprintf("DNA-assayed main-state tumors; full assay roster retained (n=%d)",
                          sum(dna_roster$SubtypeClass == "Main")),
       x = NULL, y = NULL,
       caption = "FGFR2 denotes SNV/indel only; fusion status is not represented in this table.") +
  theme_classic(base_size = 10) + theme(legend.position = "top")
save_pair(p_driver, "step4_driver_snv_indel_frequency_main_subtypes", 7.8, 7.2)

driver_forest <- driver_tests[is.finite(driver_tests$CMHOddsRatio) & driver_tests$CMHOddsRatio > 0 &
                                is.finite(driver_tests$CMHLower95) & driver_tests$CMHLower95 > 0 &
                                is.finite(driver_tests$CMHUpper95), ]
driver_forest$Gene <- factor(driver_forest$Gene, levels = rev(drivers))
p_driver_forest <- ggplot(driver_forest,
                          aes(CMHOddsRatio, Gene, xmin = CMHLower95, xmax = CMHUpper95,
                              colour = Subtype)) +
  geom_vline(xintercept = 1, linetype = 2, colour = "grey55", linewidth = 0.45) +
  geom_errorbar(orientation = "y", width = 0.15, linewidth = 0.45) +
  geom_point(size = 1.8) +
  facet_wrap(~Subtype, nrow = 1) +
  scale_x_log10() +
  scale_colour_manual(values = subtype_palette[c("C1", "C2", "C3")], guide = "none") +
  labs(title = "Study-stratified driver SNV/indel associations",
       subtitle = "Cochran–Mantel–Haenszel odds ratios and 95% confidence intervals",
       x = "Odds ratio (log scale; subtype versus other principal states)", y = NULL,
       caption = driver_fdr_caption) +
  theme_classic(base_size = 9)
save_pair(p_driver_forest, "step4_driver_snv_indel_odds_ratios", 10.5, 6.2)

burden_main <- mutation_burden[mutation_burden$SubtypeClass == "Main", ]
p_burden <- ggplot(burden_main, aes(FinalSubtype, NonsynonymousVariantCount, color = FinalSubtype)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, color = "grey35") +
  geom_jitter(width = 0.15, alpha = 0.45, size = 1) +
  facet_wrap(~Study, scales = "free_y") +
  scale_color_manual(values = subtype_palette[c("C1", "C2", "C3")], guide = "none") +
  labs(title = "Nonsynonymous variant counts within DNA-assayed studies",
       subtitle = "Counts are not tumor mutational burden because capture-size harmonization is unavailable",
       x = NULL, y = "Nonsynonymous variants per tumor") +
  theme_classic(base_size = 10)
save_pair(p_burden, "step4_nonsynonymous_variant_counts_by_subtype_and_study", 10, 5.5)

# -----------------------------------------------------------------------------
# 2022Dong continuous CNV: explicit cohort-specific denominator.
# -----------------------------------------------------------------------------
cnv_raw <- read.delim(cnv_path, check.names = FALSE)
cnv_gene <- clean_feature(sub("\\|.*", "", cnv_raw[[1]]))
cnv_numeric <- as.matrix(cnv_raw[, -(1:3), drop = FALSE])
storage.mode(cnv_numeric) <- "numeric"
cnv <- aggregate_rows_mean(cnv_numeric, cnv_gene)
colnames(cnv) <- sub("^X", "", colnames(cnv))
cnv_samples <- intersect(colnames(cnv), dna_roster$Sample[dna_roster$Study == "2022Dong"])
cnv <- cnv[, cnv_samples, drop = FALSE]
cnv_meta <- metadata[match(cnv_samples, metadata$Sample), c("Sample", "Study", "FinalSubtype", "SubtypeClass")]
stopifnot(identical(cnv_meta$Sample, colnames(cnv)))

cnv_burden <- data.frame(Sample = cnv_samples, Study = "2022Dong",
                         FinalSubtype = cnv_meta$FinalSubtype, SubtypeClass = cnv_meta$SubtypeClass,
                         MeanAbsoluteCNV = colMeans(abs(cnv), na.rm = TRUE))
write_csv(cnv_burden, file.path(table_dir, "step4_cnv_burden_2022dong.csv"))

cnv_driver_genes <- intersect(drivers, rownames(cnv))
cnv_driver <- as.data.frame(as.table(cnv[cnv_driver_genes, , drop = FALSE]), stringsAsFactors = FALSE)
names(cnv_driver) <- c("Gene", "Sample", "ContinuousCNV")
cnv_driver <- merge(cnv_driver, cnv_meta, by = "Sample", all.x = TRUE)
cnv_driver$Call <- ifelse(cnv_driver$ContinuousCNV >= 0.3, "Gain",
                          ifelse(cnv_driver$ContinuousCNV <= -0.3, "Loss", "Neutral"))
write_csv(cnv_driver, file.path(table_dir, "step4_driver_cnv_sample_long_2022dong.csv"))

cnv_freq <- aggregate(Call ~ Gene + FinalSubtype,
                      cnv_driver[cnv_driver$SubtypeClass == "Main", ],
                      function(z) c(N = length(z), Gain = sum(z == "Gain"), Loss = sum(z == "Loss")))
cnv_freq <- data.frame(Gene = cnv_freq$Gene, FinalSubtype = cnv_freq$FinalSubtype,
                       N = cnv_freq$Call[, "N"], Gain = cnv_freq$Call[, "Gain"], Loss = cnv_freq$Call[, "Loss"])
cnv_freq$GainFraction <- cnv_freq$Gain / cnv_freq$N
cnv_freq$LossFraction <- cnv_freq$Loss / cnv_freq$N
write_csv(cnv_freq, file.path(table_dir, "step4_driver_cnv_frequency_main_subtypes_2022dong.csv"))

cnv_plot <- rbind(
  data.frame(cnv_freq[, c("Gene", "FinalSubtype", "N")], Direction = "Gain", Fraction = cnv_freq$GainFraction),
  data.frame(cnv_freq[, c("Gene", "FinalSubtype", "N")], Direction = "Loss", Fraction = cnv_freq$LossFraction)
)
p_cnv <- ggplot(cnv_plot, aes(FinalSubtype, Gene)) +
  geom_point(aes(size = Fraction, fill = FinalSubtype), shape = 21, color = "grey25", stroke = 0.3) +
  facet_wrap(~Direction) +
  scale_fill_manual(values = subtype_palette[c("C1", "C2", "C3")], guide = "none") +
  scale_size_continuous(range = c(1.2, 7), labels = scales::percent, name = "Sample fraction") +
  labs(title = "Driver-gene CNV calls in the 2022Dong cohort",
       subtitle = sprintf("Continuous GISTIC2 values; exploratory gain/loss thresholds ±0.3; main-state n=%d",
                          sum(cnv_meta$SubtypeClass == "Main")),
       x = NULL, y = NULL) + theme_classic(base_size = 10) + theme(legend.position = "top")
save_pair(p_cnv, "step4_driver_cnv_frequency_2022dong", 9.2, 6.8)

# -----------------------------------------------------------------------------
# Mass-spectrometry proteomics: repair the source-scale mismatch without changing
# protected Step0 objects. The ComBat output is converted to within-study protein
# z scores, which removes any positive affine batch transformation per protein.
# Only proteins observed in at least half of both source studies are retained.
# -----------------------------------------------------------------------------
protein_env <- new.env(parent = emptyenv())
load(protein_path, envir = protein_env)
corrected_protein <- protein_env$corrected_prot_matrix
raw_dong_protein <- protein_env$pro_2022Dong
raw_deng_protein <- protein_env$pro_2023Deng
tcga_rppa <- protein_env$pro_tcga

coverage_by_gene <- function(df, gene_col, zero_is_missing = FALSE) {
  genes <- clean_feature(df[[gene_col]])
  x <- as.matrix(df[, setdiff(names(df), gene_col), drop = FALSE])
  storage.mode(x) <- "numeric"
  if (zero_is_missing) x[x <= 1e-5] <- NA_real_
  cov <- rowMeans(!is.na(x))
  z <- data.frame(Gene = genes, Coverage = cov)
  z <- z[!is.na(z$Gene), ]
  aggregate(Coverage ~ Gene, z, max)
}
cov_dong <- coverage_by_gene(raw_dong_protein, "Gene", FALSE)
cov_deng <- coverage_by_gene(raw_deng_protein, "Tumor..T..proteomic.ID", TRUE)
coverage <- merge(cov_dong, cov_deng, by = "Gene", suffixes = c("_2022Dong", "_2023Deng"))
coverage$Retained <- coverage$Coverage_2022Dong >= 0.5 & coverage$Coverage_2023Deng >= 0.5 &
  coverage$Gene %in% rownames(corrected_protein)
write_csv(coverage, file.path(table_dir, "step4_protein_feature_coverage_by_study.csv"))

protein_samples <- intersect(colnames(corrected_protein), metadata$Sample)
protein_meta <- metadata[match(protein_samples, metadata$Sample),
                         c("Sample", "Study", "FinalSubtype", "SubtypeClass")]
protein_meta <- protein_meta[protein_meta$Study %in% c("2022Dong", "2023Deng"), ]
protein_samples <- protein_meta$Sample
stopifnot(identical(protein_samples, colnames(corrected_protein)[match(protein_samples, colnames(corrected_protein))]))
reliable_proteins <- coverage$Gene[coverage$Retained]
protein <- corrected_protein[reliable_proteins, protein_samples, drop = FALSE]
protein_standardized <- protein
for (study in c("2022Dong", "2023Deng")) {
  jj <- protein_meta$Study == study
  protein_standardized[, jj] <- row_zscore(protein[, jj, drop = FALSE])
}
keep_variation <- apply(protein_standardized, 1, sd) > 0
protein_standardized <- protein_standardized[keep_variation, , drop = FALSE]
saveRDS(protein_standardized, file.path(data_dir, "step4_mass_spec_within_study_standardized.rds"))

scale_audit <- data.frame(
  Matrix = c("Step0 ComBat object", "Step4 within-study standardized object"),
  Samples = c(ncol(corrected_protein), ncol(protein_standardized)),
  Proteins = c(nrow(corrected_protein), nrow(protein_standardized)),
  Minimum = c(min(corrected_protein), min(protein_standardized)),
  Median = c(median(corrected_protein), median(protein_standardized)),
  Maximum = c(max(corrected_protein), max(protein_standardized)),
  Interpretation = c("Not used for cross-study effect-size interpretation because source scales differ",
                     "Used for Step4 standardized-effect inference")
)
write_csv(scale_audit, file.path(table_dir, "step4_protein_scale_audit.csv"))

main_protein_meta <- protein_meta[protein_meta$SubtypeClass == "Main", ]
main_protein <- protein_standardized[, main_protein_meta$Sample, drop = FALSE]
main_protein_meta$FinalSubtype <- factor(main_protein_meta$FinalSubtype, levels = c("C1", "C2", "C3"))
main_protein_meta$Study <- factor(main_protein_meta$Study)
design_pro <- model.matrix(~0 + FinalSubtype + Study, data = main_protein_meta)
colnames(design_pro) <- make.names(colnames(design_pro))
contrast_pro <- makeContrasts(
  C1 = FinalSubtypeC1 - (FinalSubtypeC2 + FinalSubtypeC3) / 2,
  C2 = FinalSubtypeC2 - (FinalSubtypeC1 + FinalSubtypeC3) / 2,
  C3 = FinalSubtypeC3 - (FinalSubtypeC1 + FinalSubtypeC2) / 2,
  levels = design_pro
)
fit_pro <- eBayes(contrasts.fit(lmFit(main_protein, design_pro), contrast_pro), robust = TRUE)
protein_deg <- do.call(rbind, lapply(c("C1", "C2", "C3"), function(st) {
  tt <- topTable(fit_pro, coef = st, number = Inf, sort.by = "none")
  data.frame(Subtype = st, Protein = rownames(tt), StandardizedDifference = tt$logFC,
             ModeratedT = tt$t, PValue = tt$P.Value, FDR = tt$adj.P.Val,
             stringsAsFactors = FALSE)
}))
protein_deg$Direction <- ifelse(protein_deg$FDR < 0.05 & protein_deg$StandardizedDifference > 0,
                                "Subtype higher",
                                ifelse(protein_deg$FDR < 0.05 & protein_deg$StandardizedDifference < 0,
                                       "Subtype lower", "Not significant"))
write_csv(protein_deg, file.path(table_dir, "step4_mass_spec_study_adjusted_deg_main_subtypes.csv"))

label_pro <- do.call(rbind, lapply(split(protein_deg, protein_deg$Subtype), function(z) {
  z <- z[z$Direction != "Not significant", ]
  head(z[order(z$FDR, -abs(z$StandardizedDifference)), ], 8)
}))
p_pro_volcano <- ggplot(protein_deg,
                        aes(StandardizedDifference, -log10(pmax(FDR, 1e-300)), color = Direction)) +
  geom_point(alpha = 0.55, size = 0.8) +
  geom_vline(xintercept = 0, color = "grey70") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey60") +
  ggrepel::geom_text_repel(data = label_pro, aes(label = Protein), size = 2.4,
                           max.overlaps = Inf, show.legend = FALSE) +
  facet_wrap(~Subtype, scales = "free") +
  scale_color_manual(values = c("Subtype higher" = "#D55E00", "Subtype lower" = "#0072B2",
                                "Not significant" = "grey80")) +
  labs(title = "Study-adjusted standardized protein differences",
       subtitle = sprintf("Main-state mass-spectrometry cohort n=%d; positive values indicate subtype higher",
                          nrow(main_protein_meta)),
       x = "Within-study standardized difference", y = expression(-log[10](FDR)), color = NULL) +
  theme_classic(base_size = 10) + theme(legend.position = "top")
save_pair(p_pro_volcano, "step4_mass_spec_deg_main_subtypes", 11, 5.6)

heat_proteins <- unique(unlist(lapply(c("C1", "C2", "C3"), function(st) {
  z <- protein_deg[protein_deg$Subtype == st & protein_deg$StandardizedDifference > 0, ]
  head(z$Protein[order(z$FDR, -z$StandardizedDifference)], 10)
})))
heat_proteins <- heat_proteins[heat_proteins %in% rownames(main_protein)]
heat_annotation <- data.frame(Subtype = as.character(main_protein_meta$FinalSubtype),
                              Study = as.character(main_protein_meta$Study),
                              row.names = main_protein_meta$Sample)
heat_colors <- list(Subtype = subtype_palette[c("C1", "C2", "C3")],
                    Study = study_palette[c("2022Dong", "2023Deng")])
for (ext in c("pdf", "png")) {
  pheatmap::pheatmap(main_protein[heat_proteins, , drop = FALSE], show_colnames = FALSE,
                     annotation_col = heat_annotation, annotation_colors = heat_colors,
                     cluster_cols = TRUE, clustering_distance_cols = "correlation",
                     border_color = NA, fontsize_row = 7,
                     main = "Study-standardized protein features across principal subtypes",
                     filename = file.path(fig_dir, paste0("step4_mass_spec_marker_heatmap.", ext)),
                     width = 10, height = 8)
}

# Ranked Hallmark/Reactome protein enrichment.
msig <- rbind(
  transform(msigdbr::msigdbr(species = "Homo sapiens", collection = "H")[, c("gs_name", "gene_symbol")],
            Collection = "Hallmark"),
  transform(msigdbr::msigdbr(species = "Homo sapiens", collection = "C2",
                             subcollection = "CP:REACTOME")[, c("gs_name", "gene_symbol")],
            Collection = "Reactome")
)
pathways <- split(msig$gene_symbol, msig$gs_name)
protein_gsea <- do.call(rbind, lapply(c("C1", "C2", "C3"), function(st) {
  z <- protein_deg[protein_deg$Subtype == st, ]
  ranks <- sort(setNames(z$ModeratedT, z$Protein), decreasing = TRUE)
  fg <- as.data.frame(fgsea::fgseaMultilevel(pathways, ranks, minSize = 15, maxSize = 500,
                                             eps = 0, nPermSimple = 10000))
  if (!nrow(fg)) return(NULL)
  fg$Subtype <- st
  fg$Collection <- ifelse(grepl("^HALLMARK_", fg$pathway), "Hallmark", "Reactome")
  fg$leadingEdge <- vapply(fg$leadingEdge, paste, collapse = ";", FUN.VALUE = character(1))
  fg
}))
protein_gsea <- protein_gsea[order(protein_gsea$Subtype, protein_gsea$padj, -abs(protein_gsea$NES)), ]
write_csv(protein_gsea, file.path(table_dir, "step4_mass_spec_ranked_gsea_main_subtypes.csv"))

# RNA-protein concordance on exactly matched samples, summarized within study.
shared_samples <- intersect(colnames(protein_standardized), colnames(expr))
shared_genes <- intersect(rownames(protein_standardized), rownames(expr))
concordance_by_study <- do.call(rbind, lapply(c("2022Dong", "2023Deng"), function(study) {
  samples <- intersect(shared_samples, protein_meta$Sample[protein_meta$Study == study])
  rho <- vapply(shared_genes, function(g) {
    suppressWarnings(cor(expr[g, samples], protein_standardized[g, samples], method = "spearman"))
  }, numeric(1))
  rho <- pmax(pmin(rho, 0.999999), -0.999999)
  data.frame(Gene = shared_genes, Study = study, N = length(samples), SpearmanRho = rho,
             FisherZ = atanh(rho), Weight = pmax(length(samples) - 3, 1))
}))
concordance_meta <- do.call(rbind, lapply(split(concordance_by_study, concordance_by_study$Gene), function(z) {
  mz <- weighted.mean(z$FisherZ, z$Weight)
  se <- 1 / sqrt(sum(z$Weight))
  data.frame(Gene = z$Gene[1], Studies = nrow(z), TotalN = sum(z$N),
             MetaSpearmanRho = tanh(mz), Lower95 = tanh(mz - 1.96 * se),
             Upper95 = tanh(mz + 1.96 * se), PValue = 2 * pnorm(-abs(mz / se)))
}))
concordance_meta$FDR <- bh(concordance_meta$PValue)
concordance_meta <- concordance_meta[order(-concordance_meta$MetaSpearmanRho), ]
write_csv(concordance_by_study, file.path(table_dir, "step4_rna_protein_gene_correlations_by_study.csv"))
write_csv(concordance_meta, file.path(table_dir, "step4_rna_protein_gene_correlation_meta.csv"))

# Matched-sample subtype contrasts in RNA and protein.
matched_main_meta <- protein_meta[protein_meta$SubtypeClass == "Main" & protein_meta$Sample %in% shared_samples, ]
matched_main_meta$FinalSubtype <- factor(matched_main_meta$FinalSubtype, levels = c("C1", "C2", "C3"))
matched_main_meta$Study <- factor(matched_main_meta$Study)
design_matched <- model.matrix(~0 + FinalSubtype + Study, data = matched_main_meta)
colnames(design_matched) <- make.names(colnames(design_matched))
contrast_matched <- makeContrasts(
  C1 = FinalSubtypeC1 - (FinalSubtypeC2 + FinalSubtypeC3) / 2,
  C2 = FinalSubtypeC2 - (FinalSubtypeC1 + FinalSubtypeC3) / 2,
  C3 = FinalSubtypeC3 - (FinalSubtypeC1 + FinalSubtypeC2) / 2,
  levels = design_matched
)
fit_rna_matched <- eBayes(contrasts.fit(lmFit(expr[shared_genes, matched_main_meta$Sample], design_matched),
                                         contrast_matched), robust = TRUE)
fit_pro_matched <- eBayes(contrasts.fit(lmFit(protein_standardized[shared_genes, matched_main_meta$Sample],
                                              design_matched), contrast_matched), robust = TRUE)
contrast_concordance <- do.call(rbind, lapply(c("C1", "C2", "C3"), function(st) {
  rna_t <- fit_rna_matched$t[, st]
  pro_t <- fit_pro_matched$t[, st]
  data.frame(Subtype = st, Gene = shared_genes, RNAModeratedT = rna_t,
             ProteinModeratedT = pro_t)
}))
contrast_summary <- do.call(rbind, lapply(split(contrast_concordance, contrast_concordance$Subtype), function(z) {
  ct <- suppressWarnings(cor.test(z$RNAModeratedT, z$ProteinModeratedT, method = "spearman", exact = FALSE))
  data.frame(Subtype = z$Subtype[1], Genes = nrow(z), SpearmanRho = unname(ct$estimate), PValue = ct$p.value)
}))
contrast_summary$FDR <- bh(contrast_summary$PValue)
write_csv(contrast_concordance, file.path(table_dir, "step4_matched_rna_protein_subtype_contrasts_long.csv"))
write_csv(contrast_summary, file.path(table_dir, "step4_matched_rna_protein_subtype_concordance.csv"))

p_concordance <- ggplot(contrast_concordance, aes(RNAModeratedT, ProteinModeratedT)) +
  geom_hline(yintercept = 0, color = "grey75") + geom_vline(xintercept = 0, color = "grey75") +
  geom_point(alpha = 0.22, size = 0.55, color = "#3C5488") +
  geom_smooth(method = "lm", se = FALSE, color = "#D55E00", linewidth = 0.7) +
  facet_wrap(~Subtype, scales = "free") +
  labs(title = "Matched RNA-protein subtype contrasts",
       subtitle = sprintf("Study-adjusted moderated t statistics; matched main-state samples n=%d",
                          nrow(matched_main_meta)),
       x = "RNA moderated t statistic", y = "Protein moderated t statistic") +
  theme_classic(base_size = 10)
save_pair(p_concordance, "step4_matched_rna_protein_contrast_concordance", 10.5, 4.5)

# Pathway-level RNA-protein concordance using the same matched-sample ranks.
matched_pathways <- do.call(rbind, lapply(c("C1", "C2", "C3"), function(st) {
  run_view <- function(stats, view) {
    fg <- as.data.frame(fgsea::fgseaMultilevel(pathways, sort(stats, decreasing = TRUE),
                                               minSize = 15, maxSize = 500, eps = 0,
                                               nPermSimple = 10000))
    fg$View <- view
    fg
  }
  rna <- run_view(setNames(fit_rna_matched$t[, st], shared_genes), "RNA")
  pro <- run_view(setNames(fit_pro_matched$t[, st], shared_genes), "Protein")
  zz <- merge(rna[, c("pathway", "NES", "padj")], pro[, c("pathway", "NES", "padj")],
              by = "pathway", suffixes = c("_RNA", "_Protein"))
  zz$Subtype <- st
  zz$Collection <- ifelse(grepl("^HALLMARK_", zz$pathway), "Hallmark", "Reactome")
  zz
}))
write_csv(matched_pathways, file.path(table_dir, "step4_matched_rna_protein_pathway_concordance.csv"))
pathway_summary <- do.call(rbind, lapply(split(matched_pathways,
                                               interaction(matched_pathways$Subtype, matched_pathways$Collection)),
                                        function(z) {
  ct <- suppressWarnings(cor.test(z$NES_RNA, z$NES_Protein, method = "spearman", exact = FALSE))
  data.frame(Subtype = z$Subtype[1], Collection = z$Collection[1], Pathways = nrow(z),
             SpearmanRho = unname(ct$estimate), PValue = ct$p.value)
}))
pathway_summary$FDR <- bh(pathway_summary$PValue)
write_csv(pathway_summary, file.path(table_dir, "step4_matched_rna_protein_pathway_concordance_summary.csv"))

plot_pathways <- matched_pathways[matched_pathways$Collection == "Hallmark", ]
p_pathway <- ggplot(plot_pathways, aes(NES_RNA, NES_Protein, color = Subtype)) +
  geom_hline(yintercept = 0, color = "grey75") + geom_vline(xintercept = 0, color = "grey75") +
  geom_point(alpha = 0.75, size = 1.5) +
  facet_wrap(~Subtype) +
  scale_color_manual(values = subtype_palette[c("C1", "C2", "C3")], guide = "none") +
  labs(title = "Matched Hallmark pathway concordance between RNA and protein",
       subtitle = "Normalized enrichment scores from the same samples and explicit subtype-minus-others contrasts",
       x = "RNA NES", y = "Protein NES") + theme_classic(base_size = 10)
save_pair(p_pathway, "step4_matched_rna_protein_hallmark_concordance", 10.5, 4.5)

# -----------------------------------------------------------------------------
# TCGA-CHOL RPPA remains a separate, small validation layer.
# -----------------------------------------------------------------------------
rppa_samples_raw <- names(tcga_rppa)[6:ncol(tcga_rppa)]
rppa_samples <- gsub("\\.", "-", rppa_samples_raw)
rppa_keep <- rppa_samples %in% metadata$Sample[metadata$Study == "TCGA-CHOL"]
rppa_x <- as.matrix(tcga_rppa[, 6:ncol(tcga_rppa), drop = FALSE])[, rppa_keep, drop = FALSE]
storage.mode(rppa_x) <- "numeric"
colnames(rppa_x) <- rppa_samples[rppa_keep]
rppa_targets <- clean_feature(tcga_rppa$peptide_target)
rppa_x <- aggregate_rows_mean(rppa_x, rppa_targets)
rppa_x <- rppa_x[rowSums(is.finite(rppa_x)) > 0, , drop = FALSE]
rppa_meta <- metadata[match(colnames(rppa_x), metadata$Sample),
                      c("Sample", "Study", "FinalSubtype", "SubtypeClass")]
stopifnot(ncol(rppa_x) == 25L, identical(rppa_meta$Sample, colnames(rppa_x)))
write_csv(rppa_meta, file.path(table_dir, "step4_tcga_rppa_matched_samples.csv"))

rppa_main_meta <- rppa_meta[rppa_meta$SubtypeClass == "Main", ]
rppa_main <- rppa_x[, rppa_main_meta$Sample, drop = FALSE]
rppa_main_meta$FinalSubtype <- factor(rppa_main_meta$FinalSubtype, levels = c("C1", "C2", "C3"))
design_rppa <- model.matrix(~0 + FinalSubtype, data = rppa_main_meta)
colnames(design_rppa) <- make.names(colnames(design_rppa))
contrast_rppa <- makeContrasts(
  C1 = FinalSubtypeC1 - (FinalSubtypeC2 + FinalSubtypeC3) / 2,
  C2 = FinalSubtypeC2 - (FinalSubtypeC1 + FinalSubtypeC3) / 2,
  C3 = FinalSubtypeC3 - (FinalSubtypeC1 + FinalSubtypeC2) / 2,
  levels = design_rppa
)
fit_rppa <- eBayes(contrasts.fit(lmFit(rppa_main, design_rppa), contrast_rppa), robust = TRUE)
rppa_results <- do.call(rbind, lapply(c("C1", "C2", "C3"), function(st) {
  tt <- topTable(fit_rppa, coef = st, number = Inf, sort.by = "none")
  data.frame(Subtype = st,
             SubtypeN = sum(rppa_main_meta$FinalSubtype == st),
             OtherN = sum(rppa_main_meta$FinalSubtype != st),
             Target = rownames(tt), Difference = tt$logFC,
             ModeratedT = tt$t, PValue = tt$P.Value, FDR = tt$adj.P.Val)
}))
write_csv(rppa_results, file.path(table_dir, "step4_tcga_rppa_subtype_associations.csv"))

rppa_scaled <- t(scale(t(rppa_x)))
rppa_scaled[!is.finite(rppa_scaled)] <- 0
rppa_scaled <- rppa_scaled[apply(rppa_scaled, 1, sd) > 0, , drop = FALSE]
rppa_pca <- stats::prcomp(t(rppa_scaled), center = FALSE, scale. = FALSE)
rppa_scores <- data.frame(Sample = rownames(rppa_pca$x), PC1 = rppa_pca$x[, 1], PC2 = rppa_pca$x[, 2])
rppa_scores <- merge(rppa_scores, rppa_meta, by = "Sample", all.x = TRUE)
write_csv(rppa_scores, file.path(table_dir, "step4_tcga_rppa_pca_scores.csv"))
rppa_scores$SubtypeLabel <- subtype_dict$DisplayLabel[match(rppa_scores$FinalSubtype, subtype_dict$FinalSubtype)]
rppa_label_palette <- setNames(subtype_dict$ColorHex, subtype_dict$DisplayLabel)
p_rppa <- ggplot(rppa_scores, aes(PC1, PC2, color = SubtypeLabel, shape = SubtypeClass)) +
  geom_point(size = 3, alpha = 0.9) +
  scale_color_manual(values = rppa_label_palette) +
  scale_shape_manual(values = c(Main = 16, Intermediate = 1)) +
  labs(title = "TCGA-CHOL RPPA sample projection",
       subtitle = sprintf("Expression-matched RPPA tumors n=%d; separate from mass spectrometry", nrow(rppa_scores)),
       x = sprintf("PC1 (%.1f%%)", 100 * summary(rppa_pca)$importance[2, 1]),
       y = sprintf("PC2 (%.1f%%)", 100 * summary(rppa_pca)$importance[2, 2]),
       color = "Subtype", shape = NULL) + theme_classic(base_size = 10) +
  theme(legend.position = "top")
save_pair(p_rppa, "step4_tcga_rppa_pca", 7.2, 5.6)

# -----------------------------------------------------------------------------
# Prepare deterministic 2022Dong RNA/protein/CNV inputs for MOFA.
# -----------------------------------------------------------------------------
rna_qc_env <- new.env(parent = emptyenv())
load(file.path(root, "0-Data/qc/0-Data-rnaseq.Rdata"), envir = rna_qc_env)
rna_dong <- as.matrix(rna_qc_env$log2tpm_2022Dong)
protein_dong_gene <- clean_feature(raw_dong_protein$Gene)
protein_dong_raw <- as.matrix(raw_dong_protein[, -1, drop = FALSE])
storage.mode(protein_dong_raw) <- "numeric"
colnames(protein_dong_raw) <- sub("^X", "", colnames(protein_dong_raw))
protein_dong <- aggregate_rows_mean(protein_dong_raw, protein_dong_gene)

triple_samples <- Reduce(intersect, list(colnames(rna_dong), colnames(protein_dong), colnames(cnv), metadata$Sample))
triple_samples <- metadata$Sample[metadata$Sample %in% triple_samples & metadata$Study == "2022Dong"]
stopifnot(length(triple_samples) == 204L)
triple_meta <- metadata[match(triple_samples, metadata$Sample),
                        c("Sample", "Study", "FinalSubtype", "SubtypeClass")]
stopifnot(identical(triple_meta$Sample, triple_samples), sum(triple_meta$SubtypeClass == "Main") == 182L)

select_view <- function(x, samples, view, top_n = 1000L, min_coverage = 0.8) {
  x <- x[, samples, drop = FALSE]
  observed <- rowMeans(is.finite(x) & !is.na(x))
  keep <- observed >= min_coverage
  x <- x[keep, , drop = FALSE]
  observed <- observed[keep]
  x[!is.finite(x)] <- NA_real_
  x <- row_median_impute(x)
  variance <- apply(x, 1, var)
  keep <- is.finite(variance) & variance > 0
  x <- x[keep, , drop = FALSE]
  variance <- variance[keep]
  observed <- observed[keep]
  ord <- order(variance, decreasing = TRUE)
  ord <- head(ord, min(top_n, length(ord)))
  x <- x[ord, , drop = FALSE]
  info <- data.frame(View = view, Feature = paste0(view, "::", rownames(x)),
                     OriginalFeature = rownames(x), Variance = variance[ord],
                     ObservedFraction = observed[ord], Rank = seq_along(ord))
  rownames(x) <- info$Feature
  list(matrix = row_zscore(x), info = info)
}

mofa_rna <- select_view(rna_dong, triple_samples, "RNA")
mofa_protein <- select_view(protein_dong, triple_samples, "Protein")
mofa_cnv <- select_view(cnv, triple_samples, "CNV")
feature_selection <- rbind(mofa_rna$info, mofa_protein$info, mofa_cnv$info)
write_csv(feature_selection, file.path(table_dir, "step4_mofa_feature_selection.csv"))
write_csv(triple_meta, file.path(input_dir, "step4_mofa_sample_metadata.csv"))

write_mofa_input <- function(obj, view) {
  out <- data.frame(Sample = colnames(obj$matrix), t(obj$matrix), check.names = FALSE)
  write_csv(out, file.path(input_dir, paste0("step4_mofa_input_", view, ".csv")))
}
write_mofa_input(mofa_rna, "RNA")
write_mofa_input(mofa_protein, "Protein")
write_mofa_input(mofa_cnv, "CNV")

assay_inclusion <- rbind(
  data.frame(Analysis = "DNA SNV/indel", Assay = "DNA-seq", Study = c("2022Dong", "2023Deng", "TCGA-CHOL"),
             TotalN = c(249, 79, 31), MainStateN = vapply(c("2022Dong", "2023Deng", "TCGA-CHOL"),
                                                         function(s) sum(dna_roster$Study == s & dna_roster$SubtypeClass == "Main"), numeric(1))),
  data.frame(Analysis = "Mass-spectrometry subtype features", Assay = "Protein mass spectrometry",
             Study = c("2022Dong", "2023Deng"),
             TotalN = vapply(c("2022Dong", "2023Deng"), function(s) sum(protein_meta$Study == s), numeric(1)),
             MainStateN = vapply(c("2022Dong", "2023Deng"),
                                 function(s) sum(protein_meta$Study == s & protein_meta$SubtypeClass == "Main"), numeric(1))),
  data.frame(Analysis = "TCGA RPPA validation", Assay = "Protein RPPA", Study = "TCGA-CHOL",
             TotalN = nrow(rppa_meta), MainStateN = sum(rppa_meta$SubtypeClass == "Main")),
  data.frame(Analysis = "MOFA matched integration", Assay = "RNA + protein + CNV", Study = "2022Dong",
             TotalN = nrow(triple_meta), MainStateN = sum(triple_meta$SubtypeClass == "Main"))
)
write_csv(assay_inclusion, file.path(table_dir, "step4_assay_specific_inclusion.csv"))

core_audit <- data.frame(
  Finding = c("Historical notebook executes out of order", "Historical mutation burden mislabeled as TMB",
              "Historical protein sample order was not identifier-matched", "Historical protein source scales differ",
              "Historical multi-omics feature selection was random", "Historical iCluster cluster-count interpretation",
              "Historical MOFA interface lacked weights, variance, and convergence outputs",
              "Historical TCGA RPPA omission"),
  Severity = c("High", "High", "Critical", "Critical", "High", "High", "High", "High"),
  Resolution = c("Rebuilt reader-facing top-to-bottom notebook", "Replaced by nonsynonymous variant count",
                 "Exact Sample matching asserted before all models", "Within-study protein standardization; protected Step0 object retained",
                 "Deterministic top-variance selection", "Removed from publication workflow; transparent PCA-k-means sensitivity follows MOFA",
                 "Python MOFA exports factors, weights, view variance, model, and training statistics",
                 "Added a separate 25-sample expression-matched RPPA validation layer")
)
write_csv(core_audit, file.path(table_dir, "step4_legacy_static_audit_findings.csv"))

cat("Step4 core analysis and MOFA inputs complete\n")
print(assay_inclusion)
