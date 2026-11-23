#!/usr/bin/env Rscript

# Rebuild Step 0 processed objects from 0-Data/matrix into an isolated staging
# directory. Canonical files in 0-Data/qc are never overwritten by this script.

rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(123)

required <- c("data.table", "dplyr", "tidyr", "readr", "readxl", "limma", "impute", "sva", "tibble", "stringr")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Missing packages: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(tibble)
  library(stringr)
})

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args)) args[[1]] else file.path(tempdir(), "icca_step0_rebuild")
project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
matrix_dir <- file.path(project_root, "0-Data", "matrix")
qc_dir <- file.path(project_root, "0-Data", "qc")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(dir.exists(matrix_dir), dir.exists(qc_dir))

sample_info_path <- file.path(qc_dir, "sample_information.xlsx")
stopifnot(file.exists(sample_info_path))
sample_info <- as.data.frame(readxl::read_excel(sample_info_path))
stopifnot(nrow(sample_info) == 1092L, !anyDuplicated(sample_info$SampleID))

deduplicate_highest_mean <- function(data, gene_col) {
  expression_cols <- setdiff(names(data), gene_col)
  means <- rowMeans(data[, expression_cols, drop = FALSE], na.rm = FALSE)
  ordered <- data[order(means, decreasing = TRUE, na.last = TRUE), , drop = FALSE]
  ordered <- ordered[!duplicated(ordered[[gene_col]]), , drop = FALSE]
  rownames(ordered) <- ordered[[gene_col]]
  ordered[, setdiff(names(ordered), gene_col), drop = FALSE]
}

filter_expressed <- function(x) {
  x[is.na(x)] <- 0
  x[rowSums(x > 1) >= ceiling(ncol(x) * 0.20), , drop = FALSE]
}

fpkm_to_tpm <- function(x) {
  column_sums <- colSums(x)
  if (any(!is.finite(column_sums) | column_sums <= 0)) stop("Invalid FPKM column sum")
  sweep(x, 2, column_sums, "/") * 1e6
}

prepare_rnaseq <- function(study) {
  path <- file.path(matrix_dir, paste0(study, "_rnaseq.csv"))
  raw <- read.csv(path, check.names = TRUE)
  gene_col <- intersect(c("X...Gene", "Gene"), names(raw))[[1]]
  x <- deduplicate_highest_mean(raw, gene_col)

  if (study == "2023Deng") {
    colnames(x) <- gsub("\\.", "_", colnames(x))
    mapping <- sample_info[sample_info$Study == study & !is.na(sample_info$RNA), c("RNA", "SampleID")]
    keep <- colnames(x) %in% mapping$RNA
    x <- x[, keep, drop = FALSE]
    colnames(x) <- mapping$SampleID[match(colnames(x), mapping$RNA)]
  } else {
    colnames(x) <- sub("^X", "", colnames(x))
    if (study == "TCGA-CHOL") colnames(x) <- gsub("\\.", "-", colnames(x))
    target <- as.character(sample_info$SampleID[sample_info$Study == study])
    x <- x[, colnames(x) %in% target, drop = FALSE]
  }

  target_order <- as.character(sample_info$SampleID[sample_info$Study == study])
  if (!setequal(colnames(x), target_order)) stop("Sample mapping failed for ", study)
  x <- x[, target_order, drop = FALSE]
  x <- filter_expressed(as.matrix(x))
  if (study == "2023Deng") x <- fpkm_to_tpm(x)
  list(tpm = as.data.frame(x), log2 = as.data.frame(log2(x + 1)))
}

rna_studies <- c("2019Ahn", "2022Dong", "2023Deng", "2023Jeon", "2023Lee", "2024Fan", "TCGA-CHOL")
rna_results <- setNames(lapply(rna_studies, prepare_rnaseq), rna_studies)
tpm_2019Ahn <- rna_results[["2019Ahn"]]$tpm
tpm_2022Dong <- rna_results[["2022Dong"]]$tpm
tpm_2023Deng <- rna_results[["2023Deng"]]$tpm
tpm_2023Jeon <- rna_results[["2023Jeon"]]$tpm
tpm_2023Lee <- rna_results[["2023Lee"]]$tpm
tpm_2024Fan <- rna_results[["2024Fan"]]$tpm
tpm_tcga <- rna_results[["TCGA-CHOL"]]$tpm
log2tpm_2019Ahn <- rna_results[["2019Ahn"]]$log2
log2tpm_2022Dong <- rna_results[["2022Dong"]]$log2
log2tpm_2023Deng <- rna_results[["2023Deng"]]$log2
log2tpm_2023Jeon <- rna_results[["2023Jeon"]]$log2
log2tpm_2023Lee <- rna_results[["2023Lee"]]$log2
log2tpm_2024Fan <- rna_results[["2024Fan"]]$log2
log2tpm_tcga <- rna_results[["TCGA-CHOL"]]$log2
save(
  log2tpm_2019Ahn, log2tpm_2022Dong, log2tpm_2023Deng, log2tpm_2023Jeon,
  log2tpm_2023Lee, log2tpm_2024Fan, log2tpm_tcga, tpm_2019Ahn,
  tpm_2022Dong, tpm_2023Deng, tpm_2023Jeon, tpm_2023Lee, tpm_2024Fan,
  tpm_tcga, sample_info, file = file.path(output_dir, "0-Data-rnaseq.Rdata")
)

prepare_illumina <- function(filename, study) {
  raw <- read.csv(file.path(matrix_dir, filename), check.names = TRUE)
  x <- deduplicate_highest_mean(raw, "Gene")
  target <- as.character(sample_info$SampleID[sample_info$Study == study])
  x <- x[, colnames(x) %in% target, drop = FALSE]
  if (!all(colnames(x) %in% target)) stop("Illumina sample mapping failed: ", filename)
  normalized <- limma::normalizeBetweenArrays(log2(as.matrix(x)), method = "quantile")
  list(raw = x, log2 = normalized)
}

ill_2012 <- prepare_illumina("2012Andersen_rnaseq.csv", "2012Andersen")
ill_2013 <- prepare_illumina("2013Sia_rnaseq.csv", "2013Sia")
ill_2017_1 <- prepare_illumina("2017Jusakual_rnaseq1.csv", "2017Jusakual")
ill_2017_2 <- prepare_illumina("2017Jusakual_rnaseq2.csv", "2017Jusakual")
dat_2012Andersen <- ill_2012$raw
dat_2013Sia <- ill_2013$raw
dat_2017Jusakual1 <- ill_2017_1$raw
dat_2017Jusakual2 <- ill_2017_2$raw
log2dat_2012Andersen <- ill_2012$log2
log2dat_2013Sia <- ill_2013$log2
log2dat_2017Jusakual1 <- ill_2017_1$log2
log2dat_2017Jusakual2 <- ill_2017_2$log2
save(
  log2dat_2012Andersen, log2dat_2013Sia, log2dat_2017Jusakual1,
  log2dat_2017Jusakual2, dat_2012Andersen, dat_2013Sia,
  dat_2017Jusakual1, dat_2017Jusakual2, sample_info,
  file = file.path(output_dir, "0-Data-rnaseq-microarray1.Rdata")
)

prepare_affymetrix <- function(filename, study) {
  raw <- fread(file.path(matrix_dir, filename), check.names = FALSE, showProgress = FALSE)
  gene_col <- names(raw)[1]
  setnames(raw, gene_col, "Gene")
  target <- as.character(sample_info$SampleID[sample_info$Study == study])
  sample_cols <- intersect(names(raw), target)
  if (!setequal(sample_cols, target)) stop("Affymetrix sample mapping failed: ", study)
  collapsed <- raw[, lapply(.SD, median, na.rm = TRUE), by = Gene, .SDcols = sample_cols]
  collapsed <- collapsed[!grepl("---", Gene) & grepl("^[A-Z0-9]+$", Gene)]
  matrix_out <- as.data.frame(collapsed)
  rownames(matrix_out) <- matrix_out$Gene
  matrix_out$Gene <- NULL
  matrix_out[, target, drop = FALSE]
}

log2dat_2020Job <- prepare_affymetrix("2020Job_rnaseq.csv", "2020Job")
log2dat_TIGER <- prepare_affymetrix("TIGER-LC_rnaseq.csv", "TIGER-LC")
save(log2dat_2020Job, log2dat_TIGER, sample_info, file = file.path(output_dir, "0-Data-rnaseq-microarray2.Rdata"))

snv_tcga <- read.csv(file.path(matrix_dir, "TCGA-CHOL_dna.csv"))
snv_tcga <- snv_tcga[, intersect(names(snv_tcga), c(
  "Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Reference_Allele",
  "Tumor_Seq_Allele2", "Variant_Classification", "Variant_Type", "Tumor_Sample_Barcode"
))]
snv_tcga$Tumor_Sample_Barcode <- substr(snv_tcga$Tumor_Sample_Barcode, 1, 12)
snv_tcga$Dataset <- "TCGA"

snv_deng_wide <- read.csv(file.path(matrix_dir, "2023Deng_dna.csv"))
patient_col <- intersect(c("Patient.No.", "X...Patient.No."), names(snv_deng_wide))[[1]]
snv_deng_wide$.patient_id <- as.character(snv_deng_wide[[patient_col]])
gene_columns <- intersect(c(
  "ARID1A", "TP53", "KMT2C", "DPCR1", "BAP1", "MGA", "KMT2A", "KMT2D", "APC", "IDH1",
  "SMAD4", "ATM", "NCOR1", "EPHA2", "KRAS", "PML", "TET2", "NF1", "ARID2", "MSH6",
  "PBRM1", "RBM47", "FGFR2", "IDH2"
), names(snv_deng_wide))
snv_2023Deng <- snv_deng_wide |>
  pivot_longer(all_of(gene_columns), names_to = "Hugo_Symbol", values_to = "Variant_Classification") |>
  filter(!is.na(Variant_Classification)) |>
  transmute(
    Hugo_Symbol, Chromosome=NA, Start_Position=NA, End_Position=NA,
    Reference_Allele=NA, Tumor_Seq_Allele2=NA, Variant_Classification,
    Variant_Type=case_when(
      grepl("Missense|Nonsense|Splice_Site|Nonstop", Variant_Classification, ignore.case=TRUE) ~ "SNP",
      grepl("Del|Deletion", Variant_Classification, ignore.case=TRUE) ~ "DEL",
      grepl("Ins|Insertion", Variant_Classification, ignore.case=TRUE) ~ "INS",
      grepl("Fusion", Variant_Classification, ignore.case=TRUE) ~ "Fusion", TRUE ~ "COMPLEX"
    ),
    Tumor_Sample_Barcode=.patient_id, Dataset="Deng2023"
  )

snv_2022Dong <- read.csv(file.path(matrix_dir, "2022Dong_dna.csv")) |>
  transmute(
    Hugo_Symbol=Gene, Chromosome=CHROM, Start_Position=POS, End_Position=POS,
    Reference_Allele=REF, Tumor_Seq_Allele2=ALT, Variant_Classification=Mutation_Type,
    Variant_Type=case_when(
      nchar(REF)==1 & nchar(ALT)==1 ~ "SNP", nchar(REF)>nchar(ALT) ~ "DEL",
      nchar(REF)<nchar(ALT) ~ "INS", TRUE ~ "COMPLEX"
    ),
    Tumor_Sample_Barcode=as.character(Sample_ID), Dataset="Dong2022"
  )
combined_maf_data <- bind_rows(snv_tcga, snv_2023Deng, snv_2022Dong)
filter_out <- c("intron_variant", "Intron", "3_prime_UTR_variant", "3'UTR", "3'Flank", "5_prime_UTR_variant", "5'UTR", "5'Flank", "intergenic_variant", "RNA", "upstream_gene_variant", "downstream_gene_variant", "sequence_feature")
combined_maf_data <- combined_maf_data |>
  filter(!is.na(Variant_Classification), !Variant_Classification %in% filter_out) |>
  mutate(Variant_Classification=case_when(
    Variant_Classification %in% c("missense_variant", "Missense_Mutation") ~ "Missense_Mutation",
    Variant_Classification %in% c("stop_gained", "Nonsense_Mutation") ~ "Nonsense_Mutation",
    Variant_Classification %in% c("stop_lost", "Nonstop_Mutation") ~ "Nonstop_Mutation",
    Variant_Classification %in% c("start_lost", "Translation_Start_Site") ~ "Translation_Start_Site",
    Variant_Classification %in% c("splice_region_variant", "splice_donor_variant", "splice_acceptor_variant", "Splice_Site") ~ "Splice_Site",
    Variant_Classification %in% c("disruptive_inframe_deletion", "conservative_inframe_deletion", "inframe_deletion") ~ "In_Frame_Del",
    Variant_Classification %in% c("disruptive_inframe_insertion", "conservative_inframe_insertion", "inframe_insertion") ~ "In_Frame_Ins",
    Variant_Classification == "synonymous_variant" ~ "Silent",
    Variant_Classification == "frameshift_variant" & nchar(Reference_Allele)>nchar(Tumor_Seq_Allele2) ~ "Frame_Shift_Del",
    Variant_Classification == "frameshift_variant" & nchar(Reference_Allele)<nchar(Tumor_Seq_Allele2) ~ "Frame_Shift_Ins",
    TRUE ~ as.character(Variant_Classification)
  ))
save(combined_maf_data, file=file.path(output_dir, "0-Data-dnaseq.Rdata"))
write_csv(combined_maf_data, file.path(output_dir, "0-Data-combined_maf_data.csv"))

pro_2022Dong <- read.csv(file.path(matrix_dir, "2022Dong_protein.csv"))
dong_long <- pro_2022Dong |>
  rename(Protein=Gene) |>
  pivot_longer(-Protein, names_to="Sample", values_to="Abundance") |>
  mutate(Sample=str_remove(Sample, "^X"), Batch="2022Dong")
pro_2023Deng <- read.csv(file.path(matrix_dir, "2023Deng_protein.csv"))
deng_long <- pro_2023Deng |>
  rename(Protein=Tumor..T..proteomic.ID) |>
  pivot_longer(-Protein, names_to="Sample", values_to="Abundance") |>
  mutate(Sample=gsub("^X", "", Sample)) |>
  left_join(sample_info |> select(`Protein profiling`, SampleID), by=c("Sample"="Protein profiling")) |>
  transmute(Protein, Sample=SampleID, Abundance, Batch="2023Deng") |>
  filter(!is.na(Sample))
pro_tcga <- read.csv(file.path(matrix_dir, "TCGA-CHOL_protein.csv"))
all_prot_long <- bind_rows(dong_long, deng_long) |>
  filter(!is.na(Abundance), !is.na(Protein)) |>
  group_by(Protein, Sample, Batch) |>
  summarise(Abundance=mean(Abundance, na.rm=TRUE), .groups="drop")
prot_matrix <- all_prot_long |>
  select(Protein, Sample, Abundance) |>
  group_by(Protein, Sample) |>
  summarise(Abundance=mean(Abundance, na.rm=TRUE), .groups="drop") |>
  pivot_wider(names_from=Sample, values_from=Abundance) |>
  column_to_rownames("Protein") |>
  as.matrix()
prot_matrix <- prot_matrix[rowSums(is.na(prot_matrix)) <= ncol(prot_matrix)*0.70, , drop=FALSE]
prot_matrix_imputed <- impute::impute.knn(prot_matrix)$data
protein_batch_info <- all_prot_long |> distinct(Sample, Batch)
protein_batch_info <- protein_batch_info[match(colnames(prot_matrix_imputed), protein_batch_info$Sample), ]
stopifnot(identical(as.character(protein_batch_info$Sample), colnames(prot_matrix_imputed)))
stopifnot(setequal(unique(as.character(protein_batch_info$Batch)), c("2022Dong", "2023Deng")))
corrected_prot_matrix <- sva::ComBat(dat=prot_matrix_imputed, batch=protein_batch_info$Batch)
save(pro_2022Dong, pro_2023Deng, pro_tcga, corrected_prot_matrix, file=file.path(output_dir, "0-Data-protein.Rdata"))
write.csv(corrected_prot_matrix, file.path(output_dir, "0-Data-corrected-protein.csv"), row.names=TRUE)

writeLines(capture.output(sessionInfo()), file.path(output_dir, "step0_rebuild_session_info.txt"))
message("Step 0 rebuild completed in: ", normalizePath(output_dir, winslash="/", mustWork=TRUE))
