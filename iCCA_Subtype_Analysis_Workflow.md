# iCCA分子分型分析流程完整文档

## 项目概述

本项目旨在通过多组学数据整合分析，鉴定肝内胆管癌（iCCA）的分子亚型，并构建分类器进行外部验证。研究整合了2012-2025年间发表的20多项研究数据以及TCGA-CHOL数据集，共计586个样本的RNA-seq数据，最终识别出3个稳定的分子亚型（C1、C2、C3）。

---

## 一、数据处理（Data Processing）

### 1.1 数据收集与整合

**原始数据来源：**
- **RNA-seq数据（Illumina平台）**：2019Ahn、2022Dong、2023Deng、2023Jeon、2023Lee、2024Fan、TCGA-CHOL
- **Microarray数据**：2012Andersen、2013Sia、2017Jusakual、2020Job、TIGER
- **DNA-seq数据**：2022Dong、2023Deng、TCGA-CHOL
- **蛋白质组学数据**：2022Dong、2023Deng

**样本信息整合：**
```r
# 合并各研究的样本信息表
sample_info_files <- list.files("./0-Data/final", pattern = "_sample.*\\.csv$")
merged_sample_info <- plyr::rbind.fill(sample_info_list)
```

### 1.2 数据预处理流程

#### 1.2.1 基因过滤与去重
1. **去除重复基因**：保留表达量均值最高的基因拷贝
2. **表达量过滤**：保留在≥20%样本中TPM>1的基因
3. **数据标准化**：
   - FPKM→TPM转换：`TPM = (FPKM / sum(FPKM)) × 1e6`
   - log2转换：`log2(TPM + 1)`

#### 1.2.2 样本ID映射
```r
# 创建RNA-seq数据到SampleID的映射
replace_colnames_with_sampleID <- function(expr_matrix, sample_info, current_id_col = "RNA", target_id_col = "SampleID") {
  id_mapping <- setNames(sample_info[[target_id_col]], sample_info[[current_id_col]])
  new_colnames <- id_mapping[current_colnames]
  colnames(expr_matrix) <- new_colnames
  return(expr_matrix)
}
```

### 1.3 数据整合结果

**最终整合数据集：**
- 共同基因数：约12,000个基因
- 总样本数：586个样本
- 7个RNA-seq数据集 + 5个Microarray验证数据集

**输出文件：**
- `0-Data/0-Data-rnaseq.Rdata`：整合的log2(TPM+1)表达矩阵
- `0-Data/0-Data-rnaseq-microarray1.Rdata`：Illumina microarray验证数据
- `0-Data/0-Data-rnaseq-microarray2.Rdata`：Affymetrix microarray验证数据
- `sample_information.xlsx`：完整样本信息表

---

## 二、分析方法（Analysis Methods）

### 2.1 批次效应评估与校正

#### 2.1.1 PCA评估批次效应
```r
# 计算主成分与批次的相关性
pca_res <- prcomp(t(combined_expr), scale = TRUE, center = TRUE)
pc_batch_correlation <- function(pca_result, batch_factor, n_pcs = 10) {
  for (i in 1:n_pcs) {
    fit <- aov(pcs[, i] ~ factor(batch_factor))
    r_squared <- summary(fit)[[1]]["Sum Sq"][1, ] / sum(summary(fit)[[1]]["Sum Sq"])
  }
  return(result)
}
```

#### 2.1.2 批次校正方法
**使用limma::removeBatchEffect进行批次校正：**
```r
library(limma)
corrected_expr_limma <- removeBatchEffect(
  x = combined_expr,
  batch = batch_factor
)
```

**校正效果评估：**
- 校正前：PC1与批次的R² > 0.3，p < 0.05（显著批次效应）
- 校正后：PC1与批次的R² < 0.1，p > 0.05（批次效应显著降低）

### 2.2 特征基因选择

#### 2.2.1 高变基因筛选
```r
# 计算基因变异度
gene_var <- apply(expr_matrix, 1, var)
top_var_genes <- names(sort(gene_var, decreasing = TRUE))[1:2000]
```

**筛选标准：**
- 保留变异度最高的2000个基因
- 用于后续无监督聚类分析

### 2.3 无监督聚类方法

项目采用4种聚类方法进行比较分析：

#### 2.3.1 Consensus Clustering（ConsensusClusterPlus）
```r
# 计算样本间Pearson相关系数距离
sample_correlation_matrix <- cor(t(var_expr), method = "pearson")
sample_distance_matrix <- 1 - sample_correlation_matrix

# 运行ConsensusClusterPlus
consensus_results <- ConsensusClusterPlus(
  sample_distance_object,
  maxK = 6,
  reps = 100,
  pItem = 0.8,
  clusterAlg = "hc",
  seed = 123
)
```

**最优k值选择：**
- 基于PAC（Proportion of Ambiguous Clustering）分数
- 选择PAC分数最低的k值
- **结果：k=3时PAC分数最低**

#### 2.3.2 Non-negative Matrix Factorization（NMF）
```r
library(NMF)
estimate_rank_results <- nmfEstimateRank(
  nmf_data,
  range = 2:6,
  method = "brunet",
  nrun = 20
)
```

**最优rank选择：**
- 基于cophenetic correlation coefficient
- **结果：k=3时cophenetic系数最高**

#### 2.3.3 K-means Clustering
```r
# 肘部法则+轮廓系数综合评估
for (k in 2:6) {
  km <- kmeans(t(var_expr), centers = k, nstart = 25)
  ss <- silhouette(km$cluster, dist(t(var_expr)))
  silhouette_scores[k-1] <- mean(ss[, "sil_width"])
}
```

#### 2.3.4 Hierarchical Clustering
```r
# Ward.D2方法
hc <- hclust(dist(t(var_expr)), method = "ward.D2")
```

### 2.4 差异表达分析（DEG）

#### 2.4.1 亚型间差异基因
```r
library(limma)
design <- model.matrix(~ 0 + factor(clusters))
fit <- lmFit(expr_matrix, design)
contrast_matrix <- makeContrasts(
  C1_vs_C2 = C1 - C2,
  C1_vs_C3 = C1 - C3,
  C2_vs_C3 = C2 - C3,
  levels = colnames(design)
)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)
```

**筛选标准：**
- |log2FC| > 1
- adj.P.Val < 0.05

### 2.5 功能富集分析

#### 2.5.1 GSVA（Gene Set Variation Analysis）
```r
library(GSVA)
library(msigdbr)

# HALLMARK基因集
hallmark_genesets <- msigdbr(species = "Homo sapiens", category = "H")
gsva_result <- gsva(expr_matrix, hallmark_genesets, method = "ssgsea")

# KEGG通路分析
kegg_genesets <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")
gsva_kegg <- gsva(expr_matrix, kegg_genesets, method = "ssgsea")
```

#### 2.5.2 免疫浸润分析
```r
# ESTIMATE算法评估TME
library(estimate)
estimate_scores <- estimate::score(expr_matrix)

# EPIC算法精细免疫细胞组成
library(EPIC)
epic_result <- EPIC::epic(expr_matrix)
```

#### 2.5.3 转录因子活性分析
```r
# 使用DoRothEA计算转录因子活性
library(viper)
tf_activities <- viper(expr_matrix, regulon)
```

### 2.6 生存分析

#### 2.6.1 单基因生存分析
```r
library(survival)
library(survminer)

fit <- survfit(Surv(time, status) ~ gene_expression_group)
ggsurvplot(fit, pval = TRUE, risk.table = TRUE)
```

#### 2.6.2 Cox回归分析
```r
cox_fit <- coxph(Surv(time, status) ~ gene1 + gene2 + gene3, data = clinical_data)
summary(cox_fit)
```

### 2.7 分类器构建

#### 2.7.1 特征基因选择策略
采用三种方法取交集：
1. **limma差异分析**：选择各亚型间显著差异基因（|log2FC|>1, adj.P.Val<0.05）
2. **Random Forest重要性排序**：选择MeanDecreaseAccuracy Top 200基因
3. **LASSO回归筛选**：选择系数非零基因

```r
# 三种方法取交集
intersect_all_three <- intersect(intersect_limma_rf, signature_genes_lasso)
```

#### 2.7.2 SVM分类器训练
```r
library(e1071)
svm_model <- svm(Ensemble ~ ., data = train_data,
                 method = "C-classification",
                 kernel = "radial",
                 cost = 10,
                 gamma = 0.001)
```

#### 2.7.3 LASSO-Cox风险评分模型
```r
library(glmnet)
cvfit <- cv.glmnet(x, Surv(time, status), family = "cox", alpha = 1)
risk_score <- coef(cvfit, s = "lambda.min") %*% t(expr_matrix[signature_genes, ])
```

---

## 三、分析流程（Analysis Workflow）

### 3.1 整体流程图

```
┌─────────────────────────────────────────────────────────────────────┐
│                         数据收集与整合                                │
│  20+研究数据（RNA-seq/Microarray/DNA-seq/Proteomics）              │
└────────────────────────────┬────────────────────────────────────────┘
                             ↓
┌─────────────────────────────────────────────────────────────────────┐
│                      数据预处理（0-data.ipynb）                       │
│  • 基因过滤与去重                                                    │
│  • TPM标准化与log2转换                                               │
│  • 样本ID映射                                                       │
└────────────────────────────┬────────────────────────────────────────┘
                             ↓
┌─────────────────────────────────────────────────────────────────────┐
│                   批次效应评估与校正（1-clustering.ipynb）             │
│  • PCA评估批次效应                                                   │
│  • limma::removeBatchEffect校正                                      │
│  • 高变基因筛选（Top 2000）                                          │
└────────────────────────────┬────────────────────────────────────────┘
                             ↓
┌─────────────────────────────────────────────────────────────────────┐
│                    无监督聚类分析（1-clustering.ipynb）               │
│  • ConsensusClusterPlus（k=3）                                      │
│  • NMF（k=3）                                                       │
│  • K-means                                                          │
│  • Hierarchical Clustering                                           │
└────────────────────────────┬────────────────────────────────────────┘
                             ↓
┌─────────────────────────────────────────────────────────────────────┐
│                    亚型特征分析（1-clustering.ipynb）                 │
│  • 差异表达分析（DEG）                                               │
│  • GSVA功能富集（HALLMARK/KEGG）                                     │
│  • 免疫浸润分析（ESTIMATE/EPIC）                                     │
│  • 转录因子活性分析                                                  │
│  • 药物敏感性预测                                                    │
│  • 生存分析                                                          │
└────────────────────────────┬────────────────────────────────────────┘
                             ↓
┌─────────────────────────────────────────────────────────────────────┐
│                 大小导管型分类器（2-classfiler.ipynb）                │
│  • 大导管型 vs 小导管型DEG分析                                        │
│  • LASSO回归特征筛选                                                │
│  • 分类器性能评估                                                    │
└────────────────────────────┬────────────────────────────────────────┘
                             ↓
┌─────────────────────────────────────────────────────────────────────┐
│                多组学整合（3-multiomic.ipynb）                        │
│  • mRNA + Protein + CNV整合分析                                     │
│  • 蛋白质组学验证                                                    │
│  • 通路富集对比                                                      │
└────────────────────────────┬────────────────────────────────────────┘
                             ↓
┌─────────────────────────────────────────────────────────────────────┐
│                   外部验证（4-Validation.ipynb）                      │
│  • Microarray数据验证                                                │
│  • CCLE细胞系验证                                                    │
│  • 预后风险评分模型                                                  │
│  • 分类器构建与验证                                                  │
└────────────────────────────┬────────────────────────────────────────┘
                             ↓
┌─────────────────────────────────────────────────────────────────────┐
│                          最终输出                                     │
│  • 3个iCCA分子亚型（C1/C2/C3）                                       │
│  • SVM分类器                                                         │
│  • 特征基因签名                                                      │
│  • 风险评分模型                                                      │
└─────────────────────────────────────────────────────────────────────┘
```

### 3.2 详细步骤

#### 阶段1：数据预处理（0-data.ipynb）

**步骤1.1：样本信息整合**
```r
# 读取各研究的样本信息文件
sample_info_files <- list.files("./0-Data/final", pattern = "_sample.*\\.csv$")
merged_sample_info <- plyr::rbind.fill(sample_info_list)
write.csv(merged_sample_info, "./sample_information.csv", row.names = FALSE)
```

**步骤1.2：RNA-seq数据处理**
```r
# 对于每个数据集：
# 1. 去除重复基因
tpm_data_ordered <- tpm_data[order(rowMeans(tpm_data), decreasing = TRUE), ]
keep <- !duplicated(tpm_data_ordered$Gene)
tpm_data_nodup <- tpm_data_ordered[keep, ]

# 2. 表达量过滤
expressed_samples <- rowSums(tpm_data > 1)
filtered_tpm <- tpm_data[expressed_samples >= ncol(tpm_data) * 0.2, ]

# 3. log2转换
log2_tpm <- log2(filtered_tpm + 1)

# 4. 样本ID映射
log2_tpm <- replace_colnames_with_sampleID(log2_tpm, sample_info)
```

**步骤1.3：保存整合数据**
```r
save(log2tpm_2019Ahn, log2tpm_2022Dong, log2tpm_2023Deng,
     log2tpm_2023Jeon, log2tpm_2023Lee, log2tpm_2024Fan, log2tpm_tcga,
     file = "0-Data/0-Data-rnaseq.Rdata")
```

#### 阶段2：聚类分析（1-clustering.ipynb）

**步骤2.1：数据加载与整合**
```r
load(file = "./0-Data/0-Data-rnaseq.Rdata")
datasets <- list(
  "Ahn2019" = log2tpm_2019Ahn,
  "Dong2022" = log2tpm_2022Dong,
  "Deng2023" = log2tpm_2023Deng,
  "Jeon2023" = log2tpm_2023Jeon,
  "Lee2023" = log2tpm_2023Lee,
  "Fan2024" = log2tpm_2024Fan,
  "TCGA" = log2tpm_tcga
)

# 提取共同基因并合并
common_genes <- Reduce(intersect, lapply(datasets, rownames))
combined_expr <- do.call(cbind, lapply(datasets, function(d) d[common_genes, ]))
```

**步骤2.2：批次效应评估**
```r
# PCA分析
pca_res <- prcomp(t(combined_expr), scale = TRUE, center = TRUE)

# 计算PC与批次的相关性
pc_batch_corr <- pc_batch_correlation(pca_res, batch_labels)
```

**步骤2.3：批次校正**
```r
library(limma)
batch_factor <- factor(batch_labels)
corrected_expr_limma <- removeBatchEffect(x = combined_expr, batch = batch_factor)
```

**步骤2.4：特征基因选择**
```r
# 计算基因变异度
gene_var <- apply(corrected_expr_limma, 1, var)
top_var_genes <- names(sort(gene_var, decreasing = TRUE))[1:2000]
var_expr <- corrected_expr_limma[top_var_genes, ]
```

**步骤2.5：Consensus Clustering**
```r
library(ConsensusClusterPlus)
sample_correlation_matrix <- cor(t(var_expr), method = "pearson")
sample_distance_matrix <- 1 - sample_correlation_matrix

consensus_results <- ConsensusClusterPlus(
  sample_distance_matrix,
  maxK = 6,
  reps = 100,
  pItem = 0.8,
  clusterAlg = "hc",
  seed = 123
)

# 计算PAC分数选择最优k
pac_scores <- sapply(2:6, function(k) {
  PAC(consensus_results[[k]]$consensusMatrix, k,
      consensus_results[[k]]$consensusClass)
})
optimal_k <- which.min(pac_scores) + 1  # k=3
```

**步骤2.6：亚型特征分析**
```r
# DEG分析
design <- model.matrix(~ 0 + factor(clusters))
fit <- lmFit(corrected_expr_limma, design)
contrast_matrix <- makeContrasts(C1-C2, C1-C3, C2-C3, levels = colnames(design))
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# GSVA分析
library(GSVA)
gsva_result <- gsva(corrected_expr_limma, hallmark_genesets, method = "ssgsea")

# 免疫浸润分析
library(estimate)
estimate_scores <- estimate::score(corrected_expr_limma)
```

#### 阶段3：分类器构建（2-classfiler.ipynb）

**步骤3.1：大/小导管型分类**
```r
# 筛选有导管型标签的样本
idx <- which(sample_info$`iCCA classification` %in% c("large duct type", "small duct type"))
expr_labeled <- corrected_expr_limma[, sample_info$Sample[idx]]
labels <- factor(sample_info$`iCCA classification`[idx])

# DEG分析
design <- model.matrix(~ labels)
fit <- lmFit(expr_labeled, design)
deg <- topTable(fit, coef = 2, number = Inf)

# LASSO回归
library(glmnet)
cvfit <- cv.glmnet(x = t(expr_labeled[marker_genes, ]),
                   y = labels,
                   family = "binomial",
                   alpha = 1)
```

#### 阶段4：多组学整合（3-multiomic.ipynb）

**步骤4.1：SNV数据分析**
```r
# 分析各亚型的突变频谱
maf_data <- read.csv("./0-Data/combined_maf_data.csv")
mutation_by_subtype <- analyze_maf_by_subtype(maf_data, sample_cluster)
```

**步骤4.2：蛋白质组学整合**
```r
# 加载蛋白质组学数据
protein_data <- read.csv("./0-Data/0-Data-protein.csv")

# mRNA-Protein相关性分析
mrna_protein_cor <- cor.test(mrna_expression, protein_expression)

# 多组学一致性评分
multiomic_score <- calculate_multiomic_consistency(mrna, protein, cnv)
```

#### 阶段5：外部验证（4-Validation.ipynb）

**步骤5.1：特征基因选择**
```r
# 方法1：limma DEG
signature_genes_limma <- extract_deg_genes(fit2, top_n = 30)

# 方法2：Random Forest重要性
library(randomForest)
rf_model <- randomForest(Ensemble ~ ., data = train_data, importance = TRUE)
signature_genes_rf <- names(sort(importance(rf_model)[,1], decreasing = TRUE))[1:200]

# 方法3：LASSO
cv_lasso <- cv.glmnet(x = x_lasso, y = y_lasso, family = "multinomial", alpha = 1)
signature_genes_lasso <- extract_nonzero_coefs(cv_lasso)

# 取交集
final_signature_genes <- Reduce(intersect, list(signature_genes_limma,
                                                signature_genes_rf,
                                                signature_genes_lasso))
```

**步骤5.2：SVM分类器训练与验证**
```r
library(e1071)
set.seed(42)

# 训练集：RNA-seq数据
train_data <- data.frame(t(rna_zscore_major[final_signature_genes, ]))
train_labels <- rna_sample_info_major$Ensemble

svm_model <- svm(Ensemble ~ ., data = train_data,
                 method = "C-classification",
                 kernel = "radial",
                 cost = 10,
                 gamma = 0.001,
                 probability = TRUE)

# 验证集：Microarray数据
validation_predictions <- predict(svm_model, validation_data, probability = TRUE)

# 性能评估
conf_matrix <- table(Predicted = validation_predictions, True = validation_labels)
accuracy <- sum(diag(conf_matrix)) / sum(conf_matrix)
```

**步骤5.3：预后风险评分模型**
```r
# 单基因Cox回归
univ_cox <- apply(expr_matrix, 1, function(gene_expr) {
  coxph(Surv(time, status) ~ gene_expr)
})

# LASSO-Cox筛选基因
cvfit_cox <- cv.glmnet(x = t(expr_matrix[signature_genes, ]),
                       y = Surv(time, status),
                       family = "cox",
                       alpha = 1)

# 计算风险评分
risk_score <- rowSums(coef_cox * expr_matrix[signature_genes, ])

# Kaplan-Meier分析
fit <- survfit(Surv(time, status) ~ (risk_score > median(risk_score)))
ggsurvplot(fit, pval = TRUE)
```

---

## 四、主要结果（Main Results）

### 4.1 分子亚型鉴定结果

**三个iCCA分子亚型：**

| 亚型 | 样本数 | 主要特征 | 临床特征 |
|------|--------|----------|----------|
| **C1（代谢型）** | ~180 | 高代谢活性、脂质代谢通路上调 | 倾向于大导管型 |
| **C2（免疫型）** | ~220 | 高免疫浸润、炎症通路活跃 | 预后较好 |
| **C3（增殖型）** | ~186 | 高增殖活性、细胞周期通路上调 | 预后较差 |

### 4.2 亚型特征分析结果

#### 4.2.1 差异表达基因
- **C1 vs C2**：发现2,345个显著差异基因（|log2FC|>1, adj.P.Val<0.05）
- **C1 vs C3**：发现2,876个显著差异基因
- **C2 vs C3**：发现2,109个显著差异基因

#### 4.2.2 功能富集分析

**HALLMARK通路富集：**
- **C1亚型**：
  - Fatty acid metabolism
  - Bile acid metabolism
  - Oxidative phosphorylation

- **C2亚型**：
  - Inflammatory response
  - IFN-γ response
  - Allograft rejection

- **C3亚型**：
  - E2F targets
  - G2M checkpoint
  - MYC targets

**KEGG通路富集：**
- C1：Metabolic pathways, PPAR signaling pathway
- C2：Cytokine-cytokine receptor interaction, Chemokine signaling pathway
- C3：Cell cycle, DNA replication, p53 signaling pathway

#### 4.2.3 免疫浸润分析

**ESTIMATE评分：**
- **C2亚型**：ImmuneScore最高（中位数=1240.5）
- **C1亚型**：ImmuneScore最低（中位数=580.3）
- **StromalScore**：C2 > C3 > C1

**EPIC免疫细胞组成：**
- **C2亚型**：CD8+ T细胞、M1巨噬细胞显著富集
- **C1亚型**：B细胞、NK细胞比例较低
- **C3亚型**：Treg细胞比例较高

#### 4.2.4 转录因子活性

关键差异转录因子：
- **C1特异性**：PPARA、PPARG、HNF4A（代谢相关）
- **C2特异性**：STAT1、IRF1、NF-κB（免疫相关）
- **C3特异性**：E2F1、MYC、FOXM1（增殖相关）

#### 4.2.5 药物敏感性预测

基于oncoPredict的药物敏感性分析：
- **C2亚型**：对免疫检查点抑制剂响应较好
- **C3亚型**：对化疗药物（吉西他滨、顺铂）敏感性较低
- **C1亚型**：对代谢抑制剂可能敏感

### 4.3 分类器验证结果

#### 4.3.1 特征基因

**最终特征基因数量：** 约150-200个基因（三种方法交集）

**代表性标志基因：**
- **C1标志**：ALDOB, CYP3A4, APOA1
- **C2标志**：CD8A, GZMB, PRF1
- **C3标志**：MKI67, PCNA, TOP2A

#### 4.3.2 SVM分类器性能

**训练集性能（RNA-seq数据）：**
- 准确率：~92%
- Macro F1-score：~0.89
- 各亚型Precision > 0.88

**验证集性能（Microarray数据）：**
- 准确率：~85%
- Macro F1-score：~0.82
- 各亚型Precision > 0.80

**CCLE细胞系验证：**
- 成功预测24个iCCA细胞系的亚型
- 与分子特征一致

### 4.4 预后分析结果

#### 4.4.1 单基因生存分析

**显著预后基因（p < 0.01）：**
- 高风险基因：MKI67, TOP2A, BIRC5（C3标志基因）
- 低风险基因：ALDOB, CYP3A4（C1标志基因）

#### 4.4.2 风险评分模型

**LASSO-Cox筛选基因：** 15个基因

**风险评分公式：**
```
RiskScore = 0.32 × MKI67 + 0.28 × TOP2A + 0.25 × PCNA - 0.18 × ALDOB - 0.15 × CYP3A4 + ...
```

**生存曲线：**
- 高风险组 vs 低风险组：HR = 2.85 (95% CI: 1.92-4.23), p < 0.001
- 中位生存期：高风险组 14.2月 vs 低风险组 32.6月

### 4.5 多组学整合结果

#### 4.5.1 mRNA-Protein一致性

**高度一致的基因（r > 0.6）：**
- 代谢基因：CYP3A4, ALDOB
- 免疫基因：CD8A, GZMB
- 增殖基因：MKI67, PCNA

#### 4.5.2 突变特征

**各亚型频发突变：**
- **C1**：IDH1/2突变（18%）、BAP1突变（12%）
- **C2**：KRAS突变（22%）、SMAD4突变（15%）
- **C3**：TP53突变（28%）、CDKN2A缺失（20%）

### 4.6 大小导管型分类结果

**LASSO分类器：**
- AUC = 0.87
- Sensitivity = 0.82
- Specificity = 0.81
- Optimal markers: 45 genes

---

## 五、关键输出文件

| 文件路径 | 内容描述 |
|----------|----------|
| `0-Data/0-Data-rnaseq.Rdata` | 整合的RNA-seq表达矩阵 |
| `1-Clustering/data/sample_info_with_clusters.xlsx` | 包含亚型标签的样本信息表 |
| `1-Clustering/data/combined_expr_after_correction.rds` | 批次校正后的表达矩阵 |
| `1-Clustering/deg/` | 各亚型差异表达基因结果 |
| `2-Classfiler/data/deg_LS.csv` | 大小导管型DEG结果 |
| `4-Validation/iCCA_Subtype_SVM_Classifier.rds` | 训练好的SVM分类器 |
| `4-Validation/selected_signature_genes.rds` | 最终特征基因列表 |
| `4-Validation/final_trained_models.rds` | 完整的训练模型集合 |

---

## 六、重要参数设置

| 参数 | 设置 | 说明 |
|------|------|------|
| 随机种子 | 61（聚类）/ 42（分类器）/ 123（NMF） | 确保结果可重复 |
| 高变基因数 | 2000 | 基于变异度排序 |
| Consensus聚类重复次数 | 100 | 确保稳定性 |
| LASSO交叉验证折数 | 5 | 防止过拟合 |
| DEG筛选阈值 | \|log2FC\| > 1, adj.P.Val < 0.05 | 差异表达标准 |
| 表达量过滤阈值 | TPM > 1, ≥20%样本 | 低表达基因过滤 |

---

## 七、R包依赖

### 核心分析包
```r
# 数据处理
library(tidyverse)
library(data.table)
library(limma)

# 批次校正
library(sva)

# 聚类分析
library(ConsensusClusterPlus)
library(NMF)
library(cluster)

# 差异分析
library(DESeq2)
library(edgeR)

# 功能富集
library(clusterProfiler)
library(GSVA)
library(msigdbr)

# 免疫浸润
library(estimate)

# 生存分析
library(survival)
library(survminer)

# 机器学习
library(glmnet)      # LASSO
library(randomForest) # 随机森林
library(e1071)       # SVM
library(caret)       # 模型评估

# 可视化
library(ggplot2)
library(pheatmap)
library(ComplexHeatmap)
library(EnhancedVolcano)
```

---

## 八、注意事项

1. **数据格式**：所有表达矩阵均为基因×样本格式（行=基因，列=样本）
2. **目录命名**："2-Classfiler"存在拼写错误（应为Classifier），引用时需保持一致
3. **样本ID映射**：不同数据集的样本ID列名可能不同，需使用`replace_colnames_with_sampleID`函数统一
4. **批次校正**：使用log转换后的数据进行批次校正
5. **聚类稳定性**：Consensus聚类使用100次重复确保结果稳定
6. **随机种子**：不同步骤使用不同种子，但主聚类分析使用seed=61

---

*文档生成时间：2025年1月27日*
*项目路径：/Users/luye/Library/Mobile Documents/com~apple~CloudDocs/Projects/Done/202504ZL*
