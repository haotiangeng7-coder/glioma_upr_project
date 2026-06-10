###############################################################################
# download_tcga.R
# 下载TCGA GBM + LGG数据（RNA-seq, 临床, 突变, CNV）
# 使用TCGAbiolinks从GDC Portal下载
###############################################################################

source("00_setup/config.R")

library(TCGAbiolinks)
library(SummarizedExperiment)
library(dplyr)

tcga_dir <- file.path(DATA_RAW, "TCGA")
dir.create(tcga_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. 下载RNA-seq表达数据 (TCGA-GBM)
# =============================================================================
message("=== Downloading TCGA-GBM RNA-seq ===")

query_gbm <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(query_gbm, directory = file.path(tcga_dir, "GDCdata"))
data_gbm <- GDCprepare(query_gbm, directory = file.path(tcga_dir, "GDCdata"))

# =============================================================================
# 2. 下载RNA-seq表达数据 (TCGA-LGG)
# =============================================================================
message("=== Downloading TCGA-LGG RNA-seq ===")

query_lgg <- GDCquery(
  project = "TCGA-LGG",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)
GDCdownload(query_lgg, directory = file.path(tcga_dir, "GDCdata"))
data_lgg <- GDCprepare(query_lgg, directory = file.path(tcga_dir, "GDCdata"))

# =============================================================================
# 3. 合并GBM + LGG数据
# =============================================================================
message("=== Merging GBM and LGG data ===")

# 取共同基因
common_genes <- intersect(rownames(data_gbm), rownames(data_lgg))
data_gbm_sub <- data_gbm[common_genes, ]
data_lgg_sub <- data_lgg[common_genes, ]

# 合并 — 对齐colData列名（GBM和LGG的临床列可能不同）
common_cols <- intersect(colnames(colData(data_gbm_sub)), colnames(colData(data_lgg_sub)))
colData(data_gbm_sub) <- colData(data_gbm_sub)[, common_cols]
colData(data_lgg_sub) <- colData(data_lgg_sub)[, common_cols]
tcga_glioma <- cbind(data_gbm_sub, data_lgg_sub)

# 提取表达矩阵（TPM）
expr_counts <- assay(tcga_glioma, "unstranded")
expr_tpm    <- assay(tcga_glioma, "tpm_unstrand")
expr_fpkm   <- assay(tcga_glioma, "fpkm_unstrand")

# 提取临床数据
clinical <- as.data.frame(colData(tcga_glioma))

# 基因注释
gene_info <- as.data.frame(rowData(tcga_glioma))

# 转换行名为gene symbol（去重）
gene_symbols <- gene_info$gene_name
dup_idx <- duplicated(gene_symbols) | is.na(gene_symbols) | gene_symbols == ""
expr_tpm_symbol <- expr_tpm[!dup_idx, ]
rownames(expr_tpm_symbol) <- gene_symbols[!dup_idx]

expr_counts_symbol <- expr_counts[!dup_idx, ]
rownames(expr_counts_symbol) <- gene_symbols[!dup_idx]

message(sprintf("TCGA glioma: %d genes x %d samples", nrow(expr_tpm_symbol), ncol(expr_tpm_symbol)))

# =============================================================================
# 4. 整理临床数据
# =============================================================================

clinical_clean <- clinical %>%
  dplyr::select(
    barcode, patient,
    project_id,
    age_at_index, gender,
    vital_status, days_to_death, days_to_last_follow_up,
    primary_diagnosis,
    paper_IDH.status,
    paper_Grade,
    paper_MGMT.promoter.status,
    paper_IDH.codel.subtype,
    paper_Telomere.Maintenance,
    paper_Transcriptome.Subtype
  ) %>%
  dplyr::mutate(
    # 计算OS
    OS.time = ifelse(vital_status == "Dead",
                     as.numeric(days_to_death),
                     as.numeric(days_to_last_follow_up)),
    OS = ifelse(vital_status == "Dead", 1, 0),
    # 简化
    IDH_status = paper_IDH.status,
    Grade = paper_Grade,
    MGMT_status = paper_MGMT.promoter.status,
    Subtype = paper_Transcriptome.Subtype
  )

# 过滤有效样本（有生存信息）
clinical_valid <- clinical_clean %>%
  dplyr::filter(!is.na(OS.time) & OS.time > 0)

message(sprintf("Valid clinical samples: %d", nrow(clinical_valid)))

# =============================================================================
# 5. 下载突变数据
# =============================================================================
message("=== Downloading mutation data ===")

query_mut_gbm <- GDCquery(
  project = "TCGA-GBM",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access = "open"
)
GDCdownload(query_mut_gbm, directory = file.path(tcga_dir, "GDCdata"))
maf_gbm <- GDCprepare(query_mut_gbm, directory = file.path(tcga_dir, "GDCdata"))

query_mut_lgg <- GDCquery(
  project = "TCGA-LGG",
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access = "open"
)
GDCdownload(query_mut_lgg, directory = file.path(tcga_dir, "GDCdata"))
maf_lgg <- GDCprepare(query_mut_lgg, directory = file.path(tcga_dir, "GDCdata"))

# 合并MAF
maf_glioma <- rbind(maf_gbm, maf_lgg)

# =============================================================================
# 6. 先保存表达和突变数据（确保不因CNV失败而丢失）
# =============================================================================
message("=== Saving expression and mutation data ===")

save(
  expr_tpm_symbol,
  expr_counts_symbol,
  clinical_clean,
  clinical_valid,
  gene_info,
  file = file.path(DATA_PROC, "tcga_glioma_expression.RData")
)
message("Expression saved: ", file.path(DATA_PROC, "tcga_glioma_expression.RData"))

save(
  maf_glioma,
  file = file.path(DATA_PROC, "tcga_glioma_mutation.RData")
)
message("Mutation saved: ", file.path(DATA_PROC, "tcga_glioma_mutation.RData"))

# =============================================================================
# 7. 下载CNV数据（可选，失败不阻塞）
# =============================================================================
message("=== Downloading CNV data ===")

tryCatch({
  query_cnv_gbm <- GDCquery(
    project = "TCGA-GBM",
    data.category = "Copy Number Variation",
    data.type = "Gene Level Copy Number"
  )
  GDCdownload(query_cnv_gbm, directory = file.path(tcga_dir, "GDCdata"))
  cnv_gbm <- GDCprepare(query_cnv_gbm, directory = file.path(tcga_dir, "GDCdata"))

  query_cnv_lgg <- GDCquery(
    project = "TCGA-LGG",
    data.category = "Copy Number Variation",
    data.type = "Gene Level Copy Number"
  )
  GDCdownload(query_cnv_lgg, directory = file.path(tcga_dir, "GDCdata"))
  cnv_lgg <- GDCprepare(query_cnv_lgg, directory = file.path(tcga_dir, "GDCdata"))

  save(cnv_gbm, cnv_lgg,
       file = file.path(DATA_PROC, "tcga_glioma_cnv.RData"))
  message("CNV saved: ", file.path(DATA_PROC, "tcga_glioma_cnv.RData"))
}, error = function(e) {
  message("CNV download/processing failed: ", e$message)
  message("This is non-blocking. Expression and mutation data are already saved.")
})

message("=== TCGA data download and processing completed ===")
message("Expression: ", file.path(DATA_PROC, "tcga_glioma_expression.RData"))
message("Mutation:   ", file.path(DATA_PROC, "tcga_glioma_mutation.RData"))
