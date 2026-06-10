###############################################################################
# download_cgga.R
# 下载CGGA (Chinese Glioma Genome Atlas) 数据
# batch1 (325例) + batch2 (693例)
# 来源：http://www.cgga.org.cn/
###############################################################################

source("00_setup/config.R")
library(data.table)
library(dplyr)

cgga_dir <- file.path(DATA_RAW, "CGGA")
dir.create(cgga_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# CGGA数据需要手动下载或通过API获取
# 下载地址: http://www.cgga.org.cn/download.jsp
# 需要下载:
#   1. CGGA.mRNAseq_693.RSEM-genes.20200506.txt (batch2 表达矩阵)
#   2. CGGA.mRNAseq_693_clinical.20200506.txt   (batch2 临床数据)
#   3. CGGA.mRNAseq_325.RSEM-genes.20200506.txt (batch1 表达矩阵)
#   4. CGGA.mRNAseq_325_clinical.20200506.txt   (batch1 临床数据)
# =============================================================================

message("=== CGGA Data Processing ===")
message("Please download CGGA data from http://www.cgga.org.cn/download.jsp")
message("Place files in: ", cgga_dir)
message("")
message("Required files:")
message("  - CGGA.mRNAseq_693.RSEM-genes.20200506.txt")
message("  - CGGA.mRNAseq_693_clinical.20200506.txt")
message("  - CGGA.mRNAseq_325.RSEM-genes.20200506.txt")
message("  - CGGA.mRNAseq_325_clinical.20200506.txt")

# 文件已手动下载，跳过自动下载（CGGA网站经常不可用）
# 如果文件是zip格式，先解压
zip_files <- list.files(cgga_dir, pattern = "\\.zip$", full.names = TRUE)
if (length(zip_files) > 0) {
  for (zf in zip_files) {
    message("Unzipping: ", basename(zf))
    unzip(zf, exdir = cgga_dir, overwrite = TRUE)
  }
}

# =============================================================================
# 处理CGGA数据（在文件存在的情况下）
# =============================================================================

process_cgga <- function(expr_file, clin_file, batch_name) {
  if (!file.exists(expr_file)) {
    message(sprintf("File not found: %s. Skipping %s.", expr_file, batch_name))
    return(NULL)
  }

  message(sprintf("Processing %s...", batch_name))

  # 读取表达矩阵
  expr <- fread(expr_file, data.table = FALSE)
  rownames(expr) <- expr[, 1]
  expr <- expr[, -1]
  expr <- as.matrix(expr)

  # 读取临床数据
  clin <- fread(clin_file, data.table = FALSE)

  # 整理临床数据 — 自动检测列名（不同版本列名可能不同）
  message("  Clinical columns: ", paste(colnames(clin), collapse = ", "))

  # 检测Censor列（可能有多种命名方式）
  censor_col <- grep("Censor|censor|status.*dead|vital", colnames(clin), value = TRUE)[1]
  os_col <- grep("^OS$", colnames(clin), value = TRUE)[1]
  idh_col <- grep("IDH_mutation|IDH.mutation|IDH_status", colnames(clin), value = TRUE)[1]
  codel_col <- grep("1p19q|codeletion", colnames(clin), value = TRUE)[1]
  mgmt_col <- grep("MGMT|mgmt", colnames(clin), value = TRUE)[1]

  clin_clean <- clin
  clin_clean$Sample_ID <- clin$CGGA_ID
  clin_clean$OS.time <- as.numeric(clin[[os_col]])
  clin_clean$OS_status <- as.numeric(clin[[censor_col]])
  clin_clean$Age <- as.numeric(clin$Age)
  if (!is.null(idh_col)) clin_clean$IDH_status <- clin[[idh_col]]
  if (!is.null(codel_col)) clin_clean$Codel_status <- clin[[codel_col]]
  if (!is.null(mgmt_col)) clin_clean$MGMT_status <- clin[[mgmt_col]]

  clin_clean <- clin_clean %>%
    dplyr::filter(!is.na(OS.time) & OS.time > 0)

  message(sprintf("  %s: %d genes x %d samples, %d with valid survival",
                   batch_name, nrow(expr), ncol(expr), nrow(clin_clean)))

  return(list(expr = expr, clinical = clin_clean))
}

# 处理batch2
cgga_b2 <- process_cgga(
  expr_file = file.path(cgga_dir, "CGGA.mRNAseq_693.RSEM-genes.20200506.txt"),
  clin_file = file.path(cgga_dir, "CGGA.mRNAseq_693_clinical.20200506.txt"),
  batch_name = "CGGA_batch2"
)

# 处理batch1
cgga_b1 <- process_cgga(
  expr_file = file.path(cgga_dir, "CGGA.mRNAseq_325.RSEM-genes.20200506.txt"),
  clin_file = file.path(cgga_dir, "CGGA.mRNAseq_325_clinical.20200506.txt"),
  batch_name = "CGGA_batch1"
)

# 保存
if (!is.null(cgga_b1) || !is.null(cgga_b2)) {
  save(cgga_b1, cgga_b2,
       file = file.path(DATA_PROC, "cgga_data.RData"))
  message("CGGA data saved to: ", file.path(DATA_PROC, "cgga_data.RData"))
} else {
  message("No CGGA data processed. Please download files first.")
}
