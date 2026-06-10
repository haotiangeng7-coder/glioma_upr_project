###############################################################################
# download_scrna.R
# 下载单细胞RNA-seq数据
# GSE131928: Neftel et al. 2019 Cell — 28例胶质瘤 ~24,000+细胞
###############################################################################

source("00_setup/config.R")
library(GEOquery)

scrna_dir <- file.path(DATA_RAW, "scRNA")
dir.create(scrna_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. GSE131928 — 单细胞数据
# =============================================================================
message("=== Downloading GSE131928 (scRNA-seq) ===")
message("Note: This dataset may be large. Download time depends on network speed.")

# 方法1: 通过GEO下载supplementary files
# GSE131928包含处理后的count matrix
tryCatch({
  # 获取GEO元数据
  gse <- getGEO("GSE131928", destdir = scrna_dir, GSEMatrix = TRUE)

  # 下载supplementary files（包含count matrix）
  supp_files <- getGEOSuppFiles("GSE131928", baseDir = scrna_dir)
  message("Downloaded supplementary files:")
  print(rownames(supp_files))

}, error = function(e) {
  message("GEO download error: ", e$message)
  message("\nAlternative: Download manually from GEO")
  message("URL: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131928")
})

# =============================================================================
# 2. 处理下载的数据
# =============================================================================

# GSE131928通常包含以下supplementary文件:
# - 10X格式的count matrix (barcodes.tsv, genes.tsv, matrix.mtx)
# - 或者处理后的TPM/count文本文件

# 检查下载的文件
downloaded_files <- list.files(
  path = file.path(scrna_dir, "GSE131928"),
  pattern = "\\.(gz|txt|csv|tsv|h5|mtx)$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(downloaded_files) > 0) {
  message("\nDownloaded files:")
  for (f in downloaded_files) {
    message("  ", basename(f), " (", round(file.size(f) / 1e6, 1), " MB)")
  }
} else {
  message("\nNo supplementary files found. Trying alternative download...")

  # 方法2: 直接从NCBI FTP下载
  # Neftel 2019数据也可能在Broad Single Cell Portal
  message("
=== MANUAL DOWNLOAD INSTRUCTIONS ===

Option A: GEO
  1. Go to: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131928
  2. Download supplementary files (count matrix)
  3. Place in: ", scrna_dir, "/GSE131928/

Option B: Broad Single Cell Portal
  1. Search for 'Neftel glioma' at https://singlecell.broadinstitute.org/
  2. Download expression matrix and metadata

Option C: Use processed data
  The Smart-seq2 data (GSE131928_10X) provides:
  - IDH-wildtype GBM samples
  - IDH-mutant glioma samples
  - Already annotated with cell states (MES/AC/OPC/NPC)

Place downloaded files in: ", scrna_dir, "
Then re-run this script.
")
}

# =============================================================================
# 3. 解压并预处理（如果文件存在）
# =============================================================================

# 解压gz文件
gz_files <- list.files(
  path = file.path(scrna_dir, "GSE131928"),
  pattern = "\\.gz$",
  recursive = TRUE,
  full.names = TRUE
)

for (gz in gz_files) {
  out_file <- sub("\\.gz$", "", gz)
  if (!file.exists(out_file)) {
    message("Decompressing: ", basename(gz))
    R.utils::gunzip(gz, destname = out_file, remove = FALSE)
  }
}

# 检查10X格式文件
mtx_files <- list.files(
  path = file.path(scrna_dir, "GSE131928"),
  pattern = "matrix\\.mtx",
  recursive = TRUE,
  full.names = TRUE
)

if (length(mtx_files) > 0) {
  message("\n10X format data found. Will be loaded in 03_single_cell/01_preprocessing.R")
} else {
  # 检查文本格式的count matrix
  txt_files <- list.files(
    path = file.path(scrna_dir, "GSE131928"),
    pattern = "\\.(txt|csv|tsv)$",
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(txt_files) > 0) {
    message("\nText format expression matrix found:")
    for (f in txt_files) {
      message("  ", basename(f))
    }
  }
}

message("\n=== scRNA-seq data download script completed ===")
message("Next step: Run 03_single_cell/01_preprocessing.R")
