###############################################################################
# install_packages.R
# 安装所有必需的R包
###############################################################################

# --- CRAN packages ---
cran_pkgs <- c(
  "survminer", "timeROC", "rms", "randomForestSRC",
  "CoxBoost", "superpc", "plsRcox",
  "survival", "glmnet", "caret", "gbm",
  "ggplot2", "ggpubr", "pheatmap", "circlize", "cowplot", "patchwork",
  "dplyr", "tidyverse", "data.table", "reshape2",
  "msigdbr", "ROCR", "pROC",
  "forestplot", "DT",
  # 评审修订新增
  "survRM2",         # RMST分析（PH假设违反时的替代）
  "mclust",          # GMM分组（UPR-high/low敏感性分析）
  "rcompanion",      # Cramér's V（Verhaak亚型混杂评估）
  "SoupX"            # Ambient RNA校正（10X数据）
)

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing CRAN package: ", pkg)
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# --- Bioconductor packages ---
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

bioc_pkgs <- c(
  "ConsensusClusterPlus", "maftools", "AUCell",
  "SingleR", "celldex", "TCGAbiolinks",
  "GEOquery", "Seurat", "harmony",
  "clusterProfiler", "GSVA", "fgsea",
  "ComplexHeatmap", "edgeR", "limma", "DESeq2", "sva",
  "AnnotationDbi", "org.Hs.eg.db", "preprocessCore",
  "MCPcounter", "EPIC", "scater", "UCell",
  "SummarizedExperiment", "GenomicRanges",
  # 评审修订新增
  "scDblFinder",     # Doublet检测（10X数据）
  "scuttle",         # scDblFinder依赖
  "STRINGdb"         # PPI网络分析
)

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing Bioconductor package: ", pkg)
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
  }
}

# --- GitHub packages ---
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

# CellChat
if (!requireNamespace("CellChat", quietly = TRUE)) {
  devtools::install_github("jinworks/CellChat")
}

# oncoPredict
if (!requireNamespace("oncoPredict", quietly = TRUE)) {
  devtools::install_github("danieldanciu/oncoPredict")
}

# immunedeconv
if (!requireNamespace("immunedeconv", quietly = TRUE)) {
  devtools::install_github("omnideconv/immunedeconv")
}

# estimate
if (!requireNamespace("estimate", quietly = TRUE)) {
  # ESTIMATE需要特殊安装
  utils::install.packages("estimate",
    repos = "https://r-forge.r-project.org",
    type = "source"
  )
}

# xCell
if (!requireNamespace("xCell", quietly = TRUE)) {
  devtools::install_github("dviraran/xCell")
}

message("=== All packages installation completed ===")
message("Run 00_setup/config.R to verify installation")
