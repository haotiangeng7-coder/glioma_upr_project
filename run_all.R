###############################################################################
# run_all.R
# 主运行脚本 — 按顺序执行所有分析
#
# 脑胶质瘤UPR驱动免疫微环境重塑研究
# Unfolded Protein Response-Driven Immune Remodeling in Glioma
###############################################################################
#
# IMPORTANT: Run this script from the repository root directory,
# e.g.: Rscript run_all.R
# All relative paths in source() calls assume working directory = repo root.

# 使用说明:
# 1. 先运行 00_setup/install_packages.R 安装所有必需包
# 2. 按顺序运行以下脚本，或直接运行本脚本
# 3. 部分脚本需要先下载数据（见01_data_download/）
# 4. 外部工具结果（TIDE, CMap）需要手动上传和下载
#
# 项目结构:
#   00_setup/           - 环境配置
#   01_data_download/   - 数据下载
#   02_gene_sets/       - UPR基因集定义
#   03_single_cell/     - Part 1: 单细胞分析
#   04_bulk_subtyping/  - Part 2: Bulk分子分型
#   05_ml_model/        - Part 3: 机器学习预后模型
#   06_clinical_translation/ - Part 4: 临床转化分析
#   figures/            - 输出图表
#   results/            - 输出结果表格
#   data/raw/           - 原始数据
#   data/processed/     - 处理后数据
#
###############################################################################

cat("
╔══════════════════════════════════════════════════════════════════╗
║     Glioma UPR-Immune Remodeling Analysis Pipeline              ║
║     脑胶质瘤UPR驱动免疫微环境重塑研究                            ║
╚══════════════════════════════════════════════════════════════════╝
\n")

# ============ 计时开始 ============
start_time <- Sys.time()

# ============ Step 0: 环境配置 ============
cat("\n=== Step 0: Configuration ===\n")
source("00_setup/config.R")

# ============ Step 1: 基因集定义 ============
cat("\n=== Step 1: UPR Gene Sets ===\n")
source("02_gene_sets/upr_gene_sets.R")

# ============ Step 2: 数据下载 ============
cat("\n=== Step 2: Data Download ===\n")
cat("  2.1 TCGA data...\n")
source("01_data_download/download_tcga.R")

cat("  2.2 CGGA data...\n")
source("01_data_download/download_cgga.R")

cat("  2.3 GEO data...\n")
source("01_data_download/download_geo.R")

cat("  2.4 scRNA-seq data...\n")
source("01_data_download/download_scrna.R")

# ============ Part 1: 单细胞分析 (Figures 1-2) ============
cat("\n╔═══ Part 1: Single-Cell Analysis ═══╗\n")

cat("  1.1 Preprocessing...\n")
source("03_single_cell/01_preprocessing.R")

cat("  1.2 Cell annotation...\n")
source("03_single_cell/02_cell_annotation.R")

cat("  1.3 UPR scoring...\n")
source("03_single_cell/03_upr_scoring.R")

cat("  1.4 Differential analysis...\n")
source("03_single_cell/04_differential_analysis.R")

cat("  1.5 CellChat analysis...\n")
source("03_single_cell/05_cellchat.R")

# ============ Part 2: Bulk分子分型 (Figures 3-4) ============
cat("\n╔═══ Part 2: Bulk Molecular Subtyping ═══╗\n")

cat("  2.1 UPR landscape...\n")
source("04_bulk_subtyping/01_upr_landscape.R")

cat("  2.2 Consensus clustering...\n")
source("04_bulk_subtyping/02_consensus_clustering.R")

cat("  2.3 Clinical characterization...\n")
source("04_bulk_subtyping/03_clinical_characterization.R")

cat("  2.4 Immune analysis...\n")
source("04_bulk_subtyping/04_immune_analysis.R")

cat("  2.5 Genomic analysis...\n")
source("04_bulk_subtyping/05_genomic_analysis.R")

cat("  2.6 Pathway analysis...\n")
source("04_bulk_subtyping/06_pathway_analysis.R")

# ============ Part 3: 机器学习模型 (Figures 5-6) ============
cat("\n╔═══ Part 3: Machine Learning Model ═══╗\n")

cat("  3.1 Feature selection...\n")
source("05_ml_model/01_feature_selection.R")

cat("  3.2 ML combinations...\n")
source("05_ml_model/02_ml_combinations.R")

cat("  3.3 Risk score...\n")
source("05_ml_model/03_risk_score.R")

cat("  3.4 Independent prognosis...\n")
source("05_ml_model/04_independent_prognosis.R")

cat("  3.5 Nomogram...\n")
source("05_ml_model/05_nomogram.R")

# ============ Part 4: 临床转化 (Figures 7-8) ============
cat("\n╔═══ Part 4: Clinical Translation ═══╗\n")

cat("  4.1 Immunotherapy prediction...\n")
source("06_clinical_translation/01_immunotherapy_prediction.R")

cat("  4.2 Drug sensitivity...\n")
source("06_clinical_translation/02_drug_sensitivity.R")

cat("  4.3 CMap analysis...\n")
source("06_clinical_translation/03_cmap_analysis.R")

cat("  4.4 Hub gene analysis...\n")
source("06_clinical_translation/04_hub_gene_analysis.R")

# ============ 完成 ============
end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "hours")

cat(sprintf("\n
╔══════════════════════════════════════════════════════════════════╗
║                    ANALYSIS COMPLETED                           ║
║  Total time: %.1f hours                                         ║
║  Figures: %s                                                    ║
║  Results: %s                                                    ║
╚══════════════════════════════════════════════════════════════════╝
\n",
as.numeric(elapsed),
FIG_DIR,
RES_DIR))

# 列出生成的文件
cat("Generated figures:\n")
fig_files <- list.files(FIG_DIR, pattern = "\\.pdf$")
for (f in fig_files) cat(sprintf("  - %s\n", f))

cat("\nGenerated result files:\n")
res_files <- list.files(RES_DIR, pattern = "\\.(csv|RData)$")
for (f in res_files) cat(sprintf("  - %s\n", f))
