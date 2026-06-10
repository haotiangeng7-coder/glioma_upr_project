###############################################################################
# config.R
# 项目全局配置：路径、参数、颜色方案
###############################################################################

# === 项目根目录 ===
PROJECT_DIR <- getwd()

# === 数据目录 ===
DATA_RAW    <- file.path(PROJECT_DIR, "data", "raw")
DATA_PROC   <- file.path(PROJECT_DIR, "data", "processed")
FIG_DIR     <- file.path(PROJECT_DIR, "figures")
RES_DIR     <- file.path(PROJECT_DIR, "results")

# 创建子目录
for (d in c(DATA_RAW, DATA_PROC, FIG_DIR, RES_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# === 随机种子（保证可重复性）===
SEED <- 42

# === 可视化参数 ===
# 论文级别图形参数
THEME_PUBLICATION <- ggplot2::theme_bw() +
  ggplot2::theme(
    plot.title   = ggplot2::element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title   = ggplot2::element_text(size = 12),
    axis.text    = ggplot2::element_text(size = 10),
    legend.title = ggplot2::element_text(size = 11),
    legend.text  = ggplot2::element_text(size = 10),
    panel.grid.minor = ggplot2::element_blank()
  )

# 颜色方案
# 评审修订B2: 亚型按预后命名，不按UPR表达水平
# UPR-high-risk(预后最差) = 红色, UPR-intermediate = teal, UPR-favorable(预后最好) = 蓝色
COLORS_SUBTYPE <- c(
  "UPR-high-risk"    = "#E64B35",
  "UPR-intermediate" = "#00A087",
  "UPR-favorable"    = "#4DBBD5"
)

COLORS_RISK <- c(
  "High" = "#E64B35",
  "Low"  = "#4DBBD5"
)

COLORS_CELLTYPE <- c(
  "Malignant"      = "#E64B35",
  "TAM"            = "#4DBBD5",
  "Microglia"      = "#00A087",
  "CD4 T"          = "#3C5488",
  "CD8 T"          = "#F39B7F",
  "Treg"           = "#8491B4",
  "DC"             = "#91D1C2",
  "Oligodendrocyte" = "#B09C85",
  "Astrocyte"      = "#7E6148",
  "Endothelial"    = "#FFDC91",
  "MDSC"           = "#DC0000"
)

# === 分析参数 ===
# 单细胞
SC_MIN_FEATURES <- 200
SC_MAX_FEATURES <- 6000
SC_MAX_MT_PCT   <- 20
SC_N_PCS        <- 30
SC_RESOLUTION   <- 0.8
SC_LOGFC_CUTOFF <- 0.5       # 评审修订：从0.25提高到0.5
SC_MIN_CELLS_EXPLORATORY <- 50  # 低于此值的细胞类型标注exploratory
SC_MIN_MALIGNANT_FOR_UPR <- 100 # CellChat样本级UPR评分所需最少恶性细胞数

# SCTransform参数
SCT_VST_FLAVOR  <- "v2"     # Choudhary & Satija 2022
SCT_NFEATURES   <- 3000

# Harmony参数
HARMONY_THETA   <- 2
HARMONY_LAMBDA  <- 1
HARMONY_MAXITER <- 20

# 细胞周期混杂阈值
CC_COR_THRESHOLD <- 0.3     # |r|>0.3时回归掉细胞周期效应

# Bulk分析
PVALUE_CUTOFF   <- 0.05
FDR_CUTOFF      <- 0.05
LOGFC_CUTOFF    <- 1.0
CORR_EFFECT_THRESHOLD <- 0.3  # 评审修订：|r|>=0.3才视为有生物学意义的相关性

# 共识聚类
CC_MAX_K        <- 6
CC_REPS         <- 1000
CC_PITEM        <- 0.8
CC_MAD_PERCENTILE <- 0.5    # 评审修订：MAD前50%筛选聚类基因

# 机器学习
ML_NFOLDS       <- 10
ML_NREPEATS     <- 100
ML_NESTED_OUTER <- 5        # 评审修订：嵌套CV外层fold数
ML_NESTED_INNER <- 5        # 评审修订：嵌套CV内层fold数
ML_NPERM        <- 1000     # 评审修订：置换检验次数

# === 基因集名称（用于引用） ===
UPR_GENESETS <- list(
  IRE1_XBP1 = "IRE1a_XBP1_pathway",
  PERK_ATF4 = "PERK_ATF4_pathway",
  ATF6      = "ATF6_pathway",
  UPR_BROAD = "UPR_broad_geneset"
)

# === 数据集信息 ===
DATASETS <- list(
  training = list(
    name = "TCGA-GBM/LGG",
    source = "GDC Portal"
  ),
  tuning_validation = list(
    CGGA_batch2 = list(name = "CGGA mRNA-seq batch2", n = 693, role = "tuning")
  ),
  held_out_test = list(
    CGGA_batch1 = list(name = "CGGA mRNA-seq batch1", n = 325, role = "blind_test"),
    GSE16011    = list(name = "GSE16011", n = 284, role = "blind_test")
  ),
  exploratory = list(
    GSE43378    = list(name = "GSE43378", n = 50, role = "exploratory")
  ),
  single_cell = list(
    GSE131928 = list(name = "GSE131928", desc = "Neftel et al. 2019 Cell")
  ),
  immunotherapy = list(
    Zhao2019 = list(name = "Zhao et al. 2019", n = 66, treatment = "anti-PD1")
  )
)

message("=== Project config loaded successfully ===")
message("PROJECT_DIR: ", PROJECT_DIR)
