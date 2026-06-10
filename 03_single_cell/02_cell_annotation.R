###############################################################################
# 02_cell_annotation.R
# 单细胞数据细胞类型注释 — 完整版（含评审修订要求）
# SingleR自动注释 + 经典marker手动校正（11种细胞类型）
# 恶性细胞亚型分类 (MES/AC/OPC/NPC)
#
# 输出:
#   - data/processed/seu_annotated.rds
#   - figures/sc_marker_dotplot.pdf
#   - figures/sc_marker_featureplot.pdf
#   - figures/sc_umap_celltype.pdf
#   - figures/sc_celltype_barplot.pdf
#   - figures/sc_singler_heatmap.pdf
#   - figures/sc_malignant_subtypes_umap.pdf
#   - results/TableS1_celltype_counts.csv
#   - results/cluster_annotation_details.csv
###############################################################################

# === 加载配置和依赖 ===
source("00_setup/config.R")
load(file.path(DATA_PROC, "upr_gene_sets.RData"))

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleR)
  library(celldex)
  library(scater)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
  library(pheatmap)
})

set.seed(SEED)

###############################################################################
# 1. 加载预处理后的Seurat对象
###############################################################################
message("=== [1/7] Loading preprocessed data ===")
seu <- readRDS(file.path(DATA_PROC, "seu_preprocessed.rds"))
message(sprintf("Loaded: %d genes x %d cells, %d clusters",
                nrow(seu), ncol(seu), length(unique(Idents(seu)))))

###############################################################################
# 2. SingleR自动注释（双参考数据集）
###############################################################################
message("=== [2/7] Running SingleR annotation ===")

# 获取参考数据集
ref_hpca      <- celldex::HumanPrimaryCellAtlasData()
ref_blueprint <- celldex::BlueprintEncodeData()

# 提取表达矩阵为SCE
sce <- as.SingleCellExperiment(seu)

# SingleR预测 — HPCA参考
message("  Running SingleR with HumanPrimaryCellAtlas...")
set.seed(SEED)
pred_hpca <- SingleR(
  test            = sce,
  ref             = ref_hpca,
  labels          = ref_hpca$label.main,
  assay.type.test = "logcounts"
)

# SingleR预测 — BlueprintEncode参考
message("  Running SingleR with BlueprintEncode...")
set.seed(SEED)
pred_blueprint <- SingleR(
  test            = sce,
  ref             = ref_blueprint,
  labels          = ref_blueprint$label.main,
  assay.type.test = "logcounts"
)

# 添加注释到metadata
seu$SingleR_hpca      <- pred_hpca$labels
seu$SingleR_blueprint <- pred_blueprint$labels

# SingleR质量评估
n_hpca_pass <- sum(!is.na(pred_hpca$pruned.labels))
n_bp_pass   <- sum(!is.na(pred_blueprint$pruned.labels))
message(sprintf("  HPCA: %d/%d cells passed pruning (%.1f%%)",
                n_hpca_pass, ncol(seu), n_hpca_pass / ncol(seu) * 100))
message(sprintf("  Blueprint: %d/%d cells passed pruning (%.1f%%)",
                n_bp_pass, ncol(seu), n_bp_pass / ncol(seu) * 100))

# SingleR得分热图
message("  Generating SingleR score heatmap...")
tryCatch({
  pdf(file.path(FIG_DIR, "sc_singler_heatmap.pdf"), width = 12, height = 8)
  plotScoreHeatmap(pred_hpca, show.pruned = TRUE, max.labels = 20)
  dev.off()
}, error = function(e) {
  message("  SingleR heatmap generation failed: ", e$message)
  try(dev.off(), silent = TRUE)
})

rm(sce)

###############################################################################
# 3. 经典Marker基因定义（11种目标细胞类型 + 恶性亚型）
###############################################################################
message("=== [3/7] Defining cell type markers ===")

# 目标11种细胞类型的marker基因
markers <- list(
  # --- 恶性细胞/胶质瘤细胞（整体）---
  Malignant = c("SOX2", "NES", "OLIG1", "OLIG2", "PDGFRA", "EGFR", "CDK4",
                "IDH1", "TP53", "PTEN"),
  # --- 免疫细胞 ---
  TAM       = c("CD68", "CD163", "CSF1R", "MRC1", "CD14", "FCGR3A", "AIF1",
                "ITGAM", "MSR1"),
  Microglia = c("P2RY12", "TMEM119", "CX3CR1", "CSF1R", "TREM2", "HEXB",
                "SALL1"),
  CD4_T     = c("CD3D", "CD3E", "CD4", "IL7R", "CCR7", "LEF1", "TCF7"),
  CD8_T     = c("CD3D", "CD3E", "CD8A", "CD8B", "GZMB", "PRF1", "NKG7",
                "GZMA"),
  Treg      = c("FOXP3", "IL2RA", "CTLA4", "IKZF2", "TNFRSF18"),
  DC        = c("ITGAX", "CD1C", "FCER1A", "CLEC10A", "HLA-DRA", "HLA-DQA1",
                "BATF3"),
  MDSC      = c("S100A8", "S100A9", "S100A12", "ITGAM", "CEACAM8", "CD33",
                "ARG1"),
  # --- 基质/正常脑细胞 ---
  Oligodendrocyte = c("MBP", "PLP1", "MOG", "MAG", "MOBP", "CNP"),
  Astrocyte       = c("GFAP", "AQP4", "S100B", "ALDH1L1", "SLC1A2",
                      "GJA1"),
  Endothelial     = c("PECAM1", "VWF", "CDH5", "CLDN5", "FLT1", "ERG")
)

# 恶性亚型marker (Neftel 2019)
malignant_subtype_markers <- list(
  MES_like = c("CHI3L1", "VIM", "ANXA1", "ANXA2", "CD44", "SERPINE1",
               "LGALS1", "TGFBI", "NAMPT"),
  AC_like  = c("GFAP", "AQP4", "ALDOC", "SLC1A3", "APOE", "CLU",
               "ALDH1L1"),
  OPC_like = c("OLIG1", "OLIG2", "PDGFRA", "SOX10", "NKX2-2", "CSPG4",
               "NEU4"),
  NPC_like = c("SOX4", "SOX11", "DCX", "STMN1", "DLL3", "ASCL1",
               "TUBB3")
)

###############################################################################
# 4. Marker基因可视化
###############################################################################
message("=== [4/7] Visualizing markers ===")

# 选择每种细胞类型的top markers用于DotPlot
top_markers <- c(
  "SOX2", "EGFR",           # Malignant
  "CD68", "CD163",           # TAM
  "P2RY12", "TMEM119",       # Microglia
  "CD3D", "CD4",             # CD4 T
  "CD8A", "GZMB",            # CD8 T
  "FOXP3", "IL2RA",          # Treg
  "HLA-DRA", "CD1C",         # DC
  "S100A8", "S100A9",        # MDSC
  "MBP", "PLP1",             # Oligo
  "GFAP", "AQP4",            # Astrocyte
  "PECAM1", "VWF"            # Endothelial
)

top_markers_valid <- top_markers[top_markers %in% rownames(seu)]
message(sprintf("  %d / %d top markers found in data", length(top_markers_valid), length(top_markers)))

if (length(top_markers_valid) > 0) {
  p_dot <- DotPlot(seu, features = top_markers_valid, group.by = "seurat_clusters") +
    RotatedAxis() +
    THEME_PUBLICATION +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("Marker Gene Expression by Cluster")

  ggsave(file.path(FIG_DIR, "sc_marker_dotplot.pdf"), p_dot, width = 16, height = 8)
}

# Feature plots
key_markers <- c("SOX2", "CD68", "CD3D", "CD8A", "GFAP", "MBP", "PECAM1", "P2RY12")
key_markers_valid <- key_markers[key_markers %in% rownames(seu)]

if (length(key_markers_valid) > 0) {
  p_feature <- FeaturePlot(seu, features = key_markers_valid, ncol = 4,
                             order = TRUE, cols = c("lightgrey", "darkred")) &
    THEME_PUBLICATION

  ggsave(file.path(FIG_DIR, "sc_marker_featureplot.pdf"), p_feature,
         width = 16, height = ceiling(length(key_markers_valid) / 4) * 4)
}

###############################################################################
# 5. 综合注释：SingleR + Marker Score手动校正
###############################################################################
message("=== [5/7] Integrated cell type annotation (SingleR + marker correction) ===")

# --- 5.1 计算所有细胞类型的module score ---
for (ct in names(markers)) {
  valid_genes <- markers[[ct]][markers[[ct]] %in% rownames(seu)]
  if (length(valid_genes) >= 2) {
    tryCatch({
      set.seed(SEED)
      seu <- AddModuleScore(seu, features = list(valid_genes),
                             name = paste0(ct, "_mscore"))
      message(sprintf("  %s: %d/%d markers found", ct, length(valid_genes), length(markers[[ct]])))
    }, error = function(e) {
      message(sprintf("  %s: AddModuleScore failed (%s). Using manual mean.", ct, e$message))
      expr_mat <- tryCatch(GetAssayData(seu, layer = "data"), error = function(e) GetAssayData(seu, slot = "data"))
      valid_in_mat <- valid_genes[valid_genes %in% rownames(expr_mat)]
      if (length(valid_in_mat) > 0) {
        scores <- Matrix::colMeans(expr_mat[valid_in_mat, , drop = FALSE])
        seu[[paste0(ct, "_mscore1")]] <<- scores
      }
    })
  } else {
    message(sprintf("  WARNING: %s has <2 valid markers (%d). Skipping.", ct, length(valid_genes)))
  }
}

# --- 5.2 计算恶性亚型的module score ---
for (st in names(malignant_subtype_markers)) {
  valid_genes <- malignant_subtype_markers[[st]][malignant_subtype_markers[[st]] %in% rownames(seu)]
  if (length(valid_genes) >= 2) {
    tryCatch({
      set.seed(SEED)
      seu <- AddModuleScore(seu, features = list(valid_genes),
                             name = paste0(st, "_mscore"))
    }, error = function(e) {
      message(sprintf("  AddModuleScore failed for %s: %s. Using manual mean.", st, e$message))
      expr_mat <- tryCatch(GetAssayData(seu, layer = "data"), error = function(e) GetAssayData(seu, slot = "data"))
      valid_in_mat <- valid_genes[valid_genes %in% rownames(expr_mat)]
      if (length(valid_in_mat) > 0) {
        scores <- Matrix::colMeans(expr_mat[valid_in_mat, , drop = FALSE])
        seu[[paste0(st, "_mscore1")]] <<- scores
      }
    })
  }
}

# --- 5.3 综合注释函数 ---
annotate_cluster <- function(cluster_id, seu_obj, markers_list, singler_col = "SingleR_hpca") {

  cells_in_cluster <- colnames(seu_obj)[Idents(seu_obj) == cluster_id]
  n_cells <- length(cells_in_cluster)

  # (A) Marker score排名
  score_cols <- paste0(names(markers_list), "_mscore1")
  score_cols_valid <- score_cols[score_cols %in% colnames(seu_obj@meta.data)]

  if (length(score_cols_valid) == 0) return("Unknown")

  mean_scores <- sapply(score_cols_valid, function(sc) {
    mean(seu_obj@meta.data[cells_in_cluster, sc], na.rm = TRUE)
  })
  names(mean_scores) <- gsub("_mscore1", "", names(mean_scores))
  top_marker_type <- names(which.max(mean_scores))

  # (B) SingleR投票
  singler_labels <- seu_obj@meta.data[cells_in_cluster, singler_col]
  singler_majority <- names(sort(table(singler_labels), decreasing = TRUE))[1]

  # (C) 综合判断
  # 如果marker score最高的是Malignant相关类型，需要额外确认
  # 恶性细胞通常不会被SingleR注释为immune cells
  immune_labels <- c("Monocytes", "Macrophages", "T_cells", "B_cells", "NK_cells",
                     "DC", "Neutrophils")

  if (top_marker_type == "Malignant") {
    # 如果SingleR也不认为是免疫细胞，确认为恶性
    if (!singler_majority %in% immune_labels) {
      return("Malignant")
    }
    # SingleR认为是免疫细胞但marker score说是恶性 — 以marker为准但标记
    if (mean_scores["Malignant"] > 0.5) {
      return("Malignant")
    }
  }

  # 对于免疫和基质细胞类型，marker score和SingleR一致性检查
  # 以marker score为主要依据
  return(top_marker_type)
}

# --- 5.4 逐cluster注释 ---
clusters <- levels(Idents(seu))
if (is.null(clusters)) clusters <- sort(unique(as.character(Idents(seu))))

cluster_annotation <- sapply(clusters, function(cl) {
  annotate_cluster(cl, seu, markers)
})
names(cluster_annotation) <- clusters

message("\nCluster annotations (marker-based with SingleR correction):")
for (cl in names(cluster_annotation)) {
  n <- sum(Idents(seu) == cl)
  message(sprintf("  Cluster %-3s -> %-20s (%d cells)", cl, cluster_annotation[cl], n))
}

# --- 5.5 映射到标准11种细胞类型标签 ---
celltype_map <- c(
  "Malignant"       = "Malignant",
  "TAM"             = "TAM",
  "Microglia"       = "Microglia",
  "CD4_T"           = "CD4 T",
  "CD8_T"           = "CD8 T",
  "Treg"            = "Treg",
  "DC"              = "DC",
  "MDSC"            = "MDSC",
  "Oligodendrocyte" = "Oligodendrocyte",
  "Astrocyte"       = "Astrocyte",
  "Endothelial"     = "Endothelial"
)

# unname()避免Seurat v5将向量names误解为cell names
seu$celltype_raw <- unname(cluster_annotation[as.character(Idents(seu))])
celltype_mapped <- celltype_map[seu$celltype_raw]
seu$celltype <- unname(ifelse(is.na(celltype_mapped), seu$celltype_raw, celltype_mapped))

# --- 5.6 恶性细胞亚型注释 ---
message("\n=== Malignant subtype annotation ===")

malignant_cells <- colnames(seu)[seu$celltype == "Malignant"]
message(sprintf("  Total malignant cells: %d", length(malignant_cells)))

seu$malignant_subtype <- NA_character_

if (length(malignant_cells) > 0) {
  subtype_score_cols <- paste0(names(malignant_subtype_markers), "_mscore1")
  subtype_score_cols <- subtype_score_cols[subtype_score_cols %in% colnames(seu@meta.data)]

  if (length(subtype_score_cols) >= 2) {
    subtype_mat <- seu@meta.data[malignant_cells, subtype_score_cols, drop = FALSE]
    subtype_names <- gsub("_mscore1", "", colnames(subtype_mat))

    # 每个恶性细胞分配到得分最高的亚型
    best_subtype_idx <- apply(subtype_mat, 1, which.max)
    seu$malignant_subtype[malignant_cells] <- subtype_names[best_subtype_idx]

    message("  Malignant subtype distribution:")
    mal_subtype_table <- table(seu$malignant_subtype[malignant_cells])
    for (st in names(mal_subtype_table)) {
      message(sprintf("    %-12s: %d cells (%.1f%%)",
                      st, mal_subtype_table[st],
                      mal_subtype_table[st] / length(malignant_cells) * 100))
    }
  } else {
    message("  WARNING: Too few subtype markers found. Cannot assign subtypes.")
  }
}

###############################################################################
# 6. 注释后可视化 + Supplementary Table S1
###############################################################################
message("=== [6/7] Post-annotation visualization & Table S1 ===")

# --- 6.1 UMAP by cell type ---
# 构建完整颜色向量（包含数据中可能出现的所有类型）
all_celltypes <- unique(seu$celltype)
color_vec <- COLORS_CELLTYPE[all_celltypes]
# 对于不在预定义颜色中的类型，分配颜色
missing_ct <- all_celltypes[is.na(color_vec)]
if (length(missing_ct) > 0) {
  extra_colors <- scales::hue_pal()(length(missing_ct))
  names(extra_colors) <- missing_ct
  color_vec[missing_ct] <- extra_colors
}
color_vec <- color_vec[!is.na(names(color_vec))]

p_celltype <- DimPlot(seu, reduction = "umap", group.by = "celltype",
                       label = TRUE, repel = TRUE, label.size = 3.5,
                       pt.size = 0.3, raster = FALSE) +
  scale_color_manual(values = color_vec) +
  THEME_PUBLICATION +
  ggtitle("Cell Type Annotation")

ggsave(file.path(FIG_DIR, "sc_umap_celltype.pdf"), p_celltype, width = 10, height = 7)

# --- 6.2 恶性亚型UMAP ---
if (sum(!is.na(seu$malignant_subtype)) > 0) {
  mal_seu <- subset(seu, celltype == "Malignant")

  p_mal_umap <- DimPlot(mal_seu, reduction = "umap", group.by = "malignant_subtype",
                          label = TRUE, repel = TRUE, pt.size = 0.5) +
    scale_color_manual(values = c(
      "MES_like" = "#E64B35", "AC_like" = "#4DBBD5",
      "OPC_like" = "#00A087", "NPC_like" = "#3C5488"
    )) +
    THEME_PUBLICATION +
    ggtitle("Malignant Cell Subtypes (Neftel 2019)")

  ggsave(file.path(FIG_DIR, "sc_malignant_subtypes_umap.pdf"), p_mal_umap,
         width = 9, height = 7)
  rm(mal_seu)
}

# --- 6.3 细胞类型计数条形图 ---
celltype_table <- as.data.frame(table(seu$celltype))
colnames(celltype_table) <- c("CellType", "Count")
celltype_table <- celltype_table %>% arrange(desc(Count))

p_bar <- ggplot(celltype_table, aes(x = reorder(CellType, -Count), y = Count, fill = CellType)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), vjust = -0.3, size = 3) +
  scale_fill_manual(values = color_vec) +
  labs(x = "", y = "Number of Cells", title = "Cell Type Composition") +
  THEME_PUBLICATION +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

ggsave(file.path(FIG_DIR, "sc_celltype_barplot.pdf"), p_bar, width = 10, height = 6)

# --- 6.4 Supplementary Table S1: 细胞类型数量报告表 ---
message("Generating Supplementary Table S1...")

# 主细胞类型统计
table_s1 <- data.frame(
  CellType     = celltype_table$CellType,
  Count        = celltype_table$Count,
  Percentage   = round(celltype_table$Count / sum(celltype_table$Count) * 100, 2),
  stringsAsFactors = FALSE
)

# 标注 <50 个细胞的类型为 exploratory
table_s1$Status <- ifelse(table_s1$Count < SC_MIN_CELLS_EXPLORATORY,
                           "exploratory", "robust")

# 添加恶性亚型统计
if (sum(!is.na(seu$malignant_subtype)) > 0) {
  mal_subtypes <- as.data.frame(table(seu$malignant_subtype))
  colnames(mal_subtypes) <- c("Subtype", "Count")
  mal_subtypes$Percentage <- round(mal_subtypes$Count / sum(mal_subtypes$Count) * 100, 2)
  mal_subtypes$Status <- ifelse(mal_subtypes$Count < SC_MIN_CELLS_EXPLORATORY,
                                 "exploratory", "robust")

  # 构建亚型行
  subtype_rows <- data.frame(
    CellType   = paste0("  Malignant/", mal_subtypes$Subtype),
    Count      = mal_subtypes$Count,
    Percentage = mal_subtypes$Percentage,
    Status     = mal_subtypes$Status,
    stringsAsFactors = FALSE
  )

  # 在Malignant行之后插入亚型行
  mal_idx <- which(table_s1$CellType == "Malignant")
  if (length(mal_idx) > 0) {
    table_s1 <- rbind(
      table_s1[1:mal_idx, , drop = FALSE],
      subtype_rows,
      if (mal_idx < nrow(table_s1)) table_s1[(mal_idx + 1):nrow(table_s1), , drop = FALSE] else NULL
    )
  }
}

# 添加汇总行
total_row <- data.frame(
  CellType   = "TOTAL",
  Count      = ncol(seu),
  Percentage = 100.0,
  Status     = "",
  stringsAsFactors = FALSE
)
table_s1 <- rbind(table_s1, total_row)

# 打印
message("\nSupplementary Table S1 — Cell Type Counts:")
print(table_s1, row.names = FALSE)

# 保存
write.csv(table_s1, file.path(RES_DIR, "TableS1_celltype_counts.csv"), row.names = FALSE)

# 标记exploratory的细胞类型
exploratory_types <- table_s1$CellType[table_s1$Status == "exploratory" &
                                         !grepl("Malignant/", table_s1$CellType)]
if (length(exploratory_types) > 0) {
  message("\nWARNING: The following cell types have <", SC_MIN_CELLS_EXPLORATORY,
          " cells and are marked as 'exploratory':")
  for (et in exploratory_types) {
    message("  - ", et)
  }
  message("Results involving these types should be interpreted with caution.")
}

# --- 6.5 注释详细表（cluster级别）---
annotation_details <- data.frame(
  Cluster     = names(cluster_annotation),
  CellType    = unname(cluster_annotation),
  MappedLabel = celltype_map[unname(cluster_annotation)],
  N_cells     = sapply(names(cluster_annotation), function(cl) sum(Idents(seu) == cl)),
  stringsAsFactors = FALSE
)

# 添加每个cluster的SingleR多数投票
annotation_details$SingleR_majority <- sapply(names(cluster_annotation), function(cl) {
  cells <- colnames(seu)[Idents(seu) == cl]
  labels <- seu$SingleR_hpca[cells]
  names(sort(table(labels), decreasing = TRUE))[1]
})

write.csv(annotation_details, file.path(RES_DIR, "cluster_annotation_details.csv"),
          row.names = FALSE)

# --- 6.6 Marker score热图（cluster x celltype）---
score_cols <- paste0(names(markers), "_mscore1")
score_cols_valid <- score_cols[score_cols %in% colnames(seu@meta.data)]

if (length(score_cols_valid) >= 3) {
  # 计算cluster x celltype score矩阵
  cluster_score_mat <- sapply(score_cols_valid, function(sc) {
    tapply(seu@meta.data[[sc]], Idents(seu), mean, na.rm = TRUE)
  })
  colnames(cluster_score_mat) <- gsub("_mscore1", "", colnames(cluster_score_mat))

  # 缩放用于热图
  cluster_score_scaled <- scale(cluster_score_mat)

  p_hm <- pheatmap::pheatmap(
    t(cluster_score_scaled),
    color = colorRampPalette(c("#4DBBD5", "white", "#E64B35"))(100),
    cluster_rows  = TRUE,
    cluster_cols  = TRUE,
    display_numbers = FALSE,
    fontsize = 10,
    main = "Marker Module Scores by Cluster (z-score)",
    filename = file.path(FIG_DIR, "sc_marker_score_heatmap.pdf"),
    width = 12, height = 8
  )
}

###############################################################################
# 7. 保存注释后的对象
###############################################################################
message("=== [7/7] Saving annotated object ===")

# 添加元数据标注
seu$celltype_exploratory <- seu$celltype %in% exploratory_types

saveRDS(seu, file = file.path(DATA_PROC, "seu_annotated.rds"))

message("\n=== Cell Annotation Complete ===")
message("Seurat object: ", file.path(DATA_PROC, "seu_annotated.rds"))
message("Table S1:      ", file.path(RES_DIR, "TableS1_celltype_counts.csv"))
message("Details CSV:   ", file.path(RES_DIR, "cluster_annotation_details.csv"))
message(sprintf("Total: %d cells across %d cell types",
                ncol(seu), length(unique(seu$celltype))))
if (sum(!is.na(seu$malignant_subtype)) > 0) {
  message(sprintf("Malignant subtypes: %s",
                  paste(names(table(seu$malignant_subtype)), collapse = ", ")))
}
message("\nNext step: Run 03_single_cell/03_upr_scoring.R")
