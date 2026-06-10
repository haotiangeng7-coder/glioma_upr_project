###############################################################################
# 03_upr_scoring.R
# 单细胞水平UPR通路激活评分 — 完整版（含评审修订要求）
# 三种方法: AUCell(主) + AddModuleScore + UCell
# AUCell自动阈值 + GMM/中位数敏感性分析
#
# 输出:
#   - data/processed/seu_upr_scored.rds
#   - figures/Fig1a_upr_score_umap.pdf
#   - figures/Fig1b_three_arms_umap.pdf
#   - figures/Fig1c_upr_violin_celltype.pdf
#   - figures/Fig1d_upr_dotplot.pdf
#   - figures/Fig1e_upr_malignant_subtypes.pdf
#   - figures/Fig1f_upr_proportion.pdf
#   - figures/Fig1_composite.pdf
#   - figures/sc_upr_method_comparison.pdf
#   - figures/sc_upr_sensitivity_analysis.pdf
#   - results/upr_scores_by_celltype.csv
#   - results/upr_pairwise_wilcoxon.csv
#   - results/upr_malignant_subtype_stats.csv
#   - results/upr_sensitivity_grouping.csv
#   - results/upr_method_correlation.csv
###############################################################################

# === 加载配置和依赖 ===
source("00_setup/config.R")
load(file.path(DATA_PROC, "upr_gene_sets.RData"))

suppressPackageStartupMessages({
  library(Seurat)
  library(AUCell)
  library(UCell)
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(ComplexHeatmap)
  library(circlize)
  library(mclust)        # GMM分组
})

set.seed(SEED)

###############################################################################
# 1. 加载数据和基因集
###############################################################################
message("=== [1/9] Loading data ===")
seu <- readRDS(file.path(DATA_PROC, "seu_annotated.rds"))
message(sprintf("Loaded: %d genes x %d cells", nrow(seu), ncol(seu)))
message(sprintf("Cell types: %s", paste(unique(seu$celltype), collapse = ", ")))

# 定义UPR通路
upr_pathways <- list(
  IRE1_XBP1 = IRE1_XBP1_genes,
  PERK_ATF4 = PERK_ATF4_genes,
  ATF6      = ATF6_genes,
  UPR_broad = UPR_broad_genes
)

# 过滤基因集中不在数据中的基因
upr_pathways_valid <- lapply(upr_pathways, function(genes) {
  genes[genes %in% rownames(seu)]
})

for (pname in names(upr_pathways_valid)) {
  message(sprintf("  %s: %d/%d genes found in data",
                  pname, length(upr_pathways_valid[[pname]]), length(upr_pathways[[pname]])))
}

# 移除基因数<5的通路
upr_pathways_valid <- upr_pathways_valid[sapply(upr_pathways_valid, length) >= 5]
if (length(upr_pathways_valid) == 0) {
  stop("No UPR pathways have >=5 genes in the data. Check gene names.")
}

###############################################################################
# 2. 方法一: AUCell评分（主分析方法）
###############################################################################
message("=== [2/9] AUCell Scoring (primary method) ===")

# 获取表达矩阵
expr_mat <- tryCatch(
  GetAssayData(seu, layer = "data"),
  error = function(e) tryCatch(GetAssayData(seu, layer = "data"), error = function(e) GetAssayData(seu, slot = "data"))
)

# 构建排名矩阵
set.seed(SEED)
cells_rankings <- AUCell_buildRankings(expr_mat, plotStats = FALSE, nCores = 1)

# 计算AUC值
cells_AUC <- AUCell_calcAUC(upr_pathways_valid, cells_rankings,
                              aucMaxRank = nrow(cells_rankings) * 0.05)

# 提取AUC分数并添加到Seurat
auc_matrix <- getAUC(cells_AUC)
for (gs_name in rownames(auc_matrix)) {
  seu[[paste0("AUCell_", gs_name)]] <- auc_matrix[gs_name, colnames(seu)]
}

message("  AUCell scores added for: ", paste(rownames(auc_matrix), collapse = ", "))

# --- AUCell自动阈值 ---
message("  Computing AUCell automatic thresholds...")
aucell_thresholds <- list()

for (gs_name in rownames(auc_matrix)) {
  tryCatch({
    set.seed(SEED)
    thres <- AUCell_exploreThresholds(cells_AUC[gs_name, ],
                                       plotHist = FALSE,
                                       assign = TRUE)
    auto_thres <- thres[[gs_name]]$aucThr$selected
    aucell_thresholds[[gs_name]] <- auto_thres
    message(sprintf("    %s: auto threshold = %.4f", gs_name, auto_thres))
  }, error = function(e) {
    # 备选: 使用75th percentile
    fallback <- quantile(auc_matrix[gs_name, ], 0.75)
    aucell_thresholds[[gs_name]] <<- fallback
    message(sprintf("    %s: auto threshold failed, using Q75 = %.4f", gs_name, fallback))
  })
}

###############################################################################
# 3. 方法二: AddModuleScore (Seurat)
###############################################################################
message("=== [3/9] AddModuleScore Scoring ===")

for (pathway_name in names(upr_pathways_valid)) {
  set.seed(SEED)
  seu <- AddModuleScore(seu, features = list(upr_pathways_valid[[pathway_name]]),
                         name = paste0("ModuleScore_", pathway_name))
}

message("  AddModuleScore scores added.")

###############################################################################
# 4. 方法三: UCell评分
###############################################################################
message("=== [4/9] UCell Scoring ===")

set.seed(SEED)
seu <- AddModuleScore_UCell(seu, features = upr_pathways_valid, name = "_UCell")

message("  UCell scores added.")

###############################################################################
# 5. 主分析评分与分组（AUCell + 自动阈值）
###############################################################################
message("=== [5/9] Primary UPR scoring and grouping ===")

# 使用AUCell_UPR_broad作为主要UPR综合评分
primary_score_col <- "AUCell_UPR_broad"
if (!primary_score_col %in% colnames(seu@meta.data)) {
  # 备选：使用所有可用AUCell分数的均值
  aucell_cols <- grep("^AUCell_", colnames(seu@meta.data), value = TRUE)
  seu$UPR_score <- rowMeans(seu@meta.data[, aucell_cols, drop = FALSE], na.rm = TRUE)
  message("  Using mean of available AUCell scores as UPR_score")
} else {
  seu$UPR_score <- seu[[primary_score_col]]
  message("  Using AUCell_UPR_broad as primary UPR_score")
}

# --- 5.1 主分组方法: AUCell自动阈值 ---
if ("UPR_broad" %in% names(aucell_thresholds)) {
  auto_threshold <- aucell_thresholds[["UPR_broad"]]
} else {
  auto_threshold <- quantile(seu$UPR_score, 0.75, na.rm = TRUE)
}

seu$UPR_group_aucell <- ifelse(seu$UPR_score >= auto_threshold, "UPR-high", "UPR-low")
message(sprintf("  AUCell auto-threshold grouping: high=%d, low=%d (threshold=%.4f)",
                sum(seu$UPR_group_aucell == "UPR-high"),
                sum(seu$UPR_group_aucell == "UPR-low"),
                auto_threshold))

# --- 5.2 细胞类型内部独立分组（AUCell） ---
message("  Cell-type-internal grouping...")
seu$UPR_group_internal <- NA_character_

for (ct in unique(seu$celltype)) {
  ct_cells <- colnames(seu)[seu$celltype == ct]
  ct_scores <- seu$UPR_score[ct_cells]

  # 对每种细胞类型内部使用AUCell阈值或中位数分组
  ct_threshold <- median(ct_scores, na.rm = TRUE)
  seu$UPR_group_internal[ct_cells] <- ifelse(ct_scores >= ct_threshold,
                                               "UPR-high", "UPR-low")
}

# 主分析使用AUCell自动阈值
seu$UPR_group <- seu$UPR_group_aucell

###############################################################################
# 6. 敏感性分析分组（GMM + 中位数）
###############################################################################
message("=== [6/9] Sensitivity analysis grouping (GMM + median) ===")

# --- 6.1 GMM分组 ---
message("  Fitting GMM (2 components)...")
set.seed(SEED)
gmm_fit <- Mclust(seu$UPR_score, G = 2, verbose = FALSE)

if (!is.null(gmm_fit)) {
  seu$UPR_group_gmm <- ifelse(gmm_fit$classification == which.max(gmm_fit$parameters$mean),
                                "UPR-high", "UPR-low")
  message(sprintf("    GMM grouping: high=%d, low=%d",
                  sum(seu$UPR_group_gmm == "UPR-high"),
                  sum(seu$UPR_group_gmm == "UPR-low")))
  message(sprintf("    GMM means: %.4f, %.4f",
                  gmm_fit$parameters$mean[1], gmm_fit$parameters$mean[2]))
} else {
  message("    GMM fitting failed. Using tertile-based grouping as fallback.")
  q33 <- quantile(seu$UPR_score, 1/3, na.rm = TRUE)
  q66 <- quantile(seu$UPR_score, 2/3, na.rm = TRUE)
  seu$UPR_group_gmm <- ifelse(seu$UPR_score >= q66, "UPR-high", "UPR-low")
}

# --- 6.2 中位数分组 ---
median_threshold <- median(seu$UPR_score, na.rm = TRUE)
seu$UPR_group_median <- ifelse(seu$UPR_score >= median_threshold, "UPR-high", "UPR-low")
message(sprintf("  Median grouping: high=%d, low=%d (median=%.4f)",
                sum(seu$UPR_group_median == "UPR-high"),
                sum(seu$UPR_group_median == "UPR-low"),
                median_threshold))

# --- 6.3 分组方法一致性分析 ---
agreement_aucell_gmm <- mean(seu$UPR_group_aucell == seu$UPR_group_gmm, na.rm = TRUE)
agreement_aucell_median <- mean(seu$UPR_group_aucell == seu$UPR_group_median, na.rm = TRUE)
agreement_gmm_median <- mean(seu$UPR_group_gmm == seu$UPR_group_median, na.rm = TRUE)

sensitivity_df <- data.frame(
  Comparison = c("AUCell vs GMM", "AUCell vs Median", "GMM vs Median"),
  Agreement  = c(agreement_aucell_gmm, agreement_aucell_median, agreement_gmm_median),
  N_high_method1 = c(sum(seu$UPR_group_aucell == "UPR-high"),
                     sum(seu$UPR_group_aucell == "UPR-high"),
                     sum(seu$UPR_group_gmm == "UPR-high")),
  N_high_method2 = c(sum(seu$UPR_group_gmm == "UPR-high"),
                     sum(seu$UPR_group_median == "UPR-high"),
                     sum(seu$UPR_group_median == "UPR-high"))
)

message("\n  Grouping method agreement:")
print(sensitivity_df)
write.csv(sensitivity_df, file.path(RES_DIR, "upr_sensitivity_grouping.csv"), row.names = FALSE)

# --- 6.4 敏感性分析可视化 ---
sens_plot_data <- data.frame(
  Cell     = colnames(seu),
  CellType = seu$celltype,
  UPR_score = seu$UPR_score,
  AUCell_group = seu$UPR_group_aucell,
  GMM_group    = seu$UPR_group_gmm,
  Median_group = seu$UPR_group_median
)

p_sens1 <- ggplot(sens_plot_data, aes(x = UPR_score, fill = AUCell_group)) +
  geom_histogram(bins = 50, alpha = 0.7) +
  geom_vline(xintercept = auto_threshold, linetype = "dashed", color = "red", linewidth = 1) +
  scale_fill_manual(values = c("UPR-high" = "#E64B35", "UPR-low" = "#4DBBD5")) +
  labs(title = "AUCell Auto-threshold", x = "UPR Score", y = "Count") +
  THEME_PUBLICATION

p_sens2 <- ggplot(sens_plot_data, aes(x = UPR_score, fill = GMM_group)) +
  geom_histogram(bins = 50, alpha = 0.7) +
  scale_fill_manual(values = c("UPR-high" = "#E64B35", "UPR-low" = "#4DBBD5")) +
  labs(title = "GMM Grouping", x = "UPR Score", y = "Count") +
  THEME_PUBLICATION

p_sens3 <- ggplot(sens_plot_data, aes(x = UPR_score, fill = Median_group)) +
  geom_histogram(bins = 50, alpha = 0.7) +
  geom_vline(xintercept = median_threshold, linetype = "dashed", color = "blue", linewidth = 1) +
  scale_fill_manual(values = c("UPR-high" = "#E64B35", "UPR-low" = "#4DBBD5")) +
  labs(title = "Median Grouping", x = "UPR Score", y = "Count") +
  THEME_PUBLICATION

p_sensitivity <- (p_sens1 | p_sens2 | p_sens3) +
  plot_annotation(
    title = "Sensitivity Analysis: UPR Grouping Methods",
    subtitle = sprintf("Agreement: AUCell-GMM=%.1f%%, AUCell-Median=%.1f%%, GMM-Median=%.1f%%",
                       agreement_aucell_gmm * 100,
                       agreement_aucell_median * 100,
                       agreement_gmm_median * 100),
    theme = theme(plot.title    = element_text(size = 14, face = "bold", hjust = 0.5),
                  plot.subtitle = element_text(size = 11, hjust = 0.5))
  )

ggsave(file.path(FIG_DIR, "sc_upr_sensitivity_analysis.pdf"), p_sensitivity,
       width = 18, height = 6)

###############################################################################
# 7. 三种评分方法一致性比较
###############################################################################
message("=== [7/9] Method comparison (AUCell vs AddModuleScore vs UCell) ===")

# 提取三种方法的UPR_broad分数
method_cols <- list(
  AUCell      = "AUCell_UPR_broad",
  ModuleScore = "ModuleScore_UPR_broad1",
  UCell       = "UPR_broad_UCell"
)

method_cols_valid <- method_cols[sapply(method_cols, function(x) x %in% colnames(seu@meta.data))]

if (length(method_cols_valid) >= 2) {
  # 相关性矩阵
  method_scores <- seu@meta.data[, unlist(method_cols_valid), drop = FALSE]
  colnames(method_scores) <- names(method_cols_valid)

  cor_mat <- cor(method_scores, method = "spearman", use = "complete.obs")

  message("\n  Scoring method correlation (Spearman):")
  print(round(cor_mat, 3))

  cor_df <- as.data.frame(as.table(cor_mat))
  colnames(cor_df) <- c("Method1", "Method2", "Spearman_rho")
  cor_df <- cor_df[cor_df$Method1 != cor_df$Method2, ]
  cor_df <- cor_df[!duplicated(paste(pmin(as.character(cor_df$Method1), as.character(cor_df$Method2)),
                                      pmax(as.character(cor_df$Method1), as.character(cor_df$Method2)))), ]
  write.csv(cor_df, file.path(RES_DIR, "upr_method_correlation.csv"), row.names = FALSE)

  # 方法对比散点图
  plot_pairs <- combn(names(method_cols_valid), 2, simplify = FALSE)
  pair_plots <- list()

  for (i in seq_along(plot_pairs)) {
    m1 <- plot_pairs[[i]][1]
    m2 <- plot_pairs[[i]][2]
    rho <- cor(method_scores[[m1]], method_scores[[m2]],
               method = "spearman", use = "complete.obs")

    pair_plots[[i]] <- ggplot(method_scores, aes(x = .data[[m1]], y = .data[[m2]])) +
      geom_point(alpha = 0.1, size = 0.3) +
      geom_smooth(method = "lm", color = "#E64B35", linewidth = 0.8) +
      labs(x = m1, y = m2,
           title = sprintf("%s vs %s (rho=%.3f)", m1, m2, rho)) +
      THEME_PUBLICATION
  }

  p_method_comp <- wrap_plots(pair_plots, ncol = length(pair_plots)) +
    plot_annotation(title = "UPR Scoring Method Comparison",
                    theme = theme(plot.title = element_text(size = 14, face = "bold",
                                                             hjust = 0.5)))

  ggsave(file.path(FIG_DIR, "sc_upr_method_comparison.pdf"), p_method_comp,
         width = 6 * length(pair_plots), height = 6)
}

###############################################################################
# 8. 统计检验
###############################################################################
message("=== [8/9] Statistical Tests ===")

# --- 8.1 Kruskal-Wallis: 各细胞类型间UPR评分差异 ---
kw_test <- kruskal.test(UPR_score ~ celltype, data = seu@meta.data)
message(sprintf("Kruskal-Wallis (UPR ~ cell type): chi-sq = %.2f, df = %d, p = %.2e",
                kw_test$statistic, kw_test$parameter, kw_test$p.value))

# --- 8.2 Wilcoxon pairwise比较（BH校正）---
message("  Running pairwise Wilcoxon tests...")
celltypes <- unique(seu$celltype)
n_ct <- length(celltypes)

if (n_ct >= 2) {
  pw_results <- list()
  pair_idx <- 1

  for (i in 1:(n_ct - 1)) {
    for (j in (i + 1):n_ct) {
      ct1 <- celltypes[i]
      ct2 <- celltypes[j]
      scores1 <- seu$UPR_score[seu$celltype == ct1]
      scores2 <- seu$UPR_score[seu$celltype == ct2]

      wt <- wilcox.test(scores1, scores2)

      # 计算effect size (rank-biserial correlation)
      n1 <- length(scores1)
      n2 <- length(scores2)
      r_effect <- 1 - (2 * wt$statistic) / (n1 * n2)

      pw_results[[pair_idx]] <- data.frame(
        CellType1     = ct1,
        CellType2     = ct2,
        W_statistic   = wt$statistic,
        p_value       = wt$p.value,
        effect_size_r = r_effect,
        n1            = n1,
        n2            = n2,
        stringsAsFactors = FALSE
      )
      pair_idx <- pair_idx + 1
    }
  }

  pw_df <- do.call(rbind, pw_results)
  pw_df$p_adjusted <- p.adjust(pw_df$p_value, method = "BH")
  pw_df$significance <- ifelse(pw_df$p_adjusted < 0.001, "***",
                                ifelse(pw_df$p_adjusted < 0.01, "**",
                                       ifelse(pw_df$p_adjusted < 0.05, "*", "ns")))
  pw_df <- pw_df %>% arrange(p_adjusted)

  write.csv(pw_df, file.path(RES_DIR, "upr_pairwise_wilcoxon.csv"), row.names = FALSE)

  message(sprintf("  %d pairwise comparisons; %d significant (FDR<0.05)",
                  nrow(pw_df), sum(pw_df$p_adjusted < 0.05)))
}

# --- 8.3 各细胞类型UPR评分描述统计 ---
upr_stats <- seu@meta.data %>%
  dplyr::group_by(celltype) %>%
  dplyr::summarise(
    n         = n(),
    mean_UPR  = mean(UPR_score, na.rm = TRUE),
    median_UPR = median(UPR_score, na.rm = TRUE),
    sd_UPR    = sd(UPR_score, na.rm = TRUE),
    Q25       = quantile(UPR_score, 0.25, na.rm = TRUE),
    Q75       = quantile(UPR_score, 0.75, na.rm = TRUE),
    pct_high  = mean(UPR_group == "UPR-high", na.rm = TRUE) * 100,
    .groups   = "drop"
  ) %>%
  dplyr::arrange(desc(mean_UPR))

message("\nUPR Score by Cell Type:")
print(as.data.frame(upr_stats))

write.csv(upr_stats, file.path(RES_DIR, "upr_scores_by_celltype.csv"), row.names = FALSE)

# --- 8.4 恶性细胞亚型UPR差异 ---
if ("malignant_subtype" %in% colnames(seu@meta.data) &&
    sum(!is.na(seu$malignant_subtype)) > 0) {

  message("\n  Malignant subtype UPR analysis...")
  mal_data <- seu@meta.data %>%
    dplyr::filter(celltype == "Malignant" & !is.na(malignant_subtype))

  if (nrow(mal_data) >= 20) {
    # Kruskal-Wallis
    kw_mal <- kruskal.test(UPR_score ~ malignant_subtype, data = mal_data)
    message(sprintf("  KW (UPR ~ subtype): chi-sq=%.2f, p=%.2e",
                    kw_mal$statistic, kw_mal$p.value))

    # 描述统计
    mal_stats <- mal_data %>%
      dplyr::group_by(malignant_subtype) %>%
      dplyr::summarise(
        n          = n(),
        mean_UPR   = mean(UPR_score, na.rm = TRUE),
        median_UPR = median(UPR_score, na.rm = TRUE),
        sd_UPR     = sd(UPR_score, na.rm = TRUE),
        pct_high   = mean(UPR_group == "UPR-high", na.rm = TRUE) * 100,
        .groups    = "drop"
      ) %>%
      dplyr::arrange(desc(mean_UPR))

    message("\n  UPR by Malignant Subtype:")
    print(as.data.frame(mal_stats))

    # Pairwise Wilcoxon for malignant subtypes
    subtypes <- unique(mal_data$malignant_subtype)
    if (length(subtypes) >= 2) {
      pw_mal_list <- list()
      idx <- 1
      for (i in 1:(length(subtypes) - 1)) {
        for (j in (i + 1):length(subtypes)) {
          s1 <- mal_data$UPR_score[mal_data$malignant_subtype == subtypes[i]]
          s2 <- mal_data$UPR_score[mal_data$malignant_subtype == subtypes[j]]
          wt <- wilcox.test(s1, s2)
          pw_mal_list[[idx]] <- data.frame(
            Subtype1  = subtypes[i],
            Subtype2  = subtypes[j],
            p_value   = wt$p.value,
            n1 = length(s1), n2 = length(s2),
            stringsAsFactors = FALSE
          )
          idx <- idx + 1
        }
      }
      pw_mal_df <- do.call(rbind, pw_mal_list)
      pw_mal_df$p_adjusted <- p.adjust(pw_mal_df$p_value, method = "BH")

      mal_stats_out <- list(summary = mal_stats, kw_pvalue = kw_mal$p.value,
                             pairwise = pw_mal_df)
    } else {
      mal_stats_out <- list(summary = mal_stats)
    }

    write.csv(mal_stats, file.path(RES_DIR, "upr_malignant_subtype_stats.csv"),
              row.names = FALSE)
  }
}

###############################################################################
# 9. Figure 1 可视化
###############################################################################
message("=== [9/9] Generating Figure 1 ===")

# 构建颜色向量
all_celltypes <- unique(seu$celltype)
color_vec <- COLORS_CELLTYPE[all_celltypes]
missing_ct <- all_celltypes[is.na(color_vec)]
if (length(missing_ct) > 0) {
  extra_colors <- scales::hue_pal()(length(missing_ct))
  names(extra_colors) <- missing_ct
  color_vec[missing_ct] <- extra_colors
}
color_vec <- color_vec[!is.na(names(color_vec))]

# --- Fig 1a: UPR评分UMAP ---
p_fig1a <- FeaturePlot(seu, features = "UPR_score",
                        cols = c("lightgrey", "darkred"), order = TRUE) +
  THEME_PUBLICATION +
  ggtitle("UPR Activity Score (AUCell)")

ggsave(file.path(FIG_DIR, "Fig1a_upr_score_umap.pdf"), p_fig1a, width = 8, height = 6)

# --- Fig 1b: 三条UPR通路 ---
arm_features <- c("AUCell_IRE1_XBP1", "AUCell_PERK_ATF4", "AUCell_ATF6")
arm_features <- arm_features[arm_features %in% colnames(seu@meta.data)]

if (length(arm_features) > 0) {
  p_fig1b <- FeaturePlot(seu, features = arm_features,
                           cols = c("lightgrey", "darkred"),
                           order = TRUE, ncol = 3) &
    THEME_PUBLICATION

  ggsave(file.path(FIG_DIR, "Fig1b_three_arms_umap.pdf"), p_fig1b,
         width = 5 * length(arm_features), height = 5)
}

# --- Fig 1c: Violin plot 各细胞类型UPR评分 ---
# 使用ggpubr添加Kruskal-Wallis统计
# 使用ggplot2替代VlnPlot（避免Seurat v5 S4兼容性问题）
vln_df <- data.frame(UPR_score = seu$UPR_score, celltype = seu$celltype)
p_fig1c <- ggplot(vln_df, aes(x = celltype, y = UPR_score, fill = celltype)) +
  geom_violin(scale = "width", trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.3) +
  scale_fill_manual(values = color_vec) +
  THEME_PUBLICATION +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
  labs(x = "Cell Type", y = "UPR Score",
       title = sprintf("UPR Score Across Cell Types (KW p=%s)",
                       format.pval(kw_test$p.value, digits = 2)))

ggsave(file.path(FIG_DIR, "Fig1c_upr_violin_celltype.pdf"), p_fig1c, width = 10, height = 6)

# --- Fig 1d: Dot Plot 三条通路 x 各细胞类型 ---
aucell_cols <- grep("^AUCell_", colnames(seu@meta.data), value = TRUE)

if (length(aucell_cols) >= 3) {
  dot_data <- seu@meta.data %>%
    dplyr::select(celltype, all_of(aucell_cols)) %>%
    tidyr::pivot_longer(cols = -celltype, names_to = "Pathway", values_to = "Score") %>%
    dplyr::mutate(Pathway = gsub("AUCell_", "", Pathway))

  # 计算每个细胞类型-通路组合的均值和激活比例
  # 激活比例基于各通路的auto-threshold
  dot_summary <- dot_data %>%
    dplyr::group_by(celltype, Pathway) %>%
    dplyr::summarise(
      mean_score = mean(Score, na.rm = TRUE),
      pct_active = {
        gs <- unique(Pathway)
        thr <- if (gs %in% names(aucell_thresholds)) aucell_thresholds[[gs]]
               else median(Score, na.rm = TRUE)
        mean(Score >= thr, na.rm = TRUE) * 100
      },
      .groups = "drop"
    )

  p_fig1d <- ggplot(dot_summary, aes(x = Pathway, y = celltype)) +
    geom_point(aes(size = pct_active, color = mean_score)) +
    scale_color_gradient2(low = "#4DBBD5", mid = "white", high = "#E64B35",
                           midpoint = median(dot_summary$mean_score)) +
    scale_size_continuous(range = c(1, 8), name = "% Active") +
    labs(x = "UPR Pathway", y = "Cell Type", color = "Mean AUC") +
    THEME_PUBLICATION +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle("UPR Pathway Activity by Cell Type")

  ggsave(file.path(FIG_DIR, "Fig1d_upr_dotplot.pdf"), p_fig1d, width = 8, height = 6)
}

# --- Fig 1e: 恶性细胞亚型UPR差异 ---
if ("malignant_subtype" %in% colnames(seu@meta.data) &&
    sum(!is.na(seu$malignant_subtype)) > 0) {

  mal_cells <- subset(seu, celltype == "Malignant" & !is.na(malignant_subtype))

  if (ncol(mal_cells) >= 20) {
    # 构建pairwise比较列表
    subtypes <- unique(mal_cells$malignant_subtype)
    if (length(subtypes) >= 2) {
      comparisons <- combn(subtypes, 2, simplify = FALSE)
    } else {
      comparisons <- NULL
    }

    mal_vln_df <- data.frame(UPR_score = mal_cells$UPR_score, subtype = mal_cells$malignant_subtype)
    p_fig1e <- ggplot(mal_vln_df, aes(x = subtype, y = UPR_score, fill = subtype)) +
      geom_violin(scale = "width", trim = FALSE) +
      geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.3) +
      scale_fill_manual(values = c(
        "MES_like" = "#E64B35", "AC_like" = "#4DBBD5",
        "OPC_like" = "#00A087", "NPC_like" = "#3C5488"
      )) +
      THEME_PUBLICATION +
      labs(x = "Malignant Subtype", y = "UPR Score") +
      ggtitle("UPR Score in Malignant Subtypes")

    if (!is.null(comparisons) && length(comparisons) > 0) {
      p_fig1e <- p_fig1e +
        stat_compare_means(method = "kruskal.test", label.y.npc = 0.95) +
        stat_compare_means(comparisons = comparisons, method = "wilcox.test",
                            label = "p.signif", step.increase = 0.05)
    }

    ggsave(file.path(FIG_DIR, "Fig1e_upr_malignant_subtypes.pdf"), p_fig1e,
           width = 8, height = 6)
    rm(mal_cells)
  }
}

# --- Fig 1f: UPR-high/low 比例堆叠图 ---
prop_data <- seu@meta.data %>%
  dplyr::group_by(celltype, UPR_group) %>%
  dplyr::summarise(n = n(), .groups = "drop") %>%
  dplyr::group_by(celltype) %>%
  dplyr::mutate(pct = n / sum(n) * 100)

# 按UPR-high比例排序
ct_order <- prop_data %>%
  dplyr::filter(UPR_group == "UPR-high") %>%
  dplyr::arrange(desc(pct)) %>%
  dplyr::pull(celltype)

prop_data$celltype <- factor(prop_data$celltype, levels = ct_order)

p_fig1f <- ggplot(prop_data, aes(x = celltype, y = pct, fill = UPR_group)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("UPR-high" = "#E64B35", "UPR-low" = "#4DBBD5")) +
  labs(x = "", y = "Percentage (%)", fill = "UPR Group",
       title = "UPR-high/low Proportion by Cell Type") +
  THEME_PUBLICATION +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(FIG_DIR, "Fig1f_upr_proportion.pdf"), p_fig1f, width = 10, height = 5)

# --- 组合Figure 1 ---
message("  Assembling composite Figure 1...")

# 检查各子图是否存在
fig1_panels <- list()
if (exists("p_fig1a")) fig1_panels[["a"]] <- p_fig1a
if (exists("p_fig1b") && length(arm_features) == 3) {
  # 单独的三通路UMAP太宽，改为缩略
  fig1_panels[["b"]] <- p_fig1b
}
if (exists("p_fig1c")) fig1_panels[["c"]] <- p_fig1c
if (exists("p_fig1d")) fig1_panels[["d"]] <- p_fig1d
if (exists("p_fig1e")) fig1_panels[["e"]] <- p_fig1e
if (exists("p_fig1f")) fig1_panels[["f"]] <- p_fig1f

# 构建复合图（根据可用面板数量调整布局）
if (length(fig1_panels) >= 4) {
  # 标准6面板布局
  tryCatch({
    p_composite <- (fig1_panels[["a"]] | fig1_panels[["c"]]) /
                   (fig1_panels[["d"]] | fig1_panels[["f"]]) +
      plot_annotation(
        title = "Figure 1: UPR Pathway Activity in Single-Cell Glioma Landscape",
        tag_levels = "A",
        theme = theme(
          plot.title = element_text(size = 16, face = "bold", hjust = 0.5)
        )
      )

    ggsave(file.path(FIG_DIR, "Fig1_composite.pdf"), p_composite,
           width = 18, height = 14)
    message("  Composite Figure 1 saved.")
  }, error = function(e) {
    message("  Composite figure assembly failed: ", e$message)
    message("  Individual panels have been saved separately.")
  })
}

###############################################################################
# 保存最终对象
###############################################################################
message("=== Saving final scored object ===")

saveRDS(seu, file = file.path(DATA_PROC, "seu_upr_scored.rds"))

message("\n=== UPR Scoring Complete ===")
message("Seurat object:  ", file.path(DATA_PROC, "seu_upr_scored.rds"))
message("Results dir:    ", RES_DIR)
message("Figures dir:    ", FIG_DIR)
message(sprintf("\nScoring summary:"))
message(sprintf("  Methods: AUCell (primary), AddModuleScore, UCell"))
message(sprintf("  Primary grouping: AUCell auto-threshold (%.4f)", auto_threshold))
message(sprintf("  UPR-high: %d cells, UPR-low: %d cells",
                sum(seu$UPR_group == "UPR-high"),
                sum(seu$UPR_group == "UPR-low")))
message(sprintf("  Sensitivity: GMM + Median grouping concordance computed"))
message("\nNext step: Run 03_single_cell/04_differential_analysis.R")
