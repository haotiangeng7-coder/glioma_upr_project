###############################################################################
# 02_consensus_clustering.R
# 基于UPR基因的共识聚类分子分型
#
# 输出:
#   - results/consensus_clustering/  (CCP内置PDF输出)
#   - figures/Fig3_km_K*.pdf          (各K值KM曲线)
#   - figures/Fig3_km_final.pdf       (最终KM曲线)
#   - figures/Fig3_km_idhwt.pdf       (IDH-WT亚组KM)
#   - figures/Fig3_cluster_metrics.pdf (K选择指标汇总)
#   - figures/FigS_sensitivity_*.pdf  (敏感性分析)
#   - results/upr_subtype_assignments.csv
#   - results/consensus_clustering_metrics.csv
#   - results/sensitivity_analysis_summary.csv
#   - data/processed/consensus_clustering_results.RData
#
# 评审修订（关键）:
#   - 聚类基因选择：全部UPR基因中MAD前50%（CC_MAD_PERCENTILE），不使用Cox筛选
#   - 敏感性分析：(1)全部78个UPR基因；(2)MSigDB Hallmark UPR基因集
#   - IDH分层：在IDH-WT亚组独立重复聚类
#   - K选择：CDF, Delta area, PAC, Silhouette综合评估
###############################################################################

source("00_setup/config.R")

library(ConsensusClusterPlus)
library(survival)
library(survminer)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(cluster)

set.seed(SEED)

# --- 加载数据 ---
load(file.path(DATA_PROC, "upr_gene_sets.RData"))
load(file.path(DATA_PROC, "tcga_glioma_expression.RData"))
load(file.path(DATA_PROC, "upr_landscape_results.RData"))

# =============================================================================
# 辅助函数
# =============================================================================

#' 运行共识聚类并返回结果
#' @param mat 表达矩阵 (genes x samples), 已z-score标准化
#' @param out_dir 输出目录
#' @param max_k 最大聚类数
#' @param reps 重采样次数
#' @param p_item 样本重采样比例
#' @param seed 随机种子
#' @return ConsensusClusterPlus结果列表
run_consensus_cluster <- function(mat, out_dir, max_k = CC_MAX_K,
                                  reps = CC_REPS, p_item = CC_PITEM,
                                  seed = SEED) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  cc_res <- ConsensusClusterPlus(
    d          = mat,
    maxK       = max_k,
    reps       = reps,
    pItem      = p_item,
    pFeature   = 1,
    clusterAlg = "hc",
    distance   = "pearson",
    innerLinkage  = "ward.D2",
    finalLinkage  = "ward.D2",
    seed       = seed,
    title      = out_dir,
    plot       = "pdf"
  )
  return(cc_res)
}

#' 计算CDF曲线下面积
#' @param cons_mat 共识矩阵
#' @return AUC值
calc_cdf_auc <- function(cons_mat) {
  vals <- cons_mat[lower.tri(cons_mat)]
  h <- hist(vals, breaks = seq(0, 1, by = 0.01), plot = FALSE)
  cdf <- cumsum(h$counts) / sum(h$counts)
  auc <- sum(cdf * diff(seq(0, 1, length.out = length(cdf) + 1)))
  return(auc)
}

#' 计算K选择综合指标（PAC, CDF AUC, Delta Area, Silhouette）
#' @param cc_res ConsensusClusterPlus结果
#' @param mat 原始标准化矩阵
#' @param max_k 最大K值
#' @return data.frame 各K值指标
calc_cluster_metrics <- function(cc_res, mat, max_k = CC_MAX_K) {
  metrics <- data.frame(
    K         = 2:max_k,
    PAC       = numeric(max_k - 1),
    CDF_AUC   = numeric(max_k - 1),
    Delta_AUC = numeric(max_k - 1),
    Silhouette_mean = numeric(max_k - 1),
    stringsAsFactors = FALSE
  )

  prev_auc <- 0
  for (i in seq_along(2:max_k)) {
    k <- i + 1
    cons_mat <- cc_res[[k]]$consensusMatrix

    # PAC: proportion of ambiguous clustering (0.1 < c_ij < 0.9)
    metrics$PAC[i] <- sum(cons_mat > 0.1 & cons_mat < 0.9) / (nrow(cons_mat)^2)

    # CDF AUC
    auc <- calc_cdf_auc(cons_mat)
    metrics$CDF_AUC[i] <- auc

    # Delta Area (relative change)
    if (i == 1) {
      metrics$Delta_AUC[i] <- auc
    } else {
      metrics$Delta_AUC[i] <- (auc - prev_auc) / prev_auc
    }
    prev_auc <- auc

    # Silhouette
    clusters <- cc_res[[k]]$consensusClass
    # 使用 1 - 共识矩阵作为距离
    dist_mat <- as.dist(1 - cons_mat)
    sil <- silhouette(clusters, dist_mat)
    metrics$Silhouette_mean[i] <- mean(sil[, "sil_width"])
  }

  return(metrics)
}

#' 对指定样本子集运行KM生存分析
#' @param cluster_assignments 命名向量 (barcode -> cluster)
#' @param clin 临床数据data.frame
#' @param pdf_path PDF输出路径
#' @param title 图标题
#' @param palette 颜色向量
#' @return log-rank p值
run_km_analysis <- function(cluster_assignments, clin, pdf_path, title,
                            palette = NULL) {
  surv_df <- data.frame(
    barcode = names(cluster_assignments),
    cluster = cluster_assignments,
    stringsAsFactors = FALSE
  )
  surv_df <- merge(surv_df, clin, by = "barcode")
  surv_df <- surv_df[!is.na(surv_df$OS.time) & surv_df$OS.time > 0, ]

  if (nrow(surv_df) < 20) {
    message("  Too few samples for KM analysis: ", nrow(surv_df))
    return(NA)
  }

  n_groups <- length(unique(surv_df$cluster))
  if (is.null(palette)) {
    palette <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F", "#8491B4")[1:n_groups]
  }

  fit_km <- survfit(Surv(OS.time / 30, OS) ~ cluster, data = surv_df)

  # log-rank检验
  lr_test <- survdiff(Surv(OS.time / 30, OS) ~ cluster, data = surv_df)
  lr_pval <- 1 - pchisq(lr_test$chisq, df = n_groups - 1)

  p_km <- ggsurvplot(
    fit_km,
    data     = surv_df,
    pval     = TRUE,
    pval.method = TRUE,
    risk.table  = TRUE,
    palette     = palette,
    xlab        = "Time (months)",
    title       = title,
    ggtheme     = THEME_PUBLICATION,
    risk.table.height = 0.25,
    conf.int    = TRUE,
    conf.int.alpha = 0.15
  )

  pdf(pdf_path, width = 8, height = 7)
  print(p_km)
  dev.off()

  return(lr_pval)
}

# =============================================================================
# 1. 聚类基因选择（评审修订: MAD筛选, 消除循环论证）
# =============================================================================
message("=== Selecting genes for clustering (MAD-based, reviewer revision) ===")
message("  NOTE: Gene selection is based on expression variability (MAD), NOT Cox regression.")
message("  This eliminates circular reasoning of using outcome to select clustering features.")

# 获取数据中存在的UPR基因
upr_in_data <- UPR_broad_genes[UPR_broad_genes %in% rownames(expr_tpm_symbol)]
expr_upr_all <- log2(expr_tpm_symbol[upr_in_data, common_samples] + 1)

message(sprintf("  Total UPR genes in data: %d", length(upr_in_data)))

# MAD (Median Absolute Deviation) 筛选
gene_mad <- apply(expr_upr_all, 1, mad)
mad_threshold <- quantile(gene_mad, 1 - CC_MAD_PERCENTILE)  # 前50%
cluster_genes <- names(gene_mad[gene_mad >= mad_threshold])

message(sprintf("  MAD percentile threshold: top %.0f%%", CC_MAD_PERCENTILE * 100))
message(sprintf("  MAD cutoff value: %.4f", mad_threshold))
message(sprintf("  Clustering genes selected: %d", length(cluster_genes)))

# 准备聚类矩阵
cluster_mat <- as.matrix(expr_upr_all[cluster_genes, ])
cluster_mat_scaled <- t(scale(t(cluster_mat)))

# =============================================================================
# 2. 主聚类分析：MAD筛选后的UPR基因
# =============================================================================
message("=== Running primary consensus clustering ===")

cc_dir <- file.path(RES_DIR, "consensus_clustering")
cc_results <- run_consensus_cluster(cluster_mat_scaled, cc_dir)

# ICL
icl <- calcICL(cc_results, title = cc_dir, plot = "pdf")

# 计算K选择综合指标
metrics <- calc_cluster_metrics(cc_results, cluster_mat_scaled)
message("\n  K selection metrics (primary analysis):")
print(metrics)

write.csv(metrics, file.path(RES_DIR, "consensus_clustering_metrics.csv"), row.names = FALSE)

# 选择最优K
optimal_k_pac <- metrics$K[which.min(metrics$PAC)]
optimal_k_sil <- metrics$K[which.max(metrics$Silhouette_mean)]
message(sprintf("  Optimal K by PAC: %d", optimal_k_pac))
message(sprintf("  Optimal K by Silhouette: %d", optimal_k_sil))

# --- K选择指标可视化 ---
metrics_long <- tidyr::pivot_longer(metrics,
  cols = c("PAC", "CDF_AUC", "Delta_AUC", "Silhouette_mean"),
  names_to = "Metric", values_to = "Value"
)

p_metrics <- ggplot(metrics_long, aes(x = K, y = Value)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.5) +
  facet_wrap(~Metric, scales = "free_y", nrow = 2) +
  scale_x_continuous(breaks = 2:CC_MAX_K) +
  labs(title = "Consensus Clustering: K Selection Metrics",
       x = "Number of Clusters (K)", y = "Value") +
  THEME_PUBLICATION +
  theme(strip.text = element_text(face = "bold"))

ggsave(file.path(FIG_DIR, "Fig3_cluster_metrics.pdf"), p_metrics,
       width = 8, height = 6, useDingbats = FALSE)

# =============================================================================
# 3. 各K值KM生存分析（全队列）
# =============================================================================
message("=== KM survival analysis for each K (full cohort) ===")

km_pvalues <- data.frame(K = integer(), logrank_p = numeric(), stringsAsFactors = FALSE)

for (k in 2:CC_MAX_K) {
  clusters <- cc_results[[k]]$consensusClass
  names(clusters) <- colnames(cluster_mat_scaled)
  cluster_named <- paste0("C", clusters)
  names(cluster_named) <- names(clusters)

  pval <- run_km_analysis(
    cluster_named, clin_matched,
    pdf_path = file.path(FIG_DIR, sprintf("Fig3_km_K%d.pdf", k)),
    title = sprintf("UPR Molecular Subtypes (K=%d, Full Cohort)", k)
  )
  km_pvalues <- rbind(km_pvalues, data.frame(K = k, logrank_p = pval))
  message(sprintf("  K=%d: log-rank p = %.4e", k, pval))
}

# =============================================================================
# 4. 选择最终K值并命名亚型
# =============================================================================
# K selection uses PRE-SPECIFIED STRUCTURAL metrics only; survival/outcome is NOT used
# to choose K (avoids circular, outcome-informed selection). Average silhouette width is
# maximized at K=2 (0.905 vs 0.859 for K=3) and the consensus CDF delta-area shows its
# largest increment at K=2. PAC marginally favors K=3 (0.230 vs 0.248) and is reported
# only as a sensitivity analysis (see k2_vs_k3_decision.R). K=2 is the primary classification.
final_k <- optimal_k_sil   # = 2 (silhouette-optimal; structural, outcome-independent)
message(sprintf("\n=== Using K=%d as final classification ===", final_k))

final_clusters <- cc_results[[final_k]]$consensusClass
names(final_clusters) <- colnames(cluster_mat_scaled)

cluster_df <- data.frame(
  barcode     = names(final_clusters),
  cluster_num = final_clusters,
  stringsAsFactors = FALSE
)

# 添加UPR均值信息（仅用于描述，不用于命名）
cluster_df <- merge(cluster_df, data.frame(
  barcode  = common_samples,
  UPR_mean = colMeans(expr_upr_all[, common_samples], na.rm = TRUE)
), by = "barcode")

# 评审修订B2: 按预后（中位OS）命名亚型，而非按UPR评分
# 原命名按UPR表达排序导致"UPR-quiescent"预后最差，名称暗示良性，具有误导性
# 新策略: 合并临床生存数据，按中位OS从好到差命名
cluster_clin <- merge(cluster_df, clin_matched[, c("barcode", "OS.time", "OS")],
                       by = "barcode", all.x = TRUE)
cluster_median_os <- tapply(cluster_clin$OS.time[!is.na(cluster_clin$OS.time)],
                             cluster_clin$cluster_num[!is.na(cluster_clin$OS.time)],
                             median)
# 按中位OS降序排列（预后最好的排第一）
rank_order_prognosis <- order(cluster_median_os, decreasing = TRUE)

# Name by prognosis position: BEST=UPR-favorable, WORST=UPR-high-risk.
# At K=2 this must be favorable + high-risk (NOT favorable + intermediate).
subtype_names <- if (final_k == 2) {
  c("UPR-favorable", "UPR-high-risk")
} else if (final_k == 3) {
  c("UPR-favorable", "UPR-intermediate", "UPR-high-risk")
} else {
  c("UPR-favorable", paste0("UPR-intermediate", seq_len(final_k - 2)), "UPR-high-risk")
}
name_map <- setNames(subtype_names, rank_order_prognosis)
cluster_df$UPR_subtype <- name_map[as.character(cluster_df$cluster_num)]

message("\n  Subtype naming based on prognosis (median OS):")
for (i in seq_along(rank_order_prognosis)) {
  cnum <- rank_order_prognosis[i]
  message(sprintf("    Cluster %d -> %s (median OS = %.0f days, mean UPR = %.3f)",
                  cnum, subtype_names[i], cluster_median_os[cnum],
                  mean(cluster_df$UPR_mean[cluster_df$cluster_num == cnum])))
}

message("\nSubtype composition (final):")
print(table(cluster_df$UPR_subtype))

# =============================================================================
# 5. 最终生存分析（全队列，命名后亚型）
# =============================================================================
message("=== Final KM survival analysis (full cohort) ===")

surv_final <- merge(cluster_df, clin_matched, by = "barcode")
surv_final <- surv_final[!is.na(surv_final$OS.time) & surv_final$OS.time > 0, ]

final_named <- setNames(surv_final$UPR_subtype, surv_final$barcode)
avail_subtypes <- intersect(names(COLORS_SUBTYPE), unique(surv_final$UPR_subtype))
palette_final <- COLORS_SUBTYPE[avail_subtypes]

pval_final <- run_km_analysis(
  final_named, clin_matched,
  pdf_path = file.path(FIG_DIR, "Fig3_km_final.pdf"),
  title = sprintf("Overall Survival by UPR Subtype (K=%d, Full Cohort)", final_k),
  palette = palette_final
)
message(sprintf("  Final KM log-rank p = %.4e", pval_final))

# =============================================================================
# 6. IDH-WT亚组独立聚类分析
# =============================================================================
message("=== IDH-WT stratified analysis ===")

if ("IDH_status" %in% colnames(clin_matched)) {
  # 获取IDH-WT样本
  idhwt_samples <- clin_matched$barcode[clin_matched$IDH_status == "WT" &
                                          !is.na(clin_matched$IDH_status)]
  idhwt_samples <- intersect(idhwt_samples, common_samples)
  message(sprintf("  IDH-WT samples: %d", length(idhwt_samples)))

  if (length(idhwt_samples) >= 30) {
    # 在IDH-WT中重新MAD筛选
    expr_upr_idhwt <- as.matrix(expr_upr_all[, idhwt_samples])
    gene_mad_idhwt <- apply(expr_upr_idhwt, 1, mad)
    mad_thresh_idhwt <- quantile(gene_mad_idhwt, 1 - CC_MAD_PERCENTILE)
    cluster_genes_idhwt <- names(gene_mad_idhwt[gene_mad_idhwt >= mad_thresh_idhwt])

    cluster_mat_idhwt <- as.matrix(expr_upr_idhwt[cluster_genes_idhwt, ])
    cluster_mat_idhwt_scaled <- t(scale(t(cluster_mat_idhwt)))

    # 运行共识聚类
    cc_dir_idhwt <- file.path(RES_DIR, "consensus_clustering_idhwt")
    cc_results_idhwt <- run_consensus_cluster(cluster_mat_idhwt_scaled, cc_dir_idhwt)

    # 指标计算
    metrics_idhwt <- calc_cluster_metrics(cc_results_idhwt, cluster_mat_idhwt_scaled)
    message("\n  IDH-WT K selection metrics:")
    print(metrics_idhwt)

    # 选择最优K
    optimal_k_idhwt <- metrics_idhwt$K[which.min(metrics_idhwt$PAC)]
    message(sprintf("  IDH-WT optimal K (PAC): %d", optimal_k_idhwt))

    # KM分析
    for (k in 2:min(4, CC_MAX_K)) {
      clusters_idhwt <- cc_results_idhwt[[k]]$consensusClass
      names(clusters_idhwt) <- colnames(cluster_mat_idhwt_scaled)
      cluster_named_idhwt <- paste0("C", clusters_idhwt)
      names(cluster_named_idhwt) <- names(clusters_idhwt)

      pval_idhwt <- run_km_analysis(
        cluster_named_idhwt, clin_matched,
        pdf_path = file.path(FIG_DIR, sprintf("Fig3_km_idhwt_K%d.pdf", k)),
        title = sprintf("IDH-WT Subgroup: UPR Subtypes (K=%d, n=%d)", k, length(idhwt_samples))
      )
      message(sprintf("  IDH-WT K=%d: log-rank p = %.4e", k, pval_idhwt))
    }

    # IDH-WT最终结果使用与全队列相同的K
    if (final_k <= CC_MAX_K) {
      clusters_idhwt_final <- cc_results_idhwt[[min(final_k, CC_MAX_K)]]$consensusClass
      names(clusters_idhwt_final) <- colnames(cluster_mat_idhwt_scaled)

      idhwt_cluster_df <- data.frame(
        barcode     = names(clusters_idhwt_final),
        idhwt_cluster_num = clusters_idhwt_final,
        stringsAsFactors = FALSE
      )

      # 评审修订B2: IDH-WT亚组同样按预后命名
      idhwt_upr <- colMeans(expr_upr_idhwt[, names(clusters_idhwt_final)], na.rm = TRUE)
      idhwt_cluster_df$UPR_mean_idhwt <- idhwt_upr[idhwt_cluster_df$barcode]
      idhwt_clin <- merge(idhwt_cluster_df, clin_matched[, c("barcode", "OS.time", "OS")],
                           by = "barcode", all.x = TRUE)
      idhwt_median_os <- tapply(idhwt_clin$OS.time[!is.na(idhwt_clin$OS.time)],
                                 idhwt_clin$idhwt_cluster_num[!is.na(idhwt_clin$OS.time)],
                                 median)
      idhwt_rank <- order(idhwt_median_os, decreasing = TRUE)
      idhwt_names <- if (final_k == 2) c("UPR-favorable", "UPR-high-risk") else
                     if (final_k == 3) c("UPR-favorable", "UPR-intermediate", "UPR-high-risk") else
                     c("UPR-favorable", paste0("UPR-intermediate", seq_len(final_k - 2)), "UPR-high-risk")
      idhwt_name_map <- setNames(idhwt_names, idhwt_rank)
      idhwt_cluster_df$UPR_subtype_idhwt <- idhwt_name_map[as.character(idhwt_cluster_df$idhwt_cluster_num)]

      message("\nIDH-WT subtype composition:")
      print(table(idhwt_cluster_df$UPR_subtype_idhwt))
    }
  } else {
    message("  Too few IDH-WT samples for independent clustering.")
    cc_results_idhwt <- NULL
    idhwt_cluster_df <- NULL
  }
} else {
  message("  IDH_status not available in clinical data.")
  cc_results_idhwt <- NULL
  idhwt_cluster_df <- NULL
}

# =============================================================================
# 7. 敏感性分析（评审修订）
# =============================================================================
message("=== Sensitivity analyses (reviewer revision) ===")

sensitivity_results <- data.frame(
  Analysis         = character(),
  N_genes          = integer(),
  Optimal_K_PAC    = integer(),
  Optimal_K_Sil    = integer(),
  LogRank_p_K2     = numeric(),
  LogRank_p_K3     = numeric(),
  ARI_vs_primary   = numeric(),
  stringsAsFactors = FALSE
)

# --- 7.1 全部UPR基因聚类 ---
message("  --- Sensitivity 1: All UPR genes (no MAD filtering) ---")

cluster_mat_all <- as.matrix(expr_upr_all)
cluster_mat_all_scaled <- t(scale(t(cluster_mat_all)))

cc_dir_all <- file.path(RES_DIR, "consensus_clustering_all_genes")
cc_results_all <- run_consensus_cluster(cluster_mat_all_scaled, cc_dir_all)

metrics_all <- calc_cluster_metrics(cc_results_all, cluster_mat_all_scaled)

# KM for K=2,3
pvals_all <- numeric(2)
for (ki in 1:2) {
  k <- ki + 1
  cls <- cc_results_all[[k]]$consensusClass
  names(cls) <- colnames(cluster_mat_all_scaled)
  cls_named <- paste0("C", cls)
  names(cls_named) <- names(cls)
  pvals_all[ki] <- run_km_analysis(
    cls_named, clin_matched,
    pdf_path = file.path(FIG_DIR, sprintf("FigS_sensitivity_all_genes_K%d.pdf", k)),
    title = sprintf("Sensitivity: All %d UPR Genes (K=%d)", nrow(cluster_mat_all), k)
  )
}

# Adjusted Rand Index vs primary
if (requireNamespace("mclust", quietly = TRUE)) {
  primary_cls <- cc_results[[final_k]]$consensusClass
  all_cls <- cc_results_all[[final_k]]$consensusClass
  shared <- intersect(names(primary_cls), names(all_cls))
  ari_all <- mclust::adjustedRandIndex(primary_cls[shared], all_cls[shared])
} else {
  ari_all <- NA
}

sensitivity_results <- rbind(sensitivity_results, data.frame(
  Analysis       = "All UPR genes (no MAD filter)",
  N_genes        = nrow(cluster_mat_all),
  Optimal_K_PAC  = metrics_all$K[which.min(metrics_all$PAC)],
  Optimal_K_Sil  = metrics_all$K[which.max(metrics_all$Silhouette_mean)],
  LogRank_p_K2   = pvals_all[1],
  LogRank_p_K3   = pvals_all[2],
  ARI_vs_primary = ari_all,
  stringsAsFactors = FALSE
))

# --- 7.2 MSigDB Hallmark UPR基因集独立聚类 ---
message("  --- Sensitivity 2: MSigDB Hallmark UPR gene set ---")

hallmark_genes <- upr_gene_list$Hallmark_UPR
hallmark_in_data <- hallmark_genes[hallmark_genes %in% rownames(expr_tpm_symbol)]
message(sprintf("    Hallmark UPR genes in data: %d", length(hallmark_in_data)))

if (length(hallmark_in_data) >= 10) {
  expr_hallmark <- log2(expr_tpm_symbol[hallmark_in_data, common_samples] + 1)
  cluster_mat_hallmark <- as.matrix(expr_hallmark)
  cluster_mat_hallmark_scaled <- t(scale(t(cluster_mat_hallmark)))

  cc_dir_hallmark <- file.path(RES_DIR, "consensus_clustering_hallmark")
  cc_results_hallmark <- run_consensus_cluster(cluster_mat_hallmark_scaled, cc_dir_hallmark)

  metrics_hallmark <- calc_cluster_metrics(cc_results_hallmark, cluster_mat_hallmark_scaled)

  # KM for K=2,3
  pvals_hallmark <- numeric(2)
  for (ki in 1:2) {
    k <- ki + 1
    cls <- cc_results_hallmark[[k]]$consensusClass
    names(cls) <- colnames(cluster_mat_hallmark_scaled)
    cls_named <- paste0("C", cls)
    names(cls_named) <- names(cls)
    pvals_hallmark[ki] <- run_km_analysis(
      cls_named, clin_matched,
      pdf_path = file.path(FIG_DIR, sprintf("FigS_sensitivity_hallmark_K%d.pdf", k)),
      title = sprintf("Sensitivity: Hallmark UPR (%d genes, K=%d)",
                       length(hallmark_in_data), k)
    )
  }

  # ARI vs primary
  if (requireNamespace("mclust", quietly = TRUE)) {
    hallmark_cls <- cc_results_hallmark[[final_k]]$consensusClass
    shared_h <- intersect(names(primary_cls), names(hallmark_cls))
    ari_hallmark <- mclust::adjustedRandIndex(primary_cls[shared_h], hallmark_cls[shared_h])
  } else {
    ari_hallmark <- NA
  }

  sensitivity_results <- rbind(sensitivity_results, data.frame(
    Analysis       = "MSigDB Hallmark UPR",
    N_genes        = length(hallmark_in_data),
    Optimal_K_PAC  = metrics_hallmark$K[which.min(metrics_hallmark$PAC)],
    Optimal_K_Sil  = metrics_hallmark$K[which.max(metrics_hallmark$Silhouette_mean)],
    LogRank_p_K2   = pvals_hallmark[1],
    LogRank_p_K3   = pvals_hallmark[2],
    ARI_vs_primary = ari_hallmark,
    stringsAsFactors = FALSE
  ))
} else {
  message("    Too few Hallmark UPR genes in data for clustering.")
}

# 保存敏感性分析结果
write.csv(sensitivity_results, file.path(RES_DIR, "sensitivity_analysis_summary.csv"),
          row.names = FALSE)
message("\nSensitivity analysis summary:")
print(sensitivity_results)

# =============================================================================
# 8. 保存所有结果
# =============================================================================
message("=== Saving results ===")

save(
  cc_results, cluster_df, surv_final, final_k, cluster_genes,
  metrics, km_pvalues,
  cc_results_idhwt, idhwt_cluster_df,
  cc_results_all, cc_results_hallmark,
  sensitivity_results,
  file = file.path(DATA_PROC, "consensus_clustering_results.RData")
)

write.csv(
  cluster_df[, c("barcode", "cluster_num", "UPR_subtype", "UPR_mean")],
  file.path(RES_DIR, "upr_subtype_assignments.csv"),
  row.names = FALSE
)

message("\n=== Consensus clustering completed ===")
message(sprintf("  Final K = %d", final_k))
message(sprintf("  Clustering genes: %d (MAD top %.0f%%)", length(cluster_genes),
                CC_MAD_PERCENTILE * 100))
message(sprintf("  Sensitivity analyses: %d completed", nrow(sensitivity_results)))
message("Next step: Run 04_bulk_subtyping/03_clinical_characterization.R")
