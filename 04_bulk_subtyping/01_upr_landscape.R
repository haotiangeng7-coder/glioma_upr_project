###############################################################################
# 01_upr_landscape.R
# UPR基因在TCGA胶质瘤中的表达景观
#
# 输出:
#   - figures/Fig3_upr_by_grade.pdf
#   - figures/Fig3_upr_by_idh.pdf
#   - figures/Fig3_upr_correlation_heatmap.pdf
#   - figures/Fig3_upr_by_subtype.pdf
#   - figures/Fig3_upr_forest_plot.pdf
#   - results/upr_ssgsea_scores.csv
#   - results/upr_genes_univariate_cox.csv
#   - results/upr_cox_ph_diagnostics.csv
#   - results/upr_gene_correlations.csv
#   - results/upr_pathway_correlations.csv
#   - data/processed/upr_landscape_results.RData
#
# 评审修订:
#   - Cox回归使用strata(IDH)调整, 并运行cox.zph()检验PH假设
#   - Cox结果仅用于描述性分析, 不用于聚类基因选择(消除循环论证)
#   - 相关性分析设|r|>=0.3效应量阈值
###############################################################################

source("00_setup/config.R")

library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)
library(reshape2)
library(GSVA)
library(survival)

set.seed(SEED)

# --- 加载数据 ---
load(file.path(DATA_PROC, "upr_gene_sets.RData"))
load(file.path(DATA_PROC, "tcga_glioma_expression.RData"))

# =============================================================================
# 1. UPR基因表达概览
# =============================================================================
message("=== UPR gene expression landscape ===")

# 检查UPR基因在数据集中的覆盖率
upr_in_data <- UPR_broad_genes[UPR_broad_genes %in% rownames(expr_tpm_symbol)]
upr_missing <- setdiff(UPR_broad_genes, rownames(expr_tpm_symbol))
message(sprintf("UPR genes in TCGA: %d/%d", length(upr_in_data), length(UPR_broad_genes)))
if (length(upr_missing) > 0) {
  message("Missing genes: ", paste(upr_missing, collapse = ", "))
}

# 提取UPR基因表达矩阵（log2 TPM+1）
expr_upr <- log2(expr_tpm_symbol[upr_in_data, ] + 1)

# =============================================================================
# 2. UPR基因与临床特征的关联
# =============================================================================
message("=== UPR genes vs clinical features ===")

# 匹配表达数据和临床数据
common_samples <- intersect(colnames(expr_upr), clinical_valid$barcode)
expr_upr_matched <- expr_upr[, common_samples]
clin_matched <- clinical_valid[match(common_samples, clinical_valid$barcode), ]

message(sprintf("Matched samples: %d", length(common_samples)))

# 计算每个样本的UPR总评分（mean log2 TPM+1）
upr_score_bulk <- colMeans(expr_upr_matched, na.rm = TRUE)

# --- 2.1 WHO Grade关联 ---
if ("Grade" %in% colnames(clin_matched) && sum(!is.na(clin_matched$Grade)) > 10) {
  plot_data_grade <- data.frame(
    UPR_score = upr_score_bulk,
    Grade     = clin_matched$Grade,
    stringsAsFactors = FALSE
  ) %>% dplyr::filter(!is.na(Grade))

  # Kruskal-Wallis检验
  kw_grade <- kruskal.test(UPR_score ~ Grade, data = plot_data_grade)
  message(sprintf("  Grade Kruskal-Wallis p = %.4e", kw_grade$p.value))

  p_grade <- ggplot(plot_data_grade, aes(x = Grade, y = UPR_score, fill = Grade)) +
    geom_boxplot(outlier.size = 0.5, width = 0.6) +
    geom_jitter(width = 0.15, size = 0.3, alpha = 0.2) +
    stat_compare_means(method = "kruskal.test", label = "p.format",
                       label.x.npc = 0.5, label.y.npc = 0.95) +
    labs(y = "Mean UPR Expression Score (log2 TPM+1)",
         x = "WHO Grade",
         title = "UPR Score by WHO Grade") +
    scale_fill_manual(values = c("G2" = "#FFF5EB", "G3" = "#FDAE6B", "G4" = "#D94701")) +
    THEME_PUBLICATION +
    theme(legend.position = "none")

  ggsave(file.path(FIG_DIR, "Fig3_upr_by_grade.pdf"), p_grade,
         width = 6, height = 5, useDingbats = FALSE)
}

# --- 2.2 IDH状态关联 ---
if ("IDH_status" %in% colnames(clin_matched) && sum(!is.na(clin_matched$IDH_status)) > 10) {
  plot_data_idh <- data.frame(
    UPR_score = upr_score_bulk,
    IDH       = clin_matched$IDH_status,
    stringsAsFactors = FALSE
  ) %>% dplyr::filter(!is.na(IDH))

  # Wilcoxon检验
  wt_idh <- wilcox.test(UPR_score ~ IDH, data = plot_data_idh)
  message(sprintf("  IDH Wilcoxon p = %.4e", wt_idh$p.value))

  p_idh <- ggplot(plot_data_idh, aes(x = IDH, y = UPR_score, fill = IDH)) +
    geom_boxplot(outlier.size = 0.5, width = 0.6) +
    geom_jitter(width = 0.15, size = 0.3, alpha = 0.2) +
    stat_compare_means(method = "wilcox.test", label = "p.format",
                       label.x.npc = 0.5, label.y.npc = 0.95) +
    labs(y = "Mean UPR Expression Score (log2 TPM+1)",
         x = "IDH Status",
         title = "UPR Score by IDH Status") +
    scale_fill_manual(values = c("Mutant" = "#00A087", "WT" = "#E64B35")) +
    THEME_PUBLICATION +
    theme(legend.position = "none")

  ggsave(file.path(FIG_DIR, "Fig3_upr_by_idh.pdf"), p_idh,
         width = 5, height = 5, useDingbats = FALSE)
}

# =============================================================================
# 3. UPR基因间Pearson相关性热图（按三条通路标注）
# =============================================================================
message("=== UPR gene co-expression correlation heatmap ===")

# 计算UPR基因间的Pearson相关性（基因 x 基因，在样本方向计算）
cor_mat <- cor(t(as.matrix(expr_upr_matched)), method = "pearson")

# 按通路标记基因（共有基因按赋值顺序优先归类）
pathway_labels <- rep("Other/Shared", nrow(cor_mat))
names(pathway_labels) <- rownames(cor_mat)
pathway_labels[rownames(cor_mat) %in% ATF6_genes]      <- "ATF6"
pathway_labels[rownames(cor_mat) %in% PERK_ATF4_genes] <- "PERK-ATF4"
pathway_labels[rownames(cor_mat) %in% IRE1_XBP1_genes] <- "IRE1-XBP1"

pathway_colors <- c(
  "IRE1-XBP1"    = "#E64B35",
  "PERK-ATF4"    = "#4DBBD5",
  "ATF6"         = "#00A087",
  "Other/Shared" = "grey80"
)

# 行注释和列注释
ha_row <- rowAnnotation(
  Pathway = pathway_labels,
  col = list(Pathway = pathway_colors),
  show_annotation_name = TRUE,
  annotation_name_side = "top"
)

ha_col <- HeatmapAnnotation(
  Pathway = pathway_labels,
  col = list(Pathway = pathway_colors),
  show_annotation_name = TRUE
)

# 统计|r| >= CORR_EFFECT_THRESHOLD的比例
n_significant_corr <- sum(abs(cor_mat[upper.tri(cor_mat)]) >= CORR_EFFECT_THRESHOLD)
n_total_pairs <- sum(upper.tri(cor_mat))
message(sprintf("  Correlation pairs with |r| >= %.1f: %d/%d (%.1f%%)",
                CORR_EFFECT_THRESHOLD, n_significant_corr, n_total_pairs,
                100 * n_significant_corr / n_total_pairs))

# 绘制相关性热图
pdf(file.path(FIG_DIR, "Fig3_upr_correlation_heatmap.pdf"), width = 14, height = 12)
ht <- Heatmap(
  cor_mat,
  name = "Pearson r",
  col = colorRamp2(c(-1, -0.3, 0, 0.3, 1),
                   c("#2166AC", "#92C5DE", "white", "#F4A582", "#B2182B")),
  left_annotation = ha_row,
  top_annotation = ha_col,
  show_row_names = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 6),
  clustering_method_rows = "ward.D2",
  clustering_method_columns = "ward.D2",
  column_title = sprintf("UPR Gene Co-expression (TCGA Glioma, n=%d)", length(common_samples)),
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  heatmap_legend_param = list(
    title = "Pearson r",
    at = c(-1, -0.5, 0, 0.5, 1),
    labels = c("-1.0", "-0.5", "0", "0.5", "1.0"),
    legend_height = unit(4, "cm")
  )
)
draw(ht)
dev.off()

# 保存相关性矩阵（标注|r|>=阈值的有生物学意义对）
cor_long <- reshape2::melt(cor_mat, varnames = c("Gene1", "Gene2"), value.name = "Pearson_r")
cor_long <- cor_long %>%
  dplyr::filter(as.character(Gene1) < as.character(Gene2)) %>%
  dplyr::mutate(
    biologically_meaningful = abs(Pearson_r) >= CORR_EFFECT_THRESHOLD,
    Gene1_pathway = pathway_labels[as.character(Gene1)],
    Gene2_pathway = pathway_labels[as.character(Gene2)]
  ) %>%
  dplyr::arrange(desc(abs(Pearson_r)))

write.csv(cor_long, file.path(RES_DIR, "upr_gene_correlations.csv"), row.names = FALSE)

# =============================================================================
# 4. ssGSEA计算三条通路活性评分
# =============================================================================
message("=== Computing pathway-level ssGSEA scores ===")

# 构建基因集（仅使用数据中存在的基因）
upr_pathway_list <- list(
  IRE1_XBP1 = IRE1_XBP1_genes[IRE1_XBP1_genes %in% rownames(expr_tpm_symbol)],
  PERK_ATF4 = PERK_ATF4_genes[PERK_ATF4_genes %in% rownames(expr_tpm_symbol)],
  ATF6      = ATF6_genes[ATF6_genes %in% rownames(expr_tpm_symbol)]
)

message("  Gene set sizes after filtering:")
for (nm in names(upr_pathway_list)) {
  message(sprintf("    %s: %d genes", nm, length(upr_pathway_list[[nm]])))
}

# ssGSEA (GSVA包新接口)
ssgsea_params <- ssgseaParam(
  exprData = as.matrix(log2(expr_tpm_symbol[, common_samples] + 1)),
  geneSets = upr_pathway_list,
  normalize = TRUE
)
ssgsea_scores <- gsva(ssgsea_params)

# 整理结果
ssgsea_df <- as.data.frame(t(ssgsea_scores))
ssgsea_df$barcode <- rownames(ssgsea_df)

write.csv(ssgsea_df, file.path(RES_DIR, "upr_ssgsea_scores.csv"), row.names = FALSE)

# 通路间相关性
cor_pathway <- cor(t(ssgsea_scores), method = "pearson")
message("\nUPR pathway-level correlations:")
print(round(cor_pathway, 3))

cor_pathway_df <- reshape2::melt(cor_pathway, varnames = c("Pathway1", "Pathway2"),
                                  value.name = "Pearson_r")
cor_pathway_df <- cor_pathway_df %>%
  dplyr::filter(as.character(Pathway1) < as.character(Pathway2)) %>%
  dplyr::mutate(meaningful = abs(Pearson_r) >= CORR_EFFECT_THRESHOLD)
write.csv(cor_pathway_df, file.path(RES_DIR, "upr_pathway_correlations.csv"), row.names = FALSE)

# =============================================================================
# 5. UPR与已知转录组亚型（Verhaak）的关联
# =============================================================================
message("=== UPR vs transcriptome subtypes ===")

if ("Subtype" %in% colnames(clin_matched) && sum(!is.na(clin_matched$Subtype)) > 10) {
  plot_data_sub <- data.frame(
    ssgsea_df[common_samples, c("IRE1_XBP1", "PERK_ATF4", "ATF6")],
    Subtype = clin_matched$Subtype,
    stringsAsFactors = FALSE
  ) %>%
    dplyr::filter(!is.na(Subtype)) %>%
    tidyr::pivot_longer(
      cols = c("IRE1_XBP1", "PERK_ATF4", "ATF6"),
      names_to = "Pathway", values_to = "Score"
    ) %>%
    dplyr::mutate(
      Pathway = factor(Pathway, levels = c("IRE1_XBP1", "PERK_ATF4", "ATF6"),
                        labels = c("IRE1a-XBP1", "PERK-ATF4", "ATF6"))
    )

  p_subtype <- ggplot(plot_data_sub, aes(x = Subtype, y = Score, fill = Subtype)) +
    geom_boxplot(outlier.size = 0.5, width = 0.7) +
    facet_wrap(~Pathway, scales = "free_y") +
    stat_compare_means(method = "kruskal.test", label = "p.format",
                       label.y.npc = 0.95, size = 3.5) +
    THEME_PUBLICATION +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none",
          strip.text = element_text(face = "bold", size = 11)) +
    labs(y = "ssGSEA Enrichment Score",
         x = "Verhaak Transcriptome Subtype",
         title = "UPR Pathway Activity by Transcriptome Subtype")

  ggsave(file.path(FIG_DIR, "Fig3_upr_by_subtype.pdf"), p_subtype,
         width = 12, height = 5, useDingbats = FALSE)
}

# =============================================================================
# 6. 单基因Cox回归（评审修订: 仅描述性, strata(IDH), cox.zph()）
# =============================================================================
message("=== Single-gene univariate Cox regression (descriptive only) ===")
message("  NOTE: Cox results are for descriptive landscape analysis ONLY.")
message("  Clustering gene selection uses MAD-based filtering (see 02_consensus_clustering.R).")

cox_results <- data.frame(
  gene        = character(),
  HR          = numeric(),
  HR_lower    = numeric(),
  HR_upper    = numeric(),
  pvalue      = numeric(),
  pathway     = character(),
  n_samples   = integer(),
  ph_global_p = numeric(),
  ph_passed   = logical(),
  stringsAsFactors = FALSE
)

ph_diagnostics <- list()

for (gene in upr_in_data) {
  expr_gene <- as.numeric(expr_upr_matched[gene, ])

  surv_data <- data.frame(
    time   = clin_matched$OS.time,
    status = clin_matched$OS,
    expr   = expr_gene,
    IDH    = clin_matched$IDH_status,
    stringsAsFactors = FALSE
  )
  surv_data <- surv_data[complete.cases(surv_data) & surv_data$time > 0, ]

  if (nrow(surv_data) < 50) next

  # Cox回归使用strata(IDH)（评审修订）
  tryCatch({
    fit <- coxph(Surv(time, status) ~ expr + strata(IDH), data = surv_data)
    s   <- summary(fit)

    # PH假设检验
    ph_test <- cox.zph(fit)
    ph_global_p <- ph_test$table["GLOBAL", "p"]
    ph_passed <- ph_global_p > 0.05

    if (!ph_passed) {
      message(sprintf("  WARNING: PH assumption violated for %s (p=%.4f)", gene, ph_global_p))
    }

    ph_diagnostics[[gene]] <- ph_test

    cox_results <- rbind(cox_results, data.frame(
      gene        = gene,
      HR          = s$conf.int[1, 1],
      HR_lower    = s$conf.int[1, 3],
      HR_upper    = s$conf.int[1, 4],
      pvalue      = s$coefficients[1, "Pr(>|z|)"],
      pathway     = dplyr::case_when(
        gene %in% IRE1_XBP1_genes ~ "IRE1-XBP1",
        gene %in% PERK_ATF4_genes ~ "PERK-ATF4",
        gene %in% ATF6_genes      ~ "ATF6",
        TRUE                      ~ "Other/Shared"
      ),
      n_samples   = nrow(surv_data),
      ph_global_p = ph_global_p,
      ph_passed   = ph_passed,
      stringsAsFactors = FALSE
    ))
  }, error = function(e) {
    message(sprintf("  Cox failed for %s: %s", gene, conditionMessage(e)))
  })
}

# 多重检验校正
cox_results$padj <- p.adjust(cox_results$pvalue, method = "BH")
cox_results <- cox_results %>% dplyr::arrange(pvalue)

# 保存结果
write.csv(cox_results, file.path(RES_DIR, "upr_genes_univariate_cox.csv"), row.names = FALSE)

# PH诊断摘要
ph_summary <- cox_results %>%
  dplyr::select(gene, pathway, HR, pvalue, padj, ph_global_p, ph_passed)
write.csv(ph_summary, file.path(RES_DIR, "upr_cox_ph_diagnostics.csv"), row.names = FALSE)

message(sprintf("Prognostic UPR genes (raw p<0.05): %d/%d",
                sum(cox_results$pvalue < 0.05), nrow(cox_results)))
message(sprintf("Prognostic UPR genes (FDR<0.05):   %d/%d",
                sum(cox_results$padj < 0.05), nrow(cox_results)))
message(sprintf("PH assumption violated:             %d/%d",
                sum(!cox_results$ph_passed), nrow(cox_results)))

# --- Forest plot ---
# 展示所有基因（按HR排序），标注显著性和PH假设状态
plot_genes <- cox_results %>%
  dplyr::mutate(
    sig_label = dplyr::case_when(
      padj < 0.001 ~ "***",
      padj < 0.01  ~ "**",
      padj < 0.05  ~ "*",
      TRUE         ~ "ns"
    ),
    ph_label = ifelse(ph_passed, "", " [PH!]"),
    gene_label = paste0(gene, ph_label)
  )

# 取top 40（如果基因太多）
if (nrow(plot_genes) > 40) {
  plot_genes_top <- plot_genes %>%
    dplyr::arrange(pvalue) %>%
    head(40)
} else {
  plot_genes_top <- plot_genes
}

if (nrow(plot_genes_top) > 0) {
  p_forest <- ggplot(plot_genes_top,
                     aes(x = HR, y = reorder(gene_label, HR))) +
    geom_point(aes(color = pathway, shape = ph_passed), size = 2.5) +
    geom_errorbarh(aes(xmin = HR_lower, xmax = HR_upper, color = pathway),
                   height = 0.25, linewidth = 0.5) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.5) +
    scale_color_manual(
      values = c("IRE1-XBP1" = "#E64B35", "PERK-ATF4" = "#4DBBD5",
                 "ATF6" = "#00A087", "Other/Shared" = "grey50"),
      name = "UPR Pathway"
    ) +
    scale_shape_manual(
      values = c("TRUE" = 16, "FALSE" = 1),
      labels = c("TRUE" = "Passed", "FALSE" = "Violated"),
      name = "PH Assumption"
    ) +
    labs(
      x = "Hazard Ratio (95% CI)",
      y = "",
      title = "UPR Genes: Univariate Cox Regression (strata IDH)",
      subtitle = sprintf("Top %d genes by p-value; [PH!] = PH assumption violated",
                          nrow(plot_genes_top))
    ) +
    THEME_PUBLICATION +
    theme(
      axis.text.y = element_text(size = 7),
      plot.subtitle = element_text(size = 9, color = "grey30")
    )

  fig_height <- max(5, nrow(plot_genes_top) * 0.28 + 2)
  ggsave(file.path(FIG_DIR, "Fig3_upr_forest_plot.pdf"), p_forest,
         width = 9, height = fig_height, useDingbats = FALSE)
}

# =============================================================================
# 7. 保存中间结果
# =============================================================================
save(
  expr_upr_matched, clin_matched, common_samples,
  ssgsea_scores, ssgsea_df, cox_results, ph_diagnostics,
  upr_in_data, upr_score_bulk, cor_mat, pathway_labels,
  file = file.path(DATA_PROC, "upr_landscape_results.RData")
)

message("\n=== UPR landscape analysis completed ===")
message(sprintf("  Figures saved to: %s", FIG_DIR))
message(sprintf("  Results saved to: %s", RES_DIR))
message("  NOTE: Cox results are DESCRIPTIVE ONLY - clustering uses MAD-based gene selection.")
message("Next step: Run 04_bulk_subtyping/02_consensus_clustering.R")
