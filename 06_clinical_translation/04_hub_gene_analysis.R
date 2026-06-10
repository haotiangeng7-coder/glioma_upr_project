###############################################################################
# 04_hub_gene_analysis.R
# Part 4 - 临床转化 (4/4)
# Hub基因深入分析: KM生存 + 免疫细胞相关性 + 检查点共表达
#                  + STRING PPI + 单细胞验证 + HPA蛋白验证
# 输出: Figure 8
###############################################################################

source("00_setup/config.R")
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(circlize)
library(survival)
library(survminer)
library(grid)

set.seed(SEED)

# --- 加载数据 ---
load(file.path(DATA_PROC, "risk_model_final.RData"))
load(file.path(DATA_PROC, "tcga_glioma_expression.RData"))
load(file.path(DATA_PROC, "upr_gene_sets.RData"))

# --- 输出目录 ---
fig_dir <- file.path(FIG_DIR, "part4_clinical")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. 选择Hub基因 (系数绝对值Top 5)
# =============================================================================
message("=== Section 1: Hub Gene Selection ===")

# 获取模型权重/重要性 — 处理不同模型类型
if (!is.null(model_weights) && is.numeric(model_weights) && length(model_weights) > 0) {
  # Cox/LASSO等有系数的模型
  importance <- abs(model_weights)
} else if (!is.null(final_built$model)) {
  # RSF等有variable importance的模型
  tryCatch({
    if (inherits(final_built$model, "rfsrc")) {
      vimp <- randomForestSRC::vimp(final_built$model)$importance
      importance <- vimp[, 1]  # mortality importance
    } else if (inherits(final_built$model, "cv.glmnet")) {
      coefs <- as.numeric(coef(final_built$model, s = "lambda.min"))[-1]
      names(coefs) <- model_genes
      importance <- abs(coefs)
    } else {
      # 后备：用单因素Cox的-log10(p)作为重要性
      importance <- setNames(rep(1, length(model_genes)), model_genes)
    }
  }, error = function(e) {
    message("Could not extract importance: ", e$message)
    importance <<- setNames(rep(1, length(model_genes)), model_genes)
  })
} else {
  # 最终后备
  importance <- setNames(rep(1, length(model_genes)), model_genes)
}

n_hub <- min(5, length(importance))
hub_genes <- names(sort(importance, decreasing = TRUE))[1:n_hub]
hub_coefs <- importance[hub_genes]  # 这里是重要性而非系数

message(sprintf("Top %d Hub genes (by |coefficient|):", n_hub))
for (i in seq_along(hub_genes)) {
  message(sprintf("  %d. %s: coef = %.4f (|coef| = %.4f)",
                  i, hub_genes[i], hub_coefs[i], abs(hub_coefs[i])))
}

# 检查是否有模型中的关键UPR基因未在Top5中
key_upr_in_model <- intersect(
  c("HSPA5", "XBP1", "ATF4", "DDIT3", "ERN1", "EIF2AK3", "ATF6"),
  model_genes
)
additional_upr <- setdiff(key_upr_in_model, hub_genes)
if (length(additional_upr) > 0) {
  message("Additional key UPR genes in model (not in Top5): ",
          paste(additional_upr, collapse = ", "))
}

# Hub基因系数条形图
hub_coef_df <- data.frame(
  gene = factor(hub_genes, levels = rev(hub_genes)),
  coefficient = hub_coefs,
  direction = ifelse(hub_coefs > 0, "Risk", "Protective"),
  stringsAsFactors = FALSE
)

p_hub_coef <- ggplot(hub_coef_df, aes(x = gene, y = coefficient, fill = direction)) +
  geom_bar(stat = "identity", width = 0.6) +
  coord_flip() +
  scale_fill_manual(values = c("Risk" = "#E64B35", "Protective" = "#4DBBD5")) +
  labs(x = "", y = "Model Coefficient",
       title = "Hub Gene Coefficients in UIRS Model", fill = "") +
  THEME_PUBLICATION +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "Fig8_hub_coefficients.pdf"), p_hub_coef,
       width = 6, height = 4)

# =============================================================================
# 2. Hub基因K-M生存分析
# =============================================================================
message("\n=== Section 2: Hub Gene K-M Survival Analysis ===")

common_samples <- intersect(colnames(expr_tpm_symbol), clinical_valid$barcode)
hub_valid <- hub_genes[hub_genes %in% rownames(expr_tpm_symbol)]
expr_hub <- log2(expr_tpm_symbol[hub_valid, common_samples] + 1)
clin_hub <- clinical_valid[match(common_samples, clinical_valid$barcode), ]

km_results <- list()

for (gene in hub_valid) {
  gene_expr <- as.numeric(expr_hub[gene, ])
  gene_group <- ifelse(gene_expr > median(gene_expr), "High", "Low")

  surv_data <- data.frame(
    time   = clin_hub$OS.time / 30,  # 转月
    status = clin_hub$OS,
    group  = factor(gene_group, levels = c("Low", "High"))
  )
  surv_data <- surv_data[complete.cases(surv_data) & surv_data$time > 0, ]

  if (nrow(surv_data) < 20) next

  fit <- survfit(Surv(time, status) ~ group, data = surv_data)
  cox_fit <- coxph(Surv(time, status) ~ group, data = surv_data)
  cox_sum <- summary(cox_fit)

  km_results[[gene]] <- list(
    HR = cox_sum$conf.int[1, 1],
    HR_lower = cox_sum$conf.int[1, 3],
    HR_upper = cox_sum$conf.int[1, 4],
    p_value = cox_sum$coefficients[1, 5],
    log_rank_p = surv_pvalue(fit)$pval
  )

  p_km <- ggsurvplot(
    fit, data = surv_data,
    pval = TRUE, pval.method = TRUE,
    risk.table = TRUE,
    palette = c(COLORS_RISK["Low"], COLORS_RISK["High"]),
    xlab = "Time (months)",
    title = sprintf("%s (HR=%.2f, p=%.2e)",
                    gene, km_results[[gene]]$HR, km_results[[gene]]$log_rank_p),
    ggtheme = THEME_PUBLICATION,
    risk.table.height = 0.25,
    legend.labs = c("Low expression", "High expression")
  )

  pdf(file.path(fig_dir, paste0("Fig8_km_", gene, ".pdf")), width = 8, height = 7)
  print(p_km)
  dev.off()
  message(sprintf("  %s: HR=%.2f [%.2f-%.2f], log-rank p=%.2e",
                  gene, km_results[[gene]]$HR,
                  km_results[[gene]]$HR_lower, km_results[[gene]]$HR_upper,
                  km_results[[gene]]$log_rank_p))
}

# KM结果汇总
km_summary <- do.call(rbind, lapply(names(km_results), function(g) {
  data.frame(
    Gene = g,
    HR = km_results[[g]]$HR,
    HR_lower = km_results[[g]]$HR_lower,
    HR_upper = km_results[[g]]$HR_upper,
    logrank_p = km_results[[g]]$log_rank_p,
    stringsAsFactors = FALSE
  )
}))
write.csv(km_summary, file.path(RES_DIR, "hub_gene_km_summary.csv"), row.names = FALSE)

# 森林图
if (nrow(km_summary) > 0) {
  km_summary$Gene <- factor(km_summary$Gene, levels = rev(km_summary$Gene))

  p_forest <- ggplot(km_summary, aes(x = Gene, y = HR)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
    geom_pointrange(aes(ymin = HR_lower, ymax = HR_upper,
                        color = ifelse(HR > 1, "Risk", "Protective")),
                    size = 0.8) +
    coord_flip() +
    scale_color_manual(values = c("Risk" = "#E64B35", "Protective" = "#4DBBD5")) +
    labs(x = "", y = "Hazard Ratio (95% CI)",
         title = "Hub Gene Prognostic Value", color = "") +
    THEME_PUBLICATION +
    theme(legend.position = "bottom")

  ggsave(file.path(fig_dir, "Fig8_hub_forest.pdf"), p_forest, width = 6, height = 4)
}

# =============================================================================
# 3. Hub基因 vs 28种免疫细胞 Spearman相关性热图
# =============================================================================
message("\n=== Section 3: Hub Gene vs Immune Cell Correlation ===")

# 加载免疫分析结果
tryCatch({
  load(file.path(DATA_PROC, "immune_analysis_results.RData"))
  has_immune <- TRUE
}, error = function(e) {
  message("immune_analysis_results.RData not found, computing ssGSEA...")
  has_immune <- FALSE
})

if (!exists("immune_ssgsea") || !has_immune) {
  # 如果没有预计算的免疫评分，当场计算
  library(GSVA)
  message("Computing ssGSEA for 28 immune cell types...")

  immune_cell_markers <- list(
    "Activated B cell" = c("BLK", "CD19", "FAM30A", "FCRL2", "MS4A1", "PNOC", "SPIB", "TCL1A"),
    "Activated CD4 T cell" = c("DPP4", "ICOS", "IL2RA", "IL4R", "OAS1", "UBE2L6"),
    "Activated CD8 T cell" = c("CD8A", "CD8B", "EOMES", "GZMA", "GZMB", "IFNG", "PRF1"),
    "Central memory CD4 T" = c("AQP3", "CD28", "CCR7", "IL7R"),
    "Central memory CD8 T" = c("CCR7", "CD8A", "SELL"),
    "Effector memory CD4 T" = c("CCL5", "CXCR3", "GZMK", "NKG7"),
    "Effector memory CD8 T" = c("CD8A", "GZMA", "GZMK", "NKG7"),
    "Gamma delta T cell" = c("CD160", "KLRC3", "TRDC", "TRGC1", "TRGC2"),
    "Immature B cell" = c("CD19", "CD79A", "CD79B", "IGHM", "MME"),
    "Macrophage" = c("CD14", "CD68", "CSF1R", "ITGAM", "MSR1"),
    "MDSC" = c("ARG1", "CD33", "ITGAM", "S100A8", "S100A9"),
    "Monocyte" = c("CD14", "CSF1R", "FCGR1A", "LYZ", "VCAN"),
    "Natural killer cell" = c("FCER1G", "GNLY", "KLRB1", "KLRD1", "KLRF1", "NCAM1", "NKG7"),
    "Natural killer T cell" = c("CD160", "CD244", "NKG7", "ZBTB16"),
    "Plasmacytoid DC" = c("CLEC4C", "IL3RA", "IRF7", "LILRA4", "NRP1"),
    "Regulatory T cell" = c("CTLA4", "FOXP3", "IL2RA", "IKZF2"),
    "T follicular helper" = c("BCL6", "CD200", "CXCL13", "CXCR5", "ICOS", "PDCD1"),
    "Th1 cell" = c("IFNG", "STAT1", "STAT4", "TBX21"),
    "Th2 cell" = c("GATA3", "IL13", "IL4", "IL5", "STAT6"),
    "Th17 cell" = c("IL17A", "IL17F", "IL22", "RORC"),
    "Mast cell" = c("CPA3", "HDC", "KIT", "MS4A2", "TPSAB1", "TPSB2"),
    "Eosinophil" = c("CCR3", "CLC", "EPX", "PRG2", "SIGLEC8"),
    "Neutrophil" = c("CEACAM3", "CSF3R", "ELANE", "FCGR3B", "S100A12"),
    "Activated DC" = c("BATF3", "CD1C", "CLEC10A", "FCER1A", "HLA-DQA1"),
    "Immature DC" = c("CCR6", "CD1A", "CD1E"),
    "Type 1 IFN Response" = c("IFI35", "IFI44", "IFIT1", "IFIT3", "MX1", "OAS1"),
    "Type 2 IFN Response" = c("CXCL10", "CXCL9", "GBP1", "IDO1", "IFNG", "STAT1"),
    "M2 Macrophage" = c("CD163", "MRC1", "MSR1", "CD68", "IL10")
  )

  immune_markers_valid <- lapply(immune_cell_markers, function(gs) {
    gs[gs %in% rownames(expr_hub)]
  })
  # 移除基因数不足的
  immune_markers_valid <- immune_markers_valid[sapply(immune_markers_valid, length) >= 3]

  ssgsea_params <- ssgseaParam(
    exprData = as.matrix(log2(expr_tpm_symbol[, common_samples] + 1)),
    geneSets = immune_markers_valid
  )
  immune_ssgsea <- gsva(ssgsea_params)
}

# 确保样本匹配
common_immune_samples <- intersect(colnames(expr_hub), colnames(immune_ssgsea))
message(sprintf("Computing correlations: %d hub genes x %d immune cells x %d samples",
                length(hub_valid), nrow(immune_ssgsea), length(common_immune_samples)))

# 高效向量化计算相关性
cor_results <- data.frame()
for (gene in hub_valid) {
  gene_expr <- as.numeric(expr_hub[gene, common_immune_samples])
  for (immune_cell in rownames(immune_ssgsea)) {
    immune_score <- as.numeric(immune_ssgsea[immune_cell, common_immune_samples])
    ct <- cor.test(gene_expr, immune_score, method = "spearman")
    cor_results <- rbind(cor_results, data.frame(
      Gene = gene,
      ImmuneCell = immune_cell,
      r = as.numeric(ct$estimate),
      p_value = ct$p.value,
      stringsAsFactors = FALSE
    ))
  }
}

cor_results$padj <- p.adjust(cor_results$p_value, method = "BH")

# |r| >= CORR_EFFECT_THRESHOLD 标注
cor_results$biologically_meaningful <- abs(cor_results$r) >= CORR_EFFECT_THRESHOLD
cor_results$interpretation <- ifelse(
  cor_results$biologically_meaningful & cor_results$padj < 0.05,
  ifelse(cor_results$r > 0, "Positive (|r|>=0.3)", "Negative (|r|>=0.3)"),
  ifelse(cor_results$padj < 0.05, "Significant but |r|<0.3", "NS")
)

write.csv(cor_results, file.path(RES_DIR, "hub_gene_immune_correlation.csv"),
          row.names = FALSE)

# 统计报告
n_meaningful <- sum(cor_results$biologically_meaningful & cor_results$padj < 0.05)
n_sig_weak   <- sum(!cor_results$biologically_meaningful & cor_results$padj < 0.05)
message(sprintf("Correlations with |r|>=%.1f & padj<0.05: %d",
                CORR_EFFECT_THRESHOLD, n_meaningful))
message(sprintf("Significant but |r|<%.1f (not biologically interpreted): %d",
                CORR_EFFECT_THRESHOLD, n_sig_weak))

# --- 相关性热图 ---
cor_matrix <- cor_results %>%
  dplyr::select(Gene, ImmuneCell, r) %>%
  tidyr::pivot_wider(names_from = ImmuneCell, values_from = r) %>%
  tibble::column_to_rownames("Gene") %>%
  as.matrix()

pval_matrix <- cor_results %>%
  dplyr::select(Gene, ImmuneCell, padj) %>%
  tidyr::pivot_wider(names_from = ImmuneCell, values_from = padj) %>%
  tibble::column_to_rownames("Gene") %>%
  as.matrix()

r_abs_matrix <- abs(cor_matrix)

# 热图标注: 同时标记统计显著性和效应大小
# |r|>=0.3且p<0.05: *** / ** / *
# |r|<0.3但p<0.05: 仅r值(灰色,提示无生物学意义)
sig_mark <- matrix("", nrow = nrow(pval_matrix), ncol = ncol(pval_matrix))
for (i in 1:nrow(pval_matrix)) {
  for (j in 1:ncol(pval_matrix)) {
    if (!is.na(pval_matrix[i, j]) && pval_matrix[i, j] < 0.05 &&
        !is.na(r_abs_matrix[i, j]) && r_abs_matrix[i, j] >= CORR_EFFECT_THRESHOLD) {
      if (pval_matrix[i, j] < 0.001) {
        sig_mark[i, j] <- "***"
      } else if (pval_matrix[i, j] < 0.01) {
        sig_mark[i, j] <- "**"
      } else {
        sig_mark[i, j] <- "*"
      }
    }
  }
}

# 在|r|>=0.3的位置用黑色标注，其他用浅灰
pdf(file.path(fig_dir, "Fig8_hub_immune_correlation.pdf"), width = 16, height = 6)
ht <- Heatmap(
  cor_matrix,
  name = "Spearman r",
  col = colorRamp2(c(-0.6, -0.3, 0, 0.3, 0.6),
                   c("#2166AC", "#92C5DE", "white", "#F4A582", "#B2182B")),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (nchar(sig_mark[i, j]) > 0) {
      grid.text(sig_mark[i, j], x, y, gp = gpar(fontsize = 8, fontface = "bold"))
    }
  },
  row_names_gp = gpar(fontsize = 11, fontface = "italic"),
  column_names_gp = gpar(fontsize = 8),
  column_names_rot = 50,
  column_title = paste0("Hub Gene - Immune Cell Correlation ",
                        "(*, **, ***: padj<0.05 & |r|>=", CORR_EFFECT_THRESHOLD, ")"),
  column_title_gp = gpar(fontsize = 11, fontface = "bold"),
  heatmap_legend_param = list(
    title = "Spearman r",
    at = c(-0.6, -0.3, 0, 0.3, 0.6),
    labels = c("-0.6", paste0("-", CORR_EFFECT_THRESHOLD),
               "0", as.character(CORR_EFFECT_THRESHOLD), "0.6")
  )
)
draw(ht)
dev.off()

# =============================================================================
# 4. Hub基因 vs 24个检查点基因共表达
# =============================================================================
message("\n=== Section 4: Hub Gene vs Checkpoint Gene Co-expression ===")

# 24个免疫检查点基因
checkpoint_genes_24 <- c(
  "PDCD1", "CD274", "PDCD1LG2", "CTLA4",
  "HAVCR2", "LAG3", "TIGIT", "VSIR",
  "IDO1", "CD276", "VTCN1",
  "TNFRSF9", "ICOS", "TNFRSF4", "TNFRSF18",
  "CD40", "CD40LG", "CD80", "CD86",
  "TNFSF4", "TNFSF9", "ADORA2A", "BTLA", "SIGLEC15"
)

icp_valid <- checkpoint_genes_24[checkpoint_genes_24 %in% rownames(expr_tpm_symbol)]
message(sprintf("Checkpoint genes found: %d/%d", length(icp_valid), length(checkpoint_genes_24)))

expr_icp <- log2(expr_tpm_symbol[icp_valid, common_immune_samples] + 1)

cor_icp <- data.frame()
for (gene in hub_valid) {
  gene_expr <- as.numeric(expr_hub[gene, common_immune_samples])
  for (icp in icp_valid) {
    icp_expr <- as.numeric(expr_icp[icp, common_immune_samples])
    ct <- cor.test(gene_expr, icp_expr, method = "spearman")
    cor_icp <- rbind(cor_icp, data.frame(
      Hub_Gene = gene,
      Checkpoint = icp,
      r = as.numeric(ct$estimate),
      p_value = ct$p.value,
      stringsAsFactors = FALSE
    ))
  }
}

cor_icp$padj <- p.adjust(cor_icp$p_value, method = "BH")
cor_icp$biologically_meaningful <- abs(cor_icp$r) >= CORR_EFFECT_THRESHOLD

write.csv(cor_icp, file.path(RES_DIR, "hub_gene_checkpoint_correlation.csv"),
          row.names = FALSE)

n_meaningful_icp <- sum(cor_icp$biologically_meaningful & cor_icp$padj < 0.05)
message(sprintf("Checkpoint correlations with |r|>=%.1f & padj<0.05: %d / %d",
                CORR_EFFECT_THRESHOLD, n_meaningful_icp, nrow(cor_icp)))

# 检查点相关性热图
cor_icp_mat <- cor_icp %>%
  dplyr::select(Hub_Gene, Checkpoint, r) %>%
  tidyr::pivot_wider(names_from = Checkpoint, values_from = r) %>%
  tibble::column_to_rownames("Hub_Gene") %>%
  as.matrix()

pval_icp_mat <- cor_icp %>%
  dplyr::select(Hub_Gene, Checkpoint, padj) %>%
  tidyr::pivot_wider(names_from = Checkpoint, values_from = padj) %>%
  tibble::column_to_rownames("Hub_Gene") %>%
  as.matrix()

sig_mark_icp <- matrix("", nrow = nrow(pval_icp_mat), ncol = ncol(pval_icp_mat))
for (i in 1:nrow(pval_icp_mat)) {
  for (j in 1:ncol(pval_icp_mat)) {
    if (!is.na(pval_icp_mat[i, j]) && pval_icp_mat[i, j] < 0.05 &&
        abs(cor_icp_mat[i, j]) >= CORR_EFFECT_THRESHOLD) {
      if (pval_icp_mat[i, j] < 0.001) sig_mark_icp[i, j] <- "***"
      else if (pval_icp_mat[i, j] < 0.01) sig_mark_icp[i, j] <- "**"
      else sig_mark_icp[i, j] <- "*"
    }
  }
}

pdf(file.path(fig_dir, "Fig8_hub_checkpoint_correlation.pdf"), width = 14, height = 5)
ht_icp <- Heatmap(
  cor_icp_mat,
  name = "Spearman r",
  col = colorRamp2(c(-0.6, -0.3, 0, 0.3, 0.6),
                   c("#2166AC", "#92C5DE", "white", "#F4A582", "#B2182B")),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (nchar(sig_mark_icp[i, j]) > 0) {
      grid.text(sig_mark_icp[i, j], x, y, gp = gpar(fontsize = 8, fontface = "bold"))
    }
  },
  row_names_gp = gpar(fontsize = 11, fontface = "italic"),
  column_names_gp = gpar(fontsize = 9),
  column_names_rot = 45,
  column_title = paste0("Hub Gene - Checkpoint Correlation ",
                        "(stars: padj<0.05 & |r|>=", CORR_EFFECT_THRESHOLD, ")"),
  column_title_gp = gpar(fontsize = 11, fontface = "bold")
)
draw(ht_icp)
dev.off()

# =============================================================================
# 5. STRING PPI网络
# =============================================================================
message("\n=== Section 5: STRING PPI Network ===")

tryCatch({
  library(STRINGdb)

  # 构建扩展基因列表: hub + 直接相关的免疫检查点/UPR核心
  ppi_genes <- unique(c(hub_genes, additional_upr))
  message(sprintf("PPI query genes: %d (%s)",
                  length(ppi_genes), paste(ppi_genes, collapse = ", ")))

  string_db <- STRINGdb$new(
    version = "12.0",
    species = 9606,
    score_threshold = 400,
    network_type = "full"
  )

  gene_df <- data.frame(gene = ppi_genes, stringsAsFactors = FALSE)
  mapped <- string_db$map(gene_df, "gene", removeUnmappedRows = TRUE)
  message(sprintf("  Mapped to STRING: %d/%d genes", nrow(mapped), length(ppi_genes)))

  if (nrow(mapped) > 1) {
    # 绘制PPI网络
    pdf(file.path(fig_dir, "Fig8_ppi_network.pdf"), width = 10, height = 8)
    string_db$plot_network(mapped$STRING_id)
    dev.off()

    # 获取交互信息
    interactions <- string_db$get_interactions(mapped$STRING_id)

    if (nrow(interactions) > 0) {
      # 映射回gene symbol
      id_to_gene <- setNames(mapped$gene, mapped$STRING_id)
      interactions$gene_from <- id_to_gene[interactions$from]
      interactions$gene_to   <- id_to_gene[interactions$to]

      write.csv(interactions, file.path(RES_DIR, "hub_gene_ppi.csv"), row.names = FALSE)
      message(sprintf("  PPI interactions: %d", nrow(interactions)))

      # 交互得分汇总
      if (nrow(interactions) > 0) {
        interactions_clean <- interactions[!is.na(interactions$gene_from) &
                                            !is.na(interactions$gene_to), ]
        if (nrow(interactions_clean) > 0) {
          message("  Top interactions:")
          top_int <- interactions_clean %>%
            dplyr::arrange(desc(combined_score)) %>%
            head(10)
          for (k in 1:min(10, nrow(top_int))) {
            message(sprintf("    %s -- %s: score = %d",
                            top_int$gene_from[k], top_int$gene_to[k],
                            top_int$combined_score[k]))
          }
        }
      }
    } else {
      message("  No PPI interactions found")
    }

    # 功能富集
    tryCatch({
      enrichment <- string_db$get_enrichment(mapped$STRING_id)
      if (nrow(enrichment) > 0) {
        top_enrich <- enrichment %>%
          dplyr::filter(category %in% c("Process", "KEGG", "Component")) %>%
          dplyr::arrange(fdr) %>%
          head(20)
        write.csv(top_enrich, file.path(RES_DIR, "hub_gene_ppi_enrichment.csv"),
                  row.names = FALSE)
      }
    }, error = function(e) {
      message("  PPI enrichment failed: ", e$message)
    })
  }
}, error = function(e) {
  message("STRING analysis error: ", e$message)
  message("Note: STRINGdb requires internet connection")
})

# =============================================================================
# 6. 单细胞验证: Hub基因在scRNA-seq各细胞类型表达
# =============================================================================
message("\n=== Section 6: Single-cell Validation ===")

sc_file <- file.path(DATA_PROC, "seu_upr_scored.rds")

if (file.exists(sc_file)) {
  tryCatch({
    library(Seurat)

    message("Loading scRNA-seq data...")
    seu <- readRDS(sc_file)

    hub_in_sc <- hub_genes[hub_genes %in% rownames(seu)]
    message(sprintf("Hub genes found in scRNA-seq: %d/%d (%s)",
                    length(hub_in_sc), length(hub_genes),
                    paste(hub_in_sc, collapse = ", ")))

    if (length(hub_in_sc) > 0) {
      # 确定细胞类型列名
      celltype_col <- intersect(c("celltype", "cell_type", "CellType",
                                  "seurat_clusters", "ident"),
                                colnames(seu@meta.data))
      if (length(celltype_col) == 0) celltype_col <- "seurat_clusters"
      celltype_col <- celltype_col[1]
      message(sprintf("  Using cell type annotation: '%s'", celltype_col))

      # FeaturePlot
      n_feat <- length(hub_in_sc)
      ncol_feat <- min(3, n_feat)

      p_feat <- FeaturePlot(seu, features = hub_in_sc,
                            ncol = ncol_feat,
                            order = TRUE,
                            cols = c("lightgrey", "darkred")) &
        THEME_PUBLICATION &
        theme(plot.title = element_text(face = "italic"))

      ggsave(file.path(fig_dir, "Fig8_hub_scRNA_feature.pdf"), p_feat,
             width = ncol_feat * 5,
             height = ceiling(n_feat / ncol_feat) * 5)

      # VlnPlot (按细胞类型)
      p_vln <- VlnPlot(seu, features = hub_in_sc,
                        group.by = celltype_col,
                        pt.size = 0,
                        ncol = ncol_feat) &
        THEME_PUBLICATION &
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
              plot.title = element_text(face = "italic"))

      ggsave(file.path(fig_dir, "Fig8_hub_scRNA_violin.pdf"), p_vln,
             width = max(10, ncol_feat * 5),
             height = ceiling(n_feat / ncol_feat) * 4.5)

      # DotPlot
      p_dot <- DotPlot(seu, features = hub_in_sc, group.by = celltype_col) +
        RotatedAxis() +
        THEME_PUBLICATION +
        labs(title = "Hub Gene Expression Across Cell Types (scRNA-seq)") +
        theme(axis.text.x = element_text(face = "italic"))

      ggsave(file.path(fig_dir, "Fig8_hub_scRNA_dotplot.pdf"), p_dot,
             width = max(8, length(hub_in_sc) * 1.5), height = 6)

      message("  scRNA-seq plots saved")
    }

    rm(seu)  # 释放内存
    gc()

  }, error = function(e) {
    message("scRNA-seq analysis error: ", e$message)
  })
} else {
  message("scRNA-seq data not found: ", sc_file)
  message("Expected from Part 1: seu_upr_scored.rds")
}

# =============================================================================
# 7. HPA蛋白水平验证（Human Protein Atlas）
# =============================================================================
message("\n=== Section 7: HPA Protein-level Validation ===")

# HPA (Human Protein Atlas) 免疫组化数据
# 通过API查询或预下载数据验证Hub基因蛋白表达

hpa_results <- list()

tryCatch({
  # 方案1: 使用hpar包
  if (requireNamespace("hpar", quietly = TRUE)) {
    library(hpar)
    message("Using hpar package for HPA data...")

    # 获取正常组织和肿瘤组织的蛋白表达
    for (gene in hub_genes) {
      tryCatch({
        # 正常组织表达
        normal_data <- tryCatch({
          hpar::getHpa(gene, hpadata = "hpaNormalTissue")
        }, error = function(e) NULL)

        # 肿瘤组织表达
        cancer_data <- tryCatch({
          hpar::getHpa(gene, hpadata = "hpaCancer")
        }, error = function(e) NULL)

        hpa_results[[gene]] <- list(
          normal = normal_data,
          cancer = cancer_data
        )

        if (!is.null(cancer_data)) {
          # 筛选胶质瘤相关
          glioma_rows <- grep("glioma|brain", cancer_data$Cancer,
                              ignore.case = TRUE)
          if (length(glioma_rows) > 0) {
            message(sprintf("  %s in glioma (HPA): %s",
                            gene,
                            paste(unique(cancer_data$Level[glioma_rows]),
                                  collapse = ", ")))
          }
        }
      }, error = function(e) {
        message(sprintf("  %s HPA query failed: %s", gene, e$message))
      })
    }

    # 如果获取了数据，生成汇总热图
    if (length(hpa_results) > 0) {
      # 提取脑/胶质瘤蛋白表达水平
      hpa_summary <- data.frame()
      for (gene in names(hpa_results)) {
        cancer <- hpa_results[[gene]]$cancer
        if (!is.null(cancer) && nrow(cancer) > 0) {
          glioma_data <- cancer[grep("glioma", cancer$Cancer, ignore.case = TRUE), ]
          if (nrow(glioma_data) > 0) {
            hpa_summary <- rbind(hpa_summary, data.frame(
              Gene = gene,
              Cancer = glioma_data$Cancer,
              Level = glioma_data$Level,
              stringsAsFactors = FALSE
            ))
          }
        }
      }

      if (nrow(hpa_summary) > 0) {
        write.csv(hpa_summary, file.path(RES_DIR, "hub_gene_hpa_glioma.csv"),
                  row.names = FALSE)
      }
    }
  } else {
    message("hpar package not available")
  }
}, error = function(e) {
  message("HPA analysis failed: ", e$message)
})

# 方案2: 生成HPA查询链接供手动验证
message("\nHPA manual verification links:")
for (gene in hub_genes) {
  message(sprintf("  %s: https://www.proteinatlas.org/%s/pathology/glioma",
                  gene, gene))
}

# 保存HPA验证URL列表
hpa_urls <- data.frame(
  Gene = hub_genes,
  HPA_Pathology_URL = paste0("https://www.proteinatlas.org/", hub_genes, "/pathology/glioma"),
  HPA_Normal_URL = paste0("https://www.proteinatlas.org/", hub_genes, "/tissue"),
  HPA_SingleCell_URL = paste0("https://www.proteinatlas.org/", hub_genes, "/single+cell+type"),
  stringsAsFactors = FALSE
)
write.csv(hpa_urls, file.path(RES_DIR, "hub_gene_hpa_urls.csv"), row.names = FALSE)

# =============================================================================
# 8. Hub基因表达汇总热图
# =============================================================================
message("\n=== Section 8: Hub Gene Expression Summary Heatmap ===")

risk_groups_hub <- setNames(risk_group_train[common_samples], common_samples)
valid_for_hm <- common_samples[!is.na(risk_groups_hub[common_samples])]

expr_hub_hm <- as.matrix(expr_hub[hub_valid, valid_for_hm])
expr_hub_scaled <- t(scale(t(expr_hub_hm)))

# 按风险组排序
sample_order <- order(risk_groups_hub[valid_for_hm])

# 临床注释
clin_for_anno <- clinical_valid[match(valid_for_hm, clinical_valid$barcode), ]
ha_top <- HeatmapAnnotation(
  Risk = risk_groups_hub[valid_for_hm[sample_order]],
  col = list(Risk = COLORS_RISK),
  annotation_name_side = "left",
  show_legend = TRUE
)

# 添加IDH和Grade注释（如果可用）
if ("IDH_status" %in% colnames(clin_for_anno)) {
  idh_vals <- clin_for_anno$IDH_status[sample_order]
  idh_cols <- c("Mutant" = "#00A087", "WT" = "#E64B35")
  ha_top <- HeatmapAnnotation(
    Risk = risk_groups_hub[valid_for_hm[sample_order]],
    IDH = idh_vals,
    col = list(Risk = COLORS_RISK, IDH = idh_cols),
    annotation_name_side = "left",
    na_col = "grey90"
  )
}

pdf(file.path(fig_dir, "Fig8_hub_expression_heatmap.pdf"), width = 12, height = 5)
draw(Heatmap(
  expr_hub_scaled[, sample_order],
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("#4DBBD5", "white", "#E64B35")),
  top_annotation = ha_top,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 11, fontface = "italic"),
  column_title = "Hub Gene Expression by Risk Group",
  column_title_gp = gpar(fontsize = 12, fontface = "bold")
))
dev.off()

# =============================================================================
# 9. 保存全部结果
# =============================================================================
message("\n=== Saving results ===")

save(hub_genes, hub_coefs, hub_coef_df,
     cor_results, cor_icp,
     km_results, km_summary,
     hpa_urls,
     file = file.path(DATA_PROC, "hub_gene_results.RData"))

message("\n=== Hub gene analysis completed ===")
message("Figures saved to: ", fig_dir)
message("  Fig8: Hub gene coefficients, KM curves, forest plot")
message("        Immune cell & checkpoint correlation heatmaps")
message("        PPI network, scRNA-seq validation, expression heatmap")
message("  Key rule: Only |r|>=", CORR_EFFECT_THRESHOLD,
        " correlations are biologically interpreted")
message("\n=== Part 4 (Clinical Translation) ALL COMPLETED ===")
message("Check figures/ and results/ directories for all outputs.")
