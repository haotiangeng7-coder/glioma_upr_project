###############################################################################
# 06_pathway_analysis.R
# 各UPR亚型的通路活性分析（Figure 4E）
# - GSVA计算50条Hallmark通路活性
# - limma鉴定亚型间差异通路
# - GSEA: UPR-favorable vs UPR-high-risk GO:BP富集
###############################################################################

source("00_setup/config.R")
library(GSVA)
library(msigdbr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(limma)
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)
library(pheatmap)

set.seed(SEED)
load(file.path(DATA_PROC, "tcga_glioma_expression.RData"))
load(file.path(DATA_PROC, "consensus_clustering_results.RData"))

# =============================================================================
# 1. 准备数据
# =============================================================================
message("=== Preparing data for pathway analysis ===")

common_samples <- cluster_df$barcode[cluster_df$barcode %in% colnames(expr_tpm_symbol)]
expr_log <- log2(expr_tpm_symbol[, common_samples] + 1)
subtypes <- setNames(cluster_df$UPR_subtype, cluster_df$barcode)[common_samples]

subtype_levels <- sort(unique(as.character(subtypes)))
comparisons_list <- combn(subtype_levels, 2, simplify = FALSE)

message(sprintf("Samples: %d | Subtypes: %s",
                length(common_samples),
                paste(names(table(subtypes)), table(subtypes), sep = "=", collapse = ", ")))

# =============================================================================
# 2. GSVA — 50条Hallmark通路活性评分
# =============================================================================
message("\n=== GSVA Hallmark pathway analysis (50 pathways) ===")

# 获取Hallmark基因集
hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_list <- split(hallmark$gene_symbol, hallmark$gs_name)
message(sprintf("  Hallmark gene sets loaded: %d", length(hallmark_list)))

# GSVA计算
gsva_params <- gsvaParam(
  exprData = as.matrix(expr_log),
  geneSets = hallmark_list,
  minSize = 10,
  maxSize = 500
)
gsva_hallmark <- gsva(gsva_params)

message(sprintf("  GSVA completed: %d pathways x %d samples",
                nrow(gsva_hallmark), ncol(gsva_hallmark)))

# =============================================================================
# 3. limma鉴定亚型间差异通路
# =============================================================================
message("\n=== Differential pathway activity (limma) ===")

# 设计矩阵
design <- model.matrix(~ 0 + factor(subtypes))
colnames(design) <- gsub("factor\\(subtypes\\)", "", colnames(design))
colnames(design) <- gsub("-", "_", colnames(design))

fit <- lmFit(gsva_hallmark, design)

# 根据亚型数量构建对比矩阵
subtype_levels_safe <- gsub("-", "_", sort(unique(subtypes)))

# Build all pairwise contrasts from the ACTUAL subtype levels present (works for any K).
# Names derive from the real labels (e.g. UPR_high_risk - UPR_favorable), not a fixed
# activated/quiescent/intermediate scheme.
contrast_formulas <- character()
for (i in 1:(length(subtype_levels_safe) - 1)) {
  for (j in (i + 1):length(subtype_levels_safe)) {
    contrast_formulas <- c(contrast_formulas,
      paste0(subtype_levels_safe[i], " - ", subtype_levels_safe[j]))
  }
}
contrast_mat <- makeContrasts(contrasts = contrast_formulas, levels = design)

fit2 <- contrasts.fit(fit, contrast_mat)
fit2 <- eBayes(fit2)

# 提取各对比的差异通路结果
diff_pathway_results <- list()
for (coef_name in colnames(contrast_mat)) {
  res <- topTable(fit2, coef = coef_name, number = Inf, sort.by = "t")
  res$pathway <- rownames(res)

  # 简化通路名用于展示
  res$pathway_short <- gsub("HALLMARK_", "", res$pathway)
  res$pathway_short <- gsub("_", " ", res$pathway_short)

  diff_pathway_results[[coef_name]] <- res

  write.csv(res, file.path(RES_DIR, paste0("gsva_diff_", coef_name, ".csv")),
            row.names = FALSE)

  sig_count <- sum(res$adj.P.Val < FDR_CUTOFF)
  message(sprintf("  %s: %d significant pathways (adj.P.Val < %.2f)",
                  coef_name, sig_count, FDR_CUTOFF))
}

# =============================================================================
# 4. Figure 4E: Hallmark通路热图
# =============================================================================
message("\n=== Figure 4E: Hallmark pathway heatmap ===")

sample_order <- order(subtypes)

# 简化通路名称
rownames_simple <- gsub("HALLMARK_", "", rownames(gsva_hallmark))
rownames_simple <- gsub("_", " ", rownames_simple)

gsva_plot <- gsva_hallmark
rownames(gsva_plot) <- rownames_simple

# Annotate significant pathways using the primary (first) pairwise contrast,
# whatever the subtype labels are (K=2 -> the single favorable-vs-high-risk contrast).
if (length(diff_pathway_results) >= 1) {
  aq_res <- diff_pathway_results[[1]]
  sig_pathways <- aq_res$pathway[aq_res$adj.P.Val < FDR_CUTOFF]
  sig_pathways_short <- gsub("HALLMARK_", "", sig_pathways)
  sig_pathways_short <- gsub("_", " ", sig_pathways_short)

  row_sig <- ifelse(rownames(gsva_plot) %in% sig_pathways_short, "Significant", "NS")
} else {
  row_sig <- rep("NS", nrow(gsva_plot))
}

# 行侧注释
row_anno <- rowAnnotation(
  Diff = row_sig,
  col = list(Diff = c("Significant" = "#E64B35", "NS" = "grey85")),
  annotation_name_side = "top",
  simple_anno_size = unit(3, "mm")
)

# 列注释
col_anno <- HeatmapAnnotation(
  UPR_subtype = subtypes[common_samples[sample_order]],
  col = list(UPR_subtype = COLORS_SUBTYPE),
  annotation_name_side = "left"
)

pdf(file.path(FIG_DIR, "Fig4E_gsva_hallmark_heatmap.pdf"), width = 14, height = 14)
ht <- Heatmap(
  gsva_plot[, common_samples[sample_order]],
  name = "GSVA\nScore",
  col = colorRamp2(c(-0.5, 0, 0.5), c("#2166AC", "white", "#B2182B")),
  top_annotation = col_anno,
  right_annotation = row_anno,
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 7),
  cluster_columns = FALSE,
  column_split = subtypes[common_samples[sample_order]],
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_title = "Hallmark Pathways",
  column_title = "Hallmark Pathway Activity by UPR Subtype"
)
draw(ht)
dev.off()

# =============================================================================
# 5. 关键通路聚焦箱线图
# =============================================================================
message("\n=== Key pathway boxplots ===")

key_pathways <- c(
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_HYPOXIA",
  "HALLMARK_ANGIOGENESIS",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_OXIDATIVE_PHOSPHORYLATION",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_PI3K_AKT_MTOR_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_DNA_REPAIR"
)

key_pathways_valid <- key_pathways[key_pathways %in% rownames(gsva_hallmark)]
message(sprintf("  Key pathways found: %d/%d", length(key_pathways_valid), length(key_pathways)))

if (length(key_pathways_valid) > 0) {
  key_data <- as.data.frame(t(gsva_hallmark[key_pathways_valid, common_samples]))
  key_data$UPR_subtype <- subtypes
  key_data$barcode <- rownames(key_data)

  key_long <- key_data %>%
    pivot_longer(cols = -c(barcode, UPR_subtype),
                 names_to = "Pathway", values_to = "Score") %>%
    mutate(Pathway = gsub("HALLMARK_", "", Pathway),
           Pathway = gsub("_", " ", Pathway))

  # 统计检验
  key_kw <- key_long %>%
    group_by(Pathway) %>%
    summarise(
      kw_pvalue = kruskal.test(Score ~ UPR_subtype)$p.value,
      .groups = "drop"
    ) %>%
    mutate(kw_padj = p.adjust(kw_pvalue, method = "BH"))

  p_key <- ggplot(key_long, aes(x = UPR_subtype, y = Score, fill = UPR_subtype)) +
    geom_boxplot(outlier.size = 0.3, width = 0.6) +
    facet_wrap(~Pathway, scales = "free_y", ncol = 5) +
    stat_compare_means(method = "kruskal.test", size = 2.5, label = "p.format") +
    stat_compare_means(comparisons = comparisons_list, method = "wilcox.test",
                       p.adjust.method = "BH", label = "p.signif",
                       size = 2, hide.ns = TRUE, step.increase = 0.06) +
    scale_fill_manual(values = COLORS_SUBTYPE) +
    THEME_PUBLICATION +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          legend.position = "bottom",
          strip.text = element_text(size = 7)) +
    labs(x = "UPR Subtype", y = "GSVA Enrichment Score") +
    ggtitle("Key Pathway Activity by UPR Subtype")

  ggsave(file.path(FIG_DIR, "Fig4E_key_pathways_boxplots.pdf"), p_key,
         width = 18, height = max(6, ceiling(length(key_pathways_valid) / 5) * 3.5))

  write.csv(key_kw, file.path(RES_DIR, "key_pathways_kw_test.csv"), row.names = FALSE)
}

# =============================================================================
# 6. GSEA — UPR-favorable vs UPR-high-risk GO:BP富集
# =============================================================================
message("\n=== GSEA: UPR-favorable vs UPR-high-risk (GO:BP) ===")

activated_samples <- common_samples[subtypes == "UPR-favorable"]
quiescent_samples <- common_samples[subtypes == "UPR-high-risk"]

gsea_go <- NULL
de_results <- NULL

if (length(activated_samples) > 10 && length(quiescent_samples) > 10) {
  message(sprintf("  favorable: %d samples, high-risk: %d samples",
                  length(activated_samples), length(quiescent_samples)))

  # limma差异分析 -> 获取ranked gene list
  group <- factor(c(rep("favorable", length(activated_samples)),
                     rep("high_risk", length(quiescent_samples))),
                  levels = c("favorable", "high_risk"))
  design_de <- model.matrix(~ 0 + group)
  colnames(design_de) <- levels(group)

  expr_de <- expr_log[, c(activated_samples, quiescent_samples)]

  fit_de <- lmFit(expr_de, design_de)
  contrast_de <- makeContrasts(favorable - high_risk, levels = design_de)
  fit_de2 <- contrasts.fit(fit_de, contrast_de)
  fit_de2 <- eBayes(fit_de2)

  de_results <- topTable(fit_de2, number = Inf, sort.by = "t")
  de_results$gene <- rownames(de_results)

  n_up <- sum(de_results$adj.P.Val < FDR_CUTOFF & de_results$logFC > LOGFC_CUTOFF)
  n_down <- sum(de_results$adj.P.Val < FDR_CUTOFF & de_results$logFC < -LOGFC_CUTOFF)
  message(sprintf("  DE genes: %d up, %d down (FDR<%.2f, |logFC|>%.1f)",
                  n_up, n_down, FDR_CUTOFF, LOGFC_CUTOFF))

  # ranked gene list (按t统计量排序)
  gene_list <- setNames(de_results$t, de_results$gene)
  gene_list <- sort(gene_list, decreasing = TRUE)

  # 去除NA和重复
  gene_list <- gene_list[!is.na(gene_list)]
  gene_list <- gene_list[!duplicated(names(gene_list))]

  # GO:BP GSEA
  go_bp <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP") %>%
    dplyr::select(gs_name, gene_symbol)

  gsea_go <- GSEA(
    geneList = gene_list,
    TERM2GENE = go_bp,
    pvalueCutoff = FDR_CUTOFF,
    minGSSize = 10,
    maxGSSize = 500,
    seed = SEED,
    eps = 1e-30
  )

  n_enriched <- nrow(gsea_go@result)
  n_up_terms <- sum(gsea_go@result$NES > 0)
  n_down_terms <- sum(gsea_go@result$NES < 0)
  message(sprintf("  GSEA GO:BP enriched terms: %d (up: %d, down: %d)",
                  n_enriched, n_up_terms, n_down_terms))

  if (n_enriched > 0) {
    # Dotplot: 分别展示上调和下调的top GO terms
    p_gsea_dot <- dotplot(gsea_go, showCategory = 20, orderBy = "NES",
                          split = ".sign") +
      facet_grid(~.sign) +
      THEME_PUBLICATION +
      theme(axis.text.y = element_text(size = 7)) +
      ggtitle("GO:BP GSEA (UPR-favorable vs UPR-high-risk)")

    ggsave(file.path(FIG_DIR, "Fig4E_gsea_gobp_dotplot.pdf"), p_gsea_dot,
           width = 14, height = 10)

    # 选取最显著的上调和下调各5条做enrichment plot
    top_up <- gsea_go@result %>%
      filter(NES > 0) %>%
      arrange(p.adjust) %>%
      head(5) %>%
      pull(ID)

    top_down <- gsea_go@result %>%
      filter(NES < 0) %>%
      arrange(p.adjust) %>%
      head(5) %>%
      pull(ID)

    top_terms <- c(top_up, top_down)

    if (length(top_terms) > 0) {
      # 使用enrichplot的ridgeplot展示GSEA结果
      tryCatch({
        library(enrichplot)
        p_ridge <- ridgeplot(gsea_go, showCategory = 20) +
          THEME_PUBLICATION +
          theme(axis.text.y = element_text(size = 6)) +
          ggtitle("GSEA GO:BP Ridge Plot")

        ggsave(file.path(FIG_DIR, "FigS_gsea_gobp_ridgeplot.pdf"), p_ridge,
               width = 10, height = 10)
      }, error = function(e) {
        message("  Ridge plot error: ", e$message)
      })

      # 单独GSEA running score plot（最显著的通路）
      if (length(top_up) > 0) {
        pdf(file.path(FIG_DIR, "FigS_gsea_running_score.pdf"), width = 8, height = 6)
        for (term_id in top_up[1:min(3, length(top_up))]) {
          tryCatch({
            p_running <- gseaplot2(gsea_go, geneSetID = term_id,
                                   title = gsub("GOBP_", "", term_id),
                                   pvalue_table = TRUE)
            print(p_running)
          }, error = function(e) {
            message(sprintf("    Running score plot error for %s: %s", term_id, e$message))
          })
        }
        dev.off()
      }
    }

    # GSEA NES条形图（top 30 terms）
    gsea_top <- gsea_go@result %>%
      arrange(p.adjust) %>%
      head(30) %>%
      mutate(
        pathway_short = gsub("GOBP_", "", ID),
        pathway_short = gsub("_", " ", pathway_short),
        pathway_short = ifelse(nchar(pathway_short) > 50,
                               paste0(substr(pathway_short, 1, 47), "..."),
                               pathway_short),
        Direction = ifelse(NES > 0, "Up in UPR-favorable", "Up in UPR-high-risk")
      )

    p_nes_bar <- ggplot(gsea_top, aes(x = reorder(pathway_short, NES), y = NES,
                                       fill = Direction)) +
      geom_bar(stat = "identity", width = 0.7) +
      coord_flip() +
      scale_fill_manual(values = c("Up in UPR-favorable" = "#4DBBD5",
                                    "Up in UPR-high-risk" = "#E64B35")) +
      geom_hline(yintercept = 0, linetype = "solid", color = "grey30") +
      labs(x = "", y = "Normalized Enrichment Score (NES)",
           fill = "Direction",
           title = "Top GO:BP Terms (UPR-favorable vs UPR-high-risk)") +
      THEME_PUBLICATION +
      theme(axis.text.y = element_text(size = 7),
            legend.position = "bottom")

    ggsave(file.path(FIG_DIR, "Fig4E_gsea_nes_barplot.pdf"), p_nes_bar,
           width = 12, height = max(8, nrow(gsea_top) * 0.3))
  }

  # 保存DE和GSEA结果
  write.csv(de_results, file.path(RES_DIR, "de_favorable_vs_high_risk.csv"),
            row.names = FALSE)
  write.csv(gsea_go@result,
            file.path(RES_DIR, "gsea_gobp_activated_vs_quiescent.csv"),
            row.names = FALSE)

} else {
  message("  Insufficient samples for GSEA analysis.")
  message(sprintf("    activated: %d, quiescent: %d (need >10 each)",
                  length(activated_samples), length(quiescent_samples)))
}

# =============================================================================
# 7. 亚型特异性通路活性 — 每个亚型vs其余的差异通路
# =============================================================================
message("\n=== Subtype-specific pathway activity (each vs rest) ===")

subtype_specific_pathways <- list()

for (st in subtype_levels) {
  st_safe <- gsub("-", "_", st)

  # one vs rest设计
  group_ovr <- factor(ifelse(subtypes == st, st_safe, "rest"),
                      levels = c(st_safe, "rest"))
  design_ovr <- model.matrix(~ 0 + group_ovr)
  colnames(design_ovr) <- levels(group_ovr)

  fit_ovr <- lmFit(gsva_hallmark, design_ovr)
  contrast_ovr <- makeContrasts(contrasts = paste0(st_safe, " - rest"), levels = design_ovr)
  fit_ovr2 <- contrasts.fit(fit_ovr, contrast_ovr)
  fit_ovr2 <- eBayes(fit_ovr2)

  res_ovr <- topTable(fit_ovr2, number = Inf, sort.by = "t")
  res_ovr$pathway <- rownames(res_ovr)
  res_ovr$subtype <- st

  subtype_specific_pathways[[st]] <- res_ovr

  sig_up <- sum(res_ovr$adj.P.Val < FDR_CUTOFF & res_ovr$logFC > 0)
  sig_down <- sum(res_ovr$adj.P.Val < FDR_CUTOFF & res_ovr$logFC < 0)
  message(sprintf("  %s vs rest: %d up, %d down (adj.P.Val < %.2f)",
                  st, sig_up, sig_down, FDR_CUTOFF))
}

# 合并并保存
all_specific <- do.call(rbind, subtype_specific_pathways)
write.csv(all_specific, file.path(RES_DIR, "gsva_subtype_specific_pathways.csv"),
          row.names = FALSE)

# =============================================================================
# 8. 保存通路分析结果
# =============================================================================
message("\n=== Saving pathway analysis results ===")

save(gsva_hallmark, diff_pathway_results, subtype_specific_pathways,
     gsea_go, de_results,
     file = file.path(DATA_PROC, "pathway_analysis_results.RData"))

message("\n=== Pathway analysis completed (Figure 4E) ===")
message("Part 2 (Bulk subtyping) completed.")
message("Next: Run 05_ml_model/01_feature_selection.R for Part 3.")
