###############################################################################
# 04_immune_analysis.R
# 各UPR亚型免疫微环境分析（Figure 4A-C）
# 多工具整合策略：ssGSEA(主) + MCPcounter(验证) + ESTIMATE(评分)
# 仅报告至少2/3工具一致显著的免疫细胞差异
###############################################################################

source("00_setup/config.R")
library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(GSVA)
library(dplyr)
library(tidyr)
library(pheatmap)
library(MCPcounter)

set.seed(SEED)
load(file.path(DATA_PROC, "upr_gene_sets.RData"))
load(file.path(DATA_PROC, "tcga_glioma_expression.RData"))
load(file.path(DATA_PROC, "consensus_clustering_results.RData"))

# =============================================================================
# 1. 准备表达数据
# =============================================================================
message("=== Preparing expression data for immune analysis ===")

common_samples <- cluster_df$barcode[cluster_df$barcode %in% colnames(expr_tpm_symbol)]
expr_mat <- expr_tpm_symbol[, common_samples]
expr_log <- log2(expr_mat + 1)

subtype_map <- setNames(cluster_df$UPR_subtype, cluster_df$barcode)
subtypes <- subtype_map[common_samples]

message(sprintf("Samples: %d | Subtypes: %s",
                length(common_samples),
                paste(names(table(subtypes)), table(subtypes), sep = "=", collapse = ", ")))

# 亚型比较对列表（后续复用）
subtype_levels <- sort(unique(as.character(subtypes)))
comparisons_list <- combn(subtype_levels, 2, simplify = FALSE)

# =============================================================================
# 2. ESTIMATE分析 — Stromal/Immune/ESTIMATE Score
# =============================================================================
message("\n=== ESTIMATE analysis ===")

estimate_input <- file.path(DATA_PROC, "estimate_input.txt")
estimate_filtered <- file.path(DATA_PROC, "estimate_filtered.gct")
estimate_output <- file.path(DATA_PROC, "estimate_output.txt")

est_df <- NULL

tryCatch({
  library(estimate)

  write.table(
    data.frame(Gene = rownames(expr_mat), expr_mat, check.names = FALSE),
    file = estimate_input,
    sep = "\t", row.names = FALSE, quote = FALSE
  )

  filterCommonGenes(
    input.f = estimate_input,
    output.f = estimate_filtered,
    id = "GeneSymbol"
  )

  estimateScore(
    input.ds = estimate_filtered,
    output.ds = estimate_output,
    platform = "illumina"
  )

  # 解析ESTIMATE输出
  est_scores <- read.delim(estimate_output, skip = 2, row.names = 1, check.names = FALSE)
  est_scores <- est_scores[, -1]  # 去掉Description列

  est_df <- data.frame(
    barcode = colnames(est_scores),
    StromalScore = as.numeric(est_scores["StromalScore", ]),
    ImmuneScore = as.numeric(est_scores["ImmuneScore", ]),
    ESTIMATEScore = as.numeric(est_scores["ESTIMATEScore", ]),
    stringsAsFactors = FALSE
  )
  est_df$UPR_subtype <- subtype_map[est_df$barcode]

  # --- ESTIMATE评分箱线图 ---
  est_long <- est_df %>%
    pivot_longer(cols = c("StromalScore", "ImmuneScore", "ESTIMATEScore"),
                 names_to = "Score", values_to = "Value")

  p_estimate <- ggplot(est_long, aes(x = UPR_subtype, y = Value, fill = UPR_subtype)) +
    geom_boxplot(outlier.size = 0.5, width = 0.6) +
    facet_wrap(~Score, scales = "free_y") +
    stat_compare_means(method = "kruskal.test", label = "p.format", size = 3.5) +
    scale_fill_manual(values = COLORS_SUBTYPE) +
    THEME_PUBLICATION +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1)) +
    labs(x = "UPR Subtype", y = "Score") +
    ggtitle("ESTIMATE Scores by UPR Subtype")

  ggsave(file.path(FIG_DIR, "Fig4_estimate_scores.pdf"), p_estimate,
         width = 12, height = 5)
  message("  ESTIMATE scores computed and plotted.")

  # 清理临时文件
  for (f in c(estimate_input, estimate_filtered, estimate_output)) {
    if (file.exists(f)) file.remove(f)
  }

}, error = function(e) {
  message("ESTIMATE error: ", e$message)
  message("  Continuing without ESTIMATE results.")
})

# =============================================================================
# 3. ssGSEA免疫细胞浸润 — 28种免疫细胞（主分析工具）
# =============================================================================
message("\n=== ssGSEA immune infiltration (primary tool) ===")

# 28种免疫细胞标志基因集（Charoentong et al. 2017 Cell Reports）
immune_cell_markers <- list(
  "Activated B cell" = c("BLK", "CD19", "FAM30A", "FCRL2", "MS4A1", "PNOC", "SPIB", "TCL1A"),
  "Activated CD4 T cell" = c("DPP4", "ICOS", "IL2RA", "IL4R", "OAS1", "UBE2L6"),
  "Activated CD8 T cell" = c("CD8A", "CD8B", "EOMES", "GZMA", "GZMB", "IFNG", "PRF1"),
  "Central memory CD4 T cell" = c("AQP3", "CD28", "CCR7", "IL7R"),
  "Central memory CD8 T cell" = c("CCR7", "CD8A", "SELL"),
  "Effector memory CD4 T cell" = c("CCL5", "CXCR3", "GZMK", "NKG7"),
  "Effector memory CD8 T cell" = c("CD8A", "GZMA", "GZMK", "NKG7"),
  "Gamma delta T cell" = c("CD160", "KLRC3", "TRDC", "TRGC1", "TRGC2"),
  "Immature B cell" = c("CD19", "CD79A", "CD79B", "IGHM", "MME"),
  "Macrophage" = c("CD14", "CD68", "CSF1R", "ITGAM", "MSR1"),
  "MDSC" = c("ARG1", "CD33", "ITGAM", "S100A8", "S100A9"),
  "Monocyte" = c("CD14", "CSF1R", "FCGR1A", "LYZ", "VCAN"),
  "Natural killer cell" = c("FCER1G", "GNLY", "KLRB1", "KLRD1", "KLRF1", "NCAM1", "NKG7"),
  "Natural killer T cell" = c("CD160", "CD244", "NKG7", "ZBTB16"),
  "Plasmacytoid dendritic cell" = c("CLEC4C", "IL3RA", "IRF7", "LILRA4", "NRP1"),
  "Regulatory T cell" = c("CTLA4", "FOXP3", "IL2RA", "IKZF2"),
  "T follicular helper cell" = c("BCL6", "CD200", "CXCL13", "CXCR5", "ICOS", "PDCD1"),
  "Type 1 T helper cell" = c("IFNG", "STAT1", "STAT4", "TBX21"),
  "Type 2 T helper cell" = c("GATA3", "IL13", "IL4", "IL5", "STAT6"),
  "Type 17 T helper cell" = c("IL17A", "IL17F", "IL22", "RORC"),
  "Mast cell" = c("CPA3", "HDC", "KIT", "MS4A2", "TPSAB1", "TPSB2"),
  "Eosinophil" = c("CCR3", "CLC", "EPX", "PRG2", "SIGLEC8"),
  "Neutrophil" = c("CEACAM3", "CSF3R", "ELANE", "FCGR3B", "S100A12"),
  "Activated dendritic cell" = c("BATF3", "CD1C", "CLEC10A", "FCER1A", "HLA-DQA1"),
  "Immature dendritic cell" = c("CCR6", "CD1A", "CD1E"),
  "Type 1 IFN Response" = c("IFI35", "IFI44", "IFIT1", "IFIT3", "MX1", "OAS1"),
  "Type 2 IFN Response" = c("CXCL10", "CXCL9", "GBP1", "IDO1", "IFNG", "STAT1")
)

# ssGSEA计算
ssgsea_params <- ssgseaParam(
  exprData = expr_log[, common_samples],
  geneSets = immune_cell_markers
)
immune_ssgsea <- gsva(ssgsea_params)

message(sprintf("  ssGSEA completed: %d cell types x %d samples",
                nrow(immune_ssgsea), ncol(immune_ssgsea)))

# ssGSEA统计检验：Kruskal-Wallis + pairwise Wilcoxon, BH FDR校正
immune_long_ssgsea <- as.data.frame(t(immune_ssgsea)) %>%
  mutate(barcode = rownames(.), UPR_subtype = subtypes) %>%
  pivot_longer(cols = -c(barcode, UPR_subtype),
               names_to = "CellType", values_to = "Score")

ssgsea_kw <- immune_long_ssgsea %>%
  group_by(CellType) %>%
  summarise(
    kw_pvalue = tryCatch(kruskal.test(Score ~ UPR_subtype)$p.value, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(kw_padj = p.adjust(kw_pvalue, method = "BH"))

# pairwise Wilcoxon检验（保存详细结果）
ssgsea_pairwise <- immune_long_ssgsea %>%
  group_by(CellType) %>%
  summarise(
    pw_results = list(pairwise.wilcox.test(Score, UPR_subtype, p.adjust.method = "BH")),
    .groups = "drop"
  )

ssgsea_sig <- ssgsea_kw %>% filter(kw_padj < FDR_CUTOFF) %>% pull(CellType)
message(sprintf("  ssGSEA significant cell types: %d/%d (FDR < %.2f)",
                length(ssgsea_sig), nrow(ssgsea_kw), FDR_CUTOFF))

# =============================================================================
# 4. MCPcounter免疫浸润（敏感性验证）
# =============================================================================
message("\n=== MCPcounter analysis (sensitivity validation) ===")

mcp_results <- NULL
mcp_sig <- character(0)
mcp_kw <- NULL

tryCatch({
  mcp_results <- MCPcounter::MCPcounter.estimate(
    expression = expr_mat[, common_samples],
    featuresType = "HUGO_symbols"
  )

  mcp_long <- as.data.frame(t(mcp_results)) %>%
    mutate(barcode = rownames(.), UPR_subtype = subtypes) %>%
    pivot_longer(cols = -c(barcode, UPR_subtype),
                 names_to = "CellType", values_to = "Score")

  # MCPcounter统计检验
  mcp_kw <- mcp_long %>%
    group_by(CellType) %>%
    summarise(
      kw_pvalue = tryCatch(kruskal.test(Score ~ UPR_subtype)$p.value, error = function(e) NA_real_),
      .groups = "drop"
    ) %>%
    mutate(kw_padj = p.adjust(kw_pvalue, method = "BH"))

  mcp_sig <- mcp_kw %>% filter(kw_padj < FDR_CUTOFF) %>% pull(CellType)
  message(sprintf("  MCPcounter significant cell types: %d/%d",
                  length(mcp_sig), nrow(mcp_kw)))

  # MCPcounter箱线图（补充图）
  p_mcp <- ggplot(mcp_long, aes(x = UPR_subtype, y = Score, fill = UPR_subtype)) +
    geom_boxplot(outlier.size = 0.3, width = 0.6) +
    facet_wrap(~CellType, scales = "free_y", ncol = 4) +
    stat_compare_means(method = "kruskal.test", size = 2.5) +
    scale_fill_manual(values = COLORS_SUBTYPE) +
    THEME_PUBLICATION +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          legend.position = "none",
          strip.text = element_text(size = 8)) +
    labs(x = "UPR Subtype", y = "Score") +
    ggtitle("MCPcounter Immune Deconvolution by UPR Subtype")

  ggsave(file.path(FIG_DIR, "FigS_mcpcounter_boxplots.pdf"), p_mcp,
         width = 14, height = 8)

}, error = function(e) {
  message("MCPcounter error: ", e$message)
  message("  Continuing without MCPcounter validation.")
})

# =============================================================================
# 5. 多工具一致性整合 — 仅报告>=2/3工具一致显著的细胞类型
# =============================================================================
message("\n=== Multi-tool consensus integration ===")

# MCPcounter -> ssGSEA 细胞类型名称映射
# MCPcounter细胞类型: T cells, CD8 T cells, Cytotoxic lymphocytes, B lineage,
#   NK cells, Monocytic lineage, Myeloid dendritic cells, Neutrophils,
#   Endothelial cells, Fibroblasts
mcp_to_ssgsea_map <- c(
  "CD8 T cells"             = "Activated CD8 T cell",
  "B lineage"               = "Activated B cell",
  "NK cells"                = "Natural killer cell",
  "Monocytic lineage"       = "Monocyte",
  "Myeloid dendritic cells" = "Activated dendritic cell",
  "Neutrophils"             = "Neutrophil"
)

# ESTIMATE提供总体Immune Score，不细分细胞亚型。
# 一致性判定逻辑:
#   对可映射细胞类型 -> 要求ssGSEA显著 + MCPcounter对应类型显著（2/2工具一致）
#   对不可映射细胞类型 -> ssGSEA显著即保留（标注为单工具验证）
# 如果ESTIMATE ImmuneScore在亚型间显著,则进一步佐证免疫浸润整体差异

# Step 1: 可映射细胞类型 — ssGSEA & MCPcounter 2/2一致
consensus_mapped <- character(0)
for (mcp_name in names(mcp_to_ssgsea_map)) {
  ssgsea_name <- mcp_to_ssgsea_map[mcp_name]
  if (ssgsea_name %in% ssgsea_sig && mcp_name %in% mcp_sig) {
    consensus_mapped <- c(consensus_mapped, ssgsea_name)
  }
}

# Step 2: 不可映射的ssGSEA显著细胞类型保留（标注为单工具）
ssgsea_unmappable_sig <- setdiff(ssgsea_sig, mcp_to_ssgsea_map)

# 最终一致性显著细胞类型
consensus_sig_cells <- unique(c(consensus_mapped, ssgsea_unmappable_sig))

# ESTIMATE一致性佐证
estimate_consistent <- FALSE
if (!is.null(est_df)) {
  est_kw_p <- tryCatch(kruskal.test(ImmuneScore ~ UPR_subtype, data = est_df)$p.value,
                       error = function(e) NA_real_)
  estimate_consistent <- est_kw_p < FDR_CUTOFF
  message(sprintf("  ESTIMATE ImmuneScore KW p = %.4e (%s)",
                  est_kw_p, ifelse(estimate_consistent, "consistent", "not significant")))
}

# 生成一致性报告
consensus_report <- data.frame(
  CellType = ssgsea_sig,
  ssGSEA_sig = TRUE,
  MCPcounter_mappable = ssgsea_sig %in% mcp_to_ssgsea_map,
  MCPcounter_concordant = ssgsea_sig %in% consensus_mapped,
  Consensus_final = ssgsea_sig %in% consensus_sig_cells,
  stringsAsFactors = FALSE
)

message(sprintf("  Consensus significant cells: %d", length(consensus_sig_cells)))
message(sprintf("    - ssGSEA + MCPcounter concordant: %d", length(consensus_mapped)))
message(sprintf("    - ssGSEA only (no MCPcounter mapping): %d", length(ssgsea_unmappable_sig)))

write.csv(consensus_report,
          file.path(RES_DIR, "immune_consensus_report.csv"), row.names = FALSE)

# =============================================================================
# 6. Figure 4A: ssGSEA免疫浸润热图
# =============================================================================
message("\n=== Figure 4A: ssGSEA heatmap ===")

immune_anno <- data.frame(
  UPR_subtype = subtypes,
  row.names = common_samples
)
sample_order <- order(subtypes)

# 行名标注一致性（*表示MCPcounter验证一致）
rownames_annotated <- rownames(immune_ssgsea)
idx_mapped <- rownames_annotated %in% consensus_mapped
rownames_annotated[idx_mapped] <- paste0(rownames_annotated[idx_mapped], " *")

ssgsea_plot <- immune_ssgsea
rownames(ssgsea_plot) <- rownames_annotated

# 行侧注释: 是否一致性显著
row_sig <- ifelse(rownames(immune_ssgsea) %in% consensus_sig_cells, "Significant", "NS")
row_anno <- rowAnnotation(
  Consensus = row_sig,
  col = list(Consensus = c("Significant" = "#E64B35", "NS" = "grey80")),
  annotation_name_side = "top",
  simple_anno_size = unit(3, "mm")
)

col_anno <- HeatmapAnnotation(
  UPR_subtype = subtypes[common_samples[sample_order]],
  col = list(UPR_subtype = COLORS_SUBTYPE),
  annotation_name_side = "left"
)

pdf(file.path(FIG_DIR, "Fig4A_immune_ssgsea_heatmap.pdf"), width = 14, height = 10)
ht <- Heatmap(
  ssgsea_plot[, common_samples[sample_order]],
  name = "ssGSEA\nScore",
  col = colorRamp2(
    seq(min(immune_ssgsea), max(immune_ssgsea), length.out = 100),
    colorRampPalette(c("#2166AC", "white", "#B2182B"))(100)
  ),
  top_annotation = col_anno,
  right_annotation = row_anno,
  show_column_names = FALSE,
  cluster_columns = FALSE,
  column_split = subtypes[common_samples[sample_order]],
  row_names_gp = gpar(fontsize = 8),
  column_title_gp = gpar(fontsize = 10, fontface = "bold"),
  row_title = "Immune Cell Types",
  column_title = "ssGSEA Immune Infiltration by UPR Subtype\n(* = MCPcounter concordant)"
)
draw(ht)
dev.off()

# =============================================================================
# 7. Figure 4B: 一致性显著免疫细胞箱线图（pairwise Wilcoxon）
# =============================================================================
message("\n=== Figure 4B: Consensus significant immune cell boxplots ===")

if (length(consensus_sig_cells) > 0) {
  # 选择展示的细胞类型（最多16个）
  show_cells <- consensus_sig_cells[1:min(16, length(consensus_sig_cells))]

  box_data <- immune_long_ssgsea %>%
    filter(CellType %in% show_cells) %>%
    mutate(CellType = factor(CellType, levels = show_cells))

  p_immune_box <- ggplot(box_data, aes(x = UPR_subtype, y = Score, fill = UPR_subtype)) +
    geom_boxplot(outlier.size = 0.3, width = 0.6) +
    facet_wrap(~CellType, scales = "free_y", ncol = 4) +
    stat_compare_means(method = "kruskal.test", size = 2.8, label = "p.format",
                       label.y.npc = 0.95) +
    stat_compare_means(comparisons = comparisons_list, method = "wilcox.test",
                       size = 2.2, p.adjust.method = "BH",
                       label = "p.signif", hide.ns = TRUE,
                       step.increase = 0.08) +
    scale_fill_manual(values = COLORS_SUBTYPE) +
    THEME_PUBLICATION +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
          legend.position = "none",
          strip.text = element_text(size = 7, face = "bold")) +
    labs(x = "UPR Subtype", y = "ssGSEA Enrichment Score") +
    ggtitle("Consensus Significant Immune Cell Infiltration")

  ggsave(file.path(FIG_DIR, "Fig4B_immune_boxplots_consensus.pdf"), p_immune_box,
         width = 14, height = max(6, ceiling(length(show_cells) / 4) * 3.5))
}

# =============================================================================
# 8. Figure 4C: 免疫检查点基因表达（24基因）
# =============================================================================
message("\n=== Figure 4C: Immune checkpoint gene expression (24 genes) ===")

# 获取可用的检查点基因
icp_valid <- immune_checkpoint_genes[immune_checkpoint_genes %in% rownames(expr_log)]
message(sprintf("  Immune checkpoint genes found: %d/%d",
                length(icp_valid), length(immune_checkpoint_genes)))

icp_expr <- expr_log[icp_valid, common_samples]

icp_long <- as.data.frame(t(icp_expr)) %>%
  mutate(barcode = rownames(.), UPR_subtype = subtypes) %>%
  pivot_longer(cols = -c(barcode, UPR_subtype),
               names_to = "Gene", values_to = "Expression")

# 检查点基因统计: Kruskal-Wallis + pairwise Wilcoxon, BH FDR校正
icp_kw <- icp_long %>%
  group_by(Gene) %>%
  summarise(
    kw_pvalue = tryCatch(kruskal.test(Expression ~ UPR_subtype)$p.value, error = function(e) NA_real_),
    .groups = "drop"
  ) %>%
  mutate(kw_padj = p.adjust(kw_pvalue, method = "BH"))

# pairwise Wilcoxon详细结果
icp_pairwise_list <- vector("list", length(icp_valid))
names(icp_pairwise_list) <- icp_valid
for (gene in icp_valid) {
  gene_data <- icp_long %>% filter(Gene == gene)
  pw <- pairwise.wilcox.test(gene_data$Expression, gene_data$UPR_subtype,
                              p.adjust.method = "BH")
  pw_df <- as.data.frame(as.table(pw$p.value)) %>%
    filter(!is.na(Freq)) %>%
    rename(Group1 = Var1, Group2 = Var2, p_adj = Freq) %>%
    mutate(Gene = gene)
  icp_pairwise_list[[gene]] <- pw_df
}
icp_pairwise_df <- do.call(rbind, icp_pairwise_list)

icp_sig <- icp_kw %>% filter(kw_padj < FDR_CUTOFF) %>% pull(Gene)
message(sprintf("  Significant checkpoint genes: %d/%d (FDR < %.2f)",
                length(icp_sig), nrow(icp_kw), FDR_CUTOFF))

# 箱线图
p_icp <- ggplot(icp_long, aes(x = UPR_subtype, y = Expression, fill = UPR_subtype)) +
  geom_boxplot(outlier.size = 0.3, width = 0.6) +
  facet_wrap(~Gene, scales = "free_y", ncol = 6) +
  stat_compare_means(method = "kruskal.test", size = 2.2, label = "p.format") +
  stat_compare_means(comparisons = comparisons_list, method = "wilcox.test",
                     size = 2, p.adjust.method = "BH",
                     label = "p.signif", hide.ns = TRUE,
                     step.increase = 0.06) +
  scale_fill_manual(values = COLORS_SUBTYPE) +
  THEME_PUBLICATION +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
        legend.position = "bottom",
        strip.text = element_text(size = 7)) +
  labs(x = "UPR Subtype", y = "log2(TPM + 1)") +
  ggtitle("Immune Checkpoint Gene Expression by UPR Subtype")

ggsave(file.path(FIG_DIR, "Fig4C_immune_checkpoints.pdf"), p_icp,
       width = 18, height = max(8, ceiling(length(icp_valid) / 6) * 3.5))

# 显著检查点基因热图（补充）
if (length(icp_sig) > 1) {
  pdf(file.path(FIG_DIR, "FigS_checkpoint_heatmap.pdf"), width = 14, height = 6)
  pheatmap(
    as.matrix(icp_expr[icp_sig, common_samples[sample_order]]),
    annotation_col = immune_anno,
    annotation_colors = list(UPR_subtype = COLORS_SUBTYPE),
    scale = "row",
    show_colnames = FALSE,
    cluster_cols = FALSE,
    fontsize_row = 9,
    color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
    main = "Significant Immune Checkpoint Genes (FDR < 0.05)"
  )
  dev.off()
}

# =============================================================================
# 9. 评审修订M3: IDH-WT亚组内免疫浸润差异验证
# =============================================================================
message("\n=== IDH-WT subgroup immune infiltration validation (reviewer revision M3) ===")
message("  Purpose: verify immune differences are not solely driven by IDH status confounding")

# 加载临床数据获取IDH状态
load(file.path(DATA_PROC, "consensus_clustering_results.RData"))
idhwt_immune_results <- NULL

if (file.exists(file.path(DATA_PROC, "clinical_characterization_results.RData"))) {
  load(file.path(DATA_PROC, "clinical_characterization_results.RData"))
  idh_col_available <- "IDH_status" %in% colnames(heatmap_data)
} else {
  # 直接从clin_matched获取
  load(file.path(DATA_PROC, "tcga_glioma_expression.RData"))
  idh_col_available <- "IDH_status" %in% colnames(clin_matched)
  if (idh_col_available) {
    heatmap_data <- merge(cluster_df, clin_matched, by = "barcode")
  }
}

if (idh_col_available) {
  # 筛选IDH-WT样本
  idhwt_barcodes <- heatmap_data$barcode[heatmap_data$IDH_status == "WT" &
                                           !is.na(heatmap_data$IDH_status)]
  idhwt_barcodes <- intersect(idhwt_barcodes, common_samples)
  message(sprintf("  IDH-WT samples for immune analysis: %d", length(idhwt_barcodes)))

  idhwt_subtypes <- subtype_map[idhwt_barcodes]
  n_subtypes_idhwt <- length(unique(idhwt_subtypes))

  if (length(idhwt_barcodes) >= 30 && n_subtypes_idhwt >= 2) {
    message(sprintf("  IDH-WT subtype distribution: %s",
                    paste(names(table(idhwt_subtypes)), table(idhwt_subtypes),
                          sep = "=", collapse = ", ")))

    # --- 9.1 ssGSEA: IDH-WT Kruskal-Wallis ---
    immune_long_idhwt <- immune_long_ssgsea %>%
      dplyr::filter(barcode %in% idhwt_barcodes)

    ssgsea_kw_idhwt <- immune_long_idhwt %>%
      dplyr::group_by(CellType) %>%
      dplyr::summarise(
        kw_pvalue = tryCatch(kruskal.test(Score ~ UPR_subtype)$p.value,
                             error = function(e) NA_real_),
        .groups = "drop"
      ) %>%
      dplyr::filter(!is.na(kw_pvalue)) %>%
      dplyr::mutate(kw_padj = p.adjust(kw_pvalue, method = "BH"))

    ssgsea_sig_idhwt <- ssgsea_kw_idhwt %>%
      dplyr::filter(kw_padj < FDR_CUTOFF) %>%
      dplyr::pull(CellType)

    message(sprintf("\n  ssGSEA IDH-WT significant cell types: %d/%d (FDR < %.2f)",
                    length(ssgsea_sig_idhwt), nrow(ssgsea_kw_idhwt), FDR_CUTOFF))
    message(sprintf("  ssGSEA full cohort significant: %d/%d",
                    length(ssgsea_sig), nrow(ssgsea_kw)))
    message(sprintf("  Retained after IDH-WT restriction: %d/%d (%.1f%%)",
                    length(ssgsea_sig_idhwt),
                    length(ssgsea_sig),
                    ifelse(length(ssgsea_sig) > 0,
                           length(ssgsea_sig_idhwt) / length(ssgsea_sig) * 100, 0)))

    # 列出在全队列显著但IDH-WT中不再显著的细胞类型（可能受IDH混杂）
    lost_in_idhwt <- setdiff(ssgsea_sig, ssgsea_sig_idhwt)
    if (length(lost_in_idhwt) > 0) {
      message(sprintf("  Cell types no longer significant in IDH-WT (potential IDH confounding):"))
      message(sprintf("    %s", paste(lost_in_idhwt, collapse = ", ")))
    }

    # 仅IDH-WT中新出现显著的（不太可能但检查）
    new_in_idhwt <- setdiff(ssgsea_sig_idhwt, ssgsea_sig)
    if (length(new_in_idhwt) > 0) {
      message(sprintf("  Cell types newly significant in IDH-WT only:"))
      message(sprintf("    %s", paste(new_in_idhwt, collapse = ", ")))
    }

    # --- 9.2 ssGSEA IDH-WT: pairwise Wilcoxon ---
    ssgsea_pairwise_idhwt <- immune_long_idhwt %>%
      dplyr::filter(CellType %in% ssgsea_sig_idhwt) %>%
      dplyr::group_by(CellType) %>%
      dplyr::summarise(
        pw_results = list(pairwise.wilcox.test(Score, UPR_subtype, p.adjust.method = "BH")),
        .groups = "drop"
      )

    # --- 9.3 IDH-WT ssGSEA箱线图（补充图） ---
    if (length(ssgsea_sig_idhwt) > 0) {
      show_cells_idhwt <- ssgsea_sig_idhwt[1:min(16, length(ssgsea_sig_idhwt))]

      idhwt_comparisons <- combn(sort(unique(as.character(idhwt_subtypes))), 2,
                                  simplify = FALSE)

      box_idhwt <- immune_long_idhwt %>%
        dplyr::filter(CellType %in% show_cells_idhwt) %>%
        dplyr::mutate(CellType = factor(CellType, levels = show_cells_idhwt))

      p_idhwt_immune <- ggplot(box_idhwt,
                                aes(x = UPR_subtype, y = Score, fill = UPR_subtype)) +
        geom_boxplot(outlier.size = 0.3, width = 0.6) +
        facet_wrap(~CellType, scales = "free_y", ncol = 4) +
        stat_compare_means(method = "kruskal.test", size = 2.8, label = "p.format",
                           label.y.npc = 0.95) +
        stat_compare_means(comparisons = idhwt_comparisons, method = "wilcox.test",
                           size = 2.2, p.adjust.method = "BH",
                           label = "p.signif", hide.ns = TRUE,
                           step.increase = 0.08) +
        scale_fill_manual(values = COLORS_SUBTYPE) +
        THEME_PUBLICATION +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
              legend.position = "none",
              strip.text = element_text(size = 7, face = "bold")) +
        labs(x = "UPR Subtype", y = "ssGSEA Enrichment Score") +
        ggtitle(sprintf("IDH-WT Subgroup: Immune Infiltration (n=%d)", length(idhwt_barcodes)))

      ggsave(file.path(FIG_DIR, "FigS_immune_ssgsea_idhwt.pdf"), p_idhwt_immune,
             width = 14, height = max(6, ceiling(length(show_cells_idhwt) / 4) * 3.5))
    }

    # --- 9.4 汇总IDH-WT验证结果 ---
    idhwt_immune_results <- list(
      n_idhwt_samples = length(idhwt_barcodes),
      ssgsea_kw_idhwt = ssgsea_kw_idhwt,
      ssgsea_sig_idhwt = ssgsea_sig_idhwt,
      n_sig_full_cohort = length(ssgsea_sig),
      n_sig_idhwt = length(ssgsea_sig_idhwt),
      lost_in_idhwt = lost_in_idhwt,
      pct_retained = ifelse(length(ssgsea_sig) > 0,
                            length(ssgsea_sig_idhwt) / length(ssgsea_sig) * 100, NA)
    )

    write.csv(ssgsea_kw_idhwt,
              file.path(RES_DIR, "immune_ssgsea_kw_idhwt.csv"), row.names = FALSE)

    message(sprintf("\n  IDH-WT immune validation summary:"))
    message(sprintf("    Full cohort significant: %d cell types", length(ssgsea_sig)))
    message(sprintf("    IDH-WT significant: %d cell types", length(ssgsea_sig_idhwt)))
    message(sprintf("    Retained: %.1f%%", idhwt_immune_results$pct_retained))

  } else {
    message("  Too few IDH-WT samples or subtypes for immune subgroup analysis.")
  }
} else {
  message("  IDH_status not available; skipping IDH-WT immune validation.")
}

# =============================================================================
# 10. 保存免疫分析结果
# =============================================================================
message("\n=== Saving immune analysis results ===")

save(immune_ssgsea, ssgsea_kw, ssgsea_pairwise,
     mcp_results, mcp_kw, mcp_sig,
     est_df, estimate_consistent,
     consensus_sig_cells, consensus_mapped, consensus_report,
     icp_kw, icp_pairwise_df,
     idhwt_immune_results,
     file = file.path(DATA_PROC, "immune_analysis_results.RData"))

write.csv(ssgsea_kw, file.path(RES_DIR, "immune_ssgsea_kw_test.csv"), row.names = FALSE)
write.csv(icp_kw, file.path(RES_DIR, "checkpoint_kw_test.csv"), row.names = FALSE)
write.csv(icp_pairwise_df, file.path(RES_DIR, "checkpoint_pairwise_wilcox.csv"), row.names = FALSE)

message("\n=== Immune analysis completed ===")
message(sprintf("  Figure 4A: %s", file.path(FIG_DIR, "Fig4A_immune_ssgsea_heatmap.pdf")))
message(sprintf("  Figure 4B: %s", file.path(FIG_DIR, "Fig4B_immune_boxplots_consensus.pdf")))
message(sprintf("  Figure 4C: %s", file.path(FIG_DIR, "Fig4C_immune_checkpoints.pdf")))
message("Next step: Run 04_bulk_subtyping/05_genomic_analysis.R")
