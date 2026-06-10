###############################################################################
# 04_differential_analysis.R
# UPR-high vs UPR-low 免疫细胞差异基因与功能富集
#
# 评审修订要点:
#   1. |log2FC| > 0.5 (SC_LOGFC_CUTOFF), adj.p < 0.05
#   2. Pseudobulk DESeq2 敏感性验证 (Squair et al. 2021)
#   3. TAM M1/M2 极化、CD8+ T 耗竭、DC 抗原呈递重点分析
#   4. fgsea + Hallmark + GO:BP GSEA
#   5. 细胞数 < 50 标注 exploratory
###############################################################################

source("00_setup/config.R")

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(EnhancedVolcano)
library(fgsea)
library(msigdbr)
library(DESeq2)
library(ComplexHeatmap)
library(circlize)

set.seed(SEED)

# =============================================================================
# 1. 加载数据与基因集
# =============================================================================
message("=== Loading UPR-scored data and gene sets ===")
seu <- readRDS(file.path(DATA_PROC, "seu_upr_scored.rds"))
load(file.path(DATA_PROC, "upr_gene_sets.RData"))

# =============================================================================
# 2. 准备 GSEA 基因集 (fgsea 格式: named list)
# =============================================================================
message("=== Preparing gene sets for fgsea ===")

# Hallmark
hallmark_df <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_gs <- split(hallmark_df$gene_symbol, hallmark_df$gs_name)

# GO:BP
gobp_df <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
gobp_gs <- split(gobp_df$gene_symbol, gobp_df$gs_name)

# =============================================================================
# 3. 各免疫细胞类型 UPR-high vs UPR-low 差异分析
# =============================================================================
message("=== Differential expression: UPR-high vs UPR-low per immune cell type ===")

immune_types <- c("TAM", "Microglia", "CD4 T", "CD8 T", "Treg", "DC", "MDSC")
immune_types <- immune_types[immune_types %in% unique(seu$celltype)]

deg_results   <- list()
gsea_results  <- list()
celltype_info <- list()   # track cell counts & exploratory flag

for (ct in immune_types) {
  message(sprintf("\n--- Processing: %s ---", ct))

  cells_ct <- subset(seu, celltype == ct)
  n_total  <- ncol(cells_ct)

  # Exploratory flag
  is_exploratory <- n_total < SC_MIN_CELLS_EXPLORATORY
  celltype_info[[ct]] <- list(n = n_total, exploratory = is_exploratory)

  if (is_exploratory) {
    message(sprintf("  [EXPLORATORY] %s has only %d cells (< %d). Results will be flagged.",
                    ct, n_total, SC_MIN_CELLS_EXPLORATORY))
  }

  if (n_total < 20) {
    message(sprintf("  Skipping %s: only %d cells — too few for any DE.", ct, n_total))
    next
  }

  # --- UPR-high / UPR-low split within this cell type ---
  median_upr <- median(cells_ct$UPR_score, na.rm = TRUE)
  cells_ct$UPR_group_ct <- ifelse(cells_ct$UPR_score > median_upr, "UPR_high", "UPR_low")

  n_high <- sum(cells_ct$UPR_group_ct == "UPR_high")
  n_low  <- sum(cells_ct$UPR_group_ct == "UPR_low")
  message(sprintf("  UPR-high: %d cells, UPR-low: %d cells", n_high, n_low))

  if (min(n_high, n_low) < 10) {
    message(sprintf("  Skipping %s: insufficient cells in one group (min = %d).",
                    ct, min(n_high, n_low)))
    next
  }

  # --- FindMarkers (Wilcoxon) ---
  Idents(cells_ct) <- "UPR_group_ct"
  deg <- FindMarkers(cells_ct,
                     ident.1     = "UPR_high",
                     ident.2     = "UPR_low",
                     test.use    = "wilcox",
                     logfc.threshold = 0.1,   # low threshold to capture full distribution
                     min.pct     = 0.1)

  deg$gene     <- rownames(deg)
  deg$celltype <- ct
  # Reviewer cutoff: |log2FC| > SC_LOGFC_CUTOFF (0.5) AND adj.p < 0.05
  deg$significant <- (deg$p_val_adj < FDR_CUTOFF) & (abs(deg$avg_log2FC) > SC_LOGFC_CUTOFF)
  deg$exploratory <- is_exploratory

  deg_results[[ct]] <- deg

  n_up   <- sum(deg$significant & deg$avg_log2FC > 0, na.rm = TRUE)
  n_down <- sum(deg$significant & deg$avg_log2FC < 0, na.rm = TRUE)
  message(sprintf("  DEGs (|log2FC|>%.1f, adj.p<%.2f): %d up, %d down%s",
                  SC_LOGFC_CUTOFF, FDR_CUTOFF, n_up, n_down,
                  ifelse(is_exploratory, " [EXPLORATORY]", "")))

  # --- Volcano plot ---
  volcano_title <- paste0(ct, ": UPR-high vs UPR-low",
                          ifelse(is_exploratory, " [EXPLORATORY]", ""))

  p_volcano <- EnhancedVolcano(deg,
    lab       = deg$gene,
    x         = "avg_log2FC",
    y         = "p_val_adj",
    title     = volcano_title,
    subtitle  = sprintf("n_high=%d, n_low=%d", n_high, n_low),
    pCutoff   = FDR_CUTOFF,
    FCcutoff  = SC_LOGFC_CUTOFF,
    pointSize = 1.5,
    labSize   = 3,
    drawConnectors = TRUE,
    max.overlaps   = 20)

  ggsave(file.path(FIG_DIR, paste0("deg_volcano_", gsub(" ", "_", ct), ".pdf")),
         p_volcano, width = 9, height = 7)

  # --- fgsea: Hallmark ---
  gene_ranks <- setNames(deg$avg_log2FC, deg$gene)
  gene_ranks <- sort(gene_ranks, decreasing = TRUE)

  fgsea_hallmark <- fgsea(pathways = hallmark_gs,
                          stats    = gene_ranks,
                          minSize  = 15,
                          maxSize  = 500,
                          nPermSimple = 10000)
  fgsea_hallmark$celltype <- ct
  gsea_results[[paste0(ct, "_hallmark")]] <- as.data.frame(fgsea_hallmark)

  # Plot top enriched Hallmark pathways
  sig_hallmark <- fgsea_hallmark[fgsea_hallmark$padj < 0.05, ]
  if (nrow(sig_hallmark) > 0) {
    sig_hallmark <- sig_hallmark[order(sig_hallmark$NES, decreasing = TRUE), ]
    top_n <- min(20, nrow(sig_hallmark))
    plot_data <- sig_hallmark[1:top_n, ]
    plot_data$pathway_short <- gsub("HALLMARK_", "", plot_data$pathway)
    plot_data$pathway_short <- gsub("_", " ", plot_data$pathway_short)
    plot_data$pathway_short <- factor(plot_data$pathway_short,
                                      levels = rev(plot_data$pathway_short))

    p_gsea_bar <- ggplot(plot_data, aes(x = NES, y = pathway_short, fill = NES > 0)) +
      geom_col() +
      scale_fill_manual(values = c("TRUE" = "#E64B35", "FALSE" = "#4DBBD5"),
                        labels = c("Downregulated", "Upregulated"), name = "") +
      labs(x = "Normalized Enrichment Score (NES)", y = "",
           title = paste0(ct, ": Hallmark GSEA (UPR-high vs UPR-low)")) +
      THEME_PUBLICATION

    ggsave(file.path(FIG_DIR, paste0("gsea_hallmark_", gsub(" ", "_", ct), ".pdf")),
           p_gsea_bar, width = 10, height = max(4, 0.3 * top_n + 2))
  }

  # --- fgsea: GO:BP ---
  fgsea_gobp <- fgsea(pathways = gobp_gs,
                      stats    = gene_ranks,
                      minSize  = 15,
                      maxSize  = 500,
                      nPermSimple = 10000)
  fgsea_gobp$celltype <- ct
  gsea_results[[paste0(ct, "_GOBP")]] <- as.data.frame(fgsea_gobp)
}

# --- Combine and save DEG results ---
all_degs <- do.call(rbind, deg_results)
write.csv(all_degs,
          file.path(RES_DIR, "deg_upr_high_vs_low_all_celltypes.csv"),
          row.names = FALSE)

# --- Save cell type info ---
celltype_info_df <- data.frame(
  celltype    = names(celltype_info),
  n_cells     = sapply(celltype_info, function(x) x$n),
  exploratory = sapply(celltype_info, function(x) x$exploratory),
  row.names   = NULL
)
write.csv(celltype_info_df,
          file.path(RES_DIR, "deg_celltype_info.csv"),
          row.names = FALSE)

# --- Save GSEA results (convert leadingEdge list column to character) ---
gsea_save <- lapply(gsea_results, function(df) {
  df$leadingEdge <- sapply(df$leadingEdge, function(x) paste(x, collapse = ";"))
  df
})
save(gsea_results, file = file.path(RES_DIR, "gsea_fgsea_results_by_celltype.RData"))
all_gsea <- do.call(rbind, gsea_save)
write.csv(all_gsea, file.path(RES_DIR, "gsea_fgsea_all_celltypes.csv"), row.names = FALSE)

# =============================================================================
# 4. Pseudobulk DESeq2 敏感性验证 (Squair et al. 2021)
#    按 sample x cell_type x UPR_group 聚合 counts
# =============================================================================
message("\n=== Pseudobulk DESeq2 sensitivity validation ===")

# Determine sample ID column
sample_col <- NULL
for (candidate in c("sample", "orig.ident", "patient", "Sample", "sample_id")) {
  if (candidate %in% colnames(seu@meta.data)) {
    sample_col <- candidate
    break
  }
}

if (is.null(sample_col)) {
  message("WARNING: No sample-level column found in metadata. ",
          "Using 'orig.ident' as fallback.")
  sample_col <- "orig.ident"
}
message(sprintf("  Using '%s' as sample identifier.", sample_col))

pseudobulk_results <- list()

# Run pseudobulk for TAM and CD8 T (at minimum, per reviewer request)
pb_celltypes <- intersect(c("TAM", "CD8 T"), names(deg_results))

for (ct in pb_celltypes) {
  message(sprintf("\n--- Pseudobulk DESeq2: %s ---", ct))

  cells_ct <- subset(seu, celltype == ct)

  # Assign UPR group within cell type (same as above)
  median_upr <- median(cells_ct$UPR_score, na.rm = TRUE)
  cells_ct$UPR_group_ct <- ifelse(cells_ct$UPR_score > median_upr, "UPR_high", "UPR_low")

  # Create pseudobulk: aggregate raw counts by sample x UPR_group
  cells_ct$pb_group <- paste0(cells_ct@meta.data[[sample_col]], "__", cells_ct$UPR_group_ct)

  pb_groups <- unique(cells_ct$pb_group)
  message(sprintf("  Aggregating %d cells into %d pseudobulk samples.", ncol(cells_ct), length(pb_groups)))

  # Need >=2 samples per UPR group for DESeq2
  pb_meta <- data.frame(
    pb_id     = pb_groups,
    sample    = sub("__.*", "", pb_groups),
    upr_group = sub(".*__", "", pb_groups),
    stringsAsFactors = FALSE
  )

  n_high_samples <- sum(pb_meta$upr_group == "UPR_high")
  n_low_samples  <- sum(pb_meta$upr_group == "UPR_low")

  if (n_high_samples < 2 || n_low_samples < 2) {
    message(sprintf("  Insufficient biological replicates for DESeq2 (high=%d, low=%d).",
                    n_high_samples, n_low_samples))

    # For single-sample data: create pseudo-replicates by random splitting
    # NOTE: This is a sensitivity check only; results should be interpreted with caution
    message("  Creating pseudo-replicates for sensitivity check (3 splits per group)...")

    n_reps <- 3
    new_pb_groups <- character(ncol(cells_ct))

    for (grp in c("UPR_high", "UPR_low")) {
      grp_idx <- which(cells_ct$UPR_group_ct == grp)
      if (length(grp_idx) < n_reps * 10) {
        message(sprintf("  Not enough cells in %s-%s for pseudo-replicates.", ct, grp))
        next
      }
      set.seed(SEED)
      rep_assignment <- sample(rep(seq_len(n_reps), length.out = length(grp_idx)))
      for (r in seq_len(n_reps)) {
        new_pb_groups[grp_idx[rep_assignment == r]] <- paste0("pseudorep_", grp, "_", r)
      }
    }

    cells_ct$pb_group <- new_pb_groups

    # Remove unassigned cells
    keep <- nchar(cells_ct$pb_group) > 0
    if (sum(keep) < 100) {
      message(sprintf("  Too few assigned cells. Skipping pseudobulk for %s.", ct))
      next
    }
    cells_ct <- cells_ct[, keep]

    pb_groups <- unique(cells_ct$pb_group)
    pb_meta <- data.frame(
      pb_id     = pb_groups,
      sample    = pb_groups,
      upr_group = ifelse(grepl("UPR_high", pb_groups), "UPR_high", "UPR_low"),
      stringsAsFactors = FALSE
    )

    n_high_samples <- sum(pb_meta$upr_group == "UPR_high")
    n_low_samples  <- sum(pb_meta$upr_group == "UPR_low")

    if (n_high_samples < 2 || n_low_samples < 2) {
      message(sprintf("  Still insufficient samples. Skipping pseudobulk for %s.", ct))
      next
    }
  }

  # Aggregate counts using RNA assay raw counts
  DefaultAssay(cells_ct) <- "RNA"
  raw_counts <- tryCatch(GetAssayData(cells_ct, layer = "counts"), error = function(e) GetAssayData(cells_ct, slot = "counts"))

  pb_counts <- matrix(0, nrow = nrow(raw_counts), ncol = length(pb_groups),
                       dimnames = list(rownames(raw_counts), pb_groups))

  for (pg in pb_groups) {
    cell_ids <- colnames(cells_ct)[cells_ct$pb_group == pg]
    if (length(cell_ids) == 1) {
      pb_counts[, pg] <- as.numeric(raw_counts[, cell_ids])
    } else {
      pb_counts[, pg] <- Matrix::rowSums(raw_counts[, cell_ids, drop = FALSE])
    }
  }

  # Filter low-expression genes
  keep_genes <- rowSums(pb_counts >= 5) >= 2
  pb_counts  <- pb_counts[keep_genes, ]
  message(sprintf("  Genes after filtering: %d", nrow(pb_counts)))

  if (nrow(pb_counts) < 100) {
    message(sprintf("  Too few genes for DESeq2 in %s. Skipping.", ct))
    next
  }

  # DESeq2
  rownames(pb_meta) <- pb_meta$pb_id
  pb_meta$upr_group <- factor(pb_meta$upr_group, levels = c("UPR_low", "UPR_high"))

  tryCatch({
    dds <- DESeqDataSetFromMatrix(countData = round(pb_counts[, pb_meta$pb_id]),
                                  colData   = pb_meta,
                                  design    = ~ upr_group)
    dds <- DESeq(dds)
    pb_res <- results(dds, contrast = c("upr_group", "UPR_high", "UPR_low"))
    pb_res_df <- as.data.frame(pb_res)
    pb_res_df$gene <- rownames(pb_res_df)
    pb_res_df$celltype <- ct
    pb_res_df$significant <- (!is.na(pb_res_df$padj)) &
                             (pb_res_df$padj < FDR_CUTOFF) &
                             (abs(pb_res_df$log2FoldChange) > SC_LOGFC_CUTOFF)

    pseudobulk_results[[ct]] <- pb_res_df

    n_sig <- sum(pb_res_df$significant, na.rm = TRUE)
    message(sprintf("  Pseudobulk DEGs: %d significant genes.", n_sig))

    # Compare with single-cell Wilcoxon results
    if (ct %in% names(deg_results)) {
      sc_sig_genes <- deg_results[[ct]]$gene[deg_results[[ct]]$significant]
      pb_sig_genes <- pb_res_df$gene[pb_res_df$significant]
      overlap <- length(intersect(sc_sig_genes, pb_sig_genes))
      message(sprintf("  Overlap with Wilcoxon DEGs: %d / SC=%d, PB=%d",
                      overlap, length(sc_sig_genes), length(pb_sig_genes)))

      # Concordance scatter plot: log2FC comparison
      common_genes <- intersect(deg_results[[ct]]$gene, pb_res_df$gene)
      if (length(common_genes) > 50) {
        concordance_df <- data.frame(
          gene      = common_genes,
          sc_log2FC = deg_results[[ct]]$avg_log2FC[match(common_genes, deg_results[[ct]]$gene)],
          pb_log2FC = pb_res_df$log2FoldChange[match(common_genes, pb_res_df$gene)],
          stringsAsFactors = FALSE
        )
        concordance_df <- concordance_df[complete.cases(concordance_df), ]

        r_val <- cor(concordance_df$sc_log2FC, concordance_df$pb_log2FC,
                     method = "spearman", use = "complete.obs")

        p_concordance <- ggplot(concordance_df, aes(x = sc_log2FC, y = pb_log2FC)) +
          geom_point(alpha = 0.3, size = 0.8) +
          geom_smooth(method = "lm", se = TRUE, color = "#E64B35") +
          geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
          annotate("text", x = min(concordance_df$sc_log2FC, na.rm = TRUE),
                   y = max(concordance_df$pb_log2FC, na.rm = TRUE),
                   label = sprintf("Spearman rho = %.3f", r_val),
                   hjust = 0, vjust = 1, size = 4) +
          labs(x = "Single-cell Wilcoxon log2FC",
               y = "Pseudobulk DESeq2 log2FC",
               title = paste0(ct, ": SC vs Pseudobulk Concordance")) +
          THEME_PUBLICATION

        ggsave(file.path(FIG_DIR, paste0("pseudobulk_concordance_", gsub(" ", "_", ct), ".pdf")),
               p_concordance, width = 6, height = 6)
      }
    }
  }, error = function(e) {
    message(sprintf("  DESeq2 error for %s: %s", ct, e$message))
  })
}

# Save pseudobulk results
if (length(pseudobulk_results) > 0) {
  all_pb <- do.call(rbind, pseudobulk_results)
  write.csv(all_pb,
            file.path(RES_DIR, "pseudobulk_deseq2_results.csv"),
            row.names = FALSE)
  save(pseudobulk_results,
       file = file.path(RES_DIR, "pseudobulk_deseq2_results.RData"))
}

# =============================================================================
# 5. Immune pathway NES heatmap (Hallmark across immune cell types)
# =============================================================================
message("\n=== Immune pathway NES heatmap ===")

immune_pathways <- c(
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "HALLMARK_COMPLEMENT",
  "HALLMARK_ALLOGRAFT_REJECTION",
  "HALLMARK_APOPTOSIS",
  "HALLMARK_HYPOXIA",
  "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
  "HALLMARK_MTORC1_SIGNALING",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_REACTIVE_OXYGEN_SPECIES_PATHWAY",
  "HALLMARK_ANGIOGENESIS"
)

# Cell types that have hallmark results
ct_with_hallmark <- sub("_hallmark$", "",
                        grep("_hallmark$", names(gsea_results), value = TRUE))

if (length(ct_with_hallmark) > 0) {
  nes_matrix  <- matrix(NA, nrow = length(immune_pathways), ncol = length(ct_with_hallmark),
                         dimnames = list(immune_pathways, ct_with_hallmark))
  pval_matrix <- nes_matrix

  for (ct in ct_with_hallmark) {
    res <- gsea_results[[paste0(ct, "_hallmark")]]
    for (pw in immune_pathways) {
      idx <- which(res$pathway == pw)
      if (length(idx) == 1) {
        nes_matrix[pw, ct]  <- res$NES[idx]
        pval_matrix[pw, ct] <- res$padj[idx]
      }
    }
  }

  # Simplify pathway names for display
  display_names <- gsub("HALLMARK_", "", rownames(nes_matrix))
  display_names <- gsub("_", " ", display_names)
  rownames(nes_matrix)  <- display_names
  rownames(pval_matrix) <- display_names

  # Significance stars
  sig_mark <- matrix("", nrow = nrow(pval_matrix), ncol = ncol(pval_matrix))
  sig_mark[!is.na(pval_matrix) & pval_matrix < 0.001] <- "***"
  sig_mark[!is.na(pval_matrix) & pval_matrix >= 0.001 & pval_matrix < 0.01]  <- "**"
  sig_mark[!is.na(pval_matrix) & pval_matrix >= 0.01  & pval_matrix < 0.05]  <- "*"

  # Add exploratory annotation to column names
  col_labels <- ct_with_hallmark
  for (i in seq_along(col_labels)) {
    if (col_labels[i] %in% names(celltype_info) && celltype_info[[col_labels[i]]]$exploratory) {
      col_labels[i] <- paste0(col_labels[i], "*")
    }
  }

  pdf(file.path(FIG_DIR, "Fig1_supp_immune_pathway_heatmap.pdf"), width = 10, height = 8)
  ht <- Heatmap(nes_matrix,
          name = "NES",
          col  = colorRamp2(c(-2, 0, 2), c("#4DBBD5", "white", "#E64B35")),
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid::grid.text(sig_mark[i, j], x, y, gp = grid::gpar(fontsize = 10))
          },
          cluster_rows    = TRUE,
          cluster_columns = FALSE,
          column_labels   = col_labels,
          row_names_gp    = grid::gpar(fontsize = 9),
          column_names_gp = grid::gpar(fontsize = 10),
          column_title    = "GSEA NES: UPR-high vs UPR-low (Hallmark)\n* = exploratory (n < 50 cells)",
          na_col = "grey90",
          heatmap_legend_param = list(title = "NES"))
  draw(ht)
  dev.off()
}

# =============================================================================
# 6. TAM M1/M2 极化分析
# =============================================================================
message("\n=== TAM M1/M2 polarization analysis ===")

tam_cells <- subset(seu, celltype == "TAM")

if (ncol(tam_cells) >= 30) {
  # Validate markers present in data
  m1_valid <- m1_markers[m1_markers %in% rownames(tam_cells)]
  m2_valid <- m2_markers[m2_markers %in% rownames(tam_cells)]

  message(sprintf("  M1 markers found: %d/%d (%s)",
                  length(m1_valid), length(m1_markers), paste(m1_valid, collapse = ", ")))
  message(sprintf("  M2 markers found: %d/%d (%s)",
                  length(m2_valid), length(m2_markers), paste(m2_valid, collapse = ", ")))

  if (length(m1_valid) >= 3 && length(m2_valid) >= 3) {
    tam_cells <- AddModuleScore(tam_cells, features = list(m1_valid), name = "M1_score")
    tam_cells <- AddModuleScore(tam_cells, features = list(m2_valid), name = "M2_score")

    # Re-assign UPR group within TAM
    median_upr <- median(tam_cells$UPR_score, na.rm = TRUE)
    tam_cells$UPR_group_tam <- ifelse(tam_cells$UPR_score > median_upr, "UPR-high", "UPR-low")

    # --- M1/M2 score comparison ---
    plot_df <- tam_cells@meta.data %>%
      dplyr::select(UPR_group_tam, M1_score1, M2_score1) %>%
      tidyr::pivot_longer(cols = c("M1_score1", "M2_score1"),
                          names_to = "Score_type", values_to = "Score") %>%
      dplyr::mutate(Score_type = ifelse(Score_type == "M1_score1", "M1 Score", "M2 Score"))

    p_m1m2 <- ggplot(plot_df, aes(x = UPR_group_tam, y = Score, fill = UPR_group_tam)) +
      geom_boxplot(outlier.size = 0.5, width = 0.6) +
      facet_wrap(~Score_type, scales = "free_y") +
      stat_compare_means(method = "wilcox.test", label = "p.format") +
      scale_fill_manual(values = c("UPR-high" = "#E64B35", "UPR-low" = "#4DBBD5")) +
      labs(x = "", y = "Module Score", fill = "UPR Group",
           title = "TAM M1/M2 Polarization: UPR-high vs UPR-low") +
      THEME_PUBLICATION

    ggsave(file.path(FIG_DIR, "Fig1_tam_m1m2_boxplot.pdf"), p_m1m2, width = 8, height = 5)

    # --- Individual M1/M2 marker gene expression ---
    all_m1m2 <- c(m1_valid, m2_valid)

    p_m1m2_dot <- DotPlot(tam_cells, features = all_m1m2, group.by = "UPR_group_tam") +
      RotatedAxis() +
      THEME_PUBLICATION +
      ggtitle("M1/M2 Markers in TAM (UPR-high vs UPR-low)")

    ggsave(file.path(FIG_DIR, "Fig1_tam_m1m2_dotplot.pdf"), p_m1m2_dot,
           width = max(6, length(all_m1m2) * 0.6 + 2), height = 4)

    # --- M1/M2 ratio ---
    tam_cells$M1M2_ratio <- tam_cells$M1_score1 - tam_cells$M2_score1

    p_ratio <- ggplot(tam_cells@meta.data, aes(x = UPR_group_tam, y = M1M2_ratio,
                                                fill = UPR_group_tam)) +
      geom_boxplot(outlier.size = 0.5, width = 0.6) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      stat_compare_means(method = "wilcox.test", label = "p.format") +
      scale_fill_manual(values = c("UPR-high" = "#E64B35", "UPR-low" = "#4DBBD5")) +
      labs(x = "", y = "M1 - M2 Score", fill = "UPR Group",
           title = "TAM M1/M2 Polarization Ratio") +
      THEME_PUBLICATION

    ggsave(file.path(FIG_DIR, "Fig1_tam_m1m2_ratio.pdf"), p_ratio, width = 5, height = 5)

    # Stats
    wt_m1 <- wilcox.test(M1_score1 ~ UPR_group_tam, data = tam_cells@meta.data)
    wt_m2 <- wilcox.test(M2_score1 ~ UPR_group_tam, data = tam_cells@meta.data)
    message(sprintf("  M1 score: Wilcoxon p = %.2e", wt_m1$p.value))
    message(sprintf("  M2 score: Wilcoxon p = %.2e", wt_m2$p.value))
  }
} else {
  message("  TAM cells < 30. Skipping M1/M2 analysis.")
}

# =============================================================================
# 7. CD8+ T 细胞耗竭分析
# =============================================================================
message("\n=== CD8+ T cell exhaustion analysis ===")

# Key exhaustion markers per reviewer requirement
exh_focus <- c("PDCD1", "HAVCR2", "LAG3", "TIGIT", "TOX", "ENTPD1")

cd8_cells <- subset(seu, celltype == "CD8 T")

if (ncol(cd8_cells) >= 20) {
  exh_valid <- t_exhaustion_genes[t_exhaustion_genes %in% rownames(cd8_cells)]
  exh_focus_valid <- exh_focus[exh_focus %in% rownames(cd8_cells)]

  message(sprintf("  CD8 T cells: %d", ncol(cd8_cells)))
  message(sprintf("  Exhaustion genes found: %d/%d", length(exh_valid), length(t_exhaustion_genes)))

  # UPR group within CD8
  median_upr <- median(cd8_cells$UPR_score, na.rm = TRUE)
  cd8_cells$UPR_group_cd8 <- ifelse(cd8_cells$UPR_score > median_upr, "UPR-high", "UPR-low")

  if (length(exh_valid) >= 3) {
    cd8_cells <- AddModuleScore(cd8_cells, features = list(exh_valid), name = "Exhaustion_score")

    # --- Violin plot (ggplot2替代VlnPlot避免S4兼容性问题) ---
    exh_df <- data.frame(score = cd8_cells$Exhaustion_score1, group = cd8_cells$UPR_group_cd8)
    p_exh_vln <- ggplot(exh_df, aes(x = group, y = score, fill = group)) +
      geom_violin(trim = FALSE) + geom_boxplot(width = 0.1, fill = "white") +
      scale_fill_manual(values = c("UPR-high" = "#E64B35", "UPR-low" = "#4DBBD5")) +
      stat_compare_means(method = "wilcox.test", label = "p.format") +
      labs(x = "", y = "Exhaustion Module Score",
           title = "CD8+ T Cell Exhaustion: UPR-high vs UPR-low") +
      THEME_PUBLICATION + theme(legend.position = "none")
    ggsave(file.path(FIG_DIR, "Fig1_cd8_exhaustion_violin.pdf"), p_exh_vln, width = 6, height = 5)

    wt_exh <- wilcox.test(Exhaustion_score1 ~ UPR_group_cd8, data = cd8_cells@meta.data)
    message(sprintf("  Exhaustion score: Wilcoxon p = %.2e", wt_exh$p.value))
  }

  # --- Individual exhaustion marker dot plot ---
  if (length(exh_focus_valid) > 0) {
    p_exh_dot <- DotPlot(cd8_cells, features = exh_focus_valid,
                          group.by = "UPR_group_cd8") +
      RotatedAxis() +
      THEME_PUBLICATION +
      ggtitle("Exhaustion Markers in CD8+ T Cells")

    ggsave(file.path(FIG_DIR, "Fig1_cd8_exh_markers_dotplot.pdf"), p_exh_dot,
           width = max(5, length(exh_focus_valid) * 0.8 + 2), height = 4)
  }

  # --- Boxplot per marker ---
  if (length(exh_focus_valid) >= 2) {
    exh_expr <- as.data.frame(t(as.matrix(
      tryCatch(GetAssayData(cd8_cells, layer = "data"), error = function(e) GetAssayData(cd8_cells, slot = "data"))[exh_focus_valid, , drop = FALSE]
    )))
    exh_expr$UPR_group <- cd8_cells$UPR_group_cd8

    exh_long <- exh_expr %>%
      tidyr::pivot_longer(cols = -UPR_group, names_to = "Gene", values_to = "Expression")

    p_exh_box <- ggplot(exh_long, aes(x = Gene, y = Expression, fill = UPR_group)) +
      geom_boxplot(outlier.size = 0.3, width = 0.6) +
      stat_compare_means(aes(group = UPR_group), method = "wilcox.test",
                         label = "p.signif", label.y.npc = 0.95) +
      scale_fill_manual(values = c("UPR-high" = "#E64B35", "UPR-low" = "#4DBBD5")) +
      labs(x = "", y = "Normalized Expression", fill = "UPR Group",
           title = "CD8+ T Exhaustion Markers: UPR-high vs UPR-low") +
      THEME_PUBLICATION +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    ggsave(file.path(FIG_DIR, "Fig1_cd8_exh_markers_boxplot.pdf"), p_exh_box,
           width = max(5, length(exh_focus_valid) * 0.9 + 2), height = 5)
  }
} else {
  message("  CD8 T cells < 20. Skipping exhaustion analysis.")
}

# =============================================================================
# 8. DC 抗原呈递分析 (HLA-II molecules)
# =============================================================================
message("\n=== DC antigen presentation analysis ===")

dc_cells <- subset(seu, celltype == "DC")

if (ncol(dc_cells) >= 20) {
  # HLA class II genes
  hla_ii_genes <- c("HLA-DRA", "HLA-DRB1", "HLA-DRB5", "HLA-DQA1", "HLA-DQB1",
                     "HLA-DPA1", "HLA-DPB1", "HLA-DMA", "HLA-DMB")
  hla_ii_valid <- hla_ii_genes[hla_ii_genes %in% rownames(dc_cells)]

  # Additional antigen presentation genes
  ap_genes <- c("CD80", "CD86", "CD40", "CD74", "CIITA", "B2M", "TAP1", "TAP2")
  ap_valid <- ap_genes[ap_genes %in% rownames(dc_cells)]

  is_dc_exploratory <- ncol(dc_cells) < SC_MIN_CELLS_EXPLORATORY

  message(sprintf("  DC cells: %d%s", ncol(dc_cells),
                  ifelse(is_dc_exploratory, " [EXPLORATORY]", "")))
  message(sprintf("  HLA-II genes found: %d/%d", length(hla_ii_valid), length(hla_ii_genes)))

  # UPR group within DC
  median_upr <- median(dc_cells$UPR_score, na.rm = TRUE)
  dc_cells$UPR_group_dc <- ifelse(dc_cells$UPR_score > median_upr, "UPR-high", "UPR-low")

  all_ap_genes <- unique(c(hla_ii_valid, ap_valid))

  if (length(all_ap_genes) >= 3) {
    # Module score for antigen presentation
    dc_cells <- AddModuleScore(dc_cells, features = list(all_ap_genes), name = "AP_score")

    dc_ap_df <- data.frame(score = dc_cells$AP_score1, group = dc_cells$UPR_group_dc)
    p_dc_ap <- ggplot(dc_ap_df, aes(x = group, y = score, fill = group)) +
      geom_violin(trim = FALSE) + geom_boxplot(width = 0.1, fill = "white") +
      scale_fill_manual(values = c("UPR-high" = "#E64B35", "UPR-low" = "#4DBBD5")) +
      stat_compare_means(method = "wilcox.test", label = "p.format") +
      labs(title = paste0("DC Antigen Presentation: UPR-high vs UPR-low",
                          ifelse(is_dc_exploratory, " [EXPLORATORY]", ""))) +
      THEME_PUBLICATION + theme(legend.position = "none")
    ggsave(file.path(FIG_DIR, "Fig1_dc_antigen_presentation_violin.pdf"),
           p_dc_ap, width = 6, height = 5)

    # Dot plot for individual HLA-II / AP genes
    p_dc_dot <- DotPlot(dc_cells, features = all_ap_genes,
                          group.by = "UPR_group_dc") +
      RotatedAxis() +
      THEME_PUBLICATION +
      ggtitle(paste0("HLA-II & AP Genes in DC",
                     ifelse(is_dc_exploratory, " [EXPLORATORY]", "")))

    ggsave(file.path(FIG_DIR, "Fig1_dc_hlaii_dotplot.pdf"), p_dc_dot,
           width = max(6, length(all_ap_genes) * 0.6 + 2), height = 4)

    wt_ap <- wilcox.test(AP_score1 ~ UPR_group_dc, data = dc_cells@meta.data)
    message(sprintf("  AP score: Wilcoxon p = %.2e", wt_ap$p.value))
  }
} else {
  message("  DC cells < 20. Skipping antigen presentation analysis.")
}

# =============================================================================
# 9. 综合结果摘要
# =============================================================================
message("\n=== Summary ===")
message("DEG results saved: ", file.path(RES_DIR, "deg_upr_high_vs_low_all_celltypes.csv"))
message("GSEA results saved: ", file.path(RES_DIR, "gsea_fgsea_all_celltypes.csv"))
if (length(pseudobulk_results) > 0) {
  message("Pseudobulk DESeq2 results: ", file.path(RES_DIR, "pseudobulk_deseq2_results.csv"))
}
message("Cell type info: ", file.path(RES_DIR, "deg_celltype_info.csv"))
message("Figures in: ", FIG_DIR)

message("\n=== Differential analysis completed ===")
message("Next step: Run 03_single_cell/05_cellchat.R")
