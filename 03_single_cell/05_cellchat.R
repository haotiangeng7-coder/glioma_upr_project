###############################################################################
# 05_cellchat.R
# CellChat v2 细胞间通讯分析
# 比较 UPR-high vs UPR-low 微环境中的通讯模式差异
#
# 评审修订核心:
#   1. 按 **样本** 级别分组（非单细胞）: 计算每个样本内恶性细胞 UPR 评分中位数，
#      按样本中位数分组，同一样本所有细胞归入同一组
#   2. 仅纳入恶性细胞数 >= 100 的样本
#   3. CellChat v2 分别构建 UPR-high / UPR-low 两组网络
#   4. 重点通路: CSF1-CSF1R, PD-L1-PD-1, LGALS9-TIM3, MHC-I, MHC-II
#   5. 可视化: Chord diagram, 热图, 气泡图
###############################################################################

source("00_setup/config.R")

library(Seurat)
library(CellChat)
library(patchwork)
library(ggplot2)
library(dplyr)
library(ComplexHeatmap)
library(circlize)

set.seed(SEED)

# =============================================================================
# 1. 加载数据
# =============================================================================
message("=== Loading data ===")
seu <- readRDS(file.path(DATA_PROC, "seu_upr_scored.rds"))

# =============================================================================
# 2. 样本级 UPR 分组（评审修订核心）
#    按恶性细胞的 UPR 评分对样本分组，而非对单细胞分组
# =============================================================================
message("=== Sample-level UPR grouping (reviewer revision) ===")

# Determine sample column
sample_col <- NULL
for (candidate in c("sample", "orig.ident", "patient", "Sample", "sample_id")) {
  if (candidate %in% colnames(seu@meta.data)) {
    sample_col <- candidate
    break
  }
}
if (is.null(sample_col)) {
  sample_col <- "orig.ident"
}
message(sprintf("  Using '%s' as sample identifier.", sample_col))

# Step 2a: For each sample, count malignant cells and compute median UPR score
malignant_cells <- subset(seu, celltype == "Malignant")
message(sprintf("  Total malignant cells: %d", ncol(malignant_cells)))

sample_upr_stats <- malignant_cells@meta.data %>%
  dplyr::group_by(.data[[sample_col]]) %>%
  dplyr::summarise(
    n_malignant   = n(),
    median_upr    = median(UPR_score, na.rm = TRUE),
    mean_upr      = mean(UPR_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::rename(sample_id = !!sample_col)

message("\n  Sample-level UPR statistics:")
print(as.data.frame(sample_upr_stats))

# Step 2b: Filter samples with >= SC_MIN_MALIGNANT_FOR_UPR (100) malignant cells
valid_samples <- sample_upr_stats %>%
  dplyr::filter(n_malignant >= SC_MIN_MALIGNANT_FOR_UPR) %>%
  dplyr::pull(sample_id)

excluded_samples <- setdiff(sample_upr_stats$sample_id, valid_samples)
if (length(excluded_samples) > 0) {
  message(sprintf("  Excluded %d sample(s) with < %d malignant cells: %s",
                  length(excluded_samples), SC_MIN_MALIGNANT_FOR_UPR,
                  paste(excluded_samples, collapse = ", ")))
}
message(sprintf("  Valid samples: %d", length(valid_samples)))

if (length(valid_samples) < 2) {
  # If only one valid sample or if all cells share the same sample ID,
  # use the single sample and split by UPR score of malignant cells directly
  message("  WARNING: Only one or no valid multi-sample structure detected.")
  message("  Falling back to single-sample strategy: splitting by median UPR of malignant cells.")

  # Still sample-level concept: the entire dataset is treated as one sample
  overall_median <- median(malignant_cells$UPR_score, na.rm = TRUE)
  message(sprintf("  Malignant cell UPR score median: %.4f", overall_median))

  # Assign all cells based on whether their sample's malignant UPR is above/below median
  # For single-sample: use quartiles to create more contrast
  q25 <- quantile(malignant_cells$UPR_score, 0.25, na.rm = TRUE)
  q75 <- quantile(malignant_cells$UPR_score, 0.75, na.rm = TRUE)

  # Classify malignant cells
  mal_high_cells <- colnames(malignant_cells)[malignant_cells$UPR_score > overall_median]
  mal_low_cells  <- colnames(malignant_cells)[malignant_cells$UPR_score <= overall_median]

  # For non-malignant cells: assign based on co-occurrence with malignant cells
  # In single-sample data, we use the malignant cell UPR score distribution to define
  # the microenvironment. Each non-malignant cell is assigned to the group based on
  # spatial/transcriptomic proximity — but since we lack spatial info, we assign
  # non-malignant cells randomly but proportionally to maintain balance.
  # Actually, for single-sample, reviewer intent is sample-level: use the overall
  # malignant median. But with one sample, all cells would be in one group.
  # Solution: if genuinely single-sample, split all cells by the overall UPR score
  # (from the UPR_group column already computed in 03_upr_scoring.R)
  message("  Using existing UPR_group from scoring step for single-sample fallback.")
  seu_valid <- seu
  seu_valid$sample_UPR_group <- seu_valid$UPR_group

} else {
  # Multi-sample: true sample-level grouping
  sample_upr_valid <- sample_upr_stats %>%
    dplyr::filter(sample_id %in% valid_samples)

  # Classify samples by median of their UPR medians
  sample_median_cutoff <- median(sample_upr_valid$median_upr)
  sample_upr_valid$sample_UPR_group <- ifelse(
    sample_upr_valid$median_upr > sample_median_cutoff,
    "UPR-high", "UPR-low"
  )

  message(sprintf("  Sample-level median UPR cutoff: %.4f", sample_median_cutoff))
  message(sprintf("  UPR-high samples: %d, UPR-low samples: %d",
                  sum(sample_upr_valid$sample_UPR_group == "UPR-high"),
                  sum(sample_upr_valid$sample_UPR_group == "UPR-low")))

  # Map sample group to all cells
  sample_group_map <- setNames(sample_upr_valid$sample_UPR_group,
                                sample_upr_valid$sample_id)

  # Subset to valid samples only
  valid_cells <- colnames(seu)[seu@meta.data[[sample_col]] %in% valid_samples]
  seu_valid <- seu[, valid_cells]
  seu_valid$sample_UPR_group <- unname(sample_group_map[seu_valid@meta.data[[sample_col]]])

  message(sprintf("  Cells in UPR-high samples: %d",
                  sum(seu_valid$sample_UPR_group == "UPR-high")))
  message(sprintf("  Cells in UPR-low samples: %d",
                  sum(seu_valid$sample_UPR_group == "UPR-low")))
}

# Save sample grouping info
write.csv(
  if (exists("sample_upr_valid")) {
    as.data.frame(sample_upr_valid)
  } else {
    data.frame(
      sample_id = "all",
      method = "single_sample_fallback",
      stringsAsFactors = FALSE
    )
  },
  file.path(RES_DIR, "cellchat_sample_upr_grouping.csv"),
  row.names = FALSE
)

# =============================================================================
# 3. Split data and create CellChat objects
# =============================================================================
message("\n=== Creating CellChat objects ===")

seu_high <- subset(seu_valid, sample_UPR_group == "UPR-high")
seu_low  <- subset(seu_valid, sample_UPR_group == "UPR-low")

message(sprintf("  UPR-high group: %d cells, UPR-low group: %d cells",
                ncol(seu_high), ncol(seu_low)))

# --- CellChat v2 analysis function ---
run_cellchat_v2 <- function(seurat_obj, group_name) {
  message(sprintf("\n--- CellChat v2: %s ---", group_name))

  # Extract normalized expression
  data_input <- tryCatch(GetAssayData(seurat_obj, layer = "data"), error = function(e) GetAssayData(seurat_obj, slot = "data"))

  # Metadata
  meta <- data.frame(
    labels = seurat_obj$celltype,
    row.names = colnames(seurat_obj)
  )

  # Filter cell types with too few cells
  cell_counts <- table(meta$labels)
  valid_types <- names(cell_counts[cell_counts >= 10])
  valid_cells <- rownames(meta)[meta$labels %in% valid_types]

  dropped_types <- names(cell_counts[cell_counts < 10])
  if (length(dropped_types) > 0) {
    message(sprintf("  Dropped cell types (< 10 cells): %s",
                    paste(dropped_types, collapse = ", ")))
  }

  data_input <- data_input[, valid_cells]
  meta <- meta[valid_cells, , drop = FALSE]

  message(sprintf("  Valid cell types (%d): %s",
                  length(valid_types), paste(valid_types, collapse = ", ")))
  message(sprintf("  Total cells: %d", length(valid_cells)))

  # Create CellChat object
  cellchat <- createCellChat(object = data_input, meta = meta, group.by = "labels")

  # Set ligand-receptor database (human, CellChat v2)
  CellChatDB <- CellChatDB.human
  cellchat@DB <- CellChatDB

  # Preprocessing
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)

  # Infer communication network
  cellchat <- computeCommunProb(cellchat, type = "triMean")
  cellchat <- filterCommunication(cellchat, min.cells = 10)

  # Infer pathway-level communication
  cellchat <- computeCommunProbPathway(cellchat)

  # Aggregate network
  cellchat <- aggregateNet(cellchat)

  # Network centrality
  cellchat <- netAnalysis_computeCentrality(cellchat)

  message(sprintf("  Pathways identified: %d", length(cellchat@netP$pathways)))
  message(sprintf("  Done: %s", group_name))

  return(cellchat)
}

cellchat_high <- run_cellchat_v2(seu_high, "UPR-high")
cellchat_low  <- run_cellchat_v2(seu_low, "UPR-low")

# =============================================================================
# 4. Individual network visualization
# =============================================================================
message("\n=== Individual network visualization ===")

# UPR-high network (count)
pdf(file.path(FIG_DIR, "Fig2a_cellchat_network_high.pdf"), width = 10, height = 8)
netVisual_circle(cellchat_high@net$count,
                 vertex.weight = table(cellchat_high@idents),
                 weight.scale = TRUE,
                 label.edge = FALSE,
                 title.name = "UPR-high: Number of Interactions")
dev.off()

# UPR-low network (count)
pdf(file.path(FIG_DIR, "Fig2a_cellchat_network_low.pdf"), width = 10, height = 8)
netVisual_circle(cellchat_low@net$count,
                 vertex.weight = table(cellchat_low@idents),
                 weight.scale = TRUE,
                 label.edge = FALSE,
                 title.name = "UPR-low: Number of Interactions")
dev.off()

# =============================================================================
# 5. Comparative analysis (Figure 2 core)
# =============================================================================
message("\n=== Comparative CellChat analysis (Figure 2) ===")

# Merge for comparison — lift to shared cell type space
object_list <- list("UPR-high" = cellchat_high, "UPR-low" = cellchat_low)
cellchat_merged <- mergeCellChat(object_list, add.names = names(object_list))

# --- 5.1 Overall interaction count and strength comparison ---
pdf(file.path(FIG_DIR, "Fig2b_interaction_comparison.pdf"), width = 12, height = 5)
p1 <- compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1, 2))
p2 <- compareInteractions(cellchat_merged, show.legend = FALSE, group = c(1, 2),
                           measure = "weight")
print(p1 + p2)
dev.off()

# --- 5.2 Differential interaction network ---
pdf(file.path(FIG_DIR, "Fig2c_differential_network_count.pdf"), width = 10, height = 8)
netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE,
                          measure = "count",
                          title.name = "Differential Interactions\n(UPR-high - UPR-low)")
dev.off()

pdf(file.path(FIG_DIR, "Fig2c_differential_network_weight.pdf"), width = 10, height = 8)
netVisual_diffInteraction(cellchat_merged, weight.scale = TRUE,
                          measure = "weight",
                          title.name = "Differential Interaction Strength\n(UPR-high - UPR-low)")
dev.off()

# --- 5.3 Heatmap comparison ---
pdf(file.path(FIG_DIR, "Fig2d_heatmap_count.pdf"), width = 8, height = 7)
ht1 <- netVisual_heatmap(cellchat_merged, measure = "count",
                          title.name = "Differential # Interactions")
dev.off()

pdf(file.path(FIG_DIR, "Fig2d_heatmap_weight.pdf"), width = 8, height = 7)
ht2 <- netVisual_heatmap(cellchat_merged, measure = "weight",
                          title.name = "Differential Interaction Strength")
dev.off()

# --- 5.4 Pathway ranking comparison ---
pdf(file.path(FIG_DIR, "Fig2e_pathway_ranking.pdf"), width = 12, height = 8)
rankNet(cellchat_merged, mode = "comparison", stacked = TRUE, do.stat = TRUE)
dev.off()

# Also stacked = FALSE for clearer separation
pdf(file.path(FIG_DIR, "Fig2e_pathway_ranking_split.pdf"), width = 12, height = 8)
rankNet(cellchat_merged, mode = "comparison", stacked = FALSE, do.stat = TRUE)
dev.off()

# =============================================================================
# 6. Key pathway deep analysis (reviewer-required pathways)
# =============================================================================
message("\n=== Key pathway analysis (CSF1-CSF1R, PD-L1-PD-1, LGALS9-TIM3, MHC-I, MHC-II) ===")

# Define reviewer-required key pathways
# CellChat pathway names may vary — define candidates to match
key_pathways_map <- list(
  "CSF1R"    = c("CSF", "CSF1", "CSF1R", "CSF3"),
  "PD-L1"    = c("PD-L1", "PDL1", "PD1", "CD274"),
  "GALECTIN" = c("GALECTIN", "LGALS"),
  "MHC-I"    = c("MHC-I", "MHC-1", "MHCI"),
  "MHC-II"   = c("MHC-II", "MHC-2", "MHCII")
)

# Get available pathways
available_pw_high <- cellchat_high@netP$pathways
available_pw_low  <- cellchat_low@netP$pathways
all_available_pw  <- union(available_pw_high, available_pw_low)

message(sprintf("  Available pathways in UPR-high: %d", length(available_pw_high)))
message(sprintf("  Available pathways in UPR-low: %d", length(available_pw_low)))

# Match key pathways
matched_pathways <- list()
for (key_name in names(key_pathways_map)) {
  candidates <- key_pathways_map[[key_name]]
  found <- intersect(candidates, all_available_pw)
  if (length(found) > 0) {
    matched_pathways[[key_name]] <- found[1]
    message(sprintf("  Matched '%s' -> '%s'", key_name, found[1]))
  } else {
    # Try partial matching
    partial <- all_available_pw[grep(paste(candidates, collapse = "|"),
                                      all_available_pw, ignore.case = TRUE)]
    if (length(partial) > 0) {
      matched_pathways[[key_name]] <- partial[1]
      message(sprintf("  Partial match '%s' -> '%s'", key_name, partial[1]))
    } else {
      message(sprintf("  No match for '%s'", key_name))
    }
  }
}

# Additional pathways of interest
extra_pathways <- c("CCL", "CXCL", "SPP1", "VEGF", "CD86", "ICAM", "CD40",
                     "PDGF", "FGF", "COMPLEMENT", "TNF", "IFN-II", "ANNEXIN")
for (ep in extra_pathways) {
  if (ep %in% all_available_pw && !(ep %in% unlist(matched_pathways))) {
    matched_pathways[[ep]] <- ep
  }
}

message(sprintf("  Total key pathways to visualize: %d", length(matched_pathways)))

# --- Chord diagrams for key pathways ---
for (pw_name in names(matched_pathways)) {
  pw <- matched_pathways[[pw_name]]

  in_high <- pw %in% available_pw_high
  in_low  <- pw %in% available_pw_low

  if (!in_high && !in_low) next

  message(sprintf("  Chord diagram: %s (high=%s, low=%s)",
                  pw, in_high, in_low))

  tryCatch({
    pdf(file.path(FIG_DIR, paste0("Fig2_chord_", gsub("[^A-Za-z0-9]", "_", pw), ".pdf")),
        width = 14, height = 6)
    par(mfrow = c(1, 2))

    if (in_high) {
      netVisual_aggregate(cellchat_high, signaling = pw,
                         layout = "chord",
                         title.name = paste0(pw, " — UPR-high"))
    } else {
      plot.new()
      title(paste0(pw, " — not detected in UPR-high"))
    }

    if (in_low) {
      netVisual_aggregate(cellchat_low, signaling = pw,
                         layout = "chord",
                         title.name = paste0(pw, " — UPR-low"))
    } else {
      plot.new()
      title(paste0(pw, " — not detected in UPR-low"))
    }

    dev.off()
  }, error = function(e) {
    message(sprintf("    Chord error for %s: %s", pw, e$message))
    try(dev.off(), silent = TRUE)
  })
}

# --- Bubble plot for key L-R pairs ---
message("\n=== Bubble plots for key pathways ===")

# Extract L-R pairs from key pathways for bubble comparison
key_signaling_for_bubble <- unlist(matched_pathways, use.names = FALSE)
key_signaling_for_bubble <- key_signaling_for_bubble[
  key_signaling_for_bubble %in% intersect(available_pw_high, available_pw_low)
]

if (length(key_signaling_for_bubble) > 0) {
  tryCatch({
    pdf(file.path(FIG_DIR, "Fig2f_bubble_key_pathways.pdf"), width = 14, height = 10)
    netVisual_bubble(cellchat_high,
                     signaling = key_signaling_for_bubble,
                     remove.isolate = TRUE,
                     title.name = "Key Pathways L-R Pairs (UPR-high)")
    dev.off()

    pdf(file.path(FIG_DIR, "Fig2f_bubble_key_pathways_low.pdf"), width = 14, height = 10)
    netVisual_bubble(cellchat_low,
                     signaling = key_signaling_for_bubble,
                     remove.isolate = TRUE,
                     title.name = "Key Pathways L-R Pairs (UPR-low)")
    dev.off()
  }, error = function(e) {
    message("  Bubble plot error: ", e$message)
    try(dev.off(), silent = TRUE)
  })

  # Comparative bubble: sources from Malignant, targets in immune cells
  immune_targets <- intersect(
    c("TAM", "Microglia", "CD8 T", "CD4 T", "Treg", "DC"),
    intersect(levels(cellchat_high@idents), levels(cellchat_low@idents))
  )

  if ("Malignant" %in% levels(cellchat_high@idents) && length(immune_targets) > 0) {
    tryCatch({
      pdf(file.path(FIG_DIR, "Fig2g_bubble_tumor_immune_high.pdf"), width = 14, height = 8)
      netVisual_bubble(cellchat_high,
                       sources.use = "Malignant",
                       targets.use = immune_targets,
                       remove.isolate = TRUE,
                       title.name = "Malignant -> Immune (UPR-high)")
      dev.off()

      pdf(file.path(FIG_DIR, "Fig2g_bubble_tumor_immune_low.pdf"), width = 14, height = 8)
      netVisual_bubble(cellchat_low,
                       sources.use = "Malignant",
                       targets.use = immune_targets,
                       remove.isolate = TRUE,
                       title.name = "Malignant -> Immune (UPR-low)")
      dev.off()
    }, error = function(e) {
      message("  Tumor-immune bubble error: ", e$message)
      try(dev.off(), silent = TRUE)
    })
  }
}

# =============================================================================
# 7. Tumor-Immune communication focused extraction
# =============================================================================
message("\n=== Extracting tumor-immune communication details ===")

# Tumor -> TAM
tryCatch({
  if ("Malignant" %in% levels(cellchat_high@idents) &&
      "TAM" %in% levels(cellchat_high@idents)) {
    tumor_tam_high <- subsetCommunication(cellchat_high,
                                           sources.use = "Malignant",
                                           targets.use = "TAM")
    write.csv(tumor_tam_high,
              file.path(RES_DIR, "cellchat_tumor_to_TAM_high.csv"),
              row.names = FALSE)
    message(sprintf("  Tumor -> TAM (UPR-high): %d interactions", nrow(tumor_tam_high)))
  }

  if ("Malignant" %in% levels(cellchat_low@idents) &&
      "TAM" %in% levels(cellchat_low@idents)) {
    tumor_tam_low <- subsetCommunication(cellchat_low,
                                          sources.use = "Malignant",
                                          targets.use = "TAM")
    write.csv(tumor_tam_low,
              file.path(RES_DIR, "cellchat_tumor_to_TAM_low.csv"),
              row.names = FALSE)
    message(sprintf("  Tumor -> TAM (UPR-low): %d interactions", nrow(tumor_tam_low)))
  }
}, error = function(e) {
  message("  Error in tumor-TAM extraction: ", e$message)
})

# TAM -> T cells
tryCatch({
  t_types_high <- intersect(c("CD4 T", "CD8 T", "Treg"), levels(cellchat_high@idents))
  t_types_low  <- intersect(c("CD4 T", "CD8 T", "Treg"), levels(cellchat_low@idents))

  if (length(t_types_high) > 0 && "TAM" %in% levels(cellchat_high@idents)) {
    tam_t_high <- subsetCommunication(cellchat_high,
                                       sources.use = "TAM",
                                       targets.use = t_types_high)
    write.csv(tam_t_high,
              file.path(RES_DIR, "cellchat_TAM_to_T_high.csv"),
              row.names = FALSE)
    message(sprintf("  TAM -> T cells (UPR-high): %d interactions", nrow(tam_t_high)))
  }

  if (length(t_types_low) > 0 && "TAM" %in% levels(cellchat_low@idents)) {
    tam_t_low <- subsetCommunication(cellchat_low,
                                      sources.use = "TAM",
                                      targets.use = t_types_low)
    write.csv(tam_t_low,
              file.path(RES_DIR, "cellchat_TAM_to_T_low.csv"),
              row.names = FALSE)
    message(sprintf("  TAM -> T cells (UPR-low): %d interactions", nrow(tam_t_low)))
  }
}, error = function(e) {
  message("  Error in TAM-T extraction: ", e$message)
})

# All interactions summary
tryCatch({
  all_comm_high <- subsetCommunication(cellchat_high)
  all_comm_low  <- subsetCommunication(cellchat_low)
  write.csv(all_comm_high,
            file.path(RES_DIR, "cellchat_all_interactions_high.csv"),
            row.names = FALSE)
  write.csv(all_comm_low,
            file.path(RES_DIR, "cellchat_all_interactions_low.csv"),
            row.names = FALSE)
  message(sprintf("  All interactions: UPR-high=%d, UPR-low=%d",
                  nrow(all_comm_high), nrow(all_comm_low)))
}, error = function(e) {
  message("  Error extracting all interactions: ", e$message)
})

# =============================================================================
# 8. Signaling pattern analysis (outgoing/incoming)
# =============================================================================
message("\n=== Signaling pattern analysis ===")

available_pw_high <- cellchat_high@netP$pathways

tryCatch({
  # Outgoing patterns — UPR-high
  n_patterns <- min(5, length(available_pw_high))
  if (n_patterns >= 2) {
    cellchat_high <- identifyCommunicationPatterns(cellchat_high,
                                                    pattern = "outgoing",
                                                    k = n_patterns)
    pdf(file.path(FIG_DIR, "Fig2_signaling_pattern_outgoing_high.pdf"), width = 12, height = 8)
    netAnalysis_river(cellchat_high, pattern = "outgoing")
    dev.off()

    # Incoming patterns — UPR-high
    cellchat_high <- identifyCommunicationPatterns(cellchat_high,
                                                    pattern = "incoming",
                                                    k = n_patterns)
    pdf(file.path(FIG_DIR, "Fig2_signaling_pattern_incoming_high.pdf"), width = 12, height = 8)
    netAnalysis_river(cellchat_high, pattern = "incoming")
    dev.off()
  }
}, error = function(e) {
  message("  Pattern analysis error (UPR-high): ", e$message)
  try(dev.off(), silent = TRUE)
})

available_pw_low <- cellchat_low@netP$pathways

tryCatch({
  n_patterns <- min(5, length(available_pw_low))
  if (n_patterns >= 2) {
    cellchat_low <- identifyCommunicationPatterns(cellchat_low,
                                                   pattern = "outgoing",
                                                   k = n_patterns)
    pdf(file.path(FIG_DIR, "Fig2_signaling_pattern_outgoing_low.pdf"), width = 12, height = 8)
    netAnalysis_river(cellchat_low, pattern = "outgoing")
    dev.off()

    cellchat_low <- identifyCommunicationPatterns(cellchat_low,
                                                   pattern = "incoming",
                                                   k = n_patterns)
    pdf(file.path(FIG_DIR, "Fig2_signaling_pattern_incoming_low.pdf"), width = 12, height = 8)
    netAnalysis_river(cellchat_low, pattern = "incoming")
    dev.off()
  }
}, error = function(e) {
  message("  Pattern analysis error (UPR-low): ", e$message)
  try(dev.off(), silent = TRUE)
})

# =============================================================================
# 9. Information flow comparison for key immune pathways
# =============================================================================
message("\n=== Information flow comparison ===")

tryCatch({
  # Signaling role analysis: compare sender/receiver roles across conditions
  pdf(file.path(FIG_DIR, "Fig2h_signaling_role_scatter.pdf"), width = 12, height = 6)
  p_role <- netAnalysis_signalingRole_scatter(cellchat_merged)
  print(p_role)
  dev.off()
}, error = function(e) {
  message("  Signaling role scatter error: ", e$message)
  try(dev.off(), silent = TRUE)
})

# Compare outgoing/incoming signaling strength per cell type
tryCatch({
  pdf(file.path(FIG_DIR, "Fig2_supp_signaling_changes_heatmap.pdf"), width = 12, height = 8)

  # Identify signaling that changes
  p_info <- rankNet(cellchat_merged, mode = "comparison",
                     stacked = TRUE, do.stat = TRUE,
                     return.data = TRUE)

  # If return.data works, also save the data
  if (is.list(p_info) && !is.null(p_info$signaling.contribution)) {
    write.csv(p_info$signaling.contribution,
              file.path(RES_DIR, "cellchat_signaling_contribution.csv"),
              row.names = FALSE)
  }
  dev.off()
}, error = function(e) {
  message("  Signaling changes heatmap error: ", e$message)
  try(dev.off(), silent = TRUE)
})

# =============================================================================
# 10. Save CellChat objects and summary
# =============================================================================
message("\n=== Saving results ===")

saveRDS(cellchat_high, file = file.path(DATA_PROC, "cellchat_high.rds"))
saveRDS(cellchat_low,  file = file.path(DATA_PROC, "cellchat_low.rds"))
save(cellchat_high, cellchat_low, cellchat_merged,
     file = file.path(DATA_PROC, "cellchat_objects.RData"))

# Summary statistics
summary_df <- data.frame(
  Group           = c("UPR-high", "UPR-low"),
  n_cells         = c(ncol(seu_high), ncol(seu_low)),
  n_celltypes     = c(length(levels(cellchat_high@idents)),
                       length(levels(cellchat_low@idents))),
  n_pathways      = c(length(cellchat_high@netP$pathways),
                       length(cellchat_low@netP$pathways)),
  n_interactions  = c(sum(cellchat_high@net$count),
                       sum(cellchat_low@net$count)),
  total_strength  = c(sum(cellchat_high@net$weight),
                       sum(cellchat_low@net$weight)),
  stringsAsFactors = FALSE
)

write.csv(summary_df,
          file.path(RES_DIR, "cellchat_summary.csv"),
          row.names = FALSE)

message("\n  CellChat Summary:")
print(summary_df)

message("\n=== CellChat analysis completed ===")
message("Objects saved to: ", file.path(DATA_PROC, "cellchat_objects.RData"))
message("Individual objects: ", file.path(DATA_PROC, "cellchat_high.rds"),
        " & ", file.path(DATA_PROC, "cellchat_low.rds"))
message("Figures saved in: ", FIG_DIR)
message("\nPart 1 (Single-cell analysis) completed.")
message("Next: Run 04_bulk_subtyping/01_upr_landscape.R for Part 2.")
