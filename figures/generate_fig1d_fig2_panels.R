###############################################################################
# generate_fig1d_fig2_panels.R
# Regenerate Fig1D and Fig2 chord/bubble panels per reviewer instructions:
#
#   Fig1D  : Dot plot of top 20 UPR_broad genes (by mean expression across all
#             cells) × cell type. No arm split.
#
#   Fig2A_chord_network_high.pdf  : Chord diagram of overall UPR-high network
#   Fig2C_chord_diff_count.pdf    : Chord of differential count (high − low)
#   Fig2D_chord_diff_weight.pdf   : Chord of differential weight (high − low)
#   Fig2f_bubble_key_pathways.pdf : Bubble plot Malignant+TAM (+Endothelial)
#                                   as key sender/receiver pairs in UPR-high
###############################################################################


suppressPackageStartupMessages({
  library(Seurat)
  library(CellChat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(circlize)
})

source("00_setup/config.R")
load(file.path(DATA_PROC, "upr_gene_sets.RData"))

FIG_DIR_MAIN <- file.path(FIG_DIR, "Main_Figures",
                           "Figure1_scRNA_UPR_landscape")
FIG_DIR_CC   <- file.path(FIG_DIR, "Main_Figures",
                           "Figure2_CellChat_communication")

###############################################################################
# ── PART 1 : Fig 1D ──────────────────────────────────────────────────────────
###############################################################################
message("\n=== Fig1D: top-20 UPR_broad dot plot ===")

seu <- readRDS(file.path(DATA_PROC, "seu_upr_scored.rds"))

# Genes present in the data
genes_present <- UPR_broad_genes[UPR_broad_genes %in% rownames(seu)]

# Rank by mean expression across ALL cells
expr_mat <- GetAssayData(seu, layer = "data")
gene_means <- rowMeans(expr_mat[genes_present, ])
top20_genes <- names(sort(gene_means, decreasing = TRUE))[1:20]

message(sprintf("  Top 20 UPR_broad genes: %s", paste(top20_genes, collapse = ", ")))

# Build a tidy data frame: for each gene × cell type, compute:
#   avg_exp  = mean log-normalised expression among all cells of that type
#   pct_exp  = % of cells with expression > 0
ct_levels <- sort(unique(seu$celltype))

dot_rows <- lapply(ct_levels, function(ct) {
  cells_ct <- colnames(seu)[seu$celltype == ct]
  mat_ct   <- expr_mat[top20_genes, cells_ct, drop = FALSE]
  data.frame(
    gene     = top20_genes,
    celltype = ct,
    avg_exp  = rowMeans(mat_ct),
    pct_exp  = rowMeans(mat_ct > 0) * 100,
    stringsAsFactors = FALSE
  )
})
dot_df <- do.call(rbind, dot_rows)

# Order genes by overall mean (largest first = bottom of y-axis in coord_flip)
dot_df$gene <- factor(dot_df$gene, levels = rev(top20_genes))

# Scale avg_exp within each gene for colour contrast
dot_df <- dot_df %>%
  group_by(gene) %>%
  mutate(avg_exp_scaled = scale(avg_exp)[, 1]) %>%
  ungroup()

p_fig1d <- ggplot(dot_df,
                   aes(x = celltype, y = gene,
                       size = pct_exp, colour = avg_exp_scaled)) +
  geom_point() +
  scale_colour_gradient2(
    low      = "#4DBBD5",
    mid      = "white",
    high     = "#E64B35",
    midpoint = 0,
    name     = "Scaled\nExpression"
  ) +
  scale_size_continuous(range = c(0.5, 7), name = "% Expressed") +
  THEME_PUBLICATION +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y  = element_text(size = 9),
    legend.position = "right",
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.3)
  ) +
  labs(
    x     = "Cell Type",
    y     = "UPR Gene",
    title = "Top 20 UPR Genes by Cell Type"
  )

out_1d_root <- file.path(FIG_DIR, "Fig1d_upr_dotplot.pdf")
out_1d_main <- file.path(FIG_DIR_MAIN, "Fig1d_upr_dotplot.pdf")

ggsave(out_1d_root, p_fig1d, width = 10, height = 8)
if (dir.exists(FIG_DIR_MAIN)) file.copy(out_1d_root, out_1d_main, overwrite = TRUE)

message("  Fig1D saved: ", out_1d_root)

# Release Seurat object from memory
rm(seu, expr_mat)
gc()

###############################################################################
# ── PART 2 : Fig 2 Chord / Bubble panels ─────────────────────────────────────
###############################################################################
message("\n=== Loading CellChat objects ===")

cc_high <- readRDS(file.path(DATA_PROC, "cellchat_high.rds"))
cc_low  <- readRDS(file.path(DATA_PROC, "cellchat_low.rds"))

object_list  <- list("UPR-high" = cc_high, "UPR-low" = cc_low)
cc_merged    <- mergeCellChat(object_list, add.names = names(object_list))

# ── Fig2A: chord diagram of UPR-high overall network ─────────────────────────
message("\n  Fig2A: chord network UPR-high")

out_2A_root <- file.path(FIG_DIR, "Fig2A_chord_network_high.pdf")
out_2A_main <- file.path(FIG_DIR_CC, "Fig2A_chord_network_high.pdf")

tryCatch({
  pdf(out_2A_root, width = 10, height = 9)
  netVisual_aggregate(cc_high,
                      signaling    = cc_high@netP$pathways,
                      layout       = "chord",
                      title.name   = "UPR-high: All Pathways (Chord)")
  dev.off()
  if (dir.exists(FIG_DIR_CC)) file.copy(out_2A_root, out_2A_main, overwrite = TRUE)
  message("  Fig2A_chord_network_high.pdf saved.")
}, error = function(e) {
  try(dev.off(), silent = TRUE)
  message("  Fig2A chord failed: ", e$message)
  message("  Falling back to circle network diagram...")
  tryCatch({
    pdf(out_2A_root, width = 10, height = 9)
    netVisual_circle(cc_high@net$count,
                     vertex.weight = table(cc_high@idents),
                     weight.scale  = TRUE,
                     label.edge    = FALSE,
                     title.name    = "UPR-high: Number of Interactions")
    dev.off()
    if (dir.exists(FIG_DIR_CC)) file.copy(out_2A_root, out_2A_main, overwrite = TRUE)
    message("  Fig2A fallback (circle) saved.")
  }, error = function(e2) {
    try(dev.off(), silent = TRUE)
    message("  Fig2A circle also failed: ", e2$message)
  })
})

# ── Fig2C: chord of differential count ───────────────────────────────────────
message("\n  Fig2C: chord differential count")

out_2C_root <- file.path(FIG_DIR, "Fig2C_chord_diff_count.pdf")
out_2C_main <- file.path(FIG_DIR_CC, "Fig2C_chord_diff_count.pdf")

tryCatch({
  pdf(out_2C_root, width = 10, height = 9)
  netVisual_diffInteraction(cc_merged,
                             weight.scale = TRUE,
                             measure      = "count",
                             title.name   = "Differential Interactions\n(UPR-high − UPR-low)")
  dev.off()
  if (dir.exists(FIG_DIR_CC)) file.copy(out_2C_root, out_2C_main, overwrite = TRUE)
  message("  Fig2C_chord_diff_count.pdf saved.")
}, error = function(e) {
  try(dev.off(), silent = TRUE)
  message("  Fig2C chord diff count failed: ", e$message)
})

# ── Fig2D: chord of differential weight ──────────────────────────────────────
message("\n  Fig2D: chord differential weight")

out_2D_root <- file.path(FIG_DIR, "Fig2D_chord_diff_weight.pdf")
out_2D_main <- file.path(FIG_DIR_CC, "Fig2D_chord_diff_weight.pdf")

tryCatch({
  pdf(out_2D_root, width = 10, height = 9)
  netVisual_diffInteraction(cc_merged,
                             weight.scale = TRUE,
                             measure      = "weight",
                             title.name   = "Differential Interaction Strength\n(UPR-high − UPR-low)")
  dev.off()
  if (dir.exists(FIG_DIR_CC)) file.copy(out_2D_root, out_2D_main, overwrite = TRUE)
  message("  Fig2D_chord_diff_weight.pdf saved.")
}, error = function(e) {
  try(dev.off(), silent = TRUE)
  message("  Fig2D chord diff weight failed: ", e$message)
})

# ── Fig2F: bubble plot Malignant + TAM (+ Endothelial) ───────────────────────
message("\n  Fig2F: bubble plot Malignant × TAM (+ Endothelial) — UPR-high")

# Cell types present in UPR-high
ct_high    <- levels(cc_high@idents)
senders    <- intersect(c("Malignant", "TAM"), ct_high)
receivers  <- intersect(c("TAM", "Malignant"), ct_high)

# Add Endothelial as third cell type if space permits
if ("Endothelial" %in% ct_high) {
  senders   <- union(senders,   "Endothelial")
  receivers <- union(receivers, "Endothelial")
}

message(sprintf("  Senders: %s", paste(senders, collapse = ", ")))
message(sprintf("  Receivers: %s", paste(receivers, collapse = ", ")))

# Select meaningful pathways: shared between high and low with at least one
# Malignant/TAM/Endothelial interaction in UPR-high
key_pw_candidates <- c(
  "SPP1", "MHC-I", "MHC-II", "CSF", "PD-L1", "GALECTIN",
  "CXCL", "CCL", "TNF", "VEGF", "PDGF", "FGF",
  "COMPLEMENT", "CD86", "ICAM", "ANNEXIN"
)
shared_pathways <- intersect(cc_high@netP$pathways, cc_low@netP$pathways)
key_pw_use <- intersect(key_pw_candidates, shared_pathways)

if (length(key_pw_use) == 0) {
  # fall back to all shared pathways
  key_pw_use <- shared_pathways
}
message(sprintf("  Key pathways for bubble: %s", paste(key_pw_use, collapse = ", ")))

out_2F_root     <- file.path(FIG_DIR, "Fig2f_bubble_key_pathways.pdf")
out_2F_low_root <- file.path(FIG_DIR, "Fig2f_bubble_key_pathways_low.pdf")
out_2F_main     <- file.path(FIG_DIR_CC, "Fig2f_bubble_key_pathways.pdf")
out_2F_low_main <- file.path(FIG_DIR_CC, "Fig2f_bubble_key_pathways_low.pdf")

tryCatch({
  p_bubble_high <- netVisual_bubble(cc_high,
                                     sources.use    = senders,
                                     targets.use    = receivers,
                                     signaling      = key_pw_use,
                                     remove.isolate = TRUE)
  pdf(out_2F_root, width = 14, height = 10)
  print(p_bubble_high)
  dev.off()
  if (dir.exists(FIG_DIR_CC)) file.copy(out_2F_root, out_2F_main, overwrite = TRUE)
  message("  Fig2f_bubble_key_pathways.pdf (UPR-high) saved.")
}, error = function(e) {
  try(dev.off(), silent = TRUE)
  message("  Fig2F UPR-high bubble failed: ", e$message)
})

ct_low <- levels(cc_low@idents)
senders_low   <- intersect(c("Malignant", "TAM", "Endothelial"), ct_low)
receivers_low <- intersect(c("Malignant", "TAM", "Endothelial"), ct_low)
key_pw_low <- intersect(key_pw_use, cc_low@netP$pathways)

tryCatch({
  p_bubble_low <- netVisual_bubble(cc_low,
                                    sources.use    = senders_low,
                                    targets.use    = receivers_low,
                                    signaling      = key_pw_low,
                                    remove.isolate = TRUE)
  pdf(out_2F_low_root, width = 14, height = 10)
  print(p_bubble_low)
  dev.off()
  if (dir.exists(FIG_DIR_CC)) file.copy(out_2F_low_root, out_2F_low_main, overwrite = TRUE)
  message("  Fig2f_bubble_key_pathways_low.pdf (UPR-low) saved.")
}, error = function(e) {
  try(dev.off(), silent = TRUE)
  message("  Fig2F UPR-low bubble failed: ", e$message)
})

message("\n=== All panels generated. ===")
message("Fig1D : ", out_1d_root)
message("Fig2A chord: ", out_2A_root)
message("Fig2C chord: ", out_2C_root)
message("Fig2D chord: ", out_2D_root)
message("Fig2F bubble (high): ", out_2F_root)
message("Fig2F bubble (low):  ", out_2F_low_root)
