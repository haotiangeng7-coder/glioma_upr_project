###############################################################################
# fig3_composite_v3.R
# Figure 3 composite — immune & genomic landscape  (4-panel version)
# Panels:
#   A  ssGSEA immune-cell infiltration heatmap (27 populations, full width ~85 mm)
#   B  Immune-cell infiltration by subtype — top 8, violin/box, 2×4 facet (~70 mm)
#   C  Immune-checkpoint gene expression (6 genes, violin/box, 2×3 facet, ~55 mm)
#   D  GSVA Hallmark pathway heatmap — top 20 pathways (~65 mm)
#   GSEA GO:BP dotplot removed (moves to Supplementary)
# Target canvas: 183 × 270 mm  (aspect ≈ 0.68)
###############################################################################


suppressPackageStartupMessages({
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(tidyr)
  library(ggpubr)
  library(patchwork)
  library(cowplot)
  library(png)
  library(grid)
  library(grDevices)
  library(showtext)
  library(sysfonts)
  library(stringr)
})

font_add("Times",
         regular    = "/usr/share/fonts/truetype/dejavu/DejaVuSerif.ttf",
         bold       = "/usr/share/fonts/truetype/dejavu/DejaVuSerif-Bold.ttf",
         italic     = "/usr/share/fonts/truetype/dejavu/DejaVuSerif-Italic.ttf",
         bolditalic = "/usr/share/fonts/truetype/dejavu/DejaVuSerif-BoldItalic.ttf")
showtext_auto()
showtext_opts(dpi = 300)

PROJ <- getwd()
DATA_PROC <- file.path(PROJ, "data", "processed")
OUT_DIR   <- file.path(PROJ, "manuscript", "submission", "Figures_v2")
TMP_DIR   <- file.path(OUT_DIR, "tmp_fig3v3")
dir.create(TMP_DIR, recursive=TRUE, showWarnings=FALSE)

W_MM      <- 183          # canvas width (mm)
MM_IN     <- 1/25.4       # mm -> inch conversion
FONT      <- "Times"
BASE      <- 8

# Target panel heights (mm) — must sum to ~270 minus patchwork overhead
# patchwork adds ~3 mm per gap × 3 gaps = ~9 mm overhead → data rows = ~261 mm
H_A  <- 104  # ssGSEA immune heatmap (extra bottom space for the bottom legend)
H_B  <- 80   # immune-cell violin/box 2×4 (extra top margin for strip label clearance)
H_C  <- 55   # checkpoint violin/box 2×3
H_D  <- 86   # GSVA heatmap (extra bottom space for the bottom legend)
# Total data rows: 104+80+55+86 = 325 mm

TOTAL_H_MM <- 316

COL_SUBTYPE <- c(
  "UPR-high-risk"    = "#E64B35",
  "UPR-intermediate" = "#00A087",
  "UPR-favorable"    = "#4DBBD5"
)
SUBTYPE_ORDER <- c("UPR-high-risk", "UPR-intermediate", "UPR-favorable")

theme_pub <- function(bs = BASE) {
  theme_classic(base_size = bs, base_family = FONT) +
    theme(
      axis.line        = element_line(linewidth = 0.4),
      axis.ticks       = element_line(linewidth = 0.3),
      axis.text        = element_text(size = bs - 1, color = "black"),
      axis.title       = element_text(size = bs),
      strip.text       = element_text(size = bs - 0.5, face = "bold"),
      strip.background = element_blank(),
      legend.text      = element_text(size = bs - 1),
      legend.title     = element_text(size = bs),
      legend.key.size  = unit(3, "mm"),
      plot.title       = element_text(size = bs, face = "bold", hjust = 0),
      panel.grid       = element_blank()
    )
}

sig_star <- function(p) {
  dplyr::case_when(
    p < 0.0001 ~ "****",
    p < 0.001  ~ "***",
    p < 0.01   ~ "**",
    p < 0.05   ~ "*",
    TRUE       ~ "ns"
  )
}

###############################################################################
# Load data
###############################################################################
message("Loading data ...")
load(file.path(DATA_PROC, "consensus_clustering_results.RData"))
load(file.path(DATA_PROC, "immune_analysis_results.RData"))
load(file.path(DATA_PROC, "pathway_analysis_results.RData"))
load(file.path(DATA_PROC, "tcga_glioma_expression.RData"))

subtype_map   <- setNames(cluster_df$UPR_subtype, cluster_df$barcode)
SUBTYPE_ORDER <- intersect(SUBTYPE_ORDER, unique(as.character(cluster_df$UPR_subtype)))
common_samples <- intersect(colnames(immune_ssgsea), cluster_df$barcode)
subtypes       <- factor(subtype_map[common_samples], levels = SUBTYPE_ORDER)

###############################################################################
# Panel A — ssGSEA immune heatmap (save as intermediate PNG, embed in composite)
###############################################################################
message("Panel A (immune heatmap) ...")
sample_order     <- order(subtypes)
ordered_samps    <- common_samples[sample_order]
ordered_subtypes <- subtypes[sample_order]

ssgsea_mat <- immune_ssgsea[, ordered_samps]
ssgsea_z   <- t(scale(t(ssgsea_mat)))
ssgsea_z[ssgsea_z >  3] <-  3
ssgsea_z[ssgsea_z < -3] <- -3

sig_cells   <- ssgsea_kw$CellType[!is.na(ssgsea_kw$kw_padj) & ssgsea_kw$kw_padj < 0.05]
row_sig_vec <- ifelse(rownames(ssgsea_z) %in% sig_cells, "FDR < 0.05", "ns")

broad_cat <- dplyr::case_when(
  grepl("CD4|CD8|T helper|Regulatory|T follicular|gamma delta|Natural killer T|memory CD",
        rownames(ssgsea_z)) ~ "T/NK cell",
  grepl("Natural killer cell", rownames(ssgsea_z)) ~ "T/NK cell",
  grepl("B cell|B lineage|Immature B",  rownames(ssgsea_z)) ~ "B cell",
  grepl("Macrophage|Monocyte|MDSC|dendritic|Neutrophil|Eosinophil|Mast",
        rownames(ssgsea_z)) ~ "Myeloid",
  grepl("IFN", rownames(ssgsea_z)) ~ "IFN response",
  TRUE ~ "Other"
)
cat_colors <- c(
  "T/NK cell"    = "#3C5488",
  "B cell"       = "#00A087",
  "Myeloid"      = "#F39B7F",
  "IFN response" = "#8491B4",
  "Other"        = "grey80"
)

row_anno <- rowAnnotation(
  Significant = row_sig_vec,
  Category    = broad_cat,
  col = list(
    Significant = c("FDR < 0.05" = "#E64B35", "ns" = "grey85"),
    Category    = cat_colors
  ),
  annotation_name_gp = gpar(fontsize = 7, fontfamily = FONT),
  simple_anno_size   = unit(3, "mm"),
  annotation_legend_param = list(
    Significant = list(title_gp  = gpar(fontsize = 7, fontfamily = FONT),
                       labels_gp = gpar(fontsize = 6, fontfamily = FONT)),
    Category    = list(title_gp  = gpar(fontsize = 7, fontfamily = FONT),
                       labels_gp = gpar(fontsize = 6, fontfamily = FONT),
                       nrow      = 2)
  )
)

st_n   <- table(ordered_subtypes)[SUBTYPE_ORDER]
ct_ttl <- paste0(SUBTYPE_ORDER, "\n(n=", st_n, ")")

col_anno_A <- HeatmapAnnotation(
  `UPR subtype` = ordered_subtypes,
  col = list(`UPR subtype` = COL_SUBTYPE),
  annotation_name_gp   = gpar(fontsize = 7, fontfamily = FONT),
  annotation_name_side = "left",
  simple_anno_size     = unit(3, "mm"),
  annotation_legend_param = list(
    `UPR subtype` = list(title_gp  = gpar(fontsize = 7, fontfamily = FONT),
                          labels_gp = gpar(fontsize = 6, fontfamily = FONT))
  )
)

ht_immune <- Heatmap(
  ssgsea_z,
  name      = "ssGSEA\n(z-score)",
  col       = colorRamp2(c(-3, 0, 3), c("#2166AC", "white", "#B2182B")),
  top_annotation   = col_anno_A,
  # Row annotation moved to the LEFT so the cell-type row names sit outermost on
  # the right edge and the right padding gives them full room to render without
  # being clipped by the annotation bars
  left_annotation  = row_anno,
  # Pin body width so the remaining canvas reserves room for the longest
  # row names (e.g. "Plasmacytoid dendritic cell"); prevents right-edge clipping
  width                     = unit(118, "mm"),
  cluster_columns           = FALSE,
  cluster_rows              = TRUE,
  clustering_method_rows    = "ward.D2",
  show_column_names         = FALSE,
  row_names_gp              = gpar(fontsize = 7, fontfamily = FONT),
  column_split              = factor(ordered_subtypes, levels = SUBTYPE_ORDER),
  column_title              = ct_ttl,
  column_title_gp           = gpar(fontsize = 8, fontface = "bold", fontfamily = FONT),
  heatmap_legend_param      = list(
    title_gp      = gpar(fontsize = 7, fontfamily = FONT),
    labels_gp     = gpar(fontsize = 6, fontfamily = FONT),
    direction     = "horizontal",
    legend_width  = unit(30, "mm")
  ),
  use_raster = FALSE
)

ht_A_png <- file.path(TMP_DIR, "ht_A.png")
# Render at H_A mm tall to match composite slot
# padding: extra right space so row labels don't overlap the merged legend
png(ht_A_png, width = W_MM, height = H_A, units = "mm", res = 300, family = FONT)
draw(ht_immune, merge_legend = TRUE,
     heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
     # Legends moved to the bottom so the heatmap body and long row labels use
     # the full canvas width; right gutter enlarged so the longest cell-type row
     # names (e.g. "Effector memory CD4 T cell") render fully without clipping at
     # the canvas right edge; bottom gutter hosts the merged legend.
     padding = unit(c(2, 6, 2, 2), "mm"))  # bottom, right, top, left
dev.off()
message("  Panel A saved (", W_MM, " x ", H_A, " mm)")

###############################################################################
# Panel B — Immune-cell infiltration by subtype (top 8, violin+box, 2 rows × 4)
###############################################################################
message("Panel B (immune violin/box top 8) ...")

top_cells <- ssgsea_kw %>%
  dplyr::filter(!is.na(kw_padj), kw_padj < 0.05) %>%
  dplyr::arrange(kw_padj) %>%
  dplyr::slice_head(n = 8) %>%
  dplyr::pull(CellType)

if (length(top_cells) == 0) top_cells <- sig_cells[seq_len(min(8, length(sig_cells)))]
message("  Top 8 cell types: ", paste(top_cells, collapse = ", "))

# Word-wrap strip labels at 16 chars so 4-column facet strips are fully legible
# Use \n (literal newline) via str_wrap — ggplot facet_wrap honours \n in factor levels
wrap_label <- function(x, w = 16) stringr::str_wrap(x, width = w)
top_cells_wrapped <- sapply(top_cells, wrap_label, USE.NAMES = FALSE)
message("  Wrapped labels: ", paste(top_cells_wrapped, collapse = " | "))

immune_long <- as.data.frame(t(immune_ssgsea[top_cells, common_samples])) %>%
  mutate(barcode = rownames(.), UPR_subtype = subtypes) %>%
  pivot_longer(-c(barcode, UPR_subtype), names_to = "CellType", values_to = "Score") %>%
  mutate(
    CellType_wrapped = factor(
      wrap_label(CellType),
      levels = top_cells_wrapped
    ),
    UPR_subtype = factor(UPR_subtype, levels = SUBTYPE_ORDER)
  )

comparisons_3 <- Filter(function(p) all(p %in% SUBTYPE_ORDER), list(
  c("UPR-high-risk", "UPR-intermediate"),
  c("UPR-high-risk", "UPR-favorable"),
  c("UPR-intermediate", "UPR-favorable")
))

p_B_base <- ggplot(immune_long, aes(x = UPR_subtype, y = Score, fill = UPR_subtype)) +
  geom_violin(alpha = 0.5, trim = TRUE, linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.35,
               fill = "white", color = "black") +
  stat_compare_means(
    comparisons       = comparisons_3,
    method            = "wilcox.test",
    p.adjust.method   = "BH",
    label             = "p.signif",
    hide.ns           = TRUE,
    size              = 1.8,
    tip.length        = 0.01,
    step.increase     = 0.06,
    family            = FONT
  ) +
  scale_fill_manual(values = COL_SUBTYPE, guide = "none") +
  scale_x_discrete(labels = c("High-risk", "Interm.", "Favorable")) +
  # Extra top expansion so the significance brackets and **** labels above the
  # tallest violins are not clipped by each facet's upper panel boundary
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  facet_wrap(~ CellType_wrapped, scales = "free_y", ncol = 4) +
  labs(x = "", y = "ssGSEA score") +
  theme_pub(bs = 7) +
  theme(
    axis.text.x  = element_text(angle = 35, hjust = 1, size = 5.5),
    axis.text.y  = element_text(size = 5.5),
    axis.title.y = element_text(size = 6.5),
    strip.text   = element_text(size = 7, face = "bold", lineheight = 1.1),
    strip.clip   = "off",
    # generous margins so strip labels at row 1 top and sig brackets at row 2 bottom
    # are never clipped by the PNG canvas edge
    plot.margin  = margin(4, 4, 4, 4, "mm")
  )

# Save panel B as intermediate PNG (same strategy as A and D) so patchwork tag
# placement never interferes with the facet strip labels
ht_B_png <- file.path(TMP_DIR, "ht_B.png")
ggsave(ht_B_png, plot = p_B_base,
       width = W_MM * MM_IN, height = H_B * MM_IN,
       units = "in", dpi = 300)
message("  Panel B saved (", W_MM, " x ", H_B, " mm)")

###############################################################################
# Panel C — Immune-checkpoint gene expression (6 genes, 2-row × 3-col facet)
###############################################################################
message("Panel C (checkpoint genes) ...")

checkpoint_show <- c("CD274", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "PDCD1")
checkpoint_show <- checkpoint_show[checkpoint_show %in% icp_kw$Gene]

expr_log <- log2(expr_tpm_symbol[, common_samples] + 1)
icp_expr <- expr_log[checkpoint_show[checkpoint_show %in% rownames(expr_log)], ]

gene_display <- c(
  "CD274"  = "PD-L1\n(CD274)",
  "HAVCR2" = "TIM3\n(HAVCR2)",
  "PDCD1"  = "PD-1\n(PDCD1)",
  "CTLA4"  = "CTLA4",
  "LAG3"   = "LAG3",
  "TIGIT"  = "TIGIT"
)

icp_long_show <- as.data.frame(t(icp_expr)) %>%
  mutate(barcode = rownames(.), UPR_subtype = subtypes) %>%
  pivot_longer(-c(barcode, UPR_subtype), names_to = "Gene", values_to = "Expression") %>%
  mutate(
    Gene        = factor(Gene, levels = checkpoint_show),
    UPR_subtype = factor(UPR_subtype, levels = SUBTYPE_ORDER),
    Gene_label  = dplyr::recode(as.character(Gene), !!!gene_display)
  )

icp_padj <- icp_kw[icp_kw$Gene %in% checkpoint_show, c("Gene", "kw_padj")]
icp_long_show <- left_join(icp_long_show, icp_padj, by = "Gene") %>%
  mutate(panel_label = paste0(Gene_label, "\n", sig_star(kw_padj)),
         panel_label = factor(panel_label, levels = unique(panel_label)))

p_C_base <- ggplot(icp_long_show, aes(x = UPR_subtype, y = Expression, fill = UPR_subtype)) +
  geom_violin(alpha = 0.5, trim = TRUE, linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.35,
               fill = "white", color = "black") +
  stat_compare_means(
    comparisons     = comparisons_3,
    method          = "wilcox.test",
    p.adjust.method = "BH",
    label           = "p.signif",
    hide.ns         = TRUE,
    size            = 1.8,
    tip.length      = 0.01,
    step.increase   = 0.08,
    family          = FONT
  ) +
  scale_fill_manual(values = COL_SUBTYPE, guide = "none") +
  scale_x_discrete(labels = c("High-risk", "Interm.", "Favorable")) +
  facet_wrap(~ panel_label, scales = "free_y", nrow = 2) +
  labs(x = "", y = "Expr. (log2 TPM+1)") +
  theme_pub(bs = 7) +
  theme(
    axis.text.x  = element_text(angle = 35, hjust = 1, size = 5.5),
    axis.text.y  = element_text(size = 5.5),
    axis.title.y = element_text(size = 6.5),
    strip.text   = element_text(size = 6.5, face = "bold"),
    plot.margin  = margin(4, 4, 4, 4, "mm")
  )

ht_C_png <- file.path(TMP_DIR, "ht_C.png")
ggsave(ht_C_png, plot = p_C_base,
       width = W_MM * MM_IN, height = H_C * MM_IN,
       units = "in", dpi = 300)
message("  Panel C saved (", W_MM, " x ", H_C, " mm)")

###############################################################################
# Panel D — GSVA Hallmark pathway heatmap: top 20 most significant pathways
###############################################################################
message("Panel D (GSVA top 20) ...")

gsva_samples  <- intersect(colnames(gsva_hallmark), common_samples)
gsva_subtypes <- factor(subtype_map[gsva_samples], levels = SUBTYPE_ORDER)
gsva_order    <- order(gsva_subtypes)
gsva_ordered  <- gsva_hallmark[, gsva_samples[gsva_order]]
gsva_sub_ord  <- gsva_subtypes[gsva_order]

gsva_z <- t(scale(t(gsva_ordered)))
gsva_z[gsva_z >  2] <-  2
gsva_z[gsva_z < -2] <- -2

# Clean pathway names: strip HALLMARK_ prefix, replace _ with space, title-case
clean_path <- function(x) {
  x <- gsub("HALLMARK_", "", x)
  x <- gsub("_", " ", x)
  x <- tolower(x)
  paste0(toupper(substr(x, 1, 1)), substr(x, 2, nchar(x)))
}
rownames(gsva_z) <- clean_path(rownames(gsva_z))

aq_res <- diff_pathway_results[[1]]
aq_res$pathway_short <- clean_path(aq_res$pathway)

# Top 20 by adjusted p-value
top20_paths <- aq_res %>%
  dplyr::arrange(adj.P.Val) %>%
  dplyr::slice_head(n = 20) %>%
  dplyr::pull(pathway_short)
top20_paths <- top20_paths[top20_paths %in% rownames(gsva_z)]
message(sprintf("  GSVA top %d pathways selected", length(top20_paths)))
gsva_z_top <- gsva_z[top20_paths, ]

sig_paths <- aq_res$pathway_short[aq_res$adj.P.Val < 0.05]
row_diff  <- ifelse(rownames(gsva_z_top) %in% sig_paths, "FDR < 0.05", "ns")
row_dir   <- dplyr::case_when(
  !(rownames(gsva_z_top) %in% sig_paths) ~ "ns",
  aq_res$logFC[match(rownames(gsva_z_top), aq_res$pathway_short)] > 0 ~ "Higher in Favorable",
  TRUE ~ "Higher in High-risk"
)
dir_colors <- c(
  "Higher in Favorable" = "#4DBBD5",
  "Higher in High-risk" = "#E64B35",
  "ns"                  = "grey85"
)

row_anno_gsva <- rowAnnotation(
  Direction = row_dir,
  col = list(Direction = dir_colors),
  annotation_name_gp = gpar(fontsize = 6, fontfamily = FONT),
  simple_anno_size   = unit(3, "mm"),
  annotation_legend_param = list(
    Direction = list(
      title_gp  = gpar(fontsize = 7, fontfamily = FONT),
      labels_gp = gpar(fontsize = 6, fontfamily = FONT)
    )
  )
)

st_n_gsva   <- table(gsva_sub_ord)[SUBTYPE_ORDER]
ct_ttl_gsva <- paste0(SUBTYPE_ORDER, "\n(n=", st_n_gsva, ")")

col_anno_gsva <- HeatmapAnnotation(
  `UPR subtype` = gsva_sub_ord,
  col = list(`UPR subtype` = COL_SUBTYPE),
  annotation_name_gp   = gpar(fontsize = 6, fontfamily = FONT),
  annotation_name_side = "left",
  simple_anno_size     = unit(3, "mm"),
  annotation_legend_param = list(
    `UPR subtype` = list(
      title_gp  = gpar(fontsize = 7, fontfamily = FONT),
      labels_gp = gpar(fontsize = 6, fontfamily = FONT)
    )
  )
)

ht_gsva <- Heatmap(
  gsva_z_top,
  name     = "GSVA\n(z-score)",
  col      = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
  top_annotation   = col_anno_gsva,
  right_annotation = row_anno_gsva,
  cluster_columns        = FALSE,
  cluster_rows           = TRUE,
  clustering_method_rows = "ward.D2",
  show_column_names      = FALSE,
  row_names_gp           = gpar(fontsize = 7, fontfamily = FONT),
  column_split           = factor(gsva_sub_ord, levels = SUBTYPE_ORDER),
  column_title           = ct_ttl_gsva,
  column_title_gp        = gpar(fontsize = 7, fontface = "bold", fontfamily = FONT),
  heatmap_legend_param   = list(
    title_gp      = gpar(fontsize = 7, fontfamily = FONT),
    labels_gp     = gpar(fontsize = 6, fontfamily = FONT),
    direction     = "horizontal",
    legend_width  = unit(30, "mm")
  ),
  use_raster = FALSE
)

ht_D_png <- file.path(TMP_DIR, "ht_D.png")
png(ht_D_png, width = W_MM, height = H_D, units = "mm", res = 300, family = FONT)
draw(ht_gsva, merge_legend = TRUE,
     heatmap_legend_side = "bottom", annotation_legend_side = "bottom",
     # Legends moved to the bottom so the heatmap body and pathway row labels
     # use the full canvas width; right gutter trimmed to a small margin; bottom
     # gutter enlarged to host the merged Direction/Category legend without clipping.
     padding = unit(c(2, 5, 2, 2), "mm"))  # bottom, right, top, left
dev.off()
message("  Panel D saved (", W_MM, " x ", H_D, " mm)")

###############################################################################
# Assemble composite: A / B / C / D  (4 panels, no GSEA dotplot)
# Strategy: all four panels rendered as fixed-size PNGs; assembled with
# cowplot::plot_grid so panel label placement never clips facet strip labels.
###############################################################################
message("=== Assembling 4-panel composite ===")

# Load all four intermediate PNGs as raster grobs
make_panel <- function(png_path, label) {
  img  <- png::readPNG(png_path)
  grob <- grid::rasterGrob(img, interpolate = TRUE)
  # Draw raster on a blank ggdraw canvas; add bold panel label top-left
  cowplot::ggdraw() +
    cowplot::draw_grob(grob) +
    cowplot::draw_label(label,
                        x = 0.01, y = 0.99,
                        hjust = 0, vjust = 1,
                        fontface = "bold",
                        size = 13,
                        fontfamily = FONT,
                        color = "black") +
    theme(plot.margin = margin(0, 0, 0, 0))
}

p_A_final <- make_panel(ht_A_png, "A")
p_B_final <- make_panel(ht_B_png, "B")
p_C_final <- make_panel(ht_C_png, "C")
p_D_final <- make_panel(ht_D_png, "D")

# Stack with cowplot::plot_grid; rel_heights proportional to panel mm heights
p_composite <- cowplot::plot_grid(
  p_A_final, p_B_final, p_C_final, p_D_final,
  ncol        = 1,
  rel_heights = c(H_A, H_B, H_C, H_D)
)

###############################################################################
# Save outputs
###############################################################################
out_pdf <- file.path(OUT_DIR, "Figure3_composite.pdf")
out_png <- file.path(OUT_DIR, "Figure3_composite.png")

ggsave(out_pdf, plot = p_composite,
       device = cairo_pdf,
       width  = W_MM * MM_IN,
       height = TOTAL_H_MM * MM_IN,
       units  = "in",
       dpi    = 300)
message("Saved PDF: ", out_pdf)

ggsave(out_png, plot = p_composite,
       width  = W_MM * MM_IN,
       height = TOTAL_H_MM * MM_IN,
       units  = "in",
       dpi    = 300)
message("Saved PNG: ", out_png)

# Verify file sizes
pdf_size <- file.info(out_pdf)$size
png_size <- file.info(out_png)$size
message(sprintf("PDF size: %.1f KB", pdf_size / 1024))
message(sprintf("PNG size: %.1f KB", png_size / 1024))

# Clean up temp dir
unlink(TMP_DIR, recursive = TRUE)

message("=== Figure 3 (4-panel) done ===")
message(sprintf("Canvas: %d x %d mm  (W/H aspect = %.2f)",
                W_MM, TOTAL_H_MM, W_MM / TOTAL_H_MM))
message("Panels present: A (ssGSEA heatmap), B (immune violin/box), ",
        "C (checkpoint), D (GSVA heatmap)")
message("Panel E (GSEA dotplot) removed — moved to Supplementary.")
message("Panel B strip labels word-wrapped at 16 chars, ncol=4, nrow=2.")
message(sprintf("Row heights: A=%dmm, B=%dmm, C=%dmm, D=%dmm", H_A, H_B, H_C, H_D))
