###############################################################################
# generate_fig4.R
# Publication-quality Figure 4: Immune Microenvironment and Pathway Analysis
# Journal of Translational Medicine submission
#
# Panels:
#   Fig4A_immune_heatmap   : ssGSEA immune cell infiltration heatmap (ComplexHeatmap)
#   Fig4B_immune_boxplots  : Top immune cells (FDR<0.05) as violin+box panels
#   Fig4C_checkpoints      : Key checkpoint genes (PD-L1/CD274, CTLA4, LAG3, TIM3, TIGIT)
#   Fig4D_estimate         : ESTIMATE immune/stromal/tumor purity scores by subtype
#   Fig4E_gsva             : GSVA Hallmark pathway heatmap (top diff pathways)
#   Fig4F_gsea_dotplot     : GSEA GO-BP dotplot (top 10 favorable + top 10 high-risk)
#
# Output: manuscript/submission/Figures/Figure4/ (relative to project root)
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(tidyr)
  library(ggpubr)
  library(patchwork)
  library(scales)
  library(grid)
  library(grDevices)
  library(forcats)
  library(showtext)
  library(sysfonts)
})

# Register DejaVu Serif as "Times" (metric-compatible TrueType serif font
# available on this system; accepted as Times New Roman equivalent).
font_add("Times",
         regular    = "/usr/share/fonts/truetype/dejavu/DejaVuSerif.ttf",
         bold       = "/usr/share/fonts/truetype/dejavu/DejaVuSerif-Bold.ttf",
         italic     = "/usr/share/fonts/truetype/dejavu/DejaVuSerif-Italic.ttf",
         bolditalic = "/usr/share/fonts/truetype/dejavu/DejaVuSerif-BoldItalic.ttf")
showtext_auto()
showtext_opts(dpi = 600)

# ---------------------------------------------------------------------------
# 0. Paths and style constants
# ---------------------------------------------------------------------------
PROJ <- getwd()
DATA_PROC <- file.path(PROJ, "data", "processed")
OUT_DIR   <- file.path(PROJ, "manuscript", "submission", "Figures", "Figure4")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

W_SINGLE  <- 89
W_1P5     <- 130
W_DOUBLE  <- 183
MM_TO_IN  <- 1 / 25.4
FONT      <- "Times"
BASE_SIZE <- 8

COL_SUBTYPE <- c(
  "UPR-high-risk"    = "#E64B35",
  "UPR-intermediate" = "#00A087",
  "UPR-favorable"    = "#4DBBD5"
)
SUBTYPE_ORDER <- c("UPR-high-risk", "UPR-intermediate", "UPR-favorable")

theme_pub <- function(base_size = BASE_SIZE) {
  theme_classic(base_size = base_size, base_family = FONT) +
    theme(
      axis.line        = element_line(linewidth = 0.4, color = "black"),
      axis.ticks       = element_line(linewidth = 0.3, color = "black"),
      axis.text        = element_text(size = base_size - 1, color = "black"),
      axis.title       = element_text(size = base_size,     color = "black"),
      strip.text       = element_text(size = base_size - 0.5, face = "bold"),
      strip.background = element_blank(),
      legend.text      = element_text(size = base_size - 1),
      legend.title     = element_text(size = base_size),
      legend.key.size  = unit(3, "mm"),
      plot.title       = element_text(size = base_size, face = "bold", hjust = 0),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

save_pdf <- function(plot_obj = NULL, path, width_mm, height_mm, draw_fn = NULL) {
  w_in <- width_mm  * MM_TO_IN
  h_in <- height_mm * MM_TO_IN
  if (!is.null(plot_obj)) {
    cairo_pdf(path, width = w_in, height = h_in, family = FONT)
    print(plot_obj)
    dev.off()
  } else if (!is.null(draw_fn)) {
    cairo_pdf(path, width = w_in, height = h_in, family = FONT)
    draw_fn()
    dev.off()
  }
  message("  Saved: ", path)
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

# ---------------------------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------------------------
message("=== Loading data ===")
load(file.path(DATA_PROC, "consensus_clustering_results.RData"))
load(file.path(DATA_PROC, "immune_analysis_results.RData"))
load(file.path(DATA_PROC, "pathway_analysis_results.RData"))
load(file.path(DATA_PROC, "clinical_characterization_results.RData"))

# Sample-to-subtype map from cluster_df
subtype_map <- setNames(cluster_df$UPR_subtype, cluster_df$barcode)

# Restrict the canonical order to subtypes actually present (K=2 drops "UPR-intermediate"),
# so every downstream table()/factor()/column_title built from SUBTYPE_ORDER stays length-consistent.
SUBTYPE_ORDER <- intersect(SUBTYPE_ORDER, unique(as.character(cluster_df$UPR_subtype)))

# Samples present in immune_ssgsea
common_samples <- intersect(colnames(immune_ssgsea), cluster_df$barcode)
subtypes       <- factor(subtype_map[common_samples], levels = SUBTYPE_ORDER)

# ---------------------------------------------------------------------------
# PANEL A — ssGSEA immune cell heatmap
# ---------------------------------------------------------------------------
message("=== Panel A: ssGSEA immune heatmap ===")

# Order samples by subtype
sample_order <- order(subtypes)
ordered_samples  <- common_samples[sample_order]
ordered_subtypes <- subtypes[sample_order]

# Scale each immune cell row (z-score) for visual contrast
ssgsea_mat <- immune_ssgsea[, ordered_samples]
ssgsea_z   <- t(scale(t(ssgsea_mat)))
ssgsea_z[ssgsea_z > 3]  <- 3
ssgsea_z[ssgsea_z < -3] <- -3

# Mark significant rows (KW FDR < 0.05)
sig_cells  <- ssgsea_kw$CellType[!is.na(ssgsea_kw$kw_padj) & ssgsea_kw$kw_padj < 0.05]
row_sig_vec <- ifelse(rownames(ssgsea_z) %in% sig_cells, "FDR < 0.05", "ns")

# Row annotation: significance + broader cell category
broad_cat <- dplyr::case_when(
  grepl("CD4|CD8|T helper|Regulatory|T follicular|gamma delta|Natural killer T|memory CD",
        rownames(ssgsea_z)) ~ "T/NK cell",
  grepl("Natural killer cell", rownames(ssgsea_z)) ~ "T/NK cell",
  grepl("B cell|B lineage|Immature B", rownames(ssgsea_z)) ~ "B cell",
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
  annotation_name_gp   = gpar(fontsize = 6, fontfamily = FONT),
  simple_anno_size     = unit(3, "mm"),
  annotation_legend_param = list(
    Significant = list(title_gp  = gpar(fontsize = 7, fontfamily = FONT),
                       labels_gp = gpar(fontsize = 6, fontfamily = FONT)),
    Category    = list(title_gp  = gpar(fontsize = 7, fontfamily = FONT),
                       labels_gp = gpar(fontsize = 6, fontfamily = FONT))
  )
)

# Subtype counts for column split titles
st_n   <- table(ordered_subtypes)[SUBTYPE_ORDER]
ct_ttl <- paste0(SUBTYPE_ORDER, "\n(n=", st_n, ")")

col_anno <- HeatmapAnnotation(
  `UPR subtype` = ordered_subtypes,
  col = list(`UPR subtype` = COL_SUBTYPE),
  annotation_name_gp = gpar(fontsize = 6, fontfamily = FONT),
  annotation_name_side = "left",
  simple_anno_size = unit(3, "mm"),
  annotation_legend_param = list(
    `UPR subtype` = list(title_gp  = gpar(fontsize = 7, fontfamily = FONT),
                          labels_gp = gpar(fontsize = 6, fontfamily = FONT))
  )
)

ht_immune <- Heatmap(
  ssgsea_z,
  name            = "ssGSEA\n(z-score)",
  col             = colorRamp2(c(-3, 0, 3), c("#2166AC", "white", "#B2182B")),
  top_annotation  = col_anno,
  right_annotation = row_anno,
  cluster_columns = FALSE,
  cluster_rows    = TRUE,
  clustering_method_rows = "ward.D2",
  show_column_names = FALSE,
  row_names_gp    = gpar(fontsize = 6, fontfamily = FONT),
  column_split    = factor(ordered_subtypes, levels = SUBTYPE_ORDER),
  column_title    = ct_ttl,
  column_title_gp = gpar(fontsize = 7, fontface = "bold", fontfamily = FONT),
  heatmap_legend_param = list(
    title_gp  = gpar(fontsize = 7, fontfamily = FONT),
    labels_gp = gpar(fontsize = 6, fontfamily = FONT),
    legend_height = unit(25, "mm")
  ),
  use_raster = FALSE
)

save_pdf(path = file.path(OUT_DIR, "Fig4A_immune_heatmap.pdf"),
         width_mm = W_DOUBLE, height_mm = 120,
         draw_fn = function() {
           draw(ht_immune,
                merge_legend = TRUE,
                heatmap_legend_side = "right",
                annotation_legend_side = "right")
         })

# ---------------------------------------------------------------------------
# PANEL B — Top immune cell boxplots (FDR < 0.05, top 12 by padj)
# ---------------------------------------------------------------------------
message("=== Panel B: Immune cell boxplots ===")

# Pick top immune cells by significance
top_cells <- ssgsea_kw %>%
  dplyr::filter(!is.na(kw_padj), kw_padj < 0.05) %>%
  dplyr::arrange(kw_padj) %>%
  dplyr::slice_head(n = 12) %>%
  dplyr::pull(CellType)

if (length(top_cells) == 0) top_cells <- sig_cells[1:min(12, length(sig_cells))]

immune_long <- as.data.frame(t(immune_ssgsea[top_cells, common_samples])) %>%
  mutate(barcode = rownames(.), UPR_subtype = subtypes) %>%
  pivot_longer(-c(barcode, UPR_subtype), names_to = "CellType", values_to = "Score") %>%
  mutate(
    CellType = factor(CellType, levels = top_cells),
    UPR_subtype = factor(UPR_subtype, levels = SUBTYPE_ORDER)
  )

# Keep only comparisons whose both groups are present (K=2 -> single favorable vs high-risk)
comparisons_3 <- Filter(function(p) all(p %in% SUBTYPE_ORDER), list(
  c("UPR-high-risk", "UPR-intermediate"),
  c("UPR-high-risk", "UPR-favorable"),
  c("UPR-intermediate", "UPR-favorable")
))

p_immune_box <- ggplot(immune_long,
                        aes(x = UPR_subtype, y = Score, fill = UPR_subtype)) +
  geom_violin(alpha = 0.5, trim = TRUE, linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.35,
               fill = "white", color = "black") +
  stat_compare_means(comparisons = comparisons_3,
                     method = "wilcox.test", p.adjust.method = "BH",
                     label = "p.signif", hide.ns = TRUE,
                     size = 2, tip.length = 0.01,
                     step.increase = 0.07,
                     family = FONT) +
  scale_fill_manual(values = COL_SUBTYPE, guide = "none") +
  scale_x_discrete(labels = c("High-risk", "Interm.", "Favorable")) +
  facet_wrap(~CellType, scales = "free_y", ncol = 4) +
  labs(x = "", y = "ssGSEA enrichment score",
       title = "Immune cell infiltration by UPR subtype (FDR < 0.05)") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 6),
        strip.text  = element_text(size = 6.5, face = "bold"))

n_rows_b <- ceiling(length(top_cells) / 4)
save_pdf(p_immune_box, file.path(OUT_DIR, "Fig4B_immune_boxplots.pdf"),
         W_DOUBLE, max(70, n_rows_b * 42))

# ---------------------------------------------------------------------------
# PANEL C — Checkpoint genes: CD274, CTLA4, LAG3, HAVCR2 (TIM3), TIGIT
# ---------------------------------------------------------------------------
message("=== Panel C: Checkpoint genes ===")

checkpoint_show <- c("CD274", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "PDCD1")
checkpoint_show <- checkpoint_show[checkpoint_show %in% icp_kw$Gene]

# Load expression
load(file.path(DATA_PROC, "tcga_glioma_expression.RData"))
expr_log <- log2(expr_tpm_symbol[, common_samples] + 1)

icp_expr <- expr_log[checkpoint_show[checkpoint_show %in% rownames(expr_log)], ]

icp_long_show <- as.data.frame(t(icp_expr)) %>%
  mutate(barcode = rownames(.), UPR_subtype = subtypes) %>%
  pivot_longer(-c(barcode, UPR_subtype),
               names_to = "Gene", values_to = "Expression") %>%
  mutate(
    Gene = factor(Gene, levels = checkpoint_show),
    UPR_subtype = factor(UPR_subtype, levels = SUBTYPE_ORDER),
    Gene_label = dplyr::recode(Gene,
                                "CD274"  = "PD-L1 (CD274)",
                                "HAVCR2" = "TIM3 (HAVCR2)",
                                "PDCD1"  = "PD-1 (PDCD1)")
  )

# Add padj annotation
icp_padj <- icp_kw[icp_kw$Gene %in% checkpoint_show, c("Gene", "kw_padj")]
icp_long_show <- left_join(icp_long_show, icp_padj, by = "Gene") %>%
  mutate(panel_label = paste0(Gene_label, "\n", sig_star(kw_padj)))

p_ckpt <- ggplot(icp_long_show,
                  aes(x = UPR_subtype, y = Expression, fill = UPR_subtype)) +
  geom_violin(alpha = 0.5, trim = TRUE, linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.35,
               fill = "white", color = "black") +
  stat_compare_means(comparisons = comparisons_3,
                     method = "wilcox.test", p.adjust.method = "BH",
                     label = "p.signif", hide.ns = TRUE,
                     size = 2, tip.length = 0.01, step.increase = 0.08,
                     family = FONT) +
  scale_fill_manual(values = COL_SUBTYPE, guide = "none") +
  scale_x_discrete(labels = c("High-risk", "Interm.", "Favorable")) +
  facet_wrap(~panel_label, scales = "free_y", nrow = 2) +
  labs(x = "", y = "Expression (log2 TPM + 1)",
       title = "Immune checkpoint gene expression by UPR subtype") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 6),
        strip.text  = element_text(size = 6.5, face = "bold"))

save_pdf(p_ckpt, file.path(OUT_DIR, "Fig4C_checkpoints.pdf"),
         W_DOUBLE, 80)

# ---------------------------------------------------------------------------
# PANEL D — ESTIMATE scores
# ---------------------------------------------------------------------------
message("=== Panel D: ESTIMATE scores ===")

# Fix barcode format mismatch (est_df uses dots, cluster_df uses dashes)
est_fixed <- est_df
est_fixed$barcode <- gsub("\\.", "-", est_fixed$barcode)

# Merge with subtype from cluster_df directly
est_merged <- merge(
  est_fixed[, c("barcode", "StromalScore", "ImmuneScore", "ESTIMATEScore")],
  cluster_df[, c("barcode", "UPR_subtype")],
  by = "barcode"
) %>%
  dplyr::filter(!is.na(UPR_subtype)) %>%
  mutate(UPR_subtype = factor(UPR_subtype, levels = SUBTYPE_ORDER))

if (nrow(est_merged) == 0) {
  # Fallback: use heatmap_data which already has UPR_subtype
  est_merged <- merge(
    est_fixed[, c("barcode", "StromalScore", "ImmuneScore", "ESTIMATEScore")],
    heatmap_data[, c("barcode", "UPR_subtype")],
    by = "barcode"
  ) %>%
    dplyr::filter(!is.na(UPR_subtype)) %>%
    mutate(UPR_subtype = factor(UPR_subtype, levels = SUBTYPE_ORDER))
}

message(sprintf("  ESTIMATE merged rows: %d", nrow(est_merged)))

est_long <- est_merged %>%
  pivot_longer(cols = c("StromalScore", "ImmuneScore", "ESTIMATEScore"),
               names_to = "Score_type", values_to = "Score") %>%
  mutate(Score_type = dplyr::recode(Score_type,
    "StromalScore"   = "Stromal score",
    "ImmuneScore"    = "Immune score",
    "ESTIMATEScore"  = "ESTIMATE score"
  ),
  Score_type = factor(Score_type,
                       levels = c("Immune score", "Stromal score", "ESTIMATE score")))

# KW tests per score type
est_kw <- est_long %>%
  dplyr::group_by(Score_type) %>%
  dplyr::summarise(
    kw_p = tryCatch(kruskal.test(Score ~ UPR_subtype)$p.value, error = function(e) NA),
    .groups = "drop"
  )

est_long <- left_join(est_long, est_kw, by = "Score_type") %>%
  mutate(panel_label = paste0(Score_type, "\n(", sig_star(kw_p), ")"))

p_estimate <- ggplot(est_long,
                      aes(x = UPR_subtype, y = Score, fill = UPR_subtype)) +
  geom_violin(alpha = 0.5, trim = TRUE, linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.35,
               fill = "white", color = "black") +
  stat_compare_means(comparisons = comparisons_3,
                     method = "wilcox.test", p.adjust.method = "BH",
                     label = "p.signif", hide.ns = TRUE,
                     size = 2, tip.length = 0.01, step.increase = 0.08,
                     family = FONT) +
  scale_fill_manual(values = COL_SUBTYPE, guide = "none") +
  scale_x_discrete(labels = c("High-risk", "Interm.", "Favorable")) +
  facet_wrap(~panel_label, scales = "free_y", nrow = 1) +
  labs(x = "", y = "ESTIMATE score",
       title = "Tumor microenvironment scores by UPR subtype") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, size = 6),
        strip.text  = element_text(size = 7, face = "bold"))

save_pdf(p_estimate, file.path(OUT_DIR, "Fig4D_estimate.pdf"),
         W_DOUBLE, 65)

# ---------------------------------------------------------------------------
# PANEL E — GSVA Hallmark heatmap (all 50 pathways, row-scaled)
# ---------------------------------------------------------------------------
message("=== Panel E: GSVA Hallmark heatmap ===")

# Order samples
gsva_samples <- intersect(colnames(gsva_hallmark), common_samples)
gsva_subtypes <- factor(subtype_map[gsva_samples], levels = SUBTYPE_ORDER)
gsva_order <- order(gsva_subtypes)
gsva_ordered <- gsva_hallmark[, gsva_samples[gsva_order]]
gsva_subtypes_ordered <- gsva_subtypes[gsva_order]

# Row-scale
gsva_z <- t(scale(t(gsva_ordered)))
gsva_z[gsva_z > 2]  <- 2
gsva_z[gsva_z < -2] <- -2

# Simplify pathway names
pathway_names_short <- gsub("HALLMARK_", "", rownames(gsva_z))
pathway_names_short <- gsub("_", " ", pathway_names_short)
pathway_names_short <- tolower(pathway_names_short)
pathway_names_short <- paste0(toupper(substr(pathway_names_short, 1, 1)),
                               substr(pathway_names_short, 2, nchar(pathway_names_short)))
rownames(gsva_z) <- pathway_names_short

# Mark differentially active pathways using the primary pairwise contrast present.
# At K=2 this is the single favorable-vs-high-risk contrast; positive logFC = higher
# in the first (alphabetically-earlier) group, i.e. UPR-favorable.
aq_res <- diff_pathway_results[[1]]
aq_res$pathway_short <- gsub("HALLMARK_", "", aq_res$pathway)
aq_res$pathway_short <- gsub("_", " ", aq_res$pathway_short)
aq_res$pathway_short <- tolower(aq_res$pathway_short)
aq_res$pathway_short <- paste0(toupper(substr(aq_res$pathway_short, 1, 1)),
                                substr(aq_res$pathway_short, 2, nchar(aq_res$pathway_short)))

sig_paths <- aq_res$pathway_short[aq_res$adj.P.Val < 0.05]
row_diff <- ifelse(rownames(gsva_z) %in% sig_paths, "FDR < 0.05", "ns")

# Direction (higher in favorable = positive logFC)
row_dir <- dplyr::case_when(
  !(rownames(gsva_z) %in% sig_paths) ~ "ns",
  aq_res$logFC[match(rownames(gsva_z), aq_res$pathway_short)] > 0 ~ "Higher in Favorable",
  TRUE ~ "Higher in High-risk"
)

dir_colors <- c(
  "Higher in Favorable"  = "#4DBBD5",
  "Higher in High-risk"  = "#E64B35",
  "ns"                   = "grey85"
)

row_anno_gsva <- rowAnnotation(
  Direction = row_dir,
  col = list(Direction = dir_colors),
  annotation_name_gp   = gpar(fontsize = 6, fontfamily = FONT),
  simple_anno_size     = unit(3, "mm"),
  annotation_legend_param = list(
    Direction = list(title_gp  = gpar(fontsize = 7, fontfamily = FONT),
                     labels_gp = gpar(fontsize = 6, fontfamily = FONT))
  )
)

st_n_gsva  <- table(gsva_subtypes_ordered)[SUBTYPE_ORDER]
ct_ttl_gsva <- paste0(SUBTYPE_ORDER, "\n(n=", st_n_gsva, ")")

col_anno_gsva <- HeatmapAnnotation(
  `UPR subtype` = gsva_subtypes_ordered,
  col = list(`UPR subtype` = COL_SUBTYPE),
  annotation_name_gp  = gpar(fontsize = 6, fontfamily = FONT),
  annotation_name_side = "left",
  simple_anno_size    = unit(3, "mm"),
  annotation_legend_param = list(
    `UPR subtype` = list(title_gp  = gpar(fontsize = 7, fontfamily = FONT),
                          labels_gp = gpar(fontsize = 6, fontfamily = FONT))
  )
)

ht_gsva <- Heatmap(
  gsva_z,
  name             = "GSVA\n(z-score)",
  col              = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
  top_annotation   = col_anno_gsva,
  right_annotation = row_anno_gsva,
  cluster_columns  = FALSE,
  cluster_rows     = TRUE,
  clustering_method_rows = "ward.D2",
  show_column_names = FALSE,
  row_names_gp     = gpar(fontsize = 5.5, fontfamily = FONT),
  column_split     = factor(gsva_subtypes_ordered, levels = SUBTYPE_ORDER),
  column_title     = ct_ttl_gsva,
  column_title_gp  = gpar(fontsize = 7, fontface = "bold", fontfamily = FONT),
  heatmap_legend_param = list(
    title_gp  = gpar(fontsize = 7, fontfamily = FONT),
    labels_gp = gpar(fontsize = 6, fontfamily = FONT),
    legend_height = unit(20, "mm")
  ),
  use_raster = FALSE
)

save_pdf(path = file.path(OUT_DIR, "Fig4E_gsva_heatmap.pdf"),
         width_mm = W_DOUBLE, height_mm = 130,
         draw_fn = function() {
           draw(ht_gsva,
                merge_legend = TRUE,
                heatmap_legend_side = "right",
                annotation_legend_side = "right")
         })

# ---------------------------------------------------------------------------
# PANEL F — GSEA GO-BP dotplot (top 10 favorable-enriched + top 10 high-risk-enriched)
# ---------------------------------------------------------------------------
message("=== Panel F: GSEA GO-BP dotplot ===")

gsea_res <- gsea_go@result

# Top 10 favorable-enriched (NES > 0, lowest p.adjust)
# Top 10 high-risk-enriched (NES < 0, lowest p.adjust)
top_fav <- gsea_res %>%
  dplyr::filter(NES > 0) %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::slice_head(n = 10) %>%
  mutate(Direction = "UPR-favorable enriched")

top_hr <- gsea_res %>%
  dplyr::filter(NES < 0) %>%
  dplyr::arrange(p.adjust) %>%
  dplyr::slice_head(n = 10) %>%
  mutate(Direction = "UPR-high-risk enriched")

gsea_top <- bind_rows(top_fav, top_hr) %>%
  mutate(
    term_short = gsub("GOBP_", "", ID),
    term_short = gsub("_", " ", term_short),
    term_short = tolower(term_short),
    term_short = paste0(toupper(substr(term_short, 1, 1)),
                         substr(term_short, 2, nchar(term_short))),
    term_short = ifelse(nchar(term_short) > 42,
                         paste0(substr(term_short, 1, 40), "\u2026"),
                         term_short),
    neg_log10_padj = -log10(p.adjust),
    Direction = factor(Direction,
                        levels = c("UPR-high-risk enriched", "UPR-favorable enriched"))
  )

# Within each direction, order by |NES|
gsea_top <- gsea_top %>%
  dplyr::arrange(Direction, NES) %>%
  mutate(term_short = factor(term_short, levels = unique(term_short)))

p_gsea_dot <- ggplot(gsea_top,
                      aes(x = NES, y = term_short,
                          size = neg_log10_padj,
                          color = Direction)) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "grey60", linetype = "dashed") +
  geom_point(alpha = 0.85) +
  scale_color_manual(values = c(
    "UPR-favorable enriched"  = "#4DBBD5",
    "UPR-high-risk enriched"  = "#E64B35"
  )) +
  scale_size_continuous(name = expression(-log[10]("FDR")),
                         range = c(1.5, 5),
                         breaks = c(5, 10, 20)) +
  scale_x_continuous(expand = expansion(mult = 0.1)) +
  facet_grid(Direction ~ ., scales = "free_y", space = "free_y") +
  labs(x = "Normalized enrichment score (NES)",
       y = "",
       color = "Enrichment direction",
       title = "GSEA GO:BP (UPR-favorable vs UPR-high-risk)") +
  theme_pub() +
  theme(
    legend.position  = "right",
    strip.text       = element_text(size = BASE_SIZE - 0.5, face = "bold"),
    axis.text.y      = element_text(size = 5.5),
    panel.border     = element_rect(color = "grey80", fill = NA, linewidth = 0.3)
  )

save_pdf(p_gsea_dot, file.path(OUT_DIR, "Fig4F_gsea_dotplot.pdf"),
         W_DOUBLE, 100)

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
message("\n=== Figure 4 complete ===")
cat("\nOutput directory:", OUT_DIR, "\n")
cat("Files produced:\n")
for (f in list.files(OUT_DIR, pattern = "\\.pdf$")) {
  cat(" ", f, "\n")
}

cat("\nKey statistics for caption:\n")
cat("  ssGSEA significant immune cells (FDR<0.05):", length(sig_cells), "/",
    nrow(ssgsea_kw), "\n")
cat("  ESTIMATE merged samples:", nrow(est_merged), "\n")
cat("  GSVA significant pathways:", sum(!is.na(aq_res$adj.P.Val) & aq_res$adj.P.Val < 0.05),
    "/ 50 Hallmark\n")
cat("  GSEA GO:BP enriched terms:", nrow(gsea_res),
    "(favorable:", sum(gsea_res$NES > 0), ", high-risk:", sum(gsea_res$NES < 0), ")\n")
