###############################################################################
# fig1_composite_v2.R
# Figure 1 composite — single-cell UPR landscape
# Panels: A UPR-score UMAP, B IRE1-XBP1 UMAP, C PERK-ATF4 UMAP, D ATF6 UMAP,
#         E cell-type violin, F pathway dotplot,
#         G malignant-subtype violin, H proportion stacked bar
# Fixes: arm titles use "IRE1-XBP1", "PERK-ATF4", "ATF6"; subtype axis "OPC-like"
# Output: Figures_v2/Figure1_composite.pdf/.png
###############################################################################


suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(showtext)
  library(sysfonts)
})

font_add("Times",
         regular    = "/usr/share/fonts/truetype/dejavu/DejaVuSerif.ttf",
         bold       = "/usr/share/fonts/truetype/dejavu/DejaVuSerif-Bold.ttf",
         italic     = "/usr/share/fonts/truetype/dejavu/DejaVuSerif-Italic.ttf",
         bolditalic = "/usr/share/fonts/truetype/dejavu/DejaVuSerif-BoldItalic.ttf")
showtext_auto()
showtext_opts(dpi = 300)

PROJECT_DIR <- getwd()
DATA_PROC   <- file.path(PROJECT_DIR, "data", "processed")
OUT_DIR     <- file.path(PROJECT_DIR, "manuscript", "submission", "Figures_v2")

# Publication theme (8 pt base — legible at 180 mm final width)
BASE_SIZE <- 8
THEME_PUB <- theme_classic(base_size = BASE_SIZE, base_family = "Times") +
  theme(
    plot.title   = element_text(size = BASE_SIZE, face = "bold", hjust = 0.5),
    axis.title   = element_text(size = BASE_SIZE),
    axis.text    = element_text(size = BASE_SIZE - 1, color = "black"),
    legend.title = element_text(size = BASE_SIZE - 1, face = "bold"),
    legend.text  = element_text(size = BASE_SIZE - 2),
    legend.key.size = unit(3, "mm"),
    strip.text   = element_text(size = BASE_SIZE - 1, face = "bold"),
    strip.background = element_blank(),
    panel.grid   = element_blank(),
    plot.margin  = margin(2, 2, 2, 2, "mm")
  )

COLORS_CELLTYPE <- c(
  "Malignant"       = "#E64B35", "TAM"             = "#4DBBD5",
  "Microglia"       = "#00A087", "CD4 T"           = "#3C5488",
  "CD8 T"           = "#F39B7F", "Treg"            = "#8491B4",
  "DC"              = "#91D1C2", "Oligodendrocyte" = "#B09C85",
  "Astrocyte"       = "#7E6148", "Endothelial"     = "#FFDC91",
  "MDSC"            = "#DC0000"
)

# ── Load data ─────────────────────────────────────────────────────────────────
message("Loading seu_upr_scored.rds ...")
seu <- readRDS(file.path(DATA_PROC, "seu_upr_scored.rds"))
message(sprintf("  %d genes x %d cells", nrow(seu), ncol(seu)))

all_celltypes <- unique(seu$celltype)
color_vec <- COLORS_CELLTYPE[all_celltypes]
missing_ct <- all_celltypes[is.na(color_vec)]
if (length(missing_ct) > 0) {
  extra_colors <- scales::hue_pal()(length(missing_ct))
  names(extra_colors) <- missing_ct
  color_vec[missing_ct] <- extra_colors
}
color_vec <- color_vec[!is.na(names(color_vec))]

###############################################################################
# Panel A — UPR score UMAP
###############################################################################
p_A <- FeaturePlot(seu, features = "UPR_score",
                   cols = c("lightgrey", "darkred"),
                   order = TRUE, combine = FALSE)[[1]] +
  THEME_PUB +
  theme(axis.line  = element_blank(), axis.ticks = element_blank(),
        axis.text  = element_blank(), axis.title = element_text(size = BASE_SIZE - 2)) +
  labs(title = "UPR Activity Score", x = "UMAP 1", y = "UMAP 2")

###############################################################################
# Panels B, C, D — Three UPR arm UMAPs
# FIX: use readable titles with dash, not underscore
###############################################################################
arm_features <- c("AUCell_IRE1_XBP1", "AUCell_PERK_ATF4", "AUCell_ATF6")
arm_features <- arm_features[arm_features %in% colnames(seu@meta.data)]
arm_titles   <- c("IRE1-XBP1", "PERK-ATF4", "ATF6")[seq_along(arm_features)]

arm_plots <- lapply(seq_along(arm_features), function(i) {
  f <- arm_features[i]
  t <- arm_titles[i]
  pl <- FeaturePlot(seu, features = f,
                    cols = c("lightgrey", "darkred"),
                    order = TRUE, combine = FALSE)[[1]]
  pl + THEME_PUB +
    theme(axis.line = element_blank(), axis.ticks = element_blank(),
          axis.text = element_blank(), axis.title = element_text(size = BASE_SIZE - 2)) +
    labs(title = t, x = "UMAP 1", y = "UMAP 2")
})

###############################################################################
# Panel E — Cell-type violin (UPR score)
###############################################################################
vln_df <- data.frame(UPR_score = seu$UPR_score, celltype = seu$celltype)

p_E <- ggplot(vln_df, aes(x = celltype, y = UPR_score, fill = celltype)) +
  geom_violin(scale = "width", trim = FALSE) +
  geom_boxplot(width = 0.12, fill = "white", outlier.size = 0.2) +
  scale_fill_manual(values = color_vec) +
  THEME_PUB +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = BASE_SIZE - 2),
        legend.position = "none") +
  labs(x = NULL, y = "UPR activity score",
       title = "UPR score by cell type")

###############################################################################
# Panel F — Pathway dotplot (UPR arm × cell type)
###############################################################################
aucell_cols <- grep("^AUCell_", colnames(seu@meta.data), value = TRUE)

if (length(aucell_cols) >= 3) {
  # Compute dot data
  dot_data <- seu@meta.data %>%
    dplyr::select(celltype, all_of(aucell_cols)) %>%
    tidyr::pivot_longer(cols = -celltype, names_to = "Pathway", values_to = "Score") %>%
    dplyr::mutate(
      # FIX: readable pathway labels with dashes
      Pathway = dplyr::recode(Pathway,
        "AUCell_IRE1_XBP1" = "IRE1-XBP1",
        "AUCell_PERK_ATF4" = "PERK-ATF4",
        "AUCell_ATF6"      = "ATF6",
        "AUCell_UPR_broad" = "UPR broad"
      )
    )

  dot_summary <- dot_data %>%
    dplyr::group_by(celltype, Pathway) %>%
    dplyr::summarise(
      mean_score = mean(Score, na.rm = TRUE),
      .groups = "drop"
    )

  # Heatmap: sequential fill encodes mean AUCell score (unidirectional, >= 0).
  # Single colour/intensity channel only; no size channel (pct-active removed).
  p_F <- ggplot(dot_summary, aes(x = Pathway, y = celltype, fill = mean_score)) +
    geom_tile(colour = "white", linewidth = 0.5) +
    scale_fill_gradient(low = "white", high = "#E64B35", name = "Mean AUC") +
    labs(x = NULL, y = NULL,
         title = "UPR Pathway Activity by Cell Type") +
    THEME_PUB +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = BASE_SIZE - 1),
          panel.grid.major = element_blank(),
          legend.position = "right",
          plot.margin = margin(2, 6, 2, 2, "mm")) +
    scale_x_discrete(expand = expansion(add = 0)) +
    scale_y_discrete(expand = expansion(add = 0)) +
    coord_cartesian(clip = "off")
} else {
  p_F <- NULL
}

###############################################################################
# Panel G — Malignant subtypes violin
# FIX: replace underscores in subtype labels with dashes/spaces
###############################################################################
p_G <- NULL
if ("malignant_subtype" %in% colnames(seu@meta.data) &&
    sum(!is.na(seu$malignant_subtype)) > 20) {

  mal_vln_df <- seu@meta.data %>%
    dplyr::filter(celltype == "Malignant" & !is.na(malignant_subtype)) %>%
    dplyr::select(UPR_score, subtype = malignant_subtype) %>%
    dplyr::mutate(
      # FIX: replace underscores with dashes for display
      subtype_label = gsub("_like$", "-like", subtype),
      subtype_label = gsub("_", "-", subtype_label)
    )

  subtype_colors <- c(
    "MES_like" = "#E64B35", "AC_like" = "#4DBBD5",
    "OPC_like" = "#00A087", "NPC_like" = "#3C5488",
    # also map the label versions for scale_fill
    "MES-like" = "#E64B35", "AC-like" = "#4DBBD5",
    "OPC-like" = "#00A087", "NPC-like" = "#3C5488"
  )

  if (nrow(mal_vln_df) >= 20) {
    p_G <- ggplot(mal_vln_df, aes(x = subtype_label, y = UPR_score, fill = subtype_label)) +
      geom_violin(scale = "width", trim = FALSE) +
      geom_boxplot(width = 0.12, fill = "white", outlier.size = 0.2) +
      scale_fill_manual(values = subtype_colors) +
      THEME_PUB +
      theme(legend.position = "none",
            axis.text.x = element_text(angle = 30, hjust = 1, size = BASE_SIZE - 1)) +
      labs(x = NULL, y = "UPR Score",
           title = "UPR Score in Malignant Subtypes")
  }
}

###############################################################################
# Panel H — UPR proportion stacked bar
###############################################################################
global_median <- median(seu$UPR_score, na.rm = TRUE)
seu$UPR_group <- ifelse(seu$UPR_score >= global_median, "UPR-high", "UPR-low")

prop_df <- data.frame(celltype = seu$celltype, UPR_group = seu$UPR_group) %>%
  dplyr::filter(!is.na(UPR_group), !is.na(celltype)) %>%
  dplyr::count(celltype, UPR_group, name = "n") %>%
  tidyr::complete(celltype, UPR_group, fill = list(n = 0)) %>%
  dplyr::group_by(celltype) %>%
  dplyr::mutate(pct = n / sum(n) * 100) %>%
  dplyr::ungroup()

ct_order <- prop_df %>%
  dplyr::filter(UPR_group == "UPR-high") %>%
  dplyr::arrange(dplyr::desc(pct)) %>%
  dplyr::pull(celltype) %>% as.character()
all_ct <- unique(as.character(prop_df$celltype))
ct_order <- c(ct_order, setdiff(all_ct, ct_order))
prop_df$celltype <- factor(prop_df$celltype, levels = ct_order)

p_H <- ggplot(prop_df, aes(x = celltype, y = pct, fill = UPR_group)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = c("UPR-high" = "#E64B35", "UPR-low" = "#4DBBD5"),
                    name = "UPR Group") +
  labs(x = NULL, y = "Percentage (%)", title = "UPR-high/low Proportion by Cell Type") +
  THEME_PUB +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = BASE_SIZE - 2))

###############################################################################
# Assemble composite
# Layout: 4 UMAPs (2x2) top half; E, F bottom-left; G, H bottom-right
# 4 rows x 2 cols = 8 cells at 183 mm wide, ~240 mm tall (just fits one page)
###############################################################################
message("=== Assembling composite ===")

panels_list <- list(p_A)
panels_list <- c(panels_list, arm_plots[1])         # row 1: A + B
if (length(arm_plots) >= 2) panels_list <- c(panels_list, arm_plots[2]) else panels_list <- c(panels_list, list(plot_spacer()))
if (length(arm_plots) >= 3) panels_list <- c(panels_list, arm_plots[3]) else panels_list <- c(panels_list, list(plot_spacer()))
# row 2: C + D done above — we have 4 UMAPs in 2x2
# E and F fill row 3
if (!is.null(p_F)) {
  panels_list <- c(panels_list, list(p_E, p_F))
} else {
  panels_list <- c(panels_list, list(p_E, plot_spacer()))
}
# G and H fill row 4
if (!is.null(p_G)) {
  panels_list <- c(panels_list, list(p_G, p_H))
} else {
  panels_list <- c(panels_list, list(plot_spacer(), p_H))
}

# Wrap: 4 rows x 2 cols = 8 panels
p_composite <- wrap_plots(panels_list, ncol = 2, nrow = 4) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size = 13, face = "bold", family = "Times",
                              color = "black", margin = margin(0,2,2,0,"mm"))
    )
  )

# Dimensions: 183 mm wide, 240 mm tall (4 rows × 60 mm each)
W_MM <- 183
H_MM <- 240

out_pdf <- file.path(OUT_DIR, "Figure1_composite.pdf")
out_png <- file.path(OUT_DIR, "Figure1_composite.png")

ggsave(out_pdf, plot = p_composite, device = cairo_pdf,
       width = W_MM, height = H_MM, units = "mm", dpi = 300)
message("Saved PDF: ", out_pdf)

ggsave(out_png, plot = p_composite,
       width = W_MM, height = H_MM, units = "mm", dpi = 300)
message("Saved PNG: ", out_png)

message("=== Figure 1 composite done ===")
message(sprintf("Dimensions: %d x %d mm  (W/H aspect = %.2f)", W_MM, H_MM, W_MM/H_MM))
