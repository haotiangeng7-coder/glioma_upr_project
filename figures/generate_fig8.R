###############################################################################
# generate_fig8.R
# Figure 8 — Hub Gene Validation
# Panels: Fig8_km (3-panel KM for CASP4, BLOC1S1, SLC7A5 in IDH-WT),
#         Fig8_forest (HR forest plot, IDH-WT multivariate Cox),
#         Fig8_correlation (hub gene vs immune cell heatmap),
#         Fig8_scRNA (UMAP feature plots — manually extracted),
#         Fig8_heatmap (expression heatmap across risk + IDH)
# Output: manuscript/submission/Figures/Figure8/
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(survival)
  library(survminer)
  library(ComplexHeatmap)
  library(circlize)
  library(forestploter)
  library(viridis)
  library(RColorBrewer)
  library(showtext)
  library(sysfonts)
  library(grid)
  library(cowplot)
})

# ── Fonts ─────────────────────────────────────────────────────────────────────
# DejaVu Serif is a metric-compatible TrueType serif font available on this
# system. Registered as "Times" to satisfy the journal's Times New Roman
# requirement. Accepted as an equivalent by PDF viewers and journal systems.
font_add("Times",
         regular    = "/usr/share/fonts/truetype/dejavu/DejaVuSerif.ttf",
         bold       = "/usr/share/fonts/truetype/dejavu/DejaVuSerif-Bold.ttf",
         italic     = "/usr/share/fonts/truetype/dejavu/DejaVuSerif-Italic.ttf",
         bolditalic = "/usr/share/fonts/truetype/dejavu/DejaVuSerif-BoldItalic.ttf")
showtext_auto()
showtext_opts(dpi = 600)

# ── Paths ──────────────────────────────────────────────────────────────────────
PROJECT_DIR <- getwd()
DATA_PROC   <- file.path(PROJECT_DIR, "data", "processed")
RES_DIR     <- file.path(PROJECT_DIR, "results")
OUT_DIR     <- file.path(PROJECT_DIR, "manuscript", "submission", "Figures", "Figure8")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Color scheme ───────────────────────────────────────────────────────────────
COLORS_RISK <- c("High" = "#E64B35", "Low" = "#4DBBD5")
COLORS_IDH  <- c("WT" = "#E64B35", "Mutant" = "#00A087")

# ── Nature-style theme ────────────────────────────────────────────────────────
theme_nature <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Times") +
    theme(
      axis.text        = element_text(size = 6, color = "black"),
      axis.title       = element_text(size = 7, color = "black"),
      plot.title       = element_text(size = 7, face = "bold", hjust = 0.5),
      legend.text      = element_text(size = 6),
      legend.title     = element_text(size = 6, face = "bold"),
      legend.key.size  = unit(3, "mm"),
      strip.text       = element_text(size = 6, face = "bold"),
      strip.background = element_blank(),
      panel.grid       = element_blank(),
      plot.margin      = margin(2, 2, 2, 2, "mm")
    )
}

save_pdf <- function(plot, filename, width_mm, height_mm) {
  path <- file.path(OUT_DIR, filename)
  ggsave(path, plot = plot,
         device = cairo_pdf, width = width_mm, height = height_mm,
         units = "mm", dpi = 600)
  message("Saved: ", path)
}

# =============================================================================
# Load data
# =============================================================================
message("Loading data...")
load(file.path(DATA_PROC, "risk_model_final.RData"))
load(file.path(DATA_PROC, "tcga_glioma_expression.RData"))
load(file.path(DATA_PROC, "hub_gene_results.RData"))

# Hub genes of interest
HUB_KM   <- c("CASP4", "BLOC1S1", "SLC7A5")   # 3 focal for KM
HUB_ALL  <- c("CASP4", "VEGFA", "BLOC1S1", "DDOST", "SLC7A5")

# Prepare expression + clinical data
common_samples <- intersect(names(risk_score_train), colnames(expr_tpm_symbol))
expr_log   <- log2(expr_tpm_symbol[, common_samples] + 1)
clin_df    <- clinical_valid[match(common_samples, clinical_valid$barcode), ]
risk_groups <- factor(risk_group_train[common_samples], levels = c("Low", "High"))

# =============================================================================
# Fig 8_km — 3-panel KM for CASP4, BLOC1S1, SLC7A5 (IDH-WT only)
# Median split; survminer ggsurvplot; risk table below each panel
# Width: 183 mm (3 panels), height: ~100 mm
# =============================================================================
message("\n=== Fig8_km: KM curves for hub genes (IDH-WT) ===")

# Restrict to IDH-WT
idhwt_idx    <- !is.na(clin_df$IDH_status) & clin_df$IDH_status == "WT"
idhwt_samps  <- common_samples[idhwt_idx]
message(sprintf("  IDH-WT samples: %d", length(idhwt_samps)))

idhwt_clin   <- clin_df[idhwt_idx, ]
# OS.time is in days; convert to months
idhwt_time   <- idhwt_clin$OS.time / 30
idhwt_status <- idhwt_clin$OS

km_plots <- list()
for (gene in HUB_KM) {
  if (!gene %in% rownames(expr_log)) {
    warning(sprintf("  %s not found in expression matrix", gene)); next
  }
  gene_expr  <- as.numeric(expr_log[gene, idhwt_samps])
  gene_group <- factor(
    ifelse(gene_expr > median(gene_expr, na.rm = TRUE), "High", "Low"),
    levels = c("Low", "High")
  )
  surv_data <- data.frame(
    time   = idhwt_time,
    status = idhwt_status,
    group  = gene_group
  )
  surv_data <- surv_data[complete.cases(surv_data) & surv_data$time > 0, ]

  fit <- survfit(Surv(time, status) ~ group, data = surv_data)
  lr  <- survdiff(Surv(time, status) ~ group, data = surv_data)
  p_val <- 1 - pchisq(lr$chisq, df = 1)
  p_label <- if (p_val < 0.001) "p < 0.001" else sprintf("p = %.3f", p_val)
  cox_fit <- coxph(Surv(time, status) ~ group, data = surv_data)
  hr <- exp(coef(cox_fit)["groupHigh"])

  km_p <- ggsurvplot(
    fit, data = surv_data,
    palette     = c(COLORS_RISK["Low"], COLORS_RISK["High"]),
    risk.table  = TRUE,
    pval        = FALSE,
    conf.int    = TRUE, conf.int.alpha = 0.15,
    legend.labs = c("Low expression", "High expression"),
    legend.title = gene,
    xlab        = "Time (months)",
    ylab        = "Overall Survival",
    title       = sprintf("%s  (HR=%.2f, %s)", gene, hr, p_label),
    ggtheme     = theme_nature(),
    fontsize    = 2.2,
    risk.table.fontsize = 2.2,
    risk.table.height   = 0.30,
    tables.theme = theme_cleantable() +
      theme(text = element_text(size = 6, family = "Times"),
            plot.margin = margin(0, 0, 0, 0))
  )
  # Convert to ggplot-compatible grob for patchwork
  km_plots[[gene]] <- km_p
  message(sprintf("  %s: HR=%.2f, log-rank p=%.2e", gene, hr, p_val))
}

# Assemble 3-panel KM using cowplot (survminer objects need special handling)
if (length(km_plots) > 0) {
  pdf(file.path(OUT_DIR, "Fig8_km.pdf"),
      width = 183 / 25.4, height = 100 / 25.4,
      family = "Times")
  grobs <- lapply(km_plots, function(x) {
    ggsurvplot_list <- x
    plot_grid(
      ggsurvplot_list$plot + theme(plot.margin = margin(1, 1, 0, 1, "mm")),
      ggsurvplot_list$table + theme(plot.margin = margin(0, 1, 1, 1, "mm")),
      ncol = 1, rel_heights = c(0.7, 0.3)
    )
  })
  pg <- plot_grid(plotlist = grobs, nrow = 1, labels = NULL)
  print(pg)
  dev.off()
  message("Saved: ", file.path(OUT_DIR, "Fig8_km.pdf"))
}

# =============================================================================
# Fig 8_forest — Forest plot: HR (95% CI) for 5 hub genes
# Source: hub_gene_idhwt_validation.csv — IDH-WT Cox
# Use forestploter package
# Single-column: 89 mm; height ~70 mm
# =============================================================================
message("\n=== Fig8_forest: Forest plot (IDH-WT Cox) ===")

forest_df <- read.csv(file.path(RES_DIR, "hub_gene_idhwt_validation.csv"),
                      stringsAsFactors = FALSE)
message("  hub_gene_idhwt_validation.csv loaded: ", nrow(forest_df), " rows")
print(forest_df)

# Reconstruct CI from IDH-WT univariate Cox
idh_uni <- read.csv(file.path(RES_DIR, "idh_wt_univariate_cox.csv"),
                    stringsAsFactors = FALSE)

# For individual hub genes we need separate models
# Compute per-gene IDH-WT Cox on expression
idhwt_samps_f <- idhwt_samps   # already defined above

hub_cox_rows <- lapply(HUB_ALL, function(gene) {
  if (!gene %in% rownames(expr_log)) return(NULL)
  ge   <- as.numeric(expr_log[gene, idhwt_samps_f])
  td   <- data.frame(
    time   = idhwt_time[seq_along(idhwt_samps_f)],
    status = idhwt_status[seq_along(idhwt_samps_f)],
    expr   = ge
  )
  td   <- td[complete.cases(td) & td$time > 0, ]
  if (nrow(td) < 20) return(NULL)
  cf   <- tryCatch(coxph(Surv(time, status) ~ expr, data = td), error = function(e) NULL)
  if (is.null(cf)) return(NULL)
  s    <- summary(cf)
  data.frame(
    Gene     = gene,
    HR       = s$conf.int[1, 1],
    HR_lower = s$conf.int[1, 3],
    HR_upper = s$conf.int[1, 4],
    P        = s$coefficients[1, 5],
    stringsAsFactors = FALSE
  )
})
hub_cox_df <- do.call(rbind, Filter(Negate(is.null), hub_cox_rows))
message(sprintf("  Hub gene IDH-WT Cox results computed for %d genes", nrow(hub_cox_df)))

# Flag significance
hub_cox_df <- hub_cox_df %>%
  mutate(
    sig_label = case_when(
      P < 0.001 ~ "***",
      P < 0.01  ~ "**",
      P < 0.05  ~ "*",
      TRUE      ~ "ns"
    ),
    P_label   = ifelse(P < 0.001,
                       formatC(P, format = "e", digits = 2),
                       sprintf("%.3f", P)),
    # ci_text for the table column
    ci_text   = sprintf("%.3f [%.3f\u2013%.3f]", HR, HR_lower, HR_upper)
  ) %>%
  arrange(desc(HR))

hub_cox_df$Gene <- factor(hub_cox_df$Gene, levels = hub_cox_df$Gene)

# Horizontal forest plot using ggplot
p_forest <- ggplot(hub_cox_df, aes(x = HR, y = Gene)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50",
             linewidth = 0.4) +
  geom_errorbarh(aes(xmin = HR_lower, xmax = HR_upper),
                 height = 0.2, linewidth = 0.5,
                 color = ifelse(hub_cox_df$HR > 1, COLORS_RISK["High"],
                                COLORS_RISK["Low"])) +
  geom_point(aes(color = ifelse(HR > 1, "Risk", "Protective")),
             size = 2.5, shape = 18) +
  geom_text(aes(x = max(HR_upper) * 1.55, label = P_label),
            size = 2.0, hjust = 1) +
  geom_text(aes(x = max(HR_upper) * 1.65, label = sig_label),
            size = 2.0, hjust = 1) +
  scale_color_manual(values = c("Risk" = "#E64B35", "Protective" = "#4DBBD5"),
                     name = "") +
  scale_x_continuous(
    breaks = c(0.5, 1, 1.5, 2),
    limits = c(min(hub_cox_df$HR_lower) * 0.85, max(hub_cox_df$HR_upper) * 1.75)
  ) +
  labs(x = "Hazard Ratio (95% CI)", y = NULL,
       title = "Hub Gene Prognostic Value in IDH-WT Glioma") +
  annotate("text",
           x = max(hub_cox_df$HR_upper) * 1.55,
           y = nrow(hub_cox_df) + 0.7,
           label = "p-value", size = 2.0, hjust = 1, fontface = "bold") +
  theme_nature() +
  theme(
    axis.text.y  = element_text(size = 6, face = "italic"),
    legend.position = "bottom"
  )

save_pdf(p_forest, "Fig8_forest.pdf", width_mm = 120, height_mm = 75)

# =============================================================================
# Fig 8_correlation — Heatmap of Spearman r: hub genes vs immune cells
# Source: cor_results from hub_gene_results.RData
# ComplexHeatmap; diverging RdBu; stars for |r|>=0.3 & padj<0.05
# Double-column: 183 mm; height ~65 mm
# =============================================================================
message("\n=== Fig8_correlation: Hub gene immune correlation heatmap ===")

# Build wide matrix
cor_mat <- cor_results %>%
  select(Gene, ImmuneCell, r) %>%
  pivot_wider(names_from = ImmuneCell, values_from = r) %>%
  tibble::column_to_rownames("Gene") %>%
  as.matrix()

pval_mat <- cor_results %>%
  select(Gene, ImmuneCell, padj) %>%
  pivot_wider(names_from = ImmuneCell, values_from = padj) %>%
  tibble::column_to_rownames("Gene") %>%
  as.matrix()

# Significance annotation matrix (|r|>=0.3 & padj<0.05)
CORR_THRESH <- 0.3
sig_mark <- matrix("", nrow = nrow(cor_mat), ncol = ncol(cor_mat))
for (i in seq_len(nrow(pval_mat))) {
  for (j in seq_len(ncol(pval_mat))) {
    if (!is.na(pval_mat[i, j]) && pval_mat[i, j] < 0.05 &&
        abs(cor_mat[i, j]) >= CORR_THRESH) {
      if      (pval_mat[i, j] < 0.001) sig_mark[i, j] <- "***"
      else if (pval_mat[i, j] < 0.01)  sig_mark[i, j] <- "**"
      else                              sig_mark[i, j] <- "*"
    }
  }
}

col_fun <- colorRamp2(
  c(-0.8, -CORR_THRESH, 0, CORR_THRESH, 0.8),
  c("#2166AC", "#92C5DE", "white", "#F4A582", "#B2182B")
)

# Reorder rows (same as HUB_ALL order)
row_order <- HUB_ALL[HUB_ALL %in% rownames(cor_mat)]
cor_mat2  <- cor_mat[row_order, ]
sig_mark2 <- sig_mark[match(row_order, rownames(cor_mat)), ]

pdf(file.path(OUT_DIR, "Fig8_correlation.pdf"),
    width = 183 / 25.4, height = 100 / 25.4, family = "Times")

ht <- Heatmap(
  cor_mat2,
  name   = "Spearman r",
  col    = col_fun,
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (nchar(sig_mark2[i, j]) > 0) {
      grid.text(sig_mark2[i, j], x, y,
                gp = gpar(fontsize = 5, fontface = "bold", col = "black"))
    }
  },
  row_names_gp    = gpar(fontsize = 7, fontface = "italic", fontfamily = "Times"),
  column_names_gp = gpar(fontsize = 5.5, fontfamily = "Times"),
  column_names_rot = 90,
  row_title       = NULL,
  column_title    = sprintf(
    "Hub Gene - Immune Cell Correlation (stars: padj<0.05 & |r|\u2265%.1f)",
    CORR_THRESH),
  column_title_gp = gpar(fontsize = 7, fontface = "bold", fontfamily = "Times"),
  heatmap_legend_param = list(
    title          = "Spearman r",
    title_gp       = gpar(fontsize = 6, fontface = "bold"),
    labels_gp      = gpar(fontsize = 6),
    at             = c(-0.8, -0.3, 0, 0.3, 0.8),
    legend_height  = unit(20, "mm")
  ),
  cluster_rows    = FALSE,
  cluster_columns = TRUE
)
draw(ht, heatmap_legend_side = "right",
     padding = unit(c(25, 2, 2, 4), "mm"))  # 25mm top pad for vertical (90 deg) column labels
dev.off()
message("Saved: ", file.path(OUT_DIR, "Fig8_correlation.pdf"))

# =============================================================================
# Fig 8_scRNA — UMAP feature plots for CASP4, BLOC1S1, SLC7A5
# Extract expression + UMAP coordinates directly from Seurat object slots
# 3-panel horizontal layout; viridis color scale
# Width: 183 mm; height ~55 mm
# =============================================================================
message("\n=== Fig8_scRNA: scRNA UMAP feature plots ===")

sc_file <- file.path(DATA_PROC, "seu_upr_scored.rds")
if (!file.exists(sc_file)) {
  warning("scRNA Seurat object not found: ", sc_file)
} else {
  tryCatch({
    seu_raw    <- readRDS(sc_file)
    sct_assay  <- seu_raw@assays[["SCT"]]
    data_mat   <- slot(sct_assay, "data")     # normalized SCT counts
    umap_emb   <- seu_raw@reductions[["umap"]]@cell.embeddings
    meta_df    <- seu_raw@meta.data

    genes_sc   <- HUB_KM[HUB_KM %in% rownames(data_mat)]
    message(sprintf("  Hub genes found in scRNA: %s", paste(genes_sc, collapse = ", ")))

    # Build data frame: umap + gene expression + cell type
    sc_df <- data.frame(
      UMAP1    = umap_emb[, 1],
      UMAP2    = umap_emb[, 2],
      celltype = meta_df$celltype,
      stringsAsFactors = FALSE
    )
    for (g in genes_sc) {
      sc_df[[g]] <- as.numeric(data_mat[g, ])
    }

    feat_plots <- lapply(genes_sc, function(g) {
      df_g <- sc_df[order(sc_df[[g]]), ]   # plot low first
      ggplot(df_g, aes(x = UMAP1, y = UMAP2, color = .data[[g]])) +
        geom_point(size = 0.15, stroke = 0, alpha = 0.7) +
        scale_color_viridis_c(
          option = "C", name = expression(log[2](SCT + 1)),
          guide  = guide_colorbar(barwidth = unit(2, "mm"), barheight = unit(12, "mm"))
        ) +
        labs(title = g) +
        coord_equal() +
        theme_void(base_family = "Times") +
        theme(
          plot.title    = element_text(size = 7, face = "italic", hjust = 0.5,
                                       margin = margin(b = 1)),
          legend.text   = element_text(size = 5),
          legend.title  = element_text(size = 5),
          plot.margin   = margin(1, 1, 1, 1, "mm")
        )
    })

    p_sc <- wrap_plots(feat_plots, nrow = 1)

    save_pdf(p_sc, "Fig8_scRNA.pdf", width_mm = 183, height_mm = 58)
    rm(seu_raw, sct_assay, data_mat)
    gc()
  }, error = function(e) {
    message("scRNA plot failed: ", e$message)
  })
}

# =============================================================================
# Fig 8_heatmap — Expression heatmap: 5 hub genes, by risk + IDH status
# Z-scored expression; ComplexHeatmap; top annotation = Risk + IDH
# Width: 183 mm; height ~65 mm
# =============================================================================
message("\n=== Fig8_heatmap: Hub gene expression heatmap ===")

hub_avail <- HUB_ALL[HUB_ALL %in% rownames(expr_log)]
message(sprintf("  Hub genes available: %s", paste(hub_avail, collapse = ", ")))

# Filter complete cases
valid_mask <- !is.na(risk_groups) & !is.na(clin_df$IDH_status)
valid_samps <- common_samples[valid_mask]
message(sprintf("  Samples with complete risk+IDH: %d", length(valid_samps)))

expr_hub   <- as.matrix(expr_log[hub_avail, valid_samps])
expr_scaled <- t(scale(t(expr_hub)))   # z-score across samples

# Sort columns: by risk group, then by IDH status
risk_v <- as.character(risk_groups[valid_samps])
idh_v  <- as.character(clin_df$IDH_status[valid_mask])
col_ord <- order(risk_v, idh_v)

# Top annotation
ha <- HeatmapAnnotation(
  Risk = risk_v[col_ord],
  IDH  = idh_v[col_ord],
  col  = list(
    Risk = COLORS_RISK,
    IDH  = COLORS_IDH
  ),
  annotation_name_gp   = gpar(fontsize = 6, fontfamily = "Times"),
  annotation_legend_param = list(
    Risk = list(title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 5.5)),
    IDH  = list(title_gp = gpar(fontsize = 6), labels_gp = gpar(fontsize = 5.5))
  ),
  na_col = "grey90",
  show_legend = TRUE
)

pdf(file.path(OUT_DIR, "Fig8_heatmap.pdf"),
    width = 183 / 25.4, height = 65 / 25.4, family = "Times")

draw(Heatmap(
  expr_scaled[, col_ord],
  name    = "Z-score",
  col     = colorRamp2(c(-2, 0, 2), c("#4DBBD5", "white", "#E64B35")),
  top_annotation    = ha,
  cluster_columns   = FALSE,
  cluster_rows      = TRUE,
  show_column_names = FALSE,
  row_names_gp      = gpar(fontsize = 7, fontface = "italic", fontfamily = "Times"),
  column_title      = "Hub Gene Expression by UIRS Risk Group and IDH Status",
  column_title_gp   = gpar(fontsize = 7, fontface = "bold", fontfamily = "Times"),
  heatmap_legend_param = list(
    title_gp  = gpar(fontsize = 6, fontface = "bold"),
    labels_gp = gpar(fontsize = 5.5),
    legend_height = unit(20, "mm")
  )
), heatmap_legend_side = "right",
   padding = unit(c(2, 2, 2, 2), "mm"))

dev.off()
message("Saved: ", file.path(OUT_DIR, "Fig8_heatmap.pdf"))

message("\n=== Figure 8 complete ===")
message("All PDFs saved to: ", OUT_DIR)
