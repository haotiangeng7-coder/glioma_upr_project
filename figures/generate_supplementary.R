###############################################################################
# generate_supplementary.R
# All Supplementary Figures for Journal of Translational Medicine submission
# FigS1_batch_correction, FigS2_cellcycle, FigS_ml_cindex, FigS_permutation,
# FigS_sensitivity_k2, FigS_sensitivity_hallmark_k2, FigS_immune_idhwt,
# FigS_dca_1yr, FigS_dca_5yr
# Output: manuscript/submission/Additional_file_1_Supplementary_Figures/
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(survival)
  library(survminer)
  library(viridis)
  library(RColorBrewer)
  library(showtext)
  library(sysfonts)
  library(grid)
  library(cowplot)
  library(dcurves)
  library(ComplexHeatmap)
  library(circlize)
})

# ── Fonts ─────────────────────────────────────────────────────────────────────
font_add("Helvetica",
         regular    = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
         bold       = "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
         italic     = "/usr/share/fonts/truetype/dejavu/DejaVuSans-Oblique.ttf",
         bolditalic = "/usr/share/fonts/truetype/dejavu/DejaVuSans-BoldOblique.ttf")
showtext_auto()
showtext_opts(dpi = 600)

# ── Paths ──────────────────────────────────────────────────────────────────────
PROJECT_DIR <- getwd()
DATA_PROC   <- file.path(PROJECT_DIR, "data", "processed")
RES_DIR     <- file.path(PROJECT_DIR, "results")
OUT_DIR     <- file.path(PROJECT_DIR, "manuscript", "submission",
                         "Additional_file_1_Supplementary_Figures")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Supplementary theme (8pt base, slightly smaller than main) ────────────────
theme_supp <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      axis.text        = element_text(size = 7, color = "black"),
      axis.title       = element_text(size = 8, color = "black"),
      plot.title       = element_text(size = 8, face = "bold", hjust = 0.5),
      legend.text      = element_text(size = 7),
      legend.title     = element_text(size = 7, face = "bold"),
      legend.key.size  = unit(3.5, "mm"),
      strip.text       = element_text(size = 7, face = "bold"),
      strip.background = element_blank(),
      panel.grid       = element_blank(),
      plot.margin      = margin(3, 3, 3, 3, "mm")
    )
}

COLORS_RISK <- c("High" = "#E64B35", "Low" = "#4DBBD5")
COLORS_SUBTYPE <- c(
  "UPR-high-risk"    = "#E64B35",
  "UPR-intermediate" = "#00A087",
  "UPR-favorable"    = "#4DBBD5"
)

save_pdf <- function(plot, filename, width_mm, height_mm) {
  path <- file.path(OUT_DIR, filename)
  ggsave(path, plot = plot,
         device = cairo_pdf, width = width_mm, height = height_mm,
         units = "mm", dpi = 600)
  message("Saved: ", path)
}

# =============================================================================
# FigS1_batch_correction — 4-panel: UMAP before/after Harmony +
#                           LISI distribution before/after
# Source: LISI_before / LISI_after in Seurat metadata; UMAP reduction
# Width: 183 mm; height ~140 mm (2×2 grid)
# =============================================================================
message("\n=== FigS1: Batch correction ===")

sc_file <- file.path(DATA_PROC, "seu_upr_scored.rds")
if (!file.exists(sc_file)) {
  warning("Seurat object not found for FigS1: ", sc_file)
} else {
  tryCatch({
    seu_raw  <- readRDS(sc_file)
    meta     <- seu_raw@meta.data
    umap_emb <- seu_raw@reductions[["umap"]]@cell.embeddings

    umap_df <- data.frame(
      UMAP1        = umap_emb[, 1],
      UMAP2        = umap_emb[, 2],
      platform     = meta$platform,
      celltype     = meta$celltype,
      LISI_before  = meta$LISI_before,
      LISI_after   = meta$LISI_after,
      stringsAsFactors = FALSE
    )

    message(sprintf("  Cells: %d;  platforms: %s",
                    nrow(umap_df), paste(unique(umap_df$platform), collapse = ", ")))

    # Subsample for speed (max 8000 points; keep proportional)
    set.seed(42)
    idx <- sample(nrow(umap_df), min(8000, nrow(umap_df)))
    umap_sub <- umap_df[idx, ]

    # Okabe-Ito for platform
    platform_cols <- c("10X" = "#0072B2", "Smart-seq2" = "#E69F00")

    p_before <- ggplot(umap_sub, aes(x = UMAP1, y = UMAP2, color = platform)) +
      geom_point(size = 0.25, alpha = 0.6, stroke = 0) +
      scale_color_manual(values = platform_cols, name = "Platform") +
      labs(title = "UMAP — Before Harmony") +
      coord_equal() +
      theme_void(base_family = "Helvetica") +
      theme(
        plot.title   = element_text(size = 8, face = "bold", hjust = 0.5),
        legend.text  = element_text(size = 7),
        legend.title = element_text(size = 7, face = "bold"),
        plot.margin  = margin(2, 2, 2, 2, "mm")
      ) +
      guides(color = guide_legend(override.aes = list(size = 2)))

    p_after <- ggplot(umap_sub, aes(x = UMAP1, y = UMAP2, color = celltype)) +
      geom_point(size = 0.25, alpha = 0.6, stroke = 0) +
      scale_color_brewer(palette = "Paired", name = "Cell type") +
      labs(title = "UMAP — After Harmony") +
      coord_equal() +
      theme_void(base_family = "Helvetica") +
      theme(
        plot.title   = element_text(size = 8, face = "bold", hjust = 0.5),
        legend.text  = element_text(size = 7),
        legend.title = element_text(size = 7, face = "bold"),
        plot.margin  = margin(2, 2, 2, 2, "mm")
      ) +
      guides(color = guide_legend(override.aes = list(size = 2), ncol = 1))

    # LISI distributions
    lisi_long <- umap_df %>%
      select(LISI_before, LISI_after, platform) %>%
      pivot_longer(c(LISI_before, LISI_after),
                   names_to = "Stage", values_to = "LISI") %>%
      mutate(Stage = factor(Stage,
                            levels = c("LISI_before", "LISI_after"),
                            labels = c("Before Harmony", "After Harmony")))

    lisi_stats <- lisi_long %>%
      group_by(Stage) %>%
      summarise(med = median(LISI, na.rm = TRUE), .groups = "drop")
    message(sprintf("  LISI median: before=%.2f, after=%.2f",
                    lisi_stats$med[1], lisi_stats$med[2]))

    p_lisi_density <- ggplot(lisi_long, aes(x = LISI, fill = Stage)) +
      geom_density(alpha = 0.55, color = NA) +
      geom_vline(
        data = lisi_stats,
        aes(xintercept = med, color = Stage),
        linetype = "dashed", linewidth = 0.5
      ) +
      scale_fill_manual(values  = c("Before Harmony" = "#CCBBAA",
                                    "After Harmony"  = "#4DBBD5"),
                        name = "Stage") +
      scale_color_manual(values = c("Before Harmony" = "#CCBBAA",
                                    "After Harmony"  = "#4DBBD5"),
                         guide = "none") +
      labs(x = "LISI Score", y = "Density",
           title = "LISI Distribution Before/After Harmony") +
      theme_supp() +
      theme(legend.position = "bottom")

    # Box for per-platform comparison
    p_lisi_box <- ggplot(lisi_long, aes(x = Stage, y = LISI, fill = Stage)) +
      geom_boxplot(width = 0.5, outlier.size = 0.3, outlier.alpha = 0.3) +
      stat_compare_means(method = "wilcox.test", label = "p.format",
                         size = 2.5, label.x.npc = 0.5) +
      scale_fill_manual(values = c("Before Harmony" = "#CCBBAA",
                                   "After Harmony"  = "#4DBBD5")) +
      labs(x = NULL, y = "LISI Score",
           title = "LISI Score Comparison") +
      theme_supp() +
      theme(legend.position = "none",
            axis.text.x = element_text(size = 7, angle = 15, hjust = 1))

    p_s1 <- (p_before | p_after) / (p_lisi_density | p_lisi_box) +
      plot_annotation(
        tag_levels = "a",
        theme = theme(plot.tag = element_text(size = 8, face = "bold",
                                              family = "Helvetica"))
      )

    save_pdf(p_s1, "FigS1_batch_correction.pdf", width_mm = 183, height_mm = 140)

    rm(seu_raw)
    gc()
  }, error = function(e) {
    message("FigS1 failed: ", e$message)
  })
}

# =============================================================================
# FigS2_cellcycle — Cell-cycle score (S.Score) vs UPR score per cell type
# Pearson r and p-value annotated; scatter with regression line
# Source: Seurat metadata columns S.Score and UPR_score
# Width: 183 mm; height ~120 mm (multi-facet)
# =============================================================================
message("\n=== FigS2: Cell cycle vs UPR score ===")

if (!file.exists(sc_file)) {
  warning("Seurat object not found for FigS2")
} else {
  tryCatch({
    seu_raw <- readRDS(sc_file)
    meta2   <- seu_raw@meta.data

    cc_df <- data.frame(
      S_Score  = meta2$S.Score,
      UPR_score = meta2$UPR_score,
      celltype  = meta2$celltype,
      stringsAsFactors = FALSE
    ) %>%
      filter(!is.na(S_Score) & !is.na(UPR_score) & !is.na(celltype))

    message(sprintf("  Cells for FigS2: %d", nrow(cc_df)))
    message("  Global correlation check...")

    # Per cell-type Pearson r
    cc_stats <- cc_df %>%
      group_by(celltype) %>%
      summarise(
        r    = cor(S_Score, UPR_score, method = "pearson"),
        pval = cor.test(S_Score, UPR_score, method = "pearson")$p.value,
        n    = n(),
        .groups = "drop"
      ) %>%
      mutate(
        padj  = p.adjust(pval, method = "BH"),
        label = sprintf("r=%.2f\np=%.2e", r, pval)
      )

    message("  Cell-type correlations:")
    print(cc_stats)

    # Subsample for plot (max 1000 per cell type)
    set.seed(42)
    cc_sub <- cc_df %>%
      group_by(celltype) %>%
      do({
        n_take <- min(nrow(.), 1000)
        .[sample(nrow(.), n_take), ]
      }) %>%
      ungroup()

    p_s2 <- ggplot(cc_sub, aes(x = S_Score, y = UPR_score)) +
      geom_point(size = 0.3, alpha = 0.35, color = "#4DBBD5") +
      geom_smooth(method = "lm", formula = y ~ x,
                  color = "#E64B35", linewidth = 0.5, se = TRUE, alpha = 0.2) +
      geom_text(
        data = cc_stats,
        aes(x = -Inf, y = Inf, label = label),
        hjust = -0.1, vjust = 1.3, size = 2.2, inherit.aes = FALSE
      ) +
      facet_wrap(~celltype, scales = "free", ncol = 3) +
      labs(
        x = "Cell Cycle S-phase Score",
        y = "UPR Score (UCell)",
        title = "Cell Cycle vs UPR Score per Cell Type"
      ) +
      theme_supp() +
      theme(strip.text = element_text(size = 7))

    save_pdf(p_s2, "FigS2_cellcycle_upr.pdf", width_mm = 183, height_mm = 130)

    rm(seu_raw)
    gc()
  }, error = function(e) {
    message("FigS2 failed: ", e$message)
  })
}

# =============================================================================
# FigS_ml_cindex — C-index distribution across 100 ML combinations
# Source: ml_results data.frame in ml_combination_results.RData
# Violin + box sorted by median nested-CV C-index
# Width: 183 mm; height ~90 mm
# =============================================================================
message("\n=== FigS_ml_cindex: ML C-index distribution ===")

load(file.path(DATA_PROC, "ml_combination_results.RData"))
message(sprintf("  ml_results: %d combinations", nrow(ml_results)))

# Extract per-combination nested CV C-indices from all_combo_results
nested_df <- lapply(all_combo_results, function(x) {
  data.frame(
    combo_name    = x$combo_name,
    nested_cindex = x$nested_cv_cindex,
    stringsAsFactors = FALSE
  )
})
nested_df <- bind_rows(nested_df)

# Sort by median (use ml_results ordering for consistency)
combo_order <- ml_results %>%
  arrange(desc(nested_cv_mean)) %>%
  pull(combo_name)
nested_df$combo_name <- factor(nested_df$combo_name, levels = rev(combo_order))

# Horizontal violin sorted by median
p_ml <- ggplot(nested_df, aes(x = combo_name, y = nested_cindex)) +
  geom_violin(fill = "#4DBBD5", alpha = 0.4, color = NA, trim = TRUE) +
  geom_boxplot(width = 0.25, outlier.size = 0.3, fill = "white",
               color = "#2166AC") +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey60",
             linewidth = 0.35) +
  coord_flip() +
  labs(
    x = NULL,
    y = "Harrell's C-index (nested 5-fold CV)",
    title = "C-index Distribution Across 100 ML Algorithm Combinations"
  ) +
  theme_supp() +
  theme(
    axis.text.y  = element_text(size = 5.5),
    axis.text.x  = element_text(size = 7)
  )

save_pdf(p_ml, "FigS_ml_cindex.pdf", width_mm = 183, height_mm = 180)

# =============================================================================
# FigS_permutation — Permutation test: histogram of permuted C-indices
# + observed C-index as red vertical line; p-value annotated
# Source: coxph_permutation_results.RData
# Width: 89 mm; height ~70 mm
# =============================================================================
message("\n=== FigS_permutation: Permutation test ===")

load(file.path(DATA_PROC, "coxph_permutation_results.RData"))
message(sprintf("  Observed C-index: %.4f, permutation p = %.3f, n_perm = %d",
                obs_cindex, perm_p, length(perm_cindex)))

perm_df <- data.frame(cindex = perm_cindex)

p_perm <- ggplot(perm_df, aes(x = cindex)) +
  geom_histogram(binwidth = 0.005, fill = "#92C5DE", color = "white",
                 linewidth = 0.2) +
  geom_vline(xintercept = obs_cindex, color = "#E64B35", linewidth = 0.7,
             linetype = "solid") +
  annotate("text",
           x = obs_cindex + 0.002,
           y = Inf,
           label = sprintf("Observed\nC-index = %.4f\np = %.3f",
                           obs_cindex, perm_p),
           vjust = 1.3, hjust = 0, size = 2.5, color = "#E64B35") +
  labs(
    x = "Permuted C-index",
    y = "Count",
    title = sprintf("Permutation Test (n = %d)", length(perm_cindex))
  ) +
  theme_supp()

save_pdf(p_perm, "FigS_permutation.pdf", width_mm = 89, height_mm = 70)

# =============================================================================
# FigS_sensitivity_k2 — KM using all UPR genes at K=2
# Source: cc_results_all from consensus_clustering_results.RData
# Width: 89 mm; height ~90 mm
# =============================================================================
message("\n=== FigS_sensitivity_k2: All-UPR-genes K=2 KM ===")

load(file.path(DATA_PROC, "consensus_clustering_results.RData"))

# cc_results_all[[2]] has consensusClass for K=2
cluster_k2_all <- cc_results_all[[2]]$consensusClass
# Match to surv_final
names_k2 <- names(cluster_k2_all)
surv_k2   <- surv_final %>%
  filter(barcode %in% names_k2) %>%
  mutate(
    cluster_k2 = factor(cluster_k2_all[barcode],
                        levels = c(1, 2),
                        labels = c("Cluster 1", "Cluster 2")),
    time_months = OS.time / 30
  ) %>%
  filter(!is.na(cluster_k2) & !is.na(OS.time) & OS.time > 0)

message(sprintf("  All-gene K=2: %d samples (%s)",
                nrow(surv_k2), paste(table(surv_k2$cluster_k2), collapse = " / ")))

fit_k2_all <- survfit(Surv(time_months, OS) ~ cluster_k2, data = surv_k2)
lr_k2      <- survdiff(Surv(time_months, OS) ~ cluster_k2, data = surv_k2)
p_lr_k2    <- 1 - pchisq(lr_k2$chisq, df = 1)

km_k2_all <- ggsurvplot(
  fit_k2_all, data = surv_k2,
  palette      = c("#4DBBD5", "#E64B35"),
  risk.table   = TRUE,
  pval         = sprintf("p = %.2e", p_lr_k2),
  pval.method  = FALSE,
  conf.int     = TRUE, conf.int.alpha = 0.15,
  legend.title = "Cluster",
  xlab         = "Time (months)",
  ylab         = "Overall Survival",
  title        = "Sensitivity Analysis: All UPR Genes, K=2",
  ggtheme      = theme_supp(),
  fontsize     = 2.5,
  risk.table.fontsize = 2.5,
  risk.table.height   = 0.28,
  tables.theme = theme_cleantable() +
    theme(text = element_text(size = 7, family = "Helvetica"))
)

pdf(file.path(OUT_DIR, "FigS_sensitivity_k2.pdf"),
    width = 89 / 25.4, height = 90 / 25.4, family = "Helvetica")
pg_k2 <- plot_grid(
  km_k2_all$plot  + theme(plot.margin = margin(1, 1, 0, 1, "mm")),
  km_k2_all$table + theme(plot.margin = margin(0, 1, 1, 1, "mm")),
  ncol = 1, rel_heights = c(0.72, 0.28)
)
print(pg_k2)
dev.off()
message("Saved: ", file.path(OUT_DIR, "FigS_sensitivity_k2.pdf"))

# =============================================================================
# FigS_sensitivity_hallmark_k2 — KM using MSigDB Hallmark UPR genes at K=2
# Source: cc_results_hallmark[[2]] from consensus_clustering_results.RData
# Width: 89 mm; height ~90 mm
# =============================================================================
message("\n=== FigS_sensitivity_hallmark_k2: Hallmark UPR K=2 KM ===")

cluster_k2_hall <- cc_results_hallmark[[2]]$consensusClass
names_k2h        <- names(cluster_k2_hall)
surv_k2h <- surv_final %>%
  filter(barcode %in% names_k2h) %>%
  mutate(
    cluster_k2 = factor(cluster_k2_hall[barcode],
                        levels = c(1, 2),
                        labels = c("Cluster 1", "Cluster 2")),
    time_months = OS.time / 30
  ) %>%
  filter(!is.na(cluster_k2) & !is.na(OS.time) & OS.time > 0)

message(sprintf("  Hallmark K=2: %d samples (%s)",
                nrow(surv_k2h), paste(table(surv_k2h$cluster_k2), collapse = " / ")))

fit_k2h <- survfit(Surv(time_months, OS) ~ cluster_k2, data = surv_k2h)
lr_k2h  <- survdiff(Surv(time_months, OS) ~ cluster_k2, data = surv_k2h)
p_lr_k2h <- 1 - pchisq(lr_k2h$chisq, df = 1)

km_k2h <- ggsurvplot(
  fit_k2h, data = surv_k2h,
  palette      = c("#4DBBD5", "#E64B35"),
  risk.table   = TRUE,
  pval         = sprintf("p = %.2e", p_lr_k2h),
  pval.method  = FALSE,
  conf.int     = TRUE, conf.int.alpha = 0.15,
  legend.title = "Cluster",
  xlab         = "Time (months)",
  ylab         = "Overall Survival",
  title        = "Sensitivity Analysis: Hallmark UPR, K=2",
  ggtheme      = theme_supp(),
  fontsize     = 2.5,
  risk.table.fontsize = 2.5,
  risk.table.height   = 0.28,
  tables.theme = theme_cleantable() +
    theme(text = element_text(size = 7, family = "Helvetica"))
)

pdf(file.path(OUT_DIR, "FigS_sensitivity_hallmark_k2.pdf"),
    width = 89 / 25.4, height = 90 / 25.4, family = "Helvetica")
pg_k2h <- plot_grid(
  km_k2h$plot  + theme(plot.margin = margin(1, 1, 0, 1, "mm")),
  km_k2h$table + theme(plot.margin = margin(0, 1, 1, 1, "mm")),
  ncol = 1, rel_heights = c(0.72, 0.28)
)
print(pg_k2h)
dev.off()
message("Saved: ", file.path(OUT_DIR, "FigS_sensitivity_hallmark_k2.pdf"))

# =============================================================================
# FigS_immune_idhwt — ssGSEA immune heatmap restricted to IDH-WT patients
# Source: immune_ssgsea_kw_idhwt.csv (CellType, kw_pvalue, kw_padj)
#         + immune_analysis_results for the score matrix
# Width: 183 mm; height ~80 mm
# =============================================================================
message("\n=== FigS_immune_idhwt: IDH-WT immune ssGSEA heatmap ===")

immune_kw_file <- file.path(RES_DIR, "immune_ssgsea_kw_idhwt.csv")
message("  immune_ssgsea_kw_idhwt.csv exists: ", file.exists(immune_kw_file))

# Load immune analysis results for the ssGSEA score matrix
immune_rdata <- file.path(DATA_PROC, "immune_analysis_results.RData")
if (file.exists(immune_rdata)) {
  load(immune_rdata)
  message("  immune_analysis_results.RData loaded")
  message("  Objects: ", paste(ls(), collapse = ", "))
}

# Load expression + risk model for IDH-WT sample list
load(file.path(DATA_PROC, "risk_model_final.RData"))
load(file.path(DATA_PROC, "tcga_glioma_expression.RData"))

# IDH-WT samples
common_samples <- intersect(names(risk_score_train), colnames(expr_tpm_symbol))
risk_groups    <- factor(risk_group_train[common_samples], levels = c("Low", "High"))
clin_c         <- clinical_valid[match(common_samples, clinical_valid$barcode), ]
idhwt_mask     <- !is.na(clin_c$IDH_status) & clin_c$IDH_status == "WT"
idhwt_samps    <- common_samples[idhwt_mask]

if (exists("immune_ssgsea")) {
  # Match IDH-WT samples
  common_imm <- intersect(idhwt_samps, colnames(immune_ssgsea))
  message(sprintf("  IDH-WT samples with ssGSEA: %d", length(common_imm)))

  imm_mat <- immune_ssgsea[, common_imm]
  # Scale per row (cell type)
  imm_scaled <- t(scale(t(as.matrix(imm_mat))))

  # Top annotation: risk group
  risk_idhwt <- as.character(risk_groups[common_imm])
  # Sort by risk group
  col_ord <- order(risk_idhwt)

  ha <- HeatmapAnnotation(
    Risk = risk_idhwt[col_ord],
    col  = list(Risk = c("High" = "#E64B35", "Low" = "#4DBBD5")),
    annotation_name_gp   = gpar(fontsize = 7, fontfamily = "Helvetica"),
    annotation_legend_param = list(
      Risk = list(title_gp = gpar(fontsize = 7), labels_gp = gpar(fontsize = 6.5))
    )
  )

  # Add KW significance from file if available
  if (file.exists(immune_kw_file)) {
    kw_df <- read.csv(immune_kw_file, stringsAsFactors = FALSE)
    sig_cells <- kw_df %>%
      filter(kw_padj < 0.05) %>%
      pull(CellType)
    message(sprintf("  Significant immune cells (IDH-WT): %d", length(sig_cells)))

    # Row annotation for significance
    ra_df <- data.frame(
      CellType  = rownames(imm_scaled),
      Significant = ifelse(rownames(imm_scaled) %in% sig_cells, "padj<0.05", "ns"),
      stringsAsFactors = FALSE
    )
    ra <- rowAnnotation(
      KW_sig = ra_df$Significant,
      col    = list(KW_sig = c("padj<0.05" = "#E64B35", "ns" = "grey80")),
      annotation_name_gp = gpar(fontsize = 6, fontfamily = "Helvetica"),
      show_legend = TRUE
    )
  } else {
    ra <- NULL
  }

  pdf(file.path(OUT_DIR, "FigS_immune_idhwt.pdf"),
      width = 183 / 25.4, height = 85 / 25.4, family = "Helvetica")

  ht_imm <- Heatmap(
    imm_scaled[, col_ord],
    name   = "Z-score",
    col    = colorRamp2(c(-2, 0, 2), c("#4DBBD5", "white", "#E64B35")),
    top_annotation    = ha,
    right_annotation  = ra,
    cluster_columns   = FALSE,
    cluster_rows      = TRUE,
    show_column_names = FALSE,
    row_names_gp      = gpar(fontsize = 5.5, fontfamily = "Helvetica"),
    column_title      = "Immune Cell ssGSEA Enrichment Scores — IDH-WT Glioma",
    column_title_gp   = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica"),
    heatmap_legend_param = list(
      title_gp      = gpar(fontsize = 7, fontface = "bold"),
      labels_gp     = gpar(fontsize = 6),
      legend_height = unit(20, "mm")
    )
  )
  draw(ht_imm, heatmap_legend_side = "right",
       padding = unit(c(2, 2, 2, 2), "mm"))
  dev.off()
  message("Saved: ", file.path(OUT_DIR, "FigS_immune_idhwt.pdf"))

} else {
  message("  immune_ssgsea object not found — recomputing from expression data...")
  # Compute ssGSEA on the fly using GSVA
  library(GSVA)
  load(file.path(DATA_PROC, "upr_gene_sets.RData"))

  immune_cell_markers <- list(
    "Activated CD4 T" = c("DPP4","ICOS","IL2RA","IL4R","OAS1","UBE2L6"),
    "Activated CD8 T" = c("CD8A","CD8B","EOMES","GZMA","GZMB","IFNG","PRF1"),
    "Regulatory T"    = c("CTLA4","FOXP3","IL2RA","IKZF2"),
    "NK cell"         = c("FCER1G","GNLY","KLRB1","KLRD1","KLRF1","NKG7"),
    "Macrophage"      = c("CD14","CD68","CSF1R","ITGAM","MSR1"),
    "MDSC"            = c("ARG1","CD33","ITGAM","S100A8","S100A9"),
    "Monocyte"        = c("CD14","CSF1R","FCGR1A","LYZ","VCAN"),
    "M2 Macrophage"   = c("CD163","MRC1","MSR1","CD68","IL10"),
    "Activated DC"    = c("BATF3","CD1C","CLEC10A","FCER1A","HLA-DQA1"),
    "T exhaust"       = c("PDCD1","HAVCR2","LAG3","TIGIT","CTLA4","TOX"),
    "IFN Response"    = c("IFI35","IFI44","IFIT1","IFIT3","MX1","OAS1")
  )
  expr_log <- log2(expr_tpm_symbol[, idhwt_samps] + 1)
  valid_gs  <- lapply(immune_cell_markers, function(g) g[g %in% rownames(expr_log)])
  valid_gs  <- valid_gs[sapply(valid_gs, length) >= 3]

  ss_params  <- ssgseaParam(as.matrix(expr_log), valid_gs)
  imm_mat_c  <- gsva(ss_params)

  # Sort by risk
  risk_idhwt <- as.character(risk_groups[idhwt_samps])
  col_ord2   <- order(risk_idhwt)
  imm_sc2    <- t(scale(t(imm_mat_c)))

  ha2 <- HeatmapAnnotation(
    Risk = risk_idhwt[col_ord2],
    col  = list(Risk = c("High" = "#E64B35", "Low" = "#4DBBD5")),
    annotation_name_gp = gpar(fontsize = 7, fontfamily = "Helvetica")
  )

  pdf(file.path(OUT_DIR, "FigS_immune_idhwt.pdf"),
      width = 183 / 25.4, height = 80 / 25.4, family = "Helvetica")
  draw(Heatmap(
    imm_sc2[, col_ord2],
    name = "Z-score",
    col  = colorRamp2(c(-2, 0, 2), c("#4DBBD5", "white", "#E64B35")),
    top_annotation = ha2,
    cluster_columns = FALSE, cluster_rows = TRUE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 6, fontfamily = "Helvetica"),
    column_title = "Immune Cell ssGSEA — IDH-WT Glioma",
    column_title_gp = gpar(fontsize = 7, fontface = "bold", fontfamily = "Helvetica")
  ), heatmap_legend_side = "right")
  dev.off()
  message("Saved: ", file.path(OUT_DIR, "FigS_immune_idhwt.pdf"))
}

# =============================================================================
# FigS_dca_1yr and FigS_dca_5yr — Decision Curve Analysis at 1-year and 5-year
# Source: nomo_df from nomogram_results.RData
# Use dcurves package
# Width: 89 mm each; height ~80 mm
# =============================================================================
message("\n=== FigS_dca: Decision Curve Analysis ===")

load(file.path(DATA_PROC, "nomogram_results.RData"))
message(sprintf("  nomo_df: %d rows, cols: %s",
                nrow(nomo_df), paste(colnames(nomo_df), collapse = ", ")))

# nomo_df has OS_months, OS, risk_score, Age, Grade, IDH, nomogram_lp
# OS_months is already in months
# Prepare survival object
nomo_surv <- nomo_df %>%
  filter(!is.na(OS_months) & OS_months > 0 & !is.na(OS)) %>%
  mutate(
    risk_score = as.numeric(risk_score),
    Age        = as.numeric(Age)
  )
message(sprintf("  DCA samples: %d", nrow(nomo_surv)))

# Standardize risk_score to 0-1 range for DCA (required by dcurves)
rs_range     <- range(nomo_surv$risk_score, na.rm = TRUE)
nomo_surv$risk_prob <- (nomo_surv$risk_score - rs_range[1]) /
  (rs_range[2] - rs_range[1])

# Compute predicted probability at t* using Cox fitted from nomo_df
surv_obj <- Surv(nomo_surv$OS_months, nomo_surv$OS)

# Fit individual and combined models
cox_risk <- coxph(surv_obj ~ risk_score,    data = nomo_surv, x = TRUE)
cox_age  <- coxph(surv_obj ~ Age,           data = nomo_surv, x = TRUE)
cox_full <- coxph(surv_obj ~ risk_score + Age + Grade + IDH,
                  data = nomo_surv, x = TRUE)

# Predict survival probability at 12 and 60 months
get_surv_prob <- function(cox_model, newdata, t_star) {
  sf <- survfit(cox_model, newdata = newdata)
  # Get probability at t_star (nearest time point)
  idx <- which.min(abs(summary(sf)$time - t_star))
  if (length(idx) == 0) return(rep(NA, nrow(newdata)))
  # survival probability = 1 - cumulative hazard at baseline ^ exp(LP)
  # Use survfit summary to extract
  sf_sum <- summary(sf, times = t_star)
  if (is.null(sf_sum$surv)) return(rep(NA, nrow(newdata)))
  sf_sum$surv
}

# For DCA we need event probability = 1 - surv_prob
make_event_prob <- function(cox_model, newdata, t_star) {
  bh <- basehaz(cox_model, centered = FALSE)
  h0 <- bh$hazard[which.min(abs(bh$time - t_star))]
  lp <- predict(cox_model, newdata = newdata, type = "lp")
  1 - exp(-h0 * exp(lp))
}

for (t_yr in c(12, 60)) {
  label_yr <- if (t_yr == 12) "1yr" else "5yr"
  message(sprintf("  Computing DCA at %d months (%s)", t_yr, label_yr))

  nomo_surv[[paste0("p_risk_", label_yr)]] <-
    make_event_prob(cox_risk, nomo_surv, t_yr)
  nomo_surv[[paste0("p_age_",  label_yr)]] <-
    make_event_prob(cox_age,  nomo_surv, t_yr)
  nomo_surv[[paste0("p_full_", label_yr)]] <-
    make_event_prob(cox_full, nomo_surv, t_yr)

  # Remove NA rows for this time point
  dca_df <- nomo_surv %>%
    filter(
      is.finite(.data[[paste0("p_risk_", label_yr)]]) &
      is.finite(.data[[paste0("p_age_",  label_yr)]]) &
      is.finite(.data[[paste0("p_full_", label_yr)]])
    )

  message(sprintf("  DCA samples at %s: %d", label_yr, nrow(dca_df)))

  tryCatch({
    # Build formula with the predicted probability variables
    f_dca <- as.formula(
      sprintf("Surv(OS_months, OS) ~ p_risk_%s + p_age_%s + p_full_%s",
              label_yr, label_yr, label_yr)
    )
    dca_res <- dca(f_dca, data = dca_df, time = t_yr,
                   thresholds = seq(0.05, 0.90, by = 0.02))

    # Build named color / label vectors outside c() to avoid parse error
    col_names  <- c("All", "None",
                    paste0("p_risk_", label_yr),
                    paste0("p_age_",  label_yr),
                    paste0("p_full_", label_yr))
    col_vals   <- c("grey70", "grey30", "#E64B35", "#4DBBD5", "#00A087")
    col_labels <- c("Treat all", "Treat none",
                    "UIRS risk score", "Age", "Combined model")
    names(col_vals)   <- col_names
    names(col_labels) <- col_names

    p_dca <- plot(dca_res, smooth = TRUE) +
      scale_color_manual(values = col_vals, labels = col_labels, name = "Model") +
      labs(
        x     = "Threshold Probability",
        y     = "Net Benefit",
        title = sprintf("Decision Curve Analysis — %s Survival",
                        if (t_yr == 12) "1-Year" else "5-Year")
      ) +
      theme_supp() +
      theme(legend.position = "right")

    save_pdf(p_dca, sprintf("FigS_dca_%s.pdf", label_yr),
             width_mm = 120, height_mm = 80)

  }, error = function(e) {
    message(sprintf("  DCA at %s failed: %s", label_yr, e$message))
    message("  Falling back to manual net benefit calculation...")

    # Manual calculation
    calc_net_benefit <- function(threshold, surv_data, t_star, prob_col) {
      p_event <- surv_data[[prob_col]]
      pred_pos <- p_event >= threshold
      n_total  <- nrow(surv_data)
      if (sum(pred_pos) == 0) return(NA_real_)

      # Kaplan-Meier event rate at t_star among predicted positives
      tryCatch({
        sf_pos <- survfit(Surv(OS_months, OS) ~ 1,
                          data = surv_data[pred_pos, ])
        sf_sum <- summary(sf_pos, times = t_star)
        if (is.null(sf_sum$surv) || length(sf_sum$surv) == 0)
          return(NA_real_)
        event_rate_pos <- 1 - sf_sum$surv[1]
        tp <- event_rate_pos * sum(pred_pos) / n_total
        fp <- (1 - event_rate_pos) * sum(pred_pos) / n_total
        nb <- tp - fp * threshold / (1 - threshold)
        nb
      }, error = function(e2) NA_real_)
    }

    # Overall event rate at t_star
    sf_all   <- survfit(Surv(OS_months, OS) ~ 1, data = nomo_surv)
    sf_t     <- summary(sf_all, times = t_yr)
    overall_event_rate <- if (!is.null(sf_t$surv)) 1 - sf_t$surv[1] else 0.5

    thresholds <- seq(0.05, 0.90, by = 0.02)
    p_col <- paste0("p_risk_", label_yr)

    nb_risk <- sapply(thresholds, function(thr)
      calc_net_benefit(thr, nomo_surv, t_yr, p_col))
    nb_all  <- sapply(thresholds, function(thr)
      overall_event_rate - (1 - overall_event_rate) * thr / (1 - thr))
    nb_none <- rep(0, length(thresholds))

    dca_manual <- data.frame(
      threshold = rep(thresholds, 3),
      net_benefit = c(nb_risk, nb_all, nb_none),
      model = rep(c("UIRS Risk Score", "Treat all", "Treat none"),
                  each = length(thresholds))
    ) %>% filter(!is.na(net_benefit))

    p_dca_manual <- ggplot(
      dca_manual, aes(x = threshold, y = net_benefit, color = model)
    ) +
      geom_line(linewidth = 0.7) +
      geom_hline(yintercept = 0, linetype = "solid", color = "grey60",
                 linewidth = 0.3) +
      scale_color_manual(
        values = c("UIRS Risk Score" = "#E64B35",
                   "Treat all"       = "grey50",
                   "Treat none"      = "grey20"),
        name = "Model"
      ) +
      coord_cartesian(ylim = c(-0.05, NA)) +
      labs(
        x = "Threshold Probability",
        y = "Net Benefit",
        title = sprintf("Decision Curve Analysis — %s Survival",
                        if (t_yr == 12) "1-Year" else "5-Year")
      ) +
      theme_supp() +
      theme(legend.position = "right")

    save_pdf(p_dca_manual, sprintf("FigS_dca_%s.pdf", label_yr),
             width_mm = 120, height_mm = 80)
  })
}

message("\n=== Supplementary Figures complete ===")
message("All PDFs saved to: ", OUT_DIR)
