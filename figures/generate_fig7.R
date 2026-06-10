###############################################################################
# generate_fig7.R
# Figure 7 — Immunotherapy and Drug Sensitivity
# Panels: Fig7A (checkpoint gene expression), Fig7B (IPS), Fig7D (drug pathway
#         bar chart), Fig7E (volcano DEG), Fig7E_GO (GO dotplot)
# Output: manuscript/submission/Figures/Figure7/
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggpubr)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(EnhancedVolcano)
  library(ggrepel)
  library(viridis)
  library(RColorBrewer)
  library(ggsci)
  library(showtext)
  library(sysfonts)
  library(grid)
  library(enrichplot)
  library(clusterProfiler)
})

# ── Fonts ─────────────────────────────────────────────────────────────────────
# DejaVu Serif is a metric-compatible TrueType serif font available on this
# system. Registered as "Times" to satisfy the journal's Times New Roman
# requirement.
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
OUT_DIR     <- file.path(PROJECT_DIR, "manuscript", "submission", "Figures", "Figure7")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Color scheme ───────────────────────────────────────────────────────────────
COLORS_RISK <- c("High" = "#E64B35", "Low" = "#4DBBD5")

# ── Nature-style theme ────────────────────────────────────────────────────────
theme_nature <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Times") +
    theme(
      axis.text        = element_text(size = 6,  color = "black"),
      axis.title       = element_text(size = 7,  color = "black"),
      plot.title       = element_text(size = 7,  face = "bold", hjust = 0.5),
      legend.text      = element_text(size = 6),
      legend.title     = element_text(size = 6,  face = "bold"),
      legend.key.size  = unit(3, "mm"),
      strip.text       = element_text(size = 6,  face = "bold"),
      strip.background = element_blank(),
      panel.grid       = element_blank(),
      plot.margin      = margin(2, 2, 2, 2, "mm")
    )
}

# ── Save helper ────────────────────────────────────────────────────────────────
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
load(file.path(DATA_PROC, "immunotherapy_prediction_results.RData"))
load(file.path(DATA_PROC, "pathway_analysis_results.RData"))
load(file.path(DATA_PROC, "upr_gene_sets.RData"))

# Build common data frame
common_samples <- intersect(names(risk_score_train), colnames(expr_tpm_symbol))
risk_groups    <- factor(risk_group_train[common_samples], levels = c("Low", "High"))
expr_log       <- log2(expr_tpm_symbol[, common_samples] + 1)

message(sprintf("Common samples: %d  (High=%d, Low=%d)",
                length(common_samples),
                sum(risk_groups == "High"),
                sum(risk_groups == "Low")))

# =============================================================================
# Fig 7A — Key immune checkpoint gene expression: High vs Low risk
# 6 genes: CD274, PDCD1, CTLA4, HAVCR2, LAG3, TIGIT
# Grouped violin + box, coloured by risk group
# Nature single-column width: 89 mm; height ~80 mm
# =============================================================================
message("\n=== Fig7A: Checkpoint gene expression ===")

genes_7a   <- c("CD274", "PDCD1", "CTLA4", "HAVCR2", "LAG3", "TIGIT")
genes_avail <- genes_7a[genes_7a %in% rownames(expr_log)]
missing_7a  <- setdiff(genes_7a, genes_avail)
if (length(missing_7a) > 0)
  warning("Fig7A: genes not in expression matrix — ", paste(missing_7a, collapse = ", "))
message(sprintf("  %d/%d genes found", length(genes_avail), length(genes_7a)))

df_7a <- data.frame(
  sample     = common_samples,
  risk_group = risk_groups,
  stringsAsFactors = FALSE
)
for (g in genes_avail) df_7a[[g]] <- as.numeric(expr_log[g, common_samples])

df_7a_long <- df_7a %>%
  pivot_longer(cols = all_of(genes_avail), names_to = "Gene", values_to = "Expression") %>%
  mutate(Gene = factor(Gene, levels = genes_avail))

# Wilcoxon comparisons
stat_7a <- df_7a_long %>%
  group_by(Gene) %>%
  summarise(
    pval = wilcox.test(Expression ~ risk_group)$p.value,
    padj = NA_real_,
    .groups = "drop"
  )
stat_7a$padj <- p.adjust(stat_7a$pval, method = "BH")

# Label positions (just above max)
label_pos <- df_7a_long %>%
  group_by(Gene) %>%
  summarise(y_pos = max(Expression, na.rm = TRUE) * 1.04, .groups = "drop") %>%
  left_join(stat_7a, by = "Gene") %>%
  mutate(
    sig_label = case_when(
      padj < 0.001 ~ "***",
      padj < 0.01  ~ "**",
      padj < 0.05  ~ "*",
      TRUE         ~ "ns"
    )
  )

p7a <- ggplot(df_7a_long, aes(x = risk_group, y = Expression, fill = risk_group)) +
  geom_violin(alpha = 0.35, color = NA, trim = TRUE) +
  geom_boxplot(width = 0.25, outlier.size = 0.2, outlier.alpha = 0.4, color = "grey30") +
  geom_text(
    data = label_pos,
    aes(x = 1.5, y = y_pos, label = sig_label),
    size = 2.2, inherit.aes = FALSE
  ) +
  facet_wrap(~Gene, nrow = 2, scales = "free_y") +
  scale_fill_manual(values = COLORS_RISK, name = "Risk") +
  labs(x = NULL, y = expression(log[2](TPM + 1)),
       title = "Immune Checkpoint Gene Expression") +
  theme_nature() +
  theme(legend.position = "right",
        axis.text.x = element_text(angle = 0, hjust = 0.5))

save_pdf(p7a, "Fig7A_checkpoint_genes.pdf", width_mm = 130, height_mm = 90)

# =============================================================================
# Fig 7B — IPS (Immunophenoscore): high vs low risk
# Main composite IPS box plot + 4 components in a 2-row layout
# Double-column width: 183 mm; height ~75 mm
# =============================================================================
message("\n=== Fig7B: IPS score ===")

ips_plot_df <- ips_df %>%
  dplyr::filter(!is.na(risk_group)) %>%
  mutate(risk_group = factor(risk_group, levels = c("Low", "High")))

ips_stat <- wilcox.test(IPS ~ risk_group, data = ips_plot_df)
ips_label <- if (ips_stat$p.value < 0.001) "p < 0.001" else
  sprintf("p = %.3f", ips_stat$p.value)

p7b_main <- ggplot(ips_plot_df, aes(x = risk_group, y = IPS, fill = risk_group)) +
  geom_boxplot(width = 0.5, outlier.size = 0.3, outlier.alpha = 0.4) +
  geom_jitter(width = 0.12, size = 0.15, alpha = 0.25, color = "grey40") +
  annotate("text", x = 1.5, y = max(ips_plot_df$IPS) * 1.03,
           label = ips_label, size = 2.2, fontface = "italic") +
  scale_fill_manual(values = COLORS_RISK) +
  labs(x = NULL, y = "Immunophenoscore (IPS)", title = "IPS") +
  theme_nature() +
  theme(legend.position = "none")

# 4 component sub-plots
components <- c("MHC", "Effector", "Suppressor", "Checkpoint")
comp_plots <- lapply(components, function(comp) {
  df_comp <- ips_plot_df %>% select(risk_group, val = all_of(comp))
  pval <- wilcox.test(val ~ risk_group, data = df_comp)$p.value
  plabel <- if (pval < 0.001) "***" else if (pval < 0.01) "**" else
            if (pval < 0.05) "*" else "ns"
  ggplot(df_comp, aes(x = risk_group, y = val, fill = risk_group)) +
    geom_boxplot(width = 0.5, outlier.size = 0.2, outlier.alpha = 0.3) +
    annotate("text", x = 1.5, y = max(df_comp$val) * 1.04,
             label = plabel, size = 2.0) +
    scale_fill_manual(values = COLORS_RISK) +
    labs(x = NULL, y = "ssGSEA score", title = comp) +
    theme_nature() +
    theme(legend.position = "none")
})

p7b <- p7b_main + wrap_plots(comp_plots, nrow = 1) +
  plot_layout(widths = c(1, 3)) +
  plot_annotation(theme = theme(plot.margin = margin(0, 0, 0, 0)))

save_pdf(p7b, "Fig7B_IPS.pdf", width_mm = 183, height_mm = 70)

# =============================================================================
# Fig 7D — Drug pathway enrichment bar chart
# drug_pathway_statistics.csv: Pathway, median_high, median_low, p_value, padj
# Bar chart of -log10(padj), coloured by direction (high > low vs low > high)
# Single-column: 89 mm
# =============================================================================
message("\n=== Fig7D: Drug pathway enrichment ===")

pw_df <- read.csv(file.path(RES_DIR, "drug_pathway_statistics.csv"),
                  stringsAsFactors = FALSE)
message(sprintf("  drug_pathway_statistics.csv: %d rows", nrow(pw_df)))
print(colnames(pw_df))

# Clean pathway names
pw_df <- pw_df %>%
  mutate(
    neg_log10_padj = -log10(padj + 1e-300),
    direction      = ifelse(median_high > median_low,
                            "Higher in High-risk", "Higher in Low-risk"),
    Pathway_clean  = gsub("_", " ", Pathway)
  ) %>%
  arrange(desc(neg_log10_padj))

# All pathways since there are only 8
pw_df$Pathway_clean <- factor(pw_df$Pathway_clean,
                               levels = rev(pw_df$Pathway_clean))

p7d <- ggplot(pw_df, aes(x = Pathway_clean, y = neg_log10_padj, fill = direction)) +
  geom_col(width = 0.65) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed",
             color = "grey50", linewidth = 0.4) +
  coord_flip() +
  scale_fill_manual(
    values = c("Higher in High-risk" = "#E64B35", "Higher in Low-risk" = "#4DBBD5"),
    name = "Activity"
  ) +
  labs(x = NULL, y = expression(-log[10](FDR)),
       title = "Drug-Pathway Activity by Risk Group") +
  theme_nature() +
  theme(legend.position = "bottom",
        legend.direction = "vertical")

save_pdf(p7d, "Fig7D_drug_pathways.pdf", width_mm = 89, height_mm = 80)

# =============================================================================
# Fig 7E — Volcano plot: DEGs between high vs low UIRS risk
# Source: deg_high_vs_low_risk.csv
# Thresholds: |logFC| > 1.5, FDR < 0.05
# UPR pathway genes highlighted in a distinct color
# =============================================================================
message("\n=== Fig7E: Volcano plot ===")

deg <- read.csv(file.path(RES_DIR, "deg_high_vs_low_risk.csv"),
                stringsAsFactors = FALSE)
message(sprintf("  deg_high_vs_low_risk.csv: %d rows", nrow(deg)))

# Significance thresholds
FC_CUTOFF  <- 1.5
FDR_CUTOFF <- 0.05

deg <- deg %>%
  mutate(
    sig_fc  = abs(logFC) > FC_CUTOFF,
    sig_fdr = adj.P.Val < FDR_CUTOFF,
    sig     = sig_fc & sig_fdr,
    is_upr  = gene %in% UPR_broad_genes,
    direction = case_when(
      sig & logFC >  FC_CUTOFF & is_upr  ~ "UPR Up",
      sig & logFC < -FC_CUTOFF & is_upr  ~ "UPR Down",
      sig & logFC >  FC_CUTOFF           ~ "Up",
      sig & logFC < -FC_CUTOFF           ~ "Down",
      TRUE                               ~ "NS"
    )
  )

# Count removed rows
n_inf <- sum(!is.finite(deg$adj.P.Val) | deg$adj.P.Val <= 0)
if (n_inf > 0) message(sprintf("  Removed %d rows with non-finite adj.P.Val", n_inf))
deg <- deg %>% filter(is.finite(adj.P.Val) & adj.P.Val > 0)

message(sprintf("  Significant (|FC|>%.1f & FDR<%.2f): %d up, %d down",
                FC_CUTOFF, FDR_CUTOFF,
                sum(deg$direction %in% c("Up", "UPR Up")),
                sum(deg$direction %in% c("Down", "UPR Down"))))
message(sprintf("  UPR genes significant: %d", sum(deg$is_upr & deg$sig)))

# Top 10 up and 10 down genes to label (prioritise UPR)
top_up   <- deg %>%
  filter(sig & logFC > 0) %>%
  arrange(desc(is_upr), adj.P.Val) %>%
  slice_head(n = 10)
top_down <- deg %>%
  filter(sig & logFC < 0) %>%
  arrange(desc(is_upr), adj.P.Val) %>%
  slice_head(n = 10)
label_genes <- bind_rows(top_up, top_down)

# Color palette — colorblind-safe
vol_colors <- c(
  "NS"       = "grey75",
  "Up"       = "#F4A582",
  "Down"     = "#92C5DE",
  "UPR Up"   = "#E64B35",
  "UPR Down" = "#2166AC"
)

p7e <- ggplot(deg, aes(x = logFC, y = -log10(adj.P.Val), color = direction)) +
  geom_point(size = 0.35, alpha = 0.55, stroke = 0) +
  geom_vline(xintercept = c(-FC_CUTOFF, FC_CUTOFF),
             linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_hline(yintercept = -log10(FDR_CUTOFF),
             linetype = "dashed", color = "grey50", linewidth = 0.3) +
  geom_text_repel(
    data = label_genes,
    aes(label = gene),
    size = 1.8,
    fontface = ifelse(label_genes$is_upr, "italic", "plain"),
    max.overlaps = 30,
    segment.size = 0.2,
    box.padding  = 0.3,
    point.padding = 0.2,
    min.segment.length = 0.1
  ) +
  scale_color_manual(
    values = vol_colors,
    labels = c("NS" = "Not significant",
               "Up" = "Up (non-UPR)", "Down" = "Down (non-UPR)",
               "UPR Up" = "UPR Up", "UPR Down" = "UPR Down"),
    name = NULL
  ) +
  labs(
    x     = expression(log[2]~"Fold Change (High / Low)"),
    y     = expression(-log[10]~"(FDR)"),
    title = "DEGs: UIRS High vs Low Risk"
  ) +
  guides(color = guide_legend(override.aes = list(size = 2), ncol = 1)) +
  theme_nature() +
  theme(legend.position = "right")

save_pdf(p7e, "Fig7E_volcano.pdf", width_mm = 120, height_mm = 100)

# =============================================================================
# Fig 7E_GO — GO enrichment dot plot (GSEA result from gseaResult object)
# Top 10 activated (NES > 0) + top 10 suppressed (NES < 0)
# Double-column width: 183 mm
# =============================================================================
message("\n=== Fig7E_GO: GO enrichment dotplot ===")

go_df <- as.data.frame(gsea_go)
go_sig <- go_df %>%
  filter(p.adjust < 0.05) %>%
  mutate(
    Description_clean = sub("^GOBP_", "", Description),
    Description_clean = gsub("_", " ", Description_clean),
    Description_clean = stringr::str_to_title(Description_clean)
  )

top_act  <- go_sig %>% filter(NES > 0) %>%
  arrange(p.adjust) %>% slice_head(n = 10)
top_supp <- go_sig %>% filter(NES < 0) %>%
  arrange(p.adjust) %>% slice_head(n = 10)
go_plot_df <- bind_rows(
  top_act  %>% mutate(Group = "Activated (NES > 0)"),
  top_supp %>% mutate(Group = "Suppressed (NES < 0)")
)

# Clamp gene ratio for display
go_plot_df <- go_plot_df %>%
  mutate(
    gene_ratio = setSize / max(setSize, na.rm = TRUE),
    neg_log10_padj = -log10(p.adjust),
    Description_clean = factor(
      Description_clean,
      levels = rev(go_plot_df$Description_clean[order(go_plot_df$NES)])
    )
  )

# Redesign: single-panel lollipop/bar plot with direction coloring
# Avoids facet overlap; shows all 20 terms on one y-axis ordered by NES
go_plot_df2 <- bind_rows(top_act, top_supp) %>%
  mutate(
    neg_log10_padj = -log10(p.adjust),
    Direction = ifelse(NES > 0, "Activated (UPR-high-risk)", "Suppressed (UPR-high-risk)"),
    Description_clean = factor(Description_clean,
                               levels = Description_clean[order(NES)])
  )

p7e_go <- ggplot(go_plot_df2,
                 aes(x = NES, y = Description_clean, color = Direction)) +
  geom_segment(aes(xend = 0, yend = Description_clean),
               linewidth = 0.4, color = "grey70") +
  geom_point(aes(size = setSize, alpha = neg_log10_padj)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40",
             linewidth = 0.5) +
  scale_color_manual(values = c("Activated (UPR-high-risk)" = "#B2182B",
                                "Suppressed (UPR-high-risk)" = "#2166AC"),
                     name = NULL) +
  scale_size_continuous(range = c(1.5, 5), name = "Gene set\nsize") +
  scale_alpha_continuous(range = c(0.5, 1),
                         name = expression(-log[10]~"(FDR)")) +
  labs(x = "Normalized Enrichment Score (NES)", y = NULL,
       title = "GO Biological Process GSEA") +
  theme_nature() +
  theme(
    axis.text.y     = element_text(size = 5.5),
    axis.text.x     = element_text(size = 6),
    legend.position = "right",
    legend.key.size = unit(3, "mm"),
    plot.title      = element_text(size = 7, face = "bold")
  )

save_pdf(p7e_go, "Fig7E_GO_dotplot.pdf", width_mm = 183, height_mm = 120)

message("\n=== Figure 7 complete ===")
message("All PDFs saved to: ", OUT_DIR)
