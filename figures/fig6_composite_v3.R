###############################################################################
# fig6_composite_v3.R  — Figure 6 hub-gene validation (5-panel redesign)
#
# Panel A: CoxPH coefficients, all 15 UIRS genes (risk-promoting vs protective)
# Panel B: Forest plot, hub-gene prognostic value in TCGA IDH-wildtype (median-split)
# Panel C: KM OS in TCGA IDH-WT — CASP4, BLOC1S1, SLC7A5 (3-in-a-row)
# Panel D: STRING functional-enrichment barplot (replaces correlation heatmap)
# Panel E: scRNA UMAP feature plots — CASP4, BLOC1S1, SLC7A5 (GSE182109)
#
# Output: Figures_v2/Figure6_composite.pdf  and  .png
###############################################################################


suppressPackageStartupMessages({
  library(SeuratObject)          # must load before Seurat to avoid RcppAnnoy error
  library(Seurat, warn.conflicts = FALSE)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(patchwork)
  library(survival)
  library(survminer)
  library(cowplot)
  library(viridis)
  library(png)
  library(grid)
  library(grDevices)
  library(showtext)
  library(sysfonts)
  library(scales)
})

# ── Fonts ──────────────────────────────────────────────────────────────────
font_add("Arial",
         regular    = "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
         bold       = "/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf",
         italic     = "/usr/share/fonts/truetype/dejavu/DejaVuSans-Oblique.ttf",
         bolditalic = "/usr/share/fonts/truetype/dejavu/DejaVuSans-BoldOblique.ttf")
showtext_auto()
showtext_opts(dpi = 300)
font_add("Times",
         regular    = "/usr/share/fonts/truetype/dejavu/DejaVuSerif.ttf",
         bold       = "/usr/share/fonts/truetype/dejavu/DejaVuSerif-Bold.ttf",
         italic     = "/usr/share/fonts/truetype/dejavu/DejaVuSerif-Italic.ttf",
         bolditalic = "/usr/share/fonts/truetype/dejavu/DejaVuSerif-BoldItalic.ttf")
FONT <- "Arial"
LABEL_FONT <- "Times"  # panel labels (A-E) match Fig3 (Times bold 13)

# ── Paths ──────────────────────────────────────────────────────────────────
PROJ <- getwd()
DATA_PROC <- file.path(PROJ, "data", "processed")
RES_DIR   <- file.path(PROJ, "results")
OUT_DIR   <- file.path(PROJ, "manuscript", "submission", "Figures_v2")
TMP_DIR   <- file.path(OUT_DIR, "tmp_fig6v3")
dir.create(TMP_DIR, recursive = TRUE, showWarnings = FALSE)

# ── Constants ──────────────────────────────────────────────────────────────
W_MM       <- 183          # double-column Nature width
H_MM       <- 160          # target height (compact 4-row layout, no footnote)
MM_IN      <- 1 / 25.4
BASE       <- 7.5
COLORS_DIR <- c("Risk" = "#E64B35", "Protective" = "#4DBBD5")
HUB_KM     <- c("CASP4", "BLOC1S1", "SLC7A5")
HUB_ALL    <- c("CASP4", "VEGFA", "BLOC1S1", "DDOST", "SLC7A5")

# ── Theme ──────────────────────────────────────────────────────────────────
theme_nat <- function(bs = BASE) {
  theme_classic(base_size = bs, base_family = FONT) +
    theme(
      axis.text        = element_text(size = bs - 1,   color = "black"),
      axis.title       = element_text(size = bs,        color = "black"),
      plot.title       = element_text(size = bs,        face = "bold", hjust = 0.5),
      legend.text      = element_text(size = bs - 1.5),
      legend.title     = element_text(size = bs - 1,   face = "bold"),
      legend.key.size  = unit(3, "mm"),
      strip.text       = element_text(size = bs - 1,   face = "bold"),
      strip.background = element_blank(),
      panel.grid       = element_blank(),
      plot.margin      = margin(2, 3, 2, 3, "mm")
    )
}

###############################################################################
# Panel A — CoxPH coefficients (15 UIRS genes) → save as PNG for proper sizing
###############################################################################
message("Panel A: UIRS model coefficients ...")

coef_df <- read.csv(file.path(RES_DIR, "coxph_model_coefficients.csv"),
                    stringsAsFactors = FALSE)
stopifnot(all(c("Gene", "Beta") %in% colnames(coef_df)))
message("  Loaded ", nrow(coef_df), " genes")

coef_df <- coef_df %>%
  mutate(
    direction  = ifelse(Beta > 0, "Risk", "Protective"),
    label_text = sprintf("%.3f", Beta)
  ) %>%
  arrange(Beta)
coef_df$Gene <- factor(coef_df$Gene, levels = coef_df$Gene)
slc_mask <- as.character(coef_df$Gene) == "SLC7A5"

p_A_raw <- ggplot(coef_df, aes(x = Beta, y = Gene, fill = direction)) +
  geom_col(width = 0.65) +
  geom_vline(xintercept = 0, linewidth = 0.4, color = "black") +
  geom_text(
    aes(label = label_text, hjust = ifelse(Beta > 0, -0.12, 1.12)),
    size     = 2.0, family = FONT,
    color    = ifelse(slc_mask, "#8B0000", "grey30"),
    fontface = ifelse(slc_mask, "bold",    "plain")
  ) +
  scale_fill_manual(values = COLORS_DIR, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.32, 0.32))) +
  scale_y_discrete(expand  = expansion(add  = c(0.4, 0.9))) +
  coord_cartesian(clip = "off") +
  labs(x = "LASSO-CoxPH coefficient", y = NULL) +
  theme_nat() +
  theme(
    axis.text.y      = element_text(size = 6.0, face = "italic"),
    legend.position  = c(0.97, 0.05),
    legend.justification = c("right", "bottom"),
    legend.background = element_blank(),
    plot.margin = margin(4, 8, 2, 9, "mm")
  )

# Panel A kept as vector ggplot object (p_A_raw); no PNG round-trip.
message("  Panel A built (vector)")

###############################################################################
# Panel B — Forest plot: hub-gene prognostic value in TCGA IDH-WT (median-split)
###############################################################################
message("Panel B: hub-gene forest plot (TCGA IDH-WT, median-split) ...")

load(file.path(DATA_PROC, "risk_model_final.RData"))
load(file.path(DATA_PROC, "tcga_glioma_expression.RData"))

common_samples <- intersect(names(risk_score_train), colnames(expr_tpm_symbol))
expr_log       <- log2(expr_tpm_symbol[, common_samples] + 1)
clin_df        <- clinical_valid[match(common_samples, clinical_valid$barcode), ]

idhwt_idx   <- !is.na(clin_df$IDH_status) & clin_df$IDH_status == "WT"
idhwt_samps <- common_samples[idhwt_idx]
idhwt_clin  <- clin_df[idhwt_idx, ]
idhwt_time  <- idhwt_clin$OS.time / 30   # days → months
idhwt_stat  <- idhwt_clin$OS
message("  TCGA IDH-WT samples: ", length(idhwt_samps))

hub_cox_rows <- lapply(HUB_ALL, function(gene) {
  if (!gene %in% rownames(expr_log)) return(NULL)
  ge  <- as.numeric(expr_log[gene, idhwt_samps])
  grp <- factor(ifelse(ge > median(ge, na.rm = TRUE), "High", "Low"),
                levels = c("Low", "High"))
  td  <- data.frame(time = idhwt_time, status = idhwt_stat, group = grp)
  td  <- td[complete.cases(td) & td$time > 0, ]
  if (nrow(td) < 20) return(NULL)
  cf <- tryCatch(coxph(Surv(time, status) ~ group, data = td), error = function(e) NULL)
  if (is.null(cf)) return(NULL)
  s  <- summary(cf)
  # Significance star uses the median-split log-rank p
  lr <- survdiff(Surv(time, status) ~ group, data = td)
  pv <- 1 - pchisq(lr$chisq, df = 1)
  data.frame(Gene = gene, HR = s$conf.int[1,1],
             HR_lower = s$conf.int[1,3], HR_upper = s$conf.int[1,4],
             P = pv, stringsAsFactors = FALSE)
})
hub_cox_df <- do.call(rbind, Filter(Negate(is.null), hub_cox_rows)) %>%
  mutate(
    sig_label = case_when(P < 0.001 ~ "***", P < 0.01 ~ "**", P < 0.05 ~ "*", TRUE ~ "ns"),
    P_label   = ifelse(P < 0.001, formatC(P, format="e", digits=2), sprintf("%.3f", P))
  ) %>%
  arrange(desc(HR))
hub_cox_df$Gene <- factor(hub_cox_df$Gene, levels = hub_cox_df$Gene)

x_annot <- max(hub_cox_df$HR_upper) * 1.75

p_B_raw <- ggplot(hub_cox_df, aes(x = HR, y = Gene)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50", linewidth = 0.4) +
  geom_errorbar(aes(xmin = HR_lower, xmax = HR_upper),
                width = 0.22, linewidth = 0.5, orientation = "y",
                color = ifelse(hub_cox_df$HR > 1, COLORS_DIR["Risk"], COLORS_DIR["Protective"])) +
  geom_point(aes(color = ifelse(HR > 1, "Risk", "Protective")), size = 2.5, shape = 18) +
  geom_text(aes(x = x_annot * 0.97, label = P_label),
            size = 2.0, hjust = 1, family = FONT) +
  geom_text(aes(x = x_annot * 1.03, label = sig_label),
            size = 2.0, hjust = 0, family = FONT) +
  annotate("text", x = x_annot * 0.97, y = nrow(hub_cox_df) + 0.75,
           label = "p-value", size = 2.0, hjust = 1, fontface = "bold", family = FONT) +
  scale_color_manual(values = COLORS_DIR, name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.08)),
                     limits = c(min(hub_cox_df$HR_lower) * 0.80, x_annot * 1.08)) +
  coord_cartesian(clip = "off") +
  labs(x = "Hazard Ratio (95% CI)", y = NULL, title = "Hub Genes in IDH-Wildtype") +
  theme_nat() +
  theme(axis.text.y = element_text(size = 7, face = "italic"),
        legend.position = "bottom",
        plot.margin = margin(2, 12, 2, 3, "mm"))

# Panel B kept as vector ggplot object (p_B_raw); no PNG round-trip.
message("  Panel B built (vector)  (", nrow(hub_cox_df), " genes)")

###############################################################################
# Panel C — KM curves: CASP4, BLOC1S1, SLC7A5 in TCGA IDH-WT (3-in-a-row) → PNG
###############################################################################
message("Panel C: KM curves (TCGA IDH-WT, 3 genes, median-split) ...")

km_grobs <- list()
for (gene in HUB_KM) {
  if (!gene %in% rownames(expr_log)) {
    message("  ", gene, " missing — skipping"); next
  }
  ge  <- as.numeric(expr_log[gene, idhwt_samps])
  grp <- factor(ifelse(ge > median(ge, na.rm = TRUE), "High", "Low"), levels = c("Low","High"))
  sd  <- data.frame(time = idhwt_time, status = idhwt_stat, group = grp)
  sd  <- sd[complete.cases(sd) & sd$time > 0, ]
  fit <- survfit(Surv(time, status) ~ group, data = sd)
  lr  <- survdiff(Surv(time, status) ~ group, data = sd)
  pv  <- 1 - pchisq(lr$chisq, df = 1)
  cox_hr <- exp(coef(coxph(Surv(time, status) ~ group, data = sd))["groupHigh"])
  plab   <- if (pv < 0.001) "p < 0.001" else sprintf("p = %.3f", pv)
  hrlab  <- sprintf("HR = %.2f", cox_hr)

  kp <- ggsurvplot(
    fit, data = sd,
    palette            = c(COLORS_DIR["Protective"], COLORS_DIR["Risk"]),
    risk.table         = TRUE,
    pval               = FALSE,
    conf.int           = TRUE,
    conf.int.alpha     = 0.15,
    legend.labs        = c("Low", "High"),
    legend.title       = NULL,
    xlab               = "Time (months)",
    ylab               = "OS probability",
    title              = sprintf("%s\n%s,  %s", gene, hrlab, plab),
    ggtheme            = theme_nat(bs = 6.5) +
      theme(plot.title = element_text(size = 6.5, hjust = 0.5, face = "bold.italic")),
    fontsize           = 2.4,
    risk.table.fontsize = 2.2,
    risk.table.height  = 0.28,
    tables.theme       = theme_cleantable() +
      theme(text = element_text(size = 5.5, family = FONT),
            plot.margin = margin(0, 1, 1, 1, "mm"))
  )
  km_grobs[[gene]] <- cowplot::plot_grid(
    kp$plot + theme(plot.margin    = margin(1, 1, 2, 1, "mm"),
                    axis.title.x   = element_text(margin = margin(t = 1.5, b = 1.0)),
                    legend.key.size = unit(2.5, "mm"),
                    legend.text    = element_text(size = 5.5)),
    kp$table + theme(plot.margin  = margin(1.5, 1, 1, 1, "mm")),
    ncol = 1, rel_heights = c(0.72, 0.28)
  )
}
# Panel C kept as vector cowplot grob object (km_row); no PNG round-trip.
km_row <- cowplot::plot_grid(plotlist = km_grobs, nrow = 1)
message("  Panel C built (vector)")

###############################################################################
# Panel D — STRING functional-enrichment horizontal barplot
###############################################################################
message("Panel D: STRING functional-enrichment barplot ...")

enrich_df <- read.csv(file.path(RES_DIR, "hub_gene_ppi_enrichment.csv"),
                      stringsAsFactors = FALSE)
stopifnot(all(c("category", "description", "fdr") %in% colnames(enrich_df)))
message("  Loaded ", nrow(enrich_df), " enrichment terms")
print(enrich_df[, c("category", "description", "fdr")])

enrich_top <- enrich_df %>%
  arrange(fdr) %>%
  slice_head(n = 8) %>%
  mutate(
    neg_log10_fdr = -log10(fdr),
    # Concise labels that fit within the plot area
    label_short = case_when(
      grepl("Intrinsic apoptotic", description) ~
        "Intrinsic apoptotic signaling (ER stress)",
      grepl("Positive regulation of transcription", description) ~
        "Positive reg. of RNAPII transcription (ER stress)",
      grepl("PERK-mediated", description) ~
        "PERK-mediated unfolded protein response",
      grepl("Protein processing", description) ~
        "Protein processing in ER (KEGG hsa04141)",
      grepl("CHOP-ATF4", description) ~
        "CHOP-ATF4 complex",
      grepl("Regulation of vascular", description) ~
        "Regulation of vascular SMC apoptosis",
      TRUE ~ description
    ),
    category = case_when(
      category == "Component" ~ "GO Component",
      category == "Process"   ~ "GO Process",
      category == "KEGG"      ~ "KEGG Pathway",
      TRUE                    ~ category
    )
  ) %>%
  arrange(neg_log10_fdr)

enrich_top$label_short <- factor(enrich_top$label_short, levels = enrich_top$label_short)

cat_colors <- c("GO Component" = "#0072B2", "GO Process" = "#009E73",
                "KEGG Pathway" = "#E69F00")

p_D_raw <- ggplot(enrich_top, aes(x = neg_log10_fdr, y = label_short, fill = category)) +
  geom_col(width = 0.82, color = NA) +
  geom_vline(xintercept = -log10(0.05), linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  # Place FDR labels just outside each bar's right end in dark text, so they are
  # never low-contrast against the coloured fill nor clipped at the left edge.
  geom_text(aes(x = neg_log10_fdr, label = sprintf("FDR=%.4f", fdr)),
            hjust = -0.10, size = 1.9, family = FONT, color = "grey15") +
  scale_fill_manual(values = cat_colors, name = "Category") +
  # Tight right expansion: just enough to seat the "FDR=…" label without leaving
  # a large empty strip, so bars run nearly full-width.
  scale_x_continuous(expand = expansion(mult = c(0, 0.12)),
                     breaks = pretty_breaks(n = 4)) +
  scale_y_discrete(expand = expansion(add = c(0.55, 0.55))) +
  coord_cartesian(clip = "off") +
  labs(x = expression(-log[10](FDR)), y = NULL,
       title = "STRING Functional Enrichment (Hub Genes & Interactors)") +
  theme_nat() +
  theme(
    axis.text.y  = element_text(size = 5.8, hjust = 1, lineheight = 1.05),
    # Horizontal legend along the bottom; fills former right-side whitespace
    # instead of floating as an isolated box in the empty lower-right corner.
    legend.position   = "bottom",
    legend.direction  = "horizontal",
    legend.title      = element_text(size = 6, face = "bold"),
    legend.text       = element_text(size = 5.5),
    legend.key.size   = unit(2.5, "mm"),
    legend.margin     = margin(0, 0, 0, 0),
    legend.box.margin = margin(t = -1, b = 0),
    plot.margin = margin(2, 8, 1, 4, "mm")
  )

# Panel D kept as vector ggplot object (p_D_raw); no PNG round-trip.
message("  Panel D built (vector)")
message("  Terms shown:")
for (tt in as.character(enrich_top$label_short))
  message("    - ", tt)

###############################################################################
# Panel E — scRNA UMAP feature plots (CASP4, BLOC1S1, SLC7A5) — GSE182109
###############################################################################
message("Panel E: scRNA UMAP feature plots ...")

sc_file <- file.path(PROJ, "results", "GSE182109_validation",
                     "checkpoint_step4_annotated.rds")
seu <- readRDS(sc_file)
message("  Seurat dims: ", paste(dim(seu), collapse = " x "))

sct_assay <- seu@assays[["SCT"]]
data_mat  <- slot(sct_assay, "data")
umap_emb  <- seu@reductions[["umap"]]@cell.embeddings

genes_sc <- HUB_KM[HUB_KM %in% rownames(data_mat)]
stopifnot(length(genes_sc) == 3)
message("  Hub genes in SCT: ", paste(genes_sc, collapse = ", "))

sc_df <- data.frame(UMAP1 = umap_emb[, 1], UMAP2 = umap_emb[, 2])
for (g in genes_sc) sc_df[[g]] <- as.numeric(data_mat[g, ])

feat_plots <- lapply(genes_sc, function(g) {
  df_g <- sc_df[order(sc_df[[g]]), ]
  ggplot(df_g, aes(x = UMAP1, y = UMAP2, color = .data[[g]])) +
    geom_point(size = 0.10, stroke = 0, alpha = 0.65) +
    scale_color_viridis_c(
      option = "C",
      name   = expression(log[2](SCT+1)),
      guide  = guide_colorbar(
        barwidth  = unit(2.0, "mm"), barheight = unit(12, "mm"),
        title.theme = element_text(size = 5.0, family = FONT),
        label.theme = element_text(size = 4.8, family = FONT)
      )
    ) +
    labs(title = g) +
    coord_equal() +
    theme_void(base_family = FONT) +
    theme(
      plot.title   = element_text(size = 7.5, face = "italic", hjust = 0.5,
                                   margin = margin(b = 1)),
      legend.text  = element_text(size = 5.0, family = FONT),
      legend.title = element_text(size = 5.0, family = FONT),
      plot.margin  = margin(1, 3, 1, 3, "mm")
    )
})

pE_png <- file.path(TMP_DIR, "pE.png")
# Save E at half-width so UMAP proportions look right
p_E_raw <- wrap_plots(feat_plots, nrow = 1)
ggsave(pE_png, plot = p_E_raw, width = 183 * MM_IN, height = 42 * MM_IN,
       units = "in", dpi = 400, bg = "white")
message("  Panel E saved: ", pE_png)

rm(seu, sct_assay, data_mat); gc()

###############################################################################
# Assemble — read all PNG panels, add bold 13 pt labels, composite with cowplot
###############################################################################
message("=== Assembling 5-panel composite ===")

make_labeled <- function(png_path, label, label_x = 0.008, label_y = 0.992) {
  img <- png::readPNG(png_path)
  cowplot::ggdraw() +
    cowplot::draw_grob(grid::rasterGrob(img, interpolate = TRUE)) +
    cowplot::draw_label(label, x = label_x, y = label_y,
                        hjust = 0, vjust = 1,
                        fontface = "bold", size = 13,
                        fontfamily = LABEL_FONT, color = "black")
}

# Vector-preserving label wrapper: draws the ggplot/grob object directly onto
# the cowplot canvas (no PNG round-trip), so panels stay as vector paths + text.
make_labeled_vec <- function(gg, label, label_x = 0.008, label_y = 0.992) {
  cowplot::ggdraw() +
    cowplot::draw_plot(gg) +
    cowplot::draw_label(label, x = label_x, y = label_y,
                        hjust = 0, vjust = 1,
                        fontface = "bold", size = 13,
                        fontfamily = LABEL_FONT, color = "black")
}

# A/B/C/D kept as vector ggplot/grob objects (no PNG round-trip); only E (≈200k
# scatter points) stays rasterized for file-size/render sanity.
pA_lbl <- make_labeled_vec(p_A_raw, "A")
pB_lbl <- make_labeled_vec(p_B_raw, "B")
pC_lbl <- make_labeled_vec(km_row,  "C")
pD_lbl <- make_labeled_vec(p_D_raw, "D")
pE_lbl <- make_labeled(pE_png, "E")

# Row 1: A (half) | B (half) — both at same 88 mm output width
row1 <- cowplot::plot_grid(pA_lbl, pB_lbl, nrow = 1, rel_widths = c(1, 1))
# Row 2: C full-width KM
# Row 3: D full-width enrichment
# Row 4: E full-width scRNA

# Heights in relative units (layout ratios per row):
# Row1 = 72 (A|B), Row2 = 78 (KM), Row3 = 50 (D, denser now that bars fill the
# width and the legend moved to the bottom), Row4 = 42 (E scRNA).
rel_h <- c(72, 78, 50, 42)
rel_h_norm <- rel_h / sum(rel_h)

p_composite <- cowplot::plot_grid(
  row1, pC_lbl, pD_lbl, pE_lbl,
  ncol = 1,
  rel_heights = rel_h_norm
)

out_pdf <- file.path(OUT_DIR, "Figure6_composite.pdf")
out_png <- file.path(OUT_DIR, "Figure6_composite.png")

ggsave(out_pdf, plot = p_composite, device = cairo_pdf,
       width = W_MM * MM_IN, height = H_MM * MM_IN, units = "in", dpi = 300)
message("Saved PDF: ", out_pdf)

ggsave(out_png, plot = p_composite,
       width = W_MM * MM_IN, height = H_MM * MM_IN, units = "in", dpi = 300)
message("Saved PNG: ", out_png)

# Verify output sizes
dims <- dim(png::readPNG(out_png))
message(sprintf("PNG verified: %d x %d px  (%.0f x %.0f mm at 300 dpi)",
                dims[2], dims[1],
                dims[2]/300*25.4, dims[1]/300*25.4))

unlink(TMP_DIR, recursive = TRUE)
message("=== Figure 6 (v3) complete ===")
message(sprintf("Target: %d x %d mm  (aspect = %.2f)", W_MM, H_MM, W_MM / H_MM))
message("Panel D: STRING enrichment barplot — all 6 terms from hub_gene_ppi_enrichment.csv")
message("Removed: correlation heatmap, expression heatmap")
