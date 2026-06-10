#!/usr/bin/env Rscript
# Assemble Supplementary Figures S3-S10 into single multi-panel PDFs with A/B/C tags.
# Source PDFs are rasterized to 300 dpi PNG (pdftoppm), white-margin trimmed in R,
# then composited with cowplot. PDF + PNG output per figure.
suppressPackageStartupMessages({
  library(png); library(grid); library(cowplot); library(ggplot2)
})
PROJ <- getwd()
FIG  <- file.path(PROJ, "figures")
ADD  <- file.path(PROJ, "manuscript/submission/Additional_file_1_Supplementary_Figures")
OUT  <- file.path(FIG, "supp_assembled")
TMP  <- file.path(OUT, "_raster")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(TMP, showWarnings = FALSE, recursive = TRUE)
DPI <- 300

# --- rasterize a single PDF page to trimmed PNG, return raster grob + aspect ---
# page = "auto" picks the first non-blank page (handles ggsurvplot blank p1).
raster_panel <- function(pdf_path, page = "auto", key = NULL, pad = 6) {
  stopifnot(file.exists(pdf_path))
  if (is.null(key)) key <- gsub("[^A-Za-z0-9]", "_", basename(pdf_path))
  render_page <- function(pg) {
    out_prefix <- file.path(TMP, sprintf("%s_p%d", key, pg))
    png_file <- paste0(out_prefix, ".png")
    cmd <- sprintf("pdftoppm -png -r %d -f %d -l %d -singlefile %s %s",
                   DPI, pg, pg, shQuote(pdf_path), shQuote(out_prefix))
    st <- system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
    if (st != 0 || !file.exists(png_file)) stop("pdftoppm failed for ", pdf_path, " p", pg)
    img <- readPNG(png_file)
    if (length(dim(img)) == 2) img <- array(rep(img, 3), dim = c(dim(img), 3))
    img
  }
  npages <- as.integer(system(sprintf("pdfinfo %s 2>/dev/null | grep Pages | tr -dc '0-9'",
                                       shQuote(pdf_path)), intern = TRUE))
  if (length(npages) == 0 || is.na(npages) || npages < 1) npages <- 1
  pick_pages <- if (identical(page, "auto")) seq_len(npages) else page
  img <- NULL; nonwhite <- NULL; rows <- integer(0); cols <- integer(0)
  for (pg in pick_pages) {
    img <- render_page(pg)
    rgb <- img[, , 1:3, drop = FALSE]
    nonwhite <- apply(rgb, c(1, 2), function(px) any(px < 0.97))
    rows <- which(apply(nonwhite, 1, any)); cols <- which(apply(nonwhite, 2, any))
    if (length(rows) > 5 && length(cols) > 5) break  # found content page
  }
  if (length(rows) && length(cols)) {
    r0 <- max(1, min(rows) - pad); r1 <- min(nrow(nonwhite), max(rows) + pad)
    c0 <- max(1, min(cols) - pad); c1 <- min(ncol(nonwhite), max(cols) + pad)
    img <- img[r0:r1, c0:c1, , drop = FALSE]
  }
  list(grob = rasterGrob(img, interpolate = TRUE),
       aspect = dim(img)[2] / dim(img)[1])  # width/height
}

panel_plot <- function(rp) ggdraw() + draw_grob(rp$grob)

assemble <- function(name, panels, ncol, rel_widths = NULL, rel_heights = NULL,
                     base_w = 7.2, label_size = 16) {
  plots <- lapply(panels, panel_plot)
  if (length(panels) == 1) {
    g <- plots[[1]]                       # single-panel: no tag, no grid
  } else {
    pg_args <- list(plotlist = plots, ncol = ncol, labels = LETTERS[seq_along(panels)],
                    label_size = label_size, label_fontface = "bold")
    if (!is.null(rel_widths))  pg_args$rel_widths  <- rel_widths
    if (!is.null(rel_heights)) pg_args$rel_heights <- rel_heights
    g <- do.call(plot_grid, pg_args)
  }
  # estimate height from aspect ratios
  nrow <- ceiling(length(panels) / ncol)
  asp <- sapply(panels, function(p) p$aspect)
  med_asp <- median(asp)
  per_w <- base_w / ncol
  per_h <- per_w / med_asp
  total_h <- per_h * nrow
  pdf_out <- file.path(OUT, paste0(name, ".pdf"))
  png_out <- file.path(OUT, paste0(name, ".png"))
  ggsave(pdf_out, g, width = base_w, height = total_h, limitsize = FALSE)
  ggsave(png_out, g, width = base_w, height = total_h, dpi = 200, limitsize = FALSE)
  cat(sprintf("  %s: %d panel(s), %.1f x %.1f in -> %s\n",
              name, length(panels), base_w, total_h, basename(pdf_out)))
}

cat("=== Assembling supplementary figures ===\n")

## ---------- S3: K2/K3 consensus comparison ----------
# (A) consensus matrices K=2 (page2) + K=3 (page3) side by side
# (B) K=3 KM full cohort + IDH-WT
# (C) recomputed pairwise log-rank panel
ccpdf <- file.path(PROJ, "results/consensus_clustering/consensus.pdf")
s3_cm_k2 <- raster_panel(ccpdf, page = 2, key = "cc_k2")
s3_cm_k3 <- raster_panel(ccpdf, page = 3, key = "cc_k3")
# Build composite A (two consensus matrices)
A_grid <- plot_grid(panel_plot(s3_cm_k2), panel_plot(s3_cm_k3), ncol = 2)
s3_kmK3_full  <- raster_panel(file.path(FIG, "Fig3_km_K3.pdf"), page = "auto", key = "km_k3_full")
s3_kmK3_idhwt <- raster_panel(file.path(FIG, "Fig3_km_idhwt_K3.pdf"), page = "auto", key = "km_k3_idhwt")
B_grid <- plot_grid(panel_plot(s3_kmK3_full), panel_plot(s3_kmK3_idhwt), ncol = 2)
s3_C <- raster_panel(file.path(OUT, "_panel_S3C_pairwise_logrank.pdf"), page = "auto", key = "s3c")
s3 <- plot_grid(
  A_grid, B_grid, panel_plot(s3_C),
  ncol = 1, labels = c("A", "B", "C"), label_size = 16, label_fontface = "bold",
  rel_heights = c(1.0, 1.0, 0.9))
ggsave(file.path(OUT, "FigureS3.pdf"), s3, width = 8.2, height = 11.0, limitsize = FALSE)
ggsave(file.path(OUT, "FigureS3.png"), s3, width = 8.2, height = 11.0, dpi = 200, limitsize = FALSE)
cat("  FigureS3: A(2 consensus matrices)+B(2 KM)+C(pairwise logrank)\n")

## ---------- S4: CoxPH permutation (single) ----------
s4 <- raster_panel(file.path(FIG, "SuppFig_permutation_uirs_coxph.pdf"), page = "auto", key = "s4")
assemble("FigureS4", list(s4), ncol = 1, base_w = 6.5)

## ---------- S5: sensitivity (A-D) ----------
s5a <- raster_panel(file.path(FIG, "FigS_sensitivity_all_genes_K2.pdf"), key = "s5a")
s5b <- raster_panel(file.path(FIG, "FigS_sensitivity_all_genes_K3.pdf"), key = "s5b")
s5c <- raster_panel(file.path(FIG, "FigS_sensitivity_hallmark_K2.pdf"), key = "s5c")
s5d <- raster_panel(file.path(FIG, "FigS_sensitivity_hallmark_K3.pdf"), key = "s5d")
assemble("FigureS5", list(s5a, s5b, s5c, s5d), ncol = 2, base_w = 8.0)

## ---------- S6: IDH-WT immune ssGSEA (single, heatmap of 27 cells) ----------
s6 <- raster_panel(file.path(ADD, "FigS_immune_idhwt.pdf"), page = "auto", key = "s6")
assemble("FigureS6", list(s6), ncol = 1, base_w = 7.5)

## ---------- S7: MCPcounter (single) ----------
s7 <- raster_panel(file.path(FIG, "FigS_mcpcounter_boxplots.pdf"), page = "auto", key = "s7")
assemble("FigureS7", list(s7), ncol = 1, base_w = 7.0)

## ---------- S8: ML C-index violin (A) + RSF time-ROC (B) + RSF perm (C) ----------
s8a <- raster_panel(file.path(ADD, "FigS_ml_cindex.pdf"), page = "auto", key = "s8a")
s8b <- raster_panel(file.path(OUT, "_panel_S8B_rsf_timeroc.pdf"), page = "auto", key = "s8b")
s8c <- raster_panel(file.path(FIG, "SuppFig_permutation_test.pdf"), page = "auto", key = "s8c")
# A is tall (100 combos); place A on top full width, B+C below side by side
BC <- plot_grid(panel_plot(s8b), panel_plot(s8c), ncol = 2,
                labels = c("B", "C"), label_size = 16, label_fontface = "bold")
s8 <- plot_grid(panel_plot(s8a), BC, ncol = 1,
                labels = c("A", ""), label_size = 16, label_fontface = "bold",
                rel_heights = c(1.5, 1.0))
ggsave(file.path(OUT, "FigureS8.pdf"), s8, width = 8.2, height = 10.5, limitsize = FALSE)
ggsave(file.path(OUT, "FigureS8.png"), s8, width = 8.2, height = 10.5, dpi = 200, limitsize = FALSE)
cat("  FigureS8: A(ml cindex)+B(RSF timeROC)+C(RSF permutation)\n")

## ---------- S9: DCA 1yr + 5yr ----------
s9a <- raster_panel(file.path(ADD, "FigS_dca_1yr.pdf"), page = "auto", key = "s9a")
s9b <- raster_panel(file.path(ADD, "FigS_dca_5yr.pdf"), page = "auto", key = "s9b")
assemble("FigureS9", list(s9a, s9b), ncol = 2, base_w = 9.0)

## ---------- S10: 5 hub gene KM ----------
hub <- c("CASP4","VEGFA","BLOC1S1","DDOST","SLC7A5")
s10 <- lapply(hub, function(g)
  raster_panel(file.path(FIG, "part4_clinical", sprintf("Fig8_km_%s.pdf", g)),
               page = "auto", key = paste0("s10_", g)))
# 5 panels: 2x3 grid (last cell empty)
plots <- lapply(s10, panel_plot)
labs <- LETTERS[1:5]
s10g <- plot_grid(plotlist = plots, ncol = 2, labels = labs,
                  label_size = 16, label_fontface = "bold")
ggsave(file.path(OUT, "FigureS10.pdf"), s10g, width = 9.5, height = 12.0, limitsize = FALSE)
ggsave(file.path(OUT, "FigureS10.png"), s10g, width = 9.5, height = 12.0, dpi = 200, limitsize = FALSE)
cat("  FigureS10: 5 hub gene KM (CASP4,VEGFA,BLOC1S1,DDOST,SLC7A5)\n")

cat("\n=== Done. Outputs in", OUT, "===\n")
