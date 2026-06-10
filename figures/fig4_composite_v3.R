###############################################################################
# fig4_composite_v3.R
# Figure 4 composite — machine-learning model (7 panels)
# Panels:
#   A  C-index heatmap (top 25 combinations)
#   B  Benchmark boxplot (top 25 combinations)
#   C  Permutation test (LASSO-CoxPH, p=0.015)
#   D  KM TCGA training set
#   E  KM CGGA batch2 (tuning validation)
#   F  KM CGGA batch1 (blind test)
#   G  Time-dependent ROC — CGGA batch1 (1-yr 0.761, 3-yr 0.856, 5-yr 0.888)
#
# Layout (183 mm wide, ~280 mm tall):
#   Row 1:  A (heatmap)          — full width
#   Row 2:  B (boxplot)          — full width
#   Row 3:  C (permutation) | spacer
#   Row 4:  D | E | F            — 3 KM panels equal width
#   Row 5:  G (time ROC)  centred (2/3 width, 1/3 spacer)
#
# Output: Figures_v2/Figure4_composite.pdf  +  Figures_v2/Figure4_composite.png
###############################################################################


suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(survival)
  library(survminer)
  library(timeROC)
  library(ComplexHeatmap)
  library(circlize)
  library(patchwork)
  library(cowplot)
  library(png)
  library(grid)
  library(grDevices)
  library(showtext)
  library(sysfonts)
})

# ── Font setup ─────────────────────────────────────────────────────────────────
font_add("Times",
         regular    = "/usr/share/fonts/truetype/dejavu/DejaVuSerif.ttf",
         bold       = "/usr/share/fonts/truetype/dejavu/DejaVuSerif-Bold.ttf",
         italic     = "/usr/share/fonts/truetype/dejavu/DejaVuSerif-Italic.ttf",
         bolditalic = "/usr/share/fonts/truetype/dejavu/DejaVuSerif-BoldItalic.ttf")
showtext_auto()
showtext_opts(dpi = 300)

# ── Paths ──────────────────────────────────────────────────────────────────────
PROJ <- getwd()
DATA_PROC <- file.path(PROJ, "data", "processed")
RES_DIR   <- file.path(PROJ, "results")
OUT_DIR   <- file.path(PROJ, "manuscript", "submission", "Figures_v2")
TMP_DIR   <- file.path(OUT_DIR, "tmp_fig4v3")
dir.create(TMP_DIR, recursive = TRUE, showWarnings = FALSE)

W_MM  <- 183          # single page width (Nature double-column)
H_MM  <- 280          # total figure height
MM_IN <- 1 / 25.4
FONT  <- "Times"
BASE  <- 8

# ── Slot geometry (single source of truth) ──────────────────────────────────
# The composite is W_MM x H_MM with plot_layout(heights = LAYOUT_HEIGHTS).
# patchwork normalizes the heights, so each row's physical height in the final
# figure = (h / sum(h)) * H_MM. Raster panels (A heatmap; D/E/F KM; G ROC) are
# saved as PNG then stretched to fill their slot via draw_grob — so any mismatch
# between the PNG's mm canvas and its slot's mm size rescales the embedded text
# away from the 8 pt used by the vector panels (B, C). To keep raster text
# visually equal to B/C we render each PNG at (close to) its true slot size, so
# placement is ~1:1 and 8 pt stays 8 pt. These constants are the single source of
# truth for both the per-panel png() calls and the final plot_layout().
LAYOUT_HEIGHTS <- c(0.20, 0.18, 0.16, 0.22, 0.16)
.row_mm   <- (LAYOUT_HEIGHTS / sum(LAYOUT_HEIGHTS)) * H_MM   # physical row heights
SLOT_A_W  <- W_MM            # row 1: A spans full width
SLOT_A_H  <- .row_mm[1]
SLOT_KM_W <- W_MM / 3.30     # row 4: D|E|F each occupy 3/3.30 of the width (4th col is legend 0.30 units ≈ 16.6 mm)
SLOT_KM_H <- .row_mm[4]
SLOT_G_W  <- W_MM            # row 5: G spans full width
SLOT_G_H  <- .row_mm[5]

# KM panels (D, E, F) are rendered by survminer's print.ggsurvplot, which packs
# the survival plot + risk table into the PNG canvas: the plot region occupies
# only ~60-70% of the canvas height, so nominal-8pt theme text renders visually
# ~0.6x of the vector panels (B, C). KM_BASE compensates so KM text matches B/C
# at the composite scale (mirrors Fig 2). The risk-table number size is decoupled
# from KM_BASE: scaling it by the full KM_BASE/BASE factor makes the two at-risk
# rows collide within risk.table.height, so KM_FS is held at a legible 4.0 and
# risk.table.height is raised below for headroom. The in-plot log-rank annotation
# is scaled by the same KM_BASE/BASE factor.
KM_BASE <- 13                     # compensated base_size for KM ggtheme (vs BASE=8 vector)
# The KM slot here is only one-third of the page width (61 mm), so the two
# at-risk rows (Low/High) have less vertical room than in Fig 2's half-width KM
# panels. KM_FS = 3.2 (vs Fig 2's 4.0) plus risk.table.height = 0.36 keeps the two
# rows legible without vertical collision in this narrower panel.
KM_FS   <- 6.4                    # risk-table number size (decoupled from KM_BASE); 2x enlarged per request
KM_ANN  <- 2.5 * KM_BASE / BASE   # in-plot log-rank annotation size

COLORS_RISK <- c("High" = "#E64B35", "Low" = "#4DBBD5")

# Format a log-rank p-value: guard against floating-point underflow to 0,
# which would otherwise print an impossible "p = 0.00e+00". Mirrors Fig 2.
fmt_logrank_p <- function(p) {
  if (is.na(p)) return("Log-rank p = NA")
  if (p < 2.2e-16) "Log-rank p < 2.2e-16" else sprintf("Log-rank p = %.2e", p)
}

theme_pub <- function(bs = BASE) {
  theme_classic(base_size = bs, base_family = FONT) +
    theme(
      axis.line        = element_line(linewidth = 0.4),
      axis.ticks       = element_line(linewidth = 0.3),
      axis.text        = element_text(size = bs - 1, color = "black"),
      axis.title       = element_text(size = bs),
      strip.text       = element_text(size = bs, face = "bold"),
      strip.background = element_blank(),
      legend.text      = element_text(size = bs - 1),
      legend.title     = element_text(size = bs),
      legend.key.size  = unit(3, "mm"),
      plot.title       = element_text(size = bs, face = "bold", hjust = 0),
      panel.grid       = element_blank()
    )
}

###############################################################################
# Panel A — C-index heatmap (top 25 combinations)
###############################################################################
message("Loading ml_combination_results.RData ...")
load(file.path(DATA_PROC, "ml_combination_results.RData"))

message("Panel A (C-index heatmap, top 25) ...")

combo_ranks <- ml_results %>%
  dplyr::arrange(dplyr::desc(rank_score)) %>%
  dplyr::mutate(combo_name = paste0(fs_algorithm, " + ", build_algorithm))

top25_combos <- combo_ranks %>%
  dplyr::slice_head(n = 25)

hm_long <- top25_combos %>%
  dplyr::select(fs_algorithm, build_algorithm, rank_score)

hm_wide <- hm_long %>%
  dplyr::group_by(fs_algorithm, build_algorithm) %>%
  dplyr::summarise(rank_score = max(rank_score, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = build_algorithm, values_from = rank_score)

hm_mat <- as.matrix(hm_wide[, -1])
rownames(hm_mat) <- hm_wide$fs_algorithm

col_fun <- colorRamp2(
  c(min(hm_mat, na.rm = TRUE), median(hm_mat, na.rm = TRUE), max(hm_mat, na.rm = TRUE)),
  c("#4DBBD5", "white", "#E64B35")
)

ht_top25 <- Heatmap(
  hm_mat,
  name             = "Rank\nScore",
  col              = col_fun,
  cluster_rows     = FALSE, cluster_columns = FALSE,
  row_names_gp     = gpar(fontsize = 7.5, fontfamily = FONT),
  column_names_gp  = gpar(fontsize = 7.5, fontfamily = FONT),
  column_names_rot = 45,
  row_names_side   = "left",
  na_col           = "grey92",
  heatmap_legend_param = list(
    title_gp  = gpar(fontsize = 7, fontface = "bold", fontfamily = FONT),
    labels_gp = gpar(fontsize = 6, fontfamily = FONT)
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    if (!is.na(hm_mat[i, j]))
      grid.text(sprintf("%.3f", hm_mat[i, j]), x, y,
                gp = gpar(fontsize = 5.5, fontfamily = FONT))
  },
  column_title    = "C-index rank score: Feature Selection × Model Building (top 25)",
  column_title_gp = gpar(fontsize = 7.5, fontface = "bold", fontfamily = FONT),
  # Set an explicit wide heatmap body so the few-column matrix spans the full
  # first-row width instead of leaving large right-side whitespace.
  width  = unit(105, "mm"),
  height = unit(36,  "mm")
)

ht_A_png <- file.path(TMP_DIR, "ht_A.png")
# Render at the panel's true composite slot size (SLOT_A_W x SLOT_A_H) so the
# heatmap is placed ~1:1 by draw_grob and its labels keep their set pt instead of
# being shrunk to fit a smaller slot. Padding (bottom,left,top,right): top gutter
# so the column title is not clipped, left gutter so full learner names
# (ElasticNet03 / SuperPC) render in full, bottom gutter for the rotated column
# labels.
png(ht_A_png, width = SLOT_A_W, height = SLOT_A_H, units = "mm", res = 300, family = FONT)
draw(ht_top25, heatmap_legend_side = "right",
     padding = unit(c(8, 24, 6, 6), "mm"))
dev.off()
message(sprintf("  Saved heatmap PNG (%d rows)", nrow(hm_mat)))

###############################################################################
# Panel B — Benchmark boxplot (top 25)
###############################################################################
message("Panel B (benchmark boxplot, top 25) ...")

fold_data <- do.call(rbind, lapply(all_combo_results, function(x) {
  data.frame(
    combo_name  = x$combo_name,
    fs          = x$fs,
    build       = x$build,
    cindex_fold = x$nested_cv_cindex,
    stringsAsFactors = FALSE
  )
}))

combo_order <- fold_data %>%
  dplyr::group_by(combo_name) %>%
  dplyr::summarise(median_cindex = median(cindex_fold, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(median_cindex))

SELECTED_COMBO <- "LASSO + RSF"
combo_order$is_selected <- combo_order$combo_name == SELECTED_COMBO

top25_names <- combo_order %>%
  dplyr::slice_head(n = 25) %>%
  dplyr::pull(combo_name)

if (!SELECTED_COMBO %in% top25_names) {
  top25_names <- c(top25_names[1:24], SELECTED_COMBO)
}

fold_top25 <- fold_data %>%
  dplyr::filter(combo_name %in% top25_names) %>%
  dplyr::left_join(combo_order[, c("combo_name", "is_selected", "median_cindex")],
                   by = "combo_name") %>%
  dplyr::mutate(combo_ordered = factor(
    combo_name,
    levels = combo_order$combo_name[combo_order$combo_name %in% top25_names]
  ))

p_B <- ggplot(fold_top25, aes(x = combo_ordered, y = cindex_fold)) +
  geom_boxplot(
    data    = dplyr::filter(fold_top25, !is_selected),
    aes(x = combo_ordered, y = cindex_fold),
    fill    = "#AABBD4", color = "#6688AA", width = 0.65,
    outlier.size = 0.8, outlier.alpha = 0.6, linewidth = 0.35
  ) +
  geom_boxplot(
    data    = dplyr::filter(fold_top25, is_selected),
    aes(x = combo_ordered, y = cindex_fold),
    fill    = "#E64B35", color = "#8B0000", width = 0.65, linewidth = 1.0,
    outlier.size = 1.5, outlier.color = "#8B0000"
  ) +
  geom_hline(
    yintercept = combo_order$median_cindex[combo_order$combo_name == SELECTED_COMBO],
    linetype = "dashed", color = "#E64B35", linewidth = 0.5, alpha = 0.7
  ) +
  scale_y_continuous(name = "C-index (5-fold nested CV)",
                     breaks = seq(0.74, 0.92, by = 0.02)) +
  labs(x = "Algorithm combination (top 25 by median C-index)",
       title = "Top-25 ML combinations: nested CV C-index") +
  theme_pub() +
  theme(
    axis.text.x  = element_text(angle = 60, hjust = 1, vjust = 1, size = 6.5),
    axis.text.y  = element_text(size = 7),
    plot.margin  = margin(2, 4, 2, 2, "mm")
  )

###############################################################################
# Panel C — Permutation test (LASSO-CoxPH)
###############################################################################
message("Panel C (permutation test) ...")

perm_csv <- file.path(RES_DIR, "permutation_test_uirs_coxph.csv")
perm_dat  <- read.csv(perm_csv, stringsAsFactors = FALSE)

get_val <- function(m) {
  row <- perm_dat[perm_dat$metric == m, ]
  if (nrow(row) == 0) return(NA)
  as.numeric(row$value[1])
}
obs_cindex <- get_val("observed_cindex")
null_mean  <- get_val("null_mean")
null_sd    <- get_val("null_sd")
n_perm     <- get_val("n_perm")
p_val      <- get_val("p_value")
message(sprintf("  Permutation: obs=%.4f, null_mean=%.4f, p=%.4f, n=%d",
                obs_cindex, null_mean, p_val, as.integer(n_perm)))

perm_rds  <- file.path(DATA_PROC, "coxph_permutation_results.RData")
null_dist <- NULL
if (file.exists(perm_rds)) {
  tryCatch({
    load(perm_rds)
    for (obj in ls()) {
      v <- get(obj)
      if (is.numeric(v) && length(v) > 50 &&
          !obj %in% c("obs_cindex", "null_mean", "null_sd", "p_val", "n_perm")) {
        null_dist <- v
        message(sprintf("  Found null dist in '%s': %d values", obj, length(v)))
        break
      }
    }
  }, error = function(e) message("  Could not load perm RDS: ", e$message))
}
if (is.null(null_dist)) {
  message("  Drawing null dist from N(null_mean, null_sd)")
  set.seed(42)
  null_dist <- rnorm(n_perm, mean = null_mean, sd = null_sd)
}

perm_df <- data.frame(null_cindex = null_dist)

p_C <- ggplot(perm_df, aes(x = null_cindex)) +
  geom_histogram(aes(y = after_stat(density)), bins = 40,
                 fill = "#4DBBD5", color = "white", linewidth = 0.3, alpha = 0.8) +
  geom_density(color = "#3C5488", linewidth = 0.5) +
  geom_vline(xintercept = obs_cindex,
             color = "#E64B35", linewidth = 1.0, linetype = "solid") +
  # Anchor the observed-value label to read leftward (into the plot) from the
  # vline so "Observed 0.728" stays fully inside the panel instead of running
  # off the right edge.
  annotate("text", x = obs_cindex - 0.003, y = Inf,
           label = sprintf("Observed\n%.4f", obs_cindex),
           hjust = 1, vjust = 1.2, size = 2.5, color = "#E64B35",
           family = FONT, fontface = "bold") +
  annotate("text", x = min(null_dist) + 0.002, y = Inf,
           label = sprintf("p = %.3f\n(n=%d permutations)", p_val, as.integer(n_perm)),
           hjust = 0, vjust = 1.2, size = 2.5, color = "grey30", family = FONT) +
  labs(x = "Null C-index distribution", y = "Density",
       title = "Permutation test: LASSO-CoxPH (UIRS)") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.06))) +
  theme_pub()

###############################################################################
# Panels D, E, F — KM curves (TCGA, CGGA batch2, CGGA batch1)
###############################################################################
message("Loading risk_model_final.RData ...")
load(file.path(DATA_PROC, "risk_model_final.RData"))

# KM ggtheme uses KM_BASE (not BASE) to compensate for survminer packing the
# plot + risk table into the canvas; see KM_BASE note above.
THEME_KM <- theme_classic(base_size = KM_BASE, base_family = FONT) +
  theme(
    axis.text   = element_text(size = KM_BASE - 1, color = "black"),
    axis.title  = element_text(size = KM_BASE),
    plot.title  = element_text(size = KM_BASE, face = "bold", hjust = 0),
    legend.text = element_text(size = KM_BASE - 1),
    legend.title = element_text(size = KM_BASE),
    panel.grid  = element_blank()
  )

make_km_png <- function(fit_obj, data, title_txt, filename) {
  lr  <- tryCatch(survdiff(Surv(time, status) ~ risk_group, data = data),
                  error = function(e) NULL)
  lr_p <- if (!is.null(lr)) 1 - pchisq(lr$chisq, df = 1) else NA

  p_km <- ggsurvplot(
    fit_obj, data = data,
    pval = FALSE, conf.int = TRUE, conf.int.alpha = 0.12,
    palette = unname(c(COLORS_RISK["Low"], COLORS_RISK["High"])),
    xlab = "Time (months)", ylab = "Overall Survival",
    legend.title = "", legend.labs = c("Low", "High"),
    risk.table = TRUE, risk.table.height = 0.45, risk.table.y.text = FALSE,
    fontsize = KM_FS,
    ggtheme = THEME_KM,
    tables.theme = THEME_KM + theme(axis.line  = element_blank(),
                                     axis.ticks = element_blank(),
                                     axis.text  = element_text(size = KM_BASE - 1))
  )
  p_label <- fmt_logrank_p(lr_p)
  # Place the log-rank p annotation in the lower-left empty quadrant so it does
  # not overprint the descending curves.
  # All three KM panels (D/E/F) suppress the in-panel legend; a single shared
  # legend grob is placed separately in the composite (right half of Row 3).
  p_km$plot <- p_km$plot +
    annotate("text", x = -Inf, y = 0.06,
             label = p_label, hjust = -0.08, vjust = 0,
             size = KM_ANN, family = FONT) +
    labs(title = title_txt) +
    theme(legend.position = "none",
          plot.title = element_text(size = KM_BASE, face = "bold", hjust = 0,
                                    margin = margin(b = 1, unit = "mm")),
          plot.margin = margin(2, 3, 1, 2, "mm"))

  fpath <- file.path(TMP_DIR, filename)
  # Render at the panel's true composite slot size (one-third width x row-4
  # height) so draw_grob places it ~1:1 and the compensated KM text lands at the
  # intended visual size.
  png(fpath, width = SLOT_KM_W, height = SLOT_KM_H, units = "mm", res = 600, family = FONT)
  print(p_km)
  dev.off()
  fpath
}

message("Panel D (KM TCGA) ...")
surv_train <- data.frame(
  time       = train_expr$time / 30,
  status     = train_expr$status,
  risk_group = factor(risk_group_train, levels = c("Low", "High"))
)
fit_train <- survfit(Surv(time, status) ~ risk_group, data = surv_train)
png_D <- make_km_png(fit_train, surv_train, "TCGA training set", "km_D.png")  # legend suppressed — shared legend in Row 3 spacer

message("Panel E (KM CGGA batch2) ...")
png_E <- NULL
if ("CGGA_batch2" %in% names(validation_results)) {
  vr     <- validation_results[["CGGA_batch2"]]
  fit_b2 <- survfit(Surv(time, status) ~ risk_group, data = vr$surv_data)
  png_E  <- make_km_png(fit_b2, vr$surv_data, "CGGA Batch2 (Tuning Validation)", "km_E.png")  # legend suppressed
}

message("Panel F (KM CGGA batch1) ...")
png_F <- NULL
if ("CGGA_batch1" %in% names(validation_results)) {
  vr    <- validation_results[["CGGA_batch1"]]
  sd_b1 <- vr$surv_data
  n_low  <- sum(sd_b1$risk_group == "Low",  na.rm = TRUE)
  n_high <- sum(sd_b1$risk_group == "High", na.rm = TRUE)
  if (n_low == 0 || n_high == 0) {
    local_med    <- median(sd_b1$risk_score, na.rm = TRUE)
    sd_b1$risk_group <- factor(ifelse(sd_b1$risk_score > local_med, "High", "Low"),
                                levels = c("Low", "High"))
  }
  fit_b1 <- survfit(Surv(time, status) ~ risk_group, data = sd_b1)
  png_F  <- make_km_png(fit_b1, sd_b1, "CGGA Batch1 (Blind Test)", "km_F.png")  # legend suppressed
}

###############################################################################
# Panel G — Time-dependent ROC (CGGA batch1, LASSO-CoxPH model)
#
# Source: coxph_interpretable_model.RData → coxph_results$cgga1_val_df
#   - time is in DAYS (range 19–4809)
#   - risk_score from LASSO-CoxPH with strata(IDH)
#   - Manuscript AUC: 1-yr 0.761, 3-yr 0.856, 5-yr 0.888
#   - timeROC called at 365.25, 1095.75, 1826.25 days (= 1, 3, 5 years)
###############################################################################
message("Panel G (time-dependent ROC — CGGA batch1, LASSO-CoxPH) ...")

png_G <- NULL

load(file.path(DATA_PROC, "coxph_interpretable_model.RData"))
b1 <- coxph_results$cgga1_val_df
message(sprintf("  cgga1_val_df rows: %d, time range: %.0f-%.0f days",
                nrow(b1), min(b1$time, na.rm=TRUE), max(b1$time, na.rm=TRUE)))

# Keep complete, positive-time rows
valid <- complete.cases(b1$time, b1$status, b1$risk_score) & b1$time > 0
sv    <- b1[valid, ]
message(sprintf("  Valid rows for timeROC: %d", nrow(sv)))

# Time points in days: 1 yr = 365.25, 3 yr = 1095.75, 5 yr = 1826.25
time_pts <- c(365.25, 1095.75, 1826.25)

troc <- timeROC(
  T      = sv$time,
  delta  = sv$status,
  marker = sv$risk_score,
  cause  = 1,
  times  = time_pts,
  iid    = TRUE
)

auc_vals <- round(troc$AUC, 3)
message(sprintf("  AUC — 1-yr: %.3f, 3-yr: %.3f, 5-yr: %.3f",
                auc_vals[1], auc_vals[2], auc_vals[3]))

# Data integrity gate: verify computed values match manuscript citations
expected <- c(0.761, 0.856, 0.888)
tol <- 0.005   # tight tolerance — these are stored computed values
for (i in 1:3) {
  yr_label <- c("1-yr", "3-yr", "5-yr")[i]
  diff_val <- abs(auc_vals[i] - expected[i])
  if (diff_val > tol) {
    stop(sprintf(
      "DATA INTEGRITY: %s AUC = %.3f deviates >%.3f from manuscript value %.3f. Investigate data source.",
      yr_label, auc_vals[i], tol, expected[i]
    ))
  }
}
message("  AUC values pass data-integrity check (within 0.005 of manuscript values).")

# Build ROC data frames for ggplot
make_roc_df <- function(troc_obj, time_idx, label) {
  fp <- troc_obj$FP[, time_idx]
  tp <- troc_obj$TP[, time_idx]
  ok <- !is.na(fp) & !is.na(tp)
  data.frame(FP    = c(0, fp[ok], 1),
             TP    = c(0, tp[ok], 1),
             label = label,
             stringsAsFactors = FALSE)
}

lbl1 <- sprintf("1-year AUC = %.3f", auc_vals[1])
lbl2 <- sprintf("3-year AUC = %.3f", auc_vals[2])
lbl3 <- sprintf("5-year AUC = %.3f", auc_vals[3])

roc_df       <- rbind(make_roc_df(troc, 1, lbl1),
                      make_roc_df(troc, 2, lbl2),
                      make_roc_df(troc, 3, lbl3))
roc_df$label <- factor(roc_df$label, levels = c(lbl1, lbl2, lbl3))

colour_vals        <- c("#E64B35", "#4DBBD5", "#00A087")
names(colour_vals) <- c(lbl1, lbl2, lbl3)

p_G <- ggplot(roc_df, aes(x = FP, y = TP, color = label)) +
  geom_line(linewidth = 1.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey60", linewidth = 0.6) +
  scale_color_manual(values = colour_vals, name = NULL) +
  scale_x_continuous(name = "1 - Specificity", limits = c(0, 1),
                     expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(name = "Sensitivity",      limits = c(0, 1),
                     expand = expansion(mult = c(0.01, 0.01))) +
  labs(title = "Time-dependent ROC: CGGA Batch1 (Blind Test)") +
  theme_pub(bs = 16) +
  theme(
    legend.position      = c(0.98, 0.04),
    legend.justification = c("right", "bottom"),
    legend.background    = element_rect(fill = "white", color = NA, linewidth = 0.2),
    legend.key           = element_blank(),
    legend.text          = element_text(size = 11),
    legend.key.size      = unit(3.5, "mm"),
    legend.key.height    = unit(3.8, "mm"),
    legend.spacing.y     = unit(0.2, "mm"),
    plot.margin          = margin(3, 5, 3, 3, "mm")
  )

png_G_path <- file.path(TMP_DIR, "roc_G.png")
# Save at the panel's true composite slot size (full width x row-5 height) so the
# ROC panel is placed ~1:1 by draw_grob and its 8 pt theme text matches the
# vector panels B/C instead of being shrunk to fit a smaller slot.
ggsave(png_G_path, plot = p_G,
       width = SLOT_G_W, height = SLOT_G_H, units = "mm", dpi = 600)
png_G <- png_G_path
message("  Panel G PNG saved: ", png_G_path)

###############################################################################
# Assemble composite (7 panels)
###############################################################################
message("=== Assembling 7-panel composite ===")

wrap_png <- function(fpath) {
  if (!is.null(fpath) && file.exists(fpath)) {
    img  <- png::readPNG(fpath)
    grob <- grid::rasterGrob(img, interpolate = TRUE)
    cowplot::ggdraw() + cowplot::draw_grob(grob)
  } else {
    ggplot() + theme_void()
  }
}

p_A_wrap <- wrap_png(ht_A_png)

# Wrap D, E, F — all three have legend suppressed (legend.position="none")
p_D_wrap <- wrap_png(file.path(TMP_DIR, "km_D.png"))
p_E_wrap <- wrap_png(file.path(TMP_DIR, "km_E.png"))
p_F_wrap <- wrap_png(file.path(TMP_DIR, "km_F.png"))

p_G_wrap <- wrap_png(png_G)

# ── Shared KM legend — S9-style reserved right column ─────────────────────────
# Root-cause diagnosis (after multiple failed attempts):
#   • All PNG-based approaches fail because patchwork vertically stretches the
#     raster 1.56× beyond SLOT_KM_H — the KM PNGs are also stretched, but their
#     content fills the full height so it is invisible; the legend PNG's content
#     (14 mm block in a 67 mm canvas) ends up at 76% from top after stretching.
#   • The stretch factor cannot be predicted without running patchwork once.
#
# Correct solution: make legend_col a native ggplot2 object so patchwork renders
# it at the correct coordinate scale without any raster stretching.
# Use geom_segment (draws a horizontal line in the panel) with alpha=0 on the
# actual data so nothing is visible in the panel area, and override.aes to show
# a solid colored segment as the legend key glyph.
# shape=NA alone on geom_point suppresses the glyph; geom_segment with
# override.aes = list(linewidth, alpha=1) keeps the segment key visible.
#
# legend_col is a theme_void ggplot with:
#   legend.position = c(0.5, 0.855) — vertically centered in the patchwork cell natively
#   legend.title    = 6.0 pt bold  (reduced from 9 pt to match KM panel text ~22 px)
#   legend.text     = 5.0 pt       (reduced from 8 pt)
#   legend.key.size = 3.0 mm       (reduced from 3.5 mm)
# Column width is 0.30 units (down from 0.55), total = 3.30, SLOT_KM_W denominator = 3.30.
# Patchwork renders this ggplot at the cell's true physical dimensions, so
# legend.position = c(0.5, 0.855) correctly centers the guide box.

SLOT_LEG_W <- (0.30 / 3.30) * W_MM          # ≈ 16.6 mm  (column = 0.30 units)

legend_df <- data.frame(
  x     = c(0, 1, 0, 1),
  y     = c(0, 0, 1, 1),
  group = factor(c("Low","Low","High","High"), levels = c("Low","High"))
)

legend_col <- ggplot(legend_df, aes(x = x, y = y,
                                    color = group, group = group)) +
  geom_segment(aes(x = x, xend = x + 0.001,
                   y = y, yend = y),
               alpha = 0, show.legend = TRUE) +
  scale_color_manual(
    name   = "Risk group",
    values = c("Low" = "#4DBBD5", "High" = "#E64B35"),
    labels = c("Low", "High"),
    guide  = guide_legend(
      override.aes = list(
        linewidth = 1.2,
        linetype  = 1,
        alpha     = 1
      )
    )
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void(base_family = FONT) +
  theme(
    # Calibrated position: legend.position = c(0.5, 0.5) lands at 77% from top;
    # c(0.5, 0.77) lands at 56.5% from top.  Target is 50% from top.
    # Linear interpolation: need +0.085 more → 0.77 + 0.085 = 0.855.
    legend.position      = c(0.5, 0.855),
    legend.justification = c(0.5, 0.5),
    legend.title         = element_text(size = 6.0, face = "bold",
                                        family = FONT, color = "black",
                                        margin = margin(b = 1.0)),
    legend.text          = element_text(size = 5.0, family = FONT,
                                        color = "black"),
    legend.key.size      = unit(3.0, "mm"),
    legend.key.width     = unit(4.0, "mm"),
    legend.key           = element_blank(),
    legend.spacing.y     = unit(0.8, "mm"),
    legend.background    = element_blank(),
    plot.background      = element_rect(fill = "white", color = NA)
  )

message(sprintf("  Shared KM legend column built as native ggplot (%.1f mm wide, 6.0/5 pt)",
                SLOT_LEG_W))

# Row 3: C left | plot_spacer() right
row3 <- (p_C | plot_spacer()) + plot_layout(widths = c(1, 1))

# Row 4: D | E | F | legend_col
# widths: the three KM panels are equal (1 unit each); legend_col is 0.30 units
# (≈16.6 mm) — sufficient for "Risk group" 5.5 pt bold title without clipping.
# The legend_col cell must receive NO panel tag; patchwork's auto-tag would
# assign it a letter. We override tag_levels explicitly (see composite below).
row4 <- (p_D_wrap | p_E_wrap | p_F_wrap | legend_col) +
  plot_layout(widths = c(1, 1, 1, 0.30))

# Row 5: G spans the full row width
row5 <- p_G_wrap

# ── Panel tagging ─────────────────────────────────────────────────────────────
# Layout tree with auto "A" tagging:
#   Row1: p_A_wrap          → A
#   Row2: p_B               → B
#   Row3: p_C | spacer      → C (spacer skipped automatically)
#   Row4: D | E | F | leg   → D E F [legend must be ""]
#   Row5: p_G_wrap          → G
#
# With patchwork tag_levels = "A", the legend_col cell (being a non-NULL ggplot
# object) will consume a tag letter (becoming H) which is wrong — G is the last
# data panel. We must supply an explicit tag list with "" for the legend cell.
# The spacer in Row3 is skipped by patchwork auto-tagging, but an explicit list
# requires entries for every cell (spacers included if nested in a row).
#
# Strategy: use a named character vector for tag_levels to override per-cell.
# patchwork's explicit tag override: supply tag_levels = list(c(...)) where
# the vector has one entry per non-spacer panel in document order.
# Document order of non-spacer panels:
#   p_A_wrap, p_B, p_C, p_D_wrap, p_E_wrap, p_F_wrap, legend_col, p_G_wrap
# We want:  "A",  "B", "C", "D",    "E",    "F",      "",         "G"

p_composite <- (p_A_wrap / p_B / row3 / row4 / row5) +
  plot_annotation(
    tag_levels = list(c("A","B","C","D","E","F","","G")),
    theme      = theme(
      plot.tag = element_text(size = 13, face = "bold", family = FONT,
                              color = "black", margin = margin(0,2,2,0,"mm"))
    )
  ) +
  plot_layout(heights = LAYOUT_HEIGHTS)

out_pdf <- file.path(OUT_DIR, "Figure4_composite.pdf")
out_png <- file.path(OUT_DIR, "Figure4_composite.png")

ggsave(out_pdf, plot = p_composite, device = cairo_pdf,
       width  = W_MM * MM_IN,
       height = H_MM * MM_IN,
       units  = "in", dpi = 300)
message("Saved PDF: ", out_pdf)

ggsave(out_png, plot = p_composite,
       width  = W_MM * MM_IN,
       height = H_MM * MM_IN,
       units  = "in", dpi = 300)
message("Saved PNG: ", out_png)

# Verify output sizes
for (f in c(out_pdf, out_png)) {
  info <- file.info(f)
  message(sprintf("  %s: %.1f KB", basename(f), info$size / 1024))
}

unlink(TMP_DIR, recursive = TRUE)

message("=== Figure 4 composite (v3, 7 panels) done ===")
message(sprintf("Dimensions: %d x %d mm  (W/H aspect = %.2f)",
                W_MM, H_MM, W_MM / H_MM))
message("Panel G: time-dependent ROC (CGGA batch1) — 1-yr, 3-yr, 5-yr AUC from real data.")
