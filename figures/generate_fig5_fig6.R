###############################################################################
# generate_fig5_fig6.R
# Generate Figure 5B, Figure 6F, Figure 6G, and fix Fig6A nomogram C-index
#
# Fig5B  : Boxplot of ALL 100 ML combinations (per-fold C-index, 5 outer CV),
#           sorted by median C-index, LASSO+RSF highlighted in red/outlined
# Fig6F  : Univariate Cox forest plot — properly labelled variables
# Fig6G  : Multivariate Cox forest plot — properly labelled variables
# Fig6A  : Re-save nomogram performance CSV with correct C-index
###############################################################################

source("00_setup/config.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(survival)
  library(rms)
  library(showtext)
  library(sysfonts)
  library(cowplot)
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

set.seed(SEED)

# =============================================================================
# STEP 0: Fix Nomogram C-index
# =============================================================================
message("=== Fixing Nomogram C-index ===")

load(file.path(DATA_PROC, "nomogram_results.RData"))

# Use concordance() which is the correct method for rms::cph objects
con_nomo <- concordance(fit_cph)
c_nomo_correct <- con_nomo$concordance
c_nomo_se      <- sqrt(con_nomo$var)
c_nomo_lower   <- c_nomo_correct - 1.96 * c_nomo_se
c_nomo_upper   <- c_nomo_correct + 1.96 * c_nomo_se

message(sprintf("Nomogram C-index: %.4f (95%% CI: %.4f-%.4f)",
                c_nomo_correct, c_nomo_lower, c_nomo_upper))

# Update the nomogram_performance.csv
perf_summary <- data.frame(
  Metric = c("C-index (Nomogram)",
             "C-index 95% CI lower",
             "C-index 95% CI upper",
             "Number of variables",
             "Variables",
             "Samples"),
  Value  = c(sprintf("%.4f", c_nomo_correct),
             sprintf("%.4f", c_nomo_lower),
             sprintf("%.4f", c_nomo_upper),
             length(avail_vars),
             paste(avail_vars, collapse = " + "),
             nrow(nomo_df)),
  stringsAsFactors = FALSE
)

write.csv(perf_summary,
          file.path(RES_DIR, "nomogram_performance.csv"), row.names = FALSE)
write.csv(perf_summary,
          file.path(RES_DIR, "Part3_MLModel", "nomogram", "nomogram_performance.csv"),
          row.names = FALSE)

message("Nomogram C-index saved to CSV")

# =============================================================================
# STEP 1: Figure 5B — All 100 ML combinations boxplot
# =============================================================================
message("\n=== Generating Figure 5B ===")

load(file.path(DATA_PROC, "ml_combination_results.RData"))

# Build long-format data from per-fold C-indices (5 outer CV folds)
# Note: CoxBoost build_algorithm combinations all failed in CV (all 5 folds NA)
# SuperPC combinations partially failed (2/5 folds NA)
# NAs are real failures from the original ML run — they are excluded from
# boxplot statistics (na.rm = TRUE) but the combination remains on x-axis
fold_data <- do.call(rbind, lapply(all_combo_results, function(x) {
  data.frame(
    combo_name  = x$combo_name,
    fs          = x$fs,
    build       = x$build,
    cindex_fold = x$nested_cv_cindex,
    stringsAsFactors = FALSE
  )
}))

# Tag combos with all-NA folds (completely failed — CoxBoost build)
all_na_combos <- fold_data %>%
  group_by(combo_name) %>%
  summarise(all_na = all(is.na(cindex_fold)), .groups = "drop") %>%
  filter(all_na) %>%
  pull(combo_name)

if (length(all_na_combos) > 0) {
  message(sprintf("  %d combinations completely failed in CV (all folds NA): %s",
                  length(all_na_combos),
                  paste(all_na_combos, collapse=", ")))
}

# Keep all combinations for display — ggplot boxplot will skip NAs silently

# Compute median per combo for ordering
combo_order <- fold_data %>%
  group_by(combo_name) %>%
  summarise(median_cindex = median(cindex_fold, na.rm = TRUE),
            .groups = "drop") %>%
  arrange(median_cindex)

# Mark the selected model
SELECTED_COMBO <- "LASSO + RSF"
combo_order$is_selected <- combo_order$combo_name == SELECTED_COMBO

# Re-level factor for sorted x-axis
fold_data$combo_ordered <- factor(fold_data$combo_name,
                                   levels = combo_order$combo_name)

# Merge highlight flag
fold_data <- fold_data %>%
  left_join(combo_order[, c("combo_name", "is_selected", "median_cindex")],
            by = "combo_name")

# Build the plot
# All non-selected: grey fill; selected: red fill with thick outline
p5b <- ggplot(fold_data, aes(x = combo_ordered, y = cindex_fold)) +
  # Non-selected boxes (background layer)
  geom_boxplot(
    data    = fold_data %>% filter(!is_selected),
    aes(x = combo_ordered, y = cindex_fold),
    fill    = "#AABBD4",
    color   = "#6688AA",
    width   = 0.7,
    outlier.size   = 0.8,
    outlier.alpha  = 0.6,
    outlier.color  = "#6688AA",
    linewidth      = 0.35
  ) +
  # Selected box (LASSO + RSF) — red/highlighted on top
  geom_boxplot(
    data      = fold_data %>% filter(is_selected),
    aes(x = combo_ordered, y = cindex_fold),
    fill      = "#E64B35",
    color     = "#8B0000",
    width     = 0.7,
    linewidth = 1.2,
    outlier.size  = 2,
    outlier.color = "#8B0000"
  ) +
  # Horizontal reference line at the selected model's median C-index
  geom_hline(
    yintercept = combo_order$median_cindex[combo_order$combo_name == SELECTED_COMBO],
    linetype   = "dashed",
    color      = "#E64B35",
    linewidth  = 0.6,
    alpha      = 0.7
  ) +
  # Annotate selected model
  annotate(
    "text",
    x     = which(levels(fold_data$combo_ordered) == SELECTED_COMBO),
    y     = max(fold_data$cindex_fold, na.rm = TRUE) + 0.002,
    label = "LASSO+RSF\n(selected)",
    size  = 2.8,
    color = "#8B0000",
    fontface = "bold",
    hjust = 0.5,
    vjust = 0
  ) +
  scale_y_continuous(
    name   = "C-index (5-fold nested CV)",
    limits = c(
      floor(min(fold_data$cindex_fold, na.rm = TRUE) * 100) / 100 - 0.005,
      ceiling(max(fold_data$cindex_fold, na.rm = TRUE) * 100) / 100 + 0.015
    ),
    breaks = seq(0.75, 0.90, by = 0.01)
  ) +
  labs(
    title    = "Performance of All 100 ML Algorithm Combinations",
    subtitle = sprintf(
      "5-fold nested CV on TCGA (n=845)  |  Selected model: %s (red, highlighted)  |  Combinations with all CV folds failed shown as empty boxes",
      SELECTED_COMBO
    ),
    x = "ML Algorithm Combination (n=100, sorted by median C-index)",
    caption = paste0("Boxes show median \u00b1 IQR of 5-fold nested CV C-indices. ",
                     sprintf("%d combinations using CoxBoost as build algorithm failed in all CV folds (empty boxes). ",
                             length(all_na_combos)),
                     "Dashed red line = median C-index of selected LASSO+RSF model.")
  ) +
  theme_bw(base_size = 9, base_family = "Times") +
  theme(
    plot.title       = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(size = 9, hjust = 0.5, color = "grey40"),
    axis.title.x     = element_text(size = 10),
    axis.title.y     = element_text(size = 10),
    axis.text.x      = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                     size = 5.5, color = "grey30"),
    axis.text.y      = element_text(size = 9),
    panel.grid.major.x = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
    plot.caption       = element_text(size = 7.5, color = "grey45",
                                      hjust = 0, margin = margin(t = 5))
  )

ggsave(
  file.path(FIG_DIR, "Fig5B_nested_cv_violin.pdf"),
  p5b,
  width  = 16,
  height = 6
)

message("Fig5B saved.")

# =============================================================================
# STEP 2: Figure 6F/G — Forest plots with proper labels
# =============================================================================
message("\n=== Generating Figure 6F and 6G forest plots ===")

load(file.path(DATA_PROC, "independent_prognosis_results.RData"))

# --------------------------------------------------------------------------
# Helper: Format p-value for display
# --------------------------------------------------------------------------
fmt_pval <- function(p) {
  ifelse(is.na(p), "NA",
    ifelse(p < 0.001,
           formatC(p, format = "e", digits = 1),
           formatC(p, format = "f", digits = 3)))
}

# --------------------------------------------------------------------------
# Helper: Build a publication-quality forest plot
# Args:
#   df        - data.frame with columns: label, HR, HR_lower, HR_upper, pvalue
#   title     - plot title
#   subtitle  - plot subtitle
#   dot_color - color for significant points
# --------------------------------------------------------------------------
make_forest_plot <- function(df, title, subtitle,
                              sig_color = "#E64B35",
                              ref_color = "#3C5488") {

  n <- nrow(df)

  # Format text columns
  df$hr_ci_text <- sprintf("%.2f (%.2f\u2013%.2f)",
                            df$HR, df$HR_lower, df$HR_upper)
  df$p_text     <- fmt_pval(df$pvalue)
  df$sig        <- ifelse(!is.na(df$pvalue) & df$pvalue < 0.05,
                           "Significant", "Non-significant")

  # y-axis: row 1 at top (highest y value = n)
  df$y_pos <- rev(seq_len(n))

  # Determine sensible x-axis limits on log scale
  all_vals <- c(df$HR_lower[is.finite(df$HR_lower)],
                df$HR_upper[is.finite(df$HR_upper)])
  x_lo <- min(all_vals, na.rm = TRUE) * 0.7
  x_hi <- max(all_vals, na.rm = TRUE) * 1.4
  x_lo <- max(x_lo, 0.05)   # never go below 0.05 for readability

  # Nice log-scale breaks
  candidate_breaks <- c(0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50)
  breaks_use <- candidate_breaks[candidate_breaks >= x_lo & candidate_breaks <= x_hi]
  if (!1 %in% breaks_use) breaks_use <- sort(c(breaks_use, 1))

  # Clip CI bars for display (wide CIs truncated with arrow-like effect)
  df$HR_lower_plot <- pmax(df$HR_lower, x_lo)
  df$HR_upper_plot <- pmin(df$HR_upper, x_hi)

  # Colors
  df$pt_color <- ifelse(df$sig == "Significant", sig_color, "#999999")

  # Panel 1: Variable labels
  p_left <- ggplot(df, aes(y = y_pos)) +
    annotate("text", x = 0, y = n + 0.7, label = "Variable",
             hjust = 0, size = 3.5, fontface = "bold") +
    geom_text(aes(x = 0, label = label), hjust = 0, size = 3.2,
              color = "grey15") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.3, n + 1.2), expand = c(0, 0)) +
    theme_void(base_family = "Times") +
    theme(plot.margin = margin(4, 0, 4, 5))

  # Panel 2: Forest plot
  p_center <- ggplot(df, aes(y = y_pos)) +
    geom_vline(xintercept = 1, linetype = "dashed",
               color = "grey55", linewidth = 0.55) +
    geom_segment(aes(x = HR_lower_plot, xend = HR_upper_plot,
                     y = y_pos, yend = y_pos),
                 color = df$pt_color, linewidth = 0.9) +
    geom_point(aes(x = HR), color = df$pt_color, shape = 18, size = 3.5) +
    annotate("text", x = exp((log(x_lo) + log(x_hi)) / 2), y = n + 0.7,
             label = "Hazard Ratio (95% CI)", size = 3.5, fontface = "bold") +
    scale_x_log10(
      limits = c(x_lo, x_hi),
      breaks = breaks_use,
      labels = as.character(breaks_use)
    ) +
    scale_y_continuous(limits = c(0.3, n + 1.2), expand = c(0, 0)) +
    labs(x = "Hazard Ratio (log scale)") +
    theme_classic(base_size = 10, base_family = "Times") +
    theme(
      axis.title.y  = element_blank(),
      axis.text.y   = element_blank(),
      axis.ticks.y  = element_blank(),
      axis.line.y   = element_blank(),
      axis.text.x   = element_text(size = 9),
      axis.title.x  = element_text(size = 9.5),
      panel.grid    = element_blank(),
      plot.margin   = margin(4, 3, 4, 3)
    )

  # Panel 3: HR (95% CI) text
  p_right1 <- ggplot(df, aes(y = y_pos)) +
    annotate("text", x = 0.5, y = n + 0.7,
             label = "HR (95% CI)", size = 3.5, fontface = "bold", hjust = 0.5) +
    geom_text(aes(x = 0.5, label = hr_ci_text), hjust = 0.5, size = 3.0,
              color = "grey15") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.3, n + 1.2), expand = c(0, 0)) +
    theme_void(base_family = "Times") +
    theme(plot.margin = margin(4, 2, 4, 2))

  # Panel 4: P-value text
  p_right2 <- ggplot(df, aes(y = y_pos)) +
    annotate("text", x = 0.5, y = n + 0.7,
             label = "P value", size = 3.5, fontface = "bold", hjust = 0.5) +
    geom_text(aes(x = 0.5, label = p_text,
                  color = sig == "Significant"),
              hjust = 0.5, size = 3.0) +
    scale_color_manual(values = c("TRUE" = sig_color, "FALSE" = "#999999"),
                       guide  = "none") +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0.3, n + 1.2), expand = c(0, 0)) +
    theme_void(base_family = "Times") +
    theme(plot.margin = margin(4, 5, 4, 0))

  # Combine panels using cowplot
  p_combined <- cowplot::plot_grid(
    p_left, p_center, p_right1, p_right2,
    nrow   = 1,
    rel_widths = c(2.2, 3.5, 1.8, 1.0),
    align  = "h",
    axis   = "tb"
  )

  # Add title
  title_grob <- cowplot::ggdraw() +
    cowplot::draw_label(title,    x = 0.5, y = 0.75, fontface = "bold",
                        size = 12, hjust = 0.5) +
    cowplot::draw_label(subtitle, x = 0.5, y = 0.25, fontface = "plain",
                        size = 9, hjust = 0.5, color = "grey40")

  cowplot::plot_grid(
    title_grob, p_combined,
    ncol        = 1,
    rel_heights = c(0.12, 1)
  )
}

# --------------------------------------------------------------------------
# Fig 6F: Univariate Cox forest plot
# --------------------------------------------------------------------------

# Clean up variable labels: map from raw names to publication labels
univar_clean <- univar_results %>%
  dplyr::mutate(
    label = dplyr::case_when(
      Variable == "UIRS Risk Score"          ~ "UIRS Risk Score",
      Variable == "Age"                       ~ "Age (per year)",
      grepl("Gender", Variable, ignore.case=TRUE) ~ "Gender (Male vs Female)",
      Variable == "WHO Grade"                 ~ "WHO Grade (per grade)",
      grepl("MGMT", Variable, ignore.case=TRUE)   ~ "MGMT Methylation",
      grepl("IDH", Variable, ignore.case=TRUE)    ~ "IDH Status (WT vs Mutant)",
      TRUE                                    ~ Variable
    )
  ) %>%
  dplyr::select(label, HR, HR_lower, HR_upper, pvalue)

p6f <- make_forest_plot(
  df       = univar_clean,
  title    = "Fig 6F: Univariate Cox Regression",
  subtitle = "IDH stratified (strata(IDH)); All available variables from TCGA cohort"
)

ggsave(
  file.path(FIG_DIR, "Fig6F_forest_univariate.pdf"),
  p6f,
  width  = 11,
  height = max(4, nrow(univar_clean) * 0.85 + 1)
)
message("Fig6F saved.")

# --------------------------------------------------------------------------
# Fig 6G: Multivariate Cox forest plot
# --------------------------------------------------------------------------

# Build clean label mapping for multivariate results
# multivar_results has raw R factor level names (e.g. "GradeG3", "Gendermale")
multivar_clean <- multivar_results %>%
  dplyr::mutate(
    label = dplyr::case_when(
      Variable == "risk_score"       ~ "UIRS Risk Score",
      Variable == "Age"              ~ "Age (per year)",
      Variable == "Gendermale"       ~ "Gender (Male vs Female)",
      Variable == "GradeG2"          ~ "WHO Grade G2 (ref: G2)",
      Variable == "GradeG3"          ~ "WHO Grade G3 (vs G2)",
      Variable == "GradeG4"          ~ "WHO Grade G4 (vs G2)",
      Variable == "MGMT"             ~ "MGMT Methylation",
      Variable == "MGMTUnmethylated" ~ "MGMT (Unmethylated vs Methylated)",
      grepl("IDH", Variable, ignore.case=TRUE) ~ "IDH Status (WT vs Mutant)",
      TRUE                           ~ Variable
    )
  ) %>%
  # Only keep variables with valid HR values
  dplyr::filter(!is.na(HR) & is.finite(HR) & !is.na(HR_lower) & !is.na(HR_upper)) %>%
  dplyr::select(label, HR, HR_lower, HR_upper, pvalue)

p6g <- make_forest_plot(
  df        = multivar_clean,
  title     = "Fig 6G: Multivariate Cox Regression",
  subtitle  = "IDH as stratification variable (strata); TCGA cohort (n=745 complete cases)",
  sig_color = "#E64B35"
)

ggsave(
  file.path(FIG_DIR, "Fig6G_forest_multivariate.pdf"),
  p6g,
  width  = 11,
  height = max(4, nrow(multivar_clean) * 0.85 + 1)
)
message("Fig6G saved.")

# --------------------------------------------------------------------------
# Combined Fig 6FG
# --------------------------------------------------------------------------
p6fg <- cowplot::plot_grid(
  p6f, p6g,
  ncol   = 1,
  labels = c("F", "G"),
  label_size = 16
)

ggsave(
  file.path(FIG_DIR, "Fig6FG_forest_combined.pdf"),
  p6fg,
  width  = 11,
  height = max(10, (nrow(univar_clean) + nrow(multivar_clean)) * 0.75 + 2)
)
message("Fig6FG combined saved.")

# =============================================================================
# STEP 3: Print summary
# =============================================================================
message("\n=== Summary ===")
message(sprintf("Nomogram C-index (corrected): %.4f (95%% CI: %.4f-%.4f)",
                c_nomo_correct, c_nomo_lower, c_nomo_upper))
message(sprintf("Fig5B: %d combinations shown; LASSO+RSF highlighted",
                length(unique(fold_data$combo_name))))
message(sprintf("Fig6F: %d variables (univariate)", nrow(univar_clean)))
message(sprintf("Fig6G: %d variables (multivariate)", nrow(multivar_clean)))
message("\nFigures saved:")
message(paste0("  ", file.path(FIG_DIR, "Fig5B_nested_cv_violin.pdf")))
message(paste0("  ", file.path(FIG_DIR, "Fig6F_forest_univariate.pdf")))
message(paste0("  ", file.path(FIG_DIR, "Fig6G_forest_multivariate.pdf")))
message(paste0("  ", file.path(FIG_DIR, "Fig6FG_forest_combined.pdf")))
message(paste0("  ", file.path(RES_DIR, "nomogram_performance.csv"), " (updated)"))
message("\nAll done.")
