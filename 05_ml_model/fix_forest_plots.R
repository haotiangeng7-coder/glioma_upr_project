###############################################################################
# fix_forest_plots.R
# Generate publication-quality forest plots for Fig6F (univariate) and
# Fig6G (multivariate) with numeric HR, 95% CI, and p-value text columns.
# IDH is used as strata (not covariate) due to PH assumption violation.
#
# Layout: Variable labels | Forest plot (HR + CI bar) | HR (95% CI) | P value
# Color: significant (p<0.05) red, non-significant gray
###############################################################################

source("00_setup/config.R")
library(ggplot2)
library(grid)
library(gridExtra)
library(gtable)
library(dplyr)

set.seed(SEED)

# Load saved Cox regression results
load(file.path(DATA_PROC, "independent_prognosis_results.RData"))

# =============================================================================
# Helper functions
# =============================================================================
format_pval <- function(p) {
  ifelse(p < 0.001, formatC(p, format = "e", digits = 1),
         formatC(p, format = "f", digits = 3))
}

# Build a single forest plot figure with table columns
# df must have: label, HR, HR_lower, HR_upper, pvalue
make_forest_plot <- function(df, title_text, subtitle_text,
                             clip_lower = NULL, clip_upper = NULL) {

  n <- nrow(df)
  df$y <- rev(seq_len(n))  # bottom-to-top: first row at top
  df$hr_ci <- sprintf("%.2f (%.2f-%.2f)", df$HR, df$HR_lower, df$HR_upper)
  df$p_text <- format_pval(df$pvalue)
  df$sig <- ifelse(df$pvalue < 0.05, "Significant", "Non-significant")
  df$color <- ifelse(df$pvalue < 0.05, "#E64B35", "#999999")

  # Determine x-axis range on log scale
  all_vals <- c(df$HR_lower, df$HR_upper)
  if (is.null(clip_lower)) clip_lower <- min(all_vals) * 0.6
  if (is.null(clip_upper)) clip_upper <- max(all_vals) * 1.5

  # Clip CIs for display (arrows would be ideal but we just clip)
  df$HR_lower_clipped <- pmax(df$HR_lower, clip_lower)
  df$HR_upper_clipped <- pmin(df$HR_upper, clip_upper)

  # Determine nice log-scale breaks
  log_range <- log10(c(clip_lower, clip_upper))
  possible_breaks <- c(0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50)
  breaks <- possible_breaks[possible_breaks >= clip_lower & possible_breaks <= clip_upper]
  if (!1 %in% breaks) breaks <- sort(c(breaks, 1))

  # --- Left panel: Variable labels ---
  p_left <- ggplot(df, aes(y = y)) +
    geom_text(aes(x = 0, label = label), hjust = 0, size = 3.5,
              fontface = "plain", family = "sans") +
    scale_y_continuous(limits = c(0.5, n + 1.5), expand = c(0, 0)) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    # Header
    annotate("text", x = 0, y = n + 0.8, label = "Variable",
             hjust = 0, size = 3.8, fontface = "bold", family = "sans") +
    theme_void() +
    theme(plot.margin = margin(5, 0, 5, 5))

  # --- Center panel: Forest plot ---
  p_center <- ggplot(df, aes(y = y)) +
    # Reference line at HR = 1
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    # CI bars
    geom_segment(aes(x = HR_lower_clipped, xend = HR_upper_clipped,
                     y = y, yend = y, color = sig),
                 linewidth = 0.8, show.legend = FALSE) +
    # HR point
    geom_point(aes(x = HR, color = sig), size = 3, shape = 15,
               show.legend = FALSE) +
    scale_color_manual(values = c("Significant" = "#E64B35",
                                  "Non-significant" = "#999999")) +
    scale_x_log10(breaks = breaks, labels = as.character(breaks),
                  limits = c(clip_lower, clip_upper)) +
    scale_y_continuous(limits = c(0.5, n + 1.5), expand = c(0, 0)) +
    # Header
    annotate("text", x = 1, y = n + 0.8, label = "Hazard Ratio",
             size = 3.8, fontface = "bold", family = "sans") +
    labs(x = "Hazard Ratio (log scale)") +
    theme_classic(base_size = 11) +
    theme(
      axis.title.y = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y  = element_blank(),
      axis.text.x  = element_text(size = 9),
      axis.title.x = element_text(size = 10, face = "plain"),
      panel.grid   = element_blank(),
      plot.margin  = margin(5, 2, 5, 2)
    )

  # --- Right panel 1: HR (95% CI) text ---
  p_hrci <- ggplot(df, aes(y = y)) +
    geom_text(aes(x = 0, label = hr_ci, color = sig),
              hjust = 0.5, size = 3.2, family = "sans", show.legend = FALSE) +
    scale_color_manual(values = c("Significant" = "#E64B35",
                                  "Non-significant" = "#999999")) +
    scale_y_continuous(limits = c(0.5, n + 1.5), expand = c(0, 0)) +
    scale_x_continuous(limits = c(-0.5, 0.5), expand = c(0, 0)) +
    annotate("text", x = 0, y = n + 0.8, label = "HR (95% CI)",
             hjust = 0.5, size = 3.8, fontface = "bold", family = "sans") +
    theme_void() +
    theme(plot.margin = margin(5, 2, 5, 2))

  # --- Right panel 2: P value text ---
  p_pval <- ggplot(df, aes(y = y)) +
    geom_text(aes(x = 0, label = p_text, color = sig),
              hjust = 0.5, size = 3.2, family = "sans", show.legend = FALSE) +
    scale_color_manual(values = c("Significant" = "#E64B35",
                                  "Non-significant" = "#999999")) +
    scale_y_continuous(limits = c(0.5, n + 1.5), expand = c(0, 0)) +
    scale_x_continuous(limits = c(-0.5, 0.5), expand = c(0, 0)) +
    annotate("text", x = 0, y = n + 0.8, label = "P value",
             hjust = 0.5, size = 3.8, fontface = "bold", family = "sans") +
    theme_void() +
    theme(plot.margin = margin(5, 5, 5, 2))

  # --- Combine panels ---
  # Convert to grobs
  g_left   <- ggplotGrob(p_left)
  g_center <- ggplotGrob(p_center)
  g_hrci   <- ggplotGrob(p_hrci)
  g_pval   <- ggplotGrob(p_pval)

  # Title grob
  title_grob <- textGrob(title_text,
                          gp = gpar(fontsize = 13, fontface = "bold",
                                    fontfamily = "sans"),
                          just = "center")
  subtitle_grob <- textGrob(subtitle_text,
                              gp = gpar(fontsize = 10, fontface = "italic",
                                        fontfamily = "sans", col = "gray40"),
                              just = "center")

  # Arrange: widths ~ label(3) | forest(4) | HR text(2.5) | p text(1.5)
  combined <- arrangeGrob(
    g_left, g_center, g_hrci, g_pval,
    nrow = 1,
    widths = unit(c(3.0, 4.0, 2.5, 1.5), "null")
  )

  # Stack title + subtitle + body
  final <- arrangeGrob(
    title_grob,
    subtitle_grob,
    combined,
    nrow = 3,
    heights = unit(c(0.8, 0.5, 5), "null")
  )

  final
}


# =============================================================================
# Fig 6F: Univariate Forest Plot
# =============================================================================
message("=== Generating Fig6F: Univariate Forest Plot ===")

uni_df <- univar_results %>%
  transmute(
    label    = Variable,
    HR       = HR,
    HR_lower = HR_lower,
    HR_upper = HR_upper,
    pvalue   = pvalue
  )

g_uni <- make_forest_plot(
  df            = uni_df,
  title_text    = "Univariate Cox Regression",
  subtitle_text = "IDH-stratified for non-IDH variables; IDH shown unstratified (PH violated)",
  clip_lower    = 0.5,
  clip_upper    = 20
)

pdf_f <- file.path(FIG_DIR, "Fig6F_forest_univariate.pdf")
ggsave(pdf_f, g_uni, width = 11, height = 5)
message("  Saved: ", pdf_f)


# =============================================================================
# Fig 6G: Multivariate Forest Plot
# =============================================================================
message("=== Generating Fig6G: Multivariate Forest Plot ===")

# Clean up multivariate variable labels
label_map <- c(
  "risk_score"       = "UIRS Risk Score",
  "Age"              = "Age",
  "Gendermale"       = "Gender (Male vs Female)",
  "GradeG3"          = "WHO Grade (G3 vs G2)",
  "GradeG4"          = "WHO Grade (G4 vs G2)",
  "MGMTUnmethylated" = "MGMT (Unmethylated vs Methylated)"
)

multi_df_plot <- multivar_results %>%
  transmute(
    label    = ifelse(Variable %in% names(label_map),
                      label_map[Variable], Variable),
    HR       = HR,
    HR_lower = HR_lower,
    HR_upper = HR_upper,
    pvalue   = pvalue
  )

g_multi <- make_forest_plot(
  df            = multi_df_plot,
  title_text    = "Multivariate Cox Regression",
  subtitle_text = "IDH used as stratification variable (strata) due to PH assumption violation",
  clip_lower    = 0.15,
  clip_upper    = 5
)

pdf_g <- file.path(FIG_DIR, "Fig6G_forest_multivariate.pdf")
ggsave(pdf_g, g_multi, width = 11, height = 5)
message("  Saved: ", pdf_g)


# =============================================================================
# Copy to all required locations
# =============================================================================
message("=== Copying to submission and Main_Figures directories ===")

fig6_main <- file.path(FIG_DIR, "Main_Figures", "Figure6_nomogram_prognosis")
fig6_sub  <- file.path(PROJECT_DIR, "manuscript", "submission", "Figures", "Figure6")
dir.create(fig6_main, recursive = TRUE, showWarnings = FALSE)
dir.create(fig6_sub,  recursive = TRUE, showWarnings = FALSE)

for (fname in c("Fig6F_forest_univariate.pdf", "Fig6G_forest_multivariate.pdf")) {
  src <- file.path(FIG_DIR, fname)
  file.copy(src, file.path(fig6_main, fname), overwrite = TRUE)
  file.copy(src, file.path(fig6_sub,  fname), overwrite = TRUE)
  message(sprintf("  Copied %s", fname))
}

message("\n=== Forest plot generation completed ===")
