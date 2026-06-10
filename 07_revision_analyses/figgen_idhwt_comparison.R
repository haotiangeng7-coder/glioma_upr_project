###############################################################################
# figgen_idhwt_comparison.R
# Figure: IDH-WT model comparison (C-index point-range, 3 models)
# Output: figures/SuppFig_idhwt_model_comparison.pdf
###############################################################################
set.seed(42)

library(ggplot2)

PROJECT_DIR <- getwd()

# --- Load data ----------------------------------------------------------------
csv_path <- file.path(PROJECT_DIR, "results/idhwt_model_improvement.csv")
raw <- read.csv(csv_path, stringsAsFactors = FALSE)

cat("Raw data:\n")
print(raw)

# Exclude the Delta row from the C-index panel
plot_df <- subset(raw, model != "Delta(combined-gene)")
plot_df$cindex   <- as.numeric(plot_df$cindex)
plot_df$ci_low   <- as.numeric(plot_df$ci_low)
plot_df$ci_high  <- as.numeric(plot_df$ci_high)

# Ordered bottom-to-top for horizontal forest: gene-only at top for emphasis
plot_df$model <- factor(plot_df$model,
  levels = c(
    "Clinical-only (Age+Grade+MGMT)",
    "Combined",
    "Gene-only (UIRS)"
  )
)

cat("\nPlot data:\n")
print(plot_df)

# Delta annotation values (from raw CSV)
delta_row <- subset(raw, model == "Delta(combined-gene)")
delta_val  <- as.numeric(delta_row$cindex)
delta_low  <- as.numeric(delta_row$ci_low)
delta_high <- as.numeric(delta_row$ci_high)
delta_label <- sprintf("Δ combined−gene = +%.3f\n[95%% CI %.3f, %.3f]",
                       delta_val, delta_low, delta_high)
cat("\nDelta annotation:", delta_label, "\n")

# Colour: Gene-only highlighted in project red; others in greys
model_colours <- c(
  "Gene-only (UIRS)"                 = "#E64B35",
  "Combined"                         = "#555555",
  "Clinical-only (Age+Grade+MGMT)"   = "#888888"
)

# --- Theme --------------------------------------------------------------------
pub_theme <- theme_classic(base_size = 9) +
  theme(
    text              = element_text(family = "sans"),
    plot.title        = element_text(size = 8, face = "bold", hjust = 0.5),
    axis.title        = element_text(size = 8),
    axis.text.y       = element_text(size = 7),
    axis.text.x       = element_text(size = 7),
    legend.position   = "none",
    panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    plot.margin        = margin(5, 12, 5, 5, "mm")
  )

# --- Build plot ---------------------------------------------------------------
# x range: give a bit of room around 0.5 reference line and CI whiskers
x_lo <- min(plot_df$ci_low)  - 0.02
x_hi <- max(plot_df$ci_high) + 0.06   # extra room for delta annotation

p <- ggplot(plot_df, aes(x = cindex, y = model, colour = model)) +
  # Reference line at 0.5 (chance)
  geom_vline(xintercept = 0.5,
             linetype = "solid",
             colour = "grey50",
             linewidth = 0.5) +
  # CI whiskers
  geom_errorbarh(aes(xmin = ci_low, xmax = ci_high),
                 height = 0.18,
                 linewidth = 0.65) +
  # Point estimates
  geom_point(size = 3, shape = 18) +
  # Delta annotation in top-right
  annotate("text",
           x = x_hi - 0.003,
           y = 0.55,
           label = delta_label,
           hjust = 1, vjust = 0,
           size = 2.4,
           colour = "grey30",
           family = "sans") +
  scale_colour_manual(values = model_colours) +
  scale_x_continuous(
    name   = "C-index (95% CI)",
    limits = c(x_lo, x_hi),
    breaks = seq(0.45, 0.75, by = 0.05),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_discrete(
    name   = NULL,
    expand = expansion(add = c(0.6, 0.6))
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "IDH-WT glioma: UIRS vs clinical features\n(CGGA-batch1 external, n=137)"
  ) +
  pub_theme

# --- Save PDF (120 mm wide, 70 mm tall) ---------------------------------------
out_pdf <- file.path(PROJECT_DIR, "figures/SuppFig_idhwt_model_comparison.pdf")

cairo_pdf(out_pdf, width = 120 / 25.4, height = 70 / 25.4)
print(p)
dev.off()

# Verify non-empty
fsize <- file.info(out_pdf)$size
cat("Output:", out_pdf, "\n")
cat("File size:", fsize, "bytes\n")
stopifnot(fsize > 0)
cat("DONE\n")
