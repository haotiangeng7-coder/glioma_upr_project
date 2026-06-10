###############################################################################
# figgen_permutation_uirs.R
# Figure: UIRS permutation null distribution histogram
# Output: figures/SuppFig_permutation_uirs_coxph.pdf
###############################################################################
set.seed(42)

library(ggplot2)

PROJECT_DIR <- getwd()

# --- Load data ----------------------------------------------------------------
obj <- readRDS(file.path(PROJECT_DIR,
               "data/processed/permutation_test_uirs_coxph.rds"))

stopifnot(all(c("observed", "perm", "pval") %in% names(obj)))

obs   <- obj$observed   # 0.7279
perm  <- obj$perm       # 1000 null C-indices
pval  <- obj$pval       # 0.015

cat("Observed C-index:", obs, "\n")
cat("Null mean:", mean(perm), "| Null SD:", sd(perm), "\n")
cat("N permutations:", length(perm), "\n")
cat("p-value:", pval, "\n")

perm_df <- data.frame(cindex = perm)

# --- Theme (Nature style, 89 mm single-column) --------------------------------
pub_theme <- theme_classic(base_size = 9) +
  theme(
    text              = element_text(family = "sans"),
    plot.title        = element_text(size = 8, face = "bold", hjust = 0.5),
    axis.title        = element_text(size = 8),
    axis.text         = element_text(size = 7),
    legend.position   = "none",
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    plot.margin       = margin(5, 10, 5, 5, "mm")
  )

# --- Build plot ---------------------------------------------------------------
# Sensible binwidth ~0.02 for 1000 null C-indices spanning ~0.27 range
bw <- 0.02

p <- ggplot(perm_df, aes(x = cindex)) +
  geom_histogram(binwidth = bw,
                 fill = "grey65",
                 colour = "white",
                 linewidth = 0.3) +
  geom_vline(xintercept = obs,
             colour = "#D62728",
             linetype = "dashed",
             linewidth = 0.9) +
  # Annotation placed LEFT of the observed line (hjust=1) to stay within canvas.
  # The observed line at 0.7279 is near the right edge; annotating to its left
  # avoids right-margin clipping on a 120 mm canvas.
  annotate("text",
           x = obs - 0.005,
           y = Inf,
           label = sprintf("Observed\nC-index = %.4f\np = %.3f", obs, pval),
           hjust = 1, vjust = 1.25,
           size = 2.5,
           colour = "#D62728",
           family = "sans") +
  scale_x_continuous(
    name   = "C-index on CGGA-batch1 (external)",
    expand = expansion(mult = c(0.04, 0.08)),   # extra right expand for last bar
    breaks = scales::pretty_breaks(n = 5)
  ) +
  scale_y_continuous(
    name   = "Count (permutations)",
    expand = expansion(mult = c(0, 0.12))        # head room so Inf annotation is visible
  ) +
  coord_cartesian(clip = "off") +
  labs(
    title = "Permutation test — UIRS (LASSO->CoxPH), external validation"
  ) +
  pub_theme

# --- Save PDF (89 mm wide, 80 mm tall) ----------------------------------------
out_pdf <- file.path(PROJECT_DIR, "figures/SuppFig_permutation_uirs_coxph.pdf")

# 120 mm (1.5-col): title "Permutation test — UIRS (LASSO->CoxPH), external validation"
# is ~54 characters at 8 pt bold — needs >= 110 mm to avoid truncation.
cairo_pdf(out_pdf, width = 120 / 25.4, height = 85 / 25.4)
print(p)
dev.off()

# Verify non-empty
fsize <- file.info(out_pdf)$size
cat("Output:", out_pdf, "\n")
cat("File size:", fsize, "bytes\n")
stopifnot(fsize > 0)
cat("DONE\n")
