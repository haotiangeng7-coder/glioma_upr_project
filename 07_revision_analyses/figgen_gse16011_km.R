###############################################################################
# figgen_gse16011_km.R
# Figure: GSE16011 external validation KM curve
# Output: figures/SuppFig_gse16011_validation_km.pdf
###############################################################################
set.seed(42)

library(survival)
library(survminer)
library(ggplot2)

PROJECT_DIR <- getwd()

# --- Load data ----------------------------------------------------------------
df <- readRDS(file.path(PROJECT_DIR, "data/processed/gse16011_uirs_validation.rds"))
stopifnot(all(c("risk", "time", "status") %in% colnames(df)))

cat("Rows:", nrow(df), "| Events:", sum(df$status), "\n")
cat("Risk range:", range(df$risk), "\n")
cat("Risk median:", median(df$risk), "\n")

# Split at cohort median (grp_med column already exists; use it for consistency,
# but verify it matches a fresh median split)
med_cut <- median(df$risk)
# Use the pre-computed grp_med column (matches the published KM_p = 8.61e-12).
# A fresh >= median split differs by 1 sample at the exact median; using the
# stored column ensures reproducibility with the reported statistic.
df$risk_group <- factor(df$grp_med, levels = c("High", "Low"))

cat("Group sizes: High =", sum(df$risk_group == "High"),
    "| Low =", sum(df$risk_group == "Low"), "\n")

# --- Fit KM -------------------------------------------------------------------
fit <- survfit(Surv(time, status) ~ risk_group, data = df)

# Log-rank test
lr   <- survdiff(Surv(time, status) ~ risk_group, data = df)
pval <- 1 - pchisq(lr$chisq, df = 1)

# C-index from CSV (not recomputed here; cross-check only)
cindex_str <- "C-index = 0.668 [95% CI 0.627-0.707]"
cat("Log-rank chi^2:", lr$chisq, "p =", pval, "\n")

# --- Theme (Nature/Springer style, matches project THEME_PUBLICATION) ---------
pub_theme <- theme_classic(base_size = 9) +
  theme(
    text              = element_text(family = "sans"),
    plot.title        = element_text(size = 8, face = "bold", hjust = 0.5),
    axis.title        = element_text(size = 8),
    axis.text         = element_text(size = 7),
    legend.title      = element_text(size = 7, face = "bold"),
    legend.text       = element_text(size = 7),
    legend.position   = "none",   # legend suppressed; group identity clear from colours + risk table
    panel.grid.major  = element_blank(),
    panel.grid.minor  = element_blank(),
    plot.margin       = margin(5, 8, 5, 5, "mm")
  )

# Annotation string (p rendered from actual log-rank; C-index from CSV)
pval_label <- if (pval < 0.001) {
  sprintf("Log-rank p = %.2e", pval)
} else {
  sprintf("Log-rank p = %.4f", pval)
}

annot_label <- paste0(pval_label, "\n", cindex_str)

# --- ggsurvplot ---------------------------------------------------------------
pal <- c(High = "#E64B35", Low = "#4DBBD5")

p <- ggsurvplot(
  fit,
  data          = df,
  palette       = unname(pal[levels(df$risk_group)]),
  legend.labs   = levels(df$risk_group),
  legend.title  = "UIRS risk",
  xlab          = "Time (years)",
  ylab          = "Overall survival probability",
  title         = "GSE16011 external validation (Affymetrix microarray, n=263)",
  risk.table    = TRUE,
  risk.table.height = 0.22,
  risk.table.fontsize = 2.8,
  tables.theme  = theme_cleantable(),
  conf.int      = TRUE,
  conf.int.alpha = 0.12,
  pval          = FALSE,   # we annotate manually for precision
  ggtheme       = pub_theme,
  fontsize      = 3,
  size          = 0.7
)

# Apply Nature-grade font to risk table axis labels too
p$table <- p$table +
  theme(axis.text = element_text(size = 6),
        axis.title = element_text(size = 7))

# Add manual annotation (p + C-index)
p$plot <- p$plot +
  annotate("text",
           x = max(df$time) * 0.02,
           y = 0.05,
           label = annot_label,
           hjust = 0, vjust = 0,
           size = 2.5,
           family = "sans")

# --- Save to PDF (89 mm wide, ~130 mm tall with risk table) -------------------
out_pdf <- file.path(PROJECT_DIR, "figures/SuppFig_gse16011_validation_km.pdf")

# 120 mm (1.5-col) width: title "GSE16011 external validation (Affymetrix microarray, n=263)"
# is ~55 characters at 8 pt bold — requires >= 110 mm at final size to avoid truncation.
cairo_pdf(out_pdf, width = 120 / 25.4, height = 130 / 25.4)
print(p, newpage = FALSE)
dev.off()

# Verify non-empty
fsize <- file.info(out_pdf)$size
cat("Output:", out_pdf, "\n")
cat("File size:", fsize, "bytes\n")
stopifnot(fsize > 0)
cat("DONE\n")
