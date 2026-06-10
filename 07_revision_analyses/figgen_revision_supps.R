#!/usr/bin/env Rscript
# Three revision supplementary figures (PDF only):
#  1. GSE16011 external validation KM
#  2. UIRS external permutation null distribution
#  3. IDH-WT model comparison (gene/clinical/combined C-index +/- 95% CI)
set.seed(42)
suppressPackageStartupMessages({
  library(survival); library(survminer); library(ggplot2)
})
PROJ <- getwd()
FIG  <- file.path(PROJ, "figures")
COL2 <- c(High = "#E64B35", Low = "#4DBBD5")

## ---- FIG 1: GSE16011 KM ----
g <- readRDS(file.path(PROJ, "data/processed/gse16011_uirs_validation.rds"))
g$grp <- factor(ifelse(g$risk > median(g$risk), "High", "Low"), levels = c("Low","High"))
sf <- survfit(Surv(time, status) ~ grp, data = g)
sd_ <- survdiff(Surv(time, status) ~ grp, data = g)
lp  <- 1 - pchisq(sd_$chisq, length(sd_$n) - 1)
cidx <- tryCatch(concordance(coxph(Surv(time, status) ~ risk, data = g))$concordance, error=function(e) NA)
p1 <- ggsurvplot(sf, data = g, palette = unname(COL2[c("Low","High")]),
                 risk.table = TRUE, pval = FALSE, conf.int = FALSE,
                 legend.labs = c("Low UIRS","High UIRS"), legend.title = "",
                 xlab = "Time (years)", ylab = "Overall survival",
                 title = "GSE16011 external validation (Affymetrix microarray, n=263)")
plab <- if (lp < 1e-12) "log-rank p < 1e-12" else sprintf("log-rank p = %.2g", lp)
p1$plot <- p1$plot + annotate("text", x = 0, y = 0.05, hjust = 0,
            label = sprintf("%s\nC-index = %.3f", plab, cidx), size = 4)
pdf(file.path(FIG, "SuppFig_gse16011_validation_km.pdf"), width = 7, height = 7)
print(p1, newpage = FALSE); dev.off()
cat(sprintf("FIG1 GSE16011 KM: n=%d events=%d logrank_p=%.3g C=%.3f\n",
            nrow(g), sum(g$status), lp, cidx))

## ---- FIG 2: permutation null ----
pm <- readRDS(file.path(PROJ, "data/processed/permutation_test_uirs_coxph.rds"))
nulls <- pm$perm[!is.na(pm$perm)]; obs <- pm$observed; pv <- pm$pval
dfp <- data.frame(c = nulls)
p2 <- ggplot(dfp, aes(c)) +
  geom_histogram(bins = 40, fill = "grey75", colour = "white") +
  geom_vline(xintercept = obs, colour = "#E64B35", linetype = "dashed", linewidth = 1) +
  annotate("text", x = obs, y = Inf, vjust = 2, hjust = 1.05, colour = "#E64B35",
           fontface = "bold", size = 4.2,
           label = sprintf("Observed = %.3f\np = %.3f", obs, pv)) +
  labs(x = "C-index on CGGA-batch1 (external)", y = "Permutations",
       title = sprintf("Permutation test: UIRS (LASSO-CoxPH), external (n=%d)", length(nulls))) +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank())
ggsave(file.path(FIG, "SuppFig_permutation_uirs_coxph.pdf"), p2, width = 7.5, height = 5)
cat(sprintf("FIG2 permutation: obs=%.4f null_mean=%.4f n=%d p=%.4f\n",
            obs, mean(nulls), length(nulls), pv))

## ---- FIG 3: IDH-WT model comparison ----
d <- read.csv(file.path(PROJ, "results/idhwt_model_improvement.csv"),
              stringsAsFactors = FALSE, check.names = FALSE)
names(d) <- tolower(gsub("[^a-z_]", "", tolower(names(d))))
# expected: model, cindex, ci_low, ci_high
mdl <- d[[grep("model", names(d))[1]]]
ci  <- as.numeric(d[[grep("^cindex|c_index|cindex", names(d))[1]]])
lo  <- as.numeric(d[[grep("low", names(d))[1]]])
hi  <- as.numeric(d[[grep("high", names(d))[1]]])
keep <- !grepl("delta", mdl, ignore.case = TRUE)
dd <- data.frame(model = mdl[keep], cindex = ci[keep], lo = lo[keep], hi = hi[keep])
dd$model <- factor(dd$model, levels = dd$model[order(dd$cindex)])
delta_row <- which(grepl("delta", mdl, ignore.case = TRUE))
dlt <- if (length(delta_row)) sprintf("Combined - Gene-only delta = %.3f (95%% CI %.3f, %.3f)",
                                      ci[delta_row[1]], lo[delta_row[1]], hi[delta_row[1]]) else ""
p3 <- ggplot(dd, aes(cindex, model)) +
  geom_vline(xintercept = 0.5, linetype = "dotted", colour = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi), colour = "#3C5488", linewidth = 0.7, size = 0.6) +
  geom_text(aes(label = sprintf("%.3f", cindex)), vjust = -0.9, size = 3.5) +
  labs(x = "C-index (95% CI), CGGA-batch1 IDH-WT external (n=137)", y = "",
       title = "IDH-WT glioma: UIRS vs clinical features",
       subtitle = dlt) +
  coord_cartesian(xlim = c(0.45, 0.75)) +
  theme_bw(base_size = 12) + theme(panel.grid.minor = element_blank())
ggsave(file.path(FIG, "SuppFig_idhwt_model_comparison.pdf"), p3, width = 7.5, height = 4)
cat(sprintf("FIG3 IDH-WT: %s | %s\n",
            paste(sprintf("%s=%.3f", dd$model, dd$cindex), collapse="; "), dlt))

cat("ALL_3_FIGS_DONE\n")
