#!/usr/bin/env Rscript
# ============================================================================
# Fork C: External prognostic validation of the LOCKED UIRS (LASSO->CoxPH, 15
# genes) on GSE16011 (Gravendeel 2009, Affymetrix microarray), using survival
# recovered from the paper's Supplementary Table 1.
#
# UIRS is applied as a LOCKED linear risk score: risk = sum(beta_g * expr_g),
# with the 15 frozen betas and the frozen TCGA-median cutpoint (no re-fitting,
# no GSE-derived cutoff for the primary analysis).
#
# CAVEAT (reported, not hidden): GSE16011 is Affymetrix microarray; UIRS was
# trained on TCGA RNA-seq log2(TPM+1). C-index is rank-based (scale-invariant);
# KM by the TCGA cutoff is scale-sensitive, so a GSE-median split is shown too.
# ============================================================================
set.seed(42)
suppressPackageStartupMessages({
  library(survival); library(timeROC); library(org.Hs.eg.db); library(AnnotationDbi)
})
PROJ <- getwd()

# ---- 1. Locked UIRS ----
load(file.path(PROJ, "data/processed/coxph_interpretable_model.RData"))  # coxph_results
betas  <- coxph_results$betas                       # 15 named gene coefficients
cutoff <- coxph_results$train_median_cutoff          # frozen TCGA-median risk cutpoint
cat(sprintf("UIRS: %d genes, TCGA cutoff=%.4f\n", length(betas), cutoff))

# ---- 2. GSE16011 expression: Entrez -> symbol ----
load(file.path(PROJ, "data/processed/gse16011_data.RData"))  # expr_16011_gene (Entrez x GSM)
ent <- rownames(expr_16011_gene)
map <- AnnotationDbi::select(org.Hs.eg.db, keys = ent, columns = "SYMBOL", keytype = "ENTREZID")
map <- map[!is.na(map$SYMBOL), ]
need <- names(betas)
# for each needed gene, pick the Entrez id with highest mean expression
expr_sym <- matrix(NA_real_, nrow = length(need), ncol = ncol(expr_16011_gene),
                   dimnames = list(need, colnames(expr_16011_gene)))
for (g in need) {
  es <- map$ENTREZID[map$SYMBOL == g]
  es <- es[es %in% rownames(expr_16011_gene)]
  if (length(es) == 0) next
  if (length(es) == 1) { expr_sym[g, ] <- expr_16011_gene[es, ] }
  else {
    mu <- rowMeans(expr_16011_gene[es, , drop = FALSE], na.rm = TRUE)
    expr_sym[g, ] <- expr_16011_gene[es[which.max(mu)], ]
  }
}
found <- need[rowSums(is.na(expr_sym)) < ncol(expr_sym)]
cat(sprintf("UIRS genes found in GSE16011: %d/%d (missing: %s)\n",
            length(found), length(need), paste(setdiff(need, found), collapse = ", ")))
cat(sprintf("GSE16011 expression range: [%.2f, %.2f] (microarray; likely already log2)\n",
            min(expr_16011_gene, na.rm = TRUE), max(expr_16011_gene, na.rm = TRUE)))

# ---- 3. UIRS linear risk score (locked betas, genes found) ----
X <- t(expr_sym[found, , drop = FALSE])            # samples x genes
risk <- as.numeric(X %*% betas[found])
names(risk) <- rownames(X)

# ---- 4. Attach recovered survival (matched by GSM) ----
surv <- read.csv(file.path(PROJ, "data/processed/gse16011_survival_from_supp.csv"),
                 stringsAsFactors = FALSE)
surv <- surv[surv$note == "ok" & surv$gsm != "", ]   # glioma with usable survival
rownames(surv) <- surv$gsm
common <- intersect(names(risk), rownames(surv))
df <- data.frame(gsm = common,
                 risk = risk[common],
                 time = as.numeric(surv[common, "OS_years"]),
                 status = as.integer(surv[common, "OS_event"]))
df <- df[complete.cases(df) & df$time > 0, ]
cat(sprintf("Validation samples (risk + survival): %d (events=%d)\n", nrow(df), sum(df$status)))

# ---- 5. Discrimination ----
cidx <- concordance(coxph(Surv(time, status) ~ risk, data = df))$concordance
# bootstrap 95% CI
set.seed(42)
bc <- replicate(1000, {
  i <- sample(nrow(df), replace = TRUE)
  tryCatch(concordance(coxph(Surv(df$time[i], df$status[i]) ~ df$risk[i]))$concordance,
           error = function(e) NA_real_)
})
ci <- quantile(bc, c(.025, .975), na.rm = TRUE)
cat(sprintf("C-index = %.4f (95%% CI %.4f-%.4f)\n", cidx, ci[1], ci[2]))

# time-dependent AUC (1/3/5 yr)
tr <- tryCatch(timeROC(T = df$time, delta = df$status, marker = df$risk,
                       cause = 1, times = c(1, 3, 5), iid = FALSE), error = function(e) NULL)
auc <- if (!is.null(tr)) tr$AUC else c(NA, NA, NA)
cat(sprintf("time-AUC 1/3/5yr = %.3f / %.3f / %.3f\n", auc[1], auc[2], auc[3]))

# ---- 6. KM: locked TCGA cutoff (primary) + GSE-median (sensitivity) ----
df$grp_tcga <- ifelse(df$risk > cutoff, "High", "Low")
df$grp_med  <- ifelse(df$risk > median(df$risk), "High", "Low")
p_tcga <- tryCatch(survdiff(Surv(time, status) ~ grp_tcga, data = df), error = function(e) NULL)
p_med  <- survdiff(Surv(time, status) ~ grp_med, data = df)
lp_tcga <- if (!is.null(p_tcga)) 1 - pchisq(p_tcga$chisq, length(p_tcga$n) - 1) else NA
lp_med  <- 1 - pchisq(p_med$chisq, length(p_med$n) - 1)
cat(sprintf("KM logrank p: TCGA-cutoff=%.3g (High n=%d/Low n=%d) | GSE-median=%.3g\n",
            lp_tcga, sum(df$grp_tcga == "High"), sum(df$grp_tcga == "Low"), lp_med))

# ---- 7. Save (robust: equal-length character vectors) ----
mk <- c("cohort","platform","n","events","genes_found","C_index","C_CI_low","C_CI_high",
        "AUC_1yr","AUC_3yr","AUC_5yr","KM_p_TCGAcutoff","KM_p_GSEmedian")
vv <- c("GSE16011","Affymetrix microarray",
        as.character(nrow(df)), as.character(sum(df$status)),
        sprintf("%d/%d", length(found), length(need)),
        sprintf("%.4f", as.numeric(cidx)), sprintf("%.4f", as.numeric(ci[1])),
        sprintf("%.4f", as.numeric(ci[2])),
        sprintf("%.3f", as.numeric(auc[1])), sprintf("%.3f", as.numeric(auc[2])),
        sprintf("%.3f", as.numeric(auc[3])),
        sprintf("%.3g", as.numeric(lp_tcga)), sprintf("%.3g", as.numeric(lp_med)))
stopifnot(length(mk) == length(vv))
write.csv(data.frame(metric = mk, value = vv, stringsAsFactors = FALSE),
          file.path(PROJ, "results/gse16011_uirs_validation.csv"), row.names = FALSE)
saveRDS(df, file.path(PROJ, "data/processed/gse16011_uirs_validation.rds"))
cat(sprintf("FORKC_SUMMARY n=%d events=%d genes=%d/%d C=%.4f CI=%.4f-%.4f AUC=%.3f/%.3f/%.3f KMp_tcga=%.3g KMp_med=%.3g\n",
            nrow(df), sum(df$status), length(found), length(need), cidx, ci[1], ci[2],
            auc[1], auc[2], auc[3], lp_tcga, lp_med))
