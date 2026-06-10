#!/usr/bin/env Rscript
# ============================================================================
# Fork B: Can adding routine molecular/clinical features improve UIRS
# discrimination within IDH-wildtype glioma (where gene-only C-index ~0.586)?
#
# Compare on TCGA IDH-WT (train) -> CGGA-batch1 IDH-WT (external test):
#   (1) Gene-only : UIRS 15-gene linear score (locked betas)
#   (2) Clinical  : Age + Grade + MGMT
#   (3) Combined  : UIRS score + Age + Grade + MGMT
# Metric: Harrell C on CGGA-b1 IDH-WT, bootstrap 95% CI; delta(combined-gene)
# with bootstrap CI. Honest reporting: if combined does not beat gene-only by a
# meaningful margin (and the delta CI includes 0), report no improvement.
# ============================================================================
set.seed(42)
suppressPackageStartupMessages(library(survival))
PROJ <- getwd()

load(file.path(PROJ, "data/processed/coxph_interpretable_model.RData")) # coxph_results
load(file.path(PROJ, "data/processed/tcga_glioma_expression.RData"))    # clinical_valid
load(file.path(PROJ, "data/processed/cgga_data.RData"))                 # cgga_b1
betas <- coxph_results$betas; genes <- names(betas)

norm_mgmt  <- function(x) ifelse(grepl("^methyl", tolower(x)), 1L,
                          ifelse(grepl("unmethyl|un-methyl", tolower(x)), 0L, NA_integer_))
norm_grade <- function(x) {
  x <- toupper(as.character(x))
  ifelse(grepl("IV|G4|4", x), 4L, ifelse(grepl("III|G3|3", x), 3L,
  ifelse(grepl("II|G2|2", x), 2L, NA_integer_)))
}
cidx <- function(lp, t, s) tryCatch(concordance(coxph(Surv(t, s) ~ lp))$concordance,
                                    error = function(e) NA_real_)

# ---------- TCGA IDH-WT training frame ----------
td <- coxph_results$train_df
td <- td[!is.na(td$IDH_status) & td$IDH_status == "WT", ]
cl <- clinical_valid[rownames(td), ]
td$Age   <- as.numeric(cl$age_at_index)
td$Grade <- norm_grade(cl$Grade)
td$MGMT  <- norm_mgmt(cl$MGMT_status)
td$uirs  <- as.numeric(as.matrix(td[, genes]) %*% betas)
td_full  <- td[complete.cases(td[, c("time","status","Age","Grade","MGMT", genes)]) & td$time > 0, ]
cat(sprintf("TCGA IDH-WT train: n=%d events=%d (complete clinical)\n",
            nrow(td_full), sum(td_full$status)))

# ---------- CGGA-batch1 IDH-WT test frame ----------
cb <- cgga_b1$clinical
cb$IDHwt <- !is.na(cb$IDH_status) & cb$IDH_status == "Wildtype"
expr <- cgga_b1$expr
gpres <- intersect(genes, rownames(expr))
em <- log2(t(expr[gpres, , drop = FALSE]) + 1)            # samples x genes, log2
rownames(em) <- colnames(expr)
test <- data.frame(
  sid    = cb$Sample_ID,
  time   = as.numeric(cb$OS.time),
  status = as.integer(cb$OS_status),
  Age    = as.numeric(cb$Age),
  Grade  = norm_grade(cb$Grade),
  MGMT   = norm_mgmt(cb$MGMT_status),
  IDHwt  = cb$IDHwt,
  stringsAsFactors = FALSE
)
test <- test[test$IDHwt %in% TRUE, ]
test <- test[test$sid %in% rownames(em), ]
test$uirs <- as.numeric(em[test$sid, gpres, drop = FALSE] %*% betas[gpres])
test <- test[complete.cases(test[, c("time","status","Age","Grade","MGMT","uirs")]) & test$time > 0, ]
cat(sprintf("CGGA-b1 IDH-WT test: n=%d events=%d\n", nrow(test), sum(test$status)))

# ---------- Fit on TCGA-WT, predict on CGGA-b1-WT ----------
m_gene <- coxph(Surv(time, status) ~ uirs, data = td_full)
m_clin <- coxph(Surv(time, status) ~ Age + Grade + MGMT, data = td_full)
m_comb <- coxph(Surv(time, status) ~ uirs + Age + Grade + MGMT, data = td_full)

lp_gene <- predict(m_gene, newdata = test, type = "lp")
lp_clin <- predict(m_clin, newdata = test, type = "lp")
lp_comb <- predict(m_comb, newdata = test, type = "lp")

c_gene <- cidx(lp_gene, test$time, test$status)
c_clin <- cidx(lp_clin, test$time, test$status)
c_comb <- cidx(lp_comb, test$time, test$status)

# ---------- Bootstrap CIs + delta(combined - gene) ----------
set.seed(42); B <- 1000
boot <- replicate(B, {
  i <- sample(nrow(test), replace = TRUE)
  tt <- test[i, ]
  cg <- cidx(lp_gene[i], tt$time, tt$status)
  cc <- cidx(lp_comb[i], tt$time, tt$status)
  cl_ <- cidx(lp_clin[i], tt$time, tt$status)
  c(cg, cc, cl_, cc - cg)
})
ci <- function(v) quantile(v, c(.025, .975), na.rm = TRUE)
ci_g <- ci(boot[1, ]); ci_c <- ci(boot[2, ]); ci_cl <- ci(boot[3, ]); ci_d <- ci(boot[4, ])

cat(sprintf("\n=== IDH-WT external C-index (CGGA-b1) ===\n"))
cat(sprintf("Gene-only (UIRS): %.3f (95%% CI %.3f-%.3f)\n", c_gene, ci_g[1], ci_g[2]))
cat(sprintf("Clinical-only   : %.3f (95%% CI %.3f-%.3f)\n", c_clin, ci_cl[1], ci_cl[2]))
cat(sprintf("Combined        : %.3f (95%% CI %.3f-%.3f)\n", c_comb, ci_c[1], ci_c[2]))
cat(sprintf("Delta(comb-gene): %.3f (95%% CI %.3f-%.3f)\n", c_comb - c_gene, ci_d[1], ci_d[2]))
improved <- (c_comb - c_gene) > 0.03 && ci_d[1] > 0
cat(sprintf("MEANINGFUL_IMPROVEMENT=%s\n", improved))

write.csv(data.frame(
  model = c("Gene-only (UIRS)","Clinical-only (Age+Grade+MGMT)","Combined",
            "Delta(combined-gene)"),
  cindex = sprintf("%.3f", c(c_gene, c_clin, c_comb, c_comb - c_gene)),
  ci_low = sprintf("%.3f", c(ci_g[1], ci_cl[1], ci_c[1], ci_d[1])),
  ci_high= sprintf("%.3f", c(ci_g[2], ci_cl[2], ci_c[2], ci_d[2])),
  stringsAsFactors = FALSE),
  file.path(PROJ, "results/idhwt_model_improvement.csv"), row.names = FALSE)
cat("written results/idhwt_model_improvement.csv\n")
cat(sprintf("FORKB_SUMMARY train_n=%d test_n=%d gene=%.3f clin=%.3f comb=%.3f delta=%.3f dCI=%.3f..%.3f improved=%s\n",
            nrow(td_full), nrow(test), c_gene, c_clin, c_comb, c_comb-c_gene, ci_d[1], ci_d[2], improved))
