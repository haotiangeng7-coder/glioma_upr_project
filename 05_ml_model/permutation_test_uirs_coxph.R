#!/usr/bin/env Rscript
# ============================================================================
# External permutation test for the FINAL UIRS model (LASSO -> CoxPH + strata(IDH))
# on the CGGA-batch1 blind test set. Resolves stats-reviewer F4 for the model
# that is actually used as UIRS (the benchmark's own permutation tested the
# top RSF model; this one tests CoxPH).
#
# Null: shuffle TCGA (time,status), refit the full UIRS pipeline (LASSO selection
# -> CoxPH with strata(IDH)), predict on the FIXED CGGA-batch1, compute C-index.
# p = (#perm C >= observed C + 1) / (n_valid + 1)
# ============================================================================
set.seed(42)
suppressPackageStartupMessages({library(survival); library(glmnet)})

PROJ <- getwd()
load(file.path(PROJ, "data/processed/coxph_interpretable_model.RData"))  # coxph_results
N_PERM <- 1000

td   <- coxph_results$train_df            # time,status,<candidate genes>,IDH_status
v1   <- coxph_results$cgga1_val_df        # prepared CGGA-batch1 blind test
cand <- coxph_results$candidate_genes     # 17 candidate genes (input to LASSO)
cat(sprintf("Train n=%d (events=%d) | CGGA-batch1 n=%d | candidate genes=%d\n",
            nrow(td), sum(td$status), nrow(v1), length(cand)))

# C-index for a risk/lp vector (higher = worse): fit Surv~lp, take concordance
cidx <- function(lp, time, status) {
  tryCatch(concordance(coxph(Surv(time, status) ~ lp))$concordance, error = function(e) NA_real_)
}

# One full UIRS fit given a training data.frame -> returns risk score on v1.
# Risk = sum(beta_i * expr_i) using the fitted CoxPH gene coefficients (the locked
# signature form). strata(IDH) only sets the baseline hazard; it does not enter the
# linear risk score, so v1 needs the genes only (it lacks an IDH_status column).
fit_predict <- function(train_df, seed) {
  x <- as.matrix(train_df[, cand]); y <- Surv(train_df$time, train_df$status)
  set.seed(seed)
  cv <- cv.glmnet(x, y, family = "cox", alpha = 1, nfolds = 10, type.measure = "C")
  cf <- coef(cv, s = "lambda.min"); sel <- rownames(cf)[as.numeric(cf) != 0]
  if (length(sel) < 2) sel <- cand
  fml <- as.formula(paste("Surv(time, status) ~", paste(sel, collapse = " + "), "+ strata(IDH_status)"))
  m <- coxph(fml, data = train_df)
  b <- coef(m); b <- b[!is.na(b)]
  g <- intersect(names(b), colnames(v1))
  as.numeric(as.matrix(v1[, g, drop = FALSE]) %*% b[g])
}

obs <- cidx(fit_predict(td, 42), v1$time, v1$status)
cat(sprintf("Observed UIRS C-index on CGGA-batch1: %.4f\n", obs))

perm <- numeric(N_PERM)
t0 <- proc.time()
for (i in seq_len(N_PERM)) {
  if (i %% 100 == 0) cat(sprintf("  perm %d/%d (%.0fs)\n", i, N_PERM, (proc.time()-t0)[3]))
  tp <- td; idx <- sample(nrow(td)); tp$time <- td$time[idx]; tp$status <- td$status[idx]
  perm[i] <- tryCatch(cidx(fit_predict(tp, 42 + i), v1$time, v1$status), error = function(e) NA_real_)
}
valid <- perm[!is.na(perm)]
pval  <- (sum(valid >= obs) + 1) / (length(valid) + 1)
cat(sprintf("\n=== UIRS (CoxPH) external permutation ===\n"))
cat(sprintf("observed=%.4f  null_mean=%.4f  null_sd=%.4f  n_valid=%d  p=%.4f\n",
            obs, mean(valid), sd(valid), length(valid), pval))

write.csv(data.frame(metric = c("model","test_set","observed_cindex","null_mean","null_sd","n_perm","p_value"),
                     value  = c("LASSO->CoxPH+strata(IDH) (UIRS)","CGGA_batch1",
                                sprintf("%.4f",obs), sprintf("%.4f",mean(valid)),
                                sprintf("%.4f",sd(valid)), length(valid), sprintf("%.4f",pval))),
          file.path(PROJ,"results/permutation_test_uirs_coxph.csv"), row.names = FALSE)
saveRDS(list(observed=obs, perm=perm, pval=pval),
        file.path(PROJ,"data/processed/permutation_test_uirs_coxph.rds"))
cat("written results/permutation_test_uirs_coxph.csv\n")
