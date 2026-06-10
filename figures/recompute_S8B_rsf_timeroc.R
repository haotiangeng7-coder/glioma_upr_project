#!/usr/bin/env Rscript
# Recompute S8B: time-dependent ROC for the BEST RSF model (ElasticNet07 + RSF)
# on CGGA-batch1 (blind test), reproducing the canonical blind C-index = 0.7173.
# Replays the EXACT best-combo blind-eval logic from 05_ml_model/02_ml_combinations.R
# on the SAVED preprocessed objects in ml_combination_results.RData.
# Real-data computation only. No hardcoded AUC/C-index.
suppressPackageStartupMessages({
  library(survival); library(glmnet); library(randomForestSRC)
  library(timeROC); library(ggplot2); library(scales)
})
PROJ <- getwd()
DATA_PROC <- file.path(PROJ, "data/processed")
OUT  <- file.path(PROJ, "figures", "supp_assembled")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

SEED <- 42
ML_NESTED_INNER <- 5

load(file.path(DATA_PROC, "ml_combination_results.RData"))
stopifnot(best_fs == "ElasticNet07", best_build == "RSF")
cat(sprintf("Best combo: %s + %s | canonical blind C = %.4f\n",
            best_fs, best_build, blind_results[["CGGA_batch1"]]))

# --- Exact fs/build functions copied verbatim from 02_ml_combinations.R ---
fs_enet07 <- function(x, y, seed = SEED) {
  set.seed(seed)
  cv <- cv.glmnet(x, y, family = "cox", alpha = 0.7, nfolds = ML_NESTED_INNER)
  coefs <- coef(cv, s = "lambda.min")
  selected <- rownames(coefs)[coefs[, 1] != 0]
  if (length(selected) == 0) selected <- colnames(x)
  list(selected = selected)
}
build_rsf <- function(x, y, seed = SEED) {
  set.seed(seed)
  df <- data.frame(time = y[, "time"], status = y[, "status"], x)
  gn <- colnames(x)
  fit <- rfsrc(Surv(time, status) ~ ., data = df, ntree = 1000, nodesize = 10, seed = seed)
  list(predict_fn = function(newx) {
    nd <- data.frame(time = 0, status = 0, newx[, gn, drop = FALSE])
    predict(fit, newdata = nd)$predicted
  })
}
calc_cindex <- function(pred_risk, time, status) {
  fit <- coxph(Surv(time, status) ~ pred_risk); concordance(fit)$concordance
}

# --- Replay blind-eval block verbatim ---
x_full <- as.matrix(train_expr[, genes])
y_full <- Surv(train_expr$time, train_expr$status)
fs_final  <- fs_enet07(x_full, y_full, seed = SEED)
sel_final <- intersect(fs_final$selected, colnames(x_full))
if (length(sel_final) < 2) sel_final <- colnames(x_full)
cat(sprintf("ElasticNet07 selected %d genes: %s\n",
            length(sel_final), paste(sel_final, collapse = ", ")))
built_final <- build_rsf(x_full[, sel_final, drop = FALSE], y_full, seed = SEED)

bt <- blind_test[["CGGA_batch1"]]
pred_bt <- built_final$predict_fn(
  data.frame(bt[, sel_final, drop = FALSE], check.names = FALSE))
c_bt <- calc_cindex(pred_bt, bt$time, bt$status)
cat(sprintf("Reproduced RSF C-index on CGGA-batch1: %.4f (canonical 0.7173)\n", c_bt))

# --- time-dependent ROC at 1/3/5 yr ---
cat(sprintf("CGGA time range (days): [%.1f, %.1f]\n", min(bt$time), max(bt$time)))
yr <- 365.25
times <- c(1*yr, 3*yr, 5*yr)
roc <- timeROC(T = bt$time, delta = bt$status, marker = pred_bt, cause = 1,
               times = times, iid = FALSE)
auc <- roc$AUC; names(auc) <- c("1yr","3yr","5yr")
cat("\n=== time-ROC AUC ===\n"); print(round(auc, 4))
cat(sprintf("5-year AUC = %.3f (caption 0.755)\n", auc["5yr"]))

# --- Plot ---
lvls <- c("1yr","3yr","5yr")
plotdf <- do.call(rbind, lapply(seq_along(times), function(i)
  data.frame(FPR = roc$FP[, i], TPR = roc$TP[, i],
             horizon = factor(lvls[i], levels = lvls))))
auc_lab <- sprintf("%s (AUC = %.3f)", c("1-year","3-year","5-year"), auc)
pal <- c("1yr" = "#4DBBD5", "3yr" = "#00A087", "5yr" = "#E64B35")
p <- ggplot(plotdf, aes(FPR, TPR, color = horizon)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  geom_line(linewidth = 0.9) +
  scale_color_manual(values = pal, labels = auc_lab, name = NULL) +
  coord_equal() +
  labs(x = "1 - Specificity", y = "Sensitivity",
       title = "RSF time-dependent ROC (CGGA-batch1)") +
  theme_bw(base_size = 12) +
  theme(legend.position = c(0.98, 0.02), legend.justification = c(1, 0),
        legend.background = element_rect(fill = alpha("white", 0.7), color = NA),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
        panel.grid.minor = element_blank())
ggsave(file.path(OUT, "_panel_S8B_rsf_timeroc.pdf"), p, width = 5.0, height = 5.0)
cat("\nSaved _panel_S8B_rsf_timeroc.pdf\n")
saveRDS(list(auc = auc, cindex = c_bt, n = nrow(bt), selected = sel_final),
        file.path(OUT, "_S8B_values.rds"))
