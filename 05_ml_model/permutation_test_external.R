###############################################################################
# permutation_test_external.R
# Fix: Permutation test evaluated on CGGA-batch1 (blind test set)
# instead of training set resubstitution
#
# Problem: RSF overfits -> null distribution on training data is ~0.87
#          -> observed 0.87 is not significant (p=0.758)
# Solution: Evaluate on external CGGA-batch1 where real signal generalizes
#           but permuted (noise) models do not
###############################################################################

set.seed(42)

library(survival)
library(glmnet)
library(randomForestSRC)
library(ggplot2)

# === Paths ===
PROJECT_DIR <- getwd()
DATA_PROC   <- file.path(PROJECT_DIR, "data", "processed")
FIG_DIR     <- file.path(PROJECT_DIR, "figures")
RES_DIR     <- file.path(PROJECT_DIR, "results")

N_PERM <- 500  # 500 permutations (halved for speed)

# =============================================================================
# 1. Load data
# =============================================================================
message("=== Loading data ===")

load(file.path(DATA_PROC, "feature_selection_results.RData"))
load(file.path(DATA_PROC, "tcga_glioma_expression.RData"))
load(file.path(DATA_PROC, "cgga_data.RData"))

genes <- candidate_genes
message(sprintf("Candidate genes: %d (%s)", length(genes), paste(genes, collapse=", ")))

# --- Training set (TCGA) ---
train_expr <- surv_df[, c("time", "status", genes)]
train_expr <- train_expr[complete.cases(train_expr), ]
message(sprintf("Training set (TCGA): %d samples", nrow(train_expr)))

# --- Test set (CGGA-batch1) ---
# Prepare CGGA-batch1 data
cgga_expr <- cgga_b1$expr
cgga_clin <- cgga_b1$clinical

common_genes <- intersect(genes, rownames(cgga_expr))
message(sprintf("Genes found in CGGA-batch1: %d / %d", length(common_genes), length(genes)))

# Build test matrix with all genes (fill missing with 0)
test_mat <- matrix(0, nrow = ncol(cgga_expr), ncol = length(genes),
                   dimnames = list(colnames(cgga_expr), genes))
test_mat[, common_genes] <- t(cgga_expr[common_genes, ])

# Match samples with clinical data
common_ids <- intersect(rownames(test_mat), cgga_clin$Sample_ID)

test_df <- data.frame(
  time   = cgga_clin$OS.time[match(common_ids, cgga_clin$Sample_ID)],
  status = cgga_clin$OS_status[match(common_ids, cgga_clin$Sample_ID)],
  test_mat[common_ids, , drop = FALSE],
  check.names = FALSE
)
test_df <- test_df[complete.cases(test_df) & test_df$time > 0, ]
message(sprintf("Test set (CGGA-batch1): %d samples", nrow(test_df)))

# =============================================================================
# 2. C-index helper
# =============================================================================
calc_cindex <- function(pred_risk, time, status) {
  tryCatch({
    fit <- coxph(Surv(time, status) ~ pred_risk)
    concordance(fit)$concordance
  }, error = function(e) NA_real_)
}

# =============================================================================
# 3. Train observed model (LASSO + RSF) on full TCGA, evaluate on CGGA-batch1
# =============================================================================
message("=== Training observed model (LASSO + RSF) ===")

train_observed_model <- function(x_train, y_train, seed = 42) {
  # Step 1: LASSO feature selection
  set.seed(seed)
  cv_lasso <- cv.glmnet(x_train, y_train, family = "cox", alpha = 1, nfolds = 5)
  coefs <- coef(cv_lasso, s = "lambda.min")
  selected <- rownames(coefs)[coefs[, 1] != 0]
  if (length(selected) < 2) selected <- colnames(x_train)

  # Step 2: RSF model building on selected features
  set.seed(seed)
  rsf_df <- data.frame(time = y_train[, "time"], status = y_train[, "status"],
                       x_train[, selected, drop = FALSE])
  rsf_fit <- rfsrc(Surv(time, status) ~ ., data = rsf_df,
                   ntree = 1000, nodesize = 10, seed = seed)

  list(
    selected = selected,
    lasso_model = cv_lasso,
    rsf_model = rsf_fit,
    predict_fn = function(newx) {
      nd <- data.frame(time = 0, status = 0,
                       newx[, selected, drop = FALSE])
      predict(rsf_fit, newdata = nd)$predicted
    }
  )
}

x_train <- as.matrix(train_expr[, genes])
y_train <- Surv(train_expr$time, train_expr$status)

obs_model <- train_observed_model(x_train, y_train, seed = 42)
message(sprintf("LASSO selected %d genes: %s",
                length(obs_model$selected), paste(obs_model$selected, collapse = ", ")))

# Predict on CGGA-batch1
obs_pred <- obs_model$predict_fn(
  data.frame(test_df[, genes, drop = FALSE], check.names = FALSE)
)
observed_cindex <- calc_cindex(obs_pred, test_df$time, test_df$status)
message(sprintf("Observed C-index on CGGA-batch1: %.4f", observed_cindex))

# =============================================================================
# 4. Permutation test (500 times)
# =============================================================================
message(sprintf("=== Permutation test (%d permutations) ===", N_PERM))
message("Shuffling TCGA training labels, re-training LASSO+RSF, evaluating on CGGA-batch1")

set.seed(42)
perm_cindices <- numeric(N_PERM)

t_start <- proc.time()

for (perm_i in seq_len(N_PERM)) {
  if (perm_i %% 50 == 0 || perm_i == 1) {
    elapsed <- (proc.time() - t_start)[3]
    eta <- elapsed / perm_i * (N_PERM - perm_i)
    message(sprintf("  Permutation %d/%d (elapsed: %.0fs, ETA: %.0fs)",
                    perm_i, N_PERM, elapsed, eta))
  }

  # Shuffle survival labels in training set
  perm_idx <- sample(nrow(train_expr))
  time_perm   <- train_expr$time[perm_idx]
  status_perm <- train_expr$status[perm_idx]

  tryCatch({
    y_perm <- Surv(time_perm, status_perm)

    # LASSO feature selection on permuted data
    set.seed(42 + perm_i)
    cv_perm <- cv.glmnet(x_train, y_perm, family = "cox", alpha = 1, nfolds = 5)
    coefs_perm <- coef(cv_perm, s = "lambda.min")
    sel_perm <- rownames(coefs_perm)[coefs_perm[, 1] != 0]
    if (length(sel_perm) < 2) sel_perm <- colnames(x_train)

    # RSF on permuted data
    set.seed(42 + perm_i)
    rsf_df_perm <- data.frame(time = time_perm, status = status_perm,
                              x_train[, sel_perm, drop = FALSE])
    rsf_perm <- rfsrc(Surv(time, status) ~ ., data = rsf_df_perm,
                      ntree = 1000, nodesize = 10, seed = 42 + perm_i)

    # Predict on CGGA-batch1 (NOT shuffled)
    nd_perm <- data.frame(time = 0, status = 0,
                          test_df[, sel_perm, drop = FALSE])
    pred_perm <- predict(rsf_perm, newdata = nd_perm)$predicted

    perm_cindices[perm_i] <- calc_cindex(pred_perm, test_df$time, test_df$status)
  }, error = function(e) {
    perm_cindices[perm_i] <<- NA_real_
  })
}

elapsed_total <- (proc.time() - t_start)[3]
message(sprintf("Permutation test completed in %.1f minutes", elapsed_total / 60))

# =============================================================================
# 5. Compute empirical p-value
# =============================================================================
valid_perm <- perm_cindices[!is.na(perm_cindices)]
perm_pvalue <- mean(valid_perm >= observed_cindex)

message(sprintf("\n========== RESULTS =========="))
message(sprintf("Observed C-index (CGGA-batch1): %.4f", observed_cindex))
message(sprintf("Null distribution mean:         %.4f", mean(valid_perm)))
message(sprintf("Null distribution SD:           %.4f", sd(valid_perm)))
message(sprintf("Null distribution range:        [%.4f, %.4f]", min(valid_perm), max(valid_perm)))
message(sprintf("Valid permutations:             %d / %d", length(valid_perm), N_PERM))
message(sprintf("Empirical p-value:              %.4f", perm_pvalue))
if (perm_pvalue == 0) {
  message(sprintf("  (p < %.4f, i.e., < 1/%d)", 1/length(valid_perm), length(valid_perm)))
}
message(sprintf("=================================\n"))

# =============================================================================
# 6. Visualization
# =============================================================================
message("=== Generating permutation test plot ===")

THEME_PUBLICATION <- theme_bw() +
  theme(
    plot.title   = element_text(size = 14, face = "bold", hjust = 0.5),
    axis.title   = element_text(size = 12),
    axis.text    = element_text(size = 10),
    legend.title = element_text(size = 11),
    legend.text  = element_text(size = 10),
    panel.grid.minor = element_blank()
  )

perm_df <- data.frame(cindex = valid_perm)

# Format p-value text
if (perm_pvalue == 0) {
  p_text <- sprintf("p < %.4f", 1 / length(valid_perm))
} else {
  p_text <- sprintf("p = %.4f", perm_pvalue)
}

p_perm <- ggplot(perm_df, aes(x = cindex)) +
  geom_histogram(bins = 40, fill = "grey70", color = "white", alpha = 0.9) +
  geom_vline(xintercept = observed_cindex, color = "#E64B35",
             linetype = "dashed", linewidth = 1.2) +
  annotate("text", x = observed_cindex, y = Inf, vjust = 2,
           label = sprintf("Observed C-index = %.4f\n%s", observed_cindex, p_text),
           color = "#E64B35", hjust = -0.1, size = 4.5, fontface = "bold") +
  labs(x = "C-index on CGGA-batch1 (External Validation)",
       y = "Count",
       title = sprintf("Permutation Test (n = %d)", length(valid_perm)),
       subtitle = "Null: Shuffle TCGA training labels, re-train LASSO+RSF, predict on CGGA-batch1") +
  THEME_PUBLICATION

ggsave(file.path(FIG_DIR, "SuppFig_permutation_test.pdf"), p_perm,
       width = 8, height = 5)
message(sprintf("Plot saved: %s", file.path(FIG_DIR, "SuppFig_permutation_test.pdf")))

# Also copy to Supplementary_Figures if directory exists
supp_dir <- file.path(FIG_DIR, "Supplementary_Figures")
if (dir.exists(supp_dir)) {
  file.copy(file.path(FIG_DIR, "SuppFig_permutation_test.pdf"),
            file.path(supp_dir, "SuppFig_permutation_test.pdf"),
            overwrite = TRUE)
  message(sprintf("Also saved to: %s", file.path(supp_dir, "SuppFig_permutation_test.pdf")))
}

# =============================================================================
# 7. Update saved results
# =============================================================================
message("=== Updating saved results ===")

# Load existing results to update permutation fields
load(file.path(DATA_PROC, "ml_combination_results.RData"))

# Overwrite permutation-related variables
observed_cindex_external <- observed_cindex
perm_cindices_external   <- perm_cindices
perm_pvalue_external     <- perm_pvalue

# Save back with updated permutation results
save(
  ml_results,
  all_combo_results,
  best_combo, best_fs, best_build,
  perm_cindices, perm_pvalue, observed_cindex,
  perm_cindices_external, perm_pvalue_external, observed_cindex_external,
  blind_results,
  train_expr, validation_data, tuning_data, blind_test, exploratory,
  genes,
  FS_ALGOS, BUILD_ALGOS,
  file = file.path(DATA_PROC, "ml_combination_results.RData")
)

# Update summary CSV
summary_df <- data.frame(
  metric = c("best_combo", "observed_cindex_external", "perm_pvalue_external",
             "perm_null_mean", "perm_null_sd", "n_permutations",
             names(blind_results)),
  value  = c(best_combo,
             sprintf("%.4f", observed_cindex_external),
             sprintf("%.4f", perm_pvalue_external),
             sprintf("%.4f", mean(valid_perm)),
             sprintf("%.4f", sd(valid_perm)),
             as.character(length(valid_perm)),
             sapply(blind_results, function(x) sprintf("%.4f", x)))
)
write.csv(summary_df, file.path(RES_DIR, "best_model_summary.csv"), row.names = FALSE)

message("\n=== Permutation test (external validation) completed ===")
