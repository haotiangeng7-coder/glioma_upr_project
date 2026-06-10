#!/usr/bin/env Rscript
# =============================================================================
# CoxPH Interpretable Risk Score Model for Glioma UPR Project
#
# Purpose: Build a transparent CoxPH model with explicit beta coefficients
#          to complement the black-box LASSO+RSF model.
#
# Pipeline: LASSO Cox variable selection -> CoxPH with strata(IDH) ->
#           Linear risk score: Risk = sum(beta_i * expr_i)
#
# Key note: TCGA surv_df is in log2 scale. CGGA is in raw FPKM scale,
#           so CGGA must be log2(x+1) transformed before applying the model.
#
# C-index note: survival::concordance(Surv ~ x) returns P(larger x -> longer T).
#   For risk scores (higher = worse), the C-statistic = 1 - concordance(Surv ~ risk).
#
# Author: Generated for glioma UPR project
# Date: 2026-03-19
# =============================================================================

set.seed(42)

suppressPackageStartupMessages({
  library(survival)
  library(glmnet)
  library(survminer)
  library(timeROC)
  library(ggplot2)
})

cat("========================================\n")
cat("CoxPH Interpretable Model Construction\n")
cat("========================================\n\n")

# =============================================================================
# 1. Load Data
# =============================================================================
cat("[1/8] Loading data...\n")

load("data/processed/feature_selection_results.RData")
load("data/processed/tcga_glioma_expression.RData")
load("data/processed/cgga_data.RData")
load("data/processed/risk_model_final.RData")

rsf_auc <- auc_summary
cat("  candidate_genes (n=", length(candidate_genes), "):",
    paste(candidate_genes, collapse = ", "), "\n")
cat("  TCGA training samples:", nrow(surv_df), "\n")
cat("  CGGA batch1 clinical:", nrow(cgga_b1$clinical), "\n")
cat("  CGGA batch2 clinical:", nrow(cgga_b2$clinical), "\n\n")

# =============================================================================
# 2. Prepare Training Data (TCGA) with IDH status
# =============================================================================
cat("[2/8] Preparing training data with IDH stratification...\n")

common_samples <- intersect(rownames(surv_df), rownames(clinical_valid))
train_df <- surv_df[common_samples, ]
train_df$IDH_status <- clinical_valid[common_samples, "IDH_status"]

keep <- complete.cases(train_df[, c("time", "status", "IDH_status")]) & train_df$time > 0
train_df <- train_df[keep, ]
cat("  Training samples after filtering:", nrow(train_df), "\n")
cat("  IDH distribution:\n")
print(table(train_df$IDH_status))

cat("\n  Expression scale check (log2, first 3 genes):\n")
for (g in candidate_genes[1:3]) {
  cat("    ", g, ": range =", round(range(train_df[[g]], na.rm = TRUE), 2), "\n")
}
cat("\n")

# =============================================================================
# 3. LASSO Cox Variable Selection
# =============================================================================
cat("[3/8] LASSO Cox variable selection...\n")

x_train <- as.matrix(train_df[, candidate_genes])
y_train <- Surv(train_df$time, train_df$status)

set.seed(42)
cv_lasso <- cv.glmnet(x_train, y_train,
                       family = "cox",
                       alpha = 1,
                       nfolds = 10,
                       type.measure = "C")

cat("  lambda.min:", cv_lasso$lambda.min, "\n")
cat("  lambda.1se:", cv_lasso$lambda.1se, "\n")

lasso_coef <- coef(cv_lasso, s = "lambda.min")
lasso_coef_vec <- as.numeric(lasso_coef)
names(lasso_coef_vec) <- rownames(lasso_coef)

selected_genes <- names(lasso_coef_vec[lasso_coef_vec != 0])
cat("  Selected genes (non-zero at lambda.min, n=", length(selected_genes), "):\n")
cat("   ", paste(selected_genes, collapse = ", "), "\n\n")

for (g in selected_genes) {
  cat("    ", g, ":", round(lasso_coef_vec[g], 4), "\n")
}
cat("\n")

# =============================================================================
# 4. Fit CoxPH Model with strata(IDH)
# =============================================================================
cat("[4/8] Fitting CoxPH model with strata(IDH)...\n")

gene_terms <- paste(selected_genes, collapse = " + ")
cox_formula <- as.formula(paste("Surv(time, status) ~", gene_terms, "+ strata(IDH_status)"))
cat("  Formula:", deparse(cox_formula), "\n")

cox_model <- coxph(cox_formula, data = train_df)
cat("\n  CoxPH Model Summary:\n")
print(summary(cox_model))

# =============================================================================
# 5. PH Assumption Test
# =============================================================================
cat("\n[5/8] Testing Proportional Hazards assumption...\n")
ph_test <- cox.zph(cox_model)
print(ph_test)

ph_results <- data.frame(
  variable = rownames(ph_test$table),
  chisq = ph_test$table[, "chisq"],
  df = ph_test$table[, "df"],
  p_value = ph_test$table[, "p"],
  PH_satisfied = ifelse(ph_test$table[, "p"] > 0.05, "Yes", "No")
)
cat("\n  PH assumption summary:\n")
print(ph_results)
cat("\n")

# =============================================================================
# 6. Extract Model Coefficients Table
# =============================================================================
cat("[6/8] Extracting model coefficients...\n")

cox_summary <- summary(cox_model)
coef_table <- data.frame(
  Gene = rownames(cox_summary$coefficients),
  Beta = cox_summary$coefficients[, "coef"],
  HR = cox_summary$coefficients[, "exp(coef)"],
  HR_lower = cox_summary$conf.int[, "lower .95"],
  HR_upper = cox_summary$conf.int[, "upper .95"],
  SE = cox_summary$coefficients[, "se(coef)"],
  z = cox_summary$coefficients[, "z"],
  p_value = cox_summary$coefficients[, "Pr(>|z|)"],
  stringsAsFactors = FALSE
)
rownames(coef_table) <- NULL
coef_table <- coef_table[order(-abs(coef_table$Beta)), ]

cat("\n  Model Coefficients:\n")
print(coef_table, digits = 4)

write.csv(coef_table, "results/coxph_model_coefficients.csv", row.names = FALSE)
cat("\n  Saved: results/coxph_model_coefficients.csv\n\n")

# =============================================================================
# 7. Compute Risk Scores and Validate
# =============================================================================
cat("[7/8] Computing risk scores and validation...\n")

betas <- coef(cox_model)
model_genes_final <- names(betas)

# --- Helper: compute C-index (Harrell's C-statistic) for risk score ---
# For a risk score: higher risk = shorter survival = higher hazard
# survival::concordance(Surv ~ x) returns P(larger x -> longer T)
# C-statistic for risk = 1 - concordance(Surv ~ risk_score)
compute_cindex <- function(time, status, risk_score) {
  df_tmp <- data.frame(time = time, status = status, rs = risk_score)
  raw_conc <- concordance(Surv(time, status) ~ rs, data = df_tmp)$concordance
  cstat <- 1 - raw_conc
  return(cstat)
}

# --- Helper: prepare CGGA validation dataset with log2 transform ---
prepare_cgga_validation <- function(cgga_obj, model_genes, dataset_name) {
  cat("\n  Preparing", dataset_name, "...\n")

  expr_mat <- cgga_obj$expr
  clin <- cgga_obj$clinical

  if ("CGGA_ID" %in% colnames(clin)) {
    rownames(clin) <- clin$CGGA_ID
  }
  common <- intersect(colnames(expr_mat), rownames(clin))
  cat("    Matched samples:", length(common), "\n")

  avail <- intersect(model_genes, rownames(expr_mat))
  cat("    Available genes:", length(avail), "/", length(model_genes), "\n")

  val_df <- data.frame(
    time = clin[common, "OS.time"],
    status = clin[common, "OS_status"],
    IDH_status = clin[common, "IDH_status"]
  )

  # Log2(x+1) transform: CGGA is in raw FPKM scale, model was trained on log2
  for (g in model_genes) {
    if (g %in% rownames(expr_mat)) {
      val_df[[g]] <- log2(as.numeric(expr_mat[g, common]) + 1)
    } else {
      val_df[[g]] <- 0
      cat("    WARNING: Gene", g, "not found, imputed as 0\n")
    }
  }

  cat("    Expression scale after log2 transform:\n")
  for (g in model_genes[1:min(2, length(model_genes))]) {
    cat("      ", g, ": range =", round(range(val_df[[g]], na.rm = TRUE), 2), "\n")
  }

  ok <- complete.cases(val_df[, c("time", "status", "IDH_status")]) & val_df$time > 0
  val_df <- val_df[ok, ]
  cat("    Valid samples after filtering:", nrow(val_df), "\n")

  return(val_df)
}

# --- A) TCGA Training Set ---
cat("\n  A) TCGA Training Set:\n")
train_risk <- as.numeric(predict(cox_model, newdata = train_df, type = "lp"))
names(train_risk) <- rownames(train_df)
train_median <- median(train_risk)
train_df$risk_score <- train_risk
train_df$risk_group <- ifelse(train_risk >= train_median, "High", "Low")
cat("    Risk score range:", round(range(train_risk), 4), "\n")
cat("    Median cutoff:", round(train_median, 4), "\n")
cat("    High/Low:", table(train_df$risk_group), "\n")

# C-index from model (handles strata correctly)
cindex_train <- summary(cox_model)$concordance[1]
cat("    C-index (strata-adjusted):", round(cindex_train, 4), "\n")

# --- B) CGGA Batch1 ---
cgga1_val_df <- prepare_cgga_validation(cgga_b1, model_genes_final, "CGGA Batch1 (Blind Test)")

cgga1_risk <- as.numeric(as.matrix(cgga1_val_df[, model_genes_final]) %*% betas)
cgga1_val_df$risk_score <- cgga1_risk
# Use cohort-specific median for KM grouping
cgga1_median <- median(cgga1_risk)
cgga1_val_df$risk_group <- ifelse(cgga1_risk >= cgga1_median, "High", "Low")
cat("    Risk score range:", round(range(cgga1_risk), 4), "\n")
cat("    Cohort median cutoff:", round(cgga1_median, 4), "\n")

cindex_b1 <- compute_cindex(cgga1_val_df$time, cgga1_val_df$status, cgga1_val_df$risk_score)
cat("    C-index:", round(cindex_b1, 4), "\n")

# --- C) CGGA Batch2 ---
cgga2_val_df <- prepare_cgga_validation(cgga_b2, model_genes_final, "CGGA Batch2 (Tuning)")

cgga2_risk <- as.numeric(as.matrix(cgga2_val_df[, model_genes_final]) %*% betas)
cgga2_val_df$risk_score <- cgga2_risk
cgga2_median <- median(cgga2_risk)
cgga2_val_df$risk_group <- ifelse(cgga2_risk >= cgga2_median, "High", "Low")
cat("    Risk score range:", round(range(cgga2_risk), 4), "\n")
cat("    Cohort median cutoff:", round(cgga2_median, 4), "\n")

cindex_b2 <- compute_cindex(cgga2_val_df$time, cgga2_val_df$status, cgga2_val_df$risk_score)
cat("    C-index:", round(cindex_b2, 4), "\n")

# =============================================================================
# 8. timeROC AUC (1, 3, 5 year)
# =============================================================================
cat("\n[8/8] Computing timeROC AUC...\n")

compute_timeroc <- function(df, label) {
  times <- c(365, 1095, 1825)

  roc_obj <- timeROC(
    T = df$time,
    delta = df$status,
    marker = df$risk_score,
    cause = 1,
    times = times,
    iid = TRUE
  )

  aucs <- roc_obj$AUC[match(times, roc_obj$times)]
  names(aucs) <- c("AUC_1yr", "AUC_3yr", "AUC_5yr")
  cat("    ", label, "- 1yr:", round(aucs[1], 3),
      " 3yr:", round(aucs[2], 3),
      " 5yr:", round(aucs[3], 3), "\n")

  return(list(roc = roc_obj, aucs = aucs))
}

roc_train <- compute_timeroc(train_df, "TCGA Training")
roc_b1 <- compute_timeroc(cgga1_val_df, "CGGA Batch1")
roc_b2 <- compute_timeroc(cgga2_val_df, "CGGA Batch2")

# =============================================================================
# Comparison Table: CoxPH vs RSF
# =============================================================================
cat("\n========================================\n")
cat("Model Comparison: CoxPH vs LASSO+RSF\n")
cat("========================================\n")

comparison_df <- data.frame(
  Dataset = c("TCGA (Training)", "CGGA_batch2 (Tuning)", "CGGA_batch1 (Blind Test)"),
  RSF_Cindex = c(rsf_auc$C_index[rsf_auc$Dataset == "TCGA (Training)"],
                  rsf_auc$C_index[rsf_auc$Dataset == "CGGA_batch2 (Tuning)"],
                  rsf_auc$C_index[rsf_auc$Dataset == "CGGA_batch1 (Blind Test)"]),
  CoxPH_Cindex = c(cindex_train, cindex_b2, cindex_b1),
  RSF_AUC_1yr = c(rsf_auc$AUC_1yr[rsf_auc$Dataset == "TCGA (Training)"],
                   rsf_auc$AUC_1yr[rsf_auc$Dataset == "CGGA_batch2 (Tuning)"],
                   rsf_auc$AUC_1yr[rsf_auc$Dataset == "CGGA_batch1 (Blind Test)"]),
  CoxPH_AUC_1yr = c(roc_train$aucs[1], roc_b2$aucs[1], roc_b1$aucs[1]),
  RSF_AUC_3yr = c(rsf_auc$AUC_3yr[rsf_auc$Dataset == "TCGA (Training)"],
                   rsf_auc$AUC_3yr[rsf_auc$Dataset == "CGGA_batch2 (Tuning)"],
                   rsf_auc$AUC_3yr[rsf_auc$Dataset == "CGGA_batch1 (Blind Test)"]),
  CoxPH_AUC_3yr = c(roc_train$aucs[2], roc_b2$aucs[2], roc_b1$aucs[2]),
  RSF_AUC_5yr = c(rsf_auc$AUC_5yr[rsf_auc$Dataset == "TCGA (Training)"],
                   rsf_auc$AUC_5yr[rsf_auc$Dataset == "CGGA_batch2 (Tuning)"],
                   rsf_auc$AUC_5yr[rsf_auc$Dataset == "CGGA_batch1 (Blind Test)"]),
  CoxPH_AUC_5yr = c(roc_train$aucs[3], roc_b2$aucs[3], roc_b1$aucs[3]),
  stringsAsFactors = FALSE
)

print(comparison_df, digits = 4)
write.csv(comparison_df, "results/coxph_vs_rsf_comparison.csv", row.names = FALSE)
cat("\n  Saved: results/coxph_vs_rsf_comparison.csv\n")

# =============================================================================
# Output Risk Formula
# =============================================================================
cat("\n========================================\n")
cat("Clinical Risk Score Formula\n")
cat("========================================\n\n")

formula_parts <- sprintf("%.4f * %s", betas, names(betas))
risk_formula <- paste("Risk =", paste(formula_parts, collapse = "\n     + "))
cat(risk_formula, "\n\n")

writeLines(c(
  "UPR-based Interpretable Risk Score Formula (CoxPH with IDH stratification)",
  "==========================================================================",
  "",
  "Model: Cox Proportional Hazards with strata(IDH_status)",
  paste("Number of genes:", length(betas)),
  paste("Training set: TCGA glioma (n =", nrow(train_df), ")"),
  "",
  "Risk Score = Linear Predictor = sum(beta_i * log2(expression_i + 1))",
  "",
  risk_formula,
  "",
  "Note: Expression values should be log2(TPM+1) or log2(FPKM+1) normalized.",
  "Patients are stratified by IDH mutation status.",
  paste("Risk group cutoff: use cohort-specific median of risk scores"),
  paste("Training set median:", round(train_median, 4)),
  "  Risk >= median => High Risk",
  "  Risk <  median => Low Risk"
), "results/coxph_risk_formula.txt")
cat("  Saved: results/coxph_risk_formula.txt\n")

# =============================================================================
# Generate PDF Figures
# =============================================================================
cat("\n[Plotting] Generating KM and ROC figures...\n")

# --- KM helper ---
plot_km <- function(df, title_str) {
  df$risk_group <- factor(df$risk_group, levels = c("Low", "High"))
  stopifnot(length(unique(df$risk_group)) == 2)

  fit <- survfit(Surv(time, status) ~ risk_group, data = df)

  p <- ggsurvplot(
    fit, data = df,
    pval = TRUE,
    risk.table = TRUE,
    risk.table.col = "strata",
    palette = c("#2E9FDF", "#E7B800"),
    title = title_str,
    xlab = "Time (days)",
    ylab = "Overall Survival",
    legend.title = "Risk Group",
    legend.labs = c("Low Risk", "High Risk"),
    conf.int = FALSE,
    ggtheme = theme_bw(base_size = 12)
  )
  return(p)
}

# --- Figure 1: KM Survival Curves ---
# Use pure ggplot2 approach to avoid survminer's blank page issue
# Build KM plots as standard ggplot objects without risk table
plot_km_gg <- function(df, title_str) {
  df$risk_group <- factor(df$risk_group, levels = c("Low", "High"))
  fit <- survfit(Surv(time, status) ~ risk_group, data = df)

  # Log-rank p-value
  lr <- survdiff(Surv(time, status) ~ risk_group, data = df)
  lr_p <- pchisq(lr$chisq, df = 1, lower.tail = FALSE)
  p_label <- ifelse(lr_p < 0.0001, "p < 0.0001", paste0("p = ", formatC(lr_p, format = "e", digits = 2)))

  # Extract data from survfit for ggplot
  sdata <- data.frame(
    time = fit$time,
    surv = fit$surv,
    upper = fit$upper,
    lower = fit$lower,
    strata = rep(names(fit$strata), fit$strata)
  )
  sdata$group <- ifelse(grepl("Low", sdata$strata), "Low Risk", "High Risk")
  sdata$group <- factor(sdata$group, levels = c("Low Risk", "High Risk"))

  # Add starting point (0, 1) for each group
  start_df <- data.frame(
    time = c(0, 0), surv = c(1, 1), upper = c(1, 1), lower = c(1, 1),
    strata = c("a", "b"),
    group = factor(c("Low Risk", "High Risk"), levels = c("Low Risk", "High Risk"))
  )
  sdata <- rbind(start_df, sdata)

  # Number at risk at specific time points
  time_breaks <- pretty(c(0, max(df$time)), n = 5)
  n_risk_low <- sapply(time_breaks, function(t) sum(df$time[df$risk_group == "Low"] >= t))
  n_risk_high <- sapply(time_breaks, function(t) sum(df$time[df$risk_group == "High"] >= t))
  risk_label <- paste0(
    "Low Risk: ", paste(n_risk_low, collapse = "  "),
    "\nHigh Risk: ", paste(n_risk_high, collapse = "  ")
  )

  gg <- ggplot(sdata, aes(x = time, y = surv, color = group)) +
    geom_step(linewidth = 1) +
    scale_color_manual(values = c("Low Risk" = "#2E9FDF", "High Risk" = "#E7B800")) +
    labs(title = title_str, x = "Time (days)", y = "Overall Survival",
         color = "Risk Group") +
    annotate("text", x = max(df$time) * 0.05, y = 0.15,
             label = p_label, hjust = 0, size = 4.5, fontface = "bold") +
    coord_cartesian(ylim = c(0, 1)) +
    theme_bw(base_size = 13) +
    theme(legend.position = c(0.85, 0.85),
          legend.background = element_rect(fill = alpha("white", 0.8)),
          plot.title = element_text(size = 14, face = "bold"))

  return(gg)
}

g1 <- plot_km_gg(train_df, paste0("TCGA Training (n=", nrow(train_df), "), C-index=", round(cindex_train, 3)))
g2 <- plot_km_gg(cgga1_val_df, paste0("CGGA Batch1 Blind Test (n=", nrow(cgga1_val_df), "), C-index=", round(cindex_b1, 3)))
g3 <- plot_km_gg(cgga2_val_df, paste0("CGGA Batch2 Tuning (n=", nrow(cgga2_val_df), "), C-index=", round(cindex_b2, 3)))

pdf("figures/coxph_km_survival.pdf", width = 10, height = 7)
print(g1)
print(g2)
print(g3)
dev.off()
cat("  Saved: figures/coxph_km_survival.pdf\n")

# --- Figure 2: timeROC Curves ---
pdf("figures/coxph_timeroc.pdf", width = 15, height = 5)
par(mfrow = c(1, 3))

plot_roc <- function(roc_obj, aucs, title_str) {
  colors <- c("#E41A1C", "#377EB8", "#4DAF4A")

  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
       xlab = "1 - Specificity", ylab = "Sensitivity", main = title_str)
  abline(0, 1, lty = 2, col = "gray60")

  for (i in 1:3) {
    lines(roc_obj$FP[, i], roc_obj$TP[, i], col = colors[i], lwd = 2)
  }

  legend("bottomright",
         legend = c(paste0("1-yr AUC=", round(aucs[1], 3)),
                    paste0("3-yr AUC=", round(aucs[2], 3)),
                    paste0("5-yr AUC=", round(aucs[3], 3))),
         col = colors, lwd = 2, bty = "n", cex = 0.9)
}

plot_roc(roc_train$roc, roc_train$aucs,
         paste0("TCGA Training (n=", nrow(train_df), ")"))
plot_roc(roc_b1$roc, roc_b1$aucs,
         paste0("CGGA Batch1 (n=", nrow(cgga1_val_df), ")"))
plot_roc(roc_b2$roc, roc_b2$aucs,
         paste0("CGGA Batch2 (n=", nrow(cgga2_val_df), ")"))

dev.off()
cat("  Saved: figures/coxph_timeroc.pdf\n")

# --- Figure 3: PH Assumption Diagnostic ---
pdf("figures/coxph_ph_diagnostic.pdf", width = 12, height = 8)
plot(ph_test)
dev.off()
cat("  Saved: figures/coxph_ph_diagnostic.pdf\n")

# --- Figure 4: Forest Plot of Coefficients ---
pdf("figures/coxph_forest_plot.pdf", width = 8, height = max(4, 0.5 * nrow(coef_table) + 2))

ct <- coef_table[order(coef_table$Beta), ]
n_genes <- nrow(ct)
y_pos <- seq_len(n_genes)

par(mar = c(5, 8, 4, 2))
plot(ct$HR, y_pos, xlim = range(c(ct$HR_lower, ct$HR_upper, 0.5, 2)),
     pch = 18, cex = 1.5, xlab = "Hazard Ratio (95% CI)", ylab = "",
     yaxt = "n", main = "CoxPH Model Coefficients", log = "x")

segments(ct$HR_lower, y_pos, ct$HR_upper, y_pos, lwd = 2)
abline(v = 1, lty = 2, col = "red")
axis(2, at = y_pos, labels = ct$Gene, las = 1, cex.axis = 0.9)

for (i in seq_len(n_genes)) {
  p_text <- ifelse(ct$p_value[i] < 0.001, "***",
             ifelse(ct$p_value[i] < 0.01, "**",
              ifelse(ct$p_value[i] < 0.05, "*", "ns")))
  text(ct$HR_upper[i] * 1.1, y_pos[i], p_text, cex = 0.8, col = "darkred")
}

dev.off()
cat("  Saved: figures/coxph_forest_plot.pdf\n")

# =============================================================================
# Save All Results
# =============================================================================
cat("\n[Saving] Packaging results...\n")

coxph_results <- list(
  cox_model = cox_model,
  model_genes = model_genes_final,
  betas = betas,
  cox_formula = cox_formula,
  risk_formula_text = risk_formula,

  cv_lasso = cv_lasso,
  lasso_coef = lasso_coef_vec,
  selected_genes = selected_genes,
  candidate_genes = candidate_genes,

  ph_test = ph_test,
  ph_results = ph_results,
  coef_table = coef_table,

  train_risk = train_risk,
  train_median_cutoff = train_median,

  train_df = train_df,
  cgga1_val_df = cgga1_val_df,
  cgga2_val_df = cgga2_val_df,

  cindex = c(TCGA_train = cindex_train,
             CGGA_batch1 = cindex_b1,
             CGGA_batch2 = cindex_b2),

  roc_train = roc_train,
  roc_b1 = roc_b1,
  roc_b2 = roc_b2,

  comparison_df = comparison_df
)

save(coxph_results, file = "data/processed/coxph_interpretable_model.RData")
cat("  Saved: data/processed/coxph_interpretable_model.RData\n")

# =============================================================================
# Final Summary
# =============================================================================
cat("\n")
cat("==========================================\n")
cat("        FINAL SUMMARY                     \n")
cat("==========================================\n")
cat("Model: CoxPH with strata(IDH_status)\n")
cat("Genes: ", length(model_genes_final), " (from LASSO selection of ",
    length(candidate_genes), " candidates)\n", sep = "")
cat("Selected:", paste(model_genes_final, collapse = ", "), "\n\n")

cat("Performance Comparison (C-index):\n")
cat(sprintf("%-25s %12s %12s\n", "Dataset", "LASSO+RSF", "LASSO+CoxPH"))
cat(sprintf("%-25s %12s %12s\n", "------------------------", "-----------", "-----------"))
for (i in seq_len(nrow(comparison_df))) {
  cat(sprintf("%-25s %12.4f %12.4f\n",
              comparison_df$Dataset[i],
              comparison_df$RSF_Cindex[i],
              comparison_df$CoxPH_Cindex[i]))
}

cat("\ntimeROC AUC (CoxPH):\n")
cat(sprintf("%-25s %8s %8s %8s\n", "Dataset", "1-yr", "3-yr", "5-yr"))
cat(sprintf("%-25s %8.3f %8.3f %8.3f\n", "TCGA Training",
            roc_train$aucs[1], roc_train$aucs[2], roc_train$aucs[3]))
cat(sprintf("%-25s %8.3f %8.3f %8.3f\n", "CGGA Batch1",
            roc_b1$aucs[1], roc_b1$aucs[2], roc_b1$aucs[3]))
cat(sprintf("%-25s %8.3f %8.3f %8.3f\n", "CGGA Batch2",
            roc_b2$aucs[1], roc_b2$aucs[2], roc_b2$aucs[3]))

cat("\nPH Assumption: ",
    sum(ph_results$PH_satisfied[ph_results$variable != "GLOBAL"] == "Yes"), "/",
    sum(ph_results$variable != "GLOBAL"), " genes passed (p > 0.05)\n", sep = "")
cat("Global PH test p =",
    round(ph_results$p_value[ph_results$variable == "GLOBAL"], 4), "\n")

cat("\nOutput files:\n")
cat("  data/processed/coxph_interpretable_model.RData\n")
cat("  results/coxph_model_coefficients.csv\n")
cat("  results/coxph_vs_rsf_comparison.csv\n")
cat("  results/coxph_risk_formula.txt\n")
cat("  figures/coxph_km_survival.pdf\n")
cat("  figures/coxph_timeroc.pdf\n")
cat("  figures/coxph_ph_diagnostic.pdf\n")
cat("  figures/coxph_forest_plot.pdf\n")
cat("\nDone!\n")
