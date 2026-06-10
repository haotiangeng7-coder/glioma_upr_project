###############################################################################
# 01_feature_selection.R
# 特征筛选：单因素Cox + PH假设检验 + RMST替代 + LASSO Cox回归
# 评审修订版：加入cox.zph() PH检验，违反PH的基因用survRM2::rmst2()替代筛选
###############################################################################

source("00_setup/config.R")
library(survival)
library(glmnet)
library(survRM2)
library(ggplot2)
library(dplyr)

set.seed(SEED)
load(file.path(DATA_PROC, "upr_gene_sets.RData"))
load(file.path(DATA_PROC, "tcga_glioma_expression.RData"))
load(file.path(DATA_PROC, "upr_landscape_results.RData"))

# =============================================================================
# 1. 准备训练数据
# =============================================================================
message("=== Preparing training data ===")

# UPR基因在数据中的交集
upr_in_data <- UPR_broad_genes[UPR_broad_genes %in% rownames(expr_tpm_symbol)]

# 表达矩阵（log2 TPM+1）
common_samples <- intersect(colnames(expr_tpm_symbol), clinical_valid$barcode)
expr_train <- log2(expr_tpm_symbol[upr_in_data, common_samples] + 1)

# 匹配临床数据
clin_train <- clinical_valid[match(common_samples, clinical_valid$barcode), ]

# 构建生存数据框
surv_df <- data.frame(
  time   = clin_train$OS.time,
  status = clin_train$OS,
  t(expr_train),
  check.names = FALSE
)

# 过滤有效样本
surv_df <- surv_df[complete.cases(surv_df$time, surv_df$status), ]
surv_df <- surv_df[surv_df$time > 0, ]

message(sprintf("Training set: %d samples, %d UPR genes",
                nrow(surv_df), ncol(surv_df) - 2))

# =============================================================================
# 2. 单因素Cox回归筛选 + PH假设检验
# =============================================================================
message("=== Univariate Cox regression with PH assumption testing ===")

gene_names <- colnames(surv_df)[-(1:2)]

# 预分配list避免反复rbind
uni_list <- vector("list", length(gene_names))

for (i in seq_along(gene_names)) {
  gene <- gene_names[i]
  tryCatch({
    fml <- as.formula(paste0("Surv(time, status) ~ `", gene, "`"))
    fit <- coxph(fml, data = surv_df)
    s   <- summary(fit)

    # PH假设检验
    ph_test <- cox.zph(fit)
    ph_p    <- ph_test$table[1, "p"]

    uni_list[[i]] <- data.frame(
      gene      = gene,
      HR        = s$conf.int[1, 1],
      HR_lower  = s$conf.int[1, 3],
      HR_upper  = s$conf.int[1, 4],
      coef      = s$coefficients[1, 1],
      pvalue    = s$coefficients[1, 5],
      ph_pvalue = ph_p,
      ph_passed = (ph_p >= 0.05),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    message(sprintf("  Skipping %s: %s", gene, e$message))
  })
}

univariate_results <- do.call(rbind, uni_list[!sapply(uni_list, is.null)])
univariate_results$padj <- p.adjust(univariate_results$pvalue, method = "BH")
univariate_results <- univariate_results %>% dplyr::arrange(pvalue)

n_ph_fail <- sum(!univariate_results$ph_passed, na.rm = TRUE)
message(sprintf("PH assumption: %d/%d genes violated (p<0.05)",
                n_ph_fail, nrow(univariate_results)))

# =============================================================================
# 3. 对违反PH假设的基因使用RMST差异替代筛选
# =============================================================================
message("=== RMST analysis for PH-violating genes ===")

# 确定RMST截断时间 tau
tau_rmst <- min(
  quantile(surv_df$time[surv_df$status == 1], 0.9, na.rm = TRUE),
  quantile(surv_df$time, 0.95, na.rm = TRUE)
)
message(sprintf("RMST tau (truncation time): %.1f days", tau_rmst))

# PH通过且Cox p<0.05
ph_pass_genes <- univariate_results %>%
  dplyr::filter(ph_passed & pvalue < PVALUE_CUTOFF) %>%
  dplyr::pull(gene)

# PH不通过
ph_fail_genes <- univariate_results %>%
  dplyr::filter(!ph_passed) %>%
  dplyr::pull(gene)

rmst_list <- vector("list", length(ph_fail_genes))

if (length(ph_fail_genes) > 0) {
  for (j in seq_along(ph_fail_genes)) {
    gene <- ph_fail_genes[j]
    tryCatch({
      expr_vec <- surv_df[[gene]]
      med_expr <- median(expr_vec, na.rm = TRUE)
      arm_indicator <- as.integer(expr_vec > med_expr)

      if (sum(arm_indicator == 0) >= 10 && sum(arm_indicator == 1) >= 10) {
        rmst_fit <- rmst2(
          time   = surv_df$time,
          status = surv_df$status,
          arm    = arm_indicator,
          tau    = tau_rmst
        )
        rmst_list[[j]] <- data.frame(
          gene        = gene,
          rmst_diff   = rmst_fit$unadjusted.result[1, "Est."],
          rmst_pvalue = rmst_fit$unadjusted.result[1, "p"],
          stringsAsFactors = FALSE
        )
      }
    }, error = function(e) {
      message(sprintf("  RMST failed for %s: %s", gene, e$message))
    })
  }
}

rmst_results <- do.call(rbind, rmst_list[!sapply(rmst_list, is.null)])
if (is.null(rmst_results)) rmst_results <- data.frame(
  gene = character(0), rmst_diff = numeric(0), rmst_pvalue = numeric(0)
)

rmst_sig_genes <- character(0)
if (nrow(rmst_results) > 0) {
  rmst_sig_genes <- rmst_results %>%
    dplyr::filter(rmst_pvalue < PVALUE_CUTOFF) %>%
    dplyr::pull(gene)
  message(sprintf("RMST significant genes (PH-violating): %d/%d",
                  length(rmst_sig_genes), length(ph_fail_genes)))
}

# 合并筛选结果
sig_genes_uni <- unique(c(ph_pass_genes, rmst_sig_genes))
message(sprintf("Total significant genes (Cox + RMST): %d", length(sig_genes_uni)))

# 在univariate_results中标注筛选方法
if (nrow(rmst_results) > 0) {
  univariate_results <- univariate_results %>%
    dplyr::left_join(rmst_results, by = "gene") %>%
    dplyr::mutate(
      screening_method = dplyr::case_when(
        ph_passed & pvalue < PVALUE_CUTOFF       ~ "Cox_PH_pass",
        !ph_passed & gene %in% rmst_sig_genes    ~ "RMST_alternative",
        TRUE                                       ~ "not_selected"
      )
    )
} else {
  univariate_results$rmst_diff      <- NA_real_
  univariate_results$rmst_pvalue    <- NA_real_
  univariate_results$screening_method <- ifelse(
    univariate_results$ph_passed & univariate_results$pvalue < PVALUE_CUTOFF,
    "Cox_PH_pass", "not_selected"
  )
}

write.csv(univariate_results,
          file.path(RES_DIR, "univariate_cox_results.csv"), row.names = FALSE)

# =============================================================================
# 4. LASSO Cox回归进一步筛选
# =============================================================================
message("=== LASSO Cox regression ===")

if (length(sig_genes_uni) < 5) {
  message("Too few significant genes. Using top 50 by p-value instead.")
  sig_genes_uni <- univariate_results$gene[1:min(50, nrow(univariate_results))]
}

x <- as.matrix(surv_df[, sig_genes_uni, drop = FALSE])
y <- Surv(surv_df$time, surv_df$status)

set.seed(SEED)
cv_fit <- cv.glmnet(
  x = x, y = y,
  family       = "cox",
  alpha        = 1,
  nfolds       = ML_NFOLDS,
  type.measure = "C"
)

# --- CV曲线 ---
pdf(file.path(FIG_DIR, "Fig5_lasso_cv.pdf"), width = 8, height = 6)
plot(cv_fit, main = "")
title(main = "LASSO Cox Regression: Cross-Validation",
      sub  = sprintf("lambda.min = %.4f, lambda.1se = %.4f",
                      cv_fit$lambda.min, cv_fit$lambda.1se))
dev.off()

# lambda.min
lasso_coefs <- coef(cv_fit, s = "lambda.min")
lasso_genes <- rownames(lasso_coefs)[lasso_coefs[, 1] != 0]
message(sprintf("LASSO selected genes (lambda.min): %d", length(lasso_genes)))
message("Selected genes: ", paste(lasso_genes, collapse = ", "))

# lambda.1se
lasso_coefs_1se <- coef(cv_fit, s = "lambda.1se")
lasso_genes_1se <- rownames(lasso_coefs_1se)[lasso_coefs_1se[, 1] != 0]
message(sprintf("LASSO selected genes (lambda.1se): %d", length(lasso_genes_1se)))

candidate_genes <- lasso_genes

if (length(candidate_genes) < 5) {
  message("Too few genes with lambda.min, relaxing threshold...")
  mid_lambda <- cv_fit$lambda[length(cv_fit$lambda) %/% 2]
  lasso_coefs_mid <- coef(cv_fit, s = mid_lambda)
  candidate_genes <- rownames(lasso_coefs_mid)[lasso_coefs_mid[, 1] != 0]
}

# =============================================================================
# 5. LASSO系数可视化
# =============================================================================
message("=== Visualizing LASSO coefficients ===")

coef_df <- data.frame(
  gene        = candidate_genes,
  coefficient = as.numeric(lasso_coefs[candidate_genes, 1]),
  stringsAsFactors = FALSE
) %>% dplyr::arrange(coefficient)

p_coef <- ggplot(coef_df, aes(x = coefficient, y = reorder(gene, coefficient))) +
  geom_bar(stat = "identity",
           aes(fill = ifelse(coefficient > 0, "Risk", "Protective"))) +
  scale_fill_manual(values = c("Risk" = "#E64B35", "Protective" = "#4DBBD5")) +
  labs(x = "LASSO Coefficient", y = "", fill = "") +
  THEME_PUBLICATION +
  ggtitle("LASSO Cox: Selected Gene Coefficients")

ggsave(file.path(FIG_DIR, "Fig5_lasso_coefficients.pdf"), p_coef,
       width = 8, height = max(4, length(candidate_genes) * 0.35))

# =============================================================================
# 6. LASSO路径图
# =============================================================================
message("=== LASSO coefficient path ===")

fit_path <- glmnet(x, y, family = "cox", alpha = 1)

pdf(file.path(FIG_DIR, "Fig5_lasso_path.pdf"), width = 8, height = 6)
plot(fit_path, xvar = "lambda", label = TRUE)
abline(v = log(cv_fit$lambda.min), lty = 2, col = "red")
abline(v = log(cv_fit$lambda.1se), lty = 2, col = "blue")
legend("topright",
       legend = c("lambda.min", "lambda.1se"),
       col = c("red", "blue"), lty = 2, bty = "n")
title("LASSO Coefficient Path")
dev.off()

# =============================================================================
# 7. PH检验汇总表（Supplementary Table）
# =============================================================================
message("=== Saving PH assumption test summary ===")

ph_summary <- univariate_results %>%
  dplyr::select(gene, HR, pvalue, padj, ph_pvalue, ph_passed,
                rmst_diff, rmst_pvalue, screening_method) %>%
  dplyr::arrange(ph_pvalue)

write.csv(ph_summary,
          file.path(RES_DIR, "supplementary_ph_assumption_test.csv"),
          row.names = FALSE)

# =============================================================================
# 8. 保存特征筛选结果
# =============================================================================
save(
  univariate_results,
  sig_genes_uni,
  rmst_results,
  cv_fit,
  candidate_genes,
  lasso_coefs,
  surv_df,
  tau_rmst,
  file = file.path(DATA_PROC, "feature_selection_results.RData")
)

message(sprintf("\n=== Feature selection completed ==="))
message(sprintf("Candidate genes for ML: %d", length(candidate_genes)))
message("Genes: ", paste(candidate_genes, collapse = ", "))
message(sprintf("PH violations handled: %d genes screened via RMST",
                length(ph_fail_genes)))
message("\nNext step: Run 05_ml_model/02_ml_combinations.R")
