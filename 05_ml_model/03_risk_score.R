###############################################################################
# 03_risk_score.R
# 基于最优模型计算UPR-Immune Risk Score (UIRS)
# 训练集嵌套CV + 调优验证 + 盲测试集 K-M
# timeROC 1/3/5年 + 风险评分分布图
# GSE43378探索性验证置于Supplementary
# 输出: Figure 5C-H
###############################################################################

source("00_setup/config.R")
library(survival)
library(survminer)
library(timeROC)
library(glmnet)
library(randomForestSRC)
library(CoxBoost)
library(superpc)
library(plsRcox)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(dplyr)

set.seed(SEED)
load(file.path(DATA_PROC, "ml_combination_results.RData"))
load(file.path(DATA_PROC, "feature_selection_results.RData"))
load(file.path(DATA_PROC, "tcga_glioma_expression.RData"))

# =============================================================================
# 1. 使用最优算法组合构建最终模型
# =============================================================================
message("=== Building final risk model using best combination ===")
message(sprintf("Best combination: %s", best_combo))
message(sprintf("  Feature selection: %s", best_fs))
message(sprintf("  Model builder: %s", best_build))

# 在全训练集上训练最终模型
x_full <- as.matrix(train_expr[, genes])
y_full <- Surv(train_expr$time, train_expr$status)

# 特征选择
fs_final <- FS_ALGOS[[best_fs]](x_full, y_full, seed = SEED)
model_genes <- intersect(fs_final$selected, colnames(x_full))
if (length(model_genes) < 2) model_genes <- colnames(x_full)

message(sprintf("Final model genes: %d", length(model_genes)))
message("  Genes: ", paste(model_genes, collapse = ", "))

# 模型构建
final_built <- BUILD_ALGOS[[best_build]](
  x_full[, model_genes, drop = FALSE], y_full, seed = SEED
)

# 提取模型权重（如果可用）
model_weights <- tryCatch({
  if (best_build %in% c("LASSO", "Ridge", "ElasticNet03", "ElasticNet05", "ElasticNet07")) {
    cv_obj <- final_built$model
    coefs  <- coef(cv_obj, s = "lambda.min")
    setNames(as.numeric(coefs[model_genes, 1]), model_genes)
  } else {
    NULL
  }
}, error = function(e) NULL)

if (!is.null(model_weights)) {
  message("\nModel coefficients:")
  for (i in seq_along(model_weights)) {
    message(sprintf("  %s: %.4f", names(model_weights)[i], model_weights[i]))
  }
}

# =============================================================================
# 2. 计算Risk Score
# =============================================================================
message("=== Computing Risk Scores ===")

# 训练集
risk_score_train <- final_built$predict_fn(
  data.frame(train_expr[, model_genes, drop = FALSE], check.names = FALSE)
)
names(risk_score_train) <- rownames(train_expr)

# 使用surv_cutpoint找最优截断值
cutpoint_data <- data.frame(
  time       = train_expr$time,
  status     = train_expr$status,
  risk_score = risk_score_train
)

cutpoint    <- surv_cutpoint(cutpoint_data, time = "time", event = "status",
                              variables = "risk_score")
optimal_cut <- cutpoint$cutpoint$cutpoint
message(sprintf("Optimal cutpoint: %.4f", optimal_cut))

# 中位数截断值（更通用）
median_cut <- median(risk_score_train, na.rm = TRUE)
message(sprintf("Median cutpoint: %.4f", median_cut))

# 使用中位数分组
risk_group_train <- ifelse(risk_score_train > median_cut, "High", "Low")

# =============================================================================
# 3. 训练集K-M生存分析
# =============================================================================
message("=== Training set K-M analysis ===")

surv_train <- data.frame(
  time       = train_expr$time / 30,
  status     = train_expr$status,
  risk_score = risk_score_train,
  risk_group = factor(risk_group_train, levels = c("Low", "High"))
)

fit_train <- survfit(Surv(time, status) ~ risk_group, data = surv_train)

p_km_train <- ggsurvplot(
  fit_train,
  data           = surv_train,
  pval           = TRUE,
  pval.method    = TRUE,
  risk.table     = TRUE,
  palette        = c(COLORS_RISK["Low"], COLORS_RISK["High"]),
  xlab           = "Time (months)",
  ylab           = "Overall Survival",
  title          = "Training Set (TCGA): UIRS Risk Stratification",
  ggtheme        = THEME_PUBLICATION,
  risk.table.height = 0.25,
  conf.int       = TRUE
)

pdf(file.path(FIG_DIR, "Fig5C_km_train.pdf"), width = 8, height = 7)
print(p_km_train)
dev.off()

# =============================================================================
# 4. 验证集K-M生存分析（调优验证 + 盲测试）
# =============================================================================
message("=== Validation set K-M analysis ===")

# 用于汇总的结果列表
validation_results <- list()

# 辅助函数：对一个数据集做KM + C-index
run_validation_km <- function(val_df, val_name, fig_prefix, fig_label) {
  val_genes_avail <- intersect(model_genes, colnames(val_df)[-(1:2)])

  if (length(val_genes_avail) < length(model_genes) * 0.7) {
    message(sprintf("  %s: insufficient gene overlap (%d/%d), skipping",
                    val_name, length(val_genes_avail), length(model_genes)))
    return(NULL)
  }

  # 对缺失基因填0
  val_mat <- matrix(0, nrow = nrow(val_df), ncol = length(model_genes),
                    dimnames = list(rownames(val_df), model_genes))
  for (g in val_genes_avail) val_mat[, g] <- val_df[[g]]

  risk_score_val <- final_built$predict_fn(
    data.frame(val_mat, check.names = FALSE)
  )

  val_median <- median(risk_score_val, na.rm = TRUE)
  risk_group_val <- ifelse(risk_score_val > val_median, "High", "Low")

  surv_val <- data.frame(
    time       = val_df$time / 30,
    status     = val_df$status,
    risk_score = risk_score_val,
    risk_group = factor(risk_group_val, levels = c("Low", "High"))
  )

  # K-M曲线
  fit_val <- survfit(Surv(time, status) ~ risk_group, data = surv_val)
  p_km_val <- ggsurvplot(
    fit_val,
    data           = surv_val,
    pval           = TRUE,
    pval.method    = TRUE,
    risk.table     = TRUE,
    palette        = c(COLORS_RISK["Low"], COLORS_RISK["High"]),
    xlab           = "Time (months)",
    ylab           = "Overall Survival",
    title          = paste0(fig_label, ": ", val_name),
    ggtheme        = THEME_PUBLICATION,
    risk.table.height = 0.25
  )

  pdf(file.path(FIG_DIR, paste0(fig_prefix, "_km_", val_name, ".pdf")),
      width = 8, height = 7)
  print(p_km_val)
  dev.off()

  c_val <- calc_cindex(risk_score_val, val_df$time, val_df$status)
  message(sprintf("  %s C-index: %.4f", val_name, c_val))

  return(list(
    risk_score = risk_score_val,
    risk_group = risk_group_val,
    c_index    = c_val,
    surv_data  = surv_val
  ))
}

# C-index辅助函数
calc_cindex <- function(pred_risk, time, status) {
  tryCatch({
    fit_tmp <- coxph(Surv(time, status) ~ pred_risk)
    concordance(fit_tmp)$concordance
  }, error = function(e) NA_real_)
}

# 调优验证集
for (tn in names(tuning_data)) {
  message(sprintf("\n--- Tuning validation: %s ---", tn))
  res <- run_validation_km(tuning_data[[tn]], tn, "Fig5D", "Tuning Validation")
  if (!is.null(res)) validation_results[[tn]] <- res
}

# 盲测试集
blind_fig_labels <- c("Fig5E", "Fig5F")
bi <- 1
for (bt_name in names(blind_test)) {
  message(sprintf("\n--- Blind test: %s ---", bt_name))
  fl <- ifelse(bi <= length(blind_fig_labels), blind_fig_labels[bi], "Fig5")
  res <- run_validation_km(blind_test[[bt_name]], bt_name, fl, "Blind Test")
  if (!is.null(res)) validation_results[[bt_name]] <- res
  bi <- bi + 1
}

# 探索性（GSE43378）→ Supplementary
for (ex_name in names(exploratory)) {
  message(sprintf("\n--- Exploratory: %s ---", ex_name))
  res <- run_validation_km(exploratory[[ex_name]], ex_name,
                           "SuppFig", "Exploratory Validation")
  if (!is.null(res)) validation_results[[paste0(ex_name, "_exploratory")]] <- res
}

# =============================================================================
# 5. 时间依赖ROC (timeROC)
# =============================================================================
message("\n=== Time-dependent ROC analysis ===")

plot_timeROC <- function(surv_data, dataset_name, fig_prefix) {
  tryCatch({
    # 过滤有效数据
    valid <- complete.cases(surv_data$time, surv_data$status, surv_data$risk_score) &
             surv_data$time > 0
    sv <- surv_data[valid, ]

    troc <- timeROC(
      T     = sv$time,
      delta = sv$status,
      marker = sv$risk_score,
      cause = 1,
      times = c(12, 36, 60),
      iid   = TRUE
    )

    auc_vals <- round(troc$AUC, 3)
    message(sprintf("  %s AUC: 1y=%.3f, 3y=%.3f, 5y=%.3f",
                    dataset_name, auc_vals[1], auc_vals[2], auc_vals[3]))

    pdf(file.path(FIG_DIR, paste0(fig_prefix, "_timeROC_", dataset_name, ".pdf")),
        width = 7, height = 6)
    plot(troc, time = 12, col = "#E64B35", lwd = 2, title = FALSE)
    plot(troc, time = 36, col = "#4DBBD5", lwd = 2, add = TRUE)
    plot(troc, time = 60, col = "#00A087", lwd = 2, add = TRUE)
    legend("bottomright",
           legend = c(
             sprintf("1-year AUC = %.3f", auc_vals[1]),
             sprintf("3-year AUC = %.3f", auc_vals[2]),
             sprintf("5-year AUC = %.3f", auc_vals[3])
           ),
           col = c("#E64B35", "#4DBBD5", "#00A087"),
           lwd = 2, bty = "n")
    title(sprintf("Time-dependent ROC (%s)", dataset_name))
    dev.off()

    return(list(troc = troc, auc = auc_vals))
  }, error = function(e) {
    message(sprintf("  timeROC failed for %s: %s", dataset_name, e$message))
    return(NULL)
  })
}

# 训练集
roc_train <- plot_timeROC(surv_train, "TCGA_train", "Fig5G")

# 验证集
roc_val_results <- list()
for (vn in names(validation_results)) {
  prefix <- ifelse(grepl("exploratory", vn), "SuppFig", "Fig5G")
  roc_val_results[[vn]] <- plot_timeROC(validation_results[[vn]]$surv_data, vn, prefix)
}

# =============================================================================
# 6. 风险评分分布图 + 生存状态散点图
# =============================================================================
message("=== Risk score distribution ===")

plot_risk_distribution <- function(surv_data, dataset_name, fig_prefix) {
  sorted <- surv_data[order(surv_data$risk_score), ]
  sorted$rank <- seq_len(nrow(sorted))

  # 风险评分点图
  p1 <- ggplot(sorted, aes(x = rank, y = risk_score, color = risk_group)) +
    geom_point(size = 0.8) +
    scale_color_manual(values = COLORS_RISK) +
    geom_hline(yintercept = median_cut, linetype = "dashed", color = "black") +
    labs(x = "Patient (ranked by risk score)", y = "UIRS Risk Score",
         color = "Risk Group") +
    THEME_PUBLICATION +
    ggtitle(paste0(dataset_name, ": Risk Score Distribution"))

  # 生存状态散点图
  p2 <- ggplot(sorted, aes(x = rank, y = time,
                            color = factor(status, levels = c(0, 1)))) +
    geom_point(size = 0.8) +
    scale_color_manual(values = c("0" = "#4DBBD5", "1" = "#E64B35"),
                       labels = c("Alive", "Dead")) +
    labs(x = "Patient (ranked by risk score)", y = "Survival Time (months)",
         color = "Status") +
    THEME_PUBLICATION

  p_combined <- cowplot::plot_grid(p1, p2, ncol = 1, align = "v", rel_heights = c(1, 1))
  ggsave(file.path(FIG_DIR, paste0(fig_prefix, "_risk_dist_", dataset_name, ".pdf")),
         p_combined, width = 10, height = 8)
}

# 训练集
plot_risk_distribution(surv_train, "TCGA_train", "Fig5H")

# 验证集
for (vn in names(validation_results)) {
  prefix <- ifelse(grepl("exploratory", vn), "SuppFig", "Fig5H")
  plot_risk_distribution(validation_results[[vn]]$surv_data, vn, prefix)
}

# =============================================================================
# 7. 汇总AUC表
# =============================================================================
message("=== Summarizing AUC results ===")

auc_summary <- data.frame(
  Dataset  = character(0),
  AUC_1yr  = numeric(0),
  AUC_3yr  = numeric(0),
  AUC_5yr  = numeric(0),
  C_index  = numeric(0),
  stringsAsFactors = FALSE
)

# 训练集
if (!is.null(roc_train)) {
  auc_summary <- rbind(auc_summary, data.frame(
    Dataset = "TCGA (Training)",
    AUC_1yr = roc_train$auc[1],
    AUC_3yr = roc_train$auc[2],
    AUC_5yr = roc_train$auc[3],
    C_index = calc_cindex(risk_score_train, train_expr$time, train_expr$status)
  ))
}

# 验证集
for (vn in names(validation_results)) {
  role_label <- dplyr::case_when(
    vn %in% names(tuning_data)  ~ "Tuning",
    vn %in% names(blind_test)   ~ "Blind Test",
    TRUE                         ~ "Exploratory"
  )
  auc_row <- data.frame(
    Dataset = sprintf("%s (%s)", vn, role_label),
    AUC_1yr = NA_real_, AUC_3yr = NA_real_, AUC_5yr = NA_real_,
    C_index = validation_results[[vn]]$c_index
  )
  if (vn %in% names(roc_val_results) && !is.null(roc_val_results[[vn]])) {
    auc_row$AUC_1yr <- roc_val_results[[vn]]$auc[1]
    auc_row$AUC_3yr <- roc_val_results[[vn]]$auc[2]
    auc_row$AUC_5yr <- roc_val_results[[vn]]$auc[3]
  }
  auc_summary <- rbind(auc_summary, auc_row)
}

write.csv(auc_summary,
          file.path(RES_DIR, "uirs_auc_summary.csv"), row.names = FALSE)
message("\nAUC summary:")
print(auc_summary)

# =============================================================================
# 8. 保存最终模型和结果
# =============================================================================
save(
  final_built,
  model_genes, model_weights,
  risk_score_train, risk_group_train,
  median_cut, optimal_cut,
  validation_results,
  roc_train, roc_val_results,
  auc_summary,
  file = file.path(DATA_PROC, "risk_model_final.RData")
)

# 保存模型基因和权重为CSV（方便查阅）
model_info <- data.frame(
  gene = model_genes,
  weight = if (!is.null(model_weights)) model_weights[model_genes] else NA_real_,
  stringsAsFactors = FALSE
)
write.csv(model_info, file.path(RES_DIR, "uirs_model_genes.csv"), row.names = FALSE)

message("\n=== Risk score analysis completed ===")
message(sprintf("Model genes: %d", length(model_genes)))
message(sprintf("Training C-index: %.4f",
                calc_cindex(risk_score_train, train_expr$time, train_expr$status)))
message("Next step: Run 05_ml_model/04_independent_prognosis.R")
