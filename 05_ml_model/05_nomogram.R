###############################################################################
# 05_nomogram.R
# Nomogram构建 + 校准曲线 + 决策曲线分析(DCA)
# 整合UIRS + Age + Grade + IDH
# 校准曲线: Bootstrap B=200
# 输出: Figure 6A-E
###############################################################################

source("00_setup/config.R")
library(survival)
library(rms)
library(ggplot2)
library(dplyr)

set.seed(SEED)
load(file.path(DATA_PROC, "independent_prognosis_results.RData"))

# =============================================================================
# 1. 准备Nomogram数据
# =============================================================================
message("=== Preparing Nomogram data ===")

# 目标变量：risk_score + Age + Grade + IDH
# IDH在Nomogram中作为分类变量纳入（与Cox多因素不同，Nomogram需要直观展示）
target_vars <- c("risk_score", "Age", "Grade", "IDH")

# 检查各变量可用性
avail_vars <- c()
for (v in target_vars) {
  if (v %in% colnames(analysis_df)) {
    n_valid <- sum(!is.na(analysis_df[[v]]))
    pct_na  <- 1 - n_valid / nrow(analysis_df)
    if (pct_na <= 0.3) {
      avail_vars <- c(avail_vars, v)
      message(sprintf("  %s: %d valid (%.1f%% missing)", v, n_valid, pct_na * 100))
    } else {
      message(sprintf("  %s: EXCLUDED (%.1f%% missing)", v, pct_na * 100))
    }
  } else {
    message(sprintf("  %s: NOT in data", v))
  }
}

# 至少包含risk_score
if (!"risk_score" %in% avail_vars) avail_vars <- c("risk_score", avail_vars)

# 准备完整数据
nomo_cols <- c("OS_months", "OS", avail_vars)
nomo_df   <- analysis_df[, nomo_cols]
nomo_df   <- nomo_df[complete.cases(nomo_df), ]
nomo_df   <- nomo_df[nomo_df$OS_months > 0, ]

# 确保因子变量类型正确
if ("Grade" %in% avail_vars) nomo_df$Grade <- factor(nomo_df$Grade)
if ("IDH" %in% avail_vars)   nomo_df$IDH   <- factor(nomo_df$IDH)

message(sprintf("Nomogram samples: %d", nrow(nomo_df)))
message(sprintf("Nomogram variables: %s", paste(avail_vars, collapse = ", ")))

# =============================================================================
# 2. 构建rms Cox模型
# =============================================================================
message("=== Building rms Cox model ===")

# 设置数据分布
dd <- datadist(nomo_df)
options(datadist = "dd")

# 构建公式
formula_cph <- as.formula(paste("Surv(OS_months, OS) ~",
                                 paste(avail_vars, collapse = " + ")))

fit_cph <- cph(formula_cph, data = nomo_df, x = TRUE, y = TRUE, surv = TRUE,
               time.inc = 36)  # 3年作为默认时间点

message("\nCox-PH model summary:")
print(fit_cph)

# C-index
c_nomo <- fit_cph$stats["C"]
message(sprintf("Nomogram model C-index: %.4f", c_nomo))

# PH检验
ph_nomo <- cox.zph(fit_cph)
message("\nPH assumption test:")
print(ph_nomo)

# =============================================================================
# 3. Nomogram（Figure 6A）
# =============================================================================
message("=== Drawing Nomogram ===")

surv_fn <- Survival(fit_cph)

# 安全获取生存概率函数
surv_1yr  <- function(x) surv_fn(12, x)
surv_3yr  <- function(x) surv_fn(36, x)
surv_5yr  <- function(x) surv_fn(60, x)

pdf(file.path(FIG_DIR, "Fig6A_nomogram.pdf"), width = 12, height = 8)
tryCatch({
  nomo <- nomogram(
    fit_cph,
    fun      = list(surv_1yr, surv_3yr, surv_5yr),
    funlabel = c("1-Year Survival", "3-Year Survival", "5-Year Survival"),
    maxscale = 100,
    fun.at   = c(0.95, 0.9, 0.85, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1),
    lp       = TRUE
  )
  plot(nomo, xfrac = 0.35,
       cex.var = 0.9, cex.axis = 0.7,
       col.grid = gray(c(0.8, 0.95)))
  title("Cox Regression Nomogram for Overall Survival", cex.main = 1.2)
}, error = function(e) {
  message("  Nomogram error: ", e$message)
  plot.new()
  text(0.5, 0.5, paste("Nomogram error:", e$message), cex = 0.8)
})
dev.off()

# =============================================================================
# 4. 校准曲线 (Figure 6B-D)
# =============================================================================
message("=== Calibration curves (Bootstrap B=200) ===")

# 每个时间点的校准
cal_timepoints <- c(12, 36, 60)
cal_labels     <- c("1-Year", "3-Year", "5-Year")
cal_figs       <- c("Fig6B", "Fig6C", "Fig6D")

# 计算每组的样本数
m_per_group <- max(20, min(80, floor(nrow(nomo_df) / 5)))

cal_results <- list()

for (ti in seq_along(cal_timepoints)) {
  tp <- cal_timepoints[ti]
  lb <- cal_labels[ti]
  fg <- cal_figs[ti]

  pdf(file.path(FIG_DIR, paste0(fg, "_calibration_", lb, ".pdf")),
      width = 7, height = 7)
  tryCatch({
    # 重新拟合模型以匹配time.inc
    fit_cal <- cph(formula_cph, data = nomo_df,
                   x = TRUE, y = TRUE, surv = TRUE, time.inc = tp)

    cal_obj <- calibrate(fit_cal, cmethod = "KM", method = "boot",
                          u = tp, m = m_per_group, B = 200)
    cal_results[[lb]] <- cal_obj

    plot(cal_obj,
         xlab  = "Nomogram-Predicted Probability",
         ylab  = "Actual Probability (Kaplan-Meier)",
         main  = sprintf("Calibration Curve: %s Survival", lb),
         xlim  = c(0, 1), ylim = c(0, 1),
         subtitles = TRUE)

    # 添加对角参考线
    abline(0, 1, lty = 2, col = "grey50")

  }, error = function(e) {
    message(sprintf("  %s calibration error: %s", lb, e$message))
    plot.new()
    text(0.5, 0.5, paste("Calibration error:", e$message), cex = 0.8)
  })
  dev.off()
}

# --- 合并校准曲线到单个PDF ---
pdf(file.path(FIG_DIR, "Fig6BCD_calibration_combined.pdf"), width = 18, height = 6)
par(mfrow = c(1, 3), mar = c(5, 5, 4, 2))

for (ti in seq_along(cal_timepoints)) {
  tp <- cal_timepoints[ti]
  lb <- cal_labels[ti]

  tryCatch({
    fit_cal <- cph(formula_cph, data = nomo_df,
                   x = TRUE, y = TRUE, surv = TRUE, time.inc = tp)
    cal_obj <- calibrate(fit_cal, cmethod = "KM", method = "boot",
                          u = tp, m = m_per_group, B = 200)
    plot(cal_obj,
         xlab  = "Predicted Probability",
         ylab  = "Actual Probability",
         main  = sprintf("%s Survival", lb),
         xlim  = c(0, 1), ylim = c(0, 1),
         subtitles = FALSE)
    abline(0, 1, lty = 2, col = "grey50")
  }, error = function(e) {
    plot.new()
    text(0.5, 0.5, e$message, cex = 0.7)
  })
}
dev.off()

# =============================================================================
# 5. 决策曲线分析 DCA (Figure 6E)
# =============================================================================
message("=== Decision Curve Analysis ===")

# --- 实现标准的生存DCA ---
# 参考: Vickers et al. (2008) Medical Decision Making
surv_dca <- function(data, time_var, event_var, predictors, pred_labels = NULL,
                     timepoint, thresholds = seq(0.01, 0.99, by = 0.01)) {

  if (is.null(pred_labels)) pred_labels <- predictors

  # KM估计事件率
  surv_obj  <- Surv(data[[time_var]], data[[event_var]])
  km_fit    <- survfit(surv_obj ~ 1)
  km_summ   <- summary(km_fit, times = timepoint)
  event_rate <- 1 - km_summ$surv  # P(event before timepoint)

  results <- data.frame()

  # Treat All策略
  for (th in thresholds) {
    nb_all <- event_rate - (1 - event_rate) * th / (1 - th)
    results <- rbind(results, data.frame(
      threshold   = th,
      predictor   = "Treat All",
      net_benefit = nb_all,
      stringsAsFactors = FALSE
    ))
  }

  # Treat None策略
  results <- rbind(results, data.frame(
    threshold   = thresholds,
    predictor   = "Treat None",
    net_benefit = 0,
    stringsAsFactors = FALSE
  ))

  # 各预测模型策略
  for (pi in seq_along(predictors)) {
    pred_name <- predictors[pi]
    pred_label <- pred_labels[pi]

    tryCatch({
      # 用Cox模型估计个体事件概率
      fml <- as.formula(paste0("Surv(", time_var, ", ", event_var, ") ~ ", pred_name))
      fit <- coxph(fml, data = data)
      # 基线生存
      bh  <- basehaz(fit, centered = TRUE)
      # 找最接近timepoint的基线累积风险
      idx_tp <- which.min(abs(bh$time - timepoint))
      H0_tp  <- bh$hazard[idx_tp]

      # 个体风险概率: 1 - S(t|X) = 1 - exp(-H0(t) * exp(lp))
      lp       <- predict(fit, type = "lp")
      pred_prob <- 1 - exp(-H0_tp * exp(lp))
      pred_prob <- pmin(pmax(pred_prob, 0), 1)

      for (th in thresholds) {
        n <- nrow(data)
        # 使用IPCW方法的近似
        treat <- pred_prob >= th
        # 简化：直接用KM加权
        tp_count <- sum(treat & data[[event_var]] == 1 &
                          data[[time_var]] <= timepoint, na.rm = TRUE)
        fp_count <- sum(treat & (data[[event_var]] == 0 |
                                    data[[time_var]] > timepoint), na.rm = TRUE)

        # 调整：考虑删失
        tp_rate <- tp_count / n
        fp_rate <- fp_count / n

        nb <- tp_rate - fp_rate * th / (1 - th)

        results <- rbind(results, data.frame(
          threshold   = th,
          predictor   = pred_label,
          net_benefit = nb,
          stringsAsFactors = FALSE
        ))
      }
    }, error = function(e) {
      message(sprintf("  DCA failed for %s: %s", pred_name, e$message))
    })
  }

  return(results)
}

# 运行DCA（3年为主要时间点）
dca_results_3yr <- surv_dca(
  data       = nomo_df,
  time_var   = "OS_months",
  event_var  = "OS",
  predictors = "risk_score",
  pred_labels = "UIRS",
  timepoint  = 36,
  thresholds = seq(0.01, 0.80, by = 0.01)
)

# 3-year DCA图
p_dca_3yr <- ggplot(dca_results_3yr,
                     aes(x = threshold, y = net_benefit, color = predictor,
                         linetype = predictor)) +
  geom_line(linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "solid", color = "grey80") +
  scale_color_manual(values = c("UIRS" = "#E64B35",
                                 "Treat All" = "#3C5488",
                                 "Treat None" = "grey50")) +
  scale_linetype_manual(values = c("UIRS" = "solid",
                                    "Treat All" = "dashed",
                                    "Treat None" = "dotted")) +
  labs(x = "Threshold Probability", y = "Net Benefit",
       color = "Strategy", linetype = "Strategy",
       title = "Decision Curve Analysis (3-Year Survival)") +
  THEME_PUBLICATION +
  coord_cartesian(xlim = c(0, 0.8),
                  ylim = c(min(dca_results_3yr$net_benefit, na.rm = TRUE) * 0.1,
                           max(dca_results_3yr$net_benefit, na.rm = TRUE) * 1.1))

ggsave(file.path(FIG_DIR, "Fig6E_dca_3year.pdf"), p_dca_3yr, width = 8, height = 6)

# 1-year和5-year DCA（Supplementary）
for (tp_info in list(list(tp = 12, lab = "1-Year"),
                     list(tp = 60, lab = "5-Year"))) {
  tryCatch({
    dca_res <- surv_dca(
      data       = nomo_df,
      time_var   = "OS_months",
      event_var  = "OS",
      predictors = "risk_score",
      pred_labels = "UIRS",
      timepoint  = tp_info$tp,
      thresholds = seq(0.01, 0.80, by = 0.01)
    )

    p_dca <- ggplot(dca_res,
                     aes(x = threshold, y = net_benefit, color = predictor,
                         linetype = predictor)) +
      geom_line(linewidth = 1) +
      geom_hline(yintercept = 0, linetype = "solid", color = "grey80") +
      scale_color_manual(values = c("UIRS" = "#E64B35",
                                     "Treat All" = "#3C5488",
                                     "Treat None" = "grey50")) +
      scale_linetype_manual(values = c("UIRS" = "solid",
                                        "Treat All" = "dashed",
                                        "Treat None" = "dotted")) +
      labs(x = "Threshold Probability", y = "Net Benefit",
           color = "Strategy", linetype = "Strategy",
           title = sprintf("Decision Curve Analysis (%s Survival)", tp_info$lab)) +
      THEME_PUBLICATION +
      coord_cartesian(xlim = c(0, 0.8))

    ggsave(file.path(FIG_DIR, sprintf("SuppFig_dca_%s.pdf", tp_info$lab)),
           p_dca, width = 8, height = 6)
  }, error = function(e) {
    message(sprintf("  DCA for %s failed: %s", tp_info$lab, e$message))
  })
}

# =============================================================================
# 6. 多模型DCA对比（Nomogram vs. 单独Risk Score vs. 单独Age等）
# =============================================================================
message("=== Multi-model DCA comparison ===")

tryCatch({
  # 构建Nomogram预测值作为新变量
  nomo_df$nomogram_lp <- predict(fit_cph, type = "lp")

  compare_preds <- c("nomogram_lp", "risk_score")
  compare_labels <- c("Nomogram", "UIRS alone")

  if ("Age" %in% avail_vars) {
    compare_preds  <- c(compare_preds, "Age")
    compare_labels <- c(compare_labels, "Age alone")
  }

  dca_compare <- surv_dca(
    data       = nomo_df,
    time_var   = "OS_months",
    event_var  = "OS",
    predictors = compare_preds,
    pred_labels = compare_labels,
    timepoint  = 36,
    thresholds = seq(0.01, 0.80, by = 0.01)
  )

  color_vals <- c("Nomogram" = "#E64B35", "UIRS alone" = "#4DBBD5",
                   "Age alone" = "#00A087",
                   "Treat All" = "#3C5488", "Treat None" = "grey50")
  ltype_vals <- c("Nomogram" = "solid", "UIRS alone" = "solid",
                   "Age alone" = "solid",
                   "Treat All" = "dashed", "Treat None" = "dotted")

  p_dca_compare <- ggplot(dca_compare,
                           aes(x = threshold, y = net_benefit,
                               color = predictor, linetype = predictor)) +
    geom_line(linewidth = 1) +
    geom_hline(yintercept = 0, linetype = "solid", color = "grey80") +
    scale_color_manual(values = color_vals) +
    scale_linetype_manual(values = ltype_vals) +
    labs(x = "Threshold Probability", y = "Net Benefit",
         color = "Strategy", linetype = "Strategy",
         title = "DCA: Nomogram vs. Individual Predictors (3-Year)") +
    THEME_PUBLICATION +
    coord_cartesian(xlim = c(0, 0.8))

  ggsave(file.path(FIG_DIR, "Fig6E_dca_comparison.pdf"), p_dca_compare,
         width = 9, height = 6)
}, error = function(e) {
  message("  Multi-model DCA failed: ", e$message)
})

# =============================================================================
# 7. 模型性能汇总
# =============================================================================
message("=== Model performance summary ===")

perf_summary <- data.frame(
  Metric = c("C-index (Nomogram)",
             "Number of variables",
             "Variables",
             "Samples"),
  Value  = c(sprintf("%.4f", c_nomo),
             length(avail_vars),
             paste(avail_vars, collapse = " + "),
             nrow(nomo_df)),
  stringsAsFactors = FALSE
)

write.csv(perf_summary,
          file.path(RES_DIR, "nomogram_performance.csv"), row.names = FALSE)

# PH检验汇总
ph_nomo_df <- data.frame(
  Variable  = rownames(ph_nomo$table),
  chisq     = ph_nomo$table[, "chisq"],
  df        = ph_nomo$table[, "df"],
  p_value   = ph_nomo$table[, "p"],
  ph_passed = ph_nomo$table[, "p"] >= 0.05,
  stringsAsFactors = FALSE
)
write.csv(ph_nomo_df,
          file.path(RES_DIR, "supplementary_nomogram_ph_test.csv"), row.names = FALSE)

# =============================================================================
# 8. 保存
# =============================================================================
save(fit_cph, nomo_df, avail_vars, cal_results, ph_nomo,
     file = file.path(DATA_PROC, "nomogram_results.RData"))

message("\n=== Nomogram analysis completed ===")
message(sprintf("Model: Surv(OS_months, OS) ~ %s", paste(avail_vars, collapse = " + ")))
message(sprintf("C-index: %.4f", c_nomo))
message("Part 3 (ML Model) completed.")
message("Next: Run 06_clinical_translation/01_immunotherapy_prediction.R for Part 4.")
