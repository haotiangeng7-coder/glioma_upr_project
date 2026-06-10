###############################################################################
# 04_independent_prognosis.R
# 单因素+多因素Cox回归 → 确认Risk Score独立预后价值
# 评审修订：
#   - IDH默认strata(IDH)，不作为协变量
#   - 所有Cox模型运行cox.zph()，报告PH检验于Supplementary Table
#   - 违反PH的变量使用时变协变量或RMST
#   - IDH-WT亚组独立重复
#   - Forest plot
# 输出: Figure 6F-G
###############################################################################

source("00_setup/config.R")
library(survival)
library(survminer)
library(survRM2)
library(ggplot2)
library(dplyr)

set.seed(SEED)
load(file.path(DATA_PROC, "risk_model_final.RData"))
load(file.path(DATA_PROC, "tcga_glioma_expression.RData"))

# =============================================================================
# 1. 准备分析数据
# =============================================================================
message("=== Preparing multivariate analysis data ===")

analysis_df <- data.frame(
  barcode    = names(risk_score_train),
  risk_score = risk_score_train,
  risk_group = risk_group_train,
  stringsAsFactors = FALSE
)

analysis_df <- merge(analysis_df, clinical_valid, by = "barcode")

analysis_df <- analysis_df %>%
  dplyr::mutate(
    Age       = as.numeric(age_at_index),
    Gender    = factor(gender),
    Grade     = factor(Grade),
    IDH       = factor(IDH_status),
    MGMT      = factor(MGMT_status),
    OS_months = OS.time / 30
  ) %>%
  dplyr::filter(!is.na(OS_months) & OS_months > 0 & !is.na(OS))

message(sprintf("Analysis samples: %d", nrow(analysis_df)))
message(sprintf("IDH levels: %s", paste(levels(analysis_df$IDH), collapse = ", ")))

# =============================================================================
# 辅助函数
# =============================================================================
# 安全的单因素Cox回归，带PH检验
safe_univar_cox <- function(data, var_name, var_label, use_strata_idh = TRUE) {
  tryCatch({
    # 检查变量可用性
    if (!var_name %in% colnames(data)) return(NULL)
    n_valid <- sum(!is.na(data[[var_name]]))
    if (n_valid < 50) return(NULL)
    if (is.factor(data[[var_name]]) && length(unique(na.omit(data[[var_name]]))) < 2) return(NULL)

    # 构建公式 - IDH作为分层变量
    if (use_strata_idh && "IDH" %in% colnames(data) &&
        sum(!is.na(data$IDH)) > 50 && var_name != "IDH") {
      fml <- as.formula(paste0("Surv(OS_months, OS) ~ ", var_name, " + strata(IDH)"))
    } else {
      fml <- as.formula(paste0("Surv(OS_months, OS) ~ ", var_name))
    }

    fit <- coxph(fml, data = data)
    s   <- summary(fit)

    # PH检验
    ph <- cox.zph(fit)
    # 取目标变量的PH p值（非strata行）
    ph_row <- grep(var_name, rownames(ph$table))
    ph_p   <- if (length(ph_row) > 0) ph$table[ph_row[1], "p"] else ph$table["GLOBAL", "p"]

    # 如果多水平factor，取最高级别的HR
    n_coefs <- nrow(s$coefficients)
    # 只取包含var_name的行
    var_rows <- grep(var_name, rownames(s$coefficients))
    if (length(var_rows) == 0) return(NULL)
    idx <- var_rows[length(var_rows)]  # 最高级别

    data.frame(
      Variable    = var_label,
      HR          = s$conf.int[idx, 1],
      HR_lower    = s$conf.int[idx, 3],
      HR_upper    = s$conf.int[idx, 4],
      pvalue      = ifelse(length(var_rows) > 1, s$sctest["pvalue"], s$coefficients[idx, 5]),
      ph_pvalue   = ph_p,
      ph_passed   = (ph_p >= 0.05),
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    message(sprintf("  Univariate Cox failed for %s: %s", var_name, e$message))
    NULL
  })
}

# =============================================================================
# 2. 单因素Cox回归（IDH分层）
# =============================================================================
message("=== Univariate Cox regression (IDH stratified) ===")

univar_list <- list(
  safe_univar_cox(analysis_df, "risk_score", "UIRS Risk Score"),
  safe_univar_cox(analysis_df, "Age",        "Age"),
  safe_univar_cox(analysis_df, "Gender",     "Gender (Male vs Female)"),
  safe_univar_cox(analysis_df, "Grade",      "WHO Grade"),
  safe_univar_cox(analysis_df, "MGMT",       "MGMT Methylation")
)

# IDH单独（不加strata）
idh_res <- tryCatch({
  if (sum(!is.na(analysis_df$IDH)) > 50 &&
      length(unique(na.omit(analysis_df$IDH))) >= 2) {
    fit <- coxph(Surv(OS_months, OS) ~ IDH, data = analysis_df)
    s   <- summary(fit)
    ph  <- cox.zph(fit)
    data.frame(
      Variable  = "IDH (WT vs Mutant)",
      HR        = s$conf.int[1, 1],
      HR_lower  = s$conf.int[1, 3],
      HR_upper  = s$conf.int[1, 4],
      pvalue    = s$coefficients[1, 5],
      ph_pvalue = ph$table[1, "p"],
      ph_passed = (ph$table[1, "p"] >= 0.05),
      stringsAsFactors = FALSE
    )
  } else NULL
}, error = function(e) NULL)

univar_list <- c(univar_list, list(idh_res))

univar_results <- do.call(rbind, univar_list[!sapply(univar_list, is.null)])
rownames(univar_results) <- NULL

message("Univariate Cox results (IDH stratified):")
print(univar_results)

# =============================================================================
# 3. 多因素Cox回归（IDH作为strata分层变量）
# =============================================================================
message("\n=== Multivariate Cox regression (IDH as strata) ===")

# 构建协变量列表
formula_parts <- "risk_score"
if (sum(!is.na(analysis_df$Age)) > 50) formula_parts <- c(formula_parts, "Age")
if (sum(!is.na(analysis_df$Gender)) > 50 &&
    length(unique(na.omit(analysis_df$Gender))) >= 2) {
  formula_parts <- c(formula_parts, "Gender")
}
if (sum(!is.na(analysis_df$Grade)) > 50 &&
    length(unique(na.omit(analysis_df$Grade))) >= 2) {
  formula_parts <- c(formula_parts, "Grade")
}
if (sum(!is.na(analysis_df$MGMT)) > 50 &&
    length(unique(na.omit(analysis_df$MGMT))) >= 2) {
  formula_parts <- c(formula_parts, "MGMT")
}

# IDH作为strata（不作为协变量）
has_idh_strata <- (sum(!is.na(analysis_df$IDH)) > 50 &&
                   length(unique(na.omit(analysis_df$IDH))) >= 2)

if (has_idh_strata) {
  formula_multi <- as.formula(paste(
    "Surv(OS_months, OS) ~",
    paste(formula_parts, collapse = " + "),
    "+ strata(IDH)"
  ))
} else {
  formula_multi <- as.formula(paste(
    "Surv(OS_months, OS) ~",
    paste(formula_parts, collapse = " + ")
  ))
}

# 完整数据子集
multi_vars <- c("OS_months", "OS", formula_parts)
if (has_idh_strata) multi_vars <- c(multi_vars, "IDH")
multi_df <- analysis_df[complete.cases(analysis_df[, multi_vars]), ]
message(sprintf("Multivariate model samples: %d", nrow(multi_df)))

fit_multi <- coxph(formula_multi, data = multi_df)
s_multi   <- summary(fit_multi)

message("\nMultivariate Cox model:")
print(s_multi)

# PH检验
ph_multi <- cox.zph(fit_multi)
message("\nPH assumption test (multivariate):")
print(ph_multi)

# 提取多因素结果
multivar_results <- data.frame(
  Variable  = rownames(s_multi$coefficients),
  HR        = s_multi$conf.int[, 1],
  HR_lower  = s_multi$conf.int[, 3],
  HR_upper  = s_multi$conf.int[, 4],
  pvalue    = s_multi$coefficients[, 5],
  stringsAsFactors = FALSE
)

# 添加多因素PH检验结果
ph_multi_df <- data.frame(
  Variable  = rownames(ph_multi$table),
  ph_pvalue = ph_multi$table[, "p"],
  ph_passed = ph_multi$table[, "p"] >= 0.05,
  stringsAsFactors = FALSE
)

multivar_results <- multivar_results %>%
  dplyr::left_join(ph_multi_df, by = "Variable")

rownames(multivar_results) <- NULL

# =============================================================================
# 4. 处理PH违反的变量
# =============================================================================
message("=== Handling PH violations ===")

ph_violations <- multivar_results %>%
  dplyr::filter(!is.na(ph_passed) & !ph_passed)

if (nrow(ph_violations) > 0) {
  message("Variables violating PH assumption:")
  print(ph_violations[, c("Variable", "ph_pvalue")])

  # 对违反PH的变量尝试时变协变量模型
  for (vi in seq_len(nrow(ph_violations))) {
    var_name <- ph_violations$Variable[vi]
    message(sprintf("\n  Handling PH violation for: %s", var_name))

    # 清理变量名（去除factor level后缀）
    base_var <- gsub("(Mutant|Unmethylated|male|G4|G3|G2)$", "", var_name)
    base_var <- trimws(base_var)

    tryCatch({
      # 尝试时变协变量：tt()函数
      if (base_var %in% colnames(multi_df) && base_var != "risk_score") {
        fml_tv <- as.formula(paste0(
          "Surv(OS_months, OS) ~ ", base_var, " + tt(", base_var, ")",
          if (has_idh_strata) " + strata(IDH)" else ""
        ))
        fit_tv <- coxph(fml_tv, data = multi_df,
                        tt = function(x, t, ...) x * log(t + 1))
        message(sprintf("  Time-varying model for %s:", base_var))
        print(summary(fit_tv)$coefficients)
      }
    }, error = function(e) {
      message(sprintf("  Time-varying model failed for %s: %s", var_name, e$message))
    })
  }
} else {
  message("No PH violations detected in multivariate model.")
}

# =============================================================================
# 5. IDH-WT亚组独立分析
# =============================================================================
message("\n=== IDH-WT subgroup analysis ===")

wt_df <- analysis_df %>% dplyr::filter(IDH == "WT")
message(sprintf("IDH-WT subgroup: %d samples", nrow(wt_df)))

if (nrow(wt_df) >= 50) {
  # 单因素
  wt_univar_list <- list(
    safe_univar_cox(wt_df, "risk_score", "UIRS Risk Score", use_strata_idh = FALSE),
    safe_univar_cox(wt_df, "Age",        "Age",             use_strata_idh = FALSE),
    safe_univar_cox(wt_df, "Gender",     "Gender",          use_strata_idh = FALSE),
    safe_univar_cox(wt_df, "Grade",      "WHO Grade",       use_strata_idh = FALSE),
    safe_univar_cox(wt_df, "MGMT",       "MGMT",            use_strata_idh = FALSE)
  )
  wt_univar <- do.call(rbind, wt_univar_list[!sapply(wt_univar_list, is.null)])

  if (!is.null(wt_univar) && nrow(wt_univar) > 0) {
    message("IDH-WT univariate Cox:")
    print(wt_univar)
    write.csv(wt_univar, file.path(RES_DIR, "idh_wt_univariate_cox.csv"), row.names = FALSE)
  }

  # 多因素（不含IDH）
  wt_formula_parts <- "risk_score"
  if (sum(!is.na(wt_df$Age)) > 20) wt_formula_parts <- c(wt_formula_parts, "Age")
  if (sum(!is.na(wt_df$Grade)) > 20 &&
      length(unique(na.omit(wt_df$Grade))) >= 2) {
    wt_formula_parts <- c(wt_formula_parts, "Grade")
  }
  if (sum(!is.na(wt_df$MGMT)) > 20 &&
      length(unique(na.omit(wt_df$MGMT))) >= 2) {
    wt_formula_parts <- c(wt_formula_parts, "MGMT")
  }

  wt_multi_vars <- c("OS_months", "OS", wt_formula_parts)
  wt_multi_df   <- wt_df[complete.cases(wt_df[, wt_multi_vars]), ]

  if (nrow(wt_multi_df) >= 30) {
    tryCatch({
      wt_fml <- as.formula(paste("Surv(OS_months, OS) ~",
                                  paste(wt_formula_parts, collapse = " + ")))
      wt_fit <- coxph(wt_fml, data = wt_multi_df)
      wt_s   <- summary(wt_fit)
      wt_ph  <- cox.zph(wt_fit)

      message("\nIDH-WT multivariate Cox:")
      print(wt_s)
      message("\nIDH-WT PH test:")
      print(wt_ph)

      wt_multivar <- data.frame(
        Variable  = rownames(wt_s$coefficients),
        HR        = wt_s$conf.int[, 1],
        HR_lower  = wt_s$conf.int[, 3],
        HR_upper  = wt_s$conf.int[, 4],
        pvalue    = wt_s$coefficients[, 5],
        stringsAsFactors = FALSE
      )

      write.csv(wt_multivar,
                file.path(RES_DIR, "idh_wt_multivariate_cox.csv"), row.names = FALSE)
    }, error = function(e) {
      message("  IDH-WT multivariate failed: ", e$message)
    })
  }
} else {
  message("Insufficient IDH-WT samples for subgroup analysis.")
}

# =============================================================================
# 6. Forest plot
# =============================================================================
message("=== Generating forest plots ===")

# 格式化辅助
format_forest_data <- function(df, analysis_type) {
  df$ci_text <- sprintf("%.2f (%.2f-%.2f)", df$HR, df$HR_lower, df$HR_upper)
  df$p_text  <- ifelse(df$pvalue < 0.001,
                        sprintf("%.2e", df$pvalue),
                        sprintf("%.3f", df$pvalue))
  df$sig     <- ifelse(df$pvalue < 0.05, "*", "")
  df$analysis <- analysis_type
  df
}

# --- Fig 6F: 单因素 Forest plot ---
univar_plot <- format_forest_data(univar_results, "Univariate")

p_forest_uni <- ggplot(univar_plot, aes(x = HR, y = reorder(Variable, -pvalue))) +
  geom_point(size = 3, color = "#3C5488") +
  geom_errorbarh(aes(xmin = HR_lower, xmax = HR_upper),
                 height = 0.2, color = "#3C5488") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_text(aes(label = ci_text),
            x = max(univar_plot$HR_upper, na.rm = TRUE) * 1.5,
            hjust = 0, size = 3, color = "grey30") +
  geom_text(aes(label = paste0("p=", p_text)),
            x = max(univar_plot$HR_upper, na.rm = TRUE) * 3,
            hjust = 0, size = 3, color = "grey30") +
  labs(x = "Hazard Ratio (95% CI)", y = "",
       title = "Univariate Cox Regression (IDH stratified)",
       subtitle = "strata(IDH) applied") +
  THEME_PUBLICATION +
  scale_x_log10()

ggsave(file.path(FIG_DIR, "Fig6F_forest_univariate.pdf"), p_forest_uni,
       width = 10, height = max(4, nrow(univar_plot) * 0.7))

# --- Fig 6G: 多因素 Forest plot ---
multivar_plot <- format_forest_data(multivar_results, "Multivariate")

p_forest_multi <- ggplot(multivar_plot, aes(x = HR, y = reorder(Variable, -pvalue))) +
  geom_point(size = 3, color = "#E64B35") +
  geom_errorbarh(aes(xmin = HR_lower, xmax = HR_upper),
                 height = 0.2, color = "#E64B35") +
  geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
  geom_text(aes(label = ci_text),
            x = max(multivar_plot$HR_upper, na.rm = TRUE) * 1.5,
            hjust = 0, size = 3, color = "grey30") +
  geom_text(aes(label = paste0("p=", p_text)),
            x = max(multivar_plot$HR_upper, na.rm = TRUE) * 3,
            hjust = 0, size = 3, color = "grey30") +
  labs(x = "Hazard Ratio (95% CI)", y = "",
       title = "Multivariate Cox Regression (IDH stratified)",
       subtitle = "IDH as strata, not covariate") +
  THEME_PUBLICATION +
  scale_x_log10()

ggsave(file.path(FIG_DIR, "Fig6G_forest_multivariate.pdf"), p_forest_multi,
       width = 10, height = max(4, nrow(multivar_plot) * 0.7))

# --- 合并Forest plot ---
p_forest_combined <- cowplot::plot_grid(
  p_forest_uni, p_forest_multi,
  ncol = 1, labels = c("F", "G"), label_size = 16
)
ggsave(file.path(FIG_DIR, "Fig6FG_forest_combined.pdf"), p_forest_combined,
       width = 10, height = max(8, (nrow(univar_plot) + nrow(multivar_plot)) * 0.6))

# =============================================================================
# 7. PH检验汇总表（Supplementary Table）
# =============================================================================
message("=== Saving PH assumption tests (Supplementary) ===")

# 单因素PH汇总
ph_supp_univar <- univar_results %>%
  dplyr::select(Variable, HR, pvalue, ph_pvalue, ph_passed) %>%
  dplyr::mutate(analysis = "Univariate")

# 多因素PH汇总
ph_supp_multivar <- multivar_results %>%
  dplyr::select(Variable, HR, pvalue, ph_pvalue, ph_passed) %>%
  dplyr::mutate(analysis = "Multivariate")

ph_supp_all <- rbind(ph_supp_univar, ph_supp_multivar)
write.csv(ph_supp_all,
          file.path(RES_DIR, "supplementary_cox_ph_tests.csv"), row.names = FALSE)

# =============================================================================
# 8. 保存结果
# =============================================================================
save(univar_results, multivar_results, fit_multi, ph_multi,
     analysis_df, multi_df,
     file = file.path(DATA_PROC, "independent_prognosis_results.RData"))

write.csv(univar_results,
          file.path(RES_DIR, "univariate_cox_clinical.csv"), row.names = FALSE)
write.csv(multivar_results,
          file.path(RES_DIR, "multivariate_cox_clinical.csv"), row.names = FALSE)

message("\n=== Independent prognosis analysis completed ===")
message(sprintf("Multivariate model formula: %s",
                deparse(formula(fit_multi), width.cutoff = 200)))
message(sprintf("IDH handled as: strata (not covariate)"))
message("Next step: Run 05_ml_model/05_nomogram.R")
