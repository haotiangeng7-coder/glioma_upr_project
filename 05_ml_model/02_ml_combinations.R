###############################################################################
# 02_ml_combinations.R
# 10种算法两两组合（~100+种）+ 5x5嵌套CV + 三级验证 + 置换检验
# 评审修订版：嵌套CV无偏性能估计 + 1000次置换检验 + 小提琴图
# 参考：Liu et al. 2024 Genome Medicine
###############################################################################

source("00_setup/config.R")
library(survival)
library(glmnet)
library(randomForestSRC)
library(CoxBoost)
library(superpc)
library(plsRcox)
library(ggplot2)
library(dplyr)
library(tidyr)

set.seed(SEED)
load(file.path(DATA_PROC, "feature_selection_results.RData"))
load(file.path(DATA_PROC, "tcga_glioma_expression.RData"))

# =============================================================================
# 1. 准备训练集和验证集数据
# =============================================================================
message("=== Preparing datasets ===")

genes <- candidate_genes
message(sprintf("Using %d candidate genes", length(genes)))

# --- 训练集 (TCGA) ---
train_expr <- surv_df[, c("time", "status", genes)]
train_expr <- train_expr[complete.cases(train_expr), ]
message(sprintf("Training set (TCGA): %d samples", nrow(train_expr)))

# --- 辅助函数：准备外部验证集 ---
prepare_external_dataset <- function(expr_mat, clinical_df, id_col, genes) {
  common_genes <- intersect(genes, rownames(expr_mat))
  if (length(common_genes) < length(genes) * 0.7) return(NULL)
  # 对缺失基因填0
  val_expr <- matrix(0, nrow = ncol(expr_mat), ncol = length(genes),
                     dimnames = list(colnames(expr_mat), genes))
  val_expr[, common_genes] <- log2(t(expr_mat[common_genes, ]) + 1)
  common_ids <- intersect(rownames(val_expr), clinical_df[[id_col]])
  if (length(common_ids) <= 20) return(NULL)
  # 自动检测生存状态列（OS_status优先于OS，因OS在CGGA中是生存时间）
  status_col <- if ("OS_status" %in% colnames(clinical_df)) "OS_status"
                else if ("OS" %in% colnames(clinical_df)) "OS"
                else stop("No survival status column found")
  val_df <- data.frame(
    time   = clinical_df$OS.time[match(common_ids, clinical_df[[id_col]])],
    status = clinical_df[[status_col]][match(common_ids, clinical_df[[id_col]])],
    val_expr[common_ids, , drop = FALSE],
    check.names = FALSE
  )
  val_df <- val_df[complete.cases(val_df) & val_df$time > 0, ]
  if (nrow(val_df) <= 20) return(NULL)
  return(val_df)
}

# --- 三级验证集 ---
# 调优验证集: CGGA_batch2
# 盲测试集: CGGA_batch1, GSE16011
# 探索性: GSE43378

tuning_data    <- list()  # 参与模型排名
blind_test     <- list()  # 最终才评估一次
exploratory    <- list()  # 补充

# CGGA
if (file.exists(file.path(DATA_PROC, "cgga_data.RData"))) {
  load(file.path(DATA_PROC, "cgga_data.RData"))
  if (!is.null(cgga_b2)) {
    vd <- prepare_external_dataset(cgga_b2$expr, cgga_b2$clinical, "Sample_ID", genes)
    if (!is.null(vd)) {
      tuning_data[["CGGA_batch2"]] <- vd
      message(sprintf("  CGGA_batch2 (tuning): %d samples", nrow(vd)))
    }
  }
  if (!is.null(cgga_b1)) {
    vd <- prepare_external_dataset(cgga_b1$expr, cgga_b1$clinical, "Sample_ID", genes)
    if (!is.null(vd)) {
      blind_test[["CGGA_batch1"]] <- vd
      message(sprintf("  CGGA_batch1 (blind test): %d samples", nrow(vd)))
    }
  }
}

# GSE16011
if (file.exists(file.path(DATA_PROC, "gse16011_data.RData"))) {
  load(file.path(DATA_PROC, "gse16011_data.RData"))
  if (exists("expr_16011_gene") && exists("clin_16011_clean")) {
    # 尝试提取生存信息
    if ("OS.time" %in% colnames(clin_16011_clean) && "OS" %in% colnames(clin_16011_clean)) {
      vd <- prepare_external_dataset(expr_16011_gene, clin_16011_clean, "Sample_ID", genes)
      if (!is.null(vd)) {
        blind_test[["GSE16011"]] <- vd
        message(sprintf("  GSE16011 (blind test): %d samples", nrow(vd)))
      }
    }
  }
}

# GSE43378
if (file.exists(file.path(DATA_PROC, "gse43378_data.RData"))) {
  load(file.path(DATA_PROC, "gse43378_data.RData"))
  if (exists("expr_43378_gene") && exists("clin_43378")) {
    if ("OS.time" %in% colnames(clin_43378) && "OS" %in% colnames(clin_43378)) {
      vd <- prepare_external_dataset(expr_43378_gene, clin_43378, "Sample_ID", genes)
      if (!is.null(vd)) {
        exploratory[["GSE43378"]] <- vd
        message(sprintf("  GSE43378 (exploratory): %d samples", nrow(vd)))
      }
    }
  }
}

# 合并所有验证集用于后续保存
validation_data <- c(tuning_data, blind_test, exploratory)
message(sprintf("Total validation datasets: tuning=%d, blind=%d, exploratory=%d",
                length(tuning_data), length(blind_test), length(exploratory)))

# =============================================================================
# 2. C-index 计算辅助函数
# =============================================================================
calc_cindex <- function(pred_risk, time, status) {
  tryCatch({
    # concordance(Surv ~ pred) returns discordance for risk scores
    # For risk scores (higher = worse), use 1 - concordance or coxph wrapper
    fit <- coxph(Surv(time, status) ~ pred_risk)
    return(concordance(fit)$concordance)
  }, error = function(e) return(NA_real_))
}

# =============================================================================
# 3. 定义10种基础算法（特征选择器 + 模型构建器）
# =============================================================================
message("=== Defining ML algorithms ===")

# 每种算法返回一个predict函数和训练时C-index
# 统一接口: algo(x_train, y_surv, ...) -> list(predict_fn, model)

# --- 3.1 特征选择算法 ---
# 返回选中的基因子集

fs_lasso <- function(x, y, seed = SEED) {
  set.seed(seed)
  cv <- cv.glmnet(x, y, family = "cox", alpha = 1, nfolds = ML_NESTED_INNER)
  coefs <- coef(cv, s = "lambda.min")
  selected <- rownames(coefs)[coefs[, 1] != 0]
  if (length(selected) == 0) selected <- colnames(x)
  list(selected = selected, model = cv, type = "glmnet")
}

fs_ridge <- function(x, y, seed = SEED) {
  set.seed(seed)
  cv <- cv.glmnet(x, y, family = "cox", alpha = 0, nfolds = ML_NESTED_INNER)
  list(selected = colnames(x), model = cv, type = "glmnet")
}

fs_enet03 <- function(x, y, seed = SEED) {
  set.seed(seed)
  cv <- cv.glmnet(x, y, family = "cox", alpha = 0.3, nfolds = ML_NESTED_INNER)
  coefs <- coef(cv, s = "lambda.min")
  selected <- rownames(coefs)[coefs[, 1] != 0]
  if (length(selected) == 0) selected <- colnames(x)
  list(selected = selected, model = cv, type = "glmnet")
}

fs_enet05 <- function(x, y, seed = SEED) {
  set.seed(seed)
  cv <- cv.glmnet(x, y, family = "cox", alpha = 0.5, nfolds = ML_NESTED_INNER)
  coefs <- coef(cv, s = "lambda.min")
  selected <- rownames(coefs)[coefs[, 1] != 0]
  if (length(selected) == 0) selected <- colnames(x)
  list(selected = selected, model = cv, type = "glmnet")
}

fs_enet07 <- function(x, y, seed = SEED) {
  set.seed(seed)
  cv <- cv.glmnet(x, y, family = "cox", alpha = 0.7, nfolds = ML_NESTED_INNER)
  coefs <- coef(cv, s = "lambda.min")
  selected <- rownames(coefs)[coefs[, 1] != 0]
  if (length(selected) == 0) selected <- colnames(x)
  list(selected = selected, model = cv, type = "glmnet")
}

fs_rsf <- function(x, y, seed = SEED) {
  set.seed(seed)
  df <- data.frame(time = y[, "time"], status = y[, "status"], x)
  fit <- rfsrc(Surv(time, status) ~ ., data = df, ntree = 500, nodesize = 10,
               seed = seed, importance = TRUE)
  vimp <- fit$importance
  selected <- names(sort(vimp, decreasing = TRUE))[vimp[order(vimp, decreasing = TRUE)] > 0]
  if (length(selected) < 3) selected <- names(sort(vimp, decreasing = TRUE))[1:min(ncol(x), 10)]
  list(selected = selected, model = fit, type = "rsf")
}

fs_coxboost <- function(x, y, seed = SEED) {
  set.seed(seed)
  tryCatch({
    cv_cb <- cv.CoxBoost(
      time = y[, "time"], status = y[, "status"], x = x,
      maxstepno = 200, K = ML_NESTED_INNER, type = "verweij"
    )
    fit <- CoxBoost(
      time = y[, "time"], status = y[, "status"], x = x,
      stepno = cv_cb$optimal.step
    )
    coefs <- coef(fit)
    selected <- colnames(x)[coefs != 0]
    if (length(selected) == 0) selected <- colnames(x)
    list(selected = selected, model = fit, type = "coxboost")
  }, error = function(e) {
    list(selected = colnames(x), model = NULL, type = "coxboost_fail")
  })
}

fs_superpc <- function(x, y, seed = SEED) {
  set.seed(seed)
  tryCatch({
    spc_data <- list(x = t(x), y = y[, "time"],
                     censoring.status = y[, "status"],
                     featurenames = colnames(x))
    fit <- superpc.train(spc_data, type = "survival")
    cv <- superpc.cv(fit, spc_data, n.fold = ML_NESTED_INNER)
    best_th <- cv$thresholds[which.max(cv$scor)]
    # SuperPC选择超过阈值的特征
    scores <- abs(fit$feature.scores)
    selected <- colnames(x)[scores > best_th]
    if (length(selected) < 3) selected <- colnames(x)[order(scores, decreasing = TRUE)[1:min(ncol(x), 10)]]
    list(selected = selected, model = list(fit = fit, data = spc_data, threshold = best_th),
         type = "superpc")
  }, error = function(e) {
    list(selected = colnames(x), model = NULL, type = "superpc_fail")
  })
}

fs_stepcox <- function(x, y, seed = SEED) {
  set.seed(seed)
  tryCatch({
    df <- data.frame(time = y[, "time"], status = y[, "status"], x)
    fml <- as.formula(paste0("Surv(time, status) ~ ",
                              paste0("`", colnames(x), "`", collapse = " + ")))
    fit_full <- coxph(fml, data = df)
    fit_step <- step(fit_full, direction = "both", trace = 0)
    selected <- names(coef(fit_step))
    selected <- gsub("`", "", selected)
    if (length(selected) == 0) selected <- colnames(x)
    list(selected = selected, model = fit_step, type = "stepcox")
  }, error = function(e) {
    list(selected = colnames(x), model = NULL, type = "stepcox_fail")
  })
}

fs_plsrcox <- function(x, y, seed = SEED) {
  set.seed(seed)
  tryCatch({
    fit <- plsRcox(
      Xplan = x, time = y[, "time"], event = y[, "status"],
      nt = min(5, ncol(x) - 1)
    )
    list(selected = colnames(x), model = fit, type = "plsrcox")
  }, error = function(e) {
    list(selected = colnames(x), model = NULL, type = "plsrcox_fail")
  })
}

# --- 3.2 模型构建器 ---
# 输入selected基因子集, 返回predict_fn(newx) -> risk_scores

build_lasso <- function(x, y, seed = SEED) {
  set.seed(seed)
  cv <- cv.glmnet(x, y, family = "cox", alpha = 1, nfolds = ML_NESTED_INNER)
  list(
    predict_fn = function(newx) predict(cv, newx = as.matrix(newx[, colnames(x), drop = FALSE]),
                                         s = "lambda.min", type = "link")[, 1],
    model = cv
  )
}

build_ridge <- function(x, y, seed = SEED) {
  set.seed(seed)
  cv <- cv.glmnet(x, y, family = "cox", alpha = 0, nfolds = ML_NESTED_INNER)
  list(
    predict_fn = function(newx) predict(cv, newx = as.matrix(newx[, colnames(x), drop = FALSE]),
                                         s = "lambda.min", type = "link")[, 1],
    model = cv
  )
}

build_enet03 <- function(x, y, seed = SEED) {
  set.seed(seed)
  cv <- cv.glmnet(x, y, family = "cox", alpha = 0.3, nfolds = ML_NESTED_INNER)
  list(
    predict_fn = function(newx) predict(cv, newx = as.matrix(newx[, colnames(x), drop = FALSE]),
                                         s = "lambda.min", type = "link")[, 1],
    model = cv
  )
}

build_enet05 <- function(x, y, seed = SEED) {
  set.seed(seed)
  cv <- cv.glmnet(x, y, family = "cox", alpha = 0.5, nfolds = ML_NESTED_INNER)
  list(
    predict_fn = function(newx) predict(cv, newx = as.matrix(newx[, colnames(x), drop = FALSE]),
                                         s = "lambda.min", type = "link")[, 1],
    model = cv
  )
}

build_enet07 <- function(x, y, seed = SEED) {
  set.seed(seed)
  cv <- cv.glmnet(x, y, family = "cox", alpha = 0.7, nfolds = ML_NESTED_INNER)
  list(
    predict_fn = function(newx) predict(cv, newx = as.matrix(newx[, colnames(x), drop = FALSE]),
                                         s = "lambda.min", type = "link")[, 1],
    model = cv
  )
}

build_rsf <- function(x, y, seed = SEED) {
  set.seed(seed)
  df <- data.frame(time = y[, "time"], status = y[, "status"], x)
  gn <- colnames(x)
  fit <- rfsrc(Surv(time, status) ~ ., data = df, ntree = 1000, nodesize = 10, seed = seed)
  list(
    predict_fn = function(newx) {
      nd <- data.frame(time = 0, status = 0, newx[, gn, drop = FALSE])
      predict(fit, newdata = nd)$predicted
    },
    model = fit
  )
}

build_coxboost <- function(x, y, seed = SEED) {
  set.seed(seed)
  tryCatch({
    cv_cb <- cv.CoxBoost(
      time = y[, "time"], status = y[, "status"], x = x,
      maxstepno = 200, K = ML_NESTED_INNER, type = "verweij"
    )
    fit <- CoxBoost(
      time = y[, "time"], status = y[, "status"], x = x,
      stepno = cv_cb$optimal.step
    )
    gn <- colnames(x)
    list(
      predict_fn = function(newx) predict(fit, newdata = as.matrix(newx[, gn, drop = FALSE]),
                                           type = "lp"),
      model = fit
    )
  }, error = function(e) {
    list(predict_fn = function(newx) rep(NA, nrow(newx)), model = NULL)
  })
}

build_stepcox <- function(x, y, seed = SEED) {
  set.seed(seed)
  tryCatch({
    df <- data.frame(time = y[, "time"], status = y[, "status"], x, check.names = FALSE)
    fml <- as.formula(paste0("Surv(time, status) ~ ",
                              paste0("`", colnames(x), "`", collapse = " + ")))
    fit_full <- coxph(fml, data = df)
    fit_step <- step(fit_full, direction = "both", trace = 0)
    list(
      predict_fn = function(newx) {
        nd <- data.frame(newx, check.names = FALSE)
        predict(fit_step, newdata = nd, type = "lp")
      },
      model = fit_step
    )
  }, error = function(e) {
    list(predict_fn = function(newx) rep(NA, nrow(newx)), model = NULL)
  })
}

build_superpc <- function(x, y, seed = SEED) {
  set.seed(seed)
  tryCatch({
    spc_data <- list(x = t(x), y = y[, "time"],
                     censoring.status = y[, "status"],
                     featurenames = colnames(x))
    fit <- superpc.train(spc_data, type = "survival")
    cv <- superpc.cv(fit, spc_data, n.fold = ML_NESTED_INNER)
    best_th <- cv$thresholds[which.max(cv$scor)]
    list(
      predict_fn = function(newx) {
        new_spc <- list(x = t(as.matrix(newx[, colnames(x), drop = FALSE])),
                        y = rep(0, nrow(newx)),
                        censoring.status = rep(1, nrow(newx)),
                        featurenames = colnames(x))
        superpc.predict(fit, spc_data, new_spc, threshold = best_th)$v.pred.1df
      },
      model = list(fit = fit, data = spc_data, threshold = best_th)
    )
  }, error = function(e) {
    list(predict_fn = function(newx) rep(NA, nrow(newx)), model = NULL)
  })
}

build_plsrcox <- function(x, y, seed = SEED) {
  set.seed(seed)
  tryCatch({
    fit <- plsRcox(
      Xplan = x, time = y[, "time"], event = y[, "status"],
      nt = min(5, ncol(x) - 1)
    )
    gn <- colnames(x)
    list(
      predict_fn = function(newx) {
        tryCatch(
          as.numeric(predict(fit, newdata = as.matrix(newx[, gn, drop = FALSE]), type = "lp")),
          error = function(e) rep(NA, nrow(newx))
        )
      },
      model = fit
    )
  }, error = function(e) {
    list(predict_fn = function(newx) rep(NA, nrow(newx)), model = NULL)
  })
}

# 评审修订F1/F2/F3：将CoxPH纳入benchmark作为第11个build学习器
# Option 1（保守）：genes-only、无strata，与其他学习器同等条件公平比较；
# IDH分层部署变体在下游单独报告（增量价值）。
build_coxph <- function(x, y, seed = SEED) {
  set.seed(seed)
  tryCatch({
    gn <- colnames(x)
    cn <- make.names(gn, unique = TRUE)
    df <- data.frame(time = y[, "time"], status = y[, "status"], x, check.names = FALSE)
    colnames(df)[-(1:2)] <- cn
    fml <- as.formula(paste("Surv(time, status) ~",
                            paste(sprintf("`%s`", cn), collapse = " + ")))
    fit <- coxph(fml, data = df)
    list(
      predict_fn = function(newx) {
        nd <- as.data.frame(newx[, gn, drop = FALSE], check.names = FALSE)
        colnames(nd) <- cn
        as.numeric(predict(fit, newdata = nd, type = "lp"))
      },
      model = fit
    )
  }, error = function(e) {
    list(predict_fn = function(newx) rep(NA, nrow(newx)), model = NULL)
  })
}

# 算法注册表
FS_ALGOS <- list(
  LASSO       = fs_lasso,
  Ridge       = fs_ridge,
  ElasticNet03 = fs_enet03,
  ElasticNet05 = fs_enet05,
  ElasticNet07 = fs_enet07,
  RSF         = fs_rsf,
  CoxBoost    = fs_coxboost,
  SuperPC     = fs_superpc,
  StepCox     = fs_stepcox,
  plsRcox     = fs_plsrcox
)

BUILD_ALGOS <- list(
  LASSO       = build_lasso,
  Ridge       = build_ridge,
  ElasticNet03 = build_enet03,
  ElasticNet05 = build_enet05,
  ElasticNet07 = build_enet07,
  RSF         = build_rsf,
  CoxBoost    = build_coxboost,
  SuperPC     = build_superpc,
  StepCox     = build_stepcox,
  plsRcox     = build_plsrcox,
  CoxPH       = build_coxph
)

# =============================================================================
# 4. 生成两两组合
# =============================================================================
message("=== Generating algorithm combinations ===")

fs_names    <- names(FS_ALGOS)
build_names <- names(BUILD_ALGOS)

# 生成所有两两组合（包括自身组合，即单独使用一种算法）
combos <- expand.grid(fs = fs_names, build = build_names, stringsAsFactors = FALSE)
combos$combo_name <- paste0(combos$fs, " + ", combos$build)
message(sprintf("Total combinations: %d", nrow(combos)))

# =============================================================================
# 5. 嵌套交叉验证（5x5 nested CV）
# =============================================================================
message("=== Running nested cross-validation ===")

n_outer <- ML_NESTED_OUTER
n_inner <- ML_NESTED_INNER

# 创建外层fold
set.seed(SEED)
outer_folds <- caret::createFolds(train_expr$status, k = n_outer, returnTrain = FALSE)

# 存储所有组合的结果
all_combo_results <- vector("list", nrow(combos))

for (combo_idx in seq_len(nrow(combos))) {
  fs_name    <- combos$fs[combo_idx]
  build_name <- combos$build[combo_idx]
  combo_name <- combos$combo_name[combo_idx]

  message(sprintf("  [%d/%d] %s", combo_idx, nrow(combos), combo_name))

  # 外层CV C-index集合
  outer_cindices <- numeric(n_outer)

  for (fold_i in seq_len(n_outer)) {
    tryCatch({
      test_idx  <- outer_folds[[fold_i]]
      train_idx <- setdiff(seq_len(nrow(train_expr)), test_idx)

      x_tr <- as.matrix(train_expr[train_idx, genes])
      y_tr <- Surv(train_expr$time[train_idx], train_expr$status[train_idx])
      x_te <- as.matrix(train_expr[test_idx, genes])
      y_te <- Surv(train_expr$time[test_idx], train_expr$status[test_idx])

      # 内层：特征选择（在训练fold上）
      fs_result <- FS_ALGOS[[fs_name]](x_tr, y_tr, seed = SEED + fold_i)
      sel_genes <- fs_result$selected
      sel_genes <- intersect(sel_genes, colnames(x_tr))
      if (length(sel_genes) < 2) sel_genes <- colnames(x_tr)

      # 内层：模型构建（在训练fold上，用选中的基因）
      built <- BUILD_ALGOS[[build_name]](
        x_tr[, sel_genes, drop = FALSE], y_tr, seed = SEED + fold_i
      )

      # 外层测试fold预测
      pred_risk <- built$predict_fn(
        data.frame(x_te[, sel_genes, drop = FALSE], check.names = FALSE)
      )

      outer_cindices[fold_i] <- calc_cindex(
        pred_risk, train_expr$time[test_idx], train_expr$status[test_idx]
      )
    }, error = function(e) {
      outer_cindices[fold_i] <<- NA_real_
    })
  }

  # 调优验证集评估（在全训练集上重新训练）
  tuning_cindices <- list()
  tryCatch({
    x_full <- as.matrix(train_expr[, genes])
    y_full <- Surv(train_expr$time, train_expr$status)
    fs_full <- FS_ALGOS[[fs_name]](x_full, y_full, seed = SEED)
    sel_full <- intersect(fs_full$selected, colnames(x_full))
    if (length(sel_full) < 2) sel_full <- colnames(x_full)
    built_full <- BUILD_ALGOS[[build_name]](
      x_full[, sel_full, drop = FALSE], y_full, seed = SEED
    )

    for (tn in names(tuning_data)) {
      td <- tuning_data[[tn]]
      pred_val <- built_full$predict_fn(
        data.frame(td[, sel_full, drop = FALSE], check.names = FALSE)
      )
      tuning_cindices[[tn]] <- calc_cindex(pred_val, td$time, td$status)
    }
  }, error = function(e) {})

  all_combo_results[[combo_idx]] <- list(
    combo_name       = combo_name,
    fs               = fs_name,
    build            = build_name,
    nested_cv_cindex = outer_cindices,
    cv_mean          = mean(outer_cindices, na.rm = TRUE),
    cv_sd            = sd(outer_cindices, na.rm = TRUE),
    tuning_cindex    = tuning_cindices
  )
}

# =============================================================================
# 6. 汇总结果
# =============================================================================
message("=== Summarizing results ===")

ml_results <- data.frame(
  combo_name        = sapply(all_combo_results, `[[`, "combo_name"),
  fs_algorithm      = sapply(all_combo_results, `[[`, "fs"),
  build_algorithm   = sapply(all_combo_results, `[[`, "build"),
  nested_cv_mean    = sapply(all_combo_results, `[[`, "cv_mean"),
  nested_cv_sd      = sapply(all_combo_results, `[[`, "cv_sd"),
  stringsAsFactors  = FALSE
)

# 添加调优验证集C-index
for (tn in names(tuning_data)) {
  ml_results[[paste0("tuning_", tn)]] <- sapply(all_combo_results, function(r) {
    if (tn %in% names(r$tuning_cindex)) r$tuning_cindex[[tn]] else NA_real_
  })
}

# 计算综合排名分数：嵌套CV均值 + 调优验证集均值
tuning_cols <- grep("^tuning_", colnames(ml_results), value = TRUE)
if (length(tuning_cols) > 0) {
  ml_results$tuning_mean <- rowMeans(ml_results[, tuning_cols, drop = FALSE], na.rm = TRUE)
  ml_results$rank_score  <- (ml_results$nested_cv_mean + ml_results$tuning_mean) / 2
} else {
  ml_results$tuning_mean <- NA_real_
  ml_results$rank_score  <- ml_results$nested_cv_mean
}

ml_results <- ml_results %>% dplyr::arrange(desc(rank_score))

write.csv(ml_results,
          file.path(RES_DIR, "ml_combination_results.csv"), row.names = FALSE)

# 最优模型
best_combo <- ml_results$combo_name[1]
best_fs    <- ml_results$fs_algorithm[1]
best_build <- ml_results$build_algorithm[1]
message(sprintf("\nBest combination: %s (rank score = %.4f)",
                best_combo, ml_results$rank_score[1]))

# =============================================================================
# 7. 置换检验（1000次）
# =============================================================================
message(sprintf("=== External permutation test (%d permutations) on CGGA-batch1 for best model ===", ML_NPERM))
# 评审修订F4：置换检验在盲测试集CGGA-batch1上做（外部泛化），而非训练集resubstitution。
# 训练集resubstitution会被RSF等过拟合掩盖真实信号；外部检验中，真实信号能泛化、
# 而置换（噪声）模型不能，因此能正确刻画零分布。该检验自动跟随benchmark选出的胜者。
perm_test_set <- blind_test[["CGGA_batch1"]]
if (is.null(perm_test_set) && length(blind_test) > 0) perm_test_set <- blind_test[[1]]

if (!is.null(perm_test_set)) {
  x_full_perm <- as.matrix(train_expr[, genes])
  y_full_perm <- Surv(train_expr$time, train_expr$status)

  # 观察模型：全TCGA训练，CGGA-batch1评估
  observed_cindex <- tryCatch({
    fs_o  <- FS_ALGOS[[best_fs]](x_full_perm, y_full_perm, seed = SEED)
    sel_o <- intersect(fs_o$selected, colnames(x_full_perm))
    if (length(sel_o) < 2) sel_o <- colnames(x_full_perm)
    blt_o <- BUILD_ALGOS[[best_build]](x_full_perm[, sel_o, drop = FALSE], y_full_perm, seed = SEED)
    pr_o  <- blt_o$predict_fn(data.frame(perm_test_set[, sel_o, drop = FALSE], check.names = FALSE))
    calc_cindex(pr_o, perm_test_set$time, perm_test_set$status)
  }, error = function(e) NA_real_)
  message(sprintf("Observed C-index (CGGA-batch1, external): %.4f", observed_cindex))

  set.seed(SEED)
  perm_cindices <- numeric(ML_NPERM)
  for (perm_i in seq_len(ML_NPERM)) {
    if (perm_i %% 100 == 0) message(sprintf("  Permutation %d/%d", perm_i, ML_NPERM))
    perm_idx <- sample(nrow(train_expr))
    yp <- Surv(train_expr$time[perm_idx], train_expr$status[perm_idx])
    tryCatch({
      fs_p  <- FS_ALGOS[[best_fs]](x_full_perm, yp, seed = SEED + perm_i)
      sel_p <- intersect(fs_p$selected, colnames(x_full_perm))
      if (length(sel_p) < 2) sel_p <- colnames(x_full_perm)
      blt_p <- BUILD_ALGOS[[best_build]](x_full_perm[, sel_p, drop = FALSE], yp, seed = SEED + perm_i)
      pr_p  <- blt_p$predict_fn(data.frame(perm_test_set[, sel_p, drop = FALSE], check.names = FALSE))
      perm_cindices[perm_i] <- calc_cindex(pr_p, perm_test_set$time, perm_test_set$status)
    }, error = function(e) {
      perm_cindices[perm_i] <<- NA_real_
    })
  }
  valid_perm <- perm_cindices[!is.na(perm_cindices)]
  # +1 修正避免 p=0
  perm_pvalue <- (sum(valid_perm >= observed_cindex) + 1) / (length(valid_perm) + 1)
  message(sprintf("External permutation p-value: %.4f (observed=%.4f, null mean=%.4f, n_valid=%d)",
                  perm_pvalue, observed_cindex, mean(valid_perm), length(valid_perm)))
} else {
  message("  No blind test set available for permutation; skipping.")
  observed_cindex <- NA_real_
  perm_cindices   <- rep(NA_real_, ML_NPERM)
  perm_pvalue     <- NA_real_
}

# =============================================================================
# 8. 可视化
# =============================================================================
message("=== Generating visualizations ===")

# --- Fig 5A: C-index热图（特征选择 x 模型构建） ---
heatmap_df <- ml_results %>%
  dplyr::select(fs_algorithm, build_algorithm, rank_score) %>%
  tidyr::pivot_wider(names_from = build_algorithm, values_from = rank_score)

heatmap_mat <- as.matrix(heatmap_df[, -1])
rownames(heatmap_mat) <- heatmap_df$fs_algorithm

pdf(file.path(FIG_DIR, "Fig5A_cindex_heatmap.pdf"), width = 12, height = 10)
par(mar = c(8, 8, 4, 4))
heatmap(heatmap_mat, Rowv = NA, Colv = NA,
        col = colorRampPalette(c("white", "#4DBBD5", "#E64B35"))(100),
        scale = "none", margins = c(10, 10),
        main = "C-index: Feature Selection x Model Building",
        xlab = "Model Building", ylab = "Feature Selection")
dev.off()

# --- ggplot版热图 ---
heatmap_long <- ml_results %>%
  dplyr::select(fs_algorithm, build_algorithm, rank_score)

p_heatmap <- ggplot(heatmap_long, aes(x = build_algorithm, y = fs_algorithm, fill = rank_score)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.3f", rank_score)), size = 2.5) +
  scale_fill_gradient2(low = "white", mid = "#4DBBD5", high = "#E64B35",
                       midpoint = median(heatmap_long$rank_score, na.rm = TRUE),
                       name = "Rank\nScore") +
  labs(x = "Model Building Algorithm", y = "Feature Selection Algorithm",
       title = "ML Combination C-index Heatmap") +
  THEME_PUBLICATION +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        axis.text.y = element_text(size = 8))

ggsave(file.path(FIG_DIR, "Fig5A_cindex_heatmap_gg.pdf"), p_heatmap,
       width = 12, height = 10)

# --- Fig 5B: 嵌套CV C-index小提琴图（top 20组合） ---
top20 <- ml_results$combo_name[1:min(20, nrow(ml_results))]

violin_list <- list()
for (r in all_combo_results) {
  if (r$combo_name %in% top20) {
    violin_list[[length(violin_list) + 1]] <- data.frame(
      combo = r$combo_name,
      cindex = r$nested_cv_cindex,
      stringsAsFactors = FALSE
    )
  }
}
violin_df <- do.call(rbind, violin_list)
violin_df$combo <- factor(violin_df$combo, levels = rev(top20))

p_violin <- ggplot(violin_df, aes(x = combo, y = cindex)) +
  geom_violin(fill = "#4DBBD5", alpha = 0.6) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA) +
  geom_jitter(width = 0.05, size = 1, alpha = 0.5) +
  coord_flip() +
  labs(x = "", y = "C-index (Nested CV)",
       title = "Top 20 ML Combinations: Nested CV Performance") +
  THEME_PUBLICATION +
  theme(axis.text.y = element_text(size = 7))

ggsave(file.path(FIG_DIR, "Fig5B_nested_cv_violin.pdf"), p_violin,
       width = 10, height = 8)

# --- 全部组合的小提琴图（Supplementary） ---
all_violin_list <- list()
for (r in all_combo_results) {
  all_violin_list[[length(all_violin_list) + 1]] <- data.frame(
    combo  = r$combo_name,
    cindex = r$nested_cv_cindex,
    stringsAsFactors = FALSE
  )
}
all_violin_df <- do.call(rbind, all_violin_list)

p_all_dist <- ggplot(all_violin_df, aes(x = cindex)) +
  geom_histogram(bins = 50, fill = "#4DBBD5", color = "white", alpha = 0.7) +
  geom_vline(xintercept = ml_results$nested_cv_mean[1], color = "#E64B35",
             linetype = "dashed", linewidth = 1) +
  labs(x = "C-index", y = "Count",
       title = "Distribution of All ML Combination C-indices",
       subtitle = sprintf("Best: %s (%.4f)", best_combo, ml_results$nested_cv_mean[1])) +
  THEME_PUBLICATION

ggsave(file.path(FIG_DIR, "SuppFig_ml_cindex_distribution.pdf"), p_all_dist,
       width = 8, height = 5)

# --- 置换检验分布图 ---
perm_df <- data.frame(cindex = perm_cindices[!is.na(perm_cindices)])

p_perm <- ggplot(perm_df, aes(x = cindex)) +
  geom_histogram(bins = 50, fill = "grey70", color = "white") +
  geom_vline(xintercept = observed_cindex, color = "#E64B35",
             linetype = "dashed", linewidth = 1.2) +
  annotate("text", x = observed_cindex, y = Inf, vjust = 2,
           label = sprintf("Observed = %.4f\np = %.4f", observed_cindex, perm_pvalue),
           color = "#E64B35", hjust = -0.1, size = 4) +
  labs(x = "C-index (Permuted)", y = "Count",
       title = sprintf("Permutation Test (%d permutations)", ML_NPERM)) +
  THEME_PUBLICATION

ggsave(file.path(FIG_DIR, "SuppFig_permutation_test.pdf"), p_perm,
       width = 8, height = 5)

# =============================================================================
# 9. 盲测试集最终评估（仅对最优模型执行一次）
# =============================================================================
message("=== Blind test set evaluation (best model only) ===")

blind_results <- list()

tryCatch({
  x_full <- as.matrix(train_expr[, genes])
  y_full <- Surv(train_expr$time, train_expr$status)

  fs_final <- FS_ALGOS[[best_fs]](x_full, y_full, seed = SEED)
  sel_final <- intersect(fs_final$selected, colnames(x_full))
  if (length(sel_final) < 2) sel_final <- colnames(x_full)

  built_final <- BUILD_ALGOS[[best_build]](
    x_full[, sel_final, drop = FALSE], y_full, seed = SEED
  )

  for (bt_name in names(blind_test)) {
    bt <- blind_test[[bt_name]]
    pred_bt <- built_final$predict_fn(
      data.frame(bt[, sel_final, drop = FALSE], check.names = FALSE)
    )
    c_bt <- calc_cindex(pred_bt, bt$time, bt$status)
    blind_results[[bt_name]] <- c_bt
    message(sprintf("  %s (blind test): C-index = %.4f", bt_name, c_bt))
  }

  # 探索性
  for (ex_name in names(exploratory)) {
    ex <- exploratory[[ex_name]]
    pred_ex <- built_final$predict_fn(
      data.frame(ex[, sel_final, drop = FALSE], check.names = FALSE)
    )
    c_ex <- calc_cindex(pred_ex, ex$time, ex$status)
    blind_results[[paste0(ex_name, "_exploratory")]] <- c_ex
    message(sprintf("  %s (exploratory): C-index = %.4f", ex_name, c_ex))
  }
}, error = function(e) {
  message("  Blind test evaluation failed: ", e$message)
})

# =============================================================================
# 10. 保存
# =============================================================================
save(
  ml_results,
  all_combo_results,
  best_combo, best_fs, best_build,
  perm_cindices, perm_pvalue, observed_cindex,
  blind_results,
  train_expr, validation_data, tuning_data, blind_test, exploratory,
  genes,
  FS_ALGOS, BUILD_ALGOS,
  file = file.path(DATA_PROC, "ml_combination_results.RData")
)

write.csv(
  data.frame(
    metric = c("best_combo", "nested_cv_mean", "perm_pvalue",
               names(blind_results)),
    value  = c(best_combo,
               sprintf("%.4f", ml_results$rank_score[1]),
               sprintf("%.4f", perm_pvalue),
               sapply(blind_results, function(x) sprintf("%.4f", x)))
  ),
  file.path(RES_DIR, "best_model_summary.csv"), row.names = FALSE
)

message("\n=== ML combination analysis completed ===")
message(sprintf("Best: %s (rank=%.4f, perm p=%.4f)", best_combo,
                ml_results$rank_score[1], perm_pvalue))
message("Next step: Run 05_ml_model/03_risk_score.R")
