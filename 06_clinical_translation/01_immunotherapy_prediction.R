###############################################################################
# 01_immunotherapy_prediction.R
# Part 4 - 临床转化 (1/4)
# 免疫治疗应答预测：TIDE + IPS + 多维免疫功能评分 + 真实ICB队列验证
# 输出: Figure 7A-C
###############################################################################

source("00_setup/config.R")
library(ggplot2)
library(ggpubr)
library(survival)
library(survminer)
library(dplyr)
library(tidyr)
library(GSVA)
library(cowplot)

set.seed(SEED)

# --- 加载数据 ---
load(file.path(DATA_PROC, "risk_model_final.RData"))
load(file.path(DATA_PROC, "tcga_glioma_expression.RData"))
load(file.path(DATA_PROC, "upr_gene_sets.RData"))

# --- 创建输出子目录 ---
fig_dir <- file.path(FIG_DIR, "part4_clinical")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 0. 准备公共数据
# =============================================================================
message("=== Preparing common data ===")

common_samples <- intersect(names(risk_score_train), colnames(expr_tpm_symbol))
risk_groups <- setNames(risk_group_train[common_samples], common_samples)
expr_log <- log2(expr_tpm_symbol[, common_samples] + 1)

message(sprintf("Common samples: %d (High: %d, Low: %d)",
                length(common_samples),
                sum(risk_groups == "High"), sum(risk_groups == "Low")))

# =============================================================================
# 1. TIDE评分 -- 程序化计算（tidepy via reticulate）
# =============================================================================
message("\n=== Section 1: TIDE Score ===")

# TIDE输入: log2(TPM+1) 行均值中心化
tide_expr <- expr_log[, common_samples]
tide_expr_centered <- tide_expr - rowMeans(tide_expr)

# 保存TIDE输入文件（供手动上传或Python调用）
tide_input_file <- file.path(RES_DIR, "tide_input.txt")
write.table(tide_expr_centered, file = tide_input_file,
            sep = "\t", quote = FALSE, col.names = NA)
message("TIDE input matrix saved: ", tide_input_file)

# 尝试通过reticulate调用tidepy
tide_computed <- FALSE
tryCatch({
  library(reticulate)

  # 检查tidepy是否可用
  if (py_module_available("tidepy")) {
    message("tidepy found, computing TIDE scores programmatically...")
    tidepy <- import("tidepy")
    pd <- import("pandas")

    # 转为python对象
    tide_mat_py <- r_to_py(as.data.frame(tide_expr_centered))

    # 运行TIDE
    tide_result_py <- tidepy$predict_tide(tide_mat_py, cancer_type = "GBM")
    tide_results <- py_to_r(tide_result_py)
    tide_computed <- TRUE
    message("TIDE computation completed via tidepy")

  } else {
    message("tidepy not installed. Attempting pip install...")
    tryCatch({
      py_install("tidepy", pip = TRUE)
      tidepy <- import("tidepy")
      pd <- import("pandas")
      tide_mat_py <- r_to_py(as.data.frame(tide_expr_centered))
      tide_result_py <- tidepy$predict_tide(tide_mat_py, cancer_type = "GBM")
      tide_results <- py_to_r(tide_result_py)
      tide_computed <- TRUE
    }, error = function(e2) {
      message("tidepy installation failed: ", e2$message)
    })
  }
}, error = function(e) {
  message("reticulate/tidepy not available: ", e$message)
})

# 备选方案: 生成独立Python脚本供离线运行
if (!tide_computed) {
  message("Generating standalone Python script for TIDE computation...")

  tide_py_script <- file.path(RES_DIR, "run_tide.py")
  writeLines(c(
    "#!/usr/bin/env python3",
    "\"\"\"Run TIDE prediction on glioma expression data.\"\"\"",
    "import pandas as pd",
    "try:",
    "    import tidepy",
    "except ImportError:",
    "    import subprocess, sys",
    "    subprocess.check_call([sys.executable, '-m', 'pip', 'install', 'tidepy'])",
    "    import tidepy",
    "",
    paste0("expr = pd.read_csv('", tide_input_file, "', sep='\\t', index_col=0)"),
    "result = tidepy.predict_tide(expr, cancer_type='GBM')",
    paste0("result.to_csv('", file.path(RES_DIR, "tide_result.csv"), "')"),
    "print('TIDE results saved successfully')",
    "print(result.head())"
  ), tide_py_script)
  message("Python script saved: ", tide_py_script)
  message("Run manually: python3 ", tide_py_script)

  # 尝试自动执行
  tryCatch({
    ret <- system2("python3", args = tide_py_script, timeout = 300,
                   stdout = TRUE, stderr = TRUE)
    tide_result_file <- file.path(RES_DIR, "tide_result.csv")
    if (file.exists(tide_result_file)) {
      tide_results <- read.csv(tide_result_file, row.names = 1)
      tide_computed <- TRUE
      message("TIDE computed via standalone Python script")
    }
  }, error = function(e) {
    message("Standalone Python execution failed: ", e$message)
  })
}

# 尝试从已有结果文件加载
if (!tide_computed) {
  for (tf in c(file.path(RES_DIR, "tide_result.csv"),
               file.path(DATA_RAW, "tide_result.csv"))) {
    if (file.exists(tf)) {
      tide_results <- read.csv(tf, row.names = 1)
      tide_computed <- TRUE
      message("TIDE results loaded from: ", tf)
      break
    }
  }
}

# --- TIDE可视化 ---
tide_df <- NULL
if (tide_computed) {
  tide_col <- grep("^TIDE$", colnames(tide_results), value = TRUE)
  if (length(tide_col) == 0) tide_col <- colnames(tide_results)[1]

  tide_df <- data.frame(
    barcode = rownames(tide_results),
    TIDE = as.numeric(tide_results[[tide_col]]),
    stringsAsFactors = FALSE
  )

  common_tide <- intersect(tide_df$barcode, names(risk_groups))
  tide_df <- tide_df[tide_df$barcode %in% common_tide, ]
  tide_df$risk_group <- factor(risk_groups[tide_df$barcode], levels = c("Low", "High"))

  p_tide <- ggplot(tide_df, aes(x = risk_group, y = TIDE, fill = risk_group)) +
    geom_boxplot(outlier.size = 0.5, width = 0.6) +
    geom_jitter(width = 0.15, size = 0.3, alpha = 0.3) +
    stat_compare_means(method = "wilcox.test", label = "p.format",
                       label.x.npc = 0.5, size = 4) +
    scale_fill_manual(values = COLORS_RISK) +
    labs(x = "Risk Group", y = "TIDE Score",
         title = "TIDE Score by Risk Group") +
    THEME_PUBLICATION +
    theme(legend.position = "none")

  ggsave(file.path(fig_dir, "Fig7A_tide_score.pdf"), p_tide, width = 5, height = 5)

  # TIDE应答预测
  resp_col <- grep("Responder|Response|responder", colnames(tide_results), value = TRUE)
  if (length(resp_col) > 0) {
    tide_df$Responder <- tide_results[tide_df$barcode, resp_col[1]]

    resp_summary <- tide_df %>%
      dplyr::filter(!is.na(Responder)) %>%
      dplyr::count(risk_group, Responder) %>%
      dplyr::group_by(risk_group) %>%
      dplyr::mutate(pct = n / sum(n) * 100)

    p_resp <- ggplot(resp_summary, aes(x = risk_group, y = pct, fill = Responder)) +
      geom_bar(stat = "identity", position = "stack", width = 0.6) +
      labs(x = "Risk Group", y = "Percentage (%)",
           title = "TIDE Predicted Response") +
      scale_fill_brewer(palette = "Set2") +
      THEME_PUBLICATION

    ggsave(file.path(fig_dir, "Fig7A_tide_response.pdf"), p_resp, width = 5.5, height = 5)
  }

  tide_stat <- wilcox.test(TIDE ~ risk_group, data = tide_df)
  message(sprintf("TIDE: High vs Low, Wilcoxon p = %.2e", tide_stat$p.value))
} else {
  message("WARNING: TIDE scores not available.")
  message("Please run: python3 ", file.path(RES_DIR, "run_tide.py"))
  message("Or upload tide_input.txt to http://tide.dfci.harvard.edu/")
}

# =============================================================================
# 2. IPS (Immunophenoscore) -- 四维度计算
# =============================================================================
message("\n=== Section 2: Immunophenoscore (IPS) ===")

# Charoentong et al. 2017 Cell Reports
ips_gene_sets <- list(
  MHC = c("HLA-A", "HLA-B", "HLA-C",
          "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1",
          "HLA-DRA", "HLA-DRB1",
          "B2M", "TAP1", "TAP2", "TAPBP"),
  Effector_cells = c("CD8A", "CD8B", "GZMA", "GZMB", "GZMK", "GZMH",
                     "PRF1", "IFNG", "TBX21", "EOMES",
                     "NKG7", "GNLY", "KLRK1", "KLRD1"),
  Suppressor_cells = c("FOXP3", "CTLA4", "PDCD1", "LAG3", "HAVCR2",
                       "TIGIT", "CD274", "PDCD1LG2",
                       "IDO1", "ARG1", "TGFB1", "IL10"),
  Checkpoints = c("TNFRSF4", "TNFRSF9", "CD28", "ICOS",
                  "TNFRSF18", "CD40", "CD27",
                  "CD80", "CD86", "TNFSF4", "TNFSF9")
)

# 过滤存在于表达矩阵中的基因
ips_gene_sets_valid <- lapply(ips_gene_sets, function(gs) {
  gs[gs %in% rownames(expr_log)]
})

for (nm in names(ips_gene_sets_valid)) {
  message(sprintf("  IPS %s: %d/%d genes found",
                  nm, length(ips_gene_sets_valid[[nm]]), length(ips_gene_sets[[nm]])))
}

# ssGSEA计算
ips_params <- ssgseaParam(
  exprData = as.matrix(expr_log[, common_samples]),
  geneSets = ips_gene_sets_valid
)
ips_scores <- gsva(ips_params)

# 综合IPS: MHC + Effector + Checkpoints - Suppressor
ips_composite <- as.numeric(ips_scores["MHC", ]) +
  as.numeric(ips_scores["Effector_cells", ]) +
  as.numeric(ips_scores["Checkpoints", ]) -
  as.numeric(ips_scores["Suppressor_cells", ])

ips_df <- data.frame(
  barcode = colnames(ips_scores),
  IPS = ips_composite,
  MHC = as.numeric(ips_scores["MHC", ]),
  Effector = as.numeric(ips_scores["Effector_cells", ]),
  Suppressor = as.numeric(ips_scores["Suppressor_cells", ]),
  Checkpoint = as.numeric(ips_scores["Checkpoints", ]),
  risk_group = factor(risk_groups[colnames(ips_scores)], levels = c("Low", "High")),
  stringsAsFactors = FALSE
)
ips_df <- ips_df[!is.na(ips_df$risk_group), ]

# 综合IPS箱线图
p_ips_main <- ggplot(ips_df, aes(x = risk_group, y = IPS, fill = risk_group)) +
  geom_boxplot(outlier.size = 0.5, width = 0.6) +
  geom_jitter(width = 0.15, size = 0.3, alpha = 0.3) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 4) +
  scale_fill_manual(values = COLORS_RISK) +
  labs(x = "Risk Group", y = "Immunophenoscore (IPS)",
       title = "IPS: High vs Low Risk") +
  THEME_PUBLICATION +
  theme(legend.position = "none")

# IPS四维度分面
ips_long <- ips_df %>%
  tidyr::pivot_longer(cols = c("MHC", "Effector", "Suppressor", "Checkpoint"),
                      names_to = "Component", values_to = "Score") %>%
  dplyr::mutate(Component = factor(Component,
                                   levels = c("MHC", "Effector", "Suppressor", "Checkpoint")))

p_ips_comp <- ggplot(ips_long, aes(x = risk_group, y = Score, fill = risk_group)) +
  geom_boxplot(outlier.size = 0.3, width = 0.6) +
  facet_wrap(~Component, scales = "free_y", nrow = 1) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 3.5) +
  scale_fill_manual(values = COLORS_RISK) +
  labs(x = "Risk Group", y = "ssGSEA Score") +
  THEME_PUBLICATION +
  theme(legend.position = "none",
        strip.text = element_text(size = 11, face = "bold"))

p_ips_combined <- plot_grid(p_ips_main, p_ips_comp,
                            ncol = 1, rel_heights = c(1, 0.8))
ggsave(file.path(fig_dir, "Fig7B_ips_score.pdf"), p_ips_combined, width = 12, height = 9)

ips_stat <- wilcox.test(IPS ~ risk_group, data = ips_df)
message(sprintf("IPS: High vs Low, Wilcoxon p = %.2e", ips_stat$p.value))
for (comp in c("MHC", "Effector", "Suppressor", "Checkpoint")) {
  stat_tmp <- wilcox.test(as.formula(paste0(comp, " ~ risk_group")), data = ips_df)
  message(sprintf("  %s: p = %.2e", comp, stat_tmp$p.value))
}

# =============================================================================
# 3. 多维免疫功能评分（7个维度）
# =============================================================================
message("\n=== Section 3: Multi-dimensional Immune Function Scores ===")

immune_functions <- list(
  Cytolytic_activity = c("GZMA", "GZMB", "PRF1", "GNLY", "NKG7",
                         "GZMK", "GZMH", "KLRK1"),
  Antigen_presentation = c("HLA-A", "HLA-B", "HLA-C", "B2M",
                           "TAP1", "TAP2", "TAPBP",
                           "PSMB8", "PSMB9", "PSME1", "PSME2"),
  T_cell_activation = c("CD3D", "CD3E", "CD3G", "CD28", "ICOS",
                         "LCK", "ZAP70", "ITK", "CD247"),
  IFN_gamma_response = c("IFNG", "IFNGR1", "IFNGR2", "STAT1",
                          "IRF1", "GBP1", "GBP2", "CXCL9", "CXCL10"),
  Immune_exclusion = c("TGFB1", "TGFB2", "TGFB3",
                       "VEGFA", "VEGFB", "IL10", "WNT5A"),
  T_cell_exhaustion = c("PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4",
                         "TOX", "TOX2", "ENTPD1", "LAYN"),
  Immunosuppression = c("CD274", "IDO1", "ARG1", "IL10", "TGFB1",
                         "FOXP3", "CD33", "SIGLEC15")
)

# 过滤基因
immune_functions_valid <- lapply(immune_functions, function(gs) {
  valid <- gs[gs %in% rownames(expr_log)]
  if (length(valid) < 3) {
    warning(sprintf("Gene set has only %d genes available", length(valid)))
  }
  valid
})

imm_func_params <- ssgseaParam(
  exprData = as.matrix(expr_log[, common_samples]),
  geneSets = immune_functions_valid
)
imm_func_scores <- gsva(imm_func_params)

imm_func_df <- as.data.frame(t(imm_func_scores))
imm_func_df$barcode <- rownames(imm_func_df)
imm_func_df$risk_group <- factor(risk_groups[rownames(imm_func_df)], levels = c("Low", "High"))
imm_func_df <- imm_func_df[!is.na(imm_func_df$risk_group), ]

imm_func_long <- imm_func_df %>%
  tidyr::pivot_longer(cols = -c(barcode, risk_group),
                      names_to = "Function", values_to = "Score") %>%
  dplyr::mutate(Function = gsub("_", " ", Function))

# Wilcoxon检验
imm_func_stats <- imm_func_long %>%
  dplyr::group_by(Function) %>%
  dplyr::summarise(
    median_high = median(Score[risk_group == "High"]),
    median_low  = median(Score[risk_group == "Low"]),
    p_value = wilcox.test(Score ~ risk_group)$p.value,
    .groups = "drop"
  ) %>%
  dplyr::mutate(padj = p.adjust(p_value, method = "BH"))

write.csv(imm_func_stats,
          file.path(RES_DIR, "immune_function_risk_comparison.csv"),
          row.names = FALSE)

p_imm_func <- ggplot(imm_func_long, aes(x = risk_group, y = Score, fill = risk_group)) +
  geom_boxplot(outlier.size = 0.3, width = 0.6) +
  facet_wrap(~Function, scales = "free_y", ncol = 4) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 3) +
  scale_fill_manual(values = COLORS_RISK) +
  labs(x = "Risk Group", y = "ssGSEA Enrichment Score") +
  THEME_PUBLICATION +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 9, face = "bold"))

ggsave(file.path(fig_dir, "Fig7B_immune_functions.pdf"), p_imm_func, width = 14, height = 8)

message("Immune function comparison:")
print(imm_func_stats)

# =============================================================================
# 4. 真实ICB队列验证
# =============================================================================
message("\n=== Section 4: Real ICB Cohort Validation ===")

# 辅助函数: 在外部ICB队列中计算UIRS并做生存/应答分析
compute_uirs_external <- function(expr_mat, clin_df,
                                  time_col, event_col,
                                  sample_id_col = NULL,
                                  response_col = NULL,
                                  dataset_name = "External") {

  avail_genes <- intersect(model_genes, rownames(expr_mat))
  if (length(avail_genes) < length(model_genes) * 0.5) {
    message(sprintf("  %s: only %d/%d model genes, skipping",
                    dataset_name, length(avail_genes), length(model_genes)))
    return(NULL)
  }
  if (length(avail_genes) < length(model_genes)) {
    message(sprintf("  %s: %d/%d model genes (missing: %s)",
                    dataset_name, length(avail_genes), length(model_genes),
                    paste(setdiff(model_genes, avail_genes), collapse = ", ")))
  }

  # 匹配样本
  if (!is.null(sample_id_col)) {
    common_ids <- intersect(colnames(expr_mat), clin_df[[sample_id_col]])
    clin_matched <- clin_df[match(common_ids, clin_df[[sample_id_col]]), ]
  } else {
    common_ids <- intersect(colnames(expr_mat), rownames(clin_df))
    clin_matched <- clin_df[common_ids, ]
  }

  if (length(common_ids) < 10) {
    message(sprintf("  %s: only %d matched samples, skipping", dataset_name, length(common_ids)))
    return(NULL)
  }

  # 计算风险评分
  x_ext <- t(expr_mat[avail_genes, common_ids])
  weights_avail <- model_weights[avail_genes]
  risk_score_ext <- as.numeric(x_ext %*% weights_avail)
  names(risk_score_ext) <- common_ids
  risk_group_ext <- ifelse(risk_score_ext > median(risk_score_ext), "High", "Low")

  result <- list(
    dataset = dataset_name,
    risk_score = risk_score_ext,
    risk_group = risk_group_ext,
    n_samples = length(common_ids),
    n_genes = length(avail_genes)
  )

  # KM生存分析
  os_time  <- as.numeric(clin_matched[[time_col]])
  os_event <- as.numeric(clin_matched[[event_col]])
  valid_idx <- !is.na(os_time) & !is.na(os_event) & os_time > 0

  if (sum(valid_idx) >= 10) {
    surv_ext <- data.frame(
      time = os_time[valid_idx],
      status = os_event[valid_idx],
      risk_score = risk_score_ext[common_ids[valid_idx]],
      risk_group = factor(risk_group_ext[common_ids[valid_idx]], levels = c("Low", "High"))
    )

    fit_ext <- survfit(Surv(time, status) ~ risk_group, data = surv_ext)
    p_km <- ggsurvplot(
      fit_ext, data = surv_ext,
      pval = TRUE, pval.method = TRUE,
      risk.table = TRUE,
      palette = c(COLORS_RISK["Low"], COLORS_RISK["High"]),
      xlab = "Time (months)",
      title = paste0("ICB Cohort: ", dataset_name),
      ggtheme = THEME_PUBLICATION,
      risk.table.height = 0.25
    )

    pdf(file.path(fig_dir, paste0("Fig7C_km_", gsub("[^a-zA-Z0-9]", "_", dataset_name), ".pdf")),
        width = 8, height = 7)
    print(p_km)
    dev.off()

    result$surv_data <- surv_ext
    result$log_rank_p <- surv_pvalue(fit_ext)$pval
    message(sprintf("  %s KM: log-rank p = %.3e", dataset_name, result$log_rank_p))
  }

  # 治疗应答分析
  if (!is.null(response_col) && response_col %in% colnames(clin_matched)) {
    resp <- clin_matched[[response_col]]
    valid_resp <- !is.na(resp)

    if (sum(valid_resp) >= 5) {
      resp_df <- data.frame(
        risk_group = factor(risk_group_ext[common_ids[valid_resp]], levels = c("Low", "High")),
        response = resp[valid_resp]
      )

      resp_table <- table(resp_df$risk_group, resp_df$response)
      fisher_p <- tryCatch(fisher.test(resp_table)$p.value, error = function(e) NA)
      result$response_table <- resp_table
      result$fisher_p <- fisher_p
      message(sprintf("  %s response Fisher p = %.3e", dataset_name, fisher_p))

      # 应答率条形图
      resp_pct <- resp_df %>%
        dplyr::count(risk_group, response) %>%
        dplyr::group_by(risk_group) %>%
        dplyr::mutate(pct = n / sum(n) * 100)

      p_resp_bar <- ggplot(resp_pct, aes(x = risk_group, y = pct, fill = response)) +
        geom_bar(stat = "identity", position = "stack", width = 0.6) +
        labs(x = "Risk Group", y = "Percentage (%)",
             title = paste0(dataset_name, ": Response"), fill = "Response") +
        scale_fill_brewer(palette = "Set2") +
        THEME_PUBLICATION

      ggsave(file.path(fig_dir, paste0("Fig7C_response_",
                                        gsub("[^a-zA-Z0-9]", "_", dataset_name), ".pdf")),
             p_resp_bar, width = 5.5, height = 5)

      # 风险评分 vs 应答
      score_resp <- data.frame(
        risk_score = risk_score_ext[common_ids[valid_resp]],
        response = resp[valid_resp]
      )

      p_score_resp <- ggplot(score_resp, aes(x = response, y = risk_score, fill = response)) +
        geom_boxplot(outlier.size = 0.5, width = 0.6) +
        geom_jitter(width = 0.15, size = 1, alpha = 0.5) +
        stat_compare_means(method = "wilcox.test", label = "p.format", size = 4) +
        labs(x = "Response", y = "UIRS Risk Score",
             title = paste0(dataset_name, ": Score vs Response")) +
        scale_fill_brewer(palette = "Set2") +
        THEME_PUBLICATION + theme(legend.position = "none")

      ggsave(file.path(fig_dir, paste0("Fig7C_score_resp_",
                                        gsub("[^a-zA-Z0-9]", "_", dataset_name), ".pdf")),
             p_score_resp, width = 5.5, height = 5)
    }
  }

  return(result)
}

# --- 4A. Zhao et al. 2019 (Nature Medicine) - 66 GBM, anti-PD1 ---
icb_results <- list()

zhao_file <- file.path(DATA_RAW, "zhao2019_icb_gbm.RData")
if (file.exists(zhao_file)) {
  load(zhao_file)
  message("Zhao 2019 data loaded")
  if (exists("zhao_expr") && exists("zhao_clin")) {
    icb_results[["Zhao2019"]] <- compute_uirs_external(
      expr_mat = zhao_expr, clin_df = zhao_clin,
      time_col = "OS.time", event_col = "OS",
      response_col = "Response",
      dataset_name = "Zhao2019_GBM_antiPD1_n66"
    )
  }
} else {
  message("Zhao 2019 ICB data not found: ", zhao_file)
  message("  Expected: zhao_expr (gene x sample), zhao_clin (with OS.time, OS, Response)")
}

# --- 4B. Cloughesy et al. 2019 (Nature Medicine) - 35 rGBM ---
cloughesy_file <- file.path(DATA_RAW, "cloughesy2019_icb_gbm.RData")
if (file.exists(cloughesy_file)) {
  load(cloughesy_file)
  message("Cloughesy 2019 data loaded")
  if (exists("cloughesy_expr") && exists("cloughesy_clin")) {
    icb_results[["Cloughesy2019"]] <- compute_uirs_external(
      expr_mat = cloughesy_expr, clin_df = cloughesy_clin,
      time_col = "OS.time", event_col = "OS",
      response_col = "Response",
      dataset_name = "Cloughesy2019_rGBM_n35"
    )
  }
} else {
  message("Cloughesy 2019 data not found: ", cloughesy_file)
}

# --- 4C. IMvigor210 (辅助, 膀胱癌 anti-PD-L1) ---
imvigor_loaded <- FALSE
tryCatch({
  if (requireNamespace("IMvigor210CoreBiologies", quietly = TRUE)) {
    library(IMvigor210CoreBiologies)
    data("cds")
    imvigor_counts <- Biobase::exprs(cds)
    imvigor_clin   <- Biobase::pData(cds)
    lib_sizes <- colSums(imvigor_counts)
    imvigor_cpm <- sweep(imvigor_counts, 2, lib_sizes, "/") * 1e6
    imvigor_log <- log2(imvigor_cpm + 1)

    icb_results[["IMvigor210"]] <- compute_uirs_external(
      expr_mat = imvigor_log, clin_df = imvigor_clin,
      time_col = "os", event_col = "censOS",
      response_col = "Best.Confirmed.Overall.Response",
      dataset_name = "IMvigor210_Bladder_antiPDL1"
    )
    imvigor_loaded <- TRUE
  }
}, error = function(e) {
  message("IMvigor210 package not available: ", e$message)
})

if (!imvigor_loaded) {
  imvigor_file <- file.path(DATA_RAW, "imvigor210_data.RData")
  if (file.exists(imvigor_file)) {
    load(imvigor_file)
    if (exists("imvigor_expr") && exists("imvigor_clin")) {
      icb_results[["IMvigor210"]] <- compute_uirs_external(
        expr_mat = imvigor_expr, clin_df = imvigor_clin,
        time_col = "os", event_col = "censOS",
        response_col = "Best.Confirmed.Overall.Response",
        dataset_name = "IMvigor210_Bladder_antiPDL1"
      )
    }
  } else {
    message("IMvigor210 data not found: ", imvigor_file)
  }
}

# ICB验证汇总
if (length(icb_results) > 0) {
  icb_summary <- do.call(rbind, lapply(icb_results, function(x) {
    data.frame(
      Dataset = x$dataset,
      N_samples = x$n_samples,
      N_genes = x$n_genes,
      Log_rank_p = ifelse(!is.null(x$log_rank_p), x$log_rank_p, NA),
      Fisher_p = ifelse(!is.null(x$fisher_p), x$fisher_p, NA),
      stringsAsFactors = FALSE
    )
  }))
  write.csv(icb_summary, file.path(RES_DIR, "icb_validation_summary.csv"), row.names = FALSE)
  message("\nICB validation summary:")
  print(icb_summary)
} else {
  message("\nNo ICB cohort data available for validation.")
  message("Required files:")
  message("  - ", zhao_file)
  message("  - ", cloughesy_file)
  message("  - IMvigor210CoreBiologies R package or ", file.path(DATA_RAW, "imvigor210_data.RData"))
}

# =============================================================================
# 5. 高低风险组免疫治疗关键基因比较
# =============================================================================
message("\n=== Section 5: Key Immunotherapy Gene Comparison ===")

key_immuno_genes <- c("CD274", "PDCD1", "CTLA4", "HAVCR2", "LAG3",
                      "TIGIT", "IFNG", "GZMA", "GZMB", "PRF1")
key_valid <- key_immuno_genes[key_immuno_genes %in% rownames(expr_log)]

immuno_genes_df <- data.frame(
  barcode = common_samples,
  risk_group = factor(risk_groups[common_samples], levels = c("Low", "High")),
  stringsAsFactors = FALSE
)
for (gene in key_valid) {
  immuno_genes_df[[gene]] <- as.numeric(expr_log[gene, common_samples])
}

immuno_genes_long <- immuno_genes_df %>%
  tidyr::pivot_longer(cols = -c(barcode, risk_group),
                      names_to = "Gene", values_to = "Expression")

p_immuno_genes <- ggplot(immuno_genes_long, aes(x = risk_group, y = Expression, fill = risk_group)) +
  geom_boxplot(outlier.size = 0.3, width = 0.6) +
  facet_wrap(~Gene, scales = "free_y", ncol = 5) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 3) +
  scale_fill_manual(values = COLORS_RISK) +
  labs(x = "Risk Group", y = "log2(TPM+1)") +
  THEME_PUBLICATION +
  theme(legend.position = "none",
        strip.text = element_text(size = 10, face = "italic"))

ggsave(file.path(fig_dir, "Fig7A_immunotherapy_genes.pdf"), p_immuno_genes,
       width = 14, height = 6)

# =============================================================================
# 6. 保存全部结果
# =============================================================================
message("\n=== Saving results ===")

save(ips_df, ips_scores, imm_func_scores, imm_func_stats,
     tide_df, icb_results,
     file = file.path(DATA_PROC, "immunotherapy_prediction_results.RData"))

write.csv(ips_df, file.path(RES_DIR, "ips_scores.csv"), row.names = FALSE)

message("\n=== Immunotherapy prediction analysis completed ===")
message("Figures saved to: ", fig_dir)
message("  Fig7A: TIDE score + immunotherapy genes")
message("  Fig7B: IPS (4 dimensions) + immune function (7 dimensions)")
message("  Fig7C: ICB cohort validation")
message("\nNext: Run 06_clinical_translation/02_drug_sensitivity.R")
