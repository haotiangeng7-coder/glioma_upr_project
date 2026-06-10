###############################################################################
# 02_drug_sensitivity.R
# Part 4 - 临床转化 (2/4)
# 药物敏感性分析: oncoPredict GDSC2 IC50 + 重点药物 + 靶点通路活性
# 输出: Figure 7D
###############################################################################

source("00_setup/config.R")
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(GSVA)
library(ComplexHeatmap)
library(circlize)

set.seed(SEED)

# --- 加载数据 ---
load(file.path(DATA_PROC, "risk_model_final.RData"))
load(file.path(DATA_PROC, "tcga_glioma_expression.RData"))

# --- 输出目录 ---
fig_dir <- file.path(FIG_DIR, "part4_clinical")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. 准备数据
# =============================================================================
message("=== Preparing drug sensitivity data ===")

common_samples <- intersect(names(risk_score_train), colnames(expr_tpm_symbol))
risk_groups <- setNames(risk_group_train[common_samples], common_samples)

# 表达矩阵 (TPM, 不取log, oncoPredict要求)
expr_drug <- as.matrix(expr_tpm_symbol[, common_samples])
expr_log_drug <- log2(expr_drug + 1)

message(sprintf("Samples: %d (High: %d, Low: %d)",
                length(common_samples),
                sum(risk_groups == "High"), sum(risk_groups == "Low")))

# =============================================================================
# 2. oncoPredict -- GDSC2 IC50预测
# =============================================================================
message("\n=== Section 2: oncoPredict GDSC2 IC50 Prediction ===")

# 重点关注的药物 (与UPR/胶质瘤相关)
focus_drugs <- c(
  "GSK2606414",   # PERK抑制剂 (直接靶向UPR)
  "Bortezomib",   # 蛋白酶体抑制剂 (加剧ER应激)
  "Temozolomide",  # 胶质瘤标准化疗
  "Rapamycin",     # mTOR抑制剂 (与UPR交叉调控)
  "Erlotinib",     # EGFR-TKI (胶质瘤靶向)
  "17-AAG"         # HSP90抑制剂 (蛋白稳态)
)

# oncoPredict需要GDSC2训练数据
# 数据路径配置
gdsc_dir <- file.path(DATA_RAW, "GDSC")
dir.create(gdsc_dir, recursive = TRUE, showWarnings = FALSE)

# GDSC2训练数据文件
gdsc2_expr_file <- file.path(gdsc_dir, "GDSC2_Expr.rds")
gdsc2_res_file  <- file.path(gdsc_dir, "GDSC2_Res.rds")

oncopredict_success <- FALSE
drug_ic50_results <- list()

tryCatch({
  library(oncoPredict)

  # 方案1: 使用oncoPredict自带数据
  if (exists("GDSC2_Expr") && exists("GDSC2_Res")) {
    message("Using built-in GDSC2 data from oncoPredict...")
    train_expr_gdsc <- GDSC2_Expr
    train_res_gdsc  <- GDSC2_Res

  } else if (file.exists(gdsc2_expr_file) && file.exists(gdsc2_res_file)) {
    # 方案2: 从本地文件加载
    message("Loading GDSC2 training data from local files...")
    train_expr_gdsc <- readRDS(gdsc2_expr_file)
    train_res_gdsc  <- readRDS(gdsc2_res_file)

  } else {
    # 方案3: 尝试从包数据目录获取
    pkg_data <- system.file("extdata", package = "oncoPredict")
    if (nchar(pkg_data) > 0) {
      rds_files <- list.files(pkg_data, pattern = "GDSC2.*\\.rds$", full.names = TRUE)
      if (length(rds_files) >= 2) {
        expr_f <- rds_files[grep("Expr", rds_files)]
        res_f  <- rds_files[grep("Res", rds_files)]
        if (length(expr_f) > 0 && length(res_f) > 0) {
          train_expr_gdsc <- readRDS(expr_f[1])
          train_res_gdsc  <- readRDS(res_f[1])
        }
      }
    }
  }

  if (exists("train_expr_gdsc") && exists("train_res_gdsc")) {
    message("Running calcPhenotype for GDSC2...")
    message(sprintf("  Training data: %d genes x %d cell lines",
                    nrow(train_expr_gdsc), ncol(train_expr_gdsc)))
    message(sprintf("  Test data: %d genes x %d samples",
                    nrow(expr_drug), ncol(expr_drug)))

    # 运行IC50预测
    oncopredict_outdir <- file.path(RES_DIR, "oncoPredict_output")
    dir.create(oncopredict_outdir, recursive = TRUE, showWarnings = FALSE)

    old_wd <- getwd()
    setwd(oncopredict_outdir)

    calcPhenotype(
      trainingExprData  = train_expr_gdsc,
      trainingPtype     = train_res_gdsc,
      testExprData      = expr_drug,
      batchCorrect      = "eb",
      powerTransformPhenotype = TRUE,
      removeLowVaryingGenes   = 0.2,
      minNumSamples     = 10,
      selection         = 1,
      printOutput       = TRUE,
      removeLowVaringGenesFrom = "rawData"
    )

    setwd(old_wd)

    # 读取结果
    pred_files <- list.files(oncopredict_outdir,
                             pattern = "DrugPredictions.*\\.csv$",
                             full.names = TRUE)
    if (length(pred_files) > 0) {
      drug_pred_all <- read.csv(pred_files[1], row.names = 1, check.names = FALSE)
      message(sprintf("oncoPredict completed: %d drugs x %d samples",
                      ncol(drug_pred_all), nrow(drug_pred_all)))

      # 提取各药物结果
      all_drug_cols <- colnames(drug_pred_all)
      for (drug in focus_drugs) {
        # 模糊匹配列名（GDSC药物名可能含编号后缀）
        matched_cols <- grep(drug, all_drug_cols, ignore.case = TRUE, value = TRUE)
        if (length(matched_cols) > 0) {
          drug_ic50_results[[drug]] <- setNames(
            as.numeric(drug_pred_all[, matched_cols[1]]),
            rownames(drug_pred_all)
          )
          message(sprintf("  Found %s: column '%s'", drug, matched_cols[1]))
        } else {
          message(sprintf("  %s not found in GDSC2 predictions", drug))
        }
      }

      oncopredict_success <- TRUE
    }
  } else {
    message("GDSC2 training data not available for oncoPredict")
  }
}, error = function(e) {
  message("oncoPredict analysis failed: ", e$message)
})

# =============================================================================
# 3. 备选方案: pRRophetic
# =============================================================================
if (!oncopredict_success) {
  message("\n=== Fallback: pRRophetic analysis ===")

  tryCatch({
    if (requireNamespace("pRRophetic", quietly = TRUE)) {
      library(pRRophetic)

      # pRRophetic使用的药物名可能不同
      prroph_drugs <- c(
        "Bortezomib",
        "Temozolomide",
        "Rapamycin",
        "Erlotinib",
        "17-AAG"
      )

      for (drug in prroph_drugs) {
        message(sprintf("  Predicting IC50 for %s (pRRophetic)...", drug))
        tryCatch({
          ic50 <- pRRopheticPredict(
            testMatrix   = expr_drug,
            drug         = drug,
            tissueType   = "all",
            batchCorrect = "eb",
            selection    = 1,
            dataset      = "cgp2016"
          )
          drug_ic50_results[[drug]] <- ic50
          message(sprintf("    %s: predicted for %d samples", drug, length(ic50)))
        }, error = function(e) {
          message(sprintf("    %s failed: %s", drug, e$message))
        })
      }

      if (length(drug_ic50_results) > 0) oncopredict_success <- TRUE
    } else {
      message("pRRophetic not installed")
    }
  }, error = function(e) {
    message("pRRophetic failed: ", e$message)
  })
}

# =============================================================================
# 4. IC50高低风险组比较（Wilcoxon test）
# =============================================================================
if (length(drug_ic50_results) > 0) {
  message("\n=== Section 4: IC50 Comparison by Risk Group ===")

  # 构建长格式数据框
  ic50_list <- lapply(names(drug_ic50_results), function(drug) {
    vals <- drug_ic50_results[[drug]]
    samps <- intersect(names(vals), common_samples)
    if (length(samps) == 0) return(NULL)
    data.frame(
      barcode    = samps,
      IC50       = vals[samps],
      Drug       = drug,
      risk_group = factor(risk_groups[samps], levels = c("Low", "High")),
      stringsAsFactors = FALSE
    )
  })
  ic50_df <- do.call(rbind, ic50_list[!sapply(ic50_list, is.null)])

  # Wilcoxon检验统计
  ic50_stats <- ic50_df %>%
    dplyr::group_by(Drug) %>%
    dplyr::summarise(
      median_high = median(IC50[risk_group == "High"], na.rm = TRUE),
      median_low  = median(IC50[risk_group == "Low"], na.rm = TRUE),
      p_value     = wilcox.test(IC50 ~ risk_group)$p.value,
      .groups     = "drop"
    ) %>%
    dplyr::mutate(
      padj     = p.adjust(p_value, method = "BH"),
      sig_mark = ifelse(padj < 0.001, "***",
                        ifelse(padj < 0.01, "**",
                               ifelse(padj < 0.05, "*", "ns")))
    )

  write.csv(ic50_stats, file.path(RES_DIR, "drug_ic50_statistics.csv"), row.names = FALSE)
  write.csv(ic50_df, file.path(RES_DIR, "drug_ic50_predictions.csv"), row.names = FALSE)

  message("IC50 comparison results:")
  print(ic50_stats)

  # --- IC50箱线图 ---
  n_drugs <- length(unique(ic50_df$Drug))
  ncols <- min(3, n_drugs)
  nrows <- ceiling(n_drugs / ncols)

  p_ic50 <- ggplot(ic50_df, aes(x = risk_group, y = IC50, fill = risk_group)) +
    geom_boxplot(outlier.size = 0.3, width = 0.6) +
    facet_wrap(~Drug, scales = "free_y", ncol = ncols) +
    stat_compare_means(method = "wilcox.test", label = "p.format", size = 3.5) +
    scale_fill_manual(values = COLORS_RISK) +
    labs(x = "Risk Group", y = "Predicted IC50 (log2)",
         title = "Predicted Drug IC50 by Risk Group (GDSC2)") +
    THEME_PUBLICATION +
    theme(legend.position = "none",
          strip.text = element_text(size = 10, face = "bold"))

  ggsave(file.path(fig_dir, "Fig7D_drug_ic50.pdf"), p_ic50,
         width = ncols * 4.5, height = nrows * 4)

  # --- 单独绘制每个重点药物的精细箱线图 ---
  for (drug in names(drug_ic50_results)) {
    drug_sub <- ic50_df[ic50_df$Drug == drug, ]
    if (nrow(drug_sub) < 10) next

    p_single <- ggplot(drug_sub, aes(x = risk_group, y = IC50, fill = risk_group)) +
      geom_boxplot(outlier.shape = NA, width = 0.5) +
      geom_jitter(width = 0.12, size = 0.5, alpha = 0.3) +
      stat_compare_means(method = "wilcox.test", label = "p.format", size = 4.5) +
      scale_fill_manual(values = COLORS_RISK) +
      labs(x = "Risk Group", y = "Predicted IC50",
           title = drug) +
      THEME_PUBLICATION +
      theme(legend.position = "none")

    ggsave(file.path(fig_dir, paste0("Fig7D_ic50_", gsub("[^a-zA-Z0-9]", "_", drug), ".pdf")),
           p_single, width = 4.5, height = 5)
  }

} else {
  message("\nWARNING: No drug IC50 predictions available.")
  message("Please ensure oncoPredict is installed with GDSC2 training data.")
  message("Install: BiocManager::install('oncoPredict')")
  message("Training data: https://osf.io/c6tfx/ (GDSC2_Expr.rds, GDSC2_Res.rds)")
  message("Place in: ", gdsc_dir)
}

# =============================================================================
# 5. 药物靶点通路ssGSEA活性比较
# =============================================================================
message("\n=== Section 5: Drug Target Pathway Activity ===")

# 药物靶点通路基因集
drug_pathways <- list(
  PERK_UPR_signaling = c("EIF2AK3", "ATF4", "DDIT3", "TRIB3",
                         "PPP1R15A", "ASNS", "CHAC1", "SESN2"),
  IRE1_XBP1_signaling = c("ERN1", "XBP1", "DNAJB9", "EDEM1",
                           "HERPUD1", "SYVN1", "DERL1"),
  Proteasome_activity = c("PSMA1", "PSMA2", "PSMA3", "PSMA4",
                          "PSMB1", "PSMB2", "PSMB5", "PSMB8", "PSMB9"),
  mTOR_signaling = c("MTOR", "RPTOR", "RICTOR", "RPS6KB1",
                     "EIF4EBP1", "AKT1", "TSC1", "TSC2"),
  EGFR_signaling = c("EGFR", "ERBB2", "GRB2", "SOS1",
                     "KRAS", "BRAF", "MAP2K1", "MAPK1"),
  HSP90_chaperone = c("HSP90AA1", "HSP90AB1", "HSP90B1",
                      "HSPA5", "TRAP1", "HSPH1"),
  DNA_damage_repair = c("MGMT", "MLH1", "MSH2", "MSH6",
                        "BRCA1", "BRCA2", "ATM", "ATR"),
  Autophagy = c("BECN1", "ATG5", "ATG7", "ATG12",
                "MAP1LC3B", "SQSTM1", "ULK1", "ATG16L1")
)

# 过滤基因
drug_pathways_valid <- lapply(drug_pathways, function(gs) {
  gs[gs %in% rownames(expr_log_drug)]
})

for (nm in names(drug_pathways_valid)) {
  message(sprintf("  %s: %d/%d genes", nm,
                  length(drug_pathways_valid[[nm]]), length(drug_pathways[[nm]])))
}

# ssGSEA
drug_pw_params <- ssgseaParam(
  exprData = expr_log_drug[, common_samples],
  geneSets = drug_pathways_valid
)
drug_pw_scores <- gsva(drug_pw_params)

drug_pw_df <- as.data.frame(t(drug_pw_scores))
drug_pw_df$barcode <- rownames(drug_pw_df)
drug_pw_df$risk_group <- factor(risk_groups[rownames(drug_pw_df)], levels = c("Low", "High"))
drug_pw_df <- drug_pw_df[!is.na(drug_pw_df$risk_group), ]

drug_pw_long <- drug_pw_df %>%
  tidyr::pivot_longer(cols = -c(barcode, risk_group),
                      names_to = "Pathway", values_to = "Score") %>%
  dplyr::mutate(Pathway = gsub("_", " ", Pathway))

# Wilcoxon统计
pw_stats <- drug_pw_long %>%
  dplyr::group_by(Pathway) %>%
  dplyr::summarise(
    median_high = median(Score[risk_group == "High"]),
    median_low  = median(Score[risk_group == "Low"]),
    p_value = wilcox.test(Score ~ risk_group)$p.value,
    .groups = "drop"
  ) %>%
  dplyr::mutate(padj = p.adjust(p_value, method = "BH"))

write.csv(pw_stats, file.path(RES_DIR, "drug_pathway_statistics.csv"), row.names = FALSE)

# 通路活性箱线图
p_drug_pw <- ggplot(drug_pw_long, aes(x = risk_group, y = Score, fill = risk_group)) +
  geom_boxplot(outlier.size = 0.3, width = 0.6) +
  facet_wrap(~Pathway, scales = "free_y", ncol = 4) +
  stat_compare_means(method = "wilcox.test", label = "p.format", size = 3) +
  scale_fill_manual(values = COLORS_RISK) +
  labs(x = "Risk Group", y = "ssGSEA Enrichment Score",
       title = "Drug Target Pathway Activity by Risk Group") +
  THEME_PUBLICATION +
  theme(legend.position = "bottom",
        strip.text = element_text(size = 9, face = "bold"))

ggsave(file.path(fig_dir, "Fig7D_drug_pathways.pdf"), p_drug_pw, width = 14, height = 8)

# =============================================================================
# 6. 综合热图: IC50 + 通路活性
# =============================================================================
message("\n=== Section 6: Integrated Drug Heatmap ===")

if (length(drug_ic50_results) > 0) {
  # IC50热图矩阵
  ic50_mat <- do.call(cbind, lapply(names(drug_ic50_results), function(drug) {
    vals <- drug_ic50_results[[drug]]
    vals[common_samples]
  }))
  colnames(ic50_mat) <- names(drug_ic50_results)
  ic50_mat <- ic50_mat[!apply(ic50_mat, 1, function(x) all(is.na(x))), ]

  # Z-score标准化
  ic50_scaled <- scale(ic50_mat)

  # 标注信息
  ha_col <- HeatmapAnnotation(
    Risk = risk_groups[rownames(ic50_scaled)],
    col = list(Risk = COLORS_RISK),
    annotation_name_side = "left"
  )

  # 按风险分组排序
  sample_order <- order(risk_groups[rownames(ic50_scaled)])

  pdf(file.path(fig_dir, "Fig7D_ic50_heatmap.pdf"), width = 10, height = 8)
  draw(Heatmap(
    t(ic50_scaled[sample_order, ]),
    name = "IC50\n(Z-score)",
    col = colorRamp2(c(-2, 0, 2), c("#4DBBD5", "white", "#E64B35")),
    top_annotation = ha_col,
    cluster_columns = FALSE,
    cluster_rows = TRUE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 10),
    column_title = "Predicted Drug IC50 (GDSC2)",
    column_title_gp = gpar(fontsize = 12, fontface = "bold")
  ))
  dev.off()
}

# =============================================================================
# 7. 保存结果
# =============================================================================
message("\n=== Saving results ===")

save_objs <- c("drug_pw_scores", "drug_pathways")
if (length(drug_ic50_results) > 0) {
  save_objs <- c(save_objs, "drug_ic50_results")
}
if (exists("ic50_stats")) save_objs <- c(save_objs, "ic50_stats")

save(list = save_objs,
     file = file.path(DATA_PROC, "drug_sensitivity_results.RData"))

message("\n=== Drug sensitivity analysis completed ===")
message("Figures saved to: ", fig_dir)
message("  Fig7D: Drug IC50 boxplots + pathway activity")
message("\nNext: Run 06_clinical_translation/03_cmap_analysis.R")
