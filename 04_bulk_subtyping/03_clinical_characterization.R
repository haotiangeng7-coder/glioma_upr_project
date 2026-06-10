###############################################################################
# 03_clinical_characterization.R
# 各UPR亚型的临床特征刻画
#
# 输出:
#   - figures/Fig3_clinical_heatmap.pdf      (综合热图: UPR基因+临床注释条)
#   - figures/Fig3_age_by_subtype.pdf
#   - figures/Fig3_grade_distribution.pdf
#   - figures/Fig3_idh_distribution.pdf
#   - figures/Fig3_clinical_features_combined.pdf
#   - figures/Fig3_verhaak_confounding.pdf    (UPR亚型 x Verhaak交叉分析)
#   - results/clinical_association_tests.csv
#   - results/verhaak_confounding_analysis.csv
#   - data/processed/clinical_characterization_results.RData
#
# 评审修订:
#   - GBM分子亚型混杂评估: IDH-WT中UPR亚型 x Verhaak亚型交叉列联表 + Cramer's V
#   - 增加1p/19q共缺失状态分析
#   - 综合热图加入Verhaak亚型和1p19q注释条
###############################################################################

source("00_setup/config.R")

library(ggplot2)
library(ggpubr)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)
library(survival)
library(survminer)
# rcompanion不兼容R 4.3，手动实现cramerV
cramerV <- function(x, ...) {
  chi2 <- suppressWarnings(chisq.test(x, simulate.p.value = TRUE, B = 2000))
  n <- sum(x)
  k <- min(nrow(x), ncol(x))
  v <- sqrt(chi2$statistic / (n * (k - 1)))
  names(v) <- "Cramer.V"
  return(v)
}

set.seed(SEED)

# --- 加载数据 ---
load(file.path(DATA_PROC, "upr_gene_sets.RData"))
load(file.path(DATA_PROC, "tcga_glioma_expression.RData"))
load(file.path(DATA_PROC, "upr_landscape_results.RData"))
load(file.path(DATA_PROC, "consensus_clustering_results.RData"))

# =============================================================================
# 0. 数据准备
# =============================================================================
message("=== Preparing data for clinical characterization ===")

# 合并聚类结果和临床数据
heatmap_data <- merge(cluster_df, clin_matched, by = "barcode")
message(sprintf("  Samples with subtype + clinical: %d", nrow(heatmap_data)))

# 定义颜色方案
grade_colors <- c("G2" = "#FFF5EB", "G3" = "#FDAE6B", "G4" = "#D94701")
idh_colors   <- c("Mutant" = "#00A087", "WT" = "#E64B35")
mgmt_colors  <- c("Methylated" = "#3C5488", "Unmethylated" = "#F39B7F")
verhaak_colors <- c(
  "Classical"     = "#E64B35", "CL" = "#E64B35",
  "Mesenchymal"   = "#4DBBD5", "ME" = "#4DBBD5",
  "Neural"        = "#00A087", "NE" = "#00A087",
  "Proneural"     = "#3C5488", "PN" = "#3C5488"
)
codel_colors <- c("codel" = "#7570B3", "non-codel" = "#D95F02")

# =============================================================================
# 1. 综合热图 (Figure 3 核心panel)
# =============================================================================
message("=== Building comprehensive clinical heatmap ===")

# 准备注释数据
anno_cols <- c("barcode", "UPR_subtype")
col_spec <- list(UPR_Subtype = COLORS_SUBTYPE)

# 动态添加可用的临床注释
if ("Grade" %in% colnames(heatmap_data)) {
  anno_cols <- c(anno_cols, "Grade")
  col_spec$Grade <- grade_colors
}
if ("IDH_status" %in% colnames(heatmap_data)) {
  anno_cols <- c(anno_cols, "IDH_status")
  col_spec$IDH <- idh_colors
}
if ("MGMT_status" %in% colnames(heatmap_data)) {
  anno_cols <- c(anno_cols, "MGMT_status")
  col_spec$MGMT <- mgmt_colors
}
if ("Subtype" %in% colnames(heatmap_data)) {
  anno_cols <- c(anno_cols, "Subtype")
  col_spec$Verhaak <- verhaak_colors
}

# 检查1p/19q共缺失列
codel_col <- NULL
for (candidate in c("paper_IDH.codel.subtype", "IDH.codel.subtype", "codel_status")) {
  if (candidate %in% colnames(heatmap_data)) {
    codel_col <- candidate
    break
  }
}
if (!is.null(codel_col)) {
  # 简化1p/19q状态：提取codel vs non-codel
  heatmap_data$codel_1p19q <- ifelse(
    grepl("codel", heatmap_data[[codel_col]], ignore.case = TRUE) &
      !grepl("non-codel|non_codel", heatmap_data[[codel_col]], ignore.case = TRUE),
    "codel", "non-codel"
  )
  heatmap_data$codel_1p19q[is.na(heatmap_data[[codel_col]])] <- NA
  anno_cols <- c(anno_cols, "codel_1p19q")
  col_spec$`1p/19q` <- codel_colors
}

# 排序样本（按UPR亚型分块）
sample_order <- heatmap_data %>%
  dplyr::arrange(UPR_subtype) %>%
  dplyr::pull(barcode)

# 准备表达矩阵（聚类基因 x 排序样本）
expr_upr_heatmap <- as.matrix(log2(expr_tpm_symbol[cluster_genes, sample_order] + 1))
expr_upr_heatmap <- t(scale(t(expr_upr_heatmap)))

# 限制Z-score范围
expr_upr_heatmap[expr_upr_heatmap > 2]  <- 2
expr_upr_heatmap[expr_upr_heatmap < -2] <- -2

# 构建行注释（基因通路标签）
gene_pathway <- dplyr::case_when(
  cluster_genes %in% IRE1_XBP1_genes ~ "IRE1-XBP1",
  cluster_genes %in% PERK_ATF4_genes ~ "PERK-ATF4",
  cluster_genes %in% ATF6_genes      ~ "ATF6",
  TRUE                               ~ "Other/Shared"
)
names(gene_pathway) <- cluster_genes

row_anno <- rowAnnotation(
  Pathway = gene_pathway,
  col = list(Pathway = c("IRE1-XBP1" = "#E64B35", "PERK-ATF4" = "#4DBBD5",
                          "ATF6" = "#00A087", "Other/Shared" = "grey80")),
  show_annotation_name = TRUE
)

# 构建列注释（动态）
anno_list <- list()
hm_row <- heatmap_data[match(sample_order, heatmap_data$barcode), ]

anno_list[["UPR_Subtype"]] <- hm_row$UPR_subtype

if ("Grade" %in% colnames(hm_row)) {
  anno_list[["Grade"]] <- factor(hm_row$Grade)
}
if ("IDH_status" %in% colnames(hm_row)) {
  anno_list[["IDH"]] <- factor(hm_row$IDH_status)
}
if ("MGMT_status" %in% colnames(hm_row)) {
  anno_list[["MGMT"]] <- factor(hm_row$MGMT_status)
}
if ("Subtype" %in% colnames(hm_row)) {
  anno_list[["Verhaak"]] <- factor(hm_row$Subtype)
}
if ("codel_1p19q" %in% colnames(hm_row)) {
  anno_list[["1p/19q"]] <- factor(hm_row$codel_1p19q)
}

top_anno <- HeatmapAnnotation(
  df  = as.data.frame(anno_list),
  col = col_spec,
  annotation_name_side = "left",
  na_col = "grey90",
  show_legend = TRUE
)

# 绘制热图
pdf(file.path(FIG_DIR, "Fig3_clinical_heatmap.pdf"), width = 16, height = 12)
ht <- Heatmap(
  expr_upr_heatmap,
  name = "Z-score",
  col  = colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#B2182B")),
  top_annotation  = top_anno,
  left_annotation = row_anno,
  cluster_columns = FALSE,
  cluster_rows    = TRUE,
  clustering_method_rows = "ward.D2",
  show_column_names = FALSE,
  row_names_gp = gpar(fontsize = 6),
  column_split = factor(hm_row$UPR_subtype,
                         levels = c("UPR-high-risk", "UPR-intermediate", "UPR-favorable")),
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  column_title = "UPR Molecular Subtypes - Clinical Characterization",
  heatmap_legend_param = list(
    title = "Expression\n(Z-score)",
    at = c(-2, -1, 0, 1, 2),
    legend_height = unit(4, "cm")
  )
)
draw(ht, merge_legend = TRUE)
dev.off()

# =============================================================================
# 2. 各亚型临床特征统计检验
# =============================================================================
message("=== Clinical feature association tests ===")

clinical_stats <- list()
clinical_tests_df <- data.frame(
  Variable     = character(),
  Test         = character(),
  Statistic    = numeric(),
  P_value      = numeric(),
  N_available  = integer(),
  stringsAsFactors = FALSE
)

# --- 2.1 Age ---
if ("age_at_index" %in% colnames(heatmap_data)) {
  heatmap_data$Age <- as.numeric(heatmap_data$age_at_index)
  age_data <- heatmap_data %>% dplyr::filter(!is.na(Age))

  kw_age <- kruskal.test(Age ~ UPR_subtype, data = age_data)
  clinical_stats$age_pvalue <- kw_age$p.value
  clinical_tests_df <- rbind(clinical_tests_df, data.frame(
    Variable = "Age", Test = "Kruskal-Wallis",
    Statistic = kw_age$statistic, P_value = kw_age$p.value,
    N_available = nrow(age_data), stringsAsFactors = FALSE
  ))

  p_age <- ggplot(age_data, aes(x = UPR_subtype, y = Age, fill = UPR_subtype)) +
    geom_boxplot(outlier.size = 0.5, width = 0.6) +
    geom_jitter(width = 0.15, size = 0.3, alpha = 0.2) +
    stat_compare_means(method = "kruskal.test", label = "p.format") +
    scale_fill_manual(values = COLORS_SUBTYPE) +
    THEME_PUBLICATION +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1)) +
    labs(x = "UPR Subtype", y = "Age at Diagnosis",
         title = "Age Distribution by UPR Subtype")

  ggsave(file.path(FIG_DIR, "Fig3_age_by_subtype.pdf"), p_age,
         width = 6, height = 5, useDingbats = FALSE)

  message(sprintf("  Age: Kruskal-Wallis p = %.4e", kw_age$p.value))
}

# --- 2.2 Gender ---
if ("gender" %in% colnames(heatmap_data)) {
  gender_data <- heatmap_data %>% dplyr::filter(!is.na(gender))
  gender_table <- table(gender_data$UPR_subtype, gender_data$gender)

  if (all(dim(gender_table) >= 2)) {
    gender_test <- chisq.test(gender_table)
    clinical_stats$gender_pvalue <- gender_test$p.value
    clinical_tests_df <- rbind(clinical_tests_df, data.frame(
      Variable = "Gender", Test = "Chi-squared",
      Statistic = gender_test$statistic, P_value = gender_test$p.value,
      N_available = nrow(gender_data), stringsAsFactors = FALSE
    ))
    message(sprintf("  Gender: chi-squared p = %.4f", gender_test$p.value))
  }
}

# --- 2.3 Grade ---
if ("Grade" %in% colnames(heatmap_data)) {
  grade_data <- heatmap_data %>% dplyr::filter(!is.na(Grade))
  grade_table <- table(grade_data$UPR_subtype, grade_data$Grade)

  if (all(dim(grade_table) >= 2)) {
    grade_test <- chisq.test(grade_table)
    clinical_stats$grade_pvalue <- grade_test$p.value
    clinical_tests_df <- rbind(clinical_tests_df, data.frame(
      Variable = "WHO Grade", Test = "Chi-squared",
      Statistic = grade_test$statistic, P_value = grade_test$p.value,
      N_available = nrow(grade_data), stringsAsFactors = FALSE
    ))

    # 比例条形图
    grade_prop <- as.data.frame(prop.table(grade_table, margin = 1) * 100)
    colnames(grade_prop) <- c("Subtype", "Grade", "Percentage")

    p_grade <- ggplot(grade_prop, aes(x = Subtype, y = Percentage, fill = Grade)) +
      geom_bar(stat = "identity", position = "stack", width = 0.7) +
      scale_fill_manual(values = grade_colors) +
      labs(y = "Percentage (%)", x = "UPR Subtype",
           title = sprintf("WHO Grade Distribution (p = %.2e)", grade_test$p.value)) +
      THEME_PUBLICATION +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))

    ggsave(file.path(FIG_DIR, "Fig3_grade_distribution.pdf"), p_grade,
           width = 7, height = 5, useDingbats = FALSE)
    message(sprintf("  Grade: chi-squared p = %.4e", grade_test$p.value))
  }
}

# --- 2.4 IDH ---
if ("IDH_status" %in% colnames(heatmap_data)) {
  idh_data <- heatmap_data %>% dplyr::filter(!is.na(IDH_status))
  idh_table <- table(idh_data$UPR_subtype, idh_data$IDH_status)

  if (all(dim(idh_table) >= 2)) {
    idh_test <- chisq.test(idh_table)
    clinical_stats$idh_pvalue <- idh_test$p.value
    clinical_tests_df <- rbind(clinical_tests_df, data.frame(
      Variable = "IDH status", Test = "Chi-squared",
      Statistic = idh_test$statistic, P_value = idh_test$p.value,
      N_available = nrow(idh_data), stringsAsFactors = FALSE
    ))

    idh_prop <- as.data.frame(prop.table(idh_table, margin = 1) * 100)
    colnames(idh_prop) <- c("Subtype", "IDH", "Percentage")

    p_idh <- ggplot(idh_prop, aes(x = Subtype, y = Percentage, fill = IDH)) +
      geom_bar(stat = "identity", position = "stack", width = 0.7) +
      scale_fill_manual(values = idh_colors) +
      labs(y = "Percentage (%)", x = "UPR Subtype",
           title = sprintf("IDH Status Distribution (p = %.2e)", idh_test$p.value)) +
      THEME_PUBLICATION +
      theme(axis.text.x = element_text(angle = 30, hjust = 1))

    ggsave(file.path(FIG_DIR, "Fig3_idh_distribution.pdf"), p_idh,
           width = 7, height = 5, useDingbats = FALSE)
    message(sprintf("  IDH: chi-squared p = %.4e", idh_test$p.value))
  }
}

# --- 2.5 MGMT ---
if ("MGMT_status" %in% colnames(heatmap_data)) {
  mgmt_data <- heatmap_data %>% dplyr::filter(!is.na(MGMT_status))
  if (nrow(mgmt_data) > 10) {
    mgmt_table <- table(mgmt_data$UPR_subtype, mgmt_data$MGMT_status)
    if (all(dim(mgmt_table) >= 2)) {
      mgmt_test <- chisq.test(mgmt_table)
      clinical_stats$mgmt_pvalue <- mgmt_test$p.value
      clinical_tests_df <- rbind(clinical_tests_df, data.frame(
        Variable = "MGMT methylation", Test = "Chi-squared",
        Statistic = mgmt_test$statistic, P_value = mgmt_test$p.value,
        N_available = nrow(mgmt_data), stringsAsFactors = FALSE
      ))
      message(sprintf("  MGMT: chi-squared p = %.4e", mgmt_test$p.value))
    }
  }
}

# --- 2.6 1p/19q ---
if ("codel_1p19q" %in% colnames(heatmap_data)) {
  codel_data <- heatmap_data %>% dplyr::filter(!is.na(codel_1p19q))
  if (nrow(codel_data) > 10) {
    codel_table <- table(codel_data$UPR_subtype, codel_data$codel_1p19q)
    if (all(dim(codel_table) >= 2)) {
      codel_test <- chisq.test(codel_table)
      clinical_stats$codel_pvalue <- codel_test$p.value
      clinical_tests_df <- rbind(clinical_tests_df, data.frame(
        Variable = "1p/19q codeletion", Test = "Chi-squared",
        Statistic = codel_test$statistic, P_value = codel_test$p.value,
        N_available = nrow(codel_data), stringsAsFactors = FALSE
      ))
      message(sprintf("  1p/19q: chi-squared p = %.4e", codel_test$p.value))
    }
  }
}

# --- 2.7 Verhaak subtype ---
if ("Subtype" %in% colnames(heatmap_data)) {
  verhaak_data <- heatmap_data %>% dplyr::filter(!is.na(Subtype))
  if (nrow(verhaak_data) > 10) {
    verhaak_table <- table(verhaak_data$UPR_subtype, verhaak_data$Subtype)
    if (all(dim(verhaak_table) >= 2)) {
      verhaak_test <- chisq.test(verhaak_table)
      clinical_stats$verhaak_pvalue <- verhaak_test$p.value
      clinical_tests_df <- rbind(clinical_tests_df, data.frame(
        Variable = "Verhaak subtype", Test = "Chi-squared",
        Statistic = verhaak_test$statistic, P_value = verhaak_test$p.value,
        N_available = nrow(verhaak_data), stringsAsFactors = FALSE
      ))
      message(sprintf("  Verhaak: chi-squared p = %.4e", verhaak_test$p.value))
    }
  }
}

# FDR校正
if (nrow(clinical_tests_df) > 1) {
  clinical_tests_df$P_adjusted <- p.adjust(clinical_tests_df$P_value, method = "BH")
} else {
  clinical_tests_df$P_adjusted <- clinical_tests_df$P_value
}

write.csv(clinical_tests_df, file.path(RES_DIR, "clinical_association_tests.csv"),
          row.names = FALSE)

# =============================================================================
# 3. GBM分子亚型混杂评估（评审修订: Cramer's V）
# =============================================================================
message("=== Verhaak subtype confounding assessment (reviewer revision) ===")
message("  Evaluating whether UPR subtypes simply recapitulate Verhaak classification")

confounding_results <- data.frame(
  Analysis       = character(),
  N_samples      = integer(),
  Chi2_statistic = numeric(),
  Chi2_p         = numeric(),
  Cramers_V      = numeric(),
  V_interpretation = character(),
  stringsAsFactors = FALSE
)

if ("IDH_status" %in% colnames(heatmap_data) && "Subtype" %in% colnames(heatmap_data)) {

  # --- 3.1 IDH-WT亚组: UPR亚型 x Verhaak亚型 ---
  idhwt_data <- heatmap_data %>%
    dplyr::filter(IDH_status == "WT", !is.na(Subtype), !is.na(UPR_subtype))

  message(sprintf("  IDH-WT samples with both UPR and Verhaak subtypes: %d", nrow(idhwt_data)))

  if (nrow(idhwt_data) >= 20) {
    # 交叉列联表
    cross_table <- table(idhwt_data$UPR_subtype, idhwt_data$Subtype)
    message("\n  Cross-tabulation (IDH-WT): UPR subtype x Verhaak subtype")
    print(cross_table)

    # 比例表
    cross_prop <- prop.table(cross_table, margin = 1) * 100
    message("\n  Row percentages:")
    print(round(cross_prop, 1))

    # Chi-squared test
    chi2_idhwt <- chisq.test(cross_table)

    # Cramer's V (rcompanion包)
    cramers_v_idhwt <- cramerV(cross_table)
    v_value <- cramers_v_idhwt[1]  # Cramer's V value

    # 效应量解读 (Cohen's guidelines for Cramer's V)
    # df* = min(nrow, ncol) - 1
    df_star <- min(nrow(cross_table), ncol(cross_table)) - 1
    v_interp <- dplyr::case_when(
      v_value < 0.1  ~ "negligible",
      v_value < 0.3  ~ "small",
      v_value < 0.5  ~ "medium",
      TRUE           ~ "large"
    )

    message(sprintf("\n  Chi-squared: X2=%.2f, df=%d, p=%.4e",
                    chi2_idhwt$statistic, chi2_idhwt$parameter, chi2_idhwt$p.value))
    message(sprintf("  Cramer's V = %.4f (%s effect)", v_value, v_interp))

    if (v_value >= 0.5) {
      message("  WARNING: Large Cramer's V suggests substantial overlap with Verhaak classification.")
      message("  Consider discussing this overlap in the manuscript.")
    } else {
      message("  Cramer's V < 0.5: UPR subtypes provide non-redundant information vs Verhaak.")
    }

    confounding_results <- rbind(confounding_results, data.frame(
      Analysis       = "IDH-WT: UPR x Verhaak",
      N_samples      = nrow(idhwt_data),
      Chi2_statistic = chi2_idhwt$statistic,
      Chi2_p         = chi2_idhwt$p.value,
      Cramers_V      = v_value,
      V_interpretation = v_interp,
      stringsAsFactors = FALSE
    ))

    # --- 可视化: 热力图 + 气泡图 ---
    # 气泡图显示交叉比例
    cross_df <- as.data.frame(cross_table)
    colnames(cross_df) <- c("UPR_subtype", "Verhaak_subtype", "Count")
    cross_df <- cross_df %>%
      dplyr::group_by(UPR_subtype) %>%
      dplyr::mutate(
        Total = sum(Count),
        Percentage = Count / Total * 100
      ) %>%
      dplyr::ungroup()

    p_confound <- ggplot(cross_df,
                          aes(x = Verhaak_subtype, y = UPR_subtype)) +
      geom_point(aes(size = Percentage, color = Percentage)) +
      geom_text(aes(label = sprintf("%d\n(%.0f%%)", Count, Percentage)),
                size = 3, vjust = -0.2) +
      scale_size_continuous(range = c(2, 15), name = "% within\nUPR subtype") +
      scale_color_gradient(low = "#DEEBF7", high = "#08519C",
                            name = "% within\nUPR subtype") +
      labs(
        x = "Verhaak Transcriptome Subtype",
        y = "UPR Subtype",
        title = "UPR vs Verhaak Subtype (IDH-WT)",
        subtitle = sprintf("Cramer's V = %.3f (%s); Chi2 p = %.2e",
                            v_value, v_interp, chi2_idhwt$p.value)
      ) +
      THEME_PUBLICATION +
      theme(
        plot.subtitle = element_text(size = 10, color = "grey30"),
        axis.text = element_text(size = 10)
      )

    ggsave(file.path(FIG_DIR, "Fig3_verhaak_confounding.pdf"), p_confound,
           width = 9, height = 6, useDingbats = FALSE)
  }

  # --- 3.2 全队列: UPR亚型 x Verhaak亚型 ---
  all_verhaak <- heatmap_data %>%
    dplyr::filter(!is.na(Subtype), !is.na(UPR_subtype))

  if (nrow(all_verhaak) >= 20) {
    cross_all <- table(all_verhaak$UPR_subtype, all_verhaak$Subtype)
    chi2_all <- chisq.test(cross_all)
    cramers_v_all <- cramerV(cross_all)
    v_all <- cramers_v_all[1]
    v_interp_all <- dplyr::case_when(
      v_all < 0.1  ~ "negligible",
      v_all < 0.3  ~ "small",
      v_all < 0.5  ~ "medium",
      TRUE         ~ "large"
    )

    confounding_results <- rbind(confounding_results, data.frame(
      Analysis       = "Full cohort: UPR x Verhaak",
      N_samples      = nrow(all_verhaak),
      Chi2_statistic = chi2_all$statistic,
      Chi2_p         = chi2_all$p.value,
      Cramers_V      = v_all,
      V_interpretation = v_interp_all,
      stringsAsFactors = FALSE
    ))
    message(sprintf("  Full cohort: Cramer's V = %.4f (%s)", v_all, v_interp_all))
  }
}

# --- 评审修订B3: IDH状态 x UPR亚型 Cramer's V ---
message("\n=== IDH status x UPR subtype Cramer's V (reviewer revision B3) ===")

if ("IDH_status" %in% colnames(heatmap_data)) {
  idh_cross_data <- heatmap_data %>%
    dplyr::filter(!is.na(IDH_status), !is.na(UPR_subtype))

  if (nrow(idh_cross_data) >= 20) {
    cross_idh <- table(idh_cross_data$UPR_subtype, idh_cross_data$IDH_status)
    message("\n  Cross-tabulation: UPR subtype x IDH status")
    print(cross_idh)

    cross_idh_prop <- prop.table(cross_idh, margin = 1) * 100
    message("\n  Row percentages:")
    print(round(cross_idh_prop, 1))

    chi2_idh <- chisq.test(cross_idh)
    cramers_v_idh <- cramerV(cross_idh)
    v_idh <- cramers_v_idh[1]
    v_interp_idh <- dplyr::case_when(
      v_idh < 0.1  ~ "negligible",
      v_idh < 0.3  ~ "small",
      v_idh < 0.5  ~ "medium",
      TRUE          ~ "large"
    )

    message(sprintf("\n  IDH x UPR: Chi-squared X2=%.2f, p=%.4e",
                    chi2_idh$statistic, chi2_idh$p.value))
    message(sprintf("  Cramer's V = %.4f (%s effect)", v_idh, v_interp_idh))

    if (v_idh >= 0.5) {
      message("  WARNING: Large Cramer's V suggests UPR subtypes may largely reflect IDH status.")
      message("  IDH-stratified analyses are critical for this study.")
    }

    confounding_results <- rbind(confounding_results, data.frame(
      Analysis       = "Full cohort: UPR x IDH status",
      N_samples      = nrow(idh_cross_data),
      Chi2_statistic = chi2_idh$statistic,
      Chi2_p         = chi2_idh$p.value,
      Cramers_V      = v_idh,
      V_interpretation = v_interp_idh,
      stringsAsFactors = FALSE
    ))

    # IDH-WT亚组内部: 排除IDH混杂后检验
    idhwt_only <- heatmap_data %>%
      dplyr::filter(IDH_status == "WT", !is.na(UPR_subtype))
    if (nrow(idhwt_only) >= 20 && length(unique(idhwt_only$UPR_subtype)) > 1) {
      message(sprintf("\n  IDH-WT only subtype distribution (n=%d):", nrow(idhwt_only)))
      print(table(idhwt_only$UPR_subtype))
    }
  }
}

write.csv(confounding_results, file.path(RES_DIR, "verhaak_confounding_analysis.csv"),
          row.names = FALSE)

# =============================================================================
# 4. 组合临床特征条形图 (Stacked bar plots)
# =============================================================================
message("=== Combined clinical barplots ===")

# 收集所有分类变量的比例数据
clinical_vars <- c("Grade", "IDH_status", "MGMT_status", "Subtype")
if ("codel_1p19q" %in% colnames(heatmap_data)) {
  clinical_vars <- c(clinical_vars, "codel_1p19q")
}

bar_data <- list()
for (var in clinical_vars) {
  if (var %in% colnames(heatmap_data)) {
    tmp <- heatmap_data %>%
      dplyr::filter(!is.na(!!sym(var))) %>%
      dplyr::count(UPR_subtype, !!sym(var)) %>%
      dplyr::group_by(UPR_subtype) %>%
      dplyr::mutate(pct = n / sum(n) * 100) %>%
      dplyr::ungroup()
    colnames(tmp) <- c("Subtype", "Category", "n", "pct")

    # 美化变量名
    var_label <- dplyr::case_when(
      var == "Grade"       ~ "WHO Grade",
      var == "IDH_status"  ~ "IDH Status",
      var == "MGMT_status" ~ "MGMT Methylation",
      var == "Subtype"     ~ "Verhaak Subtype",
      var == "codel_1p19q" ~ "1p/19q Codeletion",
      TRUE                 ~ var
    )
    tmp$Variable <- var_label
    bar_data[[var]] <- tmp
  }
}

if (length(bar_data) > 0) {
  all_bar <- do.call(rbind, bar_data)

  n_vars <- length(unique(all_bar$Variable))
  p_combined <- ggplot(all_bar, aes(x = Subtype, y = pct, fill = Category)) +
    geom_bar(stat = "identity", position = "stack", width = 0.7) +
    facet_wrap(~Variable, scales = "free_y", nrow = 1) +
    labs(y = "Percentage (%)", x = "UPR Subtype", fill = "Category",
         title = "Clinical Feature Distribution by UPR Subtype") +
    THEME_PUBLICATION +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      strip.text = element_text(face = "bold", size = 9),
      legend.position = "bottom",
      legend.box = "horizontal"
    ) +
    guides(fill = guide_legend(nrow = 2))

  fig_width <- max(12, n_vars * 3.5)
  ggsave(file.path(FIG_DIR, "Fig3_clinical_features_combined.pdf"),
         p_combined, width = fig_width, height = 6, useDingbats = FALSE)
}

# =============================================================================
# 5. 亚型内部特征汇总表
# =============================================================================
message("=== Generating subtype summary table ===")

summary_table <- heatmap_data %>%
  dplyr::group_by(UPR_subtype) %>%
  dplyr::summarise(
    N = dplyr::n(),
    Age_median = median(as.numeric(age_at_index), na.rm = TRUE),
    Age_IQR = paste0(
      round(quantile(as.numeric(age_at_index), 0.25, na.rm = TRUE), 1), "-",
      round(quantile(as.numeric(age_at_index), 0.75, na.rm = TRUE), 1)
    ),
    Male_pct = ifelse("gender" %in% colnames(heatmap_data),
                       round(sum(gender == "male", na.rm = TRUE) / sum(!is.na(gender)) * 100, 1),
                       NA),
    Grade_G4_pct = ifelse("Grade" %in% colnames(heatmap_data),
                           round(sum(Grade == "G4", na.rm = TRUE) / sum(!is.na(Grade)) * 100, 1),
                           NA),
    IDH_WT_pct = ifelse("IDH_status" %in% colnames(heatmap_data),
                         round(sum(IDH_status == "WT", na.rm = TRUE) / sum(!is.na(IDH_status)) * 100, 1),
                         NA),
    MGMT_methyl_pct = ifelse("MGMT_status" %in% colnames(heatmap_data),
                              round(sum(MGMT_status == "Methylated", na.rm = TRUE) / sum(!is.na(MGMT_status)) * 100, 1),
                              NA),
    Median_OS_days = median(OS.time, na.rm = TRUE),
    .groups = "drop"
  )

message("\nSubtype summary:")
print(as.data.frame(summary_table))
write.csv(summary_table, file.path(RES_DIR, "upr_subtype_summary_table.csv"), row.names = FALSE)

# =============================================================================
# 6. 保存所有结果
# =============================================================================
save(
  heatmap_data, clinical_stats, clinical_tests_df,
  confounding_results, summary_table,
  file = file.path(DATA_PROC, "clinical_characterization_results.RData")
)

message("\n=== Clinical characterization completed ===")
message(sprintf("  Association tests: %d variables tested", nrow(clinical_tests_df)))
if (nrow(confounding_results) > 0) {
  message(sprintf("  Verhaak confounding: Cramer's V = %.3f (%s)",
                  confounding_results$Cramers_V[1],
                  confounding_results$V_interpretation[1]))
}
message("Next step: Run 04_bulk_subtyping/04_immune_analysis.R")
