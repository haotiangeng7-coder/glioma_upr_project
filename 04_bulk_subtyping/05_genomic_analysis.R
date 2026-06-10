###############################################################################
# 05_genomic_analysis.R
# 各UPR亚型的体细胞突变分析（Figure 4D）
# - maftools分析体细胞突变
# - 各亚型Top20突变瀑布图(oncoplot)
# - 亚型间差异突变(mafCompare + Forest plot)
# - TMB比较（Kruskal-Wallis + pairwise Wilcoxon）
###############################################################################

source("00_setup/config.R")
library(maftools)
library(ggplot2)
library(ggpubr)
library(dplyr)

set.seed(SEED)
load(file.path(DATA_PROC, "consensus_clustering_results.RData"))

# =============================================================================
# 1. 加载突变数据 & 构建barcode映射
# =============================================================================
message("=== Loading mutation data ===")

load(file.path(DATA_PROC, "tcga_glioma_mutation.RData"))

# 提取TCGA barcode前12位（患者ID）以匹配突变数据
cluster_df$patient_id <- substr(cluster_df$barcode, 1, 12)

# 为MAF准备临床注释
clinical_for_maf <- cluster_df %>%
  dplyr::select(patient_id, UPR_subtype) %>%
  dplyr::distinct() %>%
  dplyr::rename(Tumor_Sample_Barcode = patient_id)

# 截短MAF中的barcode以匹配
maf_glioma$Tumor_Sample_Barcode_short <- substr(maf_glioma$Tumor_Sample_Barcode, 1, 12)

subtypes_unique <- sort(unique(cluster_df$UPR_subtype))
message(sprintf("Subtypes: %s", paste(subtypes_unique, collapse = ", ")))
message(sprintf("Patients with subtype annotation: %d", nrow(clinical_for_maf)))

# =============================================================================
# 2. 各亚型Top20突变瀑布图 (oncoplot)
# =============================================================================
message("\n=== Mutation landscape by subtype (oncoplot) ===")

# 存储各亚型的maftools对象和TMB，避免重复创建
maf_obj_list <- list()

for (st in subtypes_unique) {
  message(sprintf("\n--- Processing %s ---", st))

  patients_st <- cluster_df$patient_id[cluster_df$UPR_subtype == st]

  # 过滤该亚型的MAF行
  maf_st <- maf_glioma[maf_glioma$Tumor_Sample_Barcode_short %in% patients_st, ]
  # 将barcode截短以匹配clinical_for_maf中的Tumor_Sample_Barcode
  maf_st$Tumor_Sample_Barcode <- maf_st$Tumor_Sample_Barcode_short

  if (nrow(maf_st) < 10) {
    message(sprintf("  Skipping %s: insufficient mutations (%d rows)", st, nrow(maf_st)))
    next
  }

  tryCatch({
    maf_obj <- read.maf(maf = maf_st, clinicalData = clinical_for_maf, verbose = FALSE)
    maf_obj_list[[st]] <- maf_obj

    # Top20突变瀑布图
    st_safe <- gsub("-", "_", st)
    pdf(file.path(FIG_DIR, paste0("Fig4D_oncoplot_", st_safe, ".pdf")),
        width = 12, height = 8)
    oncoplot(maf = maf_obj, top = 20,
             title = paste0("Top 20 Mutated Genes: ", st),
             fontSize = 0.6)
    dev.off()

    # 记录突变摘要
    gene_summary <- getGeneSummary(maf_obj)
    n_samples <- as.numeric(maf_obj@summary[maf_obj@summary$ID == "Samples", "summary"])
    message(sprintf("  %s: %d patients, %d mutated genes",
                    st, n_samples, nrow(gene_summary)))

    # TMB预计算
    tmb_st <- tmb(maf_obj, captureSize = 50)
    message(sprintf("  Median TMB = %.2f mut/Mb", median(tmb_st$total_perMB, na.rm = TRUE)))

  }, error = function(e) {
    message(sprintf("  Error processing %s: %s", st, e$message))
  })
}

# =============================================================================
# 3. 亚型间差异突变分析 (mafCompare + Forest plot)
# =============================================================================
message("\n=== Comparative mutation analysis (mafCompare) ===")

# 构建全量MAF对象（带亚型注释），用于subsetMaf
tryCatch({
  all_patients <- unique(maf_glioma$Tumor_Sample_Barcode_short)
  all_clinical <- data.frame(
    Tumor_Sample_Barcode = all_patients,
    stringsAsFactors = FALSE
  )
  all_clinical <- merge(all_clinical, clinical_for_maf,
                        by = "Tumor_Sample_Barcode", all.x = TRUE)

  maf_all <- read.maf(maf = maf_glioma, clinicalData = all_clinical, verbose = FALSE)

  # 对所有成对亚型进行差异突变比较
  diff_mut_results <- list()

  if (length(subtypes_unique) >= 2) {
    pairs <- combn(subtypes_unique, 2, simplify = FALSE)

    for (pair in pairs) {
      st1 <- pair[1]
      st2 <- pair[2]
      message(sprintf("\n  Comparing %s vs %s", st1, st2))

      st1_patients <- cluster_df$patient_id[cluster_df$UPR_subtype == st1]
      st2_patients <- cluster_df$patient_id[cluster_df$UPR_subtype == st2]

      maf_sub1 <- subsetMaf(maf_all, tsb = st1_patients, verbose = FALSE)
      maf_sub2 <- subsetMaf(maf_all, tsb = st2_patients, verbose = FALSE)

      tryCatch({
        diff_mut <- mafCompare(m1 = maf_sub1, m2 = maf_sub2,
                                m1Name = st1, m2Name = st2,
                                minMut = 5)

        pair_name <- paste0(gsub("-", "_", st1), "_vs_", gsub("-", "_", st2))
        diff_mut_results[[pair_name]] <- diff_mut

        n_sig <- sum(diff_mut$results$pval < PVALUE_CUTOFF, na.rm = TRUE)
        n_sig_adj <- sum(diff_mut$results$adjPval < FDR_CUTOFF, na.rm = TRUE)
        message(sprintf("    Differentially mutated genes: %d (p<%.2f), %d (FDR<%.2f)",
                        n_sig, PVALUE_CUTOFF, n_sig_adj, FDR_CUTOFF))

        # Forest plot
        if (n_sig > 0) {
          pdf(file.path(FIG_DIR, paste0("Fig4D_forest_", pair_name, ".pdf")),
              width = 10, height = max(6, min(n_sig * 0.4, 12)))
          forestPlot(mafCompareRes = diff_mut, pVal = PVALUE_CUTOFF,
                     titleSize = 1.2)
          dev.off()
        }

        # 保存差异结果
        write.csv(diff_mut$results,
                  file.path(RES_DIR, paste0("diff_mutations_", pair_name, ".csv")),
                  row.names = FALSE)

      }, error = function(e) {
        message(sprintf("    mafCompare error: %s", e$message))
      })
    }
  }

  # Co-oncoplot: UPR-favorable vs UPR-high-risk（预后最好 vs 最差的对比）
  if ("UPR-favorable" %in% subtypes_unique && "UPR-high-risk" %in% subtypes_unique) {
    fav_patients <- cluster_df$patient_id[cluster_df$UPR_subtype == "UPR-favorable"]
    hr_patients <- cluster_df$patient_id[cluster_df$UPR_subtype == "UPR-high-risk"]

    maf_fav <- subsetMaf(maf_all, tsb = fav_patients, verbose = FALSE)
    maf_hr <- subsetMaf(maf_all, tsb = hr_patients, verbose = FALSE)

    pair_key <- "UPR_favorable_vs_UPR_high_risk"
    if (pair_key %in% names(diff_mut_results)) {
      top_genes <- head(diff_mut_results[[pair_key]]$results$Hugo_Symbol, 20)
    } else {
      # 合并两组的高频突变基因
      top_genes <- unique(c(
        head(getGeneSummary(maf_fav)$Hugo_Symbol, 10),
        head(getGeneSummary(maf_hr)$Hugo_Symbol, 10)
      ))
    }

    if (length(top_genes) > 0) {
      pdf(file.path(FIG_DIR, "Fig4D_co_oncoplot.pdf"), width = 16, height = 10)
      coOncoplot(m1 = maf_fav, m2 = maf_hr,
                 m1Name = "UPR-favorable", m2Name = "UPR-high-risk",
                 genes = top_genes)
      dev.off()
    }
  }

}, error = function(e) {
  message("Comparative mutation analysis error: ", e$message)
})

# =============================================================================
# 4. TMB比较（Kruskal-Wallis + pairwise Wilcoxon）
# =============================================================================
message("\n=== TMB comparison ===")

tmb_all <- NULL

tryCatch({
  tmb_list <- vector("list", length(subtypes_unique))
  names(tmb_list) <- subtypes_unique

  for (st in subtypes_unique) {
    patients_st <- cluster_df$patient_id[cluster_df$UPR_subtype == st]
    maf_st <- maf_glioma[maf_glioma$Tumor_Sample_Barcode_short %in% patients_st, ]

  maf_st$Tumor_Sample_Barcode <- maf_st$Tumor_Sample_Barcode_short
    if (nrow(maf_st) > 0) {
      maf_obj <- read.maf(maf = maf_st, verbose = FALSE)
      tmb_st <- tmb(maf_obj, captureSize = 50)
      tmb_st$UPR_subtype <- st
      tmb_list[[st]] <- tmb_st
    }
  }

  tmb_all <- do.call(rbind, tmb_list[!sapply(tmb_list, is.null)])

  if (!is.null(tmb_all) && nrow(tmb_all) > 0) {
    # 统计检验
    tmb_kw <- kruskal.test(total_perMB ~ UPR_subtype, data = tmb_all)
    message(sprintf("  TMB Kruskal-Wallis p = %.4e", tmb_kw$p.value))

    # pairwise比较
    tmb_pw <- pairwise.wilcox.test(tmb_all$total_perMB, tmb_all$UPR_subtype,
                                    p.adjust.method = "BH")
    message("  Pairwise Wilcoxon (BH-adjusted):")
    print(tmb_pw$p.value)

    # 各亚型TMB统计摘要
    tmb_summary <- tmb_all %>%
      group_by(UPR_subtype) %>%
      summarise(
        n = n(),
        median_TMB = median(total_perMB, na.rm = TRUE),
        mean_TMB = mean(total_perMB, na.rm = TRUE),
        sd_TMB = sd(total_perMB, na.rm = TRUE),
        .groups = "drop"
      )
    message("\n  TMB Summary:")
    print(as.data.frame(tmb_summary))

    # 箱线图 + pairwise比较
    comparisons_list <- combn(sort(unique(as.character(tmb_all$UPR_subtype))), 2, simplify = FALSE)

    p_tmb <- ggplot(tmb_all, aes(x = UPR_subtype, y = total_perMB, fill = UPR_subtype)) +
      geom_boxplot(outlier.size = 0.5, width = 0.6) +
      geom_jitter(width = 0.15, size = 0.3, alpha = 0.3) +
      stat_compare_means(method = "kruskal.test", label.y.npc = 0.95, size = 3.5) +
      stat_compare_means(comparisons = comparisons_list, method = "wilcox.test",
                         p.adjust.method = "BH", label = "p.signif",
                         size = 3, hide.ns = TRUE, step.increase = 0.08) +
      scale_fill_manual(values = COLORS_SUBTYPE) +
      scale_y_log10() +
      labs(x = "UPR Subtype", y = "TMB (mutations/Mb)",
           title = "Tumor Mutation Burden by UPR Subtype") +
      THEME_PUBLICATION +
      theme(legend.position = "none")

    ggsave(file.path(FIG_DIR, "Fig4D_tmb_comparison.pdf"), p_tmb, width = 6, height = 5)

    # 保存TMB结果
    write.csv(tmb_all, file.path(RES_DIR, "tmb_by_subtype.csv"), row.names = FALSE)
    write.csv(tmb_summary, file.path(RES_DIR, "tmb_summary_by_subtype.csv"), row.names = FALSE)
  }
}, error = function(e) {
  message("TMB comparison error: ", e$message)
})

# =============================================================================
# 5. 关键胶质瘤驱动基因突变频率比较
# =============================================================================
message("\n=== Key glioma driver gene mutation frequency ===")

tryCatch({
  # 胶质瘤关键驱动基因
  driver_genes <- c("IDH1", "IDH2", "TP53", "ATRX", "CIC", "FUBP1",
                    "NOTCH1", "PIK3CA", "PIK3R1", "PTEN", "EGFR",
                    "NF1", "RB1", "CDK4", "CDKN2A", "PDGFRA",
                    "TERT", "H3F3A", "HIST1H3B", "BRAF")

  # 按亚型计算各驱动基因的突变频率
  driver_freq <- data.frame()

  for (st in subtypes_unique) {
    patients_st <- cluster_df$patient_id[cluster_df$UPR_subtype == st]
    maf_st_rows <- maf_glioma[maf_glioma$Tumor_Sample_Barcode_short %in% patients_st, ]
    maf_st_rows$Tumor_Sample_Barcode <- maf_st_rows$Tumor_Sample_Barcode_short

  maf_st$Tumor_Sample_Barcode <- maf_st$Tumor_Sample_Barcode_short
    if (nrow(maf_st_rows) == 0) next

    maf_obj <- read.maf(maf = maf_st_rows, verbose = FALSE)
    n_samples <- as.numeric(maf_obj@summary[maf_obj@summary$ID == "Samples", "summary"])
    gene_summ <- getGeneSummary(maf_obj)

    for (gene in driver_genes) {
      if (gene %in% gene_summ$Hugo_Symbol) {
        n_mutated <- gene_summ$MutatedSamples[gene_summ$Hugo_Symbol == gene]
      } else {
        n_mutated <- 0
      }
      driver_freq <- rbind(driver_freq, data.frame(
        Gene = gene,
        UPR_subtype = st,
        n_mutated = n_mutated,
        n_total = n_samples,
        freq = n_mutated / n_samples * 100,
        stringsAsFactors = FALSE
      ))
    }
  }

  if (nrow(driver_freq) > 0) {
    # 热图展示驱动基因突变频率
    freq_mat <- driver_freq %>%
      dplyr::select(Gene, UPR_subtype, freq) %>%
      tidyr::pivot_wider(names_from = UPR_subtype, values_from = freq) %>%
      tibble::column_to_rownames("Gene") %>%
      as.matrix()

    # 按最大频率排序
    freq_mat <- freq_mat[order(-apply(freq_mat, 1, max)), , drop = FALSE]

    # Fisher检验各基因在亚型间的频率差异
    fisher_p <- numeric(nrow(freq_mat))
    names(fisher_p) <- rownames(freq_mat)
    for (gene in rownames(freq_mat)) {
      gene_data <- driver_freq %>% filter(Gene == gene)
      cont_table <- matrix(c(gene_data$n_mutated, gene_data$n_total - gene_data$n_mutated),
                           nrow = nrow(gene_data), ncol = 2)
      if (nrow(cont_table) >= 2 && all(rowSums(cont_table) > 0)) {
        fisher_p[gene] <- fisher.test(cont_table, simulate.p.value = TRUE, B = 10000)$p.value
      } else {
        fisher_p[gene] <- NA
      }
    }
    fisher_padj <- p.adjust(fisher_p, method = "BH")

    driver_freq_summary <- data.frame(
      Gene = names(fisher_p),
      fisher_p = fisher_p,
      fisher_padj = fisher_padj,
      stringsAsFactors = FALSE
    )
    driver_freq_summary <- merge(
      driver_freq_summary,
      driver_freq %>%
        dplyr::select(Gene, UPR_subtype, freq) %>%
        tidyr::pivot_wider(names_from = UPR_subtype, values_from = freq,
                           names_prefix = "freq_"),
      by = "Gene"
    )

    write.csv(driver_freq_summary,
              file.path(RES_DIR, "driver_gene_mutation_freq.csv"), row.names = FALSE)

    # 驱动基因频率条形图
    library(tidyr)
    p_driver <- ggplot(driver_freq, aes(x = reorder(Gene, -freq), y = freq,
                                         fill = UPR_subtype)) +
      geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
      scale_fill_manual(values = COLORS_SUBTYPE) +
      labs(x = "Driver Gene", y = "Mutation Frequency (%)",
           fill = "UPR Subtype",
           title = "Glioma Driver Gene Mutation Frequency by UPR Subtype") +
      THEME_PUBLICATION +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

    ggsave(file.path(FIG_DIR, "Fig4D_driver_gene_freq.pdf"), p_driver,
           width = 12, height = 6)

    message(sprintf("  Significantly different driver genes (FDR < %.2f): %d",
                    FDR_CUTOFF, sum(fisher_padj < FDR_CUTOFF, na.rm = TRUE)))
  }
}, error = function(e) {
  message("Driver gene analysis error: ", e$message)
})

# =============================================================================
# 6. 保存所有基因组分析结果
# =============================================================================
message("\n=== Saving genomic analysis results ===")

save_objects <- c("tmb_all", "clinical_for_maf")
if (exists("diff_mut_results")) save_objects <- c(save_objects, "diff_mut_results")
if (exists("driver_freq_summary")) save_objects <- c(save_objects, "driver_freq_summary")

save(list = save_objects,
     file = file.path(DATA_PROC, "genomic_analysis_results.RData"))

message("\n=== Genomic analysis completed (Figure 4D) ===")
message("Next step: Run 04_bulk_subtyping/06_pathway_analysis.R")
