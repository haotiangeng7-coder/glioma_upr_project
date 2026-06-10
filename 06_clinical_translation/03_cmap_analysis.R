###############################################################################
# 03_cmap_analysis.R
# Part 4 - 临床转化 (3/4)
# CMap药物重定位: limma差异基因 + Top150 up/down + clue.io查询
# 输出: Figure 7E
###############################################################################

source("00_setup/config.R")
library(limma)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyr)
library(ggrepel)

set.seed(SEED)

# --- 加载数据 ---
load(file.path(DATA_PROC, "risk_model_final.RData"))
load(file.path(DATA_PROC, "tcga_glioma_expression.RData"))

# --- 输出目录 ---
fig_dir <- file.path(FIG_DIR, "part4_clinical")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. limma差异基因分析 (High vs Low risk)
# =============================================================================
message("=== Section 1: Differential Expression Analysis (High vs Low Risk) ===")

common_samples <- intersect(names(risk_score_train), colnames(expr_tpm_symbol))
risk_groups <- setNames(risk_group_train[common_samples], common_samples)

high_samples <- common_samples[risk_groups == "High"]
low_samples  <- common_samples[risk_groups == "Low"]

message(sprintf("High risk: %d samples, Low risk: %d samples",
                length(high_samples), length(low_samples)))

# log2(TPM+1) 表达矩阵
all_de_samples <- c(high_samples, low_samples)
expr_de <- log2(expr_tpm_symbol[, all_de_samples] + 1)

# 过滤低表达基因: 至少在10%样本中表达 > 0.5
keep_genes <- rowSums(expr_de > 0.5) >= ncol(expr_de) * 0.1
expr_de <- expr_de[keep_genes, ]
message(sprintf("Genes after filtering: %d", nrow(expr_de)))

# limma设计矩阵
group <- factor(c(rep("High", length(high_samples)),
                  rep("Low", length(low_samples))),
                levels = c("Low", "High"))
design <- model.matrix(~ 0 + group)
colnames(design) <- levels(group)

# 拟合
fit <- lmFit(expr_de, design)
contrast <- makeContrasts(High - Low, levels = design)
fit2 <- contrasts.fit(fit, contrast)
fit2 <- eBayes(fit2)

de_results <- topTable(fit2, number = Inf, sort.by = "t")
de_results$gene <- rownames(de_results)

# 显著DEGs
sig_de <- de_results %>%
  dplyr::filter(adj.P.Val < FDR_CUTOFF, abs(logFC) > LOGFC_CUTOFF)

n_up   <- sum(sig_de$logFC > 0)
n_down <- sum(sig_de$logFC < 0)
message(sprintf("Significant DEGs (|logFC|>%.1f, FDR<%.2f): %d (up: %d, down: %d)",
                LOGFC_CUTOFF, FDR_CUTOFF, nrow(sig_de), n_up, n_down))

write.csv(de_results, file.path(RES_DIR, "deg_high_vs_low_risk.csv"), row.names = FALSE)

# --- 火山图 ---
de_results$significance <- "NS"
de_results$significance[de_results$adj.P.Val < FDR_CUTOFF & de_results$logFC > LOGFC_CUTOFF] <- "Up"
de_results$significance[de_results$adj.P.Val < FDR_CUTOFF & de_results$logFC < -LOGFC_CUTOFF] <- "Down"
de_results$significance <- factor(de_results$significance, levels = c("Down", "NS", "Up"))

# 标注top基因
top_label_genes <- de_results %>%
  dplyr::filter(significance != "NS") %>%
  dplyr::arrange(adj.P.Val) %>%
  head(20) %>%
  dplyr::pull(gene)

de_results$label <- ifelse(de_results$gene %in% top_label_genes, de_results$gene, "")

p_volcano <- ggplot(de_results, aes(x = logFC, y = -log10(adj.P.Val), color = significance)) +
  geom_point(size = 0.8, alpha = 0.6) +
  geom_text_repel(aes(label = label), size = 2.5, max.overlaps = 20,
                  show.legend = FALSE) +
  scale_color_manual(values = c("Down" = "#4DBBD5", "NS" = "grey70", "Up" = "#E64B35")) +
  geom_hline(yintercept = -log10(FDR_CUTOFF), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-LOGFC_CUTOFF, LOGFC_CUTOFF), linetype = "dashed", color = "grey40") +
  labs(x = "log2 Fold Change (High / Low)",
       y = "-log10(adj. P-value)",
       title = "Differential Expression: High vs Low Risk",
       color = "Significance") +
  annotate("text", x = max(de_results$logFC) * 0.7, y = max(-log10(de_results$adj.P.Val)) * 0.95,
           label = paste0("Up: ", n_up), color = "#E64B35", size = 4) +
  annotate("text", x = min(de_results$logFC) * 0.7, y = max(-log10(de_results$adj.P.Val)) * 0.95,
           label = paste0("Down: ", n_down), color = "#4DBBD5", size = 4) +
  THEME_PUBLICATION +
  theme(legend.position = "bottom")

ggsave(file.path(fig_dir, "Fig7E_volcano.pdf"), p_volcano, width = 8, height = 7)

# =============================================================================
# 2. 准备CMap查询Signature (Top 150 up + Top 150 down)
# =============================================================================
message("\n=== Section 2: Preparing CMap Signature ===")

# 取显著上调基因中logFC最大的150个
up_genes <- de_results %>%
  dplyr::filter(logFC > 0, adj.P.Val < FDR_CUTOFF) %>%
  dplyr::arrange(desc(logFC)) %>%
  head(150) %>%
  dplyr::pull(gene)

# 取显著下调基因中logFC最小的150个
down_genes <- de_results %>%
  dplyr::filter(logFC < 0, adj.P.Val < FDR_CUTOFF) %>%
  dplyr::arrange(logFC) %>%
  head(150) %>%
  dplyr::pull(gene)

message(sprintf("CMap signature: %d up genes, %d down genes",
                length(up_genes), length(down_genes)))

# 如果显著基因不足150个，放宽标准
if (length(up_genes) < 150) {
  message(sprintf("  Note: only %d significant up genes (< 150 target)", length(up_genes)))
  up_genes_relaxed <- de_results %>%
    dplyr::filter(logFC > 0) %>%
    dplyr::arrange(desc(logFC)) %>%
    head(150) %>%
    dplyr::pull(gene)
  message(sprintf("  Using top 150 by logFC (relaxed): %d genes", length(up_genes_relaxed)))
}

if (length(down_genes) < 150) {
  message(sprintf("  Note: only %d significant down genes (< 150 target)", length(down_genes)))
  down_genes_relaxed <- de_results %>%
    dplyr::filter(logFC < 0) %>%
    dplyr::arrange(logFC) %>%
    head(150) %>%
    dplyr::pull(gene)
  message(sprintf("  Using top 150 by logFC (relaxed): %d genes", length(down_genes_relaxed)))
}

# 保存CMap查询文件 (每行一个基因)
writeLines(up_genes, file.path(RES_DIR, "cmap_up_genes.txt"))
writeLines(down_genes, file.path(RES_DIR, "cmap_down_genes.txt"))

# GMT格式
cat("UP_IN_HIGH_RISK\tna\t", paste(up_genes, collapse = "\t"), "\n",
    file = file.path(RES_DIR, "cmap_signature.gmt"), sep = "")
cat("DOWN_IN_HIGH_RISK\tna\t", paste(down_genes, collapse = "\t"), "\n",
    file = file.path(RES_DIR, "cmap_signature.gmt"), sep = "", append = TRUE)

message("\nCMap query files saved:")
message("  Up genes:   ", file.path(RES_DIR, "cmap_up_genes.txt"))
message("  Down genes: ", file.path(RES_DIR, "cmap_down_genes.txt"))
message("  GMT:        ", file.path(RES_DIR, "cmap_signature.gmt"))

# =============================================================================
# 3. 尝试clue.io API程序化查询
# =============================================================================
message("\n=== Section 3: Clue.io API Query (if API key available) ===")

# clue.io API需要用户注册后获取API key
clue_api_key <- Sys.getenv("CLUE_API_KEY")

if (nchar(clue_api_key) > 0) {
  message("Clue API key found, attempting programmatic query...")

  tryCatch({
    library(httr)
    library(jsonlite)

    # 上调基因 -> Entrez ID (clue.io推荐Entrez或Symbol)
    # 直接使用gene symbol
    up_payload <- paste(up_genes, collapse = "\n")
    down_payload <- paste(down_genes, collapse = "\n")

    # clue.io L1000 Query API
    # POST https://api.clue.io/api/query
    query_body <- list(
      tool       = "sig_fastgutc_tool",
      up         = up_genes,
      down       = down_genes,
      name       = "UIRS_high_vs_low_risk",
      species    = "human"
    )

    response <- httr::POST(
      url = "https://api.clue.io/api/query",
      httr::add_headers(
        "user_key" = clue_api_key,
        "Content-Type" = "application/json"
      ),
      body = jsonlite::toJSON(query_body, auto_unbox = TRUE),
      encode = "raw"
    )

    if (httr::status_code(response) == 200) {
      query_result <- httr::content(response, as = "parsed")
      message("CMap query submitted successfully!")
      message("Job ID: ", query_result$job_id)
      message("Check results at: https://clue.io/results/", query_result$job_id)

      # 保存job信息
      writeLines(c(
        paste0("Job ID: ", query_result$job_id),
        paste0("Submitted: ", Sys.time()),
        paste0("Up genes: ", length(up_genes)),
        paste0("Down genes: ", length(down_genes))
      ), file.path(RES_DIR, "cmap_query_info.txt"))
    } else {
      message("CMap API query failed: HTTP ", httr::status_code(response))
      message("Response: ", httr::content(response, as = "text"))
    }
  }, error = function(e) {
    message("CMap API query error: ", e$message)
  })
} else {
  message("CLUE_API_KEY not set in environment.")
  message("To use programmatic CMap query:")
  message("  1. Register at https://clue.io")
  message("  2. Get API key from https://clue.io/api")
  message("  3. export CLUE_API_KEY='your_key_here'")
  message("")
  message("Manual query instructions:")
  message("  1. Go to https://clue.io/query (L1000 CMap)")
  message("  2. Upload cmap_up_genes.txt as 'Up' gene set")
  message("  3. Upload cmap_down_genes.txt as 'Down' gene set")
  message("  4. Download results and place in: ", file.path(DATA_RAW, "cmap_results/"))
}

# =============================================================================
# 4. 处理CMap结果（如果已下载）
# =============================================================================
message("\n=== Section 4: Processing CMap Results ===")

cmap_result_dir <- file.path(DATA_RAW, "cmap_results")
cmap_processed <- FALSE

if (dir.exists(cmap_result_dir)) {
  cmap_files <- list.files(cmap_result_dir,
                           pattern = "\\.csv$|\\.txt$|\\.gct$|\\.tsv$",
                           full.names = TRUE)

  if (length(cmap_files) > 0) {
    message(sprintf("Found %d CMap result files", length(cmap_files)))

    for (f in cmap_files) {
      message("  Processing: ", basename(f))

      cmap_data <- tryCatch(
        read.csv(f, stringsAsFactors = FALSE),
        error = function(e) {
          tryCatch(
            read.delim(f, stringsAsFactors = FALSE),
            error = function(e2) NULL
          )
        }
      )

      if (is.null(cmap_data)) next

      # 寻找连接性评分列
      score_col <- grep("^(score|connectivity|tau|cs)$",
                        colnames(cmap_data), ignore.case = TRUE, value = TRUE)
      name_col  <- grep("^(name|pert_iname|compound|drug)$",
                        colnames(cmap_data), ignore.case = TRUE, value = TRUE)

      if (length(score_col) > 0) {
        sc <- score_col[1]
        cmap_data <- cmap_data %>% dplyr::arrange(!!sym(sc))

        # Top 20 候选药物（连接性最负 = 最可能逆转高风险表型）
        top_candidates <- head(cmap_data, 20)

        message("Top 20 CMap candidates (most negative connectivity score):")
        cols_to_show <- intersect(c(name_col, sc, "cell_id", "dose", "moa"),
                                  colnames(cmap_data))
        if (length(cols_to_show) == 0) cols_to_show <- colnames(cmap_data)[1:min(5, ncol(cmap_data))]
        print(top_candidates[, cols_to_show])

        write.csv(top_candidates,
                  file.path(RES_DIR, "cmap_top_candidates.csv"), row.names = FALSE)

        # --- CMap结果条形图 ---
        if (length(name_col) > 0) {
          plot_data <- top_candidates
          plot_data$drug_name <- plot_data[[name_col[1]]]
          plot_data$score <- as.numeric(plot_data[[sc]])

          p_cmap <- ggplot(plot_data, aes(x = reorder(drug_name, score), y = score)) +
            geom_bar(stat = "identity",
                     fill = ifelse(plot_data$score < 0, "#4DBBD5", "#E64B35"),
                     width = 0.7) +
            coord_flip() +
            labs(x = "", y = "Connectivity Score",
                 title = "Top 20 CMap Candidates\n(negative = reverses high-risk signature)") +
            geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
            THEME_PUBLICATION +
            theme(axis.text.y = element_text(size = 9))

          ggsave(file.path(fig_dir, "Fig7E_cmap_candidates.pdf"),
                 p_cmap, width = 8, height = 7)
        }

        cmap_processed <- TRUE
        break  # 处理第一个有效文件即可
      }
    }
  }
}

if (!cmap_processed) {
  message("CMap results not available yet.")
  message("Results will be visualized after CMap query completion.")
}

# =============================================================================
# 5. UPR靶向药物文献整理
# =============================================================================
message("\n=== Section 5: Literature-based UPR-targeting Drug Candidates ===")

er_stress_drugs <- data.frame(
  Drug = c("GSK2606414", "ISRIB", "4u8C", "STF-083010",
           "Bortezomib", "Carfilzomib", "Thapsigargin",
           "Tunicamycin", "17-AAG", "Celecoxib"),
  Target = c("PERK inhibitor",
             "ISR inhibitor (eIF2B activator)",
             "IRE1a RNase inhibitor",
             "IRE1a RNase inhibitor",
             "Proteasome inhibitor",
             "Proteasome inhibitor",
             "SERCA inhibitor",
             "N-glycosylation inhibitor",
             "HSP90 inhibitor",
             "COX-2 inhibitor"),
  UPR_Arm = c("PERK-ATF4", "ISR/PERK", "IRE1-XBP1", "IRE1-XBP1",
              "Proteostasis", "Proteostasis", "ER Ca2+", "ER protein folding",
              "Chaperone", "ER stress modulator"),
  Mechanism = c("Blocks PERK-ATF4 pro-survival signaling",
                "Restores global translation despite eIF2a phosphorylation",
                "Blocks XBP1 splicing, restores immune function",
                "Blocks XBP1 splicing",
                "Overwhelms proteasome, triggers UPR-mediated apoptosis",
                "Selective proteasome inhibition, UPR overload",
                "Depletes ER calcium, triggers terminal UPR",
                "Blocks protein glycosylation, triggers UPR",
                "Disrupts chaperone-client interaction",
                "Modulates ER stress through COX-2 pathway"),
  Evidence = c("Preclinical GBM studies (Axten et al. 2012)",
               "Enhances anti-tumor immunity (Nguyen et al. 2018)",
               "Restores DC function in TME (Cubillos-Ruiz et al. 2015)",
               "Reduces tumor growth in mouse GBM models",
               "FDA-approved, Phase II glioma (Phuphanich et al. 2010)",
               "FDA-approved for myeloma, glioma preclinical",
               "Research tool compound, ER stress inducer",
               "Research tool compound, ER stress inducer",
               "Phase II recurrent GBM (Sauvageot et al. 2009)",
               "Epidemiological glioma risk reduction"),
  Clinical_Status = c("Preclinical", "Preclinical", "Preclinical", "Preclinical",
                      "Phase II", "Preclinical (glioma)", "Tool compound", "Tool compound",
                      "Phase II", "Epidemiological"),
  stringsAsFactors = FALSE
)

write.csv(er_stress_drugs, file.path(RES_DIR, "er_stress_drug_candidates.csv"),
          row.names = FALSE)

# --- 药物靶向通路概览图 ---
p_drug_overview <- ggplot(er_stress_drugs,
                          aes(x = reorder(Drug, -seq_along(Drug)),
                              y = 1, fill = UPR_Arm)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(aes(label = Drug), size = 3, fontface = "bold") +
  geom_text(aes(y = 0.7, label = Target), size = 2.3, color = "grey30") +
  scale_fill_brewer(palette = "Set3") +
  labs(x = "", y = "", fill = "UPR Pathway",
       title = "UPR-targeting Drug Candidates for Glioma") +
  coord_flip() +
  THEME_PUBLICATION +
  theme(axis.text   = element_blank(),
        axis.ticks  = element_blank(),
        panel.grid  = element_blank(),
        panel.border = element_blank())

ggsave(file.path(fig_dir, "Fig7E_drug_candidates.pdf"), p_drug_overview,
       width = 10, height = 6)

# =============================================================================
# 6. DEG通路富集辅助分析
# =============================================================================
message("\n=== Section 6: DEG Pathway Enrichment ===")

tryCatch({
  library(clusterProfiler)
  library(org.Hs.eg.db)

  # 上调基因GO富集
  up_entrez <- tryCatch({
    bitr(up_genes, fromType = "SYMBOL", toType = "ENTREZID",
         OrgDb = org.Hs.eg.db)$ENTREZID
  }, error = function(e) NULL)

  down_entrez <- tryCatch({
    bitr(down_genes, fromType = "SYMBOL", toType = "ENTREZID",
         OrgDb = org.Hs.eg.db)$ENTREZID
  }, error = function(e) NULL)

  if (!is.null(up_entrez) && length(up_entrez) > 10) {
    ego_up <- enrichGO(gene = up_entrez,
                       OrgDb = org.Hs.eg.db,
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       readable = TRUE)

    if (nrow(ego_up@result) > 0) {
      p_go_up <- dotplot(ego_up, showCategory = 15,
                         title = "Up in High Risk: GO Biological Process")
      ggsave(file.path(fig_dir, "Fig7E_GO_up_in_high_risk.pdf"),
             p_go_up, width = 10, height = 8)
    }
  }

  if (!is.null(down_entrez) && length(down_entrez) > 10) {
    ego_down <- enrichGO(gene = down_entrez,
                         OrgDb = org.Hs.eg.db,
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         readable = TRUE)

    if (nrow(ego_down@result) > 0) {
      p_go_down <- dotplot(ego_down, showCategory = 15,
                           title = "Down in High Risk: GO Biological Process")
      ggsave(file.path(fig_dir, "Fig7E_GO_down_in_high_risk.pdf"),
             p_go_down, width = 10, height = 8)
    }
  }
}, error = function(e) {
  message("GO enrichment analysis failed: ", e$message)
  message("clusterProfiler or org.Hs.eg.db may need to be installed")
})

# =============================================================================
# 7. 保存结果
# =============================================================================
message("\n=== Saving results ===")

save(de_results, up_genes, down_genes, er_stress_drugs,
     file = file.path(DATA_PROC, "cmap_analysis_results.RData"))

message("\n=== CMap analysis completed ===")
message("Figures saved to: ", fig_dir)
message("  Fig7E: Volcano plot + CMap candidates + drug overview")
message("  CMap input files ready for clue.io query")
message("\nNext: Run 06_clinical_translation/04_hub_gene_analysis.R")
