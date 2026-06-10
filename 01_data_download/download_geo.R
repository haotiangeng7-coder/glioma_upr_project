###############################################################################
# download_geo.R
# 下载GEO验证集数据
# GSE16011 (284例), GSE43378 (50例)
###############################################################################

source("00_setup/config.R")
library(GEOquery)
library(dplyr)
library(limma)

geo_dir <- file.path(DATA_RAW, "GEO")
dir.create(geo_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. GSE16011 — 284例胶质瘤, 表达芯片
# =============================================================================
message("=== Downloading GSE16011 ===")

gse16011 <- getGEO("GSE16011", destdir = geo_dir, GSEMatrix = TRUE, getGPL = TRUE)
gse16011 <- gse16011[[1]]

# 提取表达矩阵
expr_16011 <- exprs(gse16011)

# 检查是否需要log2转换
if (max(expr_16011, na.rm = TRUE) > 50) {
  message("  Log2 transforming GSE16011...")
  expr_16011 <- log2(expr_16011 + 1)
}

# 探针注释 → Gene Symbol
feat_16011 <- fData(gse16011)
# 查找Gene Symbol列
gene_col <- grep("Gene.Symbol|gene_assignment|Symbol", colnames(feat_16011), value = TRUE)[1]
if (is.na(gene_col)) {
  message("  Available annotation columns: ", paste(colnames(feat_16011), collapse = ", "))
  # Try broader column search
  gene_col <- grep("gene", colnames(feat_16011), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(gene_col)) {
    message("  Trying to extract gene symbols from GPL...")
    tryCatch({
      gpl <- getGEO(annotation(gse16011), destdir = geo_dir)
      gpl_table <- Table(gpl)
      gene_col_gpl <- grep("Gene.Symbol|Symbol|gene_assignment|GENE_SYMBOL|ORF",
                           colnames(gpl_table), ignore.case = TRUE, value = TRUE)[1]
      if (!is.na(gene_col_gpl)) {
        gene_symbols <- gpl_table[[gene_col_gpl]]
        names(gene_symbols) <- gpl_table$ID
        gene_symbols <- gene_symbols[rownames(expr_16011)]
      } else {
        # Last resort: use probe IDs as gene names
        message("  WARNING: No gene symbol column found. Using probe IDs.")
        gene_symbols <- rownames(expr_16011)
      }
    }, error = function(e) {
      message("  GPL download failed: ", e$message)
      gene_symbols <<- rownames(expr_16011)
    })
  } else {
    gene_symbols <- feat_16011[[gene_col]]
  }
} else {
  gene_symbols <- feat_16011[[gene_col]]
}

# 探针到基因映射（取每个基因最大IQR的探针）
if (exists("gene_symbols")) {
  # 处理多个基因映射（取第一个）
  gene_symbols <- gsub(" ///.*", "", gene_symbols)
  gene_symbols <- gsub("///.*", "", gene_symbols)

  valid_idx <- !is.na(gene_symbols) & gene_symbols != "" & gene_symbols != "---"
  expr_16011_valid <- expr_16011[valid_idx, ]
  genes_valid <- gene_symbols[valid_idx]

  # 对重复基因取IQR最大的探针
  iqr_vals <- apply(expr_16011_valid, 1, IQR, na.rm = TRUE)
  probe_df <- data.frame(
    probe = rownames(expr_16011_valid),
    gene = genes_valid,
    iqr = iqr_vals,
    stringsAsFactors = FALSE
  )
  probe_df <- probe_df %>%
    dplyr::group_by(gene) %>%
    dplyr::slice_max(iqr, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()

  expr_16011_gene <- expr_16011_valid[probe_df$probe, ]
  rownames(expr_16011_gene) <- probe_df$gene

  message(sprintf("  GSE16011: %d genes x %d samples", nrow(expr_16011_gene), ncol(expr_16011_gene)))
}

# 临床数据
clin_16011 <- pData(gse16011)

# 整理临床数据（GSE16011的临床列名可能不同，需要检查）
message("  GSE16011 clinical columns: ", paste(head(colnames(clin_16011), 20), collapse = ", "))

# 提取关键信息
clin_16011_clean <- clin_16011 %>%
  dplyr::mutate(
    Sample_ID = rownames(clin_16011)
  )

# 尝试提取生存信息（列名可能为characteristics_ch1等）
surv_cols <- grep("survival|os|overall|vital|death|follow", colnames(clin_16011), ignore.case = TRUE, value = TRUE)
message("  Survival-related columns: ", paste(surv_cols, collapse = ", "))

# =============================================================================
# 2. GSE43378 — 50例胶质瘤
# =============================================================================
message("\n=== Downloading GSE43378 ===")

gse43378 <- getGEO("GSE43378", destdir = geo_dir, GSEMatrix = TRUE, getGPL = TRUE)
gse43378 <- gse43378[[1]]

expr_43378 <- exprs(gse43378)

if (max(expr_43378, na.rm = TRUE) > 50) {
  message("  Log2 transforming GSE43378...")
  expr_43378 <- log2(expr_43378 + 1)
}

# 探针注释
feat_43378 <- fData(gse43378)
gene_col_43378 <- grep("Gene.Symbol|gene_assignment|Symbol", colnames(feat_43378), value = TRUE)[1]

if (!is.na(gene_col_43378)) {
  gene_symbols_43378 <- feat_43378[[gene_col_43378]]
  gene_symbols_43378 <- gsub(" ///.*", "", gene_symbols_43378)

  valid_idx2 <- !is.na(gene_symbols_43378) & gene_symbols_43378 != "" & gene_symbols_43378 != "---"
  expr_43378_valid <- expr_43378[valid_idx2, ]
  genes_valid2 <- gene_symbols_43378[valid_idx2]

  iqr_vals2 <- apply(expr_43378_valid, 1, IQR, na.rm = TRUE)
  probe_df2 <- data.frame(
    probe = rownames(expr_43378_valid),
    gene = genes_valid2,
    iqr = iqr_vals2,
    stringsAsFactors = FALSE
  )
  probe_df2 <- probe_df2 %>%
    dplyr::group_by(gene) %>%
    dplyr::slice_max(iqr, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup()

  expr_43378_gene <- expr_43378_valid[probe_df2$probe, ]
  rownames(expr_43378_gene) <- probe_df2$gene

  message(sprintf("  GSE43378: %d genes x %d samples", nrow(expr_43378_gene), ncol(expr_43378_gene)))
}

clin_43378 <- pData(gse43378)
message("  GSE43378 clinical columns: ", paste(head(colnames(clin_43378), 20), collapse = ", "))

# =============================================================================
# 3. 保存GEO数据
# =============================================================================

save(
  expr_16011_gene, clin_16011, clin_16011_clean,
  file = file.path(DATA_PROC, "gse16011_data.RData")
)

save(
  expr_43378_gene, clin_43378,
  file = file.path(DATA_PROC, "gse43378_data.RData")
)

message("\n=== GEO data download completed ===")
message("GSE16011: ", file.path(DATA_PROC, "gse16011_data.RData"))
message("GSE43378: ", file.path(DATA_PROC, "gse43378_data.RData"))
