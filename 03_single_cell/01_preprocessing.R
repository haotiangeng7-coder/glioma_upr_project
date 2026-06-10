###############################################################################
# 01_preprocessing.R
# 单细胞RNA-seq数据预处理 — 完整版（含评审修订要求）
# GSE131928 (Neftel et al. 2019 Cell) — Smart-seq2
# 支持10X数据（GSE182109等）的SoupX/scDblFinder流程
#
# 输出:
#   - data/processed/seu_preprocessed.rds
#   - figures/FigS1_batch_correction.pdf
#   - figures/FigS2_cellcycle_vs_upr.pdf
#   - figures/sc_qc_*.pdf, sc_elbow_plot.pdf, sc_clustree.pdf
#   - results/preprocessing_summary.csv
###############################################################################

# === 加载配置和依赖 ===
source("00_setup/config.R")
load(file.path(DATA_PROC, "upr_gene_sets.RData"))

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(tidyr)
  library(harmony)
  library(clustree)
  library(lisi)           # compute_lisi
  library(scDblFinder)
  library(SingleCellExperiment)
  library(SoupX)
  if (requireNamespace("celda", quietly = TRUE)) library(celda) # decontX (optional, for 10X data)
  library(scran)
})

set.seed(SEED)

# SCTransform需要更大的全局变量限制
options(future.globals.maxSize = 4 * 1024^3)  # 4GB

# === 辅助函数 ===
safe_source_genes <- function(genes, available) {
  found <- genes[genes %in% available]
  message(sprintf("  Gene overlap: %d / %d", length(found), length(genes)))
  found
}

###############################################################################
# 1. 数据读取
###############################################################################
message("=== [1/9] Loading scRNA-seq data ===")

scrna_dir <- file.path(DATA_RAW, "scRNA", "GSE131928")
is_10x <- FALSE  # 将根据数据格式自动判断

# --- 检测数据格式 ---
mtx_files <- list.files(scrna_dir, pattern = "matrix\\.mtx", recursive = TRUE, full.names = TRUE)
h5_files  <- list.files(scrna_dir, pattern = "\\.h5$", recursive = TRUE, full.names = TRUE)
txt_files <- list.files(scrna_dir, pattern = "\\.(txt|csv|tsv)(|\\.gz)$", recursive = TRUE, full.names = TRUE)

if (length(mtx_files) > 0) {
  # ---------- 10X Genomics格式 ----------
  is_10x <- TRUE
  message("Detected 10X format data.")

  # 如果有多个样本目录，逐个处理
  sample_dirs <- unique(dirname(mtx_files))
  seu_list <- list()

  for (sd in sample_dirs) {
    sample_name <- basename(dirname(sd))
    if (sample_name == basename(scrna_dir)) sample_name <- basename(sd)
    message("  Loading sample: ", sample_name)

    counts_filt <- Read10X(data.dir = sd)

    # SoupX ambient RNA校正（仅10X数据）
    raw_sd <- gsub("filtered", "raw", sd)
    if (dir.exists(raw_sd) && length(list.files(raw_sd, pattern = "matrix\\.mtx")) > 0) {
      message("    Running SoupX ambient RNA correction...")
      tryCatch({
        counts_raw <- Read10X(data.dir = raw_sd)
        # SoupX需要先做快速聚类
        tmp_seu <- CreateSeuratObject(counts = counts_filt, min.cells = 3,
                                       min.features = SC_MIN_FEATURES)
        tmp_seu <- NormalizeData(tmp_seu, verbose = FALSE)
        tmp_seu <- FindVariableFeatures(tmp_seu, verbose = FALSE)
        tmp_seu <- ScaleData(tmp_seu, verbose = FALSE)
        set.seed(SEED)
        tmp_seu <- RunPCA(tmp_seu, npcs = 20, verbose = FALSE)
        set.seed(SEED)
        tmp_seu <- RunUMAP(tmp_seu, dims = 1:20, verbose = FALSE)
        tmp_seu <- FindNeighbors(tmp_seu, dims = 1:20, verbose = FALSE)
        set.seed(SEED)
        tmp_seu <- FindClusters(tmp_seu, resolution = 0.5, verbose = FALSE)

        sc <- SoupChannel(tod = counts_raw, toc = counts_filt)
        sc <- setClusters(sc, setNames(as.character(Idents(tmp_seu)), colnames(tmp_seu)))
        sc <- setDR(sc, Embeddings(tmp_seu, "umap")[, 1:2])
        set.seed(SEED)
        sc <- autoEstCont(sc, verbose = FALSE)
        counts_corrected <- adjustCounts(sc, verbose = FALSE)
        message("    SoupX correction completed. Contamination fraction: ",
                round(sc$fit$rhoEst, 4))
        rm(tmp_seu, sc)
        counts_filt <- counts_corrected
      }, error = function(e) {
        message("    SoupX failed: ", e$message)
        message("    Falling back to DecontX...")
        tryCatch({
          sce_tmp <- SingleCellExperiment(assays = list(counts = counts_filt))
          set.seed(SEED)
          sce_tmp <- decontX(sce_tmp)
          counts_filt <<- decontXcounts(sce_tmp)
          message("    DecontX correction completed. Mean contamination: ",
                  round(mean(sce_tmp$decontX_contamination), 4))
          rm(sce_tmp)
        }, error = function(e2) {
          message("    DecontX also failed: ", e2$message,
                  ". Proceeding without ambient correction.")
        })
      })
    } else {
      message("    No raw matrix found. Attempting DecontX for ambient correction...")
      tryCatch({
        sce_tmp <- SingleCellExperiment(assays = list(counts = counts_filt))
        set.seed(SEED)
        sce_tmp <- decontX(sce_tmp)
        counts_filt <- decontXcounts(sce_tmp)
        message("    DecontX correction completed. Mean contamination: ",
                round(mean(sce_tmp$decontX_contamination), 4))
        rm(sce_tmp)
      }, error = function(e) {
        message("    DecontX failed: ", e$message,
                ". Proceeding without ambient correction.")
      })
    }

    tmp_seu <- CreateSeuratObject(counts = counts_filt, project = sample_name,
                                   min.cells = 3, min.features = SC_MIN_FEATURES)
    tmp_seu$sample <- sample_name
    seu_list[[sample_name]] <- tmp_seu
  }

  # 合并多样本
  if (length(seu_list) > 1) {
    seu <- merge(seu_list[[1]], y = seu_list[-1], add.cell.ids = names(seu_list))
  } else {
    seu <- seu_list[[1]]
  }
  rm(seu_list)

} else if (length(txt_files) > 0) {
  # ---------- 文本/CSV格式（Smart-seq2等）----------
  is_10x <- FALSE
  count_file <- txt_files[grep("count|expression|TPM|GSE|Smartseq|10X", txt_files, ignore.case = TRUE)]
  if (length(count_file) == 0) count_file <- txt_files

  # 如果有多个文件（如Smart-seq2 + 10X），合并处理
  if (length(count_file) > 1) {
    message("Found multiple expression files — merging:")
    mat_list <- list()
    platform_labels <- c()
    for (cf in count_file) {
      fname <- basename(cf)
      message("  Loading: ", fname)
      tmp <- data.table::fread(cf, data.table = FALSE)
      rownames(tmp) <- tmp[, 1]
      tmp <- tmp[, -1]
      # 检测平台
      platform <- ifelse(grepl("Smartseq|Smart", fname, ignore.case = TRUE), "Smart-seq2",
                   ifelse(grepl("10X|10x", fname, ignore.case = TRUE), "10X", "Unknown"))
      colnames(tmp) <- paste0(colnames(tmp), "_", platform)
      platform_labels <- c(platform_labels, setNames(rep(platform, ncol(tmp)), colnames(tmp)))
      mat_list[[fname]] <- as.matrix(tmp)
      message(sprintf("    %s: %d genes x %d cells", platform, nrow(tmp), ncol(tmp)))
    }
    # 取共同基因合并
    common_genes <- Reduce(intersect, lapply(mat_list, rownames))
    counts <- do.call(cbind, lapply(mat_list, function(m) m[common_genes, ]))
    message(sprintf("  Merged: %d genes x %d cells", nrow(counts), ncol(counts)))
    rm(mat_list)
  } else {
    message("Loading text format data: ", basename(count_file[1]))
    counts <- data.table::fread(count_file[1], data.table = FALSE)
    rownames(counts) <- counts[, 1]
    counts <- counts[, -1]
    counts <- as.matrix(counts)
    platform_labels <- setNames(rep("Unknown", ncol(counts)), colnames(counts))
  }

  # TPM数据处理：若值有小数且非整数，说明是TPM/FPKM而非raw counts
  if (any(counts %% 1 != 0, na.rm = TRUE)) {
    message("  Data appears to be TPM/FPKM (non-integer values detected).")
    message("  Converting to pseudo-counts: round(TPM) for Seurat compatibility.")
    message("  Note: SCTransform will use log-normalized TPM via NormalizeData as fallback.")
    counts <- round(counts)  # 转为伪counts供CreateSeuratObject
  }

  # Smart-seq2不做SoupX（非液滴方法无ambient RNA问题）
  message("  Plate-based/processed data detected. Skipping SoupX/DecontX (not applicable).")

  seu <- CreateSeuratObject(counts = counts, project = "Glioma_scRNA",
                             min.cells = 3, min.features = SC_MIN_FEATURES)

  # 添加平台信息
  if (exists("platform_labels")) {
    seu$platform <- platform_labels[colnames(seu)]
    seu$sample <- gsub("_Smart-seq2$|_10X$", "", colnames(seu))
    # 从cell名提取样本（如MGH101-P1-A04 → MGH101）
    seu$sample <- gsub("^(MGH[0-9]+|BT[0-9]+).*", "\\1", seu$sample)
  }

  # 尝试读取元数据
  meta_files <- list.files(scrna_dir, pattern = "meta|annotation|cell.?type|sample",
                            recursive = TRUE, full.names = TRUE, ignore.case = TRUE)
  if (length(meta_files) > 0) {
    message("Loading metadata: ", basename(meta_files[1]))
    meta <- data.table::fread(meta_files[1], data.table = FALSE)
    if (any(meta[, 1] %in% colnames(seu))) {
      rownames(meta) <- meta[, 1]
      common_cells <- intersect(colnames(seu), rownames(meta))
      seu <- seu[, common_cells]
      for (col in colnames(meta)[-1]) {
        seu[[col]] <- meta[colnames(seu), col]
      }
      message("  Metadata added for ", length(common_cells), " cells")
    }
  }

  # 尝试识别样本列
  sample_cols <- grep("sample|patient|donor|subject|orig", colnames(seu@meta.data),
                       ignore.case = TRUE, value = TRUE)
  if (length(sample_cols) > 0 && !"sample" %in% colnames(seu@meta.data)) {
    seu$sample <- seu@meta.data[[sample_cols[1]]]
  } else if (!"sample" %in% colnames(seu@meta.data)) {
    seu$sample <- "GSE131928"
  }

} else if (length(h5_files) > 0) {
  # ---------- h5格式 ----------
  is_10x <- TRUE
  message("Loading h5 format data: ", basename(h5_files[1]))
  counts <- Read10X_h5(h5_files[1])
  seu <- CreateSeuratObject(counts = counts, project = "Glioma_scRNA",
                             min.cells = 3, min.features = SC_MIN_FEATURES)
  seu$sample <- "GSE131928"
} else {
  stop("No recognizable data files found in ", scrna_dir,
       "\nExpected: matrix.mtx (10X), .txt/.csv/.tsv (Smart-seq2), or .h5 files.",
       "\nPlease download data first.")
}

message(sprintf("Raw data loaded: %d genes x %d cells", nrow(seu), ncol(seu)))
message(sprintf("  Samples detected: %s", paste(unique(seu$sample), collapse = ", ")))
message(sprintf("  Data type: %s", ifelse(is_10x, "10X Chromium (droplet)",
                                           "Smart-seq2/plate-based")))

###############################################################################
# 2. Doublet检测（仅10X数据，按样本独立）
###############################################################################
message("=== [2/9] Doublet Detection ===")

n_before_doublet <- ncol(seu)

if (is_10x) {
  message("Running scDblFinder on 10X data (per sample)...")
  sce <- as.SingleCellExperiment(seu)

  samples <- unique(seu$sample)
  doublet_calls  <- character(ncol(sce))
  names(doublet_calls)  <- colnames(sce)
  doublet_scores <- numeric(ncol(sce))
  names(doublet_scores) <- colnames(sce)

  for (s in samples) {
    idx <- which(seu$sample == s)
    message(sprintf("  Sample '%s': %d cells", s, length(idx)))
    set.seed(SEED)
    sce_sub <- scDblFinder(sce[, idx])
    doublet_calls[colnames(sce_sub)]  <- sce_sub$scDblFinder.class
    doublet_scores[colnames(sce_sub)] <- sce_sub$scDblFinder.score
  }

  seu$scDblFinder_class <- doublet_calls[colnames(seu)]
  seu$scDblFinder_score <- doublet_scores[colnames(seu)]

  n_doublets <- sum(seu$scDblFinder_class == "doublet", na.rm = TRUE)
  message(sprintf("  Total doublets detected: %d (%.1f%%)",
                  n_doublets, n_doublets / ncol(seu) * 100))

  # 过滤doublets
  seu <- subset(seu, scDblFinder_class == "singlet")
  message(sprintf("  After doublet removal: %d cells", ncol(seu)))
  rm(sce)
} else {
  message("Smart-seq2/plate-based data: scDblFinder not applicable. Skipping.")
  seu$scDblFinder_class <- "singlet"
  seu$scDblFinder_score <- 0
}

###############################################################################
# 3. QC过滤
###############################################################################
message("=== [3/9] Quality Control ===")

# 计算QC指标
seu[["percent.mt"]]   <- PercentageFeatureSet(seu, pattern = "^MT-")
seu[["percent.ribo"]] <- PercentageFeatureSet(seu, pattern = "^RPS|^RPL")

# QC可视化（过滤前）— 使用ggplot2避免Seurat v5 S4兼容性问题
qc_df <- seu@meta.data[, c("nFeature_RNA", "nCount_RNA", "percent.mt")]
qc_df$cell <- rownames(qc_df)
qc_long <- tidyr::pivot_longer(qc_df, cols = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                                names_to = "metric", values_to = "value")
p_qc <- ggplot(qc_long, aes(x = metric, y = value, fill = metric)) +
  geom_violin(scale = "width", trim = FALSE) +
  facet_wrap(~metric, scales = "free", ncol = 3) +
  theme_bw() + theme(legend.position = "none") +
  labs(title = "QC Metrics (Before Filtering)", x = "", y = "Value")
ggsave(file.path(FIG_DIR, "sc_qc_violin_before.pdf"), p_qc, width = 12, height = 5)

p_sc1 <- ggplot(qc_df, aes(x = nCount_RNA, y = percent.mt)) +
  geom_point(size = 0.1, alpha = 0.3) + theme_bw() + labs(title = "nCount vs MT%")
p_sc2 <- ggplot(qc_df, aes(x = nCount_RNA, y = nFeature_RNA)) +
  geom_point(size = 0.1, alpha = 0.3) + theme_bw() + labs(title = "nCount vs nFeature")
ggsave(file.path(FIG_DIR, "sc_qc_scatter.pdf"), p_sc1 + p_sc2, width = 12, height = 5)

# 记录过滤前数量
n_before_qc <- ncol(seu)

# 过滤
seu <- subset(seu,
              nFeature_RNA > SC_MIN_FEATURES &
                nFeature_RNA < SC_MAX_FEATURES &
                percent.mt < SC_MAX_MT_PCT)

n_after_qc <- ncol(seu)
message(sprintf("QC filtering: %d -> %d cells (removed %d, %.1f%%)",
                n_before_qc, n_after_qc, n_before_qc - n_after_qc,
                (n_before_qc - n_after_qc) / n_before_qc * 100))

# QC可视化（过滤后）
qc_df2 <- seu@meta.data[, c("nFeature_RNA", "nCount_RNA", "percent.mt")]
qc_long2 <- tidyr::pivot_longer(qc_df2, cols = everything(), names_to = "metric", values_to = "value")
p_qc2 <- ggplot(qc_long2, aes(x = metric, y = value, fill = metric)) +
  geom_violin(scale = "width", trim = FALSE) +
  facet_wrap(~metric, scales = "free", ncol = 3) +
  theme_bw() + theme(legend.position = "none") +
  labs(title = "QC Metrics (After Filtering)", x = "", y = "Value")
ggsave(file.path(FIG_DIR, "sc_qc_violin_after.pdf"), p_qc2, width = 12, height = 5)

###############################################################################
# 4. SCTransform v2 标准化
###############################################################################
message("=== [4/9] SCTransform v2 Normalization ===")

set.seed(SEED)
tryCatch({
  seu <- SCTransform(seu,
                      vst.flavor          = SCT_VST_FLAVOR,
                      variable.features.n = SCT_NFEATURES,
                      verbose             = TRUE)
  message(sprintf("SCTransform v2 completed: %d variable features selected",
                  length(VariableFeatures(seu))))
}, error = function(e) {
  message("SCTransform failed (likely TPM/non-count data): ", e$message)
  message("Falling back to standard NormalizeData + ScaleData pipeline.")
  seu <<- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000)
  seu <<- FindVariableFeatures(seu, selection.method = "vst", nfeatures = SCT_NFEATURES)
  seu <<- ScaleData(seu, features = rownames(seu))
  message(sprintf("NormalizeData pipeline completed: %d variable features selected",
                  length(VariableFeatures(seu))))
})

###############################################################################
# 5. 细胞周期评估与条件回归
###############################################################################
message("=== [5/9] Cell Cycle Scoring & UPR Correlation Assessment ===")

# Seurat自带的细胞周期基因
s.genes   <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

set.seed(SEED)
seu <- CellCycleScoring(seu,
                         s.features   = s.genes,
                         g2m.features = g2m.genes,
                         set.ident    = FALSE)

message("Cell cycle phase distribution:")
cc_table <- table(seu$Phase)
print(cc_table)

# 计算临时UPR评分用于评估细胞周期相关性
upr_genes_valid <- safe_source_genes(UPR_broad_genes, rownames(seu))

regressed_cc <- FALSE
cc_upr_cor <- NA

if (length(upr_genes_valid) >= 5) {
  set.seed(SEED)
  seu <- AddModuleScore(seu, features = list(upr_genes_valid), name = "UPR_temp_score")

  # 计算细胞周期综合分数与UPR临时分数的相关性
  cc_score <- seu$S.Score + seu$G2M.Score
  upr_temp <- seu$UPR_temp_score1

  cc_upr_cor  <- cor(cc_score, upr_temp, method = "spearman", use = "complete.obs")
  cc_upr_pval <- cor.test(cc_score, upr_temp, method = "spearman")$p.value

  message(sprintf("Cell cycle vs UPR correlation: rho = %.3f, p = %.2e",
                  cc_upr_cor, cc_upr_pval))
  message(sprintf("Threshold for regression: |r| > %.1f", CC_COR_THRESHOLD))

  # === Supplementary Figure S2: 细胞周期 vs UPR 相关性 ===
  cc_upr_df <- data.frame(
    CellCycle_Score = cc_score,
    UPR_Score       = upr_temp,
    Phase           = seu$Phase,
    S_Score         = seu$S.Score,
    G2M_Score       = seu$G2M.Score
  )

  cor_s   <- cor(seu$S.Score, upr_temp, method = "spearman", use = "complete.obs")
  cor_g2m <- cor(seu$G2M.Score, upr_temp, method = "spearman", use = "complete.obs")

  p_s2a <- ggplot(cc_upr_df, aes(x = CellCycle_Score, y = UPR_Score, color = Phase)) +
    geom_point(alpha = 0.3, size = 0.5) +
    geom_smooth(method = "lm", color = "black", linewidth = 0.8, se = TRUE) +
    scale_color_manual(values = c("G1" = "#00A087", "S" = "#E64B35", "G2M" = "#3C5488")) +
    labs(x = "Cell Cycle Score (S + G2M)",
         y = "UPR Module Score",
         title = sprintf("Cell Cycle vs UPR (rho=%.3f, p=%.1e)", cc_upr_cor, cc_upr_pval)) +
    THEME_PUBLICATION

  p_s2b <- ggplot(cc_upr_df, aes(x = S_Score, y = UPR_Score)) +
    geom_point(alpha = 0.3, size = 0.5, color = "#E64B35") +
    geom_smooth(method = "lm", color = "black", linewidth = 0.8) +
    labs(x = "S Phase Score", y = "UPR Module Score",
         title = sprintf("S Phase vs UPR (rho=%.3f)", cor_s)) +
    THEME_PUBLICATION

  p_s2c <- ggplot(cc_upr_df, aes(x = G2M_Score, y = UPR_Score)) +
    geom_point(alpha = 0.3, size = 0.5, color = "#3C5488") +
    geom_smooth(method = "lm", color = "black", linewidth = 0.8) +
    labs(x = "G2M Phase Score", y = "UPR Module Score",
         title = sprintf("G2M Phase vs UPR (rho=%.3f)", cor_g2m)) +
    THEME_PUBLICATION

  p_s2d <- ggplot(cc_upr_df, aes(x = Phase, y = UPR_Score, fill = Phase)) +
    geom_violin(alpha = 0.7, scale = "width") +
    geom_boxplot(width = 0.1, fill = "white", outlier.size = 0.3) +
    scale_fill_manual(values = c("G1" = "#00A087", "S" = "#E64B35", "G2M" = "#3C5488")) +
    labs(x = "Cell Cycle Phase", y = "UPR Module Score",
         title = "UPR Score by Cell Cycle Phase") +
    THEME_PUBLICATION +
    theme(legend.position = "none")

  p_figS2 <- (p_s2a | p_s2b) / (p_s2c | p_s2d) +
    plot_annotation(
      title = "Supplementary Figure S2: Cell Cycle vs UPR Correlation",
      subtitle = sprintf(
        "Decision: %s cell cycle regression (|rho|=%.3f, threshold=%.1f)",
        ifelse(abs(cc_upr_cor) > CC_COR_THRESHOLD, "REGRESS", "SKIP"),
        abs(cc_upr_cor), CC_COR_THRESHOLD
      ),
      theme = theme(
        plot.title    = element_text(size = 16, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5)
      )
    )

  ggsave(file.path(FIG_DIR, "FigS2_cellcycle_vs_upr.pdf"), p_figS2, width = 14, height = 10)

  # 条件回归
  if (abs(cc_upr_cor) > CC_COR_THRESHOLD) {
    message("*** Cell cycle significantly correlated with UPR. Regressing out S.Score and G2M.Score ***")
    set.seed(SEED)
    seu <- SCTransform(seu,
                        vst.flavor          = SCT_VST_FLAVOR,
                        variable.features.n = SCT_NFEATURES,
                        vars.to.regress     = c("S.Score", "G2M.Score"),
                        verbose             = TRUE)
    regressed_cc <- TRUE
    message("SCTransform re-run with cell cycle regression.")
  } else {
    message("Cell cycle NOT significantly correlated with UPR. No regression needed.")
  }

  # 清理临时列
  seu$UPR_temp_score1 <- NULL
} else {
  message("WARNING: Too few UPR genes found in data (<5). Skipping cell cycle correlation check.")
}

###############################################################################
# 6. PCA降维
###############################################################################
message("=== [6/9] PCA ===")

set.seed(SEED)
seu <- RunPCA(seu, npcs = 50, verbose = FALSE)

# Elbow plot
pca_stdev <- Stdev(seu, reduction = "pca")

p_elbow <- ElbowPlot(seu, ndims = 50) +
  THEME_PUBLICATION +
  geom_vline(xintercept = SC_N_PCS, linetype = "dashed", color = "red") +
  annotate("text", x = SC_N_PCS + 1, y = max(pca_stdev[1:50]) * 0.5,
           label = paste0("PC = ", SC_N_PCS), color = "red", hjust = 0) +
  ggtitle("Elbow Plot")

ggsave(file.path(FIG_DIR, "sc_elbow_plot.pdf"), p_elbow, width = 7, height = 5)

###############################################################################
# 7. 批次校正 — Harmony + LISI/kBET评估
###############################################################################
message("=== [7/9] Batch Correction (Harmony) ===")

n_samples <- length(unique(seu$sample))
batch_var <- "sample"

if (n_samples > 1) {
  message(sprintf("Multiple samples detected (%d). Running Harmony...", n_samples))

  # --- UMAP before Harmony（用于对比）---
  set.seed(SEED)
  seu <- RunUMAP(seu, reduction = "pca", dims = 1:SC_N_PCS,
                  reduction.name = "umap_uncorrected")

  # --- Harmony batch correction ---
  set.seed(SEED)
  # Harmony批次校正
  seu <- RunHarmony(seu,
                     group.by.vars     = batch_var,
                     reduction.use     = "pca",
                     dims.use          = 1:SC_N_PCS,
                     theta             = HARMONY_THETA,
                     lambda            = HARMONY_LAMBDA,
                     max.iter.harmony  = HARMONY_MAXITER,
                     verbose           = TRUE)

  reduction_use <- "harmony"

  # --- UMAP after Harmony ---
  set.seed(SEED)
  seu <- RunUMAP(seu, reduction = "harmony", dims = 1:SC_N_PCS, reduction.name = "umap")

  # === 批次校正效果评估: LISI ===
  message("Computing LISI scores...")
  lisi_before <- compute_lisi(
    Embeddings(seu, "umap_uncorrected")[, 1:2],
    seu@meta.data, batch_var
  )
  lisi_after <- compute_lisi(
    Embeddings(seu, "umap")[, 1:2],
    seu@meta.data, batch_var
  )

  seu$LISI_before <- lisi_before[[batch_var]]
  seu$LISI_after  <- lisi_after[[batch_var]]

  mean_lisi_before <- mean(lisi_before[[batch_var]], na.rm = TRUE)
  mean_lisi_after  <- mean(lisi_after[[batch_var]], na.rm = TRUE)
  message(sprintf("  LISI (before): %.3f", mean_lisi_before))
  message(sprintf("  LISI (after):  %.3f", mean_lisi_after))
  message(sprintf("  Ideal LISI for %d batches: %.1f", n_samples, n_samples))

  # === 批次校正效果评估: kBET ===
  message("Computing kBET acceptance rate...")
  kbet_before_rate <- NA
  kbet_after_rate  <- NA

  tryCatch({
    library(kBET)

    # 采样子集以加速（最多5000细胞）
    n_sub <- min(5000, ncol(seu))
    set.seed(SEED)
    sub_idx <- sample(ncol(seu), n_sub)

    set.seed(SEED)
    kbet_before <- kBET(
      df      = Embeddings(seu, "pca")[sub_idx, 1:SC_N_PCS],
      batch   = seu@meta.data[[batch_var]][sub_idx],
      plot    = FALSE,
      do.pca  = FALSE,
      verbose = FALSE
    )
    kbet_before_rate <- 1 - kbet_before$summary$kBET.observed[1]

    set.seed(SEED)
    kbet_after <- kBET(
      df      = Embeddings(seu, "harmony")[sub_idx, 1:SC_N_PCS],
      batch   = seu@meta.data[[batch_var]][sub_idx],
      plot    = FALSE,
      do.pca  = FALSE,
      verbose = FALSE
    )
    kbet_after_rate <- 1 - kbet_after$summary$kBET.observed[1]

    message(sprintf("  kBET acceptance (before): %.3f", kbet_before_rate))
    message(sprintf("  kBET acceptance (after):  %.3f", kbet_after_rate))
  }, error = function(e) {
    message("  kBET computation failed: ", e$message)
    message("  Proceeding with LISI results only.")
  })

  # === Supplementary Figure S1: 批次校正前后UMAP + LISI/kBET ===
  message("Generating Supplementary Figure S1...")

  p_s1a <- DimPlot(seu, reduction = "umap_uncorrected", group.by = batch_var,
                    pt.size = 0.3, raster = FALSE) +
    THEME_PUBLICATION + ggtitle("Before Harmony")

  p_s1b <- DimPlot(seu, reduction = "umap", group.by = batch_var,
                    pt.size = 0.3, raster = FALSE) +
    THEME_PUBLICATION + ggtitle("After Harmony")

  lisi_df <- data.frame(
    Condition = rep(c("Before Harmony", "After Harmony"), each = ncol(seu)),
    LISI      = c(seu$LISI_before, seu$LISI_after)
  )
  lisi_df$Condition <- factor(lisi_df$Condition,
                               levels = c("Before Harmony", "After Harmony"))

  p_s1c <- ggplot(lisi_df, aes(x = LISI, fill = Condition)) +
    geom_density(alpha = 0.5) +
    geom_vline(xintercept = c(mean_lisi_before, mean_lisi_after),
               linetype = "dashed", color = c("#E64B35", "#4DBBD5")) +
    scale_fill_manual(values = c("Before Harmony" = "#E64B35",
                                  "After Harmony"  = "#4DBBD5")) +
    labs(x = "LISI Score", y = "Density",
         title = sprintf("LISI Distribution (Before: %.2f, After: %.2f)",
                         mean_lisi_before, mean_lisi_after)) +
    THEME_PUBLICATION

  stats_text <- sprintf(
    paste0("Batch Correction Summary\n\n",
           "Method: Harmony (theta=%.0f, lambda=%.0f)\n",
           "Samples: %d\n\n",
           "LISI (mean):\n  Before: %.3f\n  After:  %.3f\n  Ideal:  %.1f\n\n",
           "kBET acceptance:\n  Before: %s\n  After:  %s"),
    HARMONY_THETA, HARMONY_LAMBDA, n_samples,
    mean_lisi_before, mean_lisi_after, n_samples,
    ifelse(is.na(kbet_before_rate), "N/A", sprintf("%.3f", kbet_before_rate)),
    ifelse(is.na(kbet_after_rate), "N/A", sprintf("%.3f", kbet_after_rate))
  )

  p_s1d <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = stats_text, size = 4,
             hjust = 0.5, vjust = 0.5, family = "mono") +
    theme_void() +
    ggtitle("Integration Metrics")

  p_figS1 <- (p_s1a | p_s1b) / (p_s1c | p_s1d) +
    plot_annotation(
      title = "Supplementary Figure S1: Batch Correction Evaluation",
      theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))
    )

  ggsave(file.path(FIG_DIR, "FigS1_batch_correction.pdf"), p_figS1, width = 16, height = 12)

  # 保存批次校正指标
  batch_metrics <- data.frame(
    Metric = c("LISI_before", "LISI_after",
               "kBET_acceptance_before", "kBET_acceptance_after",
               "n_samples", "harmony_theta", "harmony_lambda"),
    Value  = c(mean_lisi_before, mean_lisi_after,
               kbet_before_rate, kbet_after_rate,
               n_samples, HARMONY_THETA, HARMONY_LAMBDA)
  )
  write.csv(batch_metrics, file.path(RES_DIR, "batch_correction_metrics.csv"), row.names = FALSE)

} else {
  message("Single sample detected. Skipping batch correction.")
  reduction_use <- "pca"

  set.seed(SEED)
  seu <- RunUMAP(seu, reduction = "pca", dims = 1:SC_N_PCS)

  # 单样本也输出FigS1
  p_figS1 <- DimPlot(seu, reduction = "umap", group.by = "sample", pt.size = 0.3) +
    THEME_PUBLICATION +
    ggtitle("Supplementary Figure S1: Single Sample (No Batch Correction Needed)")

  ggsave(file.path(FIG_DIR, "FigS1_batch_correction.pdf"), p_figS1, width = 10, height = 8)
}

###############################################################################
# 8. 聚类 — clustree分辨率选择
###############################################################################
message("=== [8/9] Clustering with Resolution Selection (clustree) ===")

seu <- FindNeighbors(seu, reduction = reduction_use, dims = 1:SC_N_PCS)

# 测试多个分辨率
resolutions <- seq(0.2, 1.6, by = 0.2)
for (res in resolutions) {
  set.seed(SEED)
  seu <- FindClusters(seu, resolution = res, verbose = FALSE)
}

# Clustree可视化
p_clustree <- clustree(seu, prefix = "SCT_snn_res.") +
  THEME_PUBLICATION +
  ggtitle("Cluster Stability Across Resolutions")

ggsave(file.path(FIG_DIR, "sc_clustree.pdf"), p_clustree, width = 12, height = 10)

# 使用config中定义的分辨率作为最终分辨率
set.seed(SEED)
seu <- FindClusters(seu, resolution = SC_RESOLUTION, verbose = FALSE)
Idents(seu) <- paste0("SCT_snn_res.", SC_RESOLUTION)

message(sprintf("Final clustering: resolution=%.1f, %d clusters",
                SC_RESOLUTION, length(unique(Idents(seu)))))

# UMAP聚类可视化
p_umap_cluster <- DimPlot(seu, reduction = "umap", label = TRUE, label.size = 4,
                            pt.size = 0.3, raster = FALSE) +
  THEME_PUBLICATION +
  ggtitle(sprintf("Clustering (resolution=%.1f)", SC_RESOLUTION))

ggsave(file.path(FIG_DIR, "sc_umap_clusters.pdf"), p_umap_cluster, width = 8, height = 7)

# 按样本着色
p_umap_sample <- DimPlot(seu, reduction = "umap", group.by = "sample",
                           pt.size = 0.3, raster = FALSE) +
  THEME_PUBLICATION +
  ggtitle("Cells by Sample")

ggsave(file.path(FIG_DIR, "sc_umap_sample.pdf"), p_umap_sample, width = 8, height = 7)

###############################################################################
# 9. 保存结果
###############################################################################
message("=== [9/9] Saving Results ===")

# 预处理汇总
summary_df <- data.frame(
  Step = c(
    "Raw cells loaded",
    "After doublet removal (10X only)",
    "After QC filtering",
    "Normalization",
    "Cell cycle regression",
    "PCA dims used",
    "Batch correction",
    "Final resolution",
    "Final clusters",
    "Final cells",
    "Final genes"
  ),
  Value = c(
    as.character(n_before_doublet),
    ifelse(is_10x,
           as.character(n_before_qc),
           "N/A (Smart-seq2)"),
    as.character(n_after_qc),
    sprintf("SCTransform v2 (%d features)", SCT_NFEATURES),
    ifelse(regressed_cc,
           sprintf("Yes (|rho|=%.3f > %.1f)", abs(cc_upr_cor), CC_COR_THRESHOLD),
           sprintf("No (|rho|=%.3f <= %.1f)", abs(ifelse(is.na(cc_upr_cor), 0, cc_upr_cor)),
                   CC_COR_THRESHOLD)),
    as.character(SC_N_PCS),
    ifelse(n_samples > 1,
           sprintf("Harmony (theta=%s, lambda=%s)", HARMONY_THETA, HARMONY_LAMBDA),
           "None (single sample)"),
    as.character(SC_RESOLUTION),
    as.character(length(unique(Idents(seu)))),
    as.character(ncol(seu)),
    as.character(nrow(seu))
  )
)

write.csv(summary_df, file.path(RES_DIR, "preprocessing_summary.csv"), row.names = FALSE)

# 保存Seurat对象
saveRDS(seu, file = file.path(DATA_PROC, "seu_preprocessed.rds"))

message("\n=== Preprocessing Complete ===")
message("Seurat object:    ", file.path(DATA_PROC, "seu_preprocessed.rds"))
message("Summary CSV:      ", file.path(RES_DIR, "preprocessing_summary.csv"))
message("Figure S1:        ", file.path(FIG_DIR, "FigS1_batch_correction.pdf"))
message("Figure S2:        ", file.path(FIG_DIR, "FigS2_cellcycle_vs_upr.pdf"))
message(sprintf("Final dimensions: %d genes x %d cells in %d clusters",
                nrow(seu), ncol(seu), length(unique(Idents(seu)))))
message("\nNext step: Run 03_single_cell/02_cell_annotation.R")
