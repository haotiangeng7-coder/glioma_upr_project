#!/usr/bin/env Rscript
# GSE182109 Validation Pipeline for UPR patterns from GSE131928 primary analysis
# Author: validation pipeline
# Date: 2026-04-18

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(harmony)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(AUCell)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(viridis)
  library(future)
})

set.seed(42)
options(future.globals.maxSize = 200 * 1024^3)  # 200 GB
plan("multicore", workers = 8)

# ---- Paths ----
PROJECT_DIR <- getwd()
RAW_DIR    <- file.path(PROJECT_DIR, "data", "raw", "GEO", "GSE182109", "raw_counts")
OUT_DIR    <- file.path(PROJECT_DIR, "results", "GSE182109_validation")
FIG_DIR    <- file.path(OUT_DIR, "figures")
UPR_RDATA  <- file.path(PROJECT_DIR, "data", "processed", "upr_gene_sets.RData")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

LOG <- file.path(OUT_DIR, "pipeline.log")
log_msg <- function(...) {
  msg <- paste0("[", format(Sys.time(), "%H:%M:%S"), "] ", paste0(..., collapse = ""))
  cat(msg, "\n")
  cat(msg, "\n", file = LOG, append = TRUE)
}

log_msg("=== GSE182109 Validation Pipeline START ===")

# ---- Load UPR gene sets ----
load(UPR_RDATA)
upr_arms <- list(
  IRE1_XBP1 = upr_gene_list$IRE1_XBP1,
  PERK_ATF4 = upr_gene_list$PERK_ATF4,
  ATF6      = upr_gene_list$ATF6,
  UPR_broad = upr_gene_list$UPR_broad
)
log_msg("UPR arms loaded: ", paste(names(upr_arms), collapse=", "))

# ============================================================
# STEP 1: Per-sample 10X loading
# ============================================================
log_msg("STEP 1: Setting up per-sample directories")

all_files <- list.files(RAW_DIR, pattern = "_(barcodes|features|matrix)")
samples_all <- unique(gsub("_(barcodes|features|matrix).*", "", all_files))
samples <- samples_all[grepl("GBM", samples_all)]   # exclude LGG
log_msg("Total samples: ", length(samples_all), "; GBM samples (ndGBM+rGBM): ", length(samples))

by_sample_dir <- file.path(RAW_DIR, "by_sample")
dir.create(by_sample_dir, showWarnings = FALSE)
for (s in samples) {
  sdir <- file.path(by_sample_dir, s)
  dir.create(sdir, recursive = TRUE, showWarnings = FALSE)
  for (suf in c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")) {
    src <- file.path(RAW_DIR, paste0(s, "_", suf))
    dst <- file.path(sdir, suf)
    if (!file.exists(dst) && file.exists(src)) {
      file.symlink(src, dst)
    }
  }
}

CHECKPOINT_RAW <- file.path(OUT_DIR, "checkpoint_step2_qc.rds")

if (file.exists(CHECKPOINT_RAW)) {
  log_msg("Loading checkpoint after QC...")
  obj <- readRDS(CHECKPOINT_RAW)
} else {
  log_msg("Loading per-sample 10X counts and creating Seurat objects")
  obj_list <- list()
  qc_pre <- data.frame(sample = character(), cells_loaded = integer(),
                       cells_post_filter = integer(), cells_post_doublet = integer(),
                       stringsAsFactors = FALSE)

  for (s in samples) {
    sdir <- file.path(by_sample_dir, s)
    counts <- tryCatch(Read10X(data.dir = sdir),
                       error = function(e) { log_msg("ERROR loading ", s, ": ", e$message); NULL })
    if (is.null(counts)) next
    so <- CreateSeuratObject(counts = counts, project = s,
                             min.cells = 3, min.features = 200)
    sample_id <- sub("^GSM[0-9]+_", "", s)
    so$sample_id <- sample_id
    so$tumor_type <- ifelse(grepl("^ndGBM", sample_id), "ndGBM",
                            ifelse(grepl("^rGBM", sample_id), "rGBM", NA))
    so[["percent.mt"]] <- PercentageFeatureSet(so, pattern = "^MT-")

    n_loaded <- ncol(so)

    # ---- QC filter ----
    so <- subset(so, subset = nFeature_RNA >= 200 & nFeature_RNA <= 6000 &
                              nCount_RNA   >= 500 & nCount_RNA   <= 50000 &
                              percent.mt   <  20)
    n_post_qc <- ncol(so)

    # ---- Doublet detection ----
    n_post_dbl <- n_post_qc
    if (ncol(so) >= 100) {
      sce <- as.SingleCellExperiment(so)
      sce <- tryCatch(scDblFinder(sce, verbose = FALSE),
                      error = function(e) { log_msg("scDblFinder fail ", s, ": ", e$message); NULL })
      if (!is.null(sce)) {
        so$scDblFinder_class <- sce$scDblFinder.class
        so$scDblFinder_score <- sce$scDblFinder.score
        so <- subset(so, scDblFinder_class == "singlet")
        n_post_dbl <- ncol(so)
      }
    }
    qc_pre <- rbind(qc_pre, data.frame(sample = sample_id,
                                       cells_loaded = n_loaded,
                                       cells_post_filter = n_post_qc,
                                       cells_post_doublet = n_post_dbl))
    log_msg(sprintf("  %s: loaded=%d, post-QC=%d, post-dbl=%d",
                    sample_id, n_loaded, n_post_qc, n_post_dbl))
    obj_list[[sample_id]] <- so
  }

  write.csv(qc_pre, file.path(OUT_DIR, "qc_per_sample.csv"), row.names = FALSE)

  # ---- Merge ----
  log_msg("Merging ", length(obj_list), " sample objects")
  obj <- merge(obj_list[[1]], y = obj_list[-1],
               add.cell.ids = names(obj_list), project = "GSE182109")
  log_msg("Merged object: ", ncol(obj), " cells, ", nrow(obj), " genes")

  saveRDS(obj, CHECKPOINT_RAW)
  log_msg("Saved QC checkpoint")
}

# ============================================================
# STEP 3: Normalize, integrate, cluster
# ============================================================
CHECKPOINT_INT <- file.path(OUT_DIR, "checkpoint_step3_integrated.rds")

if (file.exists(CHECKPOINT_INT)) {
  log_msg("Loading integrated checkpoint...")
  obj <- readRDS(CHECKPOINT_INT)
} else {
  log_msg("STEP 3: Normalize + integrate (layers already split by merge)")
  # Seurat v5 merge() auto-splits layers per sample; skip manual split

  log_msg("Running SCTransform per layer (vars.to.regress=percent.mt)")
  obj <- SCTransform(obj, vars.to.regress = "percent.mt",
                     vst.flavor = "v2", verbose = FALSE,
                     conserve.memory = TRUE)

  log_msg("Running PCA (50 comps)")
  obj <- RunPCA(obj, npcs = 50, verbose = FALSE)

  log_msg("Harmony integration")
  obj <- IntegrateLayers(object = obj, method = HarmonyIntegration,
                         orig.reduction = "pca",
                         new.reduction = "harmony",
                         normalization.method = "SCT",
                         verbose = FALSE)

  log_msg("UMAP + clustering on harmony 1:30")
  obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, verbose = FALSE)
  obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.4, verbose = FALSE)

  saveRDS(obj, CHECKPOINT_INT)
  log_msg("Saved integration checkpoint")
}

log_msg("Final cell count: ", ncol(obj))

# ============================================================
# STEP 4: Cell type annotation via marker module scores
# ============================================================
CHECKPOINT_ANN <- file.path(OUT_DIR, "checkpoint_step4_annotated.rds")

if (file.exists(CHECKPOINT_ANN)) {
  log_msg("Loading annotated checkpoint...")
  obj <- readRDS(CHECKPOINT_ANN)
} else {
  log_msg("STEP 4: Marker-based cell type annotation")
  DefaultAssay(obj) <- "SCT"

  marker_sets <- list(
    Malignant       = c("EGFR","CDK4","MDM2","PDGFRA","SOX2","OLIG2"),
    TAM             = c("CD68","CD163","AIF1","CSF1R"),
    Microglia       = c("P2RY12","TMEM119","CX3CR1"),
    CD8_T           = c("CD8A","CD8B","GZMB","NKG7"),
    CD4_T           = c("CD4","IL7R","FOXP3"),
    Endothelial     = c("PECAM1","VWF","CLDN5","CDH5"),
    Oligodendrocyte = c("MBP","MOG","PLP1","OLIG1"),
    Astrocyte       = c("GFAP","AQP4","S100B","ALDH1L1")
  )

  # Filter genes present
  marker_sets <- lapply(marker_sets, function(g) intersect(g, rownames(obj)))
  log_msg("Marker set sizes after filtering: ",
          paste(names(marker_sets), sapply(marker_sets, length), sep="=", collapse=", "))

  for (ct in names(marker_sets)) {
    obj <- AddModuleScore(obj, features = list(marker_sets[[ct]]),
                          name = paste0("score_", ct), assay = "SCT")
  }
  score_cols <- paste0("score_", names(marker_sets), "1")
  meta <- obj@meta.data

  # Per-cluster mean scores -> assign cell type with highest mean
  cluster_scores <- meta %>%
    group_by(seurat_clusters) %>%
    summarise(across(all_of(score_cols), mean), .groups = "drop")
  cs_mat <- as.matrix(cluster_scores[, score_cols])
  rownames(cs_mat) <- cluster_scores$seurat_clusters
  colnames(cs_mat) <- gsub("^score_|1$", "", colnames(cs_mat))
  cluster_assignment <- apply(cs_mat, 1, function(x) colnames(cs_mat)[which.max(x)])
  log_msg("Cluster assignments:")
  for (cl in names(cluster_assignment)) log_msg("  cluster ", cl, " -> ", cluster_assignment[cl])

  obj$cell_type <- unname(cluster_assignment[as.character(obj$seurat_clusters)])

  # Write cluster annotation table
  ca_df <- data.frame(cluster = rownames(cs_mat), assigned = cluster_assignment, cs_mat)
  write.csv(ca_df, file.path(OUT_DIR, "cluster_annotation.csv"), row.names = FALSE)

  saveRDS(obj, CHECKPOINT_ANN)
  log_msg("Saved annotation checkpoint")
}

ct_table <- table(obj$cell_type)
log_msg("Cell type composition:")
for (ct in names(ct_table)) log_msg("  ", ct, ": ", ct_table[ct])
write.csv(as.data.frame(ct_table), file.path(OUT_DIR, "celltype_composition.csv"),
          row.names = FALSE)

# ============================================================
# STEP 5: AUCell UPR scoring
# ============================================================
CHECKPOINT_AUC <- file.path(OUT_DIR, "checkpoint_step5_aucell.rds")

if (file.exists(CHECKPOINT_AUC)) {
  log_msg("Loading AUCell checkpoint...")
  obj <- readRDS(CHECKPOINT_AUC)
} else {
  log_msg("STEP 5: AUCell UPR scoring")
  DefaultAssay(obj) <- "SCT"
  exprMatrix <- GetAssayData(obj, assay = "SCT", layer = "data")
  log_msg("Building AUCell rankings on ", ncol(exprMatrix), " cells, ", nrow(exprMatrix), " genes")

  cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats = FALSE,
                                         splitByBlocks = TRUE, verbose = FALSE)
  cells_AUC <- AUCell_calcAUC(upr_arms, cells_rankings, verbose = FALSE)
  auc_mat <- t(getAUC(cells_AUC))

  obj$IRE1_XBP1 <- auc_mat[, "IRE1_XBP1"]
  obj$PERK_ATF4 <- auc_mat[, "PERK_ATF4"]
  obj$ATF6      <- auc_mat[, "ATF6"]
  obj$UPR_broad <- auc_mat[, "UPR_broad"]

  saveRDS(obj, CHECKPOINT_AUC)
  log_msg("AUCell scoring complete and saved")
}

# ============================================================
# STEP 6: Validate 4 patterns
# ============================================================
log_msg("STEP 6: Pattern validation")
md <- obj@meta.data

# ---- Pattern 1: Cell-type UPR ranking ----
log_msg("Pattern 1: Cell-type UPR ranking")
kw1 <- kruskal.test(UPR_broad ~ cell_type, data = md)
ct_summary <- md %>% group_by(cell_type) %>%
  summarise(median_UPR = median(UPR_broad), mean_UPR = mean(UPR_broad), n = n()) %>%
  arrange(desc(median_UPR))
log_msg("  KW p = ", format.pval(kw1$p.value))
log_msg("  Top cell types by median UPR:")
print(ct_summary)
top3 <- head(ct_summary$cell_type, 3)
primary_top3 <- c("Endothelial","Malignant","TAM")
p1_concord <- length(intersect(top3, primary_top3)) >= 2  # at least 2 of 3
write.csv(ct_summary, file.path(OUT_DIR, "pattern1_celltype_upr.csv"), row.names = FALSE)

# ---- Pattern 2: AC-like vs MES-like in malignant ----
log_msg("Pattern 2: AC-like vs MES-like UPR")
mal <- subset(obj, cell_type == "Malignant")
log_msg("  Malignant cells: ", ncol(mal))
subtype_sets <- list(
  AC_like  = intersect(c("GFAP","S100B","AQP4","ALDH1L1","CD44"), rownames(mal)),
  MES_like = intersect(c("VIM","FN1","CHI3L1","CD44","ANXA1","ANXA2"), rownames(mal)),
  NPC_like = intersect(c("SOX2","SOX11","ASCL1","DLX5","DCX"), rownames(mal)),
  OPC_like = intersect(c("OLIG1","OLIG2","PDGFRA","SOX10"), rownames(mal))
)
for (st in names(subtype_sets)) {
  mal <- AddModuleScore(mal, features = list(subtype_sets[[st]]),
                        name = paste0("st_", st), assay = "SCT")
}
sub_cols <- paste0("st_", names(subtype_sets), "1")
sub_mat <- as.matrix(mal@meta.data[, sub_cols])
colnames(sub_mat) <- names(subtype_sets)
mal$malignant_subtype <- colnames(sub_mat)[apply(sub_mat, 1, which.max)]
sub_table <- table(mal$malignant_subtype)
log_msg("  Subtype composition:"); print(sub_table)

ac_upr  <- mal$UPR_broad[mal$malignant_subtype == "AC_like"]
mes_upr <- mal$UPR_broad[mal$malignant_subtype == "MES_like"]
wt_p2 <- wilcox.test(ac_upr, mes_upr)
ac_med  <- median(ac_upr); mes_med <- median(mes_upr)
log_msg(sprintf("  AC median=%.4f (n=%d), MES median=%.4f (n=%d), p=%s",
                ac_med, length(ac_upr), mes_med, length(mes_upr),
                format.pval(wt_p2$p.value)))
p2_concord <- ac_med > mes_med   # primary: AC > MES
sub_df <- data.frame(subtype = c("AC_like","MES_like","NPC_like","OPC_like"),
                     n = as.integer(sub_table[c("AC_like","MES_like","NPC_like","OPC_like")]),
                     median_UPR = c(ac_med, mes_med,
                                    median(mal$UPR_broad[mal$malignant_subtype=="NPC_like"]),
                                    median(mal$UPR_broad[mal$malignant_subtype=="OPC_like"])))
write.csv(sub_df, file.path(OUT_DIR, "pattern2_malignant_subtypes.csv"), row.names = FALSE)

# ---- Pattern 3: TAM M2 enrichment in UPR-high TAMs ----
log_msg("Pattern 3: TAM M1/M2 by UPR group")
tam <- subset(obj, cell_type %in% c("TAM","Microglia"))
log_msg("  TAM/Microglia cells: ", ncol(tam))
tam_med <- median(tam$UPR_broad)
tam$UPR_group <- ifelse(tam$UPR_broad >= tam_med, "UPR_high", "UPR_low")

m1_set <- intersect(c("TNF","IL6","IL1B","CD86","NOS2","CXCL10"), rownames(tam))
m2_set <- intersect(c("MRC1","ARG1","CD163","TGFB1","IL10","MERTK"), rownames(tam))
tam <- AddModuleScore(tam, features = list(m1_set), name = "M1_score", assay = "SCT")
tam <- AddModuleScore(tam, features = list(m2_set), name = "M2_score", assay = "SCT")

m2_high <- tam$M2_score1[tam$UPR_group == "UPR_high"]
m2_low  <- tam$M2_score1[tam$UPR_group == "UPR_low"]
m1_high <- tam$M1_score1[tam$UPR_group == "UPR_high"]
m1_low  <- tam$M1_score1[tam$UPR_group == "UPR_low"]
wt_p3_m2 <- wilcox.test(m2_high, m2_low)
wt_p3_m1 <- wilcox.test(m1_high, m1_low)
log_msg(sprintf("  M2 high vs low: median %.4f vs %.4f, p=%s",
                median(m2_high), median(m2_low), format.pval(wt_p3_m2$p.value)))
log_msg(sprintf("  M1 high vs low: median %.4f vs %.4f, p=%s",
                median(m1_high), median(m1_low), format.pval(wt_p3_m1$p.value)))
p3_concord <- median(m2_high) > median(m2_low)  # primary: M2 enriched in UPR-high

# ---- Pattern 4: CD8 exhaustion paradox ----
log_msg("Pattern 4: CD8 exhaustion vs sample-level UPR-high/low tumors")
mal_obj <- subset(obj, cell_type == "Malignant")
sample_upr <- mal_obj@meta.data %>% group_by(sample_id) %>%
  summarise(med_UPR_mal = median(UPR_broad), n_mal = n())
sample_med <- median(sample_upr$med_UPR_mal)
sample_upr$tumor_UPR_group <- ifelse(sample_upr$med_UPR_mal >= sample_med,
                                     "UPR_high","UPR_low")
log_msg("  Sample-level grouping (n=", nrow(sample_upr), ")")
write.csv(sample_upr, file.path(OUT_DIR, "pattern4_sample_upr.csv"), row.names = FALSE)

cd8 <- subset(obj, cell_type == "CD8_T")
log_msg("  CD8 T cells: ", ncol(cd8))
cd8$tumor_UPR_group <- sample_upr$tumor_UPR_group[match(cd8$sample_id, sample_upr$sample_id)]

exh_set <- intersect(c("LAG3","TIGIT","HAVCR2","PDCD1","ENTPD1","TOX","CTLA4"),
                     rownames(cd8))
cd8 <- AddModuleScore(cd8, features = list(exh_set), name = "Exh_score", assay = "SCT")

exh_high <- cd8$Exh_score1[cd8$tumor_UPR_group == "UPR_high"]
exh_low  <- cd8$Exh_score1[cd8$tumor_UPR_group == "UPR_low"]
wt_p4 <- wilcox.test(exh_high, exh_low)
log_msg(sprintf("  Exh in UPR-high tumors=%.4f (n=%d), in UPR-low=%.4f (n=%d), p=%s",
                median(exh_high), length(exh_high),
                median(exh_low),  length(exh_low),
                format.pval(wt_p4$p.value)))
# primary finding: UPR-high has LOWER exhaustion -> validation concord if med_high < med_low
p4_concord <- median(exh_high) < median(exh_low)

# ============================================================
# STEP 7: Figures (PDF, theme_classic, 10pt)
# ============================================================
log_msg("STEP 7: Generating PDF figures")
my_theme <- theme_classic(base_size = 10)

# UMAP by cell type
p_ct <- DimPlot(obj, group.by = "cell_type", label = TRUE, repel = TRUE,
                raster = FALSE, pt.size = 0.1) + my_theme +
  ggtitle("Cell types (Harmony-integrated)")
ggsave(file.path(FIG_DIR, "Val_UMAP_celltype.pdf"), p_ct, width = 8, height = 6)

# UMAP UPR_broad with viridis
p_upr <- FeaturePlot(obj, features = "UPR_broad", raster = FALSE, pt.size = 0.1,
                     order = TRUE) +
  scale_color_viridis_c(option = "viridis") + my_theme +
  ggtitle("UPR_broad (AUCell)")
ggsave(file.path(FIG_DIR, "Val_UMAP_UPR.pdf"), p_upr, width = 7, height = 6)

# Pattern 1 violin
ord <- ct_summary$cell_type
md_p <- md
md_p$cell_type <- factor(md_p$cell_type, levels = ord)
p_vln <- ggplot(md_p, aes(x = cell_type, y = UPR_broad, fill = cell_type)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", alpha = 0.7) +
  my_theme + theme(axis.text.x = element_text(angle = 45, hjust = 1),
                   legend.position = "none") +
  labs(x = NULL, y = "UPR_broad (AUCell)",
       title = sprintf("Pattern 1: UPR by cell type (KW p=%s)",
                       format.pval(kw1$p.value, digits = 3)))
ggsave(file.path(FIG_DIR, "Val_Pattern1_violin.pdf"), p_vln, width = 8, height = 5)

# Pattern 2 subtype
mal_md <- mal@meta.data
p_sub <- ggplot(mal_md, aes(x = malignant_subtype, y = UPR_broad,
                            fill = malignant_subtype)) +
  geom_violin(scale = "width", trim = TRUE) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white", alpha = 0.7) +
  my_theme + theme(legend.position = "none") +
  labs(x = NULL, y = "UPR_broad",
       title = sprintf("Pattern 2: Malignant subtype UPR (AC vs MES p=%s)",
                       format.pval(wt_p2$p.value, digits = 3)))
ggsave(file.path(FIG_DIR, "Val_Pattern2_subtype.pdf"), p_sub, width = 7, height = 5)

# Pattern 3 TAM M1/M2
tam_md <- tam@meta.data %>% select(UPR_group, M1_score1, M2_score1) %>%
  pivot_longer(c(M1_score1, M2_score1), names_to = "polarization", values_to = "score") %>%
  mutate(polarization = gsub("_score1","", polarization))
p_tam <- ggplot(tam_md, aes(x = polarization, y = score, fill = UPR_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85) +
  my_theme +
  labs(x = NULL, y = "Module score",
       title = sprintf("Pattern 3: TAM M1/M2 by UPR group (M2 p=%s)",
                       format.pval(wt_p3_m2$p.value, digits = 3)))
ggsave(file.path(FIG_DIR, "Val_Pattern3_TAM.pdf"), p_tam, width = 6, height = 5)

# Pattern 4 CD8 exhaustion
cd8_md <- cd8@meta.data %>% filter(!is.na(tumor_UPR_group))
p_cd8 <- ggplot(cd8_md, aes(x = tumor_UPR_group, y = Exh_score1,
                            fill = tumor_UPR_group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.85) +
  my_theme + theme(legend.position = "none") +
  labs(x = "Tumor-level UPR group", y = "CD8 exhaustion score",
       title = sprintf("Pattern 4: CD8 exhaustion (p=%s)",
                       format.pval(wt_p4$p.value, digits = 3)))
ggsave(file.path(FIG_DIR, "Val_Pattern4_CD8.pdf"), p_cd8, width = 5, height = 5)

# ============================================================
# STEP 8: Summary CSV, RDS, MD report
# ============================================================
log_msg("STEP 8: Writing summary outputs")

eff <- function(a, b) {
  # rank-biserial-like effect: median diff
  median(a) - median(b)
}

summary_df <- data.frame(
  pattern = c("P1: Cell-type UPR ranking",
              "P2: AC-like > MES-like UPR",
              "P3: TAM M2 enrichment in UPR-high",
              "P4: CD8 exhaustion lower in UPR-high tumors"),
  primary_finding = c("Endothelial > Malignant > TAM (top 3)",
                      "AC > MES",
                      "M2 score higher in UPR-high TAMs",
                      "Exhaustion LOWER in UPR-high tumors"),
  validation_finding = c(paste("Top3:", paste(top3, collapse=", ")),
                         sprintf("AC_med=%.4f, MES_med=%.4f", ac_med, mes_med),
                         sprintf("M2_high=%.4f, M2_low=%.4f",
                                 median(m2_high), median(m2_low)),
                         sprintf("Exh_high=%.4f, Exh_low=%.4f",
                                 median(exh_high), median(exh_low))),
  p_value = c(kw1$p.value, wt_p2$p.value, wt_p3_m2$p.value, wt_p4$p.value),
  effect_size = c(NA,
                  eff(ac_upr, mes_upr),
                  eff(m2_high, m2_low),
                  eff(exh_high, exh_low)),
  concordant = c(p1_concord, p2_concord, p3_concord, p4_concord)
)
write.csv(summary_df, file.path(OUT_DIR, "gse182109_validation_summary.csv"),
          row.names = FALSE)

# Save final RDS
log_msg("Saving final Seurat object")
saveRDS(obj, file.path(OUT_DIR, "gse182109_seu.rds"))

# Markdown report
n_total <- ncol(obj)
qc_df <- read.csv(file.path(OUT_DIR, "qc_per_sample.csv"))
report <- c(
  "# GSE182109 Validation Report",
  "",
  paste0("Generated: ", format(Sys.time())),
  "",
  "## Dataset",
  paste0("- Source: GSE182109 (Abdelfattah et al., glioma scRNA-seq)"),
  paste0("- Samples used (ndGBM + rGBM, IDH-WT proxy): ", length(samples)),
  paste0("- Total cells (post-QC, post-doublet): **", n_total, "**"),
  paste0("- Cells loaded total: ", sum(qc_df$cells_loaded),
         "; post-QC: ", sum(qc_df$cells_post_filter),
         "; post-doublet: ", sum(qc_df$cells_post_doublet)),
  "",
  "## Pipeline",
  "Per-sample 10X load -> QC (nFeature 200-6000, nCount 500-50000, %mt<20) -> scDblFinder -> ",
  "merge -> SCTransform v2 (regress %mt) -> PCA(50) -> Harmony -> UMAP/Leiden(res=0.4) -> ",
  "marker-based annotation (8 cell types via AddModuleScore) -> AUCell UPR scoring (4 arms).",
  "",
  "## Cell type composition",
  "| Cell type | n cells |",
  "|---|---|",
  paste0("| ", names(ct_table), " | ", as.integer(ct_table), " |"),
  "",
  "## Pattern validation",
  "| Pattern | Primary finding | Validation finding | p-value | Concordant |",
  "|---|---|---|---|---|",
  paste0("| ", summary_df$pattern, " | ", summary_df$primary_finding, " | ",
         summary_df$validation_finding, " | ",
         format.pval(summary_df$p_value, digits = 3), " | ",
         ifelse(summary_df$concordant, "YES", "NO"), " |"),
  "",
  "## Output files",
  paste0("- Final Seurat: `", file.path(OUT_DIR,"gse182109_seu.rds"),"`"),
  paste0("- Summary CSV: `", file.path(OUT_DIR,"gse182109_validation_summary.csv"),"`"),
  paste0("- Figures: `", FIG_DIR, "/Val_*.pdf`"),
  paste0("- Per-sample QC: `", file.path(OUT_DIR,"qc_per_sample.csv"),"`"),
  paste0("- Cluster annotation: `", file.path(OUT_DIR,"cluster_annotation.csv"),"`"),
  "",
  "## Notes",
  "- Cell types assigned by max-mean module score per cluster (Seurat AddModuleScore).",
  "- Microglia and TAM clusters were merged for Pattern 3 (myeloid compartment).",
  "- IDH-WT proxy: only ndGBM and rGBM samples; LGG samples excluded by name pattern."
)
writeLines(report, file.path(OUT_DIR, "gse182109_validation_report.md"))

log_msg("=== PIPELINE COMPLETE ===")
sessionInfo()
