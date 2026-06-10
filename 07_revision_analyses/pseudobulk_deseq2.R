###############################################################################
# pseudobulk_deseq2.R
# PATIENT-LEVEL pseudobulk DESeq2 for glioma UPR manuscript revision
# Replaces underpowered sample-level module-score Wilcoxon with properly-powered
# pseudobulk DESeq2 testing.
#
# Primary: TAM M2 polarization (discovery + validation)
# Secondary: Malignant AC-vs-MES UPR, CD8 T exhaustion
###############################################################################

suppressMessages({
  library(Seurat)
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(Matrix)
  library(fgsea)
  library(limma)
})

set.seed(42)
PROJECT_DIR <- getwd()
RES_DIR <- file.path(PROJECT_DIR, "results", "revision_stats")
dir.create(RES_DIR, showWarnings = FALSE, recursive = TRUE)

# Gene sets
m2_markers <- c("ARG1", "MRC1", "CD163", "IL10", "TGFB1", "CCL22", "IRF4", "STAT6")
t_exhaustion_genes <- c("PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4",
                        "ENTPD1", "LAYN", "TOX", "TOX2", "CXCL13",
                        "BATF", "IRF4", "NFATC1", "NR4A1", "NR4A2",
                        "EOMES", "TBX21", "GZMB", "PRF1", "IFNG")
upr_genes <- c("XBP1", "ATF4", "ATF6", "ERN1", "EIF2AK3", "HSPA5",
               "DDIT3", "PPP1R15A", "TRIB3", "HERPUD1", "DNAJB9",
               "PDIA4", "HYOU1", "SDF2L1", "CRELD2")

cat("========================================================================\n")
cat("PSEUDOBULK DESEQ2 ANALYSIS\n")
cat("========================================================================\n\n")

# =============================================================================
# 1. LOAD DISCOVERY DATA
# =============================================================================
cat("--- 1. Loading discovery data (Smart-seq2 subset) ---\n")
seu_full <- readRDS(file.path(PROJECT_DIR, "data", "processed", "seu_upr_scored.rds"))
cat(sprintf("  Full object: %d cells, %d features\n", ncol(seu_full), nrow(seu_full)))

valid_samples <- grep("^(MGH|BT)[0-9]+$", unique(seu_full$sample), value = TRUE)
seu_disc <- subset(seu_full, sample %in% valid_samples)
cat(sprintf("  Smart-seq2 subset: %d cells, %d samples\n",
            ncol(seu_disc), length(unique(seu_disc$sample))))

disc_counts <- GetAssayData(seu_disc, assay = "RNA", layer = "counts")
cat(sprintf("  Raw counts dim: %d genes x %d cells\n", nrow(disc_counts), ncol(disc_counts)))

disc_ct_tbl <- table(seu_disc$celltype)
cat("  Cell types:\n")
for (ct in names(disc_ct_tbl)) {
  n_samp <- length(unique(seu_disc$sample[seu_disc$celltype == ct]))
  cat(sprintf("    %s: %d cells, %d samples\n", ct, disc_ct_tbl[ct], n_samp))
}

# =============================================================================
# 2. LOAD VALIDATION DATA
# =============================================================================
cat("\n--- 2. Loading validation data (GSE182109) ---\n")
seu_val <- readRDS(file.path(PROJECT_DIR, "results", "GSE182109_validation", "gse182109_seu.rds"))
cat(sprintf("  Validation object: %d cells, %d features\n", ncol(seu_val), nrow(seu_val)))

cat("  Joining layers (40 split layers -> unified counts)...\n")
seu_val <- JoinLayers(seu_val, assay = "RNA")
val_counts <- GetAssayData(seu_val, assay = "RNA", layer = "counts")
cat(sprintf("  Unified counts dim: %d genes x %d cells\n", nrow(val_counts), ncol(val_counts)))

val_ct_tbl <- table(seu_val$cell_type)
cat("  Cell types:\n")
for (ct in names(val_ct_tbl)) {
  n_samp <- length(unique(seu_val$sample_id[seu_val$cell_type == ct]))
  cat(sprintf("    %s: %d cells, %d samples\n", ct, val_ct_tbl[ct], n_samp))
}

# =============================================================================
# 3. DEFINE UPR GROUPS (sample-level, based on malignant cell median UPR)
# =============================================================================
cat("\n--- 3. Defining UPR groups ---\n")

# Discovery
mal_disc_meta <- seu_disc@meta.data %>%
  dplyr::filter(celltype == "Malignant") %>%
  dplyr::group_by(sample) %>%
  dplyr::summarise(
    n_mal = n(),
    med_UPR = median(UPR_score, na.rm = TRUE),
    .groups = "drop"
  )
disc_upr_cutoff <- median(mal_disc_meta$med_UPR, na.rm = TRUE)
mal_disc_meta$UPR_group <- ifelse(mal_disc_meta$med_UPR > disc_upr_cutoff, "UPR-high", "UPR-low")
cat(sprintf("  Discovery cutoff: %.4f, high=%d, low=%d\n",
            disc_upr_cutoff,
            sum(mal_disc_meta$UPR_group == "UPR-high"),
            sum(mal_disc_meta$UPR_group == "UPR-low")))

# Validation
mal_val_meta <- seu_val@meta.data %>%
  dplyr::filter(cell_type == "Malignant") %>%
  dplyr::group_by(sample_id) %>%
  dplyr::summarise(
    n_mal = n(),
    med_UPR = median(UPR_broad, na.rm = TRUE),
    .groups = "drop"
  )
val_upr_cutoff <- median(mal_val_meta$med_UPR, na.rm = TRUE)
mal_val_meta$UPR_group <- ifelse(mal_val_meta$med_UPR > val_upr_cutoff, "UPR-high", "UPR-low")
cat(sprintf("  Validation cutoff: %.4f, high=%d, low=%d\n",
            val_upr_cutoff,
            sum(mal_val_meta$UPR_group == "UPR-high"),
            sum(mal_val_meta$UPR_group == "UPR-low")))

# =============================================================================
# 4. PSEUDOBULK AGGREGATION FUNCTION
# =============================================================================

pseudobulk_aggregate <- function(counts, sample_ids, cell_mask, min_cells = 10) {
  samp_sub <- sample_ids[cell_mask]
  cnt_sub <- counts[, cell_mask, drop = FALSE]

  samp_tbl <- table(samp_sub)
  valid_samps <- names(samp_tbl)[samp_tbl >= min_cells]

  if (length(valid_samps) == 0) {
    return(list(pb = NULL, n_cells = NULL, samples = character(0)))
  }

  pb_list <- list()
  for (s in valid_samps) {
    idx <- which(samp_sub == s)
    if (length(idx) == 1) {
      pb_list[[s]] <- cnt_sub[, idx]
    } else {
      pb_list[[s]] <- Matrix::rowSums(cnt_sub[, idx, drop = FALSE])
    }
  }
  pb <- do.call(cbind, pb_list)
  colnames(pb) <- valid_samps

  n_cells <- samp_tbl[valid_samps]

  list(pb = pb, n_cells = n_cells, samples = valid_samps)
}

# =============================================================================
# 5. DESeq2 + fgsea PIPELINE FUNCTION
# =============================================================================

run_deseq2_and_test <- function(pb_counts, sample_info, gene_set, gene_set_name,
                                 individual_genes = NULL, min_count_filter = 10) {

  # Match samples
  common_samps <- intersect(colnames(pb_counts), rownames(sample_info))
  if (length(common_samps) < 4) {
    cat(sprintf("  ERROR: only %d common samples (need >=4)\n", length(common_samps)))
    return(NULL)
  }
  pb <- pb_counts[, common_samps, drop = FALSE]
  si <- sample_info[common_samps, , drop = FALSE]

  # Count groups
  n_high <- sum(si$UPR_group == "UPR-high")
  n_low  <- sum(si$UPR_group == "UPR-low")
  cat(sprintf("  Samples: UPR-high=%d, UPR-low=%d\n", n_high, n_low))

  # Filter low-count genes
  pb <- round(pb)  # ensure integer
  gene_sums <- Matrix::rowSums(pb)
  keep <- gene_sums >= min_count_filter
  pb <- pb[keep, , drop = FALSE]
  pb <- as.matrix(pb)
  cat(sprintf("  Genes after count filter (>=%d): %d\n", min_count_filter, nrow(pb)))

  # Create DESeq2 object
  dds <- DESeqDataSetFromMatrix(
    countData = pb,
    colData = si,
    design = ~ UPR_group
  )
  dds$UPR_group <- factor(dds$UPR_group, levels = c("UPR-low", "UPR-high"))

  # Run DESeq2
  dds <- DESeq2::DESeq(dds, quiet = FALSE)
  cat(sprintf("  DESeq2 done. Size factors: min=%.3f max=%.3f\n",
              min(sizeFactors(dds)), max(sizeFactors(dds))))

  # Extract results — use resultsNames to get the correct coefficient name
  # DESeq2 replaces hyphens with dots in factor level names
  res_names <- resultsNames(dds)
  cat(sprintf("  resultsNames: %s\n", paste(res_names, collapse = ", ")))
  target_coef <- grep("UPR_group.*UPR.*high.*vs.*UPR.*low", res_names, value = TRUE)
  if (length(target_coef) == 0) {
    cat("  ERROR: could not find UPR_group coefficient in resultsNames\n")
    return(NULL)
  }
  cat(sprintf("  Using coefficient: %s\n", target_coef[1]))

  # Try ashr first, fall back to normal
  res <- tryCatch({
    lfcShrink(dds, coef = target_coef[1], type = "ashr", quiet = FALSE)
  }, error = function(e) {
    cat(sprintf("  WARNING: ashr failed (%s), using normal shrinkage\n", e$message))
    lfcShrink(dds, coef = target_coef[1], type = "normal")
  })

  # Build result table with explicit gene column
  res_df <- data.frame(
    gene = rownames(res),
    baseMean = res$baseMean,
    log2FoldChange = res$log2FoldChange,
    lfcSE = res$lfcSE,
    pvalue = res$pvalue,
    padj = res$padj,
    stat = NA_real_,
    stringsAsFactors = FALSE
  )
  rownames(res_df) <- NULL

  # Compute stat for fgsea ranking: signed -log10(p) * LFC
  res_df$stat <- res_df$log2FoldChange * (-log10(pmax(res_df$pvalue, 1e-300)))

  n_sig_05 <- sum(res_df$padj < 0.05, na.rm = TRUE)
  n_sig_10 <- sum(res_df$padj < 0.10, na.rm = TRUE)
  cat(sprintf("  DEGs: padj<0.05 = %d, padj<0.10 = %d\n", n_sig_05, n_sig_10))

  # --- Individual gene tests ---
  gene_results <- NULL
  if (!is.null(individual_genes)) {
    genes_found <- individual_genes[individual_genes %in% res_df$gene]
    cat(sprintf("  Individual genes tested: %d/%d\n", length(genes_found), length(individual_genes)))
    if (length(genes_found) > 0) {
      gene_results <- res_df[res_df$gene %in% genes_found, ]
      gene_results$direction <- ifelse(gene_results$log2FoldChange > 0, "UP in UPR-high", "DOWN in UPR-high")
      gene_results$verdict <- ifelse(gene_results$padj < 0.05, "significant",
                              ifelse(gene_results$padj < 0.2, "trend", "not significant"))
      gene_results <- gene_results[, c("gene","baseMean","log2FoldChange","pvalue","padj","direction","verdict")]

      for (i in seq_len(nrow(gene_results))) {
        cat(sprintf("    %s: log2FC=%.3f, p=%.2e, padj=%.2e [%s] %s\n",
                    gene_results$gene[i],
                    gene_results$log2FoldChange[i],
                    gene_results$pvalue[i],
                    gene_results$padj[i],
                    gene_results$direction[i],
                    gene_results$verdict[i]))
      }
    } else {
      cat("  WARNING: none of the individual genes found in result!\n")
      cat(sprintf("  Checking: input genes = %s\n", paste(head(individual_genes, 5), collapse=", ")))
      cat(sprintf("  First 5 result genes: %s\n", paste(head(res_df$gene, 5), collapse=", ")))
    }
  }

  # --- fgsea ---
  # Build ranked list: use stat column, named by gene
  ranks <- res_df$stat
  names(ranks) <- res_df$gene
  ranks <- sort(ranks, decreasing = TRUE)
  # Remove duplicates (should not exist, but safety)
  ranks <- ranks[!duplicated(names(ranks))]
  ranks <- ranks[is.finite(ranks)]

  # Filter gene set to genes present in ranks
  gene_set_filt <- lapply(gene_set, function(gs) intersect(gs, names(ranks)))
  gene_set_filt <- gene_set_filt[lengths(gene_set_filt) >= 3]

  fgsea_res <- NULL
  if (length(gene_set_filt) > 0) {
    fgsea_res <- fgsea::fgsea(
      pathways = gene_set_filt,
      stats = ranks,
      minSize = 3,
      maxSize = 500,
      nPermSimple = 10000
    )

    for (i in seq_len(nrow(fgsea_res))) {
      cat(sprintf("  fgsea %s: NES=%.3f, pval=%.4f, padj=%.4f, size=%d\n",
                  fgsea_res$pathway[i],
                  fgsea_res$NES[i],
                  fgsea_res$pval[i],
                  fgsea_res$padj[i],
                  fgsea_res$size[i]))
      if (length(fgsea_res$leadingEdge[[i]]) > 0) {
        cat(sprintf("    leadingEdge: %s\n",
                    paste(head(fgsea_res$leadingEdge[[i]], 6), collapse=", ")))
      }
    }
  } else {
    cat("  fgsea: no gene sets with >=3 genes present in ranks\n")
  }

  # --- Mean LFC for gene set ---
  mean_lfc <- NA_real_
  if (!is.null(gene_set) && length(gene_set) > 0) {
    gs_all <- unique(unlist(gene_set))
    gs_in_res <- gs_all[gs_all %in% res_df$gene]
    if (length(gs_in_res) > 0) {
      mean_lfc <- mean(res_df$log2FoldChange[res_df$gene %in% gs_in_res], na.rm = TRUE)
      cat(sprintf("  Mean gene-set LFC (n=%d genes): %.3f (%s)\n",
                  length(gs_in_res), mean_lfc,
                  ifelse(mean_lfc > 0, "UP in UPR-high", "DOWN in UPR-high")))
    }
  }

  list(
    res = res_df,
    gene_results = gene_results,
    fgsea = fgsea_res,
    mean_lfc = mean_lfc,
    n_samples = ncol(pb),
    n_genes = nrow(pb),
    n_high = n_high,
    n_low = n_low,
    dds = dds
  )
}

# =============================================================================
# 6. PRIMARY: TAM M2 POLARIZATION — DISCOVERY
# =============================================================================
cat("\n========================================================================\n")
cat("6. PRIMARY: TAM M2 polarization -- DISCOVERY cohort\n")
cat("========================================================================\n")

tam_mask_disc <- seu_disc$celltype == "TAM"
cat(sprintf("  TAM cells: %d across %d samples\n",
            sum(tam_mask_disc), length(unique(seu_disc$sample[tam_mask_disc]))))

pb_disc <- pseudobulk_aggregate(disc_counts, seu_disc$sample, tam_mask_disc, min_cells = 10)
cat(sprintf("  Pseudobulk samples (>=10 TAM): %d\n", length(pb_disc$samples)))
cat(sprintf("  Excluded (<10 TAM): %s\n",
            paste(setdiff(unique(seu_disc$sample[tam_mask_disc]), pb_disc$samples), collapse = ", ")))

# Build sample info
si_disc <- data.frame(
  sample = pb_disc$samples,
  stringsAsFactors = FALSE
)
si_disc$med_UPR <- mal_disc_meta$med_UPR[match(si_disc$sample, mal_disc_meta$sample)]
si_disc$UPR_group <- ifelse(si_disc$med_UPR > disc_upr_cutoff, "UPR-high", "UPR-low")
si_disc$n_tam <- as.integer(pb_disc$n_cells[si_disc$sample])
rownames(si_disc) <- si_disc$sample

cat(sprintf("  TAM cells per sample: min=%d, max=%d, median=%d\n",
            min(si_disc$n_tam), max(si_disc$n_tam), median(si_disc$n_tam)))

disc_tam_res <- run_deseq2_and_test(
  pb_counts = pb_disc$pb,
  sample_info = si_disc,
  gene_set = list(M2_polarization = m2_markers),
  gene_set_name = "M2_polarization",
  individual_genes = m2_markers,
  min_count_filter = 10
)

# =============================================================================
# 7. PRIMARY: TAM M2 POLARIZATION — VALIDATION
# =============================================================================
cat("\n========================================================================\n")
cat("7. PRIMARY: TAM M2 polarization -- VALIDATION cohort\n")
cat("========================================================================\n")

tam_mask_val <- seu_val$cell_type == "TAM"
cat(sprintf("  TAM cells: %d across %d samples\n",
            sum(tam_mask_val), length(unique(seu_val$sample_id[tam_mask_val]))))

pb_val <- pseudobulk_aggregate(val_counts, seu_val$sample_id, tam_mask_val, min_cells = 10)
cat(sprintf("  Pseudobulk samples (>=10 TAM): %d\n", length(pb_val$samples)))

si_val <- data.frame(
  sample_id = pb_val$samples,
  stringsAsFactors = FALSE
)
si_val$med_UPR <- mal_val_meta$med_UPR[match(si_val$sample_id, mal_val_meta$sample_id)]
si_val$UPR_group <- ifelse(si_val$med_UPR > val_upr_cutoff, "UPR-high", "UPR-low")
si_val$n_tam <- as.integer(pb_val$n_cells[si_val$sample_id])
rownames(si_val) <- si_val$sample_id

cat(sprintf("  TAM cells per sample: min=%d, max=%d, median=%d\n",
            min(si_val$n_tam), max(si_val$n_tam), median(si_val$n_tam)))

val_tam_res <- run_deseq2_and_test(
  pb_counts = pb_val$pb,
  sample_info = si_val,
  gene_set = list(M2_polarization = m2_markers),
  gene_set_name = "M2_polarization",
  individual_genes = m2_markers,
  min_count_filter = 10
)

# =============================================================================
# 8. SECONDARY: Malignant AC-like vs MES-like UPR (Discovery)
# =============================================================================
cat("\n========================================================================\n")
cat("8. SECONDARY: Malignant AC-like vs MES-like UPR -- DISCOVERY\n")
cat("========================================================================\n")

mal_disc <- seu_disc@meta.data %>%
  dplyr::filter(celltype == "Malignant" & malignant_subtype %in% c("AC_like", "MES_like"))

cat(sprintf("  Malignant AC/MES cells: %d\n", nrow(mal_disc)))
cat(sprintf("  AC_like: %d cells in %d samples\n",
            sum(mal_disc$malignant_subtype == "AC_like"),
            length(unique(mal_disc$sample[mal_disc$malignant_subtype == "AC_like"]))))
cat(sprintf("  MES_like: %d cells in %d samples\n",
            sum(mal_disc$malignant_subtype == "MES_like"),
            length(unique(mal_disc$sample[mal_disc$malignant_subtype == "MES_like"]))))

# Find samples with both subtypes (paired design)
mal_sample_subtype <- mal_disc %>%
  dplyr::group_by(sample, malignant_subtype) %>%
  dplyr::summarise(n = n(), .groups = "drop") %>%
  dplyr::filter(n >= 5)

paired_samples <- mal_sample_subtype %>%
  tidyr::pivot_wider(names_from = malignant_subtype, values_from = n) %>%
  tidyr::drop_na() %>%
  dplyr::pull(sample)

cat(sprintf("  Samples with both AC_like AND MES_like (>=5 cells each): %d\n",
            length(paired_samples)))

ac_mes_res <- NULL
if (length(paired_samples) >= 4) {
  # Aggregate pseudobulk per sample per subtype
  pb_mal_list <- list()
  for (s in paired_samples) {
    for (st in c("AC_like", "MES_like")) {
      mask <- seu_disc$sample == s & seu_disc$celltype == "Malignant" &
              seu_disc$malignant_subtype == st
      cnt <- disc_counts[, mask, drop = FALSE]
      pb_mal_list[[paste0(s, "_", st)]] <- Matrix::rowSums(cnt)
    }
  }
  pb_mal <- do.call(cbind, pb_mal_list)

  coldata_mal <- data.frame(
    row.names = colnames(pb_mal),
    sample = gsub("_(AC_like|MES_like)$", "", colnames(pb_mal)),
    subtype = ifelse(grepl("_AC_like$", colnames(pb_mal)), "AC_like", "MES_like"),
    stringsAsFactors = FALSE
  )

  # Filter and convert
  gene_sums <- Matrix::rowSums(pb_mal)
  pb_mal <- pb_mal[gene_sums >= 10, ]
  pb_mal <- as.matrix(round(pb_mal))

  # DESeq2 with paired design
  dds_mal <- DESeqDataSetFromMatrix(
    countData = pb_mal,
    colData = coldata_mal,
    design = ~ sample + subtype
  )
  dds_mal$subtype <- factor(dds_mal$subtype, levels = c("MES_like", "AC_like"))
  dds_mal <- DESeq2::DESeq(dds_mal, quiet = FALSE)

  res_mal <- tryCatch({
    lfcShrink(dds_mal, coef = "subtype_AC_like_vs_MES_like", type = "ashr", quiet = FALSE)
  }, error = function(e) {
    cat(sprintf("  WARNING: ashr failed, using normal: %s\n", e$message))
    lfcShrink(dds_mal, coef = "subtype_AC_like_vs_MES_like", type = "normal")
  })

  res_mal_df <- data.frame(
    gene = rownames(res_mal),
    baseMean = res_mal$baseMean,
    log2FoldChange = res_mal$log2FoldChange,
    lfcSE = res_mal$lfcSE,
    pvalue = res_mal$pvalue,
    padj = res_mal$padj,
    stringsAsFactors = FALSE
  )
  rownames(res_mal_df) <- NULL

  cat(sprintf("  DESeq2: %d genes tested, %d paired samples\n",
              nrow(res_mal_df), length(paired_samples)))
  cat(sprintf("  Significant (padj<0.05): %d, (padj<0.10): %d\n",
              sum(res_mal_df$padj < 0.05, na.rm = TRUE),
              sum(res_mal_df$padj < 0.10, na.rm = TRUE)))

  # UPR gene program
  upr_found <- upr_genes[upr_genes %in% res_mal_df$gene]
  cat(sprintf("  UPR genes found: %d/%d\n", length(upr_found), length(upr_genes)))

  ac_mes_upr <- data.frame()
  for (g in upr_found) {
    row <- res_mal_df[res_mal_df$gene == g, ]
    verdict <- ifelse(!is.na(row$padj) && row$padj < 0.05, "significant",
                      ifelse(!is.na(row$padj) && row$padj < 0.2, "trend", "not significant"))
    ac_mes_upr <- rbind(ac_mes_upr, data.frame(
      gene = g, log2FC = round(row$log2FoldChange, 3),
      pvalue = signif(row$pvalue, 3), padj = signif(row$padj, 3),
      direction = ifelse(row$log2FoldChange > 0, "UP in AC-like", "UP in MES-like"),
      verdict = verdict,
      stringsAsFactors = FALSE
    ))
    cat(sprintf("    %s: log2FC=%.3f, p=%.2e, padj=%.2e [%s] %s\n",
                g, row$log2FoldChange, row$pvalue, row$padj,
                ifelse(row$log2FoldChange > 0, "AC>MES", "MES>AC"), verdict))
  }

  # fgsea on UPR program
  ranks_mal <- res_mal_df$log2FoldChange * (-log10(pmax(res_mal_df$pvalue, 1e-300)))
  names(ranks_mal) <- res_mal_df$gene
  ranks_mal <- sort(ranks_mal, decreasing = TRUE)
  ranks_mal <- ranks_mal[!duplicated(names(ranks_mal))]
  ranks_mal <- ranks_mal[is.finite(ranks_mal)]

  upr_set <- list(UPR_program = upr_genes[upr_genes %in% names(ranks_mal)])
  if (length(upr_set$UPR_program) >= 3) {
    fgsea_upr <- fgsea::fgsea(pathways = upr_set, stats = ranks_mal,
                              minSize = 3, nPermSimple = 10000)
    if (nrow(fgsea_upr) > 0) {
      cat(sprintf("  fgsea UPR: NES=%.3f, pval=%.4f, padj=%.4f, size=%d\n",
                  fgsea_upr$NES[1], fgsea_upr$pval[1], fgsea_upr$padj[1], fgsea_upr$size[1]))
    }
  }

  mean_upr_lfc <- mean(res_mal_df$log2FoldChange[res_mal_df$gene %in% upr_genes], na.rm = TRUE)
  cat(sprintf("  Mean UPR gene LFC (AC vs MES): %.3f (%s)\n",
              mean_upr_lfc, ifelse(mean_upr_lfc > 0, "AC-like higher", "MES-like higher")))

  ac_mes_res <- list(
    genes = ac_mes_upr,
    fgsea = if (exists("fgsea_upr")) fgsea_upr else NULL,
    mean_lfc = mean_upr_lfc,
    n_pairs = length(paired_samples),
    full_res = res_mal_df
  )
} else {
  cat("  INSUFFICIENT paired samples for DESeq2 (need >=4)\n")
}

# =============================================================================
# 9. SECONDARY: CD8 T exhaustion — BOTH cohorts
# =============================================================================
cat("\n========================================================================\n")
cat("9. SECONDARY: CD8 T exhaustion\n")
cat("========================================================================\n")

# --- Discovery ---
cat("\n  --- Discovery CD8 T exhaustion ---\n")
cd8_mask_disc <- seu_disc$celltype == "CD8 T"
cat(sprintf("  CD8 T cells: %d across %d samples\n",
            sum(cd8_mask_disc), length(unique(seu_disc$sample[cd8_mask_disc]))))

pb_cd8_disc <- pseudobulk_aggregate(disc_counts, seu_disc$sample, cd8_mask_disc, min_cells = 5)
cat(sprintf("  Pseudobulk samples (>=5 CD8 T): %d\n", length(pb_cd8_disc$samples)))

si_cd8_disc <- data.frame(
  sample = pb_cd8_disc$samples,
  stringsAsFactors = FALSE
)
si_cd8_disc$med_UPR <- mal_disc_meta$med_UPR[match(si_cd8_disc$sample, mal_disc_meta$sample)]
si_cd8_disc$UPR_group <- ifelse(si_cd8_disc$med_UPR > disc_upr_cutoff, "UPR-high", "UPR-low")
si_cd8_disc$n_cd8 <- as.integer(pb_cd8_disc$n_cells[si_cd8_disc$sample])
rownames(si_cd8_disc) <- si_cd8_disc$sample

n_high_cd8d <- sum(si_cd8_disc$UPR_group == "UPR-high")
n_low_cd8d  <- sum(si_cd8_disc$UPR_group == "UPR-low")
cat(sprintf("  UPR-high=%d, UPR-low=%d\n", n_high_cd8d, n_low_cd8d))

disc_cd8_res <- NULL
if (n_high_cd8d >= 2 && n_low_cd8d >= 2) {
  disc_cd8_res <- run_deseq2_and_test(
    pb_counts = pb_cd8_disc$pb,
    sample_info = si_cd8_disc,
    gene_set = list(CD8_exhaustion = t_exhaustion_genes),
    gene_set_name = "CD8_exhaustion",
    individual_genes = t_exhaustion_genes,
    min_count_filter = 10
  )
} else {
  cat(sprintf("  INSUFFICIENT: need >=2 per group (high=%d, low=%d)\n", n_high_cd8d, n_low_cd8d))
}

# --- Validation ---
cat("\n  --- Validation CD8 T exhaustion ---\n")
cd8_mask_val <- seu_val$cell_type == "CD8_T"
cat(sprintf("  CD8 T cells: %d across %d samples\n",
            sum(cd8_mask_val), length(unique(seu_val$sample_id[cd8_mask_val]))))

pb_cd8_val <- pseudobulk_aggregate(val_counts, seu_val$sample_id, cd8_mask_val, min_cells = 5)
cat(sprintf("  Pseudobulk samples (>=5 CD8 T): %d\n", length(pb_cd8_val$samples)))

si_cd8_val <- data.frame(
  sample_id = pb_cd8_val$samples,
  stringsAsFactors = FALSE
)
si_cd8_val$med_UPR <- mal_val_meta$med_UPR[match(si_cd8_val$sample_id, mal_val_meta$sample_id)]
si_cd8_val$UPR_group <- ifelse(si_cd8_val$med_UPR > val_upr_cutoff, "UPR-high", "UPR-low")
si_cd8_val$n_cd8 <- as.integer(pb_cd8_val$n_cells[si_cd8_val$sample_id])
rownames(si_cd8_val) <- si_cd8_val$sample_id

n_high_cd8v <- sum(si_cd8_val$UPR_group == "UPR-high")
n_low_cd8v  <- sum(si_cd8_val$UPR_group == "UPR-low")
cat(sprintf("  UPR-high=%d, UPR-low=%d\n", n_high_cd8v, n_low_cd8v))

val_cd8_res <- NULL
if (n_high_cd8v >= 3 && n_low_cd8v >= 3) {
  val_cd8_res <- run_deseq2_and_test(
    pb_counts = pb_cd8_val$pb,
    sample_info = si_cd8_val,
    gene_set = list(CD8_exhaustion = t_exhaustion_genes),
    gene_set_name = "CD8_exhaustion",
    individual_genes = t_exhaustion_genes,
    min_count_filter = 10
  )
} else {
  cat(sprintf("  INSUFFICIENT: need >=3 per group (high=%d, low=%d)\n", n_high_cd8v, n_low_cd8v))
}

# =============================================================================
# 10. SAVE RESULTS
# =============================================================================
cat("\n========================================================================\n")
cat("10. Saving results\n")
cat("========================================================================\n")

# TAM M2 Discovery
if (!is.null(disc_tam_res)) {
  write.csv(disc_tam_res$gene_results,
            file.path(RES_DIR, "pseudobulk_tam_m2_deseq2_discovery.csv"),
            row.names = FALSE)
  write.csv(disc_tam_res$res,
            file.path(RES_DIR, "pseudobulk_tam_deseq2_discovery_all_genes.csv"),
            row.names = FALSE)
}

# TAM M2 Validation
if (!is.null(val_tam_res)) {
  write.csv(val_tam_res$gene_results,
            file.path(RES_DIR, "pseudobulk_tam_m2_deseq2_validation.csv"),
            row.names = FALSE)
  write.csv(val_tam_res$res,
            file.path(RES_DIR, "pseudobulk_tam_deseq2_validation_all_genes.csv"),
            row.names = FALSE)
}

# fgsea summary
fgsea_rows <- list()
add_fgsea_row <- function(cohort, ct, test, fgobj) {
  if (is.null(fgobj) || nrow(fgobj) == 0) return(NULL)
  # Convert leadingEdge list column to character
  fgobj$leadingEdge <- sapply(fgobj$leadingEdge, function(x) paste(x, collapse = ";"))
  data.frame(cohort = cohort, celltype = ct, test = test, fgobj,
             stringsAsFactors = FALSE)
}
fgsea_rows <- list()
fgsea_rows[[length(fgsea_rows) + 1]] <- add_fgsea_row("Discovery", "TAM", "M2_polarization",
  if (!is.null(disc_tam_res)) disc_tam_res$fgsea else NULL)
fgsea_rows[[length(fgsea_rows) + 1]] <- add_fgsea_row("Validation", "TAM", "M2_polarization",
  if (!is.null(val_tam_res)) val_tam_res$fgsea else NULL)
fgsea_rows[[length(fgsea_rows) + 1]] <- add_fgsea_row("Discovery", "CD8_T", "CD8_exhaustion",
  if (!is.null(disc_cd8_res)) disc_cd8_res$fgsea else NULL)
fgsea_rows[[length(fgsea_rows) + 1]] <- add_fgsea_row("Validation", "CD8_T", "CD8_exhaustion",
  if (!is.null(val_cd8_res)) val_cd8_res$fgsea else NULL)
fgsea_rows <- fgsea_rows[!sapply(fgsea_rows, is.null)]
if (length(fgsea_rows) > 0) {
  fgsea_summary <- do.call(rbind, fgsea_rows)
  write.csv(fgsea_summary,
            file.path(RES_DIR, "pseudobulk_fgsea_summary.csv"),
            row.names = FALSE)
}

# CD8 exhaustion gene results
if (!is.null(disc_cd8_res) && !is.null(disc_cd8_res$gene_results)) {
  write.csv(disc_cd8_res$gene_results,
            file.path(RES_DIR, "pseudobulk_cd8_exhaustion_deseq2_discovery.csv"),
            row.names = FALSE)
}
if (!is.null(val_cd8_res) && !is.null(val_cd8_res$gene_results)) {
  write.csv(val_cd8_res$gene_results,
            file.path(RES_DIR, "pseudobulk_cd8_exhaustion_deseq2_validation.csv"),
            row.names = FALSE)
}

# AC vs MES
if (!is.null(ac_mes_res)) {
  write.csv(ac_mes_res$genes,
            file.path(RES_DIR, "pseudobulk_ac_vs_mes_upr_deseq2_discovery.csv"),
            row.names = FALSE)
}

# =============================================================================
# 11. FINAL SUMMARY TABLE
# =============================================================================
cat("\n========================================================================\n")
cat("FINAL SUMMARY -- PSEUDOBULK DESEQ2 RESULTS\n")
cat("========================================================================\n\n")

cat(sprintf("%-45s | %-10s | %-22s | %-15s | %-15s | %-12s | %s\n",
            "Claim", "Cohort", "N high/low", "Gene-set NES", "Gene-set padj",
            "Mean LFC", "Verdict"))
cat(strrep("-", 145), "\n")

# --- TAM M2 Discovery ---
if (!is.null(disc_tam_res)) {
  if (!is.null(disc_tam_res$fgsea) && nrow(disc_tam_res$fgsea) > 0) {
    nes <- sprintf("%.3f", disc_tam_res$fgsea$NES[1])
    padj <- sprintf("%.3f", disc_tam_res$fgsea$padj[1])
    verdict <- if (disc_tam_res$fgsea$padj[1] < 0.05) "SIGNIFICANT" else
               if (disc_tam_res$fgsea$padj[1] < 0.2) "TREND" else "NOT significant"
  } else {
    nes <- "N/A"; padj <- "N/A"; verdict <- "fgsea failed"
  }
  mlfc <- ifelse(!is.na(disc_tam_res$mean_lfc),
                 sprintf("%+.3f", disc_tam_res$mean_lfc), "N/A")
  cat(sprintf("%-45s | %-10s | high=%2d low=%2d        | %-15s | %-15s | %-12s | %s\n",
              "TAM M2 polarization", "Discovery",
              disc_tam_res$n_high, disc_tam_res$n_low, nes, padj, mlfc, verdict))
} else {
  cat(sprintf("%-45s | %-10s | %-22s | %-15s | %-15s | %-12s | %s\n",
              "TAM M2 polarization", "Discovery", "FAILED", "", "", "", "EXECUTION FAILED"))
}

# --- TAM M2 Validation ---
if (!is.null(val_tam_res)) {
  if (!is.null(val_tam_res$fgsea) && nrow(val_tam_res$fgsea) > 0) {
    nes <- sprintf("%.3f", val_tam_res$fgsea$NES[1])
    padj <- sprintf("%.3f", val_tam_res$fgsea$padj[1])
    verdict <- if (val_tam_res$fgsea$padj[1] < 0.05) "SIGNIFICANT" else
               if (val_tam_res$fgsea$padj[1] < 0.2) "TREND" else "NOT significant"
  } else {
    nes <- "N/A"; padj <- "N/A"; verdict <- "fgsea failed"
  }
  mlfc <- ifelse(!is.na(val_tam_res$mean_lfc),
                 sprintf("%+.3f", val_tam_res$mean_lfc), "N/A")
  cat(sprintf("%-45s | %-10s | high=%2d low=%2d        | %-15s | %-15s | %-12s | %s\n",
              "TAM M2 polarization", "Validation",
              val_tam_res$n_high, val_tam_res$n_low, nes, padj, mlfc, verdict))
}

# --- AC vs MES ---
if (!is.null(ac_mes_res)) {
  n_sig <- sum(ac_mes_res$genes$padj < 0.05, na.rm = TRUE)
  if (!is.null(ac_mes_res$fgsea) && nrow(ac_mes_res$fgsea) > 0) {
    nes <- sprintf("%.3f", ac_mes_res$fgsea$NES[1])
    padj <- sprintf("%.3f", ac_mes_res$fgsea$padj[1])
    verdict <- if (ac_mes_res$fgsea$padj[1] < 0.05) "SIGNIFICANT" else
               if (ac_mes_res$fgsea$padj[1] < 0.2) "TREND" else "NOT significant"
  } else {
    nes <- "N/A"; padj <- "N/A"
    verdict <- if (n_sig >= 3) "SIGNIFICANT (genes)" else
               if (n_sig >= 1) "TREND (genes)" else "NOT significant"
  }
  mlfc <- sprintf("%+.3f", ac_mes_res$mean_lfc)
  cat(sprintf("%-45s | %-10s | n_pairs=%2d            | %-15s | %-15s | %-12s | %s\n",
              "Malignant AC>MES UPR", "Discovery",
              ac_mes_res$n_pairs, nes, padj, mlfc, verdict))
}

# --- CD8 Exhaustion Discovery ---
if (!is.null(disc_cd8_res)) {
  if (!is.null(disc_cd8_res$fgsea) && nrow(disc_cd8_res$fgsea) > 0) {
    nes <- sprintf("%.3f", disc_cd8_res$fgsea$NES[1])
    padj <- sprintf("%.3f", disc_cd8_res$fgsea$padj[1])
    verdict <- if (disc_cd8_res$fgsea$padj[1] < 0.05) "SIGNIFICANT" else
               if (disc_cd8_res$fgsea$padj[1] < 0.2) "TREND" else "NOT significant"
  } else {
    nes <- "N/A"; padj <- "N/A"; verdict <- "fgsea failed"
  }
  mlfc <- ifelse(!is.na(disc_cd8_res$mean_lfc),
                 sprintf("%+.3f", disc_cd8_res$mean_lfc), "N/A")
  cat(sprintf("%-45s | %-10s | high=%2d low=%2d        | %-15s | %-15s | %-12s | %s\n",
              "CD8 T exhaustion", "Discovery",
              disc_cd8_res$n_high, disc_cd8_res$n_low, nes, padj, mlfc, verdict))
} else {
  cat(sprintf("%-45s | %-10s | high=%2d low=%2d        | %-15s | %-15s | %-12s | %s\n",
              "CD8 T exhaustion", "Discovery",
              n_high_cd8d, n_low_cd8d, "", "", "", "INSUFFICIENT SAMPLES"))
}

# --- CD8 Exhaustion Validation ---
if (!is.null(val_cd8_res)) {
  if (!is.null(val_cd8_res$fgsea) && nrow(val_cd8_res$fgsea) > 0) {
    nes <- sprintf("%.3f", val_cd8_res$fgsea$NES[1])
    padj <- sprintf("%.3f", val_cd8_res$fgsea$padj[1])
    verdict <- if (val_cd8_res$fgsea$padj[1] < 0.05) "SIGNIFICANT" else
               if (val_cd8_res$fgsea$padj[1] < 0.2) "TREND" else "NOT significant"
  } else {
    nes <- "N/A"; padj <- "N/A"; verdict <- "fgsea failed"
  }
  mlfc <- ifelse(!is.na(val_cd8_res$mean_lfc),
                 sprintf("%+.3f", val_cd8_res$mean_lfc), "N/A")
  cat(sprintf("%-45s | %-10s | high=%2d low=%2d        | %-15s | %-15s | %-12s | %s\n",
              "CD8 T exhaustion", "Validation",
              val_cd8_res$n_high, val_cd8_res$n_low, nes, padj, mlfc, verdict))
} else {
  cat(sprintf("%-45s | %-10s | high=%2d low=%2d        | %-15s | %-15s | %-12s | %s\n",
              "CD8 T exhaustion", "Validation",
              n_high_cd8v, n_low_cd8v, "", "", "", "INSUFFICIENT SAMPLES"))
}

cat("\n========================================================================\n")
cat("ANALYSIS COMPLETE\n")
cat("========================================================================\n")
cat(sprintf("Output directory: %s\n", RES_DIR))
cat("Files written:\n")
for (f in list.files(RES_DIR, pattern = "pseudobulk.*\\.csv$", full.names = TRUE)) {
  cat(sprintf("  %s\n", f))
}
writeLines(capture.output(sessionInfo()),
           file.path(RES_DIR, "pseudobulk_deseq2_sessionInfo.txt"))
