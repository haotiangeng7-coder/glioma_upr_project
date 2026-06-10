###############################################################################
# k2_vs_k3_decision.R
# Data-driven decision: K=2 vs K=3 for UPR-based glioma subtyping
# Analyses: full cohort, IDH-WT subgroup, and IDH-WT independent re-clustering
###############################################################################

set.seed(42)

# === Load libraries ===
suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(ConsensusClusterPlus)
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(cluster)
})

# === Load config ===
source("00_setup/config.R")

# === Load data ===
load(file.path(DATA_PROC, "consensus_clustering_results.RData"))
load(file.path(DATA_PROC, "tcga_glioma_expression.RData"))
load(file.path(DATA_PROC, "upr_gene_sets.RData"))

cat("================================================================\n")
cat("  K=2 vs K=3 Data-Driven Decision Analysis\n")
cat("================================================================\n\n")

###############################################################################
# HELPER FUNCTIONS
###############################################################################

# Assign cluster labels: order by median survival (descending = best prognosis first)
assign_prognosis_labels <- function(cluster_vec, surv_time, surv_event, k) {
  # Compute median survival per cluster to determine prognosis ordering
  df_tmp <- data.frame(cl = cluster_vec, time = surv_time, event = surv_event)
  df_tmp <- df_tmp[complete.cases(df_tmp), ]

  med_surv <- tapply(seq_len(nrow(df_tmp)), df_tmp$cl, function(idx) {
    fit <- survfit(Surv(time, event) ~ 1, data = df_tmp[idx, ])
    sm <- surv_median(fit)
    return(sm$median)
  })

  # Replace NA medians with Inf (group so good that median not reached)
  med_surv[is.na(med_surv)] <- Inf

  # Order: best prognosis (highest median) = favorable
  rank_order <- order(med_surv, decreasing = TRUE)

  if (k == 2) {
    label_map <- setNames(c("UPR-favorable", "UPR-high-risk"), rank_order)
  } else if (k == 3) {
    label_map <- setNames(c("UPR-favorable", "UPR-intermediate", "UPR-high-risk"), rank_order)
  } else {
    label_map <- setNames(paste0("Cluster_", seq_len(k)), rank_order)
  }

  all_samples <- names(cluster_vec)
  result <- label_map[as.character(cluster_vec)]
  names(result) <- all_samples
  return(result)
}

# Compute survival analysis metrics for a given K
compute_surv_metrics <- function(surv_df, cluster_col, k, label = "") {
  valid <- complete.cases(surv_df$OS.time, surv_df$OS, surv_df[[cluster_col]])
  df <- surv_df[valid, ]
  df[[cluster_col]] <- factor(df[[cluster_col]])

  # Overall log-rank
  formula_str <- as.formula(paste0("Surv(OS.time, OS) ~ ", cluster_col))
  sd_test <- survdiff(formula_str, data = df)
  pval <- 1 - pchisq(sd_test$chisq, df = length(sd_test$n) - 1)

  # Median OS per group
  fit <- survfit(formula_str, data = df)
  med_tbl <- surv_median(fit)

  cat(sprintf("\n--- %s (K=%d) ---\n", label, k))
  cat(sprintf("Overall log-rank p-value: %.2e\n", pval))
  cat("Median OS (days) per group:\n")
  print(med_tbl)

  # Median OS gap
  medians <- med_tbl$median
  medians_valid <- medians[!is.na(medians)]
  if (length(medians_valid) >= 2) {
    med_gap <- max(medians_valid) - min(medians_valid)
    cat(sprintf("Max median OS gap: %.0f days (%.1f months)\n", med_gap, med_gap / 30.44))
  } else {
    # If some groups don't reach median, use the range of available + note
    med_gap <- NA
    if (length(medians_valid) == 1) {
      cat(sprintf("Only 1 group reached median OS (%.0f days); others have better prognosis (median not reached)\n",
                  medians_valid[1]))
    } else {
      cat("No group reached median OS\n")
    }
  }

  # Pairwise log-rank (for K >= 3)
  pairwise_p <- NULL
  if (k >= 3) {
    groups <- sort(unique(as.character(df[[cluster_col]])))
    cat("\nPairwise log-rank tests:\n")
    pw_results <- list()
    for (i in 1:(length(groups) - 1)) {
      for (j in (i + 1):length(groups)) {
        sub <- df[as.character(df[[cluster_col]]) %in% c(groups[i], groups[j]), ]
        sub$pw_group <- factor(as.character(sub[[cluster_col]]))
        pw_fit <- survdiff(Surv(OS.time, OS) ~ pw_group, data = sub)
        pw_p <- 1 - pchisq(pw_fit$chisq, df = 1)
        n_i <- sum(as.character(sub[[cluster_col]]) == groups[i])
        n_j <- sum(as.character(sub[[cluster_col]]) == groups[j])
        cat(sprintf("  %s vs %s: p = %.2e (n=%d vs n=%d)\n",
                    groups[i], groups[j], pw_p, n_i, n_j))
        pw_results[[paste(groups[i], "vs", groups[j])]] <- pw_p
      }
    }
    pairwise_p <- pw_results
  }

  return(list(
    pval = pval,
    median_tbl = med_tbl,
    median_gap = med_gap,
    pairwise_p = pairwise_p,
    n_per_group = table(df[[cluster_col]]),
    fit = fit
  ))
}

# Create KM plot (returns a ggplot object, no risk table to avoid grid issues)
create_km_plot <- function(surv_df, cluster_col, k, title_str, colors = NULL) {
  valid <- complete.cases(surv_df$OS.time, surv_df$OS, surv_df[[cluster_col]])
  df <- surv_df[valid, ]
  # Rename to standard column to avoid formula environment issues with ggsurvplot
  df$Subtype <- factor(df[[cluster_col]])

  fit <- survfit(Surv(OS.time, OS) ~ Subtype, data = df)

  if (is.null(colors)) {
    if (k == 2) {
      colors <- c("#4DBBD5", "#E64B35")
    } else if (k == 3) {
      colors <- c("#4DBBD5", "#E64B35", "#00A087")
    } else {
      colors <- scales::hue_pal()(k)
    }
  }

  p <- ggsurvplot(
    fit, data = df,
    pval = TRUE, pval.method = TRUE,
    risk.table = FALSE,
    palette = colors,
    title = title_str,
    xlab = "Time (days)",
    ylab = "Overall Survival Probability",
    legend.title = "Subtype",
    ggtheme = theme_bw() + theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
    )
  )
  return(p$plot)
}

# Create KM plot with risk table (standalone PDF)
create_km_pdf <- function(surv_df, cluster_col, k, title_str, pdf_path, colors = NULL,
                           width = 8, height = 7) {
  valid <- complete.cases(surv_df$OS.time, surv_df$OS, surv_df[[cluster_col]])
  df <- surv_df[valid, ]
  df$Subtype <- factor(df[[cluster_col]])

  fit <- survfit(Surv(OS.time, OS) ~ Subtype, data = df)

  if (is.null(colors)) {
    if (k == 2) {
      colors <- c("#4DBBD5", "#E64B35")
    } else if (k == 3) {
      colors <- c("#4DBBD5", "#E64B35", "#00A087")
    } else {
      colors <- scales::hue_pal()(k)
    }
  }

  pdf(pdf_path, width = width, height = height)
  p <- ggsurvplot(
    fit, data = df,
    pval = TRUE, pval.method = TRUE,
    risk.table = TRUE,
    risk.table.height = 0.25,
    palette = colors,
    title = title_str,
    xlab = "Time (days)",
    ylab = "Overall Survival Probability",
    legend.title = "Subtype",
    ggtheme = theme_bw() + theme(
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
    )
  )
  print(p)
  dev.off()
  cat(sprintf("  Saved: %s\n", pdf_path))
}


###############################################################################
# PART 1: Full Cohort K=2 vs K=3
###############################################################################
cat("\n================================================================\n")
cat("PART 1: FULL COHORT ANALYSIS\n")
cat("================================================================\n")

# Get cluster assignments
cls_k2_full <- cc_results[[2]]$consensusClass
cls_k3_full <- cc_results[[3]]$consensusClass

full_samples <- names(cls_k2_full)

# Match survival data
surv_match <- surv_final[surv_final$barcode %in% full_samples, ]

# Assign labels based on prognosis (survival)
labels_k2_full <- assign_prognosis_labels(
  cls_k2_full[surv_match$barcode],
  surv_match$OS.time, surv_match$OS, 2
)
# Extend to all samples using the mapping
k2_map <- tapply(labels_k2_full, cls_k2_full[names(labels_k2_full)], function(x) x[1])
labels_k2_all <- k2_map[as.character(cls_k2_full)]
names(labels_k2_all) <- names(cls_k2_full)

labels_k3_full <- assign_prognosis_labels(
  cls_k3_full[surv_match$barcode],
  surv_match$OS.time, surv_match$OS, 3
)
k3_map <- tapply(labels_k3_full, cls_k3_full[names(labels_k3_full)], function(x) x[1])
labels_k3_all <- k3_map[as.character(cls_k3_full)]
names(labels_k3_all) <- names(cls_k3_full)

# Build survival data frames
surv_k2_full <- surv_final
surv_k2_full$subtype_k2 <- labels_k2_all[surv_k2_full$barcode]

surv_k3_full <- surv_final
surv_k3_full$subtype_k3 <- labels_k3_all[surv_k3_full$barcode]

# Compute metrics
cat(sprintf("\nFull cohort samples: %d\n", nrow(surv_final)))
cat(sprintf("K=2 group sizes: %s\n", paste(table(surv_k2_full$subtype_k2), collapse=", ")))
cat(sprintf("K=3 group sizes: %s\n", paste(table(surv_k3_full$subtype_k3), collapse=", ")))

res_k2_full <- compute_surv_metrics(surv_k2_full, "subtype_k2", 2, "Full Cohort")
res_k3_full <- compute_surv_metrics(surv_k3_full, "subtype_k3", 3, "Full Cohort")


###############################################################################
# PART 2: IDH-WT Subgroup K=2 vs K=3
###############################################################################
cat("\n================================================================\n")
cat("PART 2: IDH-WT SUBGROUP ANALYSIS\n")
cat("================================================================\n")

cls_k2_idhwt <- cc_results_idhwt[[2]]$consensusClass
cls_k3_idhwt <- cc_results_idhwt[[3]]$consensusClass

idhwt_samples <- names(cls_k2_idhwt)
idhwt_barcodes_surv <- surv_final$barcode[surv_final$IDH_status == "WT" & !is.na(surv_final$IDH_status)]
common_idhwt <- intersect(idhwt_samples, idhwt_barcodes_surv)
cat(sprintf("IDH-WT samples: clustering=%d, survival=%d, overlap=%d\n",
            length(idhwt_samples), length(idhwt_barcodes_surv), length(common_idhwt)))

# Assign labels based on prognosis
surv_idhwt_base <- surv_final[surv_final$barcode %in% common_idhwt, ]

labels_k2_idhwt <- assign_prognosis_labels(
  cls_k2_idhwt[surv_idhwt_base$barcode],
  surv_idhwt_base$OS.time, surv_idhwt_base$OS, 2
)

labels_k3_idhwt <- assign_prognosis_labels(
  cls_k3_idhwt[surv_idhwt_base$barcode],
  surv_idhwt_base$OS.time, surv_idhwt_base$OS, 3
)

surv_idhwt <- surv_idhwt_base
surv_idhwt$subtype_k2 <- labels_k2_idhwt[surv_idhwt$barcode]
surv_idhwt$subtype_k3 <- labels_k3_idhwt[surv_idhwt$barcode]

res_k2_idhwt <- compute_surv_metrics(surv_idhwt, "subtype_k2", 2, "IDH-WT Subgroup")
res_k3_idhwt <- compute_surv_metrics(surv_idhwt, "subtype_k3", 3, "IDH-WT Subgroup")


###############################################################################
# PART 3: IDH-WT Independent Re-clustering
###############################################################################
cat("\n================================================================\n")
cat("PART 3: IDH-WT INDEPENDENT RE-CLUSTERING\n")
cat("================================================================\n")

# Get IDH-WT expression data
expr_mat <- expr_tpm_symbol
idhwt_expr <- expr_mat[, common_idhwt]
cat(sprintf("IDH-WT expression matrix: %d genes x %d samples\n", nrow(idhwt_expr), ncol(idhwt_expr)))

# Log2 transform
idhwt_log2 <- log2(idhwt_expr + 1)

# Get UPR broad genes available
upr_broad_avail <- intersect(UPR_broad_genes, rownames(idhwt_log2))
cat(sprintf("UPR broad genes available: %d / %d\n", length(upr_broad_avail), length(UPR_broad_genes)))

# MAD filtering
upr_expr_idhwt <- idhwt_log2[upr_broad_avail, ]
gene_mad <- apply(upr_expr_idhwt, 1, mad)
mad_threshold <- quantile(gene_mad, probs = CC_MAD_PERCENTILE)
selected_genes <- names(gene_mad[gene_mad >= mad_threshold])
cat(sprintf("MAD filtering: %d genes pass (threshold MAD >= %.4f, percentile=%.0f%%)\n",
            length(selected_genes), mad_threshold, CC_MAD_PERCENTILE * 100))
cat("Selected genes:", paste(selected_genes, collapse = ", "), "\n")

# Prepare matrix
cc_mat <- as.matrix(upr_expr_idhwt[selected_genes, ])
cc_mat_scaled <- t(scale(t(cc_mat)))

# Run ConsensusClusterPlus
cat("\nRunning ConsensusClusterPlus on IDH-WT (K=2-4, reps=1000)...\n")
cc_dir <- file.path(FIG_DIR, "cc_idhwt_recluster")
dir.create(cc_dir, showWarnings = FALSE, recursive = TRUE)

cc_idhwt_new <- ConsensusClusterPlus(
  d = cc_mat_scaled,
  maxK = 4,
  reps = CC_REPS,
  pItem = CC_PITEM,
  pFeature = 1,
  clusterAlg = "hc",
  distance = "pearson",
  innerLinkage = "ward.D2",
  finalLinkage = "ward.D2",
  seed = SEED,
  plot = "pdf",
  title = cc_dir
)

# Compute metrics for re-clustering
cat("\nRe-clustering metrics:\n")
recluster_metrics <- data.frame(K = 2:4, PAC = NA, Silhouette = NA,
                                 logrank_p = NA, stringsAsFactors = FALSE)

for (k_idx in 1:3) {
  k <- k_idx + 1
  cls <- cc_idhwt_new[[k]]$consensusClass

  # PAC
  cmat <- cc_idhwt_new[[k]]$consensusMatrix
  cdf_vals <- sort(cmat[upper.tri(cmat)])
  pac <- mean(cdf_vals > 0.1 & cdf_vals < 0.9)

  # Silhouette
  dist_mat <- as.dist(1 - cmat)
  sil <- silhouette(cls, dist_mat)
  sil_mean <- mean(sil[, "sil_width"])

  # Survival
  surv_tmp <- surv_idhwt_base
  surv_tmp$recluster <- factor(cls[match(surv_tmp$barcode, names(cls))])
  valid_tmp <- complete.cases(surv_tmp$OS.time, surv_tmp$OS, surv_tmp$recluster)
  sd_tmp <- survdiff(Surv(OS.time, OS) ~ recluster, data = surv_tmp[valid_tmp, ])
  p_tmp <- 1 - pchisq(sd_tmp$chisq, df = length(sd_tmp$n) - 1)

  recluster_metrics$PAC[k_idx] <- pac
  recluster_metrics$Silhouette[k_idx] <- sil_mean
  recluster_metrics$logrank_p[k_idx] <- p_tmp

  cat(sprintf("  K=%d: PAC=%.4f, Silhouette=%.4f, logrank p=%.2e, groups: %s\n",
              k, pac, sil_mean, p_tmp, paste(table(cls), collapse="/")))
}

cat("\nRe-clustering metrics table:\n")
print(recluster_metrics)

# Detailed survival for re-clustered K=2 and K=3
cat("\n--- IDH-WT Re-clustered Survival Analysis ---\n")

recluster_surv_results <- list()
for (k in 2:3) {
  cls_new <- cc_idhwt_new[[k]]$consensusClass
  labels_new <- assign_prognosis_labels(
    cls_new[surv_idhwt_base$barcode],
    surv_idhwt_base$OS.time, surv_idhwt_base$OS, k
  )

  surv_tmp <- surv_idhwt_base
  surv_tmp$recluster_subtype <- labels_new[surv_tmp$barcode]

  res <- compute_surv_metrics(surv_tmp, "recluster_subtype", k,
                               "IDH-WT Re-clustered")
  recluster_surv_results[[as.character(k)]] <- res
}


###############################################################################
# PART 4: Generate PDF Figures
###############################################################################
cat("\n================================================================\n")
cat("PART 4: GENERATING PDF FIGURES\n")
cat("================================================================\n")

# --- Individual KM PDFs with risk tables ---
create_km_pdf(surv_k2_full, "subtype_k2", 2, "Full Cohort: K=2 UPR Subtypes",
              file.path(FIG_DIR, "K2vsK3_fullcohort_K2_KM.pdf"),
              colors = c("#4DBBD5", "#E64B35"))

create_km_pdf(surv_k3_full, "subtype_k3", 3, "Full Cohort: K=3 UPR Subtypes",
              file.path(FIG_DIR, "K2vsK3_fullcohort_K3_KM.pdf"),
              colors = c("#4DBBD5", "#E64B35", "#00A087"))

create_km_pdf(surv_idhwt, "subtype_k2", 2, "IDH-WT: K=2 UPR Subtypes",
              file.path(FIG_DIR, "K2vsK3_idhwt_K2_KM.pdf"),
              colors = c("#4DBBD5", "#E64B35"))

create_km_pdf(surv_idhwt, "subtype_k3", 3, "IDH-WT: K=3 UPR Subtypes",
              file.path(FIG_DIR, "K2vsK3_idhwt_K3_KM.pdf"),
              colors = c("#4DBBD5", "#E64B35", "#00A087"))

# Re-clustered KM
for (k in 2:4) {
  cls_new <- cc_idhwt_new[[k]]$consensusClass
  labels_new <- assign_prognosis_labels(
    cls_new[surv_idhwt_base$barcode],
    surv_idhwt_base$OS.time, surv_idhwt_base$OS, k
  )

  surv_tmp <- surv_idhwt_base
  surv_tmp$recluster_subtype <- labels_new[surv_tmp$barcode]

  if (k == 2) cols <- c("#4DBBD5", "#E64B35")
  else if (k == 3) cols <- c("#4DBBD5", "#E64B35", "#00A087")
  else cols <- c("#4DBBD5", "#E64B35", "#00A087", "#3C5488")

  create_km_pdf(surv_tmp, "recluster_subtype", k,
                sprintf("IDH-WT Re-clustered: K=%d", k),
                file.path(FIG_DIR, sprintf("K2vsK3_idhwt_recluster_K%d_KM.pdf", k)),
                colors = cols)
}

# --- Combined comparison figure (without risk tables) ---
pdf(file.path(FIG_DIR, "K2vsK3_comparison_panel.pdf"), width = 14, height = 14)

p1 <- create_km_plot(surv_k2_full, "subtype_k2", 2, "A) Full Cohort K=2",
                      colors = c("#4DBBD5", "#E64B35"))
p2 <- create_km_plot(surv_k3_full, "subtype_k3", 3, "B) Full Cohort K=3",
                      colors = c("#4DBBD5", "#E64B35", "#00A087"))
p3 <- create_km_plot(surv_idhwt, "subtype_k2", 2, "C) IDH-WT K=2",
                      colors = c("#4DBBD5", "#E64B35"))
p4 <- create_km_plot(surv_idhwt, "subtype_k3", 3, "D) IDH-WT K=3",
                      colors = c("#4DBBD5", "#E64B35", "#00A087"))

grid.arrange(p1, p2, p3, p4, ncol = 2,
             top = textGrob("K=2 vs K=3: Survival Comparison",
                            gp = gpar(fontsize = 18, fontface = "bold")))
dev.off()
cat("  Saved: K2vsK3_comparison_panel.pdf\n")

# --- Clustering metrics figure ---
pdf(file.path(FIG_DIR, "K2vsK3_clustering_metrics.pdf"), width = 12, height = 8)

p_pac <- ggplot(metrics, aes(x = K, y = PAC)) +
  geom_line(color = "#3C5488", linewidth = 1.2) +
  geom_point(color = "#3C5488", size = 3) +
  geom_vline(xintercept = c(2, 3), linetype = "dashed", alpha = 0.5) +
  annotate("text", x = 2.1, y = max(metrics$PAC), label = "K=2", color = "gray40", size = 3.5) +
  annotate("text", x = 3.1, y = max(metrics$PAC), label = "K=3", color = "gray40", size = 3.5) +
  labs(title = "Full Cohort: PAC (lower = better)", y = "PAC", x = "K") +
  THEME_PUBLICATION

p_sil <- ggplot(metrics, aes(x = K, y = Silhouette_mean)) +
  geom_line(color = "#E64B35", linewidth = 1.2) +
  geom_point(color = "#E64B35", size = 3) +
  geom_vline(xintercept = c(2, 3), linetype = "dashed", alpha = 0.5) +
  labs(title = "Full Cohort: Mean Silhouette (higher = better)", y = "Mean Silhouette", x = "K") +
  THEME_PUBLICATION

p_delta <- ggplot(metrics, aes(x = K, y = Delta_AUC)) +
  geom_bar(stat = "identity", fill = "#4DBBD5", width = 0.6) +
  labs(title = "Full Cohort: Delta AUC (elbow method)", y = "Delta AUC", x = "K") +
  THEME_PUBLICATION

p_recluster <- ggplot(recluster_metrics, aes(x = K)) +
  geom_line(aes(y = Silhouette, color = "Silhouette"), linewidth = 1.2) +
  geom_point(aes(y = Silhouette, color = "Silhouette"), size = 3) +
  geom_line(aes(y = 1 - PAC, color = "1-PAC"), linewidth = 1.2) +
  geom_point(aes(y = 1 - PAC, color = "1-PAC"), size = 3) +
  scale_color_manual(values = c("Silhouette" = "#E64B35", "1-PAC" = "#3C5488")) +
  labs(title = "IDH-WT Re-cluster: Quality Metrics", y = "Score", x = "K", color = "") +
  THEME_PUBLICATION

grid.arrange(p_pac, p_sil, p_delta, p_recluster, ncol = 2,
             top = textGrob("Clustering Quality Metrics",
                            gp = gpar(fontsize = 16, fontface = "bold")))
dev.off()
cat("  Saved: K2vsK3_clustering_metrics.pdf\n")

# --- Group sizes figure ---
pdf(file.path(FIG_DIR, "K2vsK3_group_sizes.pdf"), width = 12, height = 6)

size_data <- rbind(
  data.frame(Analysis = "Full K=2", Group = names(res_k2_full$n_per_group),
             N = as.integer(res_k2_full$n_per_group)),
  data.frame(Analysis = "Full K=3", Group = names(res_k3_full$n_per_group),
             N = as.integer(res_k3_full$n_per_group)),
  data.frame(Analysis = "IDH-WT K=2", Group = names(res_k2_idhwt$n_per_group),
             N = as.integer(res_k2_idhwt$n_per_group)),
  data.frame(Analysis = "IDH-WT K=3", Group = names(res_k3_idhwt$n_per_group),
             N = as.integer(res_k3_idhwt$n_per_group))
)
size_data$Analysis <- factor(size_data$Analysis, levels = c("Full K=2", "Full K=3",
                                                             "IDH-WT K=2", "IDH-WT K=3"))

p_sizes <- ggplot(size_data, aes(x = Analysis, y = N, fill = Group)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = N), position = position_dodge(width = 0.9), vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("UPR-favorable" = "#4DBBD5",
                                "UPR-intermediate" = "#00A087",
                                "UPR-high-risk" = "#E64B35")) +
  labs(title = "Group Size Distribution: K=2 vs K=3", y = "Number of Samples", x = "") +
  ylim(0, max(size_data$N) * 1.15) +
  THEME_PUBLICATION +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

print(p_sizes)
dev.off()
cat("  Saved: K2vsK3_group_sizes.pdf\n")


###############################################################################
# PART 5: Decision Report
###############################################################################
cat("\n================================================================\n")
cat("PART 5: DECISION REPORT\n")
cat("================================================================\n\n")

# Helper for formatting p-values
fmt_p <- function(p) ifelse(p < 2.2e-16, "< 2.2e-16", sprintf("%.2e", p))
fmt_gap <- function(g) ifelse(is.na(g), "N/A (median not reached in some groups)", sprintf("%.0f", g))

comparison_table <- data.frame(
  Metric = c(
    "Full Cohort: log-rank p",
    "Full Cohort: Median OS gap (days)",
    "Full Cohort: PAC",
    "Full Cohort: Silhouette",
    "Full Cohort: Min group size",
    "IDH-WT: log-rank p",
    "IDH-WT: Median OS gap (days)",
    "IDH-WT: Min group size",
    "IDH-WT Recluster: log-rank p",
    "IDH-WT Recluster: PAC",
    "IDH-WT Recluster: Silhouette"
  ),
  K2 = c(
    fmt_p(res_k2_full$pval),
    fmt_gap(res_k2_full$median_gap),
    sprintf("%.4f", metrics$PAC[metrics$K == 2]),
    sprintf("%.4f", metrics$Silhouette_mean[metrics$K == 2]),
    sprintf("%d", min(res_k2_full$n_per_group)),
    fmt_p(res_k2_idhwt$pval),
    fmt_gap(res_k2_idhwt$median_gap),
    sprintf("%d", min(res_k2_idhwt$n_per_group)),
    fmt_p(recluster_metrics$logrank_p[recluster_metrics$K == 2]),
    sprintf("%.4f", recluster_metrics$PAC[recluster_metrics$K == 2]),
    sprintf("%.4f", recluster_metrics$Silhouette[recluster_metrics$K == 2])
  ),
  K3 = c(
    fmt_p(res_k3_full$pval),
    fmt_gap(res_k3_full$median_gap),
    sprintf("%.4f", metrics$PAC[metrics$K == 3]),
    sprintf("%.4f", metrics$Silhouette_mean[metrics$K == 3]),
    sprintf("%d", min(res_k3_full$n_per_group)),
    fmt_p(res_k3_idhwt$pval),
    fmt_gap(res_k3_idhwt$median_gap),
    sprintf("%d", min(res_k3_idhwt$n_per_group)),
    fmt_p(recluster_metrics$logrank_p[recluster_metrics$K == 3]),
    sprintf("%.4f", recluster_metrics$PAC[recluster_metrics$K == 3]),
    sprintf("%.4f", recluster_metrics$Silhouette[recluster_metrics$K == 3])
  ),
  stringsAsFactors = FALSE
)

cat("\n===== K=2 vs K=3 Comparison Table =====\n\n")
print(comparison_table, row.names = FALSE, right = FALSE)

# K=3 intermediate vs high-risk significance
cat("\n\n===== K=3 Pairwise Significance (Key Question) =====\n")
cat("\n[Full Cohort K=3 pairwise]:\n")
if (!is.null(res_k3_full$pairwise_p)) {
  for (name in names(res_k3_full$pairwise_p)) {
    p <- res_k3_full$pairwise_p[[name]]
    sig <- ifelse(p < 0.05, "SIGNIFICANT", "NOT significant")
    cat(sprintf("  %s: p = %.2e (%s)\n", name, p, sig))
  }
}

cat("\n[IDH-WT K=3 pairwise]:\n")
if (!is.null(res_k3_idhwt$pairwise_p)) {
  for (name in names(res_k3_idhwt$pairwise_p)) {
    p <- res_k3_idhwt$pairwise_p[[name]]
    sig <- ifelse(p < 0.05, "SIGNIFICANT", "NOT significant")
    cat(sprintf("  %s: p = %.2e (%s)\n", name, p, sig))
  }
}

cat("\n[IDH-WT Re-clustered K=3 pairwise]:\n")
if (!is.null(recluster_surv_results[["3"]]$pairwise_p)) {
  for (name in names(recluster_surv_results[["3"]]$pairwise_p)) {
    p <- recluster_surv_results[["3"]]$pairwise_p[[name]]
    sig <- ifelse(p < 0.05, "SIGNIFICANT", "NOT significant")
    cat(sprintf("  %s: p = %.2e (%s)\n", name, p, sig))
  }
}

# Decision scoring
cat("\n\n===== DECISION SCORING =====\n\n")

score_k2 <- 0
score_k3 <- 0
reasons_k2 <- character(0)
reasons_k3 <- character(0)

# 1. PAC
if (metrics$PAC[metrics$K == 2] < metrics$PAC[metrics$K == 3]) {
  reasons_k2 <- c(reasons_k2, sprintf("Lower PAC in full cohort (%.4f vs %.4f) = cleaner clusters",
                                        metrics$PAC[metrics$K == 2], metrics$PAC[metrics$K == 3]))
  score_k2 <- score_k2 + 1
} else {
  reasons_k3 <- c(reasons_k3, sprintf("Lower PAC in full cohort (%.4f vs %.4f) = cleaner clusters",
                                        metrics$PAC[metrics$K == 3], metrics$PAC[metrics$K == 2]))
  score_k3 <- score_k3 + 1
}

# 2. Silhouette
if (metrics$Silhouette_mean[metrics$K == 2] > metrics$Silhouette_mean[metrics$K == 3]) {
  reasons_k2 <- c(reasons_k2, sprintf("Higher Silhouette (%.4f vs %.4f) = better cluster separation",
                                        metrics$Silhouette_mean[metrics$K == 2],
                                        metrics$Silhouette_mean[metrics$K == 3]))
  score_k2 <- score_k2 + 1
} else {
  reasons_k3 <- c(reasons_k3, sprintf("Higher Silhouette (%.4f vs %.4f)",
                                        metrics$Silhouette_mean[metrics$K == 3],
                                        metrics$Silhouette_mean[metrics$K == 2]))
  score_k3 <- score_k3 + 1
}

# 3. IDH-WT median OS gap
if (!is.na(res_k2_idhwt$median_gap) && !is.na(res_k3_idhwt$median_gap)) {
  if (res_k3_idhwt$median_gap > res_k2_idhwt$median_gap) {
    reasons_k3 <- c(reasons_k3, sprintf("Larger median OS gap in IDH-WT (%.0f vs %.0f days)",
                                          res_k3_idhwt$median_gap, res_k2_idhwt$median_gap))
    score_k3 <- score_k3 + 1
  } else {
    reasons_k2 <- c(reasons_k2, sprintf("Larger median OS gap in IDH-WT (%.0f vs %.0f days)",
                                          res_k2_idhwt$median_gap, res_k3_idhwt$median_gap))
    score_k2 <- score_k2 + 1
  }
} else if (is.na(res_k2_idhwt$median_gap) && !is.na(res_k3_idhwt$median_gap)) {
  reasons_k3 <- c(reasons_k3, "K=3 shows computable median OS gap; K=2 does not")
  score_k3 <- score_k3 + 1
} else if (!is.na(res_k2_idhwt$median_gap) && is.na(res_k3_idhwt$median_gap)) {
  reasons_k2 <- c(reasons_k2, "K=2 shows computable median OS gap; K=3 does not")
  score_k2 <- score_k2 + 1
}

# 4. K=3 intermediate vs high-risk significance (KEY criterion, double weight)
inter_vs_high_p <- NULL
if (!is.null(res_k3_idhwt$pairwise_p)) {
  inter_high_key <- grep("intermediate.*high|high.*intermediate",
                          names(res_k3_idhwt$pairwise_p), value = TRUE)
  if (length(inter_high_key) > 0) {
    inter_vs_high_p <- res_k3_idhwt$pairwise_p[[inter_high_key[1]]]
  }
}

if (!is.null(inter_vs_high_p)) {
  if (inter_vs_high_p < 0.05) {
    reasons_k3 <- c(reasons_k3,
                     sprintf("K=3 intermediate vs high-risk IS significant in IDH-WT (p=%.2e) -> 3 distinct prognostic groups",
                             inter_vs_high_p))
    score_k3 <- score_k3 + 2
  } else {
    reasons_k2 <- c(reasons_k2,
                     sprintf("K=3 intermediate vs high-risk NOT significant in IDH-WT (p=%.2e) -> 3rd group redundant",
                             inter_vs_high_p))
    score_k2 <- score_k2 + 2
  }
}

# 5. Minimum group size
min_k3_idhwt <- min(res_k3_idhwt$n_per_group)
if (min_k3_idhwt < 50) {
  reasons_k2 <- c(reasons_k2, sprintf("K=3 smallest IDH-WT group underpowered (n=%d < 50)", min_k3_idhwt))
  score_k2 <- score_k2 + 1
} else {
  reasons_k3 <- c(reasons_k3, sprintf("K=3 all IDH-WT groups adequate size (min n=%d)", min_k3_idhwt))
}

# 6. IDH-WT re-clustering confirmation
rc_p_k2 <- recluster_metrics$logrank_p[recluster_metrics$K == 2]
rc_p_k3 <- recluster_metrics$logrank_p[recluster_metrics$K == 3]
if (rc_p_k3 < rc_p_k2) {
  reasons_k3 <- c(reasons_k3, sprintf("IDH-WT re-clustering: K=3 better log-rank (p=%.2e vs %.2e)",
                                        rc_p_k3, rc_p_k2))
  score_k3 <- score_k3 + 1
} else {
  reasons_k2 <- c(reasons_k2, sprintf("IDH-WT re-clustering: K=2 better log-rank (p=%.2e vs %.2e)",
                                        rc_p_k2, rc_p_k3))
  score_k2 <- score_k2 + 1
}

# 7. Re-clustering quality
rc_sil_k2 <- recluster_metrics$Silhouette[recluster_metrics$K == 2]
rc_sil_k3 <- recluster_metrics$Silhouette[recluster_metrics$K == 3]
if (rc_sil_k2 > rc_sil_k3) {
  reasons_k2 <- c(reasons_k2, sprintf("IDH-WT re-cluster Silhouette: K=2 better (%.4f vs %.4f)",
                                        rc_sil_k2, rc_sil_k3))
  score_k2 <- score_k2 + 1
} else {
  reasons_k3 <- c(reasons_k3, sprintf("IDH-WT re-cluster Silhouette: K=3 better (%.4f vs %.4f)",
                                        rc_sil_k3, rc_sil_k2))
  score_k3 <- score_k3 + 1
}

# Final decision
if (score_k2 > score_k3) {
  recommendation <- "K=2"
  reason_summary <- "K=2 wins on more criteria, particularly cluster quality metrics"
} else if (score_k3 > score_k2) {
  recommendation <- "K=3"
  reason_summary <- "K=3 provides finer prognostic stratification with statistically distinct groups"
} else {
  recommendation <- "K=3 (tie-break: finer prognostic resolution preferred)"
  reason_summary <- "Tied on criteria count; K=3 preferred for clinical granularity"
}

cat(sprintf("Final Score: K=2 = %d, K=3 = %d\n\n", score_k2, score_k3))

cat("Arguments FOR K=2:\n")
for (r in reasons_k2) cat(sprintf("  [+] %s\n", r))

cat("\nArguments FOR K=3:\n")
for (r in reasons_k3) cat(sprintf("  [+] %s\n", r))

cat(sprintf("\n\n========================================\n"))
cat(sprintf(">>> FINAL RECOMMENDATION: %s\n", recommendation))
cat(sprintf(">>> Rationale: %s\n", reason_summary))
cat(sprintf("========================================\n\n"))


# --- Save Decision Report PDF ---
pdf(file.path(FIG_DIR, "K2vsK3_decision_report.pdf"), width = 14, height = 10)

# Page 1: Table
grid.newpage()
pushViewport(viewport(layout = grid.layout(3, 1, heights = unit(c(0.1, 0.05, 0.85), "npc"))))

pushViewport(viewport(layout.pos.row = 1))
grid.text("K=2 vs K=3 Decision Report: Comparison Metrics",
          gp = gpar(fontsize = 18, fontface = "bold"))
popViewport()

pushViewport(viewport(layout.pos.row = 2))
grid.text(sprintf("Full cohort: n=%d | IDH-WT: n=%d | Re-clustered IDH-WT: n=%d",
                  nrow(surv_final), nrow(surv_idhwt), length(common_idhwt)),
          gp = gpar(fontsize = 11, col = "gray30"))
popViewport()

pushViewport(viewport(layout.pos.row = 3))
tbl_grob <- tableGrob(comparison_table, rows = NULL,
                       theme = ttheme_minimal(
                         core = list(fg_params = list(fontsize = 9),
                                     bg_params = list(fill = c("white", "gray95"))),
                         colhead = list(fg_params = list(fontsize = 11, fontface = "bold"),
                                        bg_params = list(fill = "#4DBBD5"))
                       ))
grid.draw(tbl_grob)
popViewport()

popViewport()

# Page 2: Recommendation
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 1)))

grid.text("FINAL RECOMMENDATION", y = 0.95,
          gp = gpar(fontsize = 20, fontface = "bold"))
grid.text(sprintf("%s", recommendation), y = 0.87,
          gp = gpar(fontsize = 24, col = "#E64B35", fontface = "bold"))
grid.text(sprintf("Score: K=2 = %d | K=3 = %d", score_k2, score_k3), y = 0.80,
          gp = gpar(fontsize = 14))
grid.text(reason_summary, y = 0.74, gp = gpar(fontsize = 12, fontface = "italic"))

y_pos <- 0.65
grid.text("Arguments FOR K=2:", x = 0.05, y = y_pos, just = "left",
          gp = gpar(fontsize = 13, fontface = "bold", col = "#4DBBD5"))
for (r in reasons_k2) {
  y_pos <- y_pos - 0.035
  grid.text(paste0("  + ", r), x = 0.07, y = y_pos, just = "left",
            gp = gpar(fontsize = 8.5))
}

y_pos <- y_pos - 0.05
grid.text("Arguments FOR K=3:", x = 0.05, y = y_pos, just = "left",
          gp = gpar(fontsize = 13, fontface = "bold", col = "#E64B35"))
for (r in reasons_k3) {
  y_pos <- y_pos - 0.035
  grid.text(paste0("  + ", r), x = 0.07, y = y_pos, just = "left",
            gp = gpar(fontsize = 8.5))
}

popViewport()
dev.off()
cat("  Saved: K2vsK3_decision_report.pdf\n")


cat("\n================================================================\n")
cat("All PDFs saved to:", FIG_DIR, "\n")
cat("Output files:\n")
cat("  - K2vsK3_fullcohort_K2_KM.pdf\n")
cat("  - K2vsK3_fullcohort_K3_KM.pdf\n")
cat("  - K2vsK3_idhwt_K2_KM.pdf\n")
cat("  - K2vsK3_idhwt_K3_KM.pdf\n")
cat("  - K2vsK3_idhwt_recluster_K2_KM.pdf\n")
cat("  - K2vsK3_idhwt_recluster_K3_KM.pdf\n")
cat("  - K2vsK3_idhwt_recluster_K4_KM.pdf\n")
cat("  - K2vsK3_comparison_panel.pdf\n")
cat("  - K2vsK3_clustering_metrics.pdf\n")
cat("  - K2vsK3_group_sizes.pdf\n")
cat("  - K2vsK3_decision_report.pdf\n")
cat("================================================================\n")
cat("ANALYSIS COMPLETE.\n")
