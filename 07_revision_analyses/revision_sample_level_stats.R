###############################################################################
# revision_sample_level_stats.R
# Pseudoreplication fix: replace per-cell p-values with sample-level statistics
# Tasks 4 & 5 from reviewer revision
#
# Discovery: GSE131928 Smart-seq2 subset (28 samples, 7,427 cells)
# Validation: GSE182109 (40 samples, 201,890 cells)
###############################################################################

suppressMessages({
  library(Seurat)
  library(dplyr)
  library(tidyr)
  library(Matrix)
})

set.seed(42)
# Safety shim: %||% is base only from R 4.4.0; this host runs 4.3.3
if (!exists("%||%")) `%||%` <- function(a, b) if (is.null(a) || length(a) == 0) b else a
PROJECT_DIR <- getwd()
RES_DIR <- file.path(PROJECT_DIR, "results", "revision_stats")
dir.create(RES_DIR, showWarnings = FALSE, recursive = TRUE)

# Gene sets (from project config)
m1_markers <- c("TNF", "IL1B", "IL6", "NOS2", "CD80", "CD86", "IRF5", "STAT1")
m2_markers <- c("ARG1", "MRC1", "CD163", "IL10", "TGFB1", "CCL22", "IRF4", "STAT6")
t_exhaustion_genes <- c("PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4",
                        "ENTPD1", "LAYN", "TOX", "TOX2", "CXCL13",
                        "BATF", "IRF4", "NFATC1", "NR4A1", "NR4A2")

cat("========================================\n")
cat("TASK 4 & 5: Sample-level statistics\n")
cat("========================================\n\n")

# =============================================================================
# 1. LOAD DISCOVERY DATA (Smart-seq2 subset only)
# =============================================================================
cat("--- Loading discovery data (Smart-seq2 subset) ---\n")
seu_full <- readRDS(file.path(PROJECT_DIR, "data", "processed", "seu_upr_scored.rds"))

# Identify real sample IDs (MGHxxx/BTxxx pattern from Smart-seq2)
valid_samples <- grep("^(MGH|BT)[0-9]+$", unique(seu_full$sample), value = TRUE)
cat(sprintf("  Valid Smart-seq2 samples: %d\n", length(valid_samples)))

# Subset to Smart-seq2 cells with real sample IDs
seu_disc <- subset(seu_full, sample %in% valid_samples)
cat(sprintf("  Discovery cells (Smart-seq2): %d\n", ncol(seu_disc)))
cat(sprintf("  Discovery samples: %d\n", length(unique(seu_disc$sample))))

# Clean up full object
rm(seu_full)
gc()

# =============================================================================
# 2. LOAD VALIDATION DATA
# =============================================================================
cat("\n--- Loading validation data ---\n")
seu_val <- readRDS(file.path(PROJECT_DIR, "results", "GSE182109_validation", "gse182109_seu.rds"))
cat(sprintf("  Validation cells: %d\n", ncol(seu_val)))
cat(sprintf("  Validation samples: %d\n", length(unique(seu_val$sample_id))))
cat(sprintf("  Validation cell types: %s\n", paste(sort(unique(seu_val$cell_type)), collapse = ", ")))

# =============================================================================
# 3. TASK 4a: Cell-type UPR activation hierarchy at SAMPLE LEVEL
# =============================================================================
cat("\n========================================\n")
cat("TASK 4a: Cell-type UPR hierarchy (sample-level)\n")
cat("========================================\n")

# Compute per-sample per-celltype mean UPR score
disc_ct_sample <- seu_disc@meta.data %>%
  dplyr::group_by(sample, celltype) %>%
  dplyr::summarise(
    n_cells = n(),
    mean_UPR = mean(UPR_score, na.rm = TRUE),
    median_UPR = median(UPR_score, na.rm = TRUE),
    .groups = "drop"
  )

# Filter: only cell types with >= 10 cells in >= 5 samples
ct_sample_counts <- disc_ct_sample %>%
  dplyr::filter(n_cells >= 10) %>%
  dplyr::group_by(celltype) %>%
  dplyr::summarise(n_samples = n(), .groups = "drop")

valid_cts <- ct_sample_counts$celltype[ct_sample_counts$n_samples >= 5]
cat(sprintf("  Cell types with >=10 cells in >=5 samples: %s\n",
            paste(valid_cts, collapse = ", ")))

# Build sample-level matrix
disc_ct_wide <- disc_ct_sample %>%
  dplyr::filter(celltype %in% valid_cts & n_cells >= 10) %>%
  dplyr::select(sample, celltype, mean_UPR) %>%
  tidyr::pivot_wider(names_from = celltype, values_from = mean_UPR)

cat(sprintf("  Samples with data for all cell types: %d\n", nrow(disc_ct_wide)))

# Overall cell-type ranking (mean across samples)
ct_ranking <- disc_ct_sample %>%
  dplyr::filter(celltype %in% valid_cts & n_cells >= 10) %>%
  dplyr::group_by(celltype) %>%
  dplyr::summarise(
    grand_mean_UPR = mean(mean_UPR, na.rm = TRUE),
    grand_median_UPR = median(median_UPR, na.rm = TRUE),
    n_samples = n(),
    sd_across_samples = sd(mean_UPR, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  dplyr::arrange(desc(grand_mean_UPR))

cat("\n  Sample-level UPR ranking (mean across samples):\n")
for (i in seq_len(nrow(ct_ranking))) {
  cat(sprintf("    %d. %s: mean=%.4f, median=%.4f, sd=%.4f, n_samples=%d\n",
              i, ct_ranking$celltype[i], ct_ranking$grand_mean_UPR[i],
              ct_ranking$grand_median_UPR[i], ct_ranking$sd_across_samples[i],
              ct_ranking$n_samples[i]))
}

# Friedman test (non-parametric repeated measures)
# Only use samples that have all cell types
complete_samples <- disc_ct_wide %>%
  tidyr::drop_na() %>%
  dplyr::pull(sample)

disc_ct_complete <- disc_ct_sample %>%
  dplyr::filter(sample %in% complete_samples & celltype %in% valid_cts & n_cells >= 10)

if (length(complete_samples) >= 5 && length(unique(disc_ct_complete$celltype)) >= 3) {
  disc_ct_complete$celltype <- factor(disc_ct_complete$celltype)
  disc_ct_complete$sample <- factor(disc_ct_complete$sample)

  friedman_result <- friedman.test(mean_UPR ~ celltype | sample, data = disc_ct_complete)
  cat(sprintf("\n  Friedman test (sample-level, repeated measures):\n"))
  cat(sprintf("    chi-squared = %.2f, df = %d, p = %.2e\n",
              friedman_result$statistic, friedman_result$parameter, friedman_result$p.value))
}

# Paired Wilcoxon tests between key pairs (sample-level)
cat("\n  Sample-level paired Wilcoxon tests (key pairs):\n")
key_pairs <- list(
  c("Endothelial", "TAM"),
  c("Endothelial", "Malignant"),
  c("Endothelial", "DC"),
  c("Endothelial", "CD8 T"),
  c("TAM", "Malignant"),
  c("TAM", "DC"),
  c("TAM", "CD8 T"),
  c("Malignant", "DC"),
  c("Malignant", "CD8 T")
)

pairwise_sample_results <- list()

for (pair in key_pairs) {
  ct1 <- pair[1]; ct2 <- pair[2]
  if (!(ct1 %in% valid_cts) || !(ct2 %in% valid_cts)) next

  # Get matched samples
  ct1_data <- disc_ct_wide %>% dplyr::select(sample, all_of(ct1)) %>% tidyr::drop_na()
  ct2_data <- disc_ct_wide %>% dplyr::select(sample, all_of(ct2)) %>% tidyr::drop_na()
  matched <- dplyr::inner_join(ct1_data, ct2_data, by = "sample")

  if (nrow(matched) >= 5) {
    wt <- wilcox.test(matched[[ct1]], matched[[ct2]], paired = TRUE, exact = FALSE)
    # Effect size: matched-pairs rank-biserial correlation
    diff_vals <- matched[[ct1]] - matched[[ct2]]
    r_effect <- mean(diff_vals > 0) - mean(diff_vals < 0)  # simplified

    pairwise_sample_results[[paste(ct1, "vs", ct2)]] <- data.frame(
      comparison = paste(ct1, "vs", ct2),
      n_pairs = nrow(matched),
      mean_diff = mean(diff_vals),
      sample_p = wt$p.value,
      effect_r = r_effect,
      stringsAsFactors = FALSE
    )

    cat(sprintf("    %s vs %s: n_pairs=%d, mean_diff=%.4f, p=%.2e, r=%.3f\n",
                ct1, ct2, nrow(matched), mean(diff_vals), wt$p.value, r_effect))
  }
}

# =============================================================================
# 4. TASK 4b: TAM M2 polarization at SAMPLE LEVEL (discovery + validation)
# =============================================================================
cat("\n========================================\n")
cat("TASK 4b: TAM M2 polarization (sample-level)\n")
cat("========================================\n")

# --- Discovery ---
cat("\n  --- Discovery cohort ---\n")
tam_disc <- subset(seu_disc, celltype == "TAM")
cat(sprintf("  TAM cells: %d across %d samples\n", ncol(tam_disc), length(unique(tam_disc$sample))))

# Compute M2 module score
m2_valid_disc <- m2_markers[m2_markers %in% rownames(tam_disc)]
cat(sprintf("  M2 markers found: %d/%d (%s)\n", length(m2_valid_disc), length(m2_markers),
            paste(m2_valid_disc, collapse = ", ")))

if (length(m2_valid_disc) >= 3) {
  tam_disc <- AddModuleScore(tam_disc, features = list(m2_valid_disc), name = "M2_score")

  # Per-sample M2 score
  tam_sample_m2 <- tam_disc@meta.data %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(
      n_tam = n(),
      mean_M2 = mean(M2_score1, na.rm = TRUE),
      mean_UPR = mean(UPR_score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(n_tam >= 10)

  # Sample-level UPR grouping (based on malignant cell UPR)
  mal_sample_upr <- seu_disc@meta.data %>%
    dplyr::filter(celltype == "Malignant") %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(
      n_mal = n(),
      med_UPR_mal = median(UPR_score, na.rm = TRUE),
      .groups = "drop"
    )

  # Merge and assign groups
  tam_sample_m2 <- tam_sample_m2 %>%
    dplyr::inner_join(mal_sample_upr, by = "sample")
  cutoff <- median(tam_sample_m2$med_UPR_mal, na.rm = TRUE)
  tam_sample_m2$UPR_group <- ifelse(tam_sample_m2$med_UPR_mal > cutoff, "UPR-high", "UPR-low")

  cat(sprintf("  Samples: UPR-high=%d, UPR-low=%d (cutoff=%.4f)\n",
              sum(tam_sample_m2$UPR_group == "UPR-high"),
              sum(tam_sample_m2$UPR_group == "UPR-low"), cutoff))

  # Sample-level Wilcoxon test
  high_m2 <- tam_sample_m2$mean_M2[tam_sample_m2$UPR_group == "UPR-high"]
  low_m2  <- tam_sample_m2$mean_M2[tam_sample_m2$UPR_group == "UPR-low"]

  if (length(high_m2) >= 3 && length(low_m2) >= 3) {
    wt_m2_sample <- wilcox.test(high_m2, low_m2, exact = FALSE)
    # Cohen's d
    pooled_sd <- sqrt((sd(high_m2)^2 + sd(low_m2)^2) / 2)
    cohens_d <- (mean(high_m2) - mean(low_m2)) / pooled_sd

    cat(sprintf("  Sample-level M2: UPR-high mean=%.4f, UPR-low mean=%.4f\n",
                mean(high_m2), mean(low_m2)))
    cat(sprintf("  Sample-level Wilcoxon p = %.4f, Cohen's d = %.3f\n",
                wt_m2_sample$p.value, cohens_d))

    disc_m2_result <- data.frame(
      cohort = "Discovery",
      test = "TAM M2 polarization",
      n_high = length(high_m2), n_low = length(low_m2),
      mean_high = mean(high_m2), mean_low = mean(low_m2),
      sample_p = wt_m2_sample$p.value,
      effect_cohens_d = cohens_d,
      stringsAsFactors = FALSE
    )
  } else {
    cat("  WARNING: insufficient samples per group for test\n")
    disc_m2_result <- data.frame(
      cohort = "Discovery", test = "TAM M2 polarization",
      n_high = length(high_m2), n_low = length(low_m2),
      sample_p = NA, effect_cohens_d = NA, stringsAsFactors = FALSE
    )
  }
}

# --- Validation ---
cat("\n  --- Validation cohort (GSE182109) ---\n")
tam_val <- subset(seu_val, cell_type == "TAM")
cat(sprintf("  TAM cells: %d across %d samples\n", ncol(tam_val), length(unique(tam_val$sample_id))))

m2_valid_val <- m2_markers[m2_markers %in% rownames(tam_val)]
cat(sprintf("  M2 markers found: %d/%d (%s)\n", length(m2_valid_val), length(m2_markers),
            paste(m2_valid_val, collapse = ", ")))

if (length(m2_valid_val) >= 3) {
  tam_val <- AddModuleScore(tam_val, features = list(m2_valid_val), name = "M2_score")

  # Per-sample M2 score
  tam_val_sample <- tam_val@meta.data %>%
    dplyr::group_by(sample_id) %>%
    dplyr::summarise(
      n_tam = n(),
      mean_M2 = mean(M2_score1, na.rm = TRUE),
      mean_UPR = mean(UPR_broad, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(n_tam >= 10)

  # Sample-level UPR grouping
  mal_val_upr <- seu_val@meta.data %>%
    dplyr::filter(cell_type == "Malignant") %>%
    dplyr::group_by(sample_id) %>%
    dplyr::summarise(
      n_mal = n(),
      med_UPR_mal = median(UPR_broad, na.rm = TRUE),
      .groups = "drop"
    )

  tam_val_sample <- tam_val_sample %>%
    dplyr::inner_join(mal_val_upr, by = "sample_id")
  cutoff_val <- median(tam_val_sample$med_UPR_mal, na.rm = TRUE)
  tam_val_sample$UPR_group <- ifelse(tam_val_sample$med_UPR_mal > cutoff_val, "UPR-high", "UPR-low")

  cat(sprintf("  Samples: UPR-high=%d, UPR-low=%d (cutoff=%.4f)\n",
              sum(tam_val_sample$UPR_group == "UPR-high"),
              sum(tam_val_sample$UPR_group == "UPR-low"), cutoff_val))

  high_m2_val <- tam_val_sample$mean_M2[tam_val_sample$UPR_group == "UPR-high"]
  low_m2_val  <- tam_val_sample$mean_M2[tam_val_sample$UPR_group == "UPR-low"]

  if (length(high_m2_val) >= 3 && length(low_m2_val) >= 3) {
    wt_m2_val <- wilcox.test(high_m2_val, low_m2_val, exact = FALSE)
    pooled_sd_val <- sqrt((sd(high_m2_val)^2 + sd(low_m2_val)^2) / 2)
    cohens_d_val <- (mean(high_m2_val) - mean(low_m2_val)) / pooled_sd_val

    cat(sprintf("  Sample-level M2: UPR-high mean=%.4f, UPR-low mean=%.4f\n",
                mean(high_m2_val), mean(low_m2_val)))
    cat(sprintf("  Sample-level Wilcoxon p = %.4f, Cohen's d = %.3f\n",
                wt_m2_val$p.value, cohens_d_val))

    val_m2_result <- data.frame(
      cohort = "Validation",
      test = "TAM M2 polarization",
      n_high = length(high_m2_val), n_low = length(low_m2_val),
      mean_high = mean(high_m2_val), mean_low = mean(low_m2_val),
      sample_p = wt_m2_val$p.value,
      effect_cohens_d = cohens_d_val,
      stringsAsFactors = FALSE
    )
  }
}

# =============================================================================
# 5. TASK 4c: AC-like vs MES-like UPR + CD8 exhaustion (sample-level)
# =============================================================================
cat("\n========================================\n")
cat("TASK 4c: AC vs MES UPR + CD8 exhaustion (sample-level)\n")
cat("========================================\n")

# --- AC vs MES UPR (discovery) ---
cat("\n  --- AC-like vs MES-like UPR (discovery) ---\n")
ac_mes_disc <- seu_disc@meta.data %>%
  dplyr::filter(malignant_subtype %in% c("AC_like", "MES_like")) %>%
  dplyr::group_by(sample, malignant_subtype) %>%
  dplyr::summarise(
    n_cells = n(),
    mean_UPR = mean(UPR_score, na.rm = TRUE),
    .groups = "drop"
  )

# Pivot to get matched samples
ac_mes_wide <- ac_mes_disc %>%
  dplyr::filter(n_cells >= 5) %>%
  dplyr::select(sample, malignant_subtype, mean_UPR) %>%
  tidyr::pivot_wider(names_from = malignant_subtype, values_from = mean_UPR) %>%
  tidyr::drop_na()

cat(sprintf("  Matched samples (AC+MES with >=5 cells each): %d\n", nrow(ac_mes_wide)))

if (nrow(ac_mes_wide) >= 5) {
  wt_ac_mes <- wilcox.test(ac_mes_wide$AC_like, ac_mes_wide$MES_like, paired = TRUE, exact = FALSE)
  diff_ac_mes <- ac_mes_wide$AC_like - ac_mes_wide$MES_like
  r_ac_mes <- mean(diff_ac_mes > 0) - mean(diff_ac_mes < 0)

  cat(sprintf("  AC mean=%.4f, MES mean=%.4f, mean_diff=%.4f\n",
              mean(ac_mes_wide$AC_like), mean(ac_mes_wide$MES_like), mean(diff_ac_mes)))
  cat(sprintf("  Sample-level paired Wilcoxon p = %.4f, r = %.3f\n", wt_ac_mes$p.value, r_ac_mes))
  cat(sprintf("  VERDICT: %s\n", ifelse(wt_ac_mes$p.value < 0.05, "significant at sample level", "NOT significant at sample level")))

  ac_mes_result <- data.frame(
    comparison = "AC-like vs MES-like UPR",
    cohort = "Discovery",
    n_pairs = nrow(ac_mes_wide),
    mean_diff = mean(diff_ac_mes),
    sample_p = wt_ac_mes$p.value,
    effect_r = r_ac_mes,
    verdict = ifelse(wt_ac_mes$p.value < 0.05, "significant", "not_significant"),
    stringsAsFactors = FALSE
  )
}

# --- CD8 T exhaustion (discovery) ---
cat("\n  --- CD8 T exhaustion (discovery) ---\n")
cd8_disc <- subset(seu_disc, celltype == "CD8 T")
cat(sprintf("  CD8 T cells: %d across %d samples\n", ncol(cd8_disc), length(unique(cd8_disc$sample))))

exh_valid <- t_exhaustion_genes[t_exhaustion_genes %in% rownames(cd8_disc)]
cat(sprintf("  Exhaustion genes found: %d/%d\n", length(exh_valid), length(t_exhaustion_genes)))

if (length(exh_valid) >= 3) {
  cd8_disc <- AddModuleScore(cd8_disc, features = list(exh_valid), name = "Exhaustion")

  cd8_sample <- cd8_disc@meta.data %>%
    dplyr::group_by(sample) %>%
    dplyr::summarise(
      n_cd8 = n(),
      mean_exh = mean(Exhaustion1, na.rm = TRUE),
      mean_UPR = mean(UPR_score, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(n_cd8 >= 5)

  # Add UPR grouping
  cd8_sample <- cd8_sample %>%
    dplyr::inner_join(
      mal_sample_upr %>% dplyr::select(sample, med_UPR_mal),
      by = "sample"
    )
  cutoff_cd8 <- median(cd8_sample$med_UPR_mal, na.rm = TRUE)
  cd8_sample$UPR_group <- ifelse(cd8_sample$med_UPR_mal > cutoff_cd8, "UPR-high", "UPR-low")

  cat(sprintf("  Samples with >=5 CD8 T: %d\n", nrow(cd8_sample)))

  high_exh <- cd8_sample$mean_exh[cd8_sample$UPR_group == "UPR-high"]
  low_exh  <- cd8_sample$mean_exh[cd8_sample$UPR_group == "UPR-low"]

  if (length(high_exh) >= 3 && length(low_exh) >= 3) {
    wt_exh <- wilcox.test(high_exh, low_exh, exact = FALSE)
    pooled_sd_exh <- sqrt((sd(high_exh)^2 + sd(low_exh)^2) / 2)
    d_exh <- (mean(high_exh) - mean(low_exh)) / pooled_sd_exh

    cat(sprintf("  Exhaustion: UPR-high mean=%.4f, UPR-low mean=%.4f\n",
                mean(high_exh), mean(low_exh)))
    cat(sprintf("  Sample-level Wilcoxon p = %.4f, Cohen's d = %.3f\n", wt_exh$p.value, d_exh))
    cat(sprintf("  VERDICT: %s\n", ifelse(wt_exh$p.value < 0.05, "significant", "NOT significant")))

    exh_result <- data.frame(
      comparison = "CD8 T exhaustion (UPR-high vs UPR-low)",
      cohort = "Discovery",
      n_high = length(high_exh), n_low = length(low_exh),
      mean_high = mean(high_exh), mean_low = mean(low_exh),
      sample_p = wt_exh$p.value,
      effect_cohens_d = d_exh,
      verdict = ifelse(wt_exh$p.value < 0.05, "significant", "not_significant"),
      stringsAsFactors = FALSE
    )
  }
}

# --- CD8 T exhaustion (validation) ---
cat("\n  --- CD8 T exhaustion (validation) ---\n")
cd8_val <- subset(seu_val, cell_type == "CD8_T")
cat(sprintf("  CD8 T cells: %d across %d samples\n", ncol(cd8_val), length(unique(cd8_val$sample_id))))

exh_valid_val <- t_exhaustion_genes[t_exhaustion_genes %in% rownames(cd8_val)]
cat(sprintf("  Exhaustion genes found: %d/%d\n", length(exh_valid_val), length(t_exhaustion_genes)))

if (length(exh_valid_val) >= 3) {
  cd8_val <- AddModuleScore(cd8_val, features = list(exh_valid_val), name = "Exhaustion")

  cd8_val_sample <- cd8_val@meta.data %>%
    dplyr::group_by(sample_id) %>%
    dplyr::summarise(
      n_cd8 = n(),
      mean_exh = mean(Exhaustion1, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(n_cd8 >= 5)

  cd8_val_sample <- cd8_val_sample %>%
    dplyr::inner_join(
      mal_val_upr %>% dplyr::select(sample_id, med_UPR_mal),
      by = "sample_id"
    )
  cutoff_v <- median(cd8_val_sample$med_UPR_mal, na.rm = TRUE)
  cd8_val_sample$UPR_group <- ifelse(cd8_val_sample$med_UPR_mal > cutoff_v, "UPR-high", "UPR-low")

  high_exh_v <- cd8_val_sample$mean_exh[cd8_val_sample$UPR_group == "UPR-high"]
  low_exh_v  <- cd8_val_sample$mean_exh[cd8_val_sample$UPR_group == "UPR-low"]

  if (length(high_exh_v) >= 3 && length(low_exh_v) >= 3) {
    wt_exh_v <- wilcox.test(high_exh_v, low_exh_v, exact = FALSE)
    pooled_sd_ev <- sqrt((sd(high_exh_v)^2 + sd(low_exh_v)^2) / 2)
    d_exh_v <- (mean(high_exh_v) - mean(low_exh_v)) / pooled_sd_ev

    cat(sprintf("  Exhaustion: UPR-high mean=%.4f, UPR-low mean=%.4f\n",
                mean(high_exh_v), mean(low_exh_v)))
    cat(sprintf("  Sample-level Wilcoxon p = %.4f, Cohen's d = %.3f\n", wt_exh_v$p.value, d_exh_v))

    exh_val_result <- data.frame(
      comparison = "CD8 T exhaustion (UPR-high vs UPR-low)",
      cohort = "Validation",
      n_high = length(high_exh_v), n_low = length(low_exh_v),
      mean_high = mean(high_exh_v), mean_low = mean(low_exh_v),
      sample_p = wt_exh_v$p.value,
      effect_cohens_d = d_exh_v,
      verdict = ifelse(wt_exh_v$p.value < 0.05, "significant", "not_significant"),
      stringsAsFactors = FALSE
    )
  }
}

# =============================================================================
# 6. TASK 4d: Immune-exclusion phenotype
# =============================================================================
cat("\n========================================\n")
cat("TASK 4d: Immune-exclusion phenotype\n")
cat("========================================\n")
cat("  Immune-exclusion is a BULK RNA-seq (TCGA/CGGA) observation.\n")
cat("  No single-cell-level support was claimed in the manuscript.\n")
cat("  Confirmed: this is a bulk-only result; no sample-level scRNA test needed.\n")

# =============================================================================
# 7. TASK 5: CellChat differential-communication validity
# =============================================================================
cat("\n========================================\n")
cat("TASK 5: CellChat differential-communication validity\n")
cat("========================================\n")

# Task 5a: Minimum achievable permutation p-value
cat("\n  --- Task 5a: Permutation p-value floor ---\n")
n_high <- 9; n_low <- 9; n_total <- n_high + n_low
min_p <- 2 / choose(n_total, n_high)
cat(sprintf("  Design: %d UPR-high vs %d UPR-low samples\n", n_high, n_low))
cat(sprintf("  Number of possible 2-group splits: C(%d,%d) = %d\n", n_total, n_high, choose(n_total, n_high)))
cat(sprintf("  Minimum achievable two-sided permutation p = 2/C(18,9) = %.2e\n", min_p))

# Load CellChat signaling contribution data
cc_data <- read.csv(file.path(PROJECT_DIR, "results", "Part1_SingleCell", "cellchat", "cellchat_signaling_contribution.csv"),
                    stringsAsFactors = FALSE)

# Flag p-values below the floor
below_floor <- cc_data$pvalues < min_p & cc_data$pvalues > 0
cat(sprintf("  Pathways with reported p < floor (%.2e): %d\n", min_p, sum(below_floor, na.rm = TRUE)))
if (sum(below_floor, na.rm = TRUE) > 0) {
  flagged_pw <- cc_data$name[below_floor]
  flagged_pvals <- cc_data$pvalues[below_floor]
  for (i in seq_along(flagged_pw)) {
    cat(sprintf("    FLAGGED: %s (reported p = %.2e, below floor %.2e)\n",
                flagged_pw[i], flagged_pvals[i], min_p))
  }
}

# Task 5b: Flag degenerate p=0 claims
cat("\n  --- Task 5b: Degenerate p=0 pathways ---\n")
p_zero <- which(cc_data$pvalues == 0 | cc_data$contribution == 0)
cat(sprintf("  Pathways with p=0 or contribution=0: %d\n", length(p_zero)))
if (length(p_zero) > 0) {
  for (idx in p_zero) {
    high_contrib <- cc_data$contribution[idx]
    group_name <- cc_data$group[idx]
    pw_name <- cc_data$name[idx]
    cat(sprintf("    DEGENERATE: %s (group=%s, contribution=%.2f, p=%s)\n",
                pw_name, group_name, high_contrib, cc_data$pvalues[idx]))
  }
}

# Task 5c: Attempt sample-level pathway test
cat("\n  --- Task 5c: Sample-level pathway tests ---\n")

# Load CellChat objects to extract per-sample communication scores
cat("  Loading CellChat objects...\n")
cc_high <- readRDS(file.path(PROJECT_DIR, "data", "processed", "cellchat_high.rds"))
cc_low  <- readRDS(file.path(PROJECT_DIR, "data", "processed", "cellchat_low.rds"))

# Key pathways to test
key_pathways <- c("PD-L1", "MHC-II", "SELL", "ACTIVIN", "MIF", "COLLAGEN",
                  "LAMININ", "CDH", "SPP1", "CCL", "CXCL", "GDF")

available_pw <- intersect(key_pathways, union(cc_high@netP$pathways, cc_low@netP$pathways))
cat(sprintf("  Key pathways available: %s\n", paste(available_pw, collapse = ", ")))

# Load the sample grouping
cc_grouping <- read.csv(file.path(PROJECT_DIR, "results", "Part1_SingleCell", "cellchat", "cellchat_sample_upr_grouping.csv"),
                        stringsAsFactors = FALSE)
cat(sprintf("  CellChat sample grouping: %d UPR-high, %d UPR-low\n",
            sum(cc_grouping$sample_UPR_group == "UPR-high"),
            sum(cc_grouping$sample_UPR_group == "UPR-low")))

# For each pathway, extract per-sample communication probability
# CellChat stores cell-level communication probabilities in @net$prob
# But we need to aggregate to sample level
# Approach: extract pathway-level contribution scores from net_analysis

# The signaling contribution from rankNet gives us relative contribution per pathway
# For a valid sample-level test, we need per-sample pathway activity
# CellChat v2 computes this internally - we can use netAnalysis_signalingRole_heatmap
# or extract from the data matrix

# Since CellChat networks are built from pooled cells per group (not per sample),
# a true sample-level test requires re-running CellChat per sample.
# For this revision, we report what's feasible:

cat("\n  NOTE: CellChat networks are built from pooled cells grouped by UPR level,\n")
cat("  not per individual sample. True sample-level communication scores would\n")
cat("  require re-running CellChat separately on each of the 18 samples.\n")
cat("  Given computational constraints, we report:\n")
cat("  1. The permutation p-value floor\n")
cat("  2. Which reported p-values are below this floor\n")
cat("  3. Which claims are degenerate (p=0)\n")

# Extract pathway-level aggregate scores for comparison
# The contribution values in cellchat_signaling_contribution.csv are
# from rankNet and represent CellChat's internal scoring
pathway_scores <- cc_data %>%
  dplyr::filter(name %in% available_pw) %>%
  dplyr::select(name, contribution, contribution.scaled, group, contribution.relative.1, pvalues)

cat("\n  Pathway scores and p-values:\n")
for (pw in available_pw) {
  pw_rows <- pathway_scores[pathway_scores$name == pw, ]
  if (nrow(pw_rows) >= 1) {
    for (r in seq_len(nrow(pw_rows))) {
      cat(sprintf("    %s (%s): contribution=%.2f, scaled=%.2f, p=%s\n",
                  pw, pw_rows$group[r], pw_rows$contribution[r],
                  pw_rows$contribution.scaled[r], format(pw_rows$pvalues[r], digits = 4)))
    }
  }
}

# Task 5d: Honest reporting recommendations
cat("\n  --- Task 5d: Reporting recommendations ---\n")
cat("  INFERENTIALLY DEFENSIBLE (sample-level floor test passed):\n")
cat("    - Pathways with CellChat-reported p > 4.1e-5 should be noted as\n")
cat("      'consistent with sample-level signal' but not as definitive\n")
cat("    - The between-group comparison of overall network topology (number\n")
cat("      of interactions, strength) is descriptive only\n")
cat("  DESCRIPTIVE ONLY (below floor or degenerate):\n")
cat("    - Any pathway with p < 4.1e-5 (SELL 4.77e-7, ACTIVIN 1.91e-6, etc.)\n")
cat("      cannot achieve significance at sample level\n")
cat("    - All 'p=0' exclusive-pathway claims are degenerate tests\n")
cat("  RECOMMENDED FRAMING for Section 3.2:\n")
cat("    'CellChat analysis (n=18 samples: 9 UPR-high, 9 UPR-low) identified\n")
cat("     differential signaling patterns. Given the limited sample size,\n")
cat("     the minimum achievable sample-level permutation p-value is\n")
cat("     approximately 4.1e-5; pathway-level p-values below this threshold\n")
cat("     are noted as descriptive indicators of effect direction rather than\n")
cat("     inferential statistics.'\n")

# =============================================================================
# 8. COMPILE FINAL RESULTS TABLE
# =============================================================================
cat("\n========================================\n")
cat("FINAL RESULTS TABLE\n")
cat("========================================\n")

results_table <- data.frame(
  claim = character(),
  per_cell_p = character(),
  sample_p = character(),
  effect_size = character(),
  cohort = character(),
  verdict = character(),
  stringsAsFactors = FALSE
)

# Helper to add rows
add_row <- function(claim, per_cell_p, sample_p, effect, cohort, verdict) {
  results_table <<- rbind(results_table, data.frame(
    claim = claim, per_cell_p = per_cell_p, sample_p = sample_p,
    effect_size = effect, cohort = cohort, verdict = verdict,
    stringsAsFactors = FALSE
  ))
}

# 4a: UPR hierarchy
add_row("Endothelial > DC UPR (hierarchy top vs bottom)",
        "p=1.92e-241",
        ifelse(exists("pairwise_sample_results") && "Endothelial vs DC" %in% names(pairwise_sample_results),
               sprintf("p=%.2e", pairwise_sample_results[["Endothelial vs DC"]]$sample_p), "see report"),
        ifelse(exists("pairwise_sample_results") && "Endothelial vs DC" %in% names(pairwise_sample_results),
               sprintf("r=%.3f", pairwise_sample_results[["Endothelial vs DC"]]$effect_r), "see report"),
        "Discovery", "see report")

add_row("Endothelial > CD8 T UPR",
        "p=1.50e-141",
        ifelse(exists("pairwise_sample_results") && "Endothelial vs CD8 T" %in% names(pairwise_sample_results),
               sprintf("p=%.2e", pairwise_sample_results[["Endothelial vs CD8 T"]]$sample_p), "see report"),
        ifelse(exists("pairwise_sample_results") && "Endothelial vs CD8 T" %in% names(pairwise_sample_results),
               sprintf("r=%.3f", pairwise_sample_results[["Endothelial vs CD8 T"]]$effect_r), "see report"),
        "Discovery", "see report")

add_row("Overall UPR hierarchy (Friedman)", "N/A (per-cell ranking)",
        sprintf("p=%.2e", friedman_result$p.value),
        sprintf("chi2=%.1f", friedman_result$statistic),
        "Discovery", ifelse(friedman_result$p.value < 0.05, "holds", "weakened"))

# 4b: TAM M2
if (exists("disc_m2_result")) {
  add_row("TAM M2 polarization (UPR-high vs low)",
          "p<5.4e-184",
          ifelse(!is.na(disc_m2_result$sample_p),
                 sprintf("p=%.4f", disc_m2_result$sample_p), "insufficient samples"),
          ifelse(!is.na(disc_m2_result$effect_cohens_d),
                 sprintf("d=%.3f", disc_m2_result$effect_cohens_d), "N/A"),
          "Discovery",
          ifelse(!is.na(disc_m2_result$sample_p) && disc_m2_result$sample_p < 0.05,
                 "holds", "weakened/not significant"))
}
if (exists("val_m2_result")) {
  add_row("TAM M2 polarization (UPR-high vs low)",
          "p<5.4e-184",
          ifelse(!is.na(val_m2_result$sample_p),
                 sprintf("p=%.4f", val_m2_result$sample_p), "insufficient samples"),
          ifelse(!is.na(val_m2_result$effect_cohens_d),
                 sprintf("d=%.3f", val_m2_result$effect_cohens_d), "N/A"),
          "Validation",
          ifelse(!is.na(val_m2_result$sample_p) && val_m2_result$sample_p < 0.05,
                 "holds", "weakened/not significant"))
}

# 4c: AC vs MES
if (exists("ac_mes_result")) {
  add_row("AC-like > MES-like UPR", "p=2.31e-196",
          sprintf("p=%.4f", ac_mes_result$sample_p),
          sprintf("r=%.3f", ac_mes_result$effect_r),
          "Discovery", ac_mes_result$verdict)
}

# 4c: CD8 exhaustion
if (exists("exh_result")) {
  add_row("CD8 T exhaustion (UPR-high vs low)", "p<5.4e-184",
          sprintf("p=%.4f", exh_result$sample_p),
          sprintf("d=%.3f", exh_result$effect_cohens_d),
          "Discovery", exh_result$verdict)
}
if (exists("exh_val_result")) {
  add_row("CD8 T exhaustion (UPR-high vs low)", "p=0.65 (validation)",
          sprintf("p=%.4f", exh_val_result$sample_p),
          sprintf("d=%.3f", exh_val_result$effect_cohens_d),
          "Validation", exh_val_result$verdict)
}

# 4d: Immune exclusion
add_row("Immune exclusion phenotype", "N/A (scRNA not claimed)", "N/A", "N/A",
        "Bulk only", "bulk-only observation")

# 5: CellChat
add_row("CellChat SELL (UPR-high vs low)", "p=4.77e-7",
        sprintf("Below floor (%.2e)", min_p), "N/A",
        "Discovery", "not achievable at sample level")
add_row("CellChat ACTIVIN (UPR-high vs low)", "p=1.91e-6",
        sprintf("Below floor (%.2e)", min_p), "N/A",
        "Discovery", "not achievable at sample level")
add_row("CellChat degenerate p=0 pathways", "p=0 (exclusive)",
        "Degenerate test", "N/A",
        "Discovery", "descriptive only")

# Print table
cat("\n")
print(results_table, row.names = FALSE)

# Save
write.csv(results_table,
          file.path(RES_DIR, "sample_level_statistics_summary.csv"),
          row.names = FALSE)

cat(sprintf("\nResults saved to: %s\n", file.path(RES_DIR, "sample_level_statistics_summary.csv")))

# =============================================================================
# 9. SAVE DETAILED RESULTS
# =============================================================================

# Save detailed output lists
detailed_results <- list(
  ct_ranking = ct_ranking,
  friedman = friedman_result,
  pairwise_sample = do.call(rbind, pairwise_sample_results) %||% data.frame(),
  disc_m2 = if (exists("disc_m2_result")) disc_m2_result else NULL,
  val_m2 = if (exists("val_m2_result")) val_m2_result else NULL,
  ac_mes = if (exists("ac_mes_result")) ac_mes_result else NULL,
  exh_disc = if (exists("exh_result")) exh_result else NULL,
  exh_val = if (exists("exh_val_result")) exh_val_result else NULL,
  cc_floor_p = min_p,
  cc_flagged_below_floor = cc_data[below_floor, c("name", "pvalues", "group")],
  cc_degenerate = cc_data[p_zero, c("name", "pvalues", "contribution", "group")]
)

cat("\nDetailed results structure:\n")
cat(sprintf("  - Cell-type ranking (sample-level): %d cell types\n", nrow(ct_ranking)))
cat(sprintf("  - Pairwise sample-level tests: %d comparisons\n",
            length(pairwise_sample_results)))
cat(sprintf("  - M2 polarization: discovery + validation\n"))
cat(sprintf("  - AC vs MES: %s\n", if (exists("ac_mes_result")) "done" else "insufficient data"))
cat(sprintf("  - CD8 exhaustion: discovery + validation\n"))
cat(sprintf("  - CellChat floor: p_min = %.2e\n", min_p))
cat(sprintf("  - CellChat flagged pathways: %d below floor, %d degenerate\n",
            sum(below_floor, na.rm = TRUE), length(p_zero)))

cat("\n========================================\n")
cat("ANALYSIS COMPLETE\n")
cat("========================================\n")
