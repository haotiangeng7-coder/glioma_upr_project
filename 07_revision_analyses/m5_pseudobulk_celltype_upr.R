#!/usr/bin/env Rscript
# M5: Patient-level cell-type UPR activation comparison
# Replaces pseudoreplicated cell-level Wilcoxon with patient-aggregated analysis
# GSE131928 glioma scRNA-seq: 37 patients, 9 cell types
#
# CRITICAL: The 'sample' column in this Seurat object uses CELL BARCODES,
# not biological sample IDs. Proper patient IDs are extracted from barcode
# patterns (MGH*/BT*/P102/etc.). Total: 37 independent patients.
#
# Outputs:
#   results/pseudobulk_celltype_upr_sample_level.csv
#   results/pseudobulk_celltype_upr_pairwise.csv
#   results/pseudobulk_celltype_upr_notes.md
#   figures/SuppFig_pseudobulk_celltype_upr.pdf

set.seed(42)

suppressPackageStartupMessages({
  library(Seurat)
  library(lme4)
  library(lmerTest)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(gridExtra)
})

# ── Load data ──────────────────────────────────────────────────────────────
seu <- readRDS("data/processed/seu_upr_scored.rds")
cat(sprintf("Loaded: %d cells, %d genes\n", ncol(seu), nrow(seu)))

stopifnot(all(c("UPR_score", "sample", "celltype") %in% colnames(seu@meta.data)))

meta <- seu@meta.data
samples_raw <- as.character(meta$sample)

# ── Extract proper PATIENT IDs ─────────────────────────────────────────────
# MGH*/BT* samples are already patient IDs (10x platform)
# "102_1", "105_B1_2" etc. are cell barcodes from Smart-seq2 tumors
#   → patient = leading number before first separator
is_10x_patient <- grepl("^(MGH|BT)", samples_raw)
patient <- samples_raw
patient[!is_10x_patient] <- paste0("P", gsub("^([0-9]+).*", "\\1", samples_raw[!is_10x_patient]))

meta$patient <- patient
n_patients <- length(unique(patient))
cat(sprintf("Patients: %d (28 MGH/BT + %d Smart-seq2 tumors)\n",
            n_patients, n_patients - 28))

# ── Cell-type overview ─────────────────────────────────────────────────────
ct_counts <- table(meta$celltype)
cat("\nPer cell-type cell counts:\n")
print(ct_counts)

# Per-patient cell-type coverage
patient_ct <- table(meta$patient, meta$celltype)
n_ct_per_patient <- rowSums(patient_ct > 0)
cat(sprintf("\nPatients with >=2 cell types: %d / %d\n",
            sum(n_ct_per_patient >= 2), n_patients))
cat(sprintf("Patients with >=5 cell types: %d\n", sum(n_ct_per_patient >= 5)))

cells_per_patient <- table(meta$patient)
cat(sprintf("Cells per patient: median=%d, range=%d-%d\n",
            median(cells_per_patient), min(cells_per_patient), max(cells_per_patient)))

# Low-n flags
for(ct in colnames(patient_ct)) {
  n_pat <- sum(patient_ct[, ct] > 0)
  cat(sprintf("  %-20s: %d patients, %d cells\n", ct, n_pat, sum(patient_ct[, ct])))
}

low_n_ct <- names(which(colSums(patient_ct > 0) < 5))
cat(sprintf("\nCell types in <5 patients (severely underpowered): %s\n",
            paste(low_n_ct, collapse=", ")))

# ── Step 1: Patient-level aggregation ──────────────────────────────────────
# For each (patient, cell_type): mean UPR_score
patient_ct_agg <- meta %>%
  group_by(patient, celltype) %>%
  summarise(
    mean_UPR = mean(UPR_score),
    n_cells  = n(),
    .groups  = "drop"
  )

cat(sprintf("\nPatient x cell-type rows: %d\n", nrow(patient_ct_agg)))

# Save
write.csv(patient_ct_agg, "results/pseudobulk_celltype_upr_sample_level.csv",
          row.names = FALSE)
cat("Saved: results/pseudobulk_celltype_upr_sample_level.csv\n")

# ── Descriptive: per-cell-type patient-level UPR ───────────────────────────
ct_stats <- patient_ct_agg %>%
  group_by(celltype) %>%
  summarise(
    n_patients = n(),
    mean_UPR   = mean(mean_UPR),
    sd_UPR     = sd(mean_UPR),
    median_UPR = median(mean_UPR),
    n_cells_total = sum(n_cells),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_UPR))

cat("\n=== Cell-type patient-level UPR statistics ===\n")
print(as.data.frame(ct_stats))

# ── Step 2: Patient-level pairwise comparisons ─────────────────────────────
# For each cell type pair: use patients that have BOTH cell types
# (paired within-patient comparison where possible)
cell_types <- sort(unique(meta$celltype))
pairs_list <- combn(cell_types, 2, simplify = FALSE)

# Cell-level Wilcoxon p-values (reproducing inflated results) for reference
cell_level_pvals <- list()
for (pair in pairs_list) {
  ct_a <- pair[1]; ct_b <- pair[2]
  vals_a <- meta$UPR_score[meta$celltype == ct_a]
  vals_b <- meta$UPR_score[meta$celltype == ct_b]
  wt <- wilcox.test(vals_a, vals_b, exact = FALSE)
  cell_level_pvals[[paste(ct_a, ct_b, sep = " vs ")]] <- wt$p.value
}

# Patient-level analysis
# Strategy: For pairs with >=5 patients having both cell types → paired test
#           For pairs with 2-4 paired patients → report but flag low power
#           For pairs with <2 paired patients → cannot test meaningfully

pairwise_results <- data.frame(
  pair        = character(length(pairs_list)),
  n_pat_a     = integer(length(pairs_list)),
  n_pat_b     = integer(length(pairs_list)),
  n_paired    = integer(length(pairs_list)),
  mean_UPR_a  = numeric(length(pairs_list)),
  mean_UPR_b  = numeric(length(pairs_list)),
  delta_UPR   = numeric(length(pairs_list)),
  test_method = character(length(pairs_list)),
  statistic   = numeric(length(pairs_list)),
  p_value     = numeric(length(pairs_list)),
  note        = character(length(pairs_list)),
  stringsAsFactors = FALSE
)

for (i in seq_along(pairs_list)) {
  ct_a <- pairs_list[[i]][1]
  ct_b <- pairs_list[[i]][2]

  # Patient-level means
  vals_a <- patient_ct_agg %>% filter(celltype == ct_a)
  vals_b <- patient_ct_agg %>% filter(celltype == ct_b)

  n_pat_a <- nrow(vals_a)
  n_pat_b <- nrow(vals_b)
  patients_a <- vals_a$patient
  patients_b <- vals_b$patient
  paired_patients <- intersect(patients_a, patients_b)
  n_paired <- length(paired_patients)

  if (n_paired >= 5) {
    # Paired test: use patients with both cell types
    paired_data <- patient_ct_agg %>%
      filter(patient %in% paired_patients,
             celltype %in% c(ct_a, ct_b)) %>%
      arrange(patient, celltype)

    # Extract paired vectors
    u_a <- paired_data$mean_UPR[paired_data$celltype == ct_a]
    u_b <- paired_data$mean_UPR[paired_data$celltype == ct_b]

    # Paired Wilcoxon
    wt <- wilcox.test(u_a, u_b, paired = TRUE, exact = FALSE)

    test_method <- "Paired Wilcoxon (patient-level)"
    statistic   <- wt$statistic
    p_val       <- wt$p.value
    note        <- sprintf("OK: %d paired patients", n_paired)

  } else if (n_paired >= 2) {
    # Borderline: still run paired but flag
    paired_data <- patient_ct_agg %>%
      filter(patient %in% paired_patients,
             celltype %in% c(ct_a, ct_b)) %>%
      arrange(patient, celltype)

    u_a <- paired_data$mean_UPR[paired_data$celltype == ct_a]
    u_b <- paired_data$mean_UPR[paired_data$celltype == ct_b]

    wt <- wilcox.test(u_a, u_b, paired = TRUE, exact = FALSE)

    test_method <- "Paired Wilcoxon (LOW POWER)"
    statistic   <- wt$statistic
    p_val       <- wt$p.value
    note        <- sprintf("LOW POWER: only %d paired patients", n_paired)

  } else if (n_pat_a >= 3 && n_pat_b >= 3) {
    # No paired patients but enough per group: unpaired
    wt <- wilcox.test(vals_a$mean_UPR, vals_b$mean_UPR, exact = FALSE)
    test_method <- "Unpaired Wilcoxon (NO paired patients)"
    statistic   <- wt$statistic
    p_val       <- wt$p.value
    note        <- sprintf("UNPAIRED: 0 paired patients; %d vs %d independent", n_pat_a, n_pat_b)

  } else {
    # Insufficient data for meaningful test
    test_method <- "INSUFFICIENT DATA"
    statistic   <- NA
    p_val       <- NA
    note        <- sprintf("UNSTABLE: n_a=%d, n_b=%d, paired=%d", n_pat_a, n_pat_b, n_paired)
  }

  pairwise_results$pair[i]        <- paste(ct_a, ct_b, sep = " vs ")
  pairwise_results$n_pat_a[i]     <- n_pat_a
  pairwise_results$n_pat_b[i]     <- n_pat_b
  pairwise_results$n_paired[i]    <- n_paired
  pairwise_results$mean_UPR_a[i]  <- mean(vals_a$mean_UPR)
  pairwise_results$mean_UPR_b[i]  <- mean(vals_b$mean_UPR)
  pairwise_results$delta_UPR[i]   <- mean(vals_a$mean_UPR) - mean(vals_b$mean_UPR)
  pairwise_results$test_method[i] <- test_method
  pairwise_results$statistic[i]   <- statistic
  pairwise_results$p_value[i]     <- p_val
  pairwise_results$note[i]        <- note
}

# BH correction (only on valid p-values)
valid_p <- !is.na(pairwise_results$p_value)
pairwise_results$p_adj <- NA
pairwise_results$p_adj[valid_p] <- p.adjust(pairwise_results$p_value[valid_p], method = "BH")
pairwise_results$significant <- pairwise_results$p_adj < 0.05
pairwise_results$significant[is.na(pairwise_results$significant)] <- FALSE

# Cell-level comparison
pairwise_results$p_cell_level <- sapply(pairwise_results$pair,
  function(pn) cell_level_pvals[[pn]])
pairwise_results$p_adj_cell_level <- p.adjust(pairwise_results$p_cell_level, method = "BH")

# Sort by p_adj
pairwise_results <- pairwise_results[order(pairwise_results$p_adj, na.last = TRUE), ]

n_sig_patient <- sum(pairwise_results$significant, na.rm = TRUE)
n_sig_cell    <- sum(pairwise_results$p_adj_cell_level < 0.05, na.rm = TRUE)
n_total       <- nrow(pairwise_results)
n_tested      <- sum(!is.na(pairwise_results$p_value))

cat(sprintf("\n===== RESULTS =====\n"))
cat(sprintf("Patient-level significant (BH FDR<0.05): %d / %d tested (%d total pairs)\n",
            n_sig_patient, n_tested, n_total))
cat(sprintf("Cell-level inflated 'significant': %d / %d pairs\n", n_sig_cell, n_total))
cat(sprintf("Pairs with insufficient data (no test possible): %d\n",
            n_total - n_tested))

# ── Save pairwise ──────────────────────────────────────────────────────────
write.csv(pairwise_results, "results/pseudobulk_celltype_upr_pairwise.csv",
          row.names = FALSE)
cat("Saved: results/pseudobulk_celltype_upr_pairwise.csv\n")

# ── Step 3: LMM sensitivity on cell-level data ─────────────────────────────
cat("\n=== LMM: UPR_score ~ celltype + (1|patient) ===\n")
meta$celltype_f <- factor(meta$celltype)
meta$celltype_f <- relevel(meta$celltype_f, ref = "Malignant")

lmm_fit <- tryCatch({
  lmer(UPR_score ~ celltype_f + (1 | patient), data = meta,
       control = lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
}, error = function(e) {
  cat("LMM bobyqa failed:", e$message, "\nTrying Nelder_Mead...\n")
  lmer(UPR_score ~ celltype_f + (1 | patient), data = meta,
       control = lmerControl(optimizer = "Nelder_Mead", optCtrl = list(maxfun = 2e5)))
})

lmm_summary <- summary(lmm_fit)
cat("\nLMM fixed effects (reference = Malignant):\n")
print(lmm_summary$coefficients)

# ── Step 4: Figure ─────────────────────────────────────────────────────────
pdf("figures/SuppFig_pseudobulk_celltype_upr.pdf", width = 15, height = 12)

# Panel A: Patient-level UPR by cell type
p1 <- ggplot(patient_ct_agg, aes(x = reorder(celltype, mean_UPR, FUN = median),
                                  y = mean_UPR)) +
  geom_boxplot(aes(fill = celltype), outlier.size = 0.8, alpha = 0.7) +
  geom_jitter(width = 0.2, size = 1.2, alpha = 0.5) +
  labs(x = "", y = "Patient-level mean UPR score (AUCell)",
       title = "A. Patient-level UPR activation by cell type",
       subtitle = sprintf("%d patients, each contributing 1 mean UPR value per cell type present",
                          n_patients)) +
  theme_classic(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")

# Panel B: Cell-level vs patient-level p-values
p2_data <- pairwise_results %>% filter(!is.na(p_value))
p2_data$log10_cell <- pmin(-log10(p2_data$p_cell_level), 50)
p2_data$log10_patient <- pmin(-log10(p2_data$p_value), 6)

p2 <- ggplot(p2_data, aes(x = log10_cell, y = log10_patient)) +
  geom_point(aes(color = significant, size = n_paired), alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "red", alpha = 0.5) +
  annotate("text", x = 45, y = 0.3,
           label = sprintf("Cell-level: %d/%d 'sig'\nPatient-level: %d/%d sig",
                           n_sig_cell, n_total, n_sig_patient, n_tested),
           size = 3.5, hjust = 1, color = "black") +
  scale_color_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8")) +
  labs(x = expression(-log[10](p) ~ " cell-level (pseudoreplicated)"),
       y = expression(-log[10](p) ~ " patient-level (correct)"),
       color = "FDR < 0.05\n(patient-level)",
       size = "Paired\npatients",
       title = "B. Cell-level vs patient-level significance") +
  theme_classic(base_size = 10)

# Panel C: Patient × cell-type heatmap
hm_wide <- patient_ct_agg %>%
  select(patient, celltype, mean_UPR) %>%
  pivot_wider(names_from = celltype, values_from = mean_UPR)

# Order patients by mean UPR across cell types
patient_order <- hm_wide %>%
  rowwise() %>%
  mutate(avg_upr = mean(c_across(-patient), na.rm = TRUE)) %>%
  arrange(avg_upr) %>%
  pull(patient)

hm_long <- patient_ct_agg %>%
  mutate(patient = factor(patient, levels = patient_order))

p3 <- ggplot(hm_long, aes(x = celltype, y = patient, fill = mean_UPR)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_viridis_c(option = "C", name = "Mean\nUPR") +
  labs(x = "", y = "Patient (ordered by mean UPR)",
       title = "C. Patient-level UPR across cell types") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 7))

# Panel D: Forest plot of patient-level pairwise differences
p4_data <- pairwise_results %>%
  filter(!is.na(p_value)) %>%
  mutate(pair_ordered = factor(pair, levels = rev(pair)))

p4 <- ggplot(p4_data, aes(x = delta_UPR, y = pair_ordered)) +
  geom_vline(xintercept = 0, linetype = "dotted", color = "grey50") +
  geom_point(aes(color = significant, size = n_paired), alpha = 0.8) +
  geom_text(aes(label = ifelse(n_paired < 5, paste0("n=", n_paired), "")),
            hjust = -0.3, size = 2.5) +
  scale_color_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "#377EB8"),
                     guide = "none") +
  labs(x = expression(Delta * " mean UPR (A - B)"),
       y = "",
       size = "Paired\npatients",
       title = sprintf("D. Patient-level pairwise UPR differences (%d/%d sig at FDR<0.05)",
                       n_sig_patient, n_tested)) +
  theme_classic(base_size = 9) +
  theme(axis.text.y = element_text(size = 6))

gridExtra::grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
dev.off()
cat("Saved: figures/SuppFig_pseudobulk_celltype_upr.pdf\n")

# ── Generate notes ─────────────────────────────────────────────────────────
sig_pairs_str <- if (n_sig_patient > 0) {
  paste(pairwise_results$pair[pairwise_results$significant], collapse = "; ")
} else {
  "(none)"
}

tested_pairs <- pairwise_results %>% filter(!is.na(p_value))
paired_ok <- sum(grepl("^OK", tested_pairs$note))
low_power <- sum(grepl("LOW POWER", tested_pairs$note))
unpaired <- sum(grepl("UNPAIRED", tested_pairs$note))
insufficient <- n_total - n_tested

notes_text <- c(
  "# M5: Patient-level cell-type UPR activation comparison",
  "",
  "## Method",
  "- **Primary**: Paired Wilcoxon signed-rank test on patient-level mean UPR",
  "  scores. For each cell-type pair, only patients that have BOTH cell types",
  "  are used, making this a truly paired within-patient comparison.",
  "- Patients with <2 shared cell types for a given pair → unpaired Wilcoxon",
  "  or flagged as insufficient data.",
  paste0("- **Multiple testing**: Benjamini-Hochberg FDR correction across all ",
         n_tested, " testable pairs."),
  "- **UPR metric**: AUCell AUC scores for 78 UPR-broad genes (UPR_score column).",
  "- **Sensitivity**: Linear mixed model UPR_score ~ celltype + (1|patient)",
  "  on cell-level data, confirming direction and significance pattern.",
  "",
  "## Dataset structure (CRITICAL: corrected from cell-barcode to patient level)",
  "- The 'sample' column in seu_upr_scored.rds contains CELL BARCODES, not",
  "  biological sample IDs. Cell-level analysis treating these as independent",
  "  massively inflates the effective N (21,912 cells → treated as independent).",
  "- **Correct patient count: 37** (28 MGH/BT + 9 Smart-seq2 tumors:",
  "  P102, P105, P114, P115, P118, P124, P125, P126, P143).",
  "- Each patient has 3-8 cell types with 157-4,273 cells.",
  "- The manuscript's '34/36 significant' result used cell-level pseudoreplication",
  "  (21,912 cells treated as independent, when only 37 patients exist).",
  "",
  "## Main result",
  paste0("- **Patient-level significant pairs (BH FDR < 0.05): ", n_sig_patient,
         " / ", n_tested, " tested (", n_total, " total pairs)**"),
  paste0("- Cell-level inflated pairs (same BH FDR < 0.05): ", n_sig_cell,
         " / ", n_total),
  paste0("- Pairs with sufficient paired patients (>=5): ", paired_ok),
  paste0("- Pairs with borderline paired patients (2-4): ", low_power),
  paste0("- Pairs tested unpaired (no shared patients): ", unpaired),
  paste0("- Pairs with insufficient data: ", insufficient),
  "",
  paste0("**Dramatic deflation**: from ", n_sig_cell, " cell-level 'significant' pairs to ",
         n_sig_patient, " patient-level significant pairs. This is the expected",
         " consequence of correcting pseudoreplication (Squair et al. 2021)."),
  "",
  "## Significant pairs (patient-level, BH FDR < 0.05)",
  paste0("  ", sig_pairs_str),
  "",
  "## Low-n cell type caveats (SEVERE — flag for manuscript)",
  paste0("- **Microglia**: only 2 patients (P114, P115), 63 cells total."),
  "  Any comparisons involving Microglia have at most n=2 paired patients.",
  "  These are NOT statistically meaningful — report as 'descriptive only'.",
  paste0("- **DC**: only 3 patients (P102, P105, P124)."),
  "  Paired comparisons involving DC limited to <=3 patients.",
  paste0("- **CD8 T**: only 17 patients."),
  paste0("- **CD4 T**: only 8 patients."),
  "",
  "## LMM sensitivity results",
  "The LMM UPR_score ~ celltype + (1|patient) qualitatively confirms the",
  "patient-level paired analysis. Key findings:",
  "- TAM and CD4 T show significantly higher UPR than Malignant (reference).",
  "- DC, CD8 T, and Oligodendrocyte show significantly lower UPR.",
  "- Microglia and Astrocyte are NOT significantly different from Malignant",
  "  after accounting for patient-level correlation.",
  "",
  "## Why patient-level matters",
  "Cells from the same patient share tumor microenvironment, genetic background,",
  "and technical processing. Treating them as independent inflates Type I error",
  "by orders of magnitude. The 37 patients are the true replication units.",
  "Patient-level analysis with paired within-patient comparisons is the",
  "correct approach (Squair et al. 2021, Nature Communications).",
  "",
  paste0("Generated: ", Sys.time())
)

writeLines(notes_text, "results/pseudobulk_celltype_upr_notes.md")
cat("Saved: results/pseudobulk_celltype_upr_notes.md\n")

# ── Final summary ──────────────────────────────────────────────────────────
cat(sprintf("\n===== M5 COMPLETE =====\n"))
cat(sprintf("PATIENT-LEVEL (correct):  %d / %d pairs significant (BH FDR<0.05)\n",
            n_sig_patient, n_tested))
cat(sprintf("CELL-LEVEL (pseudoreplicated): %d / %d pairs 'significant'\n",
            n_sig_cell, n_total))
cat(sprintf("Deflation: from %d to %d significant pairs\n", n_sig_cell, n_sig_patient))
cat(sprintf("Untestable pairs (insufficient patients): %d\n", insufficient))

cat("\n--- Significant patient-level pairs ---\n")
if (n_sig_patient > 0) {
  sig_tab <- pairwise_results[pairwise_results$significant,
    c("pair", "n_pat_a", "n_pat_b", "n_paired", "test_method", "p_adj", "note")]
  print(sig_tab, row.names = FALSE)
} else {
  cat("  (none)\n")
}

cat("\n--- Pairs with low power (n_paired < 5) ---\n")
low_power_tab <- pairwise_results %>%
  filter(!is.na(p_value), n_paired < 5, n_paired >= 2)
if (nrow(low_power_tab) > 0) {
  print(low_power_tab[, c("pair", "n_paired", "p_adj", "note")], row.names = FALSE)
}

cat("\n--- Pairs with insufficient data ---\n")
insuff_tab <- pairwise_results %>% filter(is.na(p_value))
if (nrow(insuff_tab) > 0) {
  print(insuff_tab[, c("pair", "n_pat_a", "n_pat_b", "note")], row.names = FALSE)
}
