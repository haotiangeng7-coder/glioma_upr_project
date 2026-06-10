###############################################################################
# M2: Head-to-head comparison of UIRS vs published glioma ER-stress/UPR
#     prognostic signatures on CGGA-batch1 (blind test set)
# Outputs: results/published_signature_comparison.csv,
#          results/published_signature_comparison_notes.md,
#          figures/SuppFig_published_signature_comparison.pdf
###############################################################################

suppressPackageStartupMessages({
  library(survival)
})

set.seed(42)
PROJECT_DIR <- getwd()
DATA_PROC   <- file.path(PROJECT_DIR, "data", "processed")
FIG_DIR     <- file.path(PROJECT_DIR, "figures")
RES_DIR     <- file.path(PROJECT_DIR, "results")

# =============================================================================
# 1. Load data
# =============================================================================
message("=== Loading data ===")
load(file.path(DATA_PROC, "coxph_interpretable_model.RData"))
load(file.path(DATA_PROC, "cgga_data.RData"))
load(file.path(DATA_PROC, "tcga_glioma_expression.RData"))
load(file.path(DATA_PROC, "upr_gene_sets.RData"))

uirs_betas  <- coxph_results$betas           # 15 genes, named numeric
train_df    <- coxph_results$train_df         # TCGA training
cgga1_val   <- coxph_results$cgga1_val_df     # CGGA-batch1 (has precomputed risk_score)

tcga_tpm  <- expr_tpm_symbol                   # 59427 x 925
cgga1_raw <- cgga_b1$expr                      # 24326 x 325
cgga1_clin <- cgga_b1$clinical                 # 313 rows

# =============================================================================
# 2. Bootstrap C-index (reverse=TRUE: higher risk_score = worse survival)
# =============================================================================
boot_cindex <- function(time, status, risk_score, B = 1000) {
  keep <- !is.na(time) & !is.na(status) & !is.na(risk_score) & is.finite(risk_score)
  time <- time[keep]; status <- status[keep]; risk_score <- risk_score[keep]
  n <- length(time)
  if (n < 20) return(list(cindex = NA, ci_low = NA, ci_high = NA))
  c_obs <- as.numeric(concordance(Surv(time, status) ~ risk_score, reverse = TRUE)$concordance)
  boot_c <- replicate(B, {
    idx <- sample(n, replace = TRUE)
    tryCatch(
      as.numeric(concordance(Surv(time[idx], status[idx]) ~ risk_score[idx], reverse = TRUE)$concordance),
      error = function(e) NA_real_)
  })
  boot_c <- boot_c[!is.na(boot_c)]
  if (length(boot_c) < 10) return(list(cindex = c_obs, ci_low = NA, ci_high = NA))
  ci <- quantile(boot_c, c(0.025, 0.975), na.rm = TRUE)
  list(cindex = c_obs, ci_low = ci[1], ci_high = ci[2])
}

# =============================================================================
# 3. UIRS reference C-index on CGGA-batch1 (precomputed risk_score)
# =============================================================================
message("=== UIRS reference C-index ===")
uirs_c <- boot_cindex(cgga1_val$time, cgga1_val$status, cgga1_val$risk_score, B = 1000)
message(sprintf("UIRS: %.4f (%.4f-%.4f)", uirs_c$cindex, uirs_c$ci_low, uirs_c$ci_high))

# =============================================================================
# 4. Match CGGA-batch1 val to expression via clinical sample IDs
# =============================================================================
message("=== Matching CGGA-batch1 to expression ===")

match_val_to_expr <- function(val_df, clinical, expr_mat) {
  # val_df rownames are indices; match to clinical by (time, status, IDH_status)
  v_key <- paste(val_df$time, val_df$status, val_df$IDH_status)
  c_key <- paste(clinical$OS.time, clinical$OS_status, clinical$IDH_status)

  m <- match(v_key, c_key)
  # greedy assignment for duplicate keys
  dup_keys <- unique(v_key[duplicated(v_key)])
  for (k in dup_keys) {
    v_pos <- which(v_key == k)
    c_pos <- which(c_key == k)
    for (i in seq_along(v_pos)) {
      if (i <= length(c_pos)) m[v_pos[i]] <- c_pos[i]
    }
  }

  clin_idx <- m[!is.na(m)]
  val_idx  <- which(!is.na(m))
  sample_ids <- clinical$Sample_ID[clin_idx]

  # Extract expression for matched samples
  matched_cols <- intersect(sample_ids, colnames(expr_mat))
  if (length(matched_cols) == 0) stop("No expression columns matched!")

  expr_sub <- expr_mat[, matched_cols, drop = FALSE]
  # log2 transform
  expr_log <- log2(expr_sub + 1)

  # Align val_df rows
  list(
    expr_log    = expr_log,
    sample_ids  = matched_cols,
    time        = val_df$time[val_idx][match(matched_cols, sample_ids)],
    status      = val_df$status[val_idx][match(matched_cols, sample_ids)],
    val_idx     = val_idx[match(matched_cols, sample_ids)],
    risk_score  = val_df$risk_score[val_idx][match(matched_cols, sample_ids)]
  )
}

cgga1_data <- match_val_to_expr(cgga1_val, cgga1_clin, cgga1_raw)
message(sprintf("CGGA-b1 matched expression: %d genes x %d samples",
                nrow(cgga1_data$expr_log), ncol(cgga1_data$expr_log)))

# =============================================================================
# 5. Prepare TCGA training expression (full, log2)
# =============================================================================
message("=== Preparing TCGA training expression ===")
train_samples <- intersect(rownames(train_df), colnames(tcga_tpm))
tcga_log <- log2(tcga_tpm[, train_samples, drop = FALSE] + 1)
tcga_time   <- train_df[train_samples, "time"]
tcga_status <- train_df[train_samples, "status"]
message(sprintf("TCGA training: %d genes x %d samples", nrow(tcga_log), ncol(tcga_log)))

# =============================================================================
# 6. Comparator 1: Zhang et al. 2021 (PMID 33611848) — 16-gene locked
# =============================================================================
message("=== Comparator 1: Zhang et al. 2021 ===")

zhang_genes <- c("CYP2E1","SLN","BRCA1","CISD2","LRRK2","BMP2","MYH7","HSPB1",
                 "DNM1L","SHISA5","RNF185","RCN1","SPP1","RPN2","PDIA3","ATP2A2")
zhang_coefs <- c(CYP2E1 = -0.1318, SLN = -0.1122, BRCA1 = 0.4543, CISD2 = 0.3481,
                 LRRK2 = 0.1811, BMP2 = -0.1236, MYH7 = -0.1183, HSPB1 = 0.4296,
                 DNM1L = 0.3839, SHISA5 = 0.5514, RNF185 = -0.5315, RCN1 = 0.2889,
                 SPP1 = 0.1260, RPN2 = 0.5187, PDIA3 = -0.7779, ATP2A2 = 0.2649)

z_cgga_genes <- intersect(zhang_genes, rownames(cgga1_data$expr_log))
z_score <- rep(0, ncol(cgga1_data$expr_log))
for (g in intersect(z_cgga_genes, names(zhang_coefs))) {
  z_score <- z_score + zhang_coefs[g] * cgga1_data$expr_log[g, ]
}

z_valid <- !is.na(z_score) & is.finite(z_score)
z_c <- if (sum(z_valid) >= 20) boot_cindex(
  cgga1_data$time[z_valid], cgga1_data$status[z_valid], z_score[z_valid], B = 1000
) else list(cindex = NA, ci_low = NA, ci_high = NA)
z_ovl <- length(intersect(names(uirs_betas), zhang_genes))

message(sprintf("Zhang: C=%.4f (%.4f-%.4f), genes=%d/%d in CGGA, ovl=%d",
                z_c$cindex, z_c$ci_low, z_c$ci_high, length(z_cgga_genes), 16, z_ovl))

# =============================================================================
# 7. Comparator 2: Shao et al. 2025 (PMID 40718795) — 6-gene ESURATAG, refit
# =============================================================================
message("=== Comparator 2: Shao et al. 2025 ===")

esa_genes <- c("DERL2","RPN2","SEC13","SEC61A1","SEC61B","STT3A")
esa_tcga  <- esa_genes[esa_genes %in% rownames(tcga_log)]
esa_cgga  <- esa_genes[esa_genes %in% rownames(cgga1_data$expr_log)]

# Refit Cox on TCGA
esa_train <- as.data.frame(t(tcga_log[esa_tcga, , drop = FALSE]))
colnames(esa_train) <- esa_tcga
esa_train$time <- tcga_time; esa_train$status <- tcga_status
esa_fit <- coxph(as.formula(paste0("Surv(time, status) ~ ", paste(esa_tcga, collapse = " + "))),
                  data = esa_train)
esa_betas <- coef(esa_fit)

# Score CGGA
esa_score <- rep(0, ncol(cgga1_data$expr_log))
for (g in intersect(esa_cgga, names(esa_betas))) {
  esa_score <- esa_score + esa_betas[g] * cgga1_data$expr_log[g, ]
}

esa_valid <- !is.na(esa_score) & is.finite(esa_score)
esa_c <- if (sum(esa_valid) >= 20) boot_cindex(
  cgga1_data$time[esa_valid], cgga1_data$status[esa_valid], esa_score[esa_valid], B = 1000
) else list(cindex = NA, ci_low = NA, ci_high = NA)
esa_ovl <- length(intersect(names(uirs_betas), esa_genes))

message(sprintf("ESURATAG: C=%.4f (%.4f-%.4f), genes=%d/%d in CGGA, ovl=%d",
                esa_c$cindex, esa_c$ci_low, esa_c$ci_high, length(esa_cgga), 6, esa_ovl))

# =============================================================================
# 8. Comparator 3: HALLMARK_UNFOLDED_PROTEIN_RESPONSE (mean-z baseline)
# =============================================================================
message("=== Comparator 3: HALLMARK_UPR mean-z ===")

hupr_all <- upr_gene_list$Hallmark_UPR
hupr_use <- Reduce(intersect, list(hupr_all, rownames(tcga_log), rownames(cgga1_data$expr_log)))

# z-score params from TCGA
hupr_tcga <- tcga_log[hupr_use, , drop = FALSE]
hupr_mean <- rowMeans(hupr_tcga, na.rm = TRUE)
hupr_sd   <- apply(hupr_tcga, 1, sd, na.rm = TRUE); hupr_sd[hupr_sd == 0] <- 1

# Mean-z on CGGA
hupr_cgga <- cgga1_data$expr_log[hupr_use, , drop = FALSE]
hupr_z <- sweep(sweep(hupr_cgga, 1, hupr_mean, "-"), 1, hupr_sd, "/")
hupr_score <- colMeans(hupr_z, na.rm = TRUE)

h_valid <- !is.na(hupr_score) & is.finite(hupr_score)
h_c <- if (sum(h_valid) >= 20) boot_cindex(
  cgga1_data$time[h_valid], cgga1_data$status[h_valid], hupr_score[h_valid], B = 1000
) else list(cindex = NA, ci_low = NA, ci_high = NA)
h_ovl <- length(intersect(names(uirs_betas), hupr_all))

message(sprintf("Hallmark: C=%.4f (%.4f-%.4f), genes_used=%d, ovl=%d",
                h_c$cindex, h_c$ci_low, h_c$ci_high, length(hupr_use), h_ovl))

# =============================================================================
# 9. Paired delta C-index (UIRS - comparator)
# =============================================================================
message("=== Delta C-index ===")

delta_boot <- function(time, status, r1, r2, B = 1000) {
  keep <- !is.na(r1) & !is.na(r2) & is.finite(r1) & is.finite(r2) &
          !is.na(time) & !is.na(status)
  time <- time[keep]; status <- status[keep]; r1 <- r1[keep]; r2 <- r2[keep]
  n <- length(time)
  if (n < 20) return(list(delta = NA, ci_low = NA, ci_high = NA))
  d_obs <- as.numeric(concordance(Surv(time, status) ~ r1, reverse = TRUE)$concordance -
                      concordance(Surv(time, status) ~ r2, reverse = TRUE)$concordance)
  boot_d <- replicate(B, {
    idx <- sample(n, replace = TRUE)
    c1 <- tryCatch(as.numeric(concordance(Surv(time[idx],status[idx]) ~ r1[idx], reverse = TRUE)$concordance), error = function(e) NA)
    c2 <- tryCatch(as.numeric(concordance(Surv(time[idx],status[idx]) ~ r2[idx], reverse = TRUE)$concordance), error = function(e) NA)
    if (is.na(c1) || is.na(c2)) NA else c1 - c2
  })
  boot_d <- boot_d[!is.na(boot_d)]
  if (length(boot_d) < 10) return(list(delta = d_obs, ci_low = NA, ci_high = NA))
  ci <- quantile(boot_d, c(0.025, 0.975))
  list(delta = d_obs, ci_low = ci[1], ci_high = ci[2])
}

# UIRS risk_scores aligned with CGGA expression samples
uirs_r <- cgga1_data$risk_score
names(uirs_r) <- cgga1_data$sample_ids

compute_delta <- function(c_score, c_valid) {
  if (sum(c_valid) < 20) return(list(delta = NA, ci_low = NA, ci_high = NA))
  delta_boot(cgga1_data$time[c_valid], cgga1_data$status[c_valid],
             uirs_r[c_valid], c_score[c_valid], B = 1000)
}

d_z <- compute_delta(z_score, z_valid)
d_e <- compute_delta(esa_score, esa_valid)
d_h <- compute_delta(hupr_score, h_valid)

message(sprintf("Delta Zhang:   %+.4f (%.4f-%.4f)", d_z$delta, d_z$ci_low, d_z$ci_high))
message(sprintf("Delta ESURATAG: %+.4f (%.4f-%.4f)", d_e$delta, d_e$ci_low, d_e$ci_high))
message(sprintf("Delta Hallmark: %+.4f (%.4f-%.4f)", d_h$delta, d_h$ci_low, d_h$ci_high))

# =============================================================================
# 10. Results CSV
# =============================================================================
message("=== Writing CSV ===")
comp_df <- data.frame(
  name               = c("Zhang_2021_16gene", "Shao_2025_ESURATAG", "HALLMARK_UPR_meanz"),
  source_pmid         = c("33611848", "40718795", "MSigDB_HALLMARK_UPR"),
  n_genes_total       = c(16L, 6L, length(hupr_all)),
  genes_found_cgga    = c(length(z_cgga_genes), length(esa_cgga), length(hupr_use)),
  genes_overlap_uirs  = c(z_ovl, esa_ovl, h_ovl),
  n_samples           = c(sum(z_valid), sum(esa_valid), sum(h_valid)),
  cindex             = c(z_c$cindex, esa_c$cindex, h_c$cindex),
  ci_low             = c(z_c$ci_low, esa_c$ci_low, h_c$ci_low),
  ci_high            = c(z_c$ci_high, esa_c$ci_high, h_c$ci_high),
  uirs_cindex        = rep(uirs_c$cindex, 3),
  delta              = c(d_z$delta, d_e$delta, d_h$delta),
  delta_ci_low       = c(d_z$ci_low, d_e$ci_low, d_h$ci_low),
  delta_ci_high      = c(d_z$ci_high, d_e$ci_high, d_h$ci_high),
  coeff_source       = c("published_locked", "refit_on_TCGA", "mean_z_score"),
  stringsAsFactors   = FALSE
)
write.csv(comp_df, file.path(RES_DIR, "published_signature_comparison.csv"), row.names = FALSE)

# =============================================================================
# 11. Notes
# =============================================================================
notes <- c(
  "# Published Signature Comparison -- Notes",
  "",
  "## UIRS (reference)",
  sprintf("- 15-gene LASSO->CoxPH signature with strata(IDH)"),
  sprintf("- C-index on CGGA-batch1 (blind test): %.4f (%.4f-%.4f)", uirs_c$cindex, uirs_c$ci_low, uirs_c$ci_high),
  sprintf("- N = %d, events = %d", nrow(cgga1_val), sum(cgga1_val$status == 1)),
  "",
  "## Comparator 1: Zhang et al. 2021 (PMID 33611848)",
  "- J Cell Mol Med, 16-gene ER-stress signature, published LASSO coefficients (locked)",
  sprintf("- C-index: %.4f (%.4f-%.4f)", z_c$cindex, z_c$ci_low, z_c$ci_high),
  sprintf("- Delta (UIRS - Zhang): %+.4f (%.4f-%.4f)", d_z$delta, d_z$ci_low, d_z$ci_high),
  sprintf("- Genes (%d found/%d total): %s", length(z_cgga_genes), 16, paste(zhang_genes, collapse=", ")),
  sprintf("- Overlap with UIRS: %d genes (%s)", z_ovl,
          paste(intersect(names(uirs_betas), zhang_genes), collapse=", ")),
  "",
  "## Comparator 2: Shao et al. 2025 (PMID 40718795)",
  "- Front Mol Biosci, 6-gene ESURATAG, Cox refit on TCGA (no published coefs)",
  sprintf("- C-index: %.4f (%.4f-%.4f)", esa_c$cindex, esa_c$ci_low, esa_c$ci_high),
  sprintf("- Delta (UIRS - ESURATAG): %+.4f (%.4f-%.4f)", d_e$delta, d_e$ci_low, d_e$ci_high),
  sprintf("- Genes (%d found/%d total): %s", length(esa_cgga), 6, paste(esa_genes, collapse=", ")),
  sprintf("- Overlap with UIRS: %d genes (%s)", esa_ovl,
          paste(intersect(names(uirs_betas), esa_genes), collapse=", ")),
  "",
  "## Comparator 3: HALLMARK_UPR (MSigDB geneset baseline)",
  sprintf("- %d-gene set, mean-z score (%d genes scored in CGGA)", length(hupr_all), length(hupr_use)),
  sprintf("- C-index: %.4f (%.4f-%.4f)", h_c$cindex, h_c$ci_low, h_c$ci_high),
  sprintf("- Delta (UIRS - Hallmark): %+.4f (%.4f-%.4f)", d_h$delta, d_h$ci_low, d_h$ci_high),
  "",
  "## Unreproducible signatures",
  "- Li et al. 2022 (PMID 35614390): 7-gene signature, only 4/7 confirmed from accessible text",
  "- Fan et al. 2022 (J Oncol): 8-gene list complete but no published coefficients",
  "",
  "## Methods",
  "- All C-indices reported with reverse=TRUE (higher risk/lp = worse survival)",
  "- Bootstrap 95% CI: B=1000 resamples",
  "- CGGA-batch1 matching: (OS.time, OS_status, IDH_status) composite key",
  "- TCGA expression: log2(TPM+1); CGGA expression: log2(FPKM+1)"
)
writeLines(notes, file.path(RES_DIR, "published_signature_comparison_notes.md"))

# =============================================================================
# 12. Figure (base R barplot)
# =============================================================================
message("=== Figure ===")
pdf(file.path(FIG_DIR, "SuppFig_published_signature_comparison.pdf"), width = 9, height = 5)
par(mar = c(4, 10, 3, 2))
labels <- c("UIRS (15-gene)", "Zhang 2021 (16-gene)", "Shao 2025 (6-gene)", "Hallmark UPR (mean-z)")
ypos <- 4:1
cvals <- c(uirs_c$cindex, z_c$cindex, esa_c$cindex, h_c$cindex)
clows <- c(uirs_c$ci_low, z_c$ci_low, esa_c$ci_low, h_c$ci_low)
chighs <- c(uirs_c$ci_high, z_c$ci_high, esa_c$ci_high, h_c$ci_high)
cols <- c("#E64B35", "#4DBBD5", "#00A087", "grey50")
xr <- range(c(clows, chighs), na.rm = TRUE)
plot(NA, NA, xlim = c(xr[1] - 0.03, xr[2] + 0.03), ylim = c(0.5, 4.5),
     xlab = "Harrell C-index (CGGA-batch1, bootstrap 95% CI)", ylab = "",
     yaxt = "n", main = "Prognostic Performance: UIRS vs Published Signatures")
axis(2, at = ypos, labels = labels, las = 1, cex.axis = 0.9)
for (i in 1:4) {
  if (!is.na(clows[i]) && !is.na(chighs[i]) && clows[i] < chighs[i])
    segments(clows[i], ypos[i], chighs[i], ypos[i], lwd = 2, col = cols[i])
  if (!is.na(cvals[i])) {
    points(cvals[i], ypos[i], pch = 19, cex = 2, col = cols[i])
    text(cvals[i], ypos[i] + 0.3, sprintf("%.3f", cvals[i]), cex = 0.8)
  }
}
legend("bottomright",
       legend = c("UIRS (this study)", "Published (locked coef)", "Published (refit)", "Gene-set baseline"),
       col = cols, pch = 19, bty = "n", cex = 0.8)
dev.off()

message("=== M2 completed ===")
cat(sprintf("UIRS:         %.4f (%.4f-%.4f)\n", uirs_c$cindex, uirs_c$ci_low, uirs_c$ci_high))
cat(sprintf("Zhang 2021:   %.4f (%.4f-%.4f)\n", z_c$cindex, z_c$ci_low, z_c$ci_high))
cat(sprintf("Shao 2025:    %.4f (%.4f-%.4f)\n", esa_c$cindex, esa_c$ci_low, esa_c$ci_high))
cat(sprintf("Hallmark UPR: %.4f (%.4f-%.4f)\n", h_c$cindex, h_c$ci_low, h_c$ci_high))
