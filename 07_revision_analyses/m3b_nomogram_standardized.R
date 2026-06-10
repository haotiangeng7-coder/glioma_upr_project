###############################################################################
# M3b: Nomogram with FROZEN z-standardized risk_score
# Fix: risk_score scale mismatch ~50x (TCGA~97 vs CGGA~2) causes LP reversal
# Method: z-standardize to TCGA-frozen mean/sd BEFORE model fitting
# Outputs: results/nomogram_standardized_*.csv, figures/Fig6_*_standardized.pdf,
#          results/nomogram_standardized_notes.md
###############################################################################

suppressPackageStartupMessages({
  library(survival)
  library(rms)
  library(timeROC)
})

set.seed(42)
PROJECT_DIR <- getwd()
DATA_PROC   <- file.path(PROJECT_DIR, "data", "processed")
FIG_DIR     <- file.path(PROJECT_DIR, "figures")
RES_DIR     <- file.path(PROJECT_DIR, "results")

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RES_DIR, showWarnings = FALSE, recursive = TRUE)

# =============================================================================
# 1. Load all ground-truth objects
# =============================================================================
message("=== Loading data ===")
load(file.path(DATA_PROC, "coxph_interpretable_model.RData"))      # -> coxph_results
load(file.path(DATA_PROC, "independent_prognosis_results.RData"))   # -> analysis_df
load(file.path(DATA_PROC, "cgga_data.RData"))                       # -> cgga_b1, cgga_b2

betas         <- coxph_results$betas
gene_names    <- names(betas)
train_df      <- coxph_results$train_df
cgga1_val     <- coxph_results$cgga1_val_df
cgga2_val     <- coxph_results$cgga2_val_df

# =============================================================================
# 2. Compute TCGA raw UIRS from locked betas, freeze standardization
# =============================================================================
message("=== Step 1: Frozen z-standardization of risk_score ===")

# Raw UIRS = sum(beta_g * expr_g) on TCGA training genes
expr_tcga     <- as.matrix(train_df[, gene_names])
raw_risk_tcga <- as.numeric(expr_tcga %*% betas)

mean_TCGA <- mean(raw_risk_tcga)
sd_TCGA   <- sd(raw_risk_tcga)

message(sprintf("FROZEN TCGA standardization: mean=%.6f, sd=%.6f", mean_TCGA, sd_TCGA))
message(sprintf("TCGA raw UIRS range: [%.3f, %.3f]", min(raw_risk_tcga), max(raw_risk_tcga)))

# z-transform TCGA
z_risk_tcga <- (raw_risk_tcga - mean_TCGA) / sd_TCGA

message(sprintf("TCGA z-risk: mean=%.3f, sd=%.3f, range=[%.3f, %.3f]",
                mean(z_risk_tcga), sd(z_risk_tcga),
                min(z_risk_tcga), max(z_risk_tcga)))

# =============================================================================
# 3. Merge TCGA z_risk with clinical (Age, Grade, IDH, OS) from analysis_df
# =============================================================================
message("=== Step 2: Merging TCGA z_risk with clinical ===")

# Build lookup: train_df sample_id -> z_risk
train_df$sample_id <- substr(rownames(train_df), 1, 16)
train_df$z_risk    <- z_risk_tcga

# Deduplicate: keep first occurrence per sample_id
train_uniq <- train_df[!duplicated(train_df$sample_id), ]

# Merge with analysis_df
analysis_df$sample_id <- substr(analysis_df$barcode, 1, 16)
# Deduplicate analysis_df too (keep first)
analysis_uniq <- analysis_df[!duplicated(analysis_df$sample_id), ]

merged <- merge(train_uniq[, c("sample_id", "z_risk")],
                analysis_uniq[, c("sample_id", "OS_months", "OS", "Age", "Grade", "IDH")],
                by = "sample_id")

message(sprintf("Merged TCGA: %d samples", nrow(merged)))
message(sprintf("  Events: %d (%.1f%%)", sum(merged$OS), 100*mean(merged$OS)))
message(sprintf("  Age: mean=%.1f, range=[%.0f, %.0f]", mean(merged$Age, na.rm=TRUE),
                min(merged$Age, na.rm=TRUE), max(merged$Age, na.rm=TRUE)))
message("  Grade table:")
print(table(merged$Grade))
message("  IDH table:")
print(table(merged$IDH))

# Filter to complete cases, OS_months > 0
nomo_df <- merged[, c("OS_months", "OS", "z_risk", "Age", "Grade", "IDH")]
nomo_df <- nomo_df[complete.cases(nomo_df), ]
nomo_df <- nomo_df[nomo_df$OS_months > 0, ]
nomo_df$Grade <- factor(nomo_df$Grade, levels = c("G2", "G3", "G4"))
nomo_df$IDH   <- factor(nomo_df$IDH, levels = c("Mutant", "WT"))

message(sprintf("TCGA nomogram complete cases: %d samples, %d events",
                nrow(nomo_df), sum(nomo_df$OS)))

# =============================================================================
# 4. Build rms::cph model on TCGA
# =============================================================================
message("=== Step 3: Building TCGA rms Cox model ===")

dd <- datadist(nomo_df)
options(datadist = "dd")

fit_cph <- cph(Surv(OS_months, OS) ~ z_risk + Age + Grade + IDH,
               data = nomo_df, x = TRUE, y = TRUE, surv = TRUE, time.inc = 36)

# Internal C-index with bootstrap 95% CI (B=1000, Harrell's)
B <- 1000
boot_cidx <- replicate(B, {
  idx <- sample(nrow(nomo_df), replace = TRUE)
  tryCatch({
    bfit <- cph(Surv(OS_months, OS) ~ z_risk + Age + Grade + IDH,
                data = nomo_df[idx, ], surv = TRUE)
    (bfit$stats["Dxy"] + 1) / 2
  }, error = function(e) NA_real_)
})
boot_cidx <- boot_cidx[!is.na(boot_cidx)]
c_int <- (fit_cph$stats["Dxy"] + 1) / 2
ci_int <- quantile(boot_cidx, c(0.025, 0.975), na.rm = TRUE)

message(sprintf("TCGA internal C-index: %.4f (95%% CI: %.4f-%.4f)",
                c_int, ci_int[1], ci_int[2]))

# PH diagnostics
ph_test <- cox.zph(fit_cph)
message(sprintf("Global PH p-value: %.4f", ph_test$table["GLOBAL", "p"]))
message("Per-term PH:")
print(ph_test$table[, c("chisq","df","p")])

# Model coefficients
message("Model coefficients:")
for (nm in names(coef(fit_cph))) {
  message(sprintf("  %s: %.4f", nm, coef(fit_cph)[nm]))
}

# =============================================================================
# 5. Internal bootstrap calibration (B=200) at 1/3/5yr
# =============================================================================
message("=== Step 4: Internal calibration (Bootstrap B=200) ===")

cal_timepoints <- c(12, 36, 60)
cal_labels     <- c("1-Year", "3-Year", "5-Year")
internal_cal <- data.frame(stringsAsFactors = FALSE)

for (ti in seq_along(cal_timepoints)) {
  tp <- cal_timepoints[ti]
  lb <- cal_labels[ti]

  # Re-fit with matching time.inc
  fit_cal <- cph(Surv(OS_months, OS) ~ z_risk + Age + Grade + IDH,
                 data = nomo_df, x = TRUE, y = TRUE, surv = TRUE, time.inc = tp)

  m_per_group <- max(20, min(80, floor(nrow(nomo_df) / 5)))
  cal_obj <- tryCatch(
    calibrate(fit_cal, cmethod = "KM", method = "boot", u = tp,
              m = m_per_group, B = 200),
    error = function(e) NULL
  )

  if (!is.null(cal_obj)) {
    # Extract calibration stats
    cal_data <- data.frame(
      pred   = cal_obj[, "mean.predicted"],
      obs    = cal_obj[, "KM"],
      n      = cal_obj[, "n"],
      stringsAsFactors = FALSE
    )
    # Calibration slope via internal bootstrap Cox
    slp <- tryCatch(coef(coxph(Surv(OS_months, OS) ~ predict(fit_cph, type = "lp"),
                               data = nomo_df)), error = function(e) NA_real_)
    internal_cal <- rbind(internal_cal, data.frame(
      timepoint_months = tp,
      timepoint_label  = lb,
      cal_slope        = slp,
      n                = nrow(nomo_df),
      n_events         = sum(nomo_df$OS),
      stringsAsFactors  = FALSE
    ))
  }

  # Calibration plot
  pdf(file.path(FIG_DIR, sprintf("Fig6_calibration_internal_%dmo_standardized.pdf", tp)),
      width = 7, height = 7)
  if (!is.null(cal_obj)) {
    plot(cal_obj,
         xlab = "Nomogram-Predicted Survival Probability",
         ylab = "Observed Survival (KM)",
         main = sprintf("Internal Calibration: %s Survival (Bootstrap B=200)", lb),
         xlim = c(0, 1), ylim = c(0, 1), subtitles = TRUE)
    abline(0, 1, lty = 2, col = "grey50")
  } else {
    plot.new()
    text(0.5, 0.5, sprintf("Calibration failed for %s", lb), cex = 0.8)
  }
  dev.off()
  message(sprintf("  %s calibration saved", lb))
}

# Combined 3-panel internal calibration
pdf(file.path(FIG_DIR, "Fig6_calibration_internal_standardized.pdf"), width = 18, height = 6)
par(mfrow = c(1, 3), mar = c(5, 5, 4, 2))
for (ti in seq_along(cal_timepoints)) {
  tp <- cal_timepoints[ti]
  lb <- cal_labels[ti]
  fit_cal <- cph(Surv(OS_months, OS) ~ z_risk + Age + Grade + IDH,
                 data = nomo_df, x = TRUE, y = TRUE, surv = TRUE, time.inc = tp)
  m_per_group <- max(20, min(80, floor(nrow(nomo_df) / 5)))
  cal_obj <- calibrate(fit_cal, cmethod = "KM", method = "boot", u = tp,
                        m = m_per_group, B = 200)
  plot(cal_obj,
       xlab = "Predicted Probability", ylab = "Observed (KM)",
       main = sprintf("%s Survival", lb),
       xlim = c(0, 1), ylim = c(0, 1), subtitles = FALSE)
  abline(0, 1, lty = 2, col = "grey50")
}
dev.off()
message("  Combined internal calibration figure saved")

# =============================================================================
# 6. Nomogram figure
# =============================================================================
message("=== Step 5: Nomogram figure ===")

surv_fn <- Survival(fit_cph)
surv_1yr <- function(x) surv_fn(12, x)
surv_3yr <- function(x) surv_fn(36, x)
surv_5yr <- function(x) surv_fn(60, x)

pdf(file.path(FIG_DIR, "Fig6_nomogram_standardized.pdf"), width = 12, height = 8)
nom <- nomogram(fit_cph,
                fun      = list(surv_1yr, surv_3yr, surv_5yr),
                funlabel = c("1-Year Survival", "3-Year Survival", "5-Year Survival"),
                maxscale = 100,
                fun.at   = c(0.95, 0.9, 0.85, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1),
                lp       = TRUE)
plot(nom, xfrac = 0.35,
     cex.var = 0.9, cex.axis = 0.7,
     col.grid = gray(c(0.8, 0.95)))
title("Cox Regression Nomogram for Overall Survival (z-standardized UIRS)", cex.main = 1.2)
dev.off()
message("  Nomogram figure saved")

# =============================================================================
# 7. Prepare CGGA nomogram frames with FROZEN z-standardization
# =============================================================================
message("=== Step 6: Preparing CGGA nomogram frames ===")

prepare_cgga <- function(val_df, raw_clinical, batch_name, mean_tcga, sd_tcga) {
  # Step A: Compute raw UIRS from locked betas on CGGA expression
  expr_cgga <- as.matrix(val_df[, gene_names])
  raw_risk_cgga <- as.numeric(expr_cgga %*% betas)

  # Step B: z-transform with TCGA-frozen mean/sd
  z_risk_cgga <- (raw_risk_cgga - mean_tcga) / sd_tcga

  # Step C: Match val_df to clinical by (OS.time, OS_status, IDH_status)
  # CGGA val_df time is in days (same as clinical OS.time)
  v_key <- paste(round(val_df$time, 0), val_df$status, val_df$IDH_status)
  c_key <- paste(round(raw_clinical$OS.time, 0), raw_clinical$OS_status, raw_clinical$IDH_status)

  m <- match(v_key, c_key)
  # Handle duplicate keys: assign sequentially
  dup_keys <- unique(v_key[duplicated(v_key)])
  for (k in dup_keys) {
    v_pos <- which(v_key == k)
    c_pos <- which(c_key == k)
    for (i in seq_along(v_pos)) {
      if (i <= length(c_pos)) m[v_pos[i]] <- c_pos[i]
    }
  }

  n_matched <- sum(!is.na(m))
  message(sprintf("  %s: %d/%d matched by (time,status,IDH)", batch_name, n_matched, nrow(val_df)))

  keep_v <- which(!is.na(m))
  c_idx  <- m[!is.na(m)]

  # Convert OS.time from days to months
  os_months <- val_df$time[keep_v] / 30.4375
  os_status <- val_df$status[keep_v]
  z_risk    <- z_risk_cgga[keep_v]
  age       <- raw_clinical$Age[c_idx]

  grade <- factor(
    dplyr::case_when(
      raw_clinical$Grade[c_idx] == "WHO II"  ~ "G2",
      raw_clinical$Grade[c_idx] == "WHO III" ~ "G3",
      raw_clinical$Grade[c_idx] == "WHO IV"  ~ "G4",
      TRUE ~ NA_character_
    ),
    levels = c("G2", "G3", "G4")
  )

  idh <- factor(
    dplyr::case_when(
      raw_clinical$IDH_status[c_idx] == "Mutant"   ~ "Mutant",
      raw_clinical$IDH_status[c_idx] == "Wildtype" ~ "WT",
      TRUE ~ NA_character_
    ),
    levels = c("Mutant", "WT")
  )

  df <- data.frame(OS_months = os_months, OS = os_status, z_risk = z_risk,
                   Age = age, Grade = grade, IDH = idh, stringsAsFactors = FALSE)
  df <- df[complete.cases(df) & df$OS_months > 0, ]

  message(sprintf("  %s: %d complete cases, %d events", batch_name, nrow(df), sum(df$OS)))
  message(sprintf("  %s z_risk: mean=%.3f, sd=%.3f, range=[%.3f, %.3f]",
                  batch_name, mean(df$z_risk), sd(df$z_risk),
                  min(df$z_risk), max(df$z_risk)))
  message(sprintf("  %s Age: mean=%.1f, Grade=%s, IDH=%s",
                  batch_name, mean(df$Age, na.rm=TRUE),
                  paste(table(df$Grade), collapse=","),
                  paste(table(df$IDH), collapse=",")))
  df
}

cgga_b1_nomo <- prepare_cgga(cgga1_val, cgga_b1$clinical, "CGGA-batch1", mean_TCGA, sd_TCGA)
cgga_b2_nomo <- prepare_cgga(cgga2_val, cgga_b2$clinical, "CGGA-batch2", mean_TCGA, sd_TCGA)

# =============================================================================
# 8. External C-index with bootstrap 95% CI
# =============================================================================
message("=== Step 7: External C-index ===")

boot_c <- function(time, status, predictor, B = 1000, reverse = TRUE) {
  keep <- complete.cases(time, status, predictor) & time > 0 & is.finite(predictor)
  time <- time[keep]; status <- status[keep]; predictor <- predictor[keep]
  n <- length(time)
  if (n < 20) return(list(cindex = NA, ci_low = NA, ci_high = NA, n = n))
  c_obs <- as.numeric(concordance(Surv(time, status) ~ predictor,
                                   reverse = reverse)$concordance)
  boot_vals <- replicate(B, {
    idx <- sample(n, replace = TRUE)
    tryCatch(
      as.numeric(concordance(Surv(time[idx], status[idx]) ~ predictor[idx],
                                reverse = reverse)$concordance),
      error = function(e) NA_real_)
  })
  boot_vals <- boot_vals[!is.na(boot_vals)]
  if (length(boot_vals) < 10) return(list(cindex = c_obs, ci_low = NA, ci_high = NA, n = n))
  ci <- quantile(boot_vals, c(0.025, 0.975))
  list(cindex = c_obs, ci_low = ci[1], ci_high = ci[2], n = n)
}

# Nomogram LP (higher = worse survival -> reverse=TRUE)
lp_b1 <- predict(fit_cph, newdata = cgga_b1_nomo, type = "lp")
lp_b2 <- predict(fit_cph, newdata = cgga_b2_nomo, type = "lp")

# C-index for UIRS z_risk alone (higher = worse -> reverse=TRUE)
zrisk_c_b1 <- boot_c(cgga_b1_nomo$OS_months, cgga_b1_nomo$OS,
                     cgga_b1_nomo$z_risk, B = 1000, reverse = TRUE)
zrisk_c_b2 <- boot_c(cgga_b2_nomo$OS_months, cgga_b2_nomo$OS,
                     cgga_b2_nomo$z_risk, B = 1000, reverse = TRUE)

# C-index for nomogram LP (BOTH directions)
nomo_c_b1_rev <- boot_c(cgga_b1_nomo$OS_months, cgga_b1_nomo$OS,
                        lp_b1, B = 1000, reverse = TRUE)
nomo_c_b2_rev <- boot_c(cgga_b2_nomo$OS_months, cgga_b2_nomo$OS,
                        lp_b2, B = 1000, reverse = TRUE)
nomo_c_b1_norev <- boot_c(cgga_b1_nomo$OS_months, cgga_b1_nomo$OS,
                          lp_b1, B = 1000, reverse = FALSE)
nomo_c_b2_norev <- boot_c(cgga_b2_nomo$OS_months, cgga_b2_nomo$OS,
                          lp_b2, B = 1000, reverse = FALSE)

# Check LP direction
direction_b1 <- if (nomo_c_b1_rev$cindex >= 0.5) "correct" else "REVERSED"
direction_b2 <- if (nomo_c_b2_rev$cindex >= 0.5) "correct" else "REVERSED"

for (tmp in list(
  list(label = "CGGA-b1 UIRS_z", r = zrisk_c_b1),
  list(label = "CGGA-b2 UIRS_z", r = zrisk_c_b2),
  list(label = "CGGA-b1 Nomogram (rev=T)", r = nomo_c_b1_rev),
  list(label = "CGGA-b2 Nomogram (rev=T)", r = nomo_c_b2_rev),
  list(label = "CGGA-b1 Nomogram (rev=F)", r = nomo_c_b1_norev),
  list(label = "CGGA-b2 Nomogram (rev=F)", r = nomo_c_b2_norev)
)) {
  message(sprintf("%s: C=%.4f (95%% CI: %.4f-%.4f) n=%d",
                  tmp$label, tmp$r$cindex, tmp$r$ci_low, tmp$r$ci_high, tmp$r$n))
}

message(sprintf("CGGA-b1 LP direction: %s", direction_b1))
message(sprintf("CGGA-b2 LP direction: %s", direction_b2))

# =============================================================================
# 9. External calibration (slope + calibration-in-the-large)
# =============================================================================
message("=== Step 8: External calibration ===")

external_cal <- function(fit, ext_df, timepoints = c(12, 36, 60), label = "") {
  ext_lp <- predict(fit, newdata = ext_df, type = "lp")
  surv_fn <- Survival(fit)

  results <- list()
  for (tp in timepoints) {
    ext_pred <- surv_fn(tp, lp = ext_lp)

    # Group by predicted risk
    n_grp <- min(5, max(2, floor(nrow(ext_df) / 15)))
    grp <- cut(ext_pred, breaks = quantile(ext_pred,
               probs = seq(0, 1, len = n_grp + 1), na.rm = TRUE),
               include.lowest = TRUE, labels = FALSE)

    cal_data <- data.frame(stringsAsFactors = FALSE)
    for (g in seq_len(n_grp)) {
      idx <- which(!is.na(grp) & grp == g)
      if (length(idx) < 3) next
      pred_mean <- mean(ext_pred[idx], na.rm = TRUE)
      km <- survfit(Surv(ext_df$OS_months[idx], ext_df$OS[idx]) ~ 1)
      km_s <- summary(km, times = tp, extend = TRUE)
      cal_data <- rbind(cal_data, data.frame(
        group = g, n = length(idx),
        pred_mean = pred_mean,
        obs_surv = if (length(km_s$surv) > 0) km_s$surv[1] else NA,
        obs_se = if (length(km_s$std.err) > 0) km_s$std.err[1] else NA,
        stringsAsFactors = FALSE
      ))
    }

    # Calibration slope: Cox(ext_Surv ~ ext_lp)
    cal_slope <- tryCatch({
      cc <- coxph(Surv(ext_df$OS_months, ext_df$OS) ~ ext_lp)
      as.numeric(coef(cc)["ext_lp"])
    }, error = function(e) NA_real_)

    # Calibration-in-the-large: log(O/E)
    cal_large <- tryCatch({
      obs_ev <- sum(ext_df$OS)
      exp_ev <- sum(predict(fit, newdata = ext_df, type = "expected"), na.rm = TRUE)
      if (!is.na(exp_ev) && exp_ev > 0) log(obs_ev / exp_ev) else NA_real_
    }, error = function(e) NA_real_)

    results[[as.character(tp)]] <- list(
      timepoint = tp, cal_data = cal_data,
      cal_slope = cal_slope, cal_large = cal_large,
      n = nrow(ext_df), n_events = sum(ext_df$OS),
      ext_pred = ext_pred, ext_lp = ext_lp
    )

    message(sprintf("  %s %dmo: slope=%.4f, cal-in-large=%.4f, n=%d, events=%d",
                    label, tp, cal_slope, cal_large, nrow(ext_df), sum(ext_df$OS)))
  }
  results
}

cal_b1 <- external_cal(fit_cph, cgga_b1_nomo, label = "CGGA-b1")
cal_b2 <- external_cal(fit_cph, cgga_b2_nomo, label = "CGGA-b2")

# =============================================================================
# 10. Output CSV: internal calibration
# =============================================================================
write.csv(data.frame(
  cindex = c_int, ci_low = ci_int[1], ci_high = ci_int[2],
  n = nrow(nomo_df), n_events = sum(nomo_df$OS),
  timepoint_months = "1/3/5",
  stringsAsFactors = FALSE
), file.path(RES_DIR, "nomogram_standardized_internal.csv"), row.names = FALSE)

# =============================================================================
# 11. Output CSV: external calibration
# =============================================================================
cal_summ <- data.frame(stringsAsFactors = FALSE)
for (nm in names(cal_b1)) {
  r <- cal_b1[[nm]]
  cal_summ <- rbind(cal_summ, data.frame(
    cohort = "CGGA-batch1", timepoint_months = r$timepoint,
    cal_slope = r$cal_slope, cal_in_large = r$cal_large,
    n = r$n, n_events = r$n_events, stringsAsFactors = FALSE
  ))
}
for (nm in names(cal_b2)) {
  r <- cal_b2[[nm]]
  cal_summ <- rbind(cal_summ, data.frame(
    cohort = "CGGA-batch2", timepoint_months = r$timepoint,
    cal_slope = r$cal_slope, cal_in_large = r$cal_large,
    n = r$n, n_events = r$n_events, stringsAsFactors = FALSE
  ))
}
write.csv(cal_summ, file.path(RES_DIR, "nomogram_external_calibration_standardized.csv"),
          row.names = FALSE)

# =============================================================================
# 12. Output CSV: external C-index
# =============================================================================
cindex_df <- data.frame(
  cohort     = c("TCGA", "CGGA-batch1", "CGGA-batch2",
                 "CGGA-batch1", "CGGA-batch2"),
  model      = c("Nomogram_z", "Nomogram_z", "Nomogram_z",
                 "UIRS_z_alone", "UIRS_z_alone"),
  cindex     = c(c_int,
                 nomo_c_b1_rev$cindex, nomo_c_b2_rev$cindex,
                 zrisk_c_b1$cindex, zrisk_c_b2$cindex),
  ci_low     = c(ci_int[1],
                 nomo_c_b1_rev$ci_low, nomo_c_b2_rev$ci_low,
                 zrisk_c_b1$ci_low, zrisk_c_b2$ci_low),
  ci_high    = c(ci_int[2],
                 nomo_c_b1_rev$ci_high, nomo_c_b2_rev$ci_high,
                 zrisk_c_b1$ci_high, zrisk_c_b2$ci_high),
  n          = c(nrow(nomo_df),
                 nomo_c_b1_rev$n, nomo_c_b2_rev$n,
                 zrisk_c_b1$n, zrisk_c_b2$n),
  validation = c("internal", "external", "external", "external", "external"),
  note       = c("",
                 sprintf("LP_direction=%s", direction_b1),
                 sprintf("LP_direction=%s", direction_b2),
                 "", ""),
  stringsAsFactors = FALSE
)
write.csv(cindex_df, file.path(RES_DIR, "nomogram_external_cindex_standardized.csv"),
          row.names = FALSE)

# =============================================================================
# 13. External calibration figures
# =============================================================================
message("=== Step 9: External calibration figures ===")

make_cal_fig <- function(cal_list, cohort, filename) {
  pdf(file.path(FIG_DIR, filename), width = 12, height = 5)
  par(mfrow = c(1, 3), mar = c(5, 5, 4, 2))
  labels <- c("12" = "1-Year", "36" = "3-Year", "60" = "5-Year")
  for (nm in names(cal_list)) {
    r <- cal_list[[nm]]
    cd <- r$cal_data
    if (nrow(cd) > 0 && !all(is.na(cd$obs_surv))) {
      ylim <- range(c(cd$obs_surv - 1.96 * cd$obs_se, cd$obs_surv + 1.96 * cd$obs_se),
                    na.rm = TRUE)
      ylim <- c(max(0, ylim[1] - 0.1), min(1, ylim[2] + 0.1))
      plot(cd$pred_mean, cd$obs_surv,
           xlim = c(0, 1), ylim = c(0, 1),
           xlab = "Predicted Survival (TCGA model)", ylab = "Observed (KM)",
           main = sprintf("%s -- %s", cohort, labels[[nm]]),
           pch = 19, col = if (r$cal_slope > 0) "#00A087" else "#E64B35",
           cex = 1.5)
      abline(0, 1, lty = 2, col = "grey50")
      if (!all(is.na(cd$obs_se))) {
        segments(cd$pred_mean, pmax(0, cd$obs_surv - 1.96 * cd$obs_se),
                 cd$pred_mean, pmin(1, cd$obs_surv + 1.96 * cd$obs_se),
                 col = if (r$cal_slope > 0) "#00A087" else "#E64B35", lwd = 1)
      }
      legend("topleft",
             legend = c(sprintf("Cal. slope: %.2f", r$cal_slope),
                        sprintf("Cal-in-large: %.2f", r$cal_large),
                        sprintf("n=%d", r$n)),
             bty = "n", cex = 0.85)
    } else {
      plot.new(); text(0.5, 0.5, "Insufficient data", cex = 1)
    }
  }
  dev.off()
  message(sprintf("  Saved: %s", filename))
}

make_cal_fig(cal_b1, "CGGA-batch1 (z-standardized)",
             "Fig6_calibration_external_CGGAb1_standardized.pdf")
make_cal_fig(cal_b2, "CGGA-batch2 (z-standardized)",
             "SuppFig_calibration_external_CGGAb2_standardized.pdf")

# =============================================================================
# 14. Verdict + notes
# =============================================================================
message("=== Verdict ===")

# Determine if external calibration improved
old_slope_b1 <- -0.897
old_slope_b2 <- -0.727
new_slope_b1 <- cal_b1[["12"]]$cal_slope
new_slope_b2 <- cal_b2[["12"]]$cal_slope
uirs_alone_b1 <- zrisk_c_b1$cindex
uirs_alone_b2 <- zrisk_c_b2$cindex
nomo_ext_b1   <- nomo_c_b1_rev$cindex

slope_fixed_b1 <- new_slope_b1 > 0
slope_fixed_b2 <- new_slope_b2 > 0
cindex_improved_b1 <- nomo_ext_b1 > uirs_alone_b1

# Write notes
notes <- c(
  "# Nomogram Standardization Fix -- Results Notes",
  "",
  "## Frozen standardization parameters (TCGA training distribution)",
  sprintf("- mean_TCGA = %.6f", mean_TCGA),
  sprintf("- sd_TCGA   = %.6f", sd_TCGA),
  sprintf("- Computed from train_df (n=%d) raw UIRS = sum(beta_g * expr_g)", nrow(train_df)),
  "",
  "## Internal performance (TCGA)",
  sprintf("- C-index: %.4f (95%% CI: %.4f-%.4f)", c_int, ci_int[1], ci_int[2]),
  sprintf("- Samples: %d, Events: %d", nrow(nomo_df), sum(nomo_df$OS)),
  sprintf("- Global PH: p=%.4f", ph_test$table["GLOBAL", "p"]),
  "",
  "## External calibration slopes",
  "| Cohort | Timepoint | Old Slope | New Slope | Status |",
  "|--------|-----------|-----------|-----------|--------|",
  sprintf("| CGGA-b1 | 12mo | %.3f | %.4f | %s |",
          old_slope_b1, new_slope_b1,
          if (slope_fixed_b1) "FIXED (+)" else "STILL NEGATIVE"),
  sprintf("| CGGA-b1 | 36mo | %.3f | %.4f | |",
          old_slope_b1, cal_b1[["36"]]$cal_slope),
  sprintf("| CGGA-b1 | 60mo | %.3f | %.4f | |",
          old_slope_b1, cal_b1[["60"]]$cal_slope),
  sprintf("| CGGA-b2 | 12mo | %.3f | %.4f | %s |",
          old_slope_b2, new_slope_b2,
          if (slope_fixed_b2) "FIXED (+)" else "STILL NEGATIVE"),
  sprintf("| CGGA-b2 | 36mo | %.3f | %.4f | |",
          old_slope_b2, cal_b2[["36"]]$cal_slope),
  sprintf("| CGGA-b2 | 60mo | %.3f | %.4f | |",
          old_slope_b2, cal_b2[["60"]]$cal_slope),
  "",
  "## External C-index comparison",
  "| Cohort | Model | C-index | 95% CI |",
  "|--------|-------|---------|--------|",
  sprintf("| CGGA-b1 | UIRS_z alone | %.4f | %.4f-%.4f |",
          zrisk_c_b1$cindex, zrisk_c_b1$ci_low, zrisk_c_b1$ci_high),
  sprintf("| CGGA-b1 | Nomogram_z | %.4f | %.4f-%.4f |",
          nomo_c_b1_rev$cindex, nomo_c_b1_rev$ci_low, nomo_c_b1_rev$ci_high),
  sprintf("| CGGA-b2 | UIRS_z alone | %.4f | %.4f-%.4f |",
          zrisk_c_b2$cindex, zrisk_c_b2$ci_low, zrisk_c_b2$ci_high),
  sprintf("| CGGA-b2 | Nomogram_z | %.4f | %.4f-%.4f |",
          nomo_c_b2_rev$cindex, nomo_c_b2_rev$ci_low, nomo_c_b2_rev$ci_high),
  "",
  "## VERDICT",
  sprintf("- CGGA-b1 calibration slope: %.4f (old: %.3f) -> %s",
          new_slope_b1, old_slope_b1,
          if (slope_fixed_b1) "POSITIVE -- direction restored" else "STILL NEGATIVE -- fix incomplete"),
  sprintf("- CGGA-b2 calibration slope: %.4f (old: %.3f) -> %s",
          new_slope_b2, old_slope_b2,
          if (slope_fixed_b2) "POSITIVE -- direction restored" else "STILL NEGATIVE -- fix incomplete"),
  sprintf("- LP direction CGGA-b1: %s (C_rev=%.4f)", direction_b1, nomo_c_b1_rev$cindex),
  sprintf("- LP direction CGGA-b2: %s (C_rev=%.4f)", direction_b2, nomo_c_b2_rev$cindex),
  sprintf("- Nomogram_z vs UIRS_z alone: delta-C = %.4f (b1), %.4f (b2)",
          nomo_ext_b1 - uirs_alone_b1,
          nomo_c_b2_rev$cindex - uirs_alone_b2),
  "",
  if (slope_fixed_b1 && slope_fixed_b2) {
    "OVERALL VERDICT: STANDARDIZATION FIXED EXTERNAL CALIBRATION.\nBoth CGGA cohorts now show positive calibration slopes and correct LP direction.\n"
  } else if (slope_fixed_b1 || slope_fixed_b2) {
    "OVERALL VERDICT: PARTIALLY FIXED.\nOne cohort calibrates externally, the other still shows issues.\n"
  } else {
    "OVERALL VERDICT: STANDARDIZATION ALONE DID NOT FIX EXTERNAL CALIBRATION.\nThe nomogram still does not generalize externally. The manuscript should report\nthe UIRS signature alone (external C ~0.73) with the nomogram as internal-only.\n"
  }
)

writeLines(notes, file.path(RES_DIR, "nomogram_standardized_notes.md"))
cat(paste(notes, collapse = "\n"), "\n")

# =============================================================================
# 15. Save fitted objects
# =============================================================================
save(fit_cph, nomo_df, cgga_b1_nomo, cgga_b2_nomo, mean_TCGA, sd_TCGA,
     cal_b1, cal_b2, boot_cidx, ph_test,
     file = file.path(DATA_PROC, "nomogram_standardized_results.RData"))

message("\n=== M3b completed ===")
message(sprintf("TCGA internal C-index: %.4f (%.4f-%.4f)", c_int, ci_int[1], ci_int[2]))
message(sprintf("CGGA-b1 cal slope: %.4f (was %.3f)", new_slope_b1, old_slope_b1))
message(sprintf("CGGA-b2 cal slope: %.4f (was %.3f)", new_slope_b2, old_slope_b2))
message(sprintf("FROZEN: mean=%.4f, sd=%.4f", mean_TCGA, sd_TCGA))
