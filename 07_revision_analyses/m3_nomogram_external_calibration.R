###############################################################################
# M3: External nomogram calibration on CGGA cohorts
# Outputs: results/nomogram_external_calibration.csv,
#          results/nomogram_external_cindex.csv,
#          figures/Fig6_calibration_external_CGGAb1.pdf,
#          figures/SuppFig_calibration_external_CGGAb2.pdf
###############################################################################

suppressPackageStartupMessages({
  library(survival)
  library(rms)
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
load(file.path(DATA_PROC, "independent_prognosis_results.RData"))

cgga1_val <- coxph_results$cgga1_val_df
cgga2_val <- coxph_results$cgga2_val_df

# =============================================================================
# 2. Build frozen nomogram Cox model on TCGA
# =============================================================================
message("=== Building TCGA nomogram model ===")

nomo_df <- analysis_df[, c("OS_months", "OS", "risk_score", "Age", "Grade", "IDH")]
nomo_df <- nomo_df[complete.cases(nomo_df), ]
nomo_df <- nomo_df[nomo_df$OS_months > 0, ]
nomo_df$Grade <- factor(nomo_df$Grade, levels = c("G2", "G3", "G4"))
nomo_df$IDH   <- factor(nomo_df$IDH, levels = c("Mutant", "WT"))

message(sprintf("TCGA nomogram: %d samples, %d events", nrow(nomo_df), sum(nomo_df$OS)))
message(sprintf("TCGA risk_score: mean=%.1f, sd=%.1f, range=[%.1f, %.1f]",
                mean(nomo_df$risk_score), sd(nomo_df$risk_score),
                min(nomo_df$risk_score), max(nomo_df$risk_score)))

dd <- datadist(nomo_df)
options(datadist = "dd")
fit_cph <- cph(Surv(OS_months, OS) ~ risk_score + Age + Grade + IDH,
               data = nomo_df, x = TRUE, y = TRUE, surv = TRUE, time.inc = 36)

# TCGA internal performance
c_tcga <- (fit_cph$stats["Dxy"] + 1) / 2
ph_test <- cox.zph(fit_cph)
message(sprintf("TCGA C-index: %.4f, PH global p: %.4f", c_tcga, ph_test$table["GLOBAL", "p"]))
message("TCGA coefficients:")
for (nm in names(coef(fit_cph))) {
  message(sprintf("  %s: %.4f", nm, coef(fit_cph)[nm]))
}

# =============================================================================
# 3. Prepare CGGA nomogram frames
# =============================================================================
message("=== Preparing CGGA nomogram frames ===")

prepare_cgga <- function(val_df, raw_clinical, batch_name) {
  # Match val_df to clinical by (OS.time, OS_status, IDH_status)
  v_key <- paste(val_df$time, val_df$status, val_df$IDH_status)
  c_key <- paste(raw_clinical$OS.time, raw_clinical$OS_status, raw_clinical$IDH_status)

  m <- match(v_key, c_key)
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
  c_idx <- m[!is.na(m)]

  os_months <- val_df$time[keep_v] / 30.4375
  os_status <- val_df$status[keep_v]
  risk_score <- val_df$risk_score[keep_v]
  age <- raw_clinical$Age[c_idx]

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

  df <- data.frame(OS_months = os_months, OS = os_status, risk_score = risk_score,
                   Age = age, Grade = grade, IDH = idh, stringsAsFactors = FALSE)
  df <- df[complete.cases(df) & df$OS_months > 0, ]

  message(sprintf("  %s: %d complete, %d events", batch_name, nrow(df), sum(df$OS)))
  message(sprintf("  %s: risk_score mean=%.2f, sd=%.2f", batch_name,
                  mean(df$risk_score), sd(df$risk_score)))
  df
}

cgga_b1_nomo <- prepare_cgga(cgga1_val, cgga_b1$clinical, "CGGA-batch1")
cgga_b2_nomo <- prepare_cgga(cgga2_val, cgga_b2$clinical, "CGGA-batch2")

# =============================================================================
# 4. External C-index with bootstrap CI
#    Convention: higher LP/risk_score = worse survival -> use reverse=TRUE
# =============================================================================
message("=== External C-index ===")

boot_c <- function(time, status, predictor, B = 1000, reverse = TRUE) {
  keep <- complete.cases(time, status, predictor) & time > 0 & is.finite(predictor)
  time <- time[keep]; status <- status[keep]; predictor <- predictor[keep]
  n <- length(time)
  if (n < 20) return(list(cindex = NA, ci_low = NA, ci_high = NA))
  c_obs <- as.numeric(concordance(Surv(time, status) ~ predictor, reverse = reverse)$concordance)
  boot_vals <- replicate(B, {
    idx <- sample(n, replace = TRUE)
    tryCatch(
      as.numeric(concordance(Surv(time[idx], status[idx]) ~ predictor[idx], reverse = reverse)$concordance),
      error = function(e) NA_real_)
  })
  boot_vals <- boot_vals[!is.na(boot_vals)]
  if (length(boot_vals) < 10) return(list(cindex = c_obs, ci_low = NA, ci_high = NA))
  ci <- quantile(boot_vals, c(0.025, 0.975))
  list(cindex = c_obs, ci_low = ci[1], ci_high = ci[2])
}

# Nomogram LP (from TCGA-trained model)
lp_b1 <- predict(fit_cph, newdata = cgga_b1_nomo, type = "lp")
lp_b2 <- predict(fit_cph, newdata = cgga_b2_nomo, type = "lp")

# UIRS risk_score (from locked model)
risk_b1 <- cgga_b1_nomo$risk_score
risk_b2 <- cgga_b2_nomo$risk_score

# C-index for risk_score (higher = worse, correct direction)
uirs_c_b1 <- boot_c(cgga_b1_nomo$OS_months, cgga_b1_nomo$OS, risk_b1, B = 1000, reverse = TRUE)
uirs_c_b2 <- boot_c(cgga_b2_nomo$OS_months, cgga_b2_nomo$OS, risk_b2, B = 1000, reverse = TRUE)

# C-index for nomogram LP: report BOTH directions to show reversal if present
# For LP, higher should = worse survival (same convention). If LP is reversed
# on CGGA, reverse=TRUE will give C<0.5, revealing the miscalibration.
nomo_c_b1_rev <- boot_c(cgga_b1_nomo$OS_months, cgga_b1_nomo$OS, lp_b1, B = 1000, reverse = TRUE)
nomo_c_b2_rev <- boot_c(cgga_b2_nomo$OS_months, cgga_b2_nomo$OS, lp_b2, B = 1000, reverse = TRUE)
nomo_c_b1_norev <- boot_c(cgga_b1_nomo$OS_months, cgga_b1_nomo$OS, lp_b1, B = 1000, reverse = FALSE)
nomo_c_b2_norev <- boot_c(cgga_b2_nomo$OS_months, cgga_b2_nomo$OS, lp_b2, B = 1000, reverse = FALSE)

message(sprintf("CGGA-b1 UIRS C (higher=worse): %.4f (%.4f-%.4f)",
                uirs_c_b1$cindex, uirs_c_b1$ci_low, uirs_c_b1$ci_high))
message(sprintf("CGGA-b2 UIRS C (higher=worse): %.4f (%.4f-%.4f)",
                uirs_c_b2$cindex, uirs_c_b2$ci_low, uirs_c_b2$ci_high))
message(sprintf("CGGA-b1 Nomogram LP C (reverse=TRUE): %.4f (%.4f-%.4f)",
                nomo_c_b1_rev$cindex, nomo_c_b1_rev$ci_low, nomo_c_b1_rev$ci_high))
message(sprintf("CGGA-b1 Nomogram LP C (reverse=FALSE): %.4f (%.4f-%.4f)",
                nomo_c_b1_norev$cindex, nomo_c_b1_norev$ci_low, nomo_c_b1_norev$ci_high))
message(sprintf("CGGA-b2 Nomogram LP C (reverse=TRUE): %.4f (%.4f-%.4f)",
                nomo_c_b2_rev$cindex, nomo_c_b2_rev$ci_low, nomo_c_b2_rev$ci_high))
message(sprintf("CGGA-b2 Nomogram LP C (reverse=FALSE): %.4f (%.4f-%.4f)",
                nomo_c_b2_norev$cindex, nomo_c_b2_norev$ci_low, nomo_c_b2_norev$ci_high))

# Check LP direction on CGGA
if (nomo_c_b1_rev$cindex < 0.5) {
  message("WARNING: Nomogram LP direction REVERSED on CGGA-b1. Higher LP = BETTER survival.")
  message("  This indicates the nomogram FAILS to generalize — risk predictions are inverted.")
}
if (nomo_c_b2_rev$cindex < 0.5) {
  message("WARNING: Nomogram LP direction REVERSED on CGGA-b2. Higher LP = BETTER survival.")
}

# =============================================================================
# 5. External calibration
# =============================================================================
message("=== External calibration ===")

external_cal <- function(fit, ext_df, timepoints = c(12, 36, 60), label = "") {
  ext_lp <- predict(fit, newdata = ext_df, type = "lp")
  surv_fn <- Survival(fit)

  results <- list()
  for (tp in timepoints) {
    ext_pred <- surv_fn(tp, lp = ext_lp)

    # Group by predicted risk
    n_grp <- min(5, max(2, floor(nrow(ext_df) / 15)))
    grp <- cut(ext_pred, breaks = quantile(ext_pred, probs = seq(0, 1, len = n_grp + 1), na.rm = TRUE),
               include.lowest = TRUE, labels = FALSE)

    cal_data <- data.frame(stringsAsFactors = FALSE)
    for (g in seq_len(n_grp)) {
      idx <- which(grp == g)
      if (length(idx) < 3) next
      pred_mean <- mean(ext_pred[idx], na.rm = TRUE)
      km <- survfit(Surv(ext_df$OS_months[idx], ext_df$OS[idx]) ~ 1)
      km_s <- summary(km, times = tp, extend = TRUE)
      cal_data <- rbind(cal_data, data.frame(
        group = g, n = length(idx),
        pred_mean = pred_mean,
        obs_surv = if (length(km_s$surv) > 0) km_s$surv else NA,
        obs_se = if (length(km_s$std.err) > 0) km_s$std.err else NA,
        stringsAsFactors = FALSE
      ))
    }

    # Calibration slope: coxph(Surv ~ lp) on external data
    # slope ~1 = good calibration; negative = LP direction reversed
    cal_slope <- tryCatch({
      cc <- coxph(Surv(ext_df$OS_months, ext_df$OS) ~ ext_lp)
      as.numeric(coef(cc)["ext_lp"])
    }, error = function(e) NA_real_)

    # Calibration-in-the-large: log(O/E ratio)
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
  }
  results
}

cal_b1 <- external_cal(fit_cph, cgga_b1_nomo, label = "CGGA-batch1")
cal_b2 <- external_cal(fit_cph, cgga_b2_nomo, label = "CGGA-batch2")

for (tp_chr in names(cal_b1)) {
  r <- cal_b1[[tp_chr]]
  message(sprintf("CGGA-b1 %d-mo: slope=%.3f, cal-in-large=%.3f",
                  r$timepoint, r$cal_slope, r$cal_large))
}
for (tp_chr in names(cal_b2)) {
  r <- cal_b2[[tp_chr]]
  message(sprintf("CGGA-b2 %d-mo: slope=%.3f, cal-in-large=%.3f",
                  r$timepoint, r$cal_slope, r$cal_large))
}

# =============================================================================
# 6. Output CSVs
# =============================================================================
message("=== Writing CSVs ===")

cal_summary <- data.frame(stringsAsFactors = FALSE)
for (nm in names(cal_b1)) {
  r <- cal_b1[[nm]]
  cal_summary <- rbind(cal_summary, data.frame(
    cohort = "CGGA-batch1", timepoint_months = r$timepoint,
    cal_slope = r$cal_slope, cal_in_large = r$cal_large,
    n = r$n, n_events = r$n_events, stringsAsFactors = FALSE
  ))
}
for (nm in names(cal_b2)) {
  r <- cal_b2[[nm]]
  cal_summary <- rbind(cal_summary, data.frame(
    cohort = "CGGA-batch2", timepoint_months = r$timepoint,
    cal_slope = r$cal_slope, cal_in_large = r$cal_large,
    n = r$n, n_events = r$n_events, stringsAsFactors = FALSE
  ))
}
write.csv(cal_summary, file.path(RES_DIR, "nomogram_external_calibration.csv"), row.names = FALSE)

cindex_df <- data.frame(
  cohort     = c("TCGA", "CGGA-batch1", "CGGA-batch2",
                 "CGGA-batch1", "CGGA-batch2"),
  model      = c("Nomogram", "Nomogram", "Nomogram",
                 "UIRS_alone", "UIRS_alone"),
  cindex     = c(c_tcga, nomo_c_b1_rev$cindex, nomo_c_b2_rev$cindex,
                 uirs_c_b1$cindex, uirs_c_b2$cindex),
  ci_low     = c(c_tcga, nomo_c_b1_rev$ci_low, nomo_c_b2_rev$ci_low,
                 uirs_c_b1$ci_low, uirs_c_b2$ci_low),
  ci_high    = c(c_tcga, nomo_c_b1_rev$ci_high, nomo_c_b2_rev$ci_high,
                 uirs_c_b1$ci_high, uirs_c_b2$ci_high),
  validation = c("internal", "external", "external", "external", "external"),
  note       = c("", "LP_direction_reversed", "LP_direction_reversed", "", ""),
  stringsAsFactors = FALSE
)
write.csv(cindex_df, file.path(RES_DIR, "nomogram_external_cindex.csv"), row.names = FALSE)

# =============================================================================
# 7. Calibration plots
# =============================================================================
message("=== Calibration figures ===")

make_cal_fig <- function(cal_list, cohort, filename) {
  pdf(file.path(FIG_DIR, filename), width = 12, height = 5)
  par(mfrow = c(1, 3), mar = c(5, 5, 4, 2))
  labels <- c("12" = "1-Year", "36" = "3-Year", "60" = "5-Year")
  for (nm in names(cal_list)) {
    r <- cal_list[[nm]]
    cd <- r$cal_data
    if (nrow(cd) > 0 && !all(is.na(cd$obs_surv))) {
      plot(cd$pred_mean, cd$obs_surv,
           xlim = c(0, 1), ylim = c(0, 1),
           xlab = "Predicted Survival (TCGA model)", ylab = "Observed (KM)",
           main = sprintf("%s -- %s", cohort, labels[[nm]]),
           pch = 19, col = "#E64B35", cex = 1.5)
      abline(0, 1, lty = 2, col = "grey50")
      if (!all(is.na(cd$obs_se))) {
        segments(cd$pred_mean, pmax(0, cd$obs_surv - 1.96 * cd$obs_se),
                 cd$pred_mean, pmin(1, cd$obs_surv + 1.96 * cd$obs_se),
                 col = "#E64B35", lwd = 1)
      }
      legend("topleft",
             legend = c(sprintf("Cal. slope: %.2f", r$cal_slope),
                        sprintf("Cal-in-large: %.2f", r$cal_large)),
             bty = "n", cex = 0.85)
    } else {
      plot.new(); text(0.5, 0.5, "Insufficient data", cex = 1)
    }
  }
  dev.off()
  message(sprintf("  Saved: %s", filename))
}

make_cal_fig(cal_b1, "CGGA-batch1 (external)", "Fig6_calibration_external_CGGAb1.pdf")
make_cal_fig(cal_b2, "CGGA-batch2 (external)", "SuppFig_calibration_external_CGGAb2.pdf")

# =============================================================================
# 8. Summary
# =============================================================================
message("=== M3 completed ===")
cat(sprintf("TCGA (internal):           C=%.4f\n", c_tcga))
cat(sprintf("CGGA-b1 Nomogram LP (ext): C=%.4f (%.4f-%.4f) [reversed direction!]\n",
            nomo_c_b1_rev$cindex, nomo_c_b1_rev$ci_low, nomo_c_b1_rev$ci_high))
cat(sprintf("CGGA-b1 UIRS alone (ext):  C=%.4f (%.4f-%.4f)\n",
            uirs_c_b1$cindex, uirs_c_b1$ci_low, uirs_c_b1$ci_high))
cat(sprintf("CGGA-b1 cal slopes: 12mo=%.3f, 36mo=%.3f, 60mo=%.3f\n",
            cal_b1[["12"]]$cal_slope, cal_b1[["36"]]$cal_slope, cal_b1[["60"]]$cal_slope))
cat(sprintf("CGGA-b2 cal slopes: 12mo=%.3f, 36mo=%.3f, 60mo=%.3f\n",
            cal_b2[["12"]]$cal_slope, cal_b2[["36"]]$cal_slope, cal_b2[["60"]]$cal_slope))
cat("\nCAVEAT: TCGA risk_score (mean=%.1f) and CGGA risk_score (mean=%.2f) are on",
    "different absolute scales, causing the nomogram LP to be dominated by Age/Grade/IDH",
    "rather than risk_score on CGGA. The nomogram does NOT generalize to CGGA.\n")
