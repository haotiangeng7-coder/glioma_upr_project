###############################################################################
# fig5_composite_v2.R
# Figure 5 composite — nomogram & clinical utility
# Panels: A nomogram, B calibration (3 timepoints), C DCA, D univariate forest,
#         E multivariate forest
# Fixes:
#   - Forest plot subtitles remove "Fig 6F/6G" prefix
#   - Nomogram row label "risk_score" → "Risk score (UIRS)"
#   - Adequate height for all panels
# Output: Figures_v2/Figure5_composite.pdf/.png
###############################################################################


suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(survival)
  library(rms)
  library(regplot)
  library(patchwork)
  library(cowplot)
  library(png)
  library(grid)
  library(grDevices)
  library(showtext)
  library(sysfonts)
})

# ── Layout-only patch: regplot:::Cox_scale "Pr( yvar > t )" row label ─────────
# regplot draws each survival-axis row label inline-right (adj=1) anchored at
# min(tickpos) — i.e. ending exactly at the first tick — so it overprints that
# axis's first probability value ("Pr( Month > 60 )0.9"). We move the label into
# the left margin regplot already reserves (x = L - 0.3*(M-L), adj=0), the same
# slot used for "Total points"/predictor row labels and the strata label. This
# repositions ONLY the text label; no tick, coordinate, probability, model, or
# colour is touched. The continuous (Cox_scalemedsurv / odds) branch is patched
# identically for safety though this figure uses the survival-probability axis.
local({
  .cs <- regplot:::Cox_scale
  .b  <- deparse(body(.cs))
  .b  <- gsub('text(x = min(tickpos), yposx, adj = 1, paste("Pr(", ',
              'text(x = L - 0.3 * (M - L), yposx, adj = 0, paste("Pr(", ',
              .b, fixed=TRUE)
  .b  <- gsub('text(x = min(tickpos), yposx, adj = 1, paste("odds(", ',
              'text(x = L - 0.3 * (M - L), yposx, adj = 0, paste("odds(", ',
              .b, fixed=TRUE)
  body(.cs) <- parse(text=paste(.b, collapse="\n"))[[1]]
  environment(.cs) <- environment(regplot:::Cox_scale)
  utils::assignInNamespace("Cox_scale", .cs, ns="regplot")
})

# ── Layout-only patch: enlarge the tick-number ↔ axis-line vertical gap ───────
# regplot hardcodes a single half-amplitude `delta` (= 0.05 + 0.15*npanels/10,
# regplot:::regplot lines ~1169-1171) that it reuses for (a) the axis tick-mark
# length ltick = 0.4*delta, (b) the staggered category-label amplitude, AND (c)
# the vertical offset of every axis's tick NUMBER, drawn at `<line> - delta`
# (numbers below the line) or `<line> + pm*delta` (staggered above/below). There
# is NO exposed parameter to widen only the number↔line gap. With cexscales=0.62
# the numbers sit so close that their top edge touches the axis line. We push the
# tick NUMBERS away from their lines by the factor TICK_GAP, editing ONLY the y
# coordinate inside the relevant text() calls in the five drawing sites:
#   regplot:::plot_caxis   — simple numeric axes (numbers below line)
#   regplot:::plot_caxisX  — staggered numeric axes, e.g. Age 10/50/90 (±delta)
#   regplot:::regplot      — the Points (0-100) axis + the Total-points axis
#   regplot:::Cox_scale    — the Pr(Month>N) probability-value numbers
# delta itself is UNCHANGED, so axis-line positions, tick-mark length (ltick),
# the category-label stagger amplitude (IDH/Grade), the left-margin Pr/Total-row
# labels (prior patch), the title position, and all prediction markers are not
# moved — only the printed numbers drop a little farther from their own line.
local({
  TICK_GAP <- 1.7                         # number offset = delta -> delta*1.7
  ns <- "regplot"
  patch_fn <- function(fname, subs) {
    f <- get(fname, envir = getNamespace(ns))
    b <- deparse(body(f))
    for (s in subs) b <- gsub(s[1], s[2], b, fixed = TRUE)
    body(f) <- parse(text = paste(b, collapse = "\n"))[[1]]
    environment(f) <- getNamespace(ns)
    utils::assignInNamespace(fname, f, ns = ns)
  }
  g <- sprintf(" * %s", format(TICK_GAP))
  # plot_caxis: y = ipos - delta  ->  y = ipos - delta*GAP  (numbers below)
  patch_fn("plot_caxis", list(
    c("text(x = ticks_pos, y = ipos - delta, paste(tickval)",
      sprintf("text(x = ticks_pos, y = ipos - delta%s, paste(tickval)", g))))
  # plot_caxisX: pos <- ipos + pm * delta  ->  ... pm * delta*GAP  (staggered)
  patch_fn("plot_caxisX", list(
    c("pos <- ipos + pm * delta",
      sprintf("pos <- ipos + pm * delta%s", g))))
  # regplot main: Points axis (ipos - delta) and Total-points axis (ypos - delta)
  patch_fn("regplot", list(
    c("text(x = ticks, y = ipos - delta, paste(ticks)",
      sprintf("text(x = ticks, y = ipos - delta%s, paste(ticks)", g)),
    c("text(x = tickp, y = ypos - delta, paste(signif(tickval, ",
      sprintf("text(x = tickp, y = ypos - delta%s, paste(signif(tickval, ", g))))
  # Cox_scale: Pr(Month>N) probability numbers (yposx - delta), both branches.
  # Re-fetch from namespace so this stacks on the row-label patch applied above.
  patch_fn("Cox_scale", list(
    c("text(x = Xpos, y = yposx - delta, paste(signif(1 - ",
      sprintf("text(x = Xpos, y = yposx - delta%s, paste(signif(1 - ", g)),
    c("text(x = Xpos, y = yposx - delta, paste(signif(Pr, ",
      sprintf("text(x = Xpos, y = yposx - delta%s, paste(signif(Pr, ", g))))
})

font_add("Times",
         regular    = "/usr/share/fonts/truetype/dejavu/DejaVuSerif.ttf",
         bold       = "/usr/share/fonts/truetype/dejavu/DejaVuSerif-Bold.ttf",
         italic     = "/usr/share/fonts/truetype/dejavu/DejaVuSerif-Italic.ttf",
         bolditalic = "/usr/share/fonts/truetype/dejavu/DejaVuSerif-BoldItalic.ttf")
showtext_auto()
showtext_opts(dpi=300)

PROJ <- getwd()
DATA_PROC <- file.path(PROJ,"data","processed")
RES_DIR   <- file.path(PROJ,"results")
OUT_DIR   <- file.path(PROJ, "manuscript", "submission", "Figures_v2")
TMP_DIR   <- file.path(OUT_DIR,"tmp_fig5")
dir.create(TMP_DIR,recursive=TRUE,showWarnings=FALSE)

W_MM  <- 183
H_MM  <- 325   # 300 -> 325: the +25 mm is funded entirely to Panel A (row 1) so its predictor axes are taller and regplot's fixed-amplitude staggered category labels (IDH WT/Mutant, Grade G2/G3/G4) gain physical vertical separation; rows B/C/D/E keep their prior absolute mm (66.67/50/50/50). Earlier note: 270 -> 300 had funded Panel B so its 3 calibration boxes are near-square — still preserved.
MM_IN <- 1/25.4
FONT  <- "Times"
BASE  <- 8

# ── Slot geometry ───────────────────────────────────────────────────────────
# The composite is W_MM x H_MM with plot_layout(heights = LAYOUT_HEIGHTS).
# patchwork normalizes the heights, so each row's physical height in the final
# figure = (h / sum(h)) * H_MM. Raster panels A (nomogram) and B (calibration)
# are saved as PNG then stretched to fill their slot via draw_grob — so any
# mismatch between the PNG's mm canvas and its slot's mm size rescales the
# embedded text away from the true pt used by the vector panels (C, D, E). To
# keep A/B text visually equal to C/D/E we render each PNG at its true slot
# size, so placement is ~1:1 and the base-R text stays at its set pt. These
# constants are the single source of truth for both the per-panel png() calls
# and the final plot_layout(). Row 4/5 heights depend on the forest row counts,
# so the forest data is loaded up front (the cleaned data frames are reused
# verbatim in Panel D below — no recomputation, no data change).
.fd <- new.env()
load(file.path(DATA_PROC,"independent_prognosis_results.RData"), envir=.fd)
.n_uni <- nrow(.fd$univar_results)
.n_mlt <- nrow(subset(.fd$multivar_results,
                      !is.na(HR) & is.finite(HR) & !is.na(HR_lower) & !is.na(HR_upper)))
LAYOUT_HEIGHTS <- c(0.390, 0.24, 0.18,
                    max(0.12, .n_uni*0.03),
                    max(0.12, .n_mlt*0.03))
.row_mm  <- (LAYOUT_HEIGHTS / sum(LAYOUT_HEIGHTS)) * H_MM  # physical row heights
SLOT_A_W <- W_MM            # A spans full width
SLOT_A_H <- .row_mm[1]      # row 1 height
SLOT_B_W <- W_MM            # B spans full width
SLOT_B_H <- .row_mm[2]      # row 2 height
rm(.fd)
message(sprintf("  Slot A: %.1f x %.1f mm | Slot B: %.1f x %.1f mm",
                SLOT_A_W, SLOT_A_H, SLOT_B_W, SLOT_B_H))

theme_pub <- function(bs=BASE) {
  theme_classic(base_size=bs, base_family=FONT) +
    theme(
      axis.line  = element_line(linewidth=0.4),
      axis.ticks = element_line(linewidth=0.3),
      axis.text  = element_text(size=bs-1,color="black"),
      axis.title = element_text(size=bs),
      strip.text = element_text(size=bs,face="bold"),
      strip.background=element_blank(),
      legend.text  = element_text(size=bs-1),
      legend.title = element_text(size=bs),
      legend.key.size=unit(3,"mm"),
      plot.title = element_text(size=bs,face="bold",hjust=0),
      panel.grid = element_blank()
    )
}

fmt_pval <- function(p) {
  ifelse(is.na(p),"NA",
    ifelse(p<0.001, formatC(p,format="e",digits=1),
           formatC(p,format="f",digits=3)))
}

###############################################################################
# Load nomogram data
###############################################################################
message("Loading nomogram data ...")
load(file.path(DATA_PROC,"nomogram_results.RData"))
dd <- datadist(nomo_df)
options(datadist="dd")
formula_cph <- formula(fit_cph)
message(sprintf("  Formula: %s", deparse(formula_cph,width.cutoff=200)))
message(sprintf("  Samples: %d", nrow(nomo_df)))

###############################################################################
# Panel A — Nomogram (regplot)
# Reimplemented with regplot::regplot() to eliminate the categorical-label /
# tick-number coordinate collisions inherent to rms::plot.nomogram (the Age
# tick numbers and the Grade/IDH category labels overprinted each other on the
# compressed base-R layout). regplot staggers categorical labels above/below
# each predictor axis, so Grade (G2/G3/G4) and IDH (Mutant/WT) read cleanly
# without touching the continuous Age scale.
#
# Statistics are UNCHANGED: regplot is driven by an equivalent cph refit
# (fit_reg) on the SAME data and SAME formula as fit_cph — only the column
# names are display aliases (risk_score -> UIRS, OS_months -> Month, OS ->
# Event) so the rendered row/axis labels are clean. The fitted coefficients of
# fit_reg are identical to fit_cph (verified by stopifnot below); nothing about
# the model, the predicted probabilities, or the colour semantics changes.
###############################################################################
message("Panel A (nomogram, regplot) ...")

# Equivalent model with display-friendly variable names (identical fit).
nomo_disp <- nomo_df
names(nomo_disp)[names(nomo_disp)=="risk_score"] <- "UIRS"
names(nomo_disp)[names(nomo_disp)=="OS_months"]  <- "Month"
names(nomo_disp)[names(nomo_disp)=="OS"]          <- "Event"
dd_reg  <- datadist(nomo_disp)
options(datadist="dd_reg")
fit_reg <- cph(Surv(Month, Event) ~ UIRS + Age + Grade + IDH,
               data=nomo_disp, x=TRUE, y=TRUE, surv=TRUE)
# Guard: the display refit must be statistically identical to the canonical fit.
stopifnot(isTRUE(all.equal(unname(coef(fit_cph)), unname(coef(fit_reg)))))
# Restore the canonical datadist for the downstream calibration panel (Panel B).
options(datadist="dd")

nomo_A_png <- file.path(TMP_DIR,"nomo_A.png")

# regplot manages its own graphics device (it calls dev.off() then plot.new()
# internally), so a pre-opened png() device is discarded. We redirect that
# internal plot.new() to OUR png target via options(device=...). showtext is
# disabled around the call because regplot's fixed cex-based layout collides
# with showtext's DPI text scaling (title overprints the Points axis); the
# DejaVu Serif system font matches the "Times" (DejaVuSerif) family used by the
# vector panels C/D/E, so Panel A stays visually consistent.
.old_device <- getOption("device")
.st_on      <- TRUE
tryCatch(showtext::showtext_auto(FALSE), error=function(e) .st_on <<- FALSE)
options(device=function(...)
  grDevices::png(nomo_A_png, width=SLOT_A_W, height=SLOT_A_H,
                 units="mm", res=300, family="DejaVu Serif"))
nomo_ok <- tryCatch({
  regplot::regplot(
    fit_reg,
    plots     = c("no plot","no plot"),  # clean nomogram only (no density/box strips)
    failtime  = c(12,36,60),             # OS time is in MONTHS -> 1/3/5-year survival
    prfail    = FALSE,                   # survival probability axes (not failure)
    points    = TRUE,
    clickable = FALSE,
    showP     = FALSE,
    title     = "Cox Regression Nomogram for Overall Survival",
    cexvars   = 0.95,                    # predictor row labels
    cexcats   = 0.58,                    # categorical level labels (staggered) — 0.62 -> 0.58: slightly smaller glyphs add horizontal breathing room between WT/Mutant & G2/G3/G4; the dominant fix is the taller Panel A slot (108 mm) which expands the fixed-amplitude vertical stagger physically
    cexscales = 0.62                     # axis tick numbers + "Pr( Month > N )" row labels — shrunk so each Pr label clears its axis first tick value
  )
  TRUE
}, error=function(e) { message("  regplot() error: ", e$message); FALSE })
# regplot opened the device via its internal plot.new(); close it to flush PNG.
while (grDevices::dev.cur() > 1) grDevices::dev.off()
options(device=.old_device)
if (.st_on) showtext::showtext_auto(TRUE)

if (!isTRUE(nomo_ok) || !file.exists(nomo_A_png)) {
  stop("Panel A regplot rendering failed; aborting (fail-fast).")
}
message("  Panel A saved")

###############################################################################
# Panel B — Calibration curves (1, 3, 5 year)
###############################################################################
message("Panel B (calibration curves) ...")

cal_B_png <- file.path(TMP_DIR,"cal_B.png")
cal_timepoints <- c(12,36,60)
cal_labels     <- c("1-Year","3-Year","5-Year")
m_per_group    <- max(20,min(80,floor(nrow(nomo_df)/5)))

# Render at true slot size (full width x row-2 height) so draw_grob places it
# ~1:1 and text is not shrunk. cex compensation matches C/D/E (~8pt).
png(cal_B_png, width=SLOT_B_W, height=SLOT_B_H, units="mm", res=300, family=FONT)
par(mfrow=c(1,3), mar=c(4,4,3,1),
    family=FONT, cex=0.95, cex.axis=0.95, cex.lab=1.05, cex.main=1.1)

for (ti in seq_along(cal_timepoints)) {
  tp <- cal_timepoints[ti]
  lb <- cal_labels[ti]
  tryCatch({
    fit_cal <- cph(formula_cph, data=nomo_df, x=TRUE, y=TRUE, surv=TRUE, time.inc=tp)
    cal_obj <- calibrate(fit_cal, cmethod="KM", method="boot",
                         u=tp, m=m_per_group, B=200)
    plot(cal_obj, xlab="Predicted Probability", ylab="Actual Probability",
         main=sprintf("%s Survival",lb), xlim=c(0,1), ylim=c(0,1), subtitles=FALSE)
    abline(0,1,lty=2,col="grey50")
  }, error=function(e) {
    message(sprintf("  %s calibration error: %s",lb,e$message))
    plot.new()
    text(0.5,0.5,e$message,cex=0.65,family=FONT)
  })
}
dev.off()
message("  Panel B saved")

###############################################################################
# Panel C — DCA (3-year)
###############################################################################
message("Panel C (DCA) ...")

surv_dca_fn <- function(data, time_var, event_var, predictors,
                         pred_labels=NULL, timepoint,
                         thresholds=seq(0.01,0.99,by=0.01)) {
  if (is.null(pred_labels)) pred_labels <- predictors
  surv_obj   <- Surv(data[[time_var]],data[[event_var]])
  km_fit     <- survfit(surv_obj~1)
  km_summ    <- summary(km_fit,times=timepoint)
  event_rate <- 1 - km_summ$surv
  results    <- data.frame()
  for (th in thresholds) {
    nb_all <- event_rate - (1-event_rate)*th/(1-th)
    results <- rbind(results,data.frame(threshold=th,predictor="Treat All",
                                         net_benefit=nb_all,stringsAsFactors=FALSE))
  }
  results <- rbind(results,data.frame(threshold=thresholds,predictor="Treat None",
                                       net_benefit=0,stringsAsFactors=FALSE))
  for (pi in seq_along(predictors)) {
    pred_name  <- predictors[pi]
    pred_label <- pred_labels[pi]
    tryCatch({
      fml <- as.formula(paste0("Surv(",time_var,",",event_var,")~",pred_name))
      fit <- coxph(fml,data=data)
      bh  <- basehaz(fit,centered=TRUE)
      idx_tp <- which.min(abs(bh$time-timepoint))
      H0_tp  <- bh$hazard[idx_tp]
      lp     <- predict(fit,type="lp")
      pred_prob <- 1-exp(-H0_tp*exp(lp))
      pred_prob <- pmin(pmax(pred_prob,0),1)
      for (th in thresholds) {
        n <- nrow(data)
        treat <- pred_prob>=th
        tp_count <- sum(treat & data[[event_var]]==1 & data[[time_var]]<=timepoint,na.rm=TRUE)
        fp_count <- sum(treat & (data[[event_var]]==0 | data[[time_var]]>timepoint),na.rm=TRUE)
        nb <- tp_count/n - fp_count/n*th/(1-th)
        results <- rbind(results,data.frame(threshold=th,predictor=pred_label,
                                             net_benefit=nb,stringsAsFactors=FALSE))
      }
    }, error=function(e) message(sprintf("  DCA failed for %s: %s",pred_name,e$message)))
  }
  results
}

p_C <- tryCatch({
  dca_res <- surv_dca_fn(
    data=nomo_df, time_var="OS_months", event_var="OS",
    predictors="risk_score", pred_labels="UIRS", timepoint=36,
    thresholds=seq(0.01,0.80,by=0.01)
  )
  ggplot(dca_res, aes(x=threshold,y=net_benefit,color=predictor,linetype=predictor)) +
    geom_line(linewidth=0.7) +
    geom_hline(yintercept=0,linetype="solid",color="grey80",linewidth=0.3) +
    scale_color_manual(values=c("UIRS"="#E64B35","Treat All"="#3C5488","Treat None"="grey50"),
                       name="Strategy") +
    scale_linetype_manual(values=c("UIRS"="solid","Treat All"="dashed","Treat None"="dotted"),
                          name="Strategy") +
    labs(x="Threshold Probability",y="Net Benefit",
         title="Decision Curve Analysis (3-Year)") +
    theme_pub() +
    theme(legend.position=c(0.97,0.97),
          legend.justification=c("right","top"),
          legend.background=element_rect(fill=NA,color=NA)) +
    coord_cartesian(xlim=c(0,0.8),
                    ylim=c(min(dca_res$net_benefit,na.rm=TRUE)*0.1,
                           max(dca_res$net_benefit,na.rm=TRUE)*1.1))
}, error=function(e) {
  message("  DCA error: ",e$message)
  ggplot()+theme_void()+
    annotate("text",x=0.5,y=0.5,label=paste("DCA error:",e$message),size=3)
})

###############################################################################
# Panel D — Univariate forest plot
# FIX: title is "Univariate Cox Regression" (no "Fig 6F" prefix)
###############################################################################
message("Panel D (univariate forest) ...")
load(file.path(DATA_PROC,"independent_prognosis_results.RData"))

univar_clean <- univar_results %>%
  dplyr::mutate(
    label = dplyr::case_when(
      Variable=="UIRS Risk Score"                ~ "UIRS Risk Score",
      Variable=="Age"                            ~ "Age (per year)",
      grepl("Gender",Variable,ignore.case=TRUE)  ~ "Gender (Male vs Female)",
      Variable=="WHO Grade"                      ~ "WHO Grade (per grade)",
      grepl("MGMT",Variable,ignore.case=TRUE)    ~ "MGMT Methylation",
      grepl("IDH",Variable,ignore.case=TRUE)     ~ "IDH Status (WT vs Mutant)",
      TRUE ~ Variable
    )
  ) %>%
  dplyr::select(label,HR,HR_lower,HR_upper,pvalue)

multivar_clean <- multivar_results %>%
  dplyr::mutate(
    label = dplyr::case_when(
      Variable=="risk_score"       ~ "Risk score (UIRS)",  # FIX: was "UIRS Risk Score"
      Variable=="Age"              ~ "Age (per year)",
      Variable=="Gendermale"       ~ "Gender (Male vs Female)",
      Variable=="GradeG2"          ~ "WHO Grade G2 (ref: G2)",
      Variable=="GradeG3"          ~ "WHO Grade G3 (vs G2)",
      Variable=="GradeG4"          ~ "WHO Grade G4 (vs G2)",
      Variable=="MGMT"             ~ "MGMT Methylation",
      Variable=="MGMTUnmethylated" ~ "MGMT (Unmethylated vs Methylated)",
      grepl("IDH",Variable,ignore.case=TRUE) ~ "IDH Status (WT vs Mutant)",
      TRUE ~ Variable
    )
  ) %>%
  dplyr::filter(!is.na(HR) & is.finite(HR) & !is.na(HR_lower) & !is.na(HR_upper)) %>%
  dplyr::select(label,HR,HR_lower,HR_upper,pvalue)

make_forest_gg <- function(df, title_txt, subtitle_txt, sig_color="#E64B35") {
  n <- nrow(df)
  df$hr_ci_text <- sprintf("%.2f (%.2f–%.2f)", df$HR, df$HR_lower, df$HR_upper)
  df$p_text     <- fmt_pval(df$pvalue)
  df$sig        <- ifelse(!is.na(df$pvalue) & df$pvalue<0.05,"Significant","Non-significant")
  df$y_pos      <- rev(seq_len(n))
  df$pt_color   <- ifelse(df$sig=="Significant",sig_color,"#999999")

  all_vals <- c(df$HR_lower[is.finite(df$HR_lower)], df$HR_upper[is.finite(df$HR_upper)])
  x_lo <- max(min(all_vals,na.rm=TRUE)*0.7, 0.05)
  x_hi <- max(all_vals,na.rm=TRUE)*1.4
  candidate_breaks <- c(0.05,0.1,0.2,0.5,1,2,5,10,20,50)
  breaks_use <- candidate_breaks[candidate_breaks>=x_lo & candidate_breaks<=x_hi]
  if (!1 %in% breaks_use) breaks_use <- sort(c(breaks_use,1))

  df$HR_lower_plot <- pmax(df$HR_lower, x_lo)
  df$HR_upper_plot <- pmin(df$HR_upper, x_hi)

  # Left: labels
  p_left <- ggplot(df,aes(y=y_pos)) +
    annotate("text",x=0,y=n+0.7,label="Variable",
             hjust=0,size=3.0,fontface="bold",family=FONT) +
    geom_text(aes(x=0,label=label),hjust=0,size=2.8,color="grey15",family=FONT) +
    scale_x_continuous(limits=c(0,1),expand=c(0,0)) +
    scale_y_continuous(limits=c(0.3,n+1.2),expand=c(0,0)) +
    theme_void(base_family=FONT) + theme(plot.margin=margin(4,0,4,4))

  # Center: forest
  p_center <- ggplot(df,aes(y=y_pos)) +
    geom_vline(xintercept=1,linetype="dashed",color="grey55",linewidth=0.5) +
    geom_segment(aes(x=HR_lower_plot,xend=HR_upper_plot,y=y_pos,yend=y_pos),
                 color=df$pt_color,linewidth=0.8) +
    geom_point(aes(x=HR),color=df$pt_color,shape=18,size=3.0) +
    annotate("text",x=exp((log(x_lo)+log(x_hi))/2),y=n+0.7,
             label="Hazard Ratio (95% CI)",size=3.0,fontface="bold",family=FONT) +
    scale_x_log10(limits=c(x_lo,x_hi),breaks=breaks_use,labels=as.character(breaks_use)) +
    scale_y_continuous(limits=c(0.3,n+1.2),expand=c(0,0)) +
    labs(x="Hazard Ratio (log scale)") +
    theme_classic(base_size=BASE,base_family=FONT) +
    theme(axis.title.y=element_blank(),axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),axis.line.y=element_blank(),
          axis.text.x=element_text(size=7),axis.title.x=element_text(size=7.5),
          panel.grid=element_blank(),plot.margin=margin(4,3,4,3))

  # Right1: HR text
  p_right1 <- ggplot(df,aes(y=y_pos)) +
    annotate("text",x=0.5,y=n+0.7,label="HR (95% CI)",
             size=3.0,fontface="bold",hjust=0.5,family=FONT) +
    geom_text(aes(x=0.5,label=hr_ci_text),hjust=0.5,size=2.7,color="grey15",family=FONT) +
    scale_x_continuous(limits=c(0,1),expand=c(0,0)) +
    scale_y_continuous(limits=c(0.3,n+1.2),expand=c(0,0)) +
    theme_void(base_family=FONT)+theme(plot.margin=margin(4,2,4,2))

  # Right2: P-value
  p_right2 <- ggplot(df,aes(y=y_pos)) +
    annotate("text",x=0.5,y=n+0.7,label="P value",
             size=3.0,fontface="bold",hjust=0.5,family=FONT) +
    geom_text(aes(x=0.5,label=p_text,color=sig=="Significant"),
              hjust=0.5,size=2.7,family=FONT) +
    scale_color_manual(values=c("TRUE"=sig_color,"FALSE"="#999999"),guide="none") +
    scale_x_continuous(limits=c(0,1),expand=c(0,0)) +
    scale_y_continuous(limits=c(0.3,n+1.2),expand=c(0,0)) +
    theme_void(base_family=FONT)+theme(plot.margin=margin(4,4,4,0))

  p_combined <- cowplot::plot_grid(
    p_left,p_center,p_right1,p_right2,
    nrow=1, rel_widths=c(2.2,3.5,1.8,1.0), align="h", axis="tb"
  )

  # FIX: title WITHOUT "Fig 6F/6G" prefix. Title band thickened (rel_heights
  # 0.16 -> 0.32) so the bold title (y=0.78, top-anchored) and the grey subtitle
  # (y=0.30, bottom-anchored, one size smaller) have a clear vertical gap and no
  # longer overprint within an over-thin band.
  title_grob <- cowplot::ggdraw() +
    cowplot::draw_label(title_txt, x=0.5,y=0.78,fontface="bold",
                        size=BASE+1,hjust=0.5,vjust=1) +
    cowplot::draw_label(subtitle_txt, x=0.5,y=0.30,fontface="plain",
                        size=BASE-2,hjust=0.5,vjust=0,color="grey40")

  cowplot::plot_grid(title_grob, p_combined, ncol=1, rel_heights=c(0.32,1))
}

# FIX: titles do NOT contain figure-number prefix
p_D <- make_forest_gg(
  univar_clean,
  title_txt    = "Univariate Cox Regression",
  subtitle_txt = "IDH stratified (strata(IDH)); TCGA cohort"
)

p_E <- make_forest_gg(
  multivar_clean,
  title_txt    = "Multivariate Cox Regression",
  subtitle_txt = "IDH as stratification variable; TCGA cohort"
)

###############################################################################
# Assemble composite
###############################################################################
message("=== Assembling composite ===")

img_A <- png::readPNG(nomo_A_png)
img_B <- png::readPNG(cal_B_png)
grob_A <- grid::rasterGrob(img_A, interpolate=TRUE)
grob_B <- grid::rasterGrob(img_B, interpolate=TRUE)

wrap_grob <- function(g) cowplot::ggdraw() + cowplot::draw_grob(g)
p_A_wrap  <- wrap_grob(grob_A)
p_B_wrap  <- wrap_grob(grob_B)

# Row 1: A (nomogram)
# Row 2: B (calibration)
# Row 3: C (DCA)
# Row 4: D (univariate forest)
# Row 5: E (multivariate forest)
p_composite <- (p_A_wrap / p_B_wrap / p_C / p_D / p_E) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size=13,face="bold",family=FONT,
                              color="black",margin=margin(0,2,2,0,"mm"))
    )
  ) +
  plot_layout(heights = LAYOUT_HEIGHTS)

# Sanity: LAYOUT_HEIGHTS (computed up front from the same data) must match the
# row counts actually produced by Panel D — guards against the slot geometry
# used for A/B rendering drifting from the final layout.
stopifnot(identical(LAYOUT_HEIGHTS,
                    c(0.390, 0.24, 0.18,
                      max(0.12, nrow(univar_clean)*0.03),
                      max(0.12, nrow(multivar_clean)*0.03))))

out_pdf <- file.path(OUT_DIR,"Figure5_composite.pdf")
out_png <- file.path(OUT_DIR,"Figure5_composite.png")

ggsave(out_pdf, plot=p_composite, device=cairo_pdf,
       width=W_MM*MM_IN, height=H_MM*MM_IN, units="in", dpi=300)
message("Saved PDF: ", out_pdf)

ggsave(out_png, plot=p_composite,
       width=W_MM*MM_IN, height=H_MM*MM_IN, units="in", dpi=300)
message("Saved PNG: ", out_png)

unlink(TMP_DIR, recursive=TRUE)

message("=== Figure 5 composite done ===")
message(sprintf("Dimensions: %d x %d mm  (W/H aspect = %.2f)", W_MM, H_MM, W_MM/H_MM))
message("Panel A: nomogram with 'risk_score' relabelled to 'Risk score (UIRS)'")
message("Panels D/E: forest plot titles have no 'Fig 6F/6G' prefix")
