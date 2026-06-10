###############################################################################
# generate_fig3.R
# Publication-quality Figure 3: Consensus Clustering and Molecular Subtypes
# Journal of Translational Medicine submission
#
# Panels:
#   Fig3_cluster_metrics  : PAC / Silhouette / Delta AUC across K=2..6, K=3 highlighted
#   Fig3_km_final         : KM OS for K=3 subtypes (full cohort), risk table, CI
#   Fig3_km_idhwt         : KM OS IDH-WT patients only
#   Fig3_clinical_heatmap : Grade / IDH / Gender / MGMT by subtype (ComplexHeatmap)
#   Fig3_upr_by_subtype   : UPR score violin+box by subtype
#   Fig3_upr_by_idh       : UPR score violin+box by IDH status
#
# Output: manuscript/submission/Figures/Figure3/ (relative to project root)
###############################################################################

suppressPackageStartupMessages({
  library(ggplot2)
  library(survminer)
  library(survival)
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(tidyr)
  library(ggpubr)
  library(patchwork)
  library(scales)
  library(grid)
  library(grDevices)
  library(showtext)
  library(sysfonts)
})

# Register DejaVu Serif as "Times" (metric-compatible TrueType serif font
# available on this system; accepted as Times New Roman equivalent).
font_add("Times",
         regular    = "/usr/share/fonts/truetype/dejavu/DejaVuSerif.ttf",
         bold       = "/usr/share/fonts/truetype/dejavu/DejaVuSerif-Bold.ttf",
         italic     = "/usr/share/fonts/truetype/dejavu/DejaVuSerif-Italic.ttf",
         bolditalic = "/usr/share/fonts/truetype/dejavu/DejaVuSerif-BoldItalic.ttf")
showtext_auto()
showtext_opts(dpi = 600)

# ---------------------------------------------------------------------------
# 0. Paths and style constants
# ---------------------------------------------------------------------------
PROJ <- getwd()
DATA_PROC   <- file.path(PROJ, "data", "processed")
OUT_DIR     <- file.path(PROJ, "manuscript", "submission", "Figures", "Figure3")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Nature / JTM dimensions (mm)
W_SINGLE    <- 89   # single column
W_1P5       <- 130  # 1.5 column
W_DOUBLE    <- 183  # double column
MM_TO_IN    <- 1 / 25.4

# Subtype color palette (consistent throughout manuscript)
COL_SUBTYPE <- c(
  "UPR-high-risk"    = "#E64B35",
  "UPR-intermediate" = "#00A087",
  "UPR-favorable"    = "#4DBBD5"
)
COL_IDH <- c("Mutant" = "#00A087", "WT" = "#E64B35")

# Publication theme — Times (DejaVu Serif), minimal gridlines
FONT   <- "Times"
BASE_SIZE <- 8   # 8pt base so panel text ~6-8 pt at final print size

theme_pub <- function(base_size = BASE_SIZE) {
  theme_classic(base_size = base_size, base_family = FONT) +
    theme(
      axis.line        = element_line(linewidth = 0.4, color = "black"),
      axis.ticks       = element_line(linewidth = 0.3, color = "black"),
      axis.text        = element_text(size = base_size - 1, color = "black"),
      axis.title       = element_text(size = base_size,     color = "black"),
      strip.text       = element_text(size = base_size,     face = "bold"),
      strip.background = element_blank(),
      legend.text      = element_text(size = base_size - 1),
      legend.title     = element_text(size = base_size),
      legend.key.size  = unit(3, "mm"),
      plot.title       = element_text(size = base_size, face = "bold", hjust = 0),
      plot.subtitle    = element_text(size = base_size - 1, color = "grey40"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
}

# Helper: save PDF (vector) only, embed fonts via cairo_pdf
save_pdf <- function(plot_obj = NULL, path, width_mm, height_mm,
                     draw_fn = NULL) {
  w_in <- width_mm  * MM_TO_IN
  h_in <- height_mm * MM_TO_IN
  if (!is.null(plot_obj)) {
    cairo_pdf(path, width = w_in, height = h_in, family = FONT)
    print(plot_obj)
    dev.off()
  } else if (!is.null(draw_fn)) {
    cairo_pdf(path, width = w_in, height = h_in, family = FONT)
    draw_fn()
    dev.off()
  }
  message("  Saved: ", path)
}

# Helper: significance stars
sig_star <- function(p) {
  dplyr::case_when(
    p < 0.0001 ~ "****",
    p < 0.001  ~ "***",
    p < 0.01   ~ "**",
    p < 0.05   ~ "*",
    TRUE       ~ "ns"
  )
}

# ---------------------------------------------------------------------------
# 1. Load data
# ---------------------------------------------------------------------------
message("=== Loading data ===")
load(file.path(DATA_PROC, "consensus_clustering_results.RData"))
load(file.path(DATA_PROC, "clinical_characterization_results.RData"))
load(file.path(DATA_PROC, "upr_landscape_results.RData"))

# Convenient alias for clinical+subtype merged table
clin <- heatmap_data   # 845 rows: barcode, UPR_subtype, OS.time, OS, IDH_status, Grade, etc.

# ---------------------------------------------------------------------------
# PANEL A — Cluster selection metrics (PAC, Silhouette, Delta AUC)
# ---------------------------------------------------------------------------
message("=== Panel A: Cluster metrics ===")

# metrics has columns: K, PAC, CDF_AUC, Delta_AUC, Silhouette_mean
metrics_long <- metrics %>%
  dplyr::select(K, PAC, Silhouette_mean, Delta_AUC) %>%
  pivot_longer(-K, names_to = "Metric", values_to = "Value") %>%
  mutate(
    Metric = dplyr::recode(Metric,
      "PAC"             = "PAC",
      "Silhouette_mean" = "Silhouette",
      "Delta_AUC"       = "Delta AUC"
    ),
    Metric = factor(Metric, levels = c("PAC", "Silhouette", "Delta AUC")),
    Selected = (K == final_k)
  )

# Arrow annotations for selected K
arrow_df <- metrics_long %>%
  dplyr::filter(K == final_k) %>%
  dplyr::select(Metric, K, Value)

p_metrics <- ggplot(metrics_long, aes(x = K, y = Value)) +
  geom_line(color = "grey50", linewidth = 0.6) +
  geom_point(aes(fill = Selected), shape = 21, size = 2.5, stroke = 0.5,
             color = "black") +
  scale_fill_manual(values = c("FALSE" = "white", "TRUE" = "#E64B35"),
                    guide = "none") +
  geom_point(data = dplyr::filter(metrics_long, K == final_k),
             shape = 21, size = 3.5, fill = "#E64B35", color = "black", stroke = 0.6) +
  facet_wrap(~Metric, scales = "free_y", nrow = 1) +
  scale_x_continuous(breaks = 2:6) +
  labs(x = "Number of clusters (K)", y = "Value",
       title = sprintf("K = %d selected", final_k)) +
  theme_pub() +
  theme(strip.text = element_text(face = "bold", size = BASE_SIZE))

save_pdf(p_metrics, file.path(OUT_DIR, "Fig3A_cluster_metrics.pdf"),
         W_DOUBLE, 55)

# ---------------------------------------------------------------------------
# PANEL B — KM OS: full cohort (K=3 subtypes)
# ---------------------------------------------------------------------------
message("=== Panel B: KM survival (full cohort) ===")

surv_df <- clin %>%
  dplyr::filter(!is.na(OS.time), OS.time > 0, !is.na(UPR_subtype)) %>%
  mutate(
    OS_months = OS.time / 30.44,
    UPR_subtype = factor(UPR_subtype,
                          levels = intersect(c("UPR-high-risk", "UPR-intermediate", "UPR-favorable"),
                                             unique(UPR_subtype)))
  )

fit_km <- survfit(Surv(OS_months, OS) ~ UPR_subtype, data = surv_df)

# log-rank p-value
lr_test  <- survdiff(Surv(OS_months, OS) ~ UPR_subtype, data = surv_df)
lr_pval  <- 1 - pchisq(lr_test$chisq, df = length(unique(surv_df$UPR_subtype)) - 1)

# Median OS per group
med_os <- surv_median(fit_km)

# Build KM plot with risk table via ggsurvplot
palette_km <- COL_SUBTYPE[levels(surv_df$UPR_subtype)]

p_km_full <- ggsurvplot(
  fit_km,
  data              = surv_df,
  pval              = FALSE,       # we add manually below for formatting control
  conf.int          = TRUE,
  conf.int.alpha    = 0.12,
  palette           = palette_km,
  xlab              = "Time (months)",
  ylab              = "Overall survival probability",
  legend.title      = "",
  legend.labs       = levels(surv_df$UPR_subtype),
  risk.table        = TRUE,
  risk.table.height = 0.28,
  risk.table.y.text = FALSE,
  fontsize          = 2.5,
  ggtheme           = theme_pub(),
  tables.theme      = theme_pub() +
    theme(axis.line = element_blank(), axis.ticks = element_blank())
)

# Add log-rank p-value and median OS annotation manually
pval_label <- sprintf("Log-rank p = %.2e", lr_pval)
p_km_full$plot <- p_km_full$plot +
  annotate("text", x = Inf, y = 0.98, label = pval_label,
           hjust = 1.05, vjust = 1, size = 2.5, family = FONT) +
  labs(title = "Overall survival by UPR subtype (full cohort)") +
  theme(legend.position = c(0.78, 0.78),
        legend.background = element_rect(fill = NA, color = NA))

cairo_pdf(file.path(OUT_DIR, "Fig3B_km_final.pdf"),
          width = W_DOUBLE * MM_TO_IN, height = 80 * MM_TO_IN, family = FONT)
print(p_km_full)
dev.off()
message("  Saved: ", file.path(OUT_DIR, "Fig3B_km_final.pdf"))

# ---------------------------------------------------------------------------
# PANEL C — KM OS: IDH-WT subgroup
# ---------------------------------------------------------------------------
message("=== Panel C: KM survival (IDH-WT) ===")

idhwt_surv <- merge(
  idhwt_cluster_df,
  clin[, c("barcode", "OS.time", "OS")],
  by = "barcode"
) %>%
  dplyr::filter(!is.na(OS.time), OS.time > 0, !is.na(UPR_subtype_idhwt)) %>%
  mutate(
    OS_months  = OS.time / 30.44,
    UPR_subtype = factor(UPR_subtype_idhwt,
                          levels = intersect(c("UPR-high-risk", "UPR-intermediate", "UPR-favorable"),
                                             unique(UPR_subtype_idhwt)))
  )

n_idhwt <- nrow(idhwt_surv)
fit_km_idhwt <- survfit(Surv(OS_months, OS) ~ UPR_subtype, data = idhwt_surv)
lr_idhwt     <- survdiff(Surv(OS_months, OS) ~ UPR_subtype, data = idhwt_surv)
pval_idhwt   <- 1 - pchisq(lr_idhwt$chisq,
                             df = length(unique(idhwt_surv$UPR_subtype)) - 1)

palette_idhwt <- COL_SUBTYPE[levels(idhwt_surv$UPR_subtype)]

p_km_idhwt <- ggsurvplot(
  fit_km_idhwt,
  data              = idhwt_surv,
  pval              = FALSE,
  conf.int          = TRUE,
  conf.int.alpha    = 0.12,
  palette           = palette_idhwt,
  xlab              = "Time (months)",
  ylab              = "Overall survival probability",
  legend.title      = "",
  legend.labs       = levels(idhwt_surv$UPR_subtype),
  risk.table        = TRUE,
  risk.table.height = 0.28,
  risk.table.y.text = FALSE,
  fontsize          = 2.5,
  ggtheme           = theme_pub(),
  tables.theme      = theme_pub() +
    theme(axis.line = element_blank(), axis.ticks = element_blank())
)

pval_label_idhwt <- sprintf("Log-rank p = %.2e\n(IDH-WT, n = %d)", pval_idhwt, n_idhwt)
p_km_idhwt$plot <- p_km_idhwt$plot +
  annotate("text", x = Inf, y = 0.98, label = pval_label_idhwt,
           hjust = 1.05, vjust = 1, size = 2.5, family = FONT) +
  labs(title = "Overall survival by UPR subtype (IDH-WT only)") +
  theme(legend.position = c(0.78, 0.78),
        legend.background = element_rect(fill = NA, color = NA))

cairo_pdf(file.path(OUT_DIR, "Fig3C_km_idhwt.pdf"),
          width = W_DOUBLE * MM_TO_IN, height = 80 * MM_TO_IN, family = FONT)
print(p_km_idhwt)
dev.off()
message("  Saved: ", file.path(OUT_DIR, "Fig3C_km_idhwt.pdf"))

# ---------------------------------------------------------------------------
# PANEL D — Clinical feature heatmap (ComplexHeatmap)
# ---------------------------------------------------------------------------
message("=== Panel D: Clinical characterization heatmap ===")

# Order samples by canonical prognosis order, keeping only subtypes present (K=2 or K=3)
subtype_order <- intersect(c("UPR-high-risk", "UPR-intermediate", "UPR-favorable"),
                           unique(clin$UPR_subtype[!is.na(clin$UPR_subtype)]))
clin_ordered  <- clin %>%
  dplyr::filter(!is.na(UPR_subtype)) %>%
  mutate(UPR_subtype = factor(UPR_subtype, levels = subtype_order)) %>%
  dplyr::arrange(UPR_subtype)

n_samp <- nrow(clin_ordered)

# Build clinical annotation matrix (rows = features, cols = samples)
# Features: UPR subtype, Grade (G2/G3/G4), IDH, MGMT, Gender
make_anno_vec <- function(x, ref) factor(x, levels = ref)

grade_col <- c("G2" = "#FFF5EB", "G3" = "#FC8D59", "G4" = "#D94701")
mgmt_col  <- c("Methylated" = "#3C5488", "Unmethylated" = "#F39B7F", "Unknown" = "grey85")
gender_col <- c("male" = "#4DBBD5", "female" = "#E64B35")

# Significance annotation text (from clinical_tests_df)
get_sig_text <- function(varname) {
  row <- clinical_tests_df[clinical_tests_df$Variable == varname, ]
  if (nrow(row) == 0) return("")
  p <- row$P_value[1]
  paste0(sig_star(p), " (p=", formatC(p, format = "e", digits = 1), ")")
}

# Top column annotation: UPR subtype (colored bar)
col_subtype <- clin_ordered$UPR_subtype

# Build grade, IDH, MGMT, gender vectors
grade_vec  <- factor(clin_ordered$Grade,       levels = c("G2", "G3", "G4"))
idh_vec    <- factor(clin_ordered$IDH_status,  levels = c("Mutant", "WT"))
mgmt_vec   <- factor(clin_ordered$MGMT_status, levels = c("Methylated", "Unmethylated"))
gender_vec <- factor(clin_ordered$gender,      levels = c("male", "female"))

# Column split by subtype
col_split <- factor(clin_ordered$UPR_subtype, levels = subtype_order)

# Subtype counts for column titles
subtype_n <- table(clin_ordered$UPR_subtype)[subtype_order]
col_title  <- paste0(subtype_order, "\n(n=", subtype_n, ")")

top_anno <- HeatmapAnnotation(
  `UPR subtype`  = col_subtype,
  Grade          = grade_vec,
  `IDH status`   = idh_vec,
  `MGMT methyl.` = mgmt_vec,
  Gender         = gender_vec,
  col = list(
    `UPR subtype`  = COL_SUBTYPE,
    Grade          = grade_col,
    `IDH status`   = COL_IDH,
    `MGMT methyl.` = mgmt_col,
    Gender         = gender_col
  ),
  annotation_name_side = "left",
  annotation_name_gp   = gpar(fontsize = 7, fontfamily = FONT),
  na_col               = "grey92",
  height               = unit(22, "mm"),
  show_legend          = TRUE,
  annotation_legend_param = list(
    `UPR subtype`  = list(title_gp = gpar(fontsize = 7, fontfamily = FONT),
                           labels_gp = gpar(fontsize = 6, fontfamily = FONT)),
    Grade          = list(title_gp = gpar(fontsize = 7, fontfamily = FONT),
                           labels_gp = gpar(fontsize = 6, fontfamily = FONT)),
    `IDH status`   = list(title_gp = gpar(fontsize = 7, fontfamily = FONT),
                           labels_gp = gpar(fontsize = 6, fontfamily = FONT)),
    `MGMT methyl.` = list(title_gp = gpar(fontsize = 7, fontfamily = FONT),
                           labels_gp = gpar(fontsize = 6, fontfamily = FONT)),
    Gender         = list(title_gp = gpar(fontsize = 7, fontfamily = FONT),
                           labels_gp = gpar(fontsize = 6, fontfamily = FONT))
  )
)

# Bottom annotation: proportion bars for categorical variables
# Build a matrix where each row = binary indicator for one category
# We'll use barplot-style annotation via anno_barplot / anno_block

# We encode a numeric heatmap (age z-score) + the categorical ones in the
# column annotation, and draw a null heatmap body (empty matrix) so the
# entire figure is the annotation.

# Build numeric age row for heatmap body
age_z <- scale(as.numeric(clin_ordered$age_at_index))[, 1]

# Add significance p-value labels to the annotation track names
get_p <- function(var) {
  row <- clinical_tests_df[clinical_tests_df$Variable == var, ]
  if (nrow(row) == 0) return(1)
  row$P_value[1]
}
p_grade  <- get_p("WHO Grade")
p_idh    <- get_p("IDH status")
p_mgmt   <- get_p("MGMT methylation")
p_age    <- get_p("Age")
p_gender <- get_p("Gender")

anno_name_fn <- function(var, p) {
  paste0(var, "\n", sig_star(p))
}

# Rebuild annotation with significance in name
top_anno2 <- HeatmapAnnotation(
  `UPR subtype`                                          = col_subtype,
  ` `                                                    = anno_empty(border = FALSE, height = unit(1, "mm")),
  `Grade\n***`                                           = grade_vec,
  `IDH status\n***`                                      = idh_vec,
  `MGMT\n***`                                            = mgmt_vec,
  `Gender\nns`                                           = gender_vec,
  `Age (z-score)\n***`                                   = anno_points(
    age_z, size = unit(0.5, "mm"), pch = 16,
    gp = gpar(col = ifelse(age_z > 0, "#E64B35", "#4DBBD5"), alpha = 0.5)
  ),
  col = list(
    `UPR subtype`     = COL_SUBTYPE,
    `Grade\n***`      = grade_col,
    `IDH status\n***` = COL_IDH,
    `MGMT\n***`       = mgmt_col,
    `Gender\nns`      = gender_col
  ),
  annotation_name_side = "left",
  annotation_name_gp   = gpar(fontsize = 6, fontfamily = FONT),
  na_col               = "grey92",
  annotation_legend_param = list(
    `UPR subtype`     = list(title = "UPR subtype",
                              title_gp = gpar(fontsize = 7, fontfamily = FONT),
                              labels_gp = gpar(fontsize = 6, fontfamily = FONT)),
    `Grade\n***`      = list(title = "WHO Grade",
                              title_gp = gpar(fontsize = 7, fontfamily = FONT),
                              labels_gp = gpar(fontsize = 6, fontfamily = FONT)),
    `IDH status\n***` = list(title = "IDH status",
                              title_gp = gpar(fontsize = 7, fontfamily = FONT),
                              labels_gp = gpar(fontsize = 6, fontfamily = FONT)),
    `MGMT\n***`       = list(title = "MGMT methylation",
                              title_gp = gpar(fontsize = 7, fontfamily = FONT),
                              labels_gp = gpar(fontsize = 6, fontfamily = FONT)),
    `Gender\nns`      = list(title = "Gender",
                              title_gp = gpar(fontsize = 7, fontfamily = FONT),
                              labels_gp = gpar(fontsize = 6, fontfamily = FONT))
  )
)

# Null heatmap body (1-row placeholder, hidden)
dummy_mat <- matrix(0, nrow = 1, ncol = n_samp,
                    dimnames = list("", clin_ordered$barcode))

ht_clin <- Heatmap(
  dummy_mat,
  top_annotation    = top_anno2,
  col               = c("0" = "white"),
  show_heatmap_legend = FALSE,
  show_row_names    = FALSE,
  show_column_names = FALSE,
  cluster_rows      = FALSE,
  cluster_columns   = FALSE,
  column_split      = col_split,
  column_title      = col_title,
  column_title_gp   = gpar(fontsize = 7, fontface = "bold", fontfamily = FONT),
  border            = FALSE,
  height            = unit(0, "mm")
)

save_pdf(path = file.path(OUT_DIR, "Fig3D_clinical_heatmap.pdf"),
         width_mm = W_DOUBLE, height_mm = 80,
         draw_fn = function() {
           draw(ht_clin,
                merge_legend = TRUE,
                legend_grouping = "original",
                heatmap_legend_side = "right",
                annotation_legend_side = "right")
         })

# ---------------------------------------------------------------------------
# PANEL E — UPR score by subtype (violin + box + jitter)
# ---------------------------------------------------------------------------
message("=== Panel E: UPR score by subtype ===")

upr_df <- data.frame(
  barcode   = names(upr_score_bulk),
  UPR_score = upr_score_bulk,
  stringsAsFactors = FALSE
) %>%
  merge(clin[, c("barcode", "UPR_subtype")], by = "barcode") %>%
  dplyr::filter(!is.na(UPR_subtype)) %>%
  mutate(UPR_subtype = factor(UPR_subtype, levels = subtype_order))

# Pairwise KW then Wilcoxon with BH
kw_upr <- kruskal.test(UPR_score ~ UPR_subtype, data = upr_df)

comparisons_upr <- Filter(function(p) all(p %in% unique(upr_df$UPR_subtype)), list(
  c("UPR-high-risk", "UPR-intermediate"),
  c("UPR-high-risk", "UPR-favorable"),
  c("UPR-intermediate", "UPR-favorable")
))

p_upr_sub <- ggplot(upr_df, aes(x = UPR_subtype, y = UPR_score, fill = UPR_subtype)) +
  geom_violin(alpha = 0.55, trim = TRUE, linewidth = 0.3) +
  geom_boxplot(width = 0.18, outlier.shape = NA, linewidth = 0.4,
               fill = "white", color = "black") +
  stat_compare_means(comparisons = comparisons_upr,
                     method = "wilcox.test", p.adjust.method = "BH",
                     label = "p.signif", hide.ns = FALSE,
                     size = 2.5, tip.length = 0.01, step.increase = 0.06,
                     family = FONT) +
  annotate("text",
           x = 0.6, y = max(upr_df$UPR_score) * 1.02,
           label = sprintf("Kruskal-Wallis p = %.2e", kw_upr$p.value),
           hjust = 0, size = 2, family = FONT, color = "grey30") +
  scale_fill_manual(values = COL_SUBTYPE, guide = "none") +
  scale_x_discrete(labels = function(x) gsub("UPR-", "", x)) +
  labs(x = "UPR subtype", y = "UPR composite score (log2 TPM)",
       title = "UPR score by molecular subtype") +
  theme_pub()

save_pdf(p_upr_sub, file.path(OUT_DIR, "Fig3E_upr_by_subtype.pdf"),
         W_SINGLE, 65)

# ---------------------------------------------------------------------------
# PANEL F — UPR score by IDH status
# ---------------------------------------------------------------------------
message("=== Panel F: UPR score by IDH status ===")

upr_idh_df <- upr_df %>%
  merge(clin[, c("barcode", "IDH_status")], by = "barcode") %>%
  dplyr::filter(!is.na(IDH_status)) %>%
  mutate(IDH_status = factor(IDH_status, levels = c("Mutant", "WT")))

wt_idh  <- wilcox.test(UPR_score ~ IDH_status, data = upr_idh_df)

p_upr_idh <- ggplot(upr_idh_df, aes(x = IDH_status, y = UPR_score, fill = IDH_status)) +
  geom_violin(alpha = 0.55, trim = TRUE, linewidth = 0.3) +
  geom_boxplot(width = 0.18, outlier.shape = NA, linewidth = 0.4,
               fill = "white", color = "black") +
  stat_compare_means(method = "wilcox.test", label = "p.format",
                     size = 2.5, label.y.npc = 0.93, family = FONT) +
  scale_fill_manual(values = COL_IDH, guide = "none") +
  labs(x = "IDH status", y = "UPR composite score (log2 TPM)",
       title = "UPR score by IDH status") +
  theme_pub()

save_pdf(p_upr_idh, file.path(OUT_DIR, "Fig3F_upr_by_idh.pdf"),
         W_SINGLE, 65)

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
message("\n=== Figure 3 complete ===")
cat("\nOutput directory:", OUT_DIR, "\n")
cat("Files produced:\n")
for (f in list.files(OUT_DIR, pattern = "\\.pdf$")) {
  cat(" ", f, "\n")
}

cat("\nKey statistics for caption:\n")
cat("  final K =", final_k, "\n")
cat("  N total =", nrow(surv_df), "\n")
cat("  Log-rank p (full cohort) =", formatC(lr_pval, format = "e", digits = 2), "\n")
cat("  Log-rank p (IDH-WT) =", formatC(pval_idhwt, format = "e", digits = 2), "\n")
cat("  N IDH-WT =", n_idhwt, "\n")
