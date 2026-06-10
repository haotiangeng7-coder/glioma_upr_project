###############################################################################
# fig2_composite_v2.R
# Figure 2 composite — consensus clustering / molecular subtypes
# Panels: A cluster metrics, B KM full, C KM IDH-WT, D clinical heatmap,
#         E UPR by subtype, F UPR by IDH
# Fixes: enlarged KM fonts; enlarged clinical-heatmap annotation-track labels
# Output: Figures_v2/Figure2_composite.pdf/.png
#
# Strategy: ComplexHeatmap and ggsurvplot cannot be combined via patchwork.
#           We save each panel to a temp PDF then assemble into one PDF by
#           reading back as rasterGrob at 150 dpi (sufficient for composite).
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
  library(cowplot)
  library(grid)
  library(png)
  library(grDevices)
  library(showtext)
  library(sysfonts)
})

font_add("Times",
         regular    = "/usr/share/fonts/truetype/dejavu/DejaVuSerif.ttf",
         bold       = "/usr/share/fonts/truetype/dejavu/DejaVuSerif-Bold.ttf",
         italic     = "/usr/share/fonts/truetype/dejavu/DejaVuSerif-Italic.ttf",
         bolditalic = "/usr/share/fonts/truetype/dejavu/DejaVuSerif-BoldItalic.ttf")
showtext_auto()
showtext_opts(dpi = 300)

PROJ <- getwd()
DATA_PROC <- file.path(PROJ, "data", "processed")
OUT_DIR   <- file.path(PROJ, "manuscript", "submission", "Figures_v2")
TMP_DIR   <- file.path(OUT_DIR, "tmp_fig2")
dir.create(TMP_DIR, recursive = TRUE, showWarnings = FALSE)

W_MM   <- 183
H_MM   <- 260
MM_IN  <- 1 / 25.4
FONT   <- "Times"
BASE   <- 8

# ── Slot geometry ───────────────────────────────────────────────────────────
# The composite is W_MM x H_MM with plot_layout(heights = LAYOUT_HEIGHTS).
# patchwork normalizes the heights, so each row's physical height in the final
# figure = (h / sum(h)) * H_MM. Raster panels (B, C KM; D heatmap) are saved as
# PNG then stretched to fill their slot via draw_grob — so any mismatch between
# the PNG's mm canvas and its slot's mm size rescales the embedded text away from
# the 8pt used by the vector panels (A, E, F). To keep B/C/D text visually equal
# to A/E/F we render each PNG at (close to) its true slot size, so placement is
# ~1:1 and 8pt stays 8pt. These constants are the single source of truth for both
# the per-panel png() calls and the final plot_layout().
LAYOUT_HEIGHTS <- c(0.15, 0.28, 0.22, 0.22)
.row_mm <- (LAYOUT_HEIGHTS / sum(LAYOUT_HEIGHTS)) * H_MM   # physical row heights
SLOT_BC_W <- W_MM / 2          # B|C each occupy half the figure width
SLOT_BC_H <- .row_mm[2]        # row 2 height
SLOT_D_W  <- W_MM              # D spans full width
SLOT_D_H  <- .row_mm[3]        # row 3 height

COL_SUBTYPE <- c(
  "UPR-high-risk"    = "#E64B35",
  "UPR-intermediate" = "#00A087",
  "UPR-favorable"    = "#4DBBD5"
)
COL_IDH <- c("Mutant" = "#00A087", "WT" = "#E64B35")

# KM panels (B, C) are rendered by survminer's print.ggsurvplot, which packs the
# survival plot + risk table into the PNG canvas: the plot region occupies only
# ~60% of the canvas height, so nominal-8pt theme text renders visually ~0.6x of
# the vector panels (A, E, F). KM_BASE compensates so KM text matches A/E/F at the
# composite scale; the mm-based fontsize (risk table) and annotation size are
# scaled by the same factor (see KM_FS / KM_ANN below).
KM_BASE <- 13            # compensated base_size for KM ggtheme (vs BASE=8 vector)
# Risk-table number size is decoupled from KM_BASE: scaling it by the full
# KM_BASE/BASE factor makes the two at-risk rows collide within risk.table.height.
# 4.0 keeps the numbers legible and matched to the rest of the KM text without
# vertical overlap (risk.table.height is also raised to 0.32 below for headroom).
KM_FS   <- 4.0                    # risk-table number size (was 3.5)
KM_ANN  <- 3   * KM_BASE / BASE   # in-plot log-rank annotation size (was 3)

theme_pub <- function(bs = BASE) {
  theme_classic(base_size = bs, base_family = FONT) +
    theme(
      axis.line  = element_line(linewidth = 0.4),
      axis.ticks = element_line(linewidth = 0.3),
      axis.text  = element_text(size = bs - 1, color = "black"),
      axis.title = element_text(size = bs),
      strip.text = element_text(size = bs, face = "bold"),
      strip.background = element_blank(),
      legend.text  = element_text(size = bs - 1),
      legend.title = element_text(size = bs),
      legend.key.size = unit(3, "mm"),
      plot.title   = element_text(size = bs, face = "bold", hjust = 0),
      panel.grid   = element_blank()
    )
}

# Format a log-rank p-value: guard against floating-point underflow to 0,
# which would print an impossible "p = 0.00e+00".
fmt_p <- function(p) if (p < 2.2e-16) "Log-rank p < 2.2e-16" else sprintf("Log-rank p = %.2e", p)

# ── Load data ─────────────────────────────────────────────────────────────────
message("Loading data ...")
load(file.path(DATA_PROC, "consensus_clustering_results.RData"))
load(file.path(DATA_PROC, "clinical_characterization_results.RData"))
load(file.path(DATA_PROC, "upr_landscape_results.RData"))

clin <- heatmap_data
subtype_order <- intersect(c("UPR-high-risk","UPR-intermediate","UPR-favorable"),
                           unique(clin$UPR_subtype[!is.na(clin$UPR_subtype)]))

###############################################################################
# Panel A
###############################################################################
metrics_long <- metrics %>%
  dplyr::select(K, PAC, Silhouette_mean, Delta_AUC) %>%
  pivot_longer(-K, names_to = "Metric", values_to = "Value") %>%
  mutate(
    Metric = dplyr::recode(Metric,
      "PAC"="PAC","Silhouette_mean"="Silhouette","Delta_AUC"="Delta AUC"),
    Metric   = factor(Metric, levels = c("PAC","Silhouette","Delta AUC")),
    Selected = (K == final_k)
  )

p_A <- ggplot(metrics_long, aes(x = K, y = Value)) +
  geom_line(color = "grey50", linewidth = 0.5) +
  geom_point(shape = 21, size = 2.5, stroke = 0.5, color = "black",
             aes(fill = Selected)) +
  scale_fill_manual(values = c("FALSE"="white","TRUE"="#E64B35"), guide="none") +
  geom_point(data = dplyr::filter(metrics_long, K==final_k),
             shape = 21, size = 3.5, fill="#E64B35", color="black", stroke=0.6) +
  facet_wrap(~Metric, scales="free_y", nrow=1) +
  scale_x_continuous(breaks=2:6) +
  labs(x="Number of clusters (K)", y="Value",
       title=sprintf("K = %d selected", final_k)) +
  theme_pub()

###############################################################################
# Panel B — KM full cohort  (save as PNG, read back)
###############################################################################
message("Panel B ...")
surv_df <- clin %>%
  dplyr::filter(!is.na(OS.time), OS.time > 0, !is.na(UPR_subtype)) %>%
  mutate(OS_months = OS.time/30.44,
         UPR_subtype = factor(UPR_subtype, levels=intersect(
           c("UPR-high-risk","UPR-intermediate","UPR-favorable"),unique(UPR_subtype))))
fit_km  <- survfit(Surv(OS_months, OS) ~ UPR_subtype, data=surv_df)
lr_test <- survdiff(Surv(OS_months, OS) ~ UPR_subtype, data=surv_df)
lr_pval <- 1 - pchisq(lr_test$chisq, df=length(unique(surv_df$UPR_subtype))-1)
palette_km <- COL_SUBTYPE[levels(surv_df$UPR_subtype)]

p_km_full <- ggsurvplot(
  fit_km, data=surv_df,
  pval=FALSE, conf.int=TRUE, conf.int.alpha=0.12,
  palette=palette_km,
  xlab="Time (months)", ylab="Overall survival probability",
  legend.title="", legend.labs=levels(surv_df$UPR_subtype),
  risk.table=TRUE, risk.table.height=0.32, risk.table.y.text=FALSE,
  fontsize=KM_FS,
  ggtheme=theme_pub(KM_BASE),
  tables.theme=theme_pub(KM_BASE)+theme(axis.line=element_blank(),
                                  axis.ticks=element_blank(),
                                  axis.text=element_text(size=KM_BASE-1))
)
p_km_full$plot <- p_km_full$plot +
  annotate("text",x=Inf,y=0.98,
           label=fmt_p(lr_pval),
           hjust=1.05,vjust=1,size=KM_ANN,family=FONT) +
  labs(title="OS by UPR subtype (full cohort)") +
  theme(legend.position="bottom",
        legend.background=element_rect(fill=NA,color=NA),
        legend.spacing.x=unit(4,"mm"),
        legend.key.width=unit(7,"mm"),
        legend.text=element_text(margin=margin(r=6,unit="pt")))

km_B_png <- file.path(TMP_DIR,"km_B.png")
png(km_B_png, width=SLOT_BC_W, height=SLOT_BC_H, units="mm", res=600, family=FONT)
print(p_km_full)
dev.off()

###############################################################################
# Panel C — KM IDH-WT
###############################################################################
message("Panel C ...")
idhwt_surv <- merge(
  idhwt_cluster_df,
  clin[,c("barcode","OS.time","OS")], by="barcode"
) %>%
  dplyr::filter(!is.na(OS.time), OS.time>0, !is.na(UPR_subtype_idhwt)) %>%
  mutate(OS_months=OS.time/30.44,
         UPR_subtype=factor(UPR_subtype_idhwt,
           levels=intersect(c("UPR-high-risk","UPR-intermediate","UPR-favorable"),
                            unique(UPR_subtype_idhwt))))
n_idhwt <- nrow(idhwt_surv)
fit_km_idhwt <- survfit(Surv(OS_months,OS)~UPR_subtype, data=idhwt_surv)
lr_idhwt <- survdiff(Surv(OS_months,OS)~UPR_subtype, data=idhwt_surv)
pval_idhwt <- 1 - pchisq(lr_idhwt$chisq, df=length(unique(idhwt_surv$UPR_subtype))-1)
palette_idhwt <- COL_SUBTYPE[levels(idhwt_surv$UPR_subtype)]

p_km_idhwt <- ggsurvplot(
  fit_km_idhwt, data=idhwt_surv,
  pval=FALSE, conf.int=TRUE, conf.int.alpha=0.12,
  palette=palette_idhwt,
  xlab="Time (months)", ylab="Overall survival probability",
  legend.title="", legend.labs=levels(idhwt_surv$UPR_subtype),
  risk.table=TRUE, risk.table.height=0.32, risk.table.y.text=FALSE,
  fontsize=KM_FS,
  ggtheme=theme_pub(KM_BASE),
  tables.theme=theme_pub(KM_BASE)+theme(axis.line=element_blank(),
                                  axis.ticks=element_blank(),
                                  axis.text=element_text(size=KM_BASE-1))
)
p_km_idhwt$plot <- p_km_idhwt$plot +
  annotate("text",x=Inf,y=0.98,
           label=sprintf("%s\n(IDH-WT, n=%d)",fmt_p(pval_idhwt),n_idhwt),
           hjust=1.05,vjust=1,size=KM_ANN,family=FONT) +
  labs(title="OS by UPR subtype (IDH-WT)") +
  theme(legend.position="bottom",
        legend.background=element_rect(fill=NA,color=NA),
        legend.spacing.x=unit(4,"mm"),
        legend.key.width=unit(7,"mm"),
        legend.text=element_text(margin=margin(r=6,unit="pt")))

km_C_png <- file.path(TMP_DIR,"km_C.png")
png(km_C_png, width=SLOT_BC_W, height=SLOT_BC_H, units="mm", res=600, family=FONT)
print(p_km_idhwt)
dev.off()

###############################################################################
# Panel D — Clinical heatmap (save to PNG, read back)
###############################################################################
message("Panel D ...")
clin_ordered <- clin %>%
  dplyr::filter(!is.na(UPR_subtype)) %>%
  mutate(UPR_subtype=factor(UPR_subtype,levels=subtype_order)) %>%
  dplyr::arrange(UPR_subtype)

n_samp <- nrow(clin_ordered)
col_split <- factor(clin_ordered$UPR_subtype, levels=subtype_order)
subtype_n <- table(clin_ordered$UPR_subtype)[subtype_order]
col_title <- paste0(subtype_order,"\n(n=",subtype_n,")")

grade_col  <- c("G2"="#FFF5EB","G3"="#FC8D59","G4"="#D94701")
mgmt_col   <- c("Methylated"="#3C5488","Unmethylated"="#F39B7F","Unknown"="grey85")
gender_col <- c("male"="#4DBBD5","female"="#E64B35")

grade_vec  <- factor(clin_ordered$Grade,       levels=c("G2","G3","G4"))
idh_vec    <- factor(clin_ordered$IDH_status,  levels=c("Mutant","WT"))
mgmt_vec   <- factor(clin_ordered$MGMT_status, levels=c("Methylated","Unmethylated"))
gender_vec <- factor(clin_ordered$gender,      levels=c("male","female"))

get_p <- function(var) {
  row <- clinical_tests_df[clinical_tests_df$Variable==var,]
  if(nrow(row)==0) return(1); row$P_value[1]
}

sig_star <- function(p) {
  dplyr::case_when(p<0.0001~"****",p<0.001~"***",p<0.01~"**",p<0.05~"*",TRUE~"ns")
}

top_anno <- HeatmapAnnotation(
  `UPR subtype` = factor(clin_ordered$UPR_subtype, levels=subtype_order),
  ` `           = anno_empty(border=FALSE, height=unit(1,"mm")),
  `Grade`       = grade_vec,
  `IDH status`  = idh_vec,
  `MGMT`        = mgmt_vec,
  `Gender`      = gender_vec,
  col = list(
    `UPR subtype` = COL_SUBTYPE,
    `Grade`       = grade_col,
    `IDH status`  = COL_IDH,
    `MGMT`        = mgmt_col,
    `Gender`      = gender_col
  ),
  annotation_name_side = "left",
  annotation_name_gp   = gpar(fontsize=8, fontfamily=FONT),
  na_col="grey92",
  annotation_legend_param=list(
    `UPR subtype`=list(title_gp=gpar(fontsize=8,fontfamily=FONT),
                        labels_gp=gpar(fontsize=7,fontfamily=FONT),nrow=3),
    `Grade`=list(title_gp=gpar(fontsize=8,fontfamily=FONT),
                  labels_gp=gpar(fontsize=7,fontfamily=FONT),nrow=3),
    `IDH status`=list(title_gp=gpar(fontsize=8,fontfamily=FONT),
                       labels_gp=gpar(fontsize=7,fontfamily=FONT),nrow=2),
    `MGMT`=list(title_gp=gpar(fontsize=8,fontfamily=FONT),
                 labels_gp=gpar(fontsize=7,fontfamily=FONT),nrow=2),
    `Gender`=list(title_gp=gpar(fontsize=8,fontfamily=FONT),
                   labels_gp=gpar(fontsize=7,fontfamily=FONT),nrow=2)
  )
)

dummy_mat <- matrix(0, nrow=1, ncol=n_samp,
                    dimnames=list("",clin_ordered$barcode))
ht_clin <- Heatmap(
  dummy_mat,
  top_annotation=top_anno,
  col=c("0"="white"),
  show_heatmap_legend=FALSE,
  show_row_names=FALSE,
  show_column_names=FALSE,
  cluster_rows=FALSE, cluster_columns=FALSE,
  column_split=col_split,
  column_title=col_title,
  column_title_gp=gpar(fontsize=8,fontface="bold",fontfamily=FONT),
  border=FALSE, height=unit(0,"mm")
)

ht_D_png <- file.path(TMP_DIR,"ht_D.png")
png(ht_D_png, width=SLOT_D_W, height=SLOT_D_H, units="mm", res=300, family=FONT)
draw(ht_clin, merge_legend=TRUE,
     heatmap_legend_side="right", annotation_legend_side="right",
     align_annotation_legend="heatmap_top",
     legend_grouping="original",
     padding=unit(c(2,4,2,12),"mm"))  # bottom, left(rowname gutter), top, right(legend gutter)
dev.off()

###############################################################################
# Panel E
###############################################################################
upr_df <- data.frame(barcode=names(upr_score_bulk),UPR_score=upr_score_bulk,
                     stringsAsFactors=FALSE) %>%
  merge(clin[,c("barcode","UPR_subtype")],by="barcode") %>%
  dplyr::filter(!is.na(UPR_subtype)) %>%
  mutate(UPR_subtype=factor(UPR_subtype,levels=subtype_order))

kw_upr <- kruskal.test(UPR_score~UPR_subtype, data=upr_df)
comps <- Filter(function(p) all(p%in%unique(upr_df$UPR_subtype)), list(
  c("UPR-high-risk","UPR-intermediate"),
  c("UPR-high-risk","UPR-favorable"),
  c("UPR-intermediate","UPR-favorable")
))

p_E <- ggplot(upr_df, aes(x=UPR_subtype,y=UPR_score,fill=UPR_subtype)) +
  geom_violin(alpha=0.55,trim=TRUE,linewidth=0.3) +
  geom_boxplot(width=0.18,outlier.shape=NA,fill="white",color="black") +
  stat_compare_means(comparisons=comps,method="wilcox.test",
                     p.adjust.method="BH",label="p.signif",
                     size=2.5,tip.length=0.01,step.increase=0.06,family=FONT) +
  scale_fill_manual(values=COL_SUBTYPE,guide="none") +
  scale_x_discrete(labels=function(x) gsub("UPR-","",x)) +
  labs(x="UPR subtype",y="UPR composite score",
       title="UPR score by molecular subtype") +
  theme_pub()

###############################################################################
# Panel F
###############################################################################
upr_idh_df <- upr_df %>%
  merge(clin[,c("barcode","IDH_status")],by="barcode") %>%
  dplyr::filter(!is.na(IDH_status)) %>%
  mutate(IDH_status=factor(IDH_status,levels=c("Mutant","WT")))

p_F <- ggplot(upr_idh_df, aes(x=IDH_status,y=UPR_score,fill=IDH_status)) +
  geom_violin(alpha=0.55,trim=TRUE,linewidth=0.3) +
  geom_boxplot(width=0.18,outlier.shape=NA,fill="white",color="black") +
  # Use the same bracket + significance-star notation as Panel E (was numeric
  # "p.format" overprinting the violin); the bracket sits in the headroom above.
  stat_compare_means(comparisons=list(c("Mutant","WT")),method="wilcox.test",
                     label="p.signif",size=2.5,tip.length=0.01,
                     step.increase=0.06,family=FONT) +
  scale_fill_manual(values=COL_IDH,guide="none") +
  scale_y_continuous(expand=expansion(mult=c(0.05,0.12))) +
  labs(x="IDH status",y="UPR composite score",
       title="UPR score by IDH status") +
  theme_pub()

###############################################################################
# Assemble final composite as ggplot-only (using rasterGrob for heatmap and KM)
###############################################################################
message("=== Assembling composite ===")

# Load PNG images as rasterGrobs
img_B   <- png::readPNG(km_B_png)
img_C   <- png::readPNG(km_C_png)
img_D   <- png::readPNG(ht_D_png)

grob_B  <- grid::rasterGrob(img_B, interpolate=TRUE)
grob_C  <- grid::rasterGrob(img_C, interpolate=TRUE)
grob_D  <- grid::rasterGrob(img_D, interpolate=TRUE)

# Wrap grobs as ggplot-compatible via cowplot wrap_elements
wrap_grob <- function(g) {
  cowplot::ggdraw() + cowplot::draw_grob(g)
}

p_B_wrap <- wrap_grob(grob_B)
p_C_wrap <- wrap_grob(grob_C)
p_D_wrap <- wrap_grob(grob_D)

# Build composite: 4 row layout
# Row 1: A (full width)
# Row 2: B | C  (KM side by side)
# Row 3: D (full width heatmap)
# Row 4: E | F  (violins)
row1 <- p_A
row2 <- p_B_wrap | p_C_wrap
row3 <- p_D_wrap
row4 <- p_E | p_F

p_composite <- (row1 / row2 / row3 / row4) +
  plot_annotation(
    tag_levels = "A",
    theme = theme(
      plot.tag = element_text(size=13, face="bold", family=FONT,
                              color="black", margin=margin(0,2,2,0,"mm"))
    )
  ) +
  plot_layout(heights = LAYOUT_HEIGHTS)

out_pdf <- file.path(OUT_DIR,"Figure2_composite.pdf")
out_png <- file.path(OUT_DIR,"Figure2_composite.png")

W_IN <- W_MM * MM_IN
H_IN <- H_MM * MM_IN

ggsave(out_pdf, plot=p_composite, device=cairo_pdf,
       width=W_IN, height=H_IN, units="in", dpi=300)
message("Saved PDF: ", out_pdf)

ggsave(out_png, plot=p_composite,
       width=W_IN, height=H_IN, units="in", dpi=300)
message("Saved PNG: ", out_png)

# Cleanup temp files
unlink(TMP_DIR, recursive=TRUE)

message("=== Figure 2 composite done ===")
message(sprintf("Dimensions: %d x 260 mm  (W/H aspect = %.2f)", W_MM, W_MM/260))
