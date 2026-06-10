#!/usr/bin/env Rscript
# Recompute S3C: pairwise log-rank tests in IDH-WT subgroup for K=3
# Real-data computation only. No hardcoded p-values.
suppressPackageStartupMessages({
  library(survival); library(survminer); library(ggplot2)
})
PROJ <- getwd()
OUT  <- file.path(PROJ, "figures", "supp_assembled")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

load(file.path(PROJ, "data/processed/consensus_clustering_results.RData"))

# K=3 IDH-WT consensus labels (raw cluster numbers)
k3 <- cc_results_idhwt[[3]]$consensusClass
stopifnot(length(k3) == 394)

df <- data.frame(barcode = names(k3), cl = as.integer(k3), stringsAsFactors = FALSE)
sf <- surv_final[, c("barcode", "OS.time", "OS")]
df <- merge(df, sf, by = "barcode")
df <- df[!is.na(df$OS.time) & df$OS.time > 0 & !is.na(df$OS), ]
df$time_months <- df$OS.time / 30
cat(sprintf("IDH-WT K=3 samples with survival: %d\n", nrow(df)))

# Name clusters by median OS (longest = favorable, shortest = high-risk)
med_os <- tapply(df$OS.time, df$cl, median, na.rm = TRUE)
ord <- order(med_os, decreasing = TRUE)            # cluster index by descending median OS
nm  <- setNames(c("UPR-favorable", "UPR-intermediate", "UPR-high-risk"),
                as.character(ord))
df$subtype <- factor(nm[as.character(df$cl)],
                     levels = c("UPR-favorable", "UPR-intermediate", "UPR-high-risk"))
cat("Cluster median OS (days):\n"); print(med_os)
cat("Subtype composition:\n"); print(table(df$subtype))

# Pairwise log-rank (BH not applied; report raw pairwise p as in caption)
pw <- pairwise_survdiff(Surv(time_months, OS) ~ subtype, data = df,
                        p.adjust.method = "none")
cat("\n=== Pairwise log-rank p-value matrix (raw) ===\n")
print(pw$p.value)

pm <- pw$p.value
# Extract the two values quoted in caption
get_p <- function(a, b) {
  v <- NA
  if (a %in% rownames(pm) && b %in% colnames(pm)) v <- pm[a, b]
  if (is.na(v) && b %in% rownames(pm) && a %in% colnames(pm)) v <- pm[b, a]
  v
}
p_int_fav  <- get_p("UPR-intermediate", "UPR-favorable")
p_int_high <- get_p("UPR-high-risk",    "UPR-intermediate")
cat(sprintf("\nintermediate vs favorable : p = %.3f  (caption: 0.186)\n", p_int_fav))
cat(sprintf("intermediate vs high-risk : p = %.3f  (caption: 0.108)\n", p_int_high))

# All three pairwise for completeness
p_fav_high <- get_p("UPR-high-risk", "UPR-favorable")
cat(sprintf("favorable vs high-risk    : p = %.3g\n", p_fav_high))

# --- Render S3C panel: pairwise p-value heatmap/table ---
labs <- c("UPR-favorable","UPR-intermediate","UPR-high-risk")
short <- c("Favorable","Intermediate","High-risk")
mat <- matrix(NA_real_, 3, 3, dimnames = list(labs, labs))
for (i in 1:3) for (j in 1:3) if (i != j) mat[i,j] <- get_p(labs[i], labs[j])

celldf <- expand.grid(row = factor(short, levels = rev(short)),
                      col = factor(short, levels = short))
celldf$p <- mapply(function(r, c) {
  ri <- short == as.character(r); ci <- short == as.character(c)
  if (identical(which(ri), which(ci))) return(NA_real_)
  get_p(labs[which(ri)], labs[which(ci)])
}, celldf$row, celldf$col)
celldf$lab <- ifelse(is.na(celldf$p), "", sprintf("p = %.3f", celldf$p))

p_panel <- ggplot(celldf, aes(col, row, fill = p)) +
  geom_tile(color = "white", linewidth = 1.2) +
  geom_text(aes(label = lab), size = 4) +
  scale_fill_gradient(low = "#E64B35", high = "#4DBBD5", na.value = "grey92",
                      limits = c(0, 1), name = "log-rank p") +
  labs(x = NULL, y = NULL,
       title = "IDH-WT K=3 pairwise log-rank") +
  coord_equal() +
  theme_minimal(base_size = 12) +
  theme(panel.grid = element_blank(),
        axis.text = element_text(color = "black"),
        plot.title = element_text(face = "bold", hjust = 0.5, size = 13),
        legend.position = "right")

ggsave(file.path(OUT, "_panel_S3C_pairwise_logrank.pdf"), p_panel,
       width = 5.2, height = 4.0)
cat("\nSaved _panel_S3C_pairwise_logrank.pdf\n")
saveRDS(list(p_int_fav=p_int_fav, p_int_high=p_int_high, p_fav_high=p_fav_high,
             n=nrow(df), comp=table(df$subtype)),
        file.path(OUT, "_S3C_values.rds"))
