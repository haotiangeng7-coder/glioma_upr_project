# Requirements

## R version

R v4.3.x

## Global random seed

`set.seed(42)` — used consistently throughout the analysis pipeline.

## Key R packages and tested versions

| Package | Version | Purpose |
|---|---|---|
| TCGAbiolinks | 2.28.4 | TCGA data download |
| Seurat | 5.0 | Single-cell analysis |
| Harmony | 1.0 | Batch correction |
| SingleR | 2.2 | Cell type annotation |
| SCTransform | 0.4 | Normalization (vst.flavor="v2") |
| AUCell | 1.18 | Gene set scoring |
| UCell | 2.4 | Gene set scoring |
| CellChat | 2.0 | Ligand-receptor analysis |
| scDblFinder | — | Doublet detection |
| DESeq2 | 1.38 | Differential expression |
| ConsensusClusterPlus | 1.60 | Consensus clustering |
| GSVA | 1.44 | Gene set variation analysis |
| MCPcounter | 1.2 | Immune cell deconvolution |
| survival | 3.5 | Survival analysis |
| survminer | 0.4.9 | Survival visualization |
| timeROC | 0.4 | Time-dependent ROC |
| survRM2 | 1.0 | Restricted mean survival |
| glmnet | 4.1 | LASSO / elastic net |
| randomForestSRC | 3.2 | Random survival forest |
| CoxBoost | 1.4 | Cox model boosting |
| rms | 6.7 | Nomogram / calibration |
| clusterProfiler | 4.8 | Enrichment analysis |
| fgsea | 1.22 | Fast GSEA |
| maftools | 2.16 | Mutation analysis |
| oncoPredict | 0.2 | Drug sensitivity |
| ComplexHeatmap | 2.14 | Heatmap visualization |
| ggplot2 | 3.4 | Base visualization |

## Environment variables

`06_clinical_translation/03_cmap_analysis.R` requires the environment variable
`CLUE_API_KEY` to be set (the script reads it via `Sys.getenv("CLUE_API_KEY")`).
No key is included in this repository. Obtain one from https://clue.io/.

## Python packages (for 07_revision_analyses/)

- tidepy
- pdfplumber
