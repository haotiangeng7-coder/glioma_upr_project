# A 15-Gene Unfolded Protein Response Signature for Prognostic Risk Stratification and Immunotherapy Guidance in Glioma

This repository contains the complete analysis code for the manuscript
"A 15-Gene Unfolded Protein Response Signature for Prognostic Risk
Stratification and Immunotherapy Guidance in Glioma." The analysis
encompasses: single-cell characterization of UPR activity across the
glioma tumor microenvironment; bulk transcriptomic consensus clustering
to define UPR-driven molecular subtypes; a two-stage LASSO pipeline with
110-combination benchmarking and 5x5 nested cross-validation to build a
15-gene UIRS prognostic model; and clinical translation analyses
including immunotherapy response prediction, drug sensitivity, and CMap
connectivity mapping.

## Repository structure

```
.
├── 00_setup/                  Environment config, random seeds, package list
├── 01_data_download/          TCGA, CGGA, GEO, and scRNA-seq data download
├── 02_gene_sets/              UPR gene set definitions (three arms + broad)
├── 03_single_cell/            Part 1: scRNA-seq preprocessing, annotation,
│                              UPR scoring, differential analysis, CellChat
├── 04_bulk_subtyping/         Part 2: UPR landscape, consensus clustering,
│                              clinical/immune/genomic/pathway characterization
├── 05_ml_model/               Part 3: feature selection, 110-combination
│                              benchmark, risk score, independent prognosis,
│                              nomogram, permutation tests
├── 06_clinical_translation/   Part 4: immunotherapy prediction, drug
│                              sensitivity, CMap analysis, hub gene analysis
├── 07_revision_analyses/      Revision analyses: GSE182109 scRNA validation,
│                              published signature comparison, nomogram
│                              calibration, pseudobulk DESeq2, IDHwt
│                              improvement, GSE16011 validation, permutation
│                              tests, TIDE prediction
├── figures/                   Figure generation: main figures (Fig1-8),
│                              composite assembly, supplementary figures
├── run_all.R                  Master pipeline script
├── requirements.md            Package versions and dependencies
├── LICENSE                    MIT License
├── .gitignore
└── PATH_NORMALIZATION_LOG.md  Record of path de-hardcoding edits
```

## Data availability

All datasets analyzed are publicly available. Bulk transcriptomic data:
TCGA-GBM/LGG via the GDC portal (https://portal.gdc.cancer.gov/), CGGA
(http://www.cgga.org.cn/), and GEO accessions GSE16011 and GSE43378.
Single-cell/single-nucleus RNA-seq: GEO GSE131928 and GSE182109. Raw
data are not included in this repository; download via 01_data_download/
before running.

## Requirements

See [requirements.md](requirements.md) for the full list of R and Python
package dependencies with tested versions.

## How to run

1. Download the public datasets using scripts in `01_data_download/`.
2. From the repository root, run:
   ```
   Rscript run_all.R
   ```
   The global seed is set to 42 in `00_setup/config.R` to ensure
   reproducibility.

## License

MIT. See [LICENSE](LICENSE).

## Citation

Geng H, et al. A 15-Gene Unfolded Protein Response Signature Enables Prognostic Stratification Across Glioma Molecular Subtypes and Reveals an Immune-Excluded Tumor Microenvironment. Journal of Neuro-Oncology. 2026 (under review). DOI: TBD
