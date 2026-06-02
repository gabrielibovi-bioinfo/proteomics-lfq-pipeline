# 🔬 proteomics-lfq-pipeline

> Reproducible R workflow for label-free quantification (LFQ) proteomics analysis — from MaxQuant output to differential expression and pathway enrichment.

![Platform](https://img.shields.io/badge/platform-Linux%20%7C%20macOS%20%7C%20Windows-blue)
![R](https://img.shields.io/badge/R-%3E%3D4.1-276DC3?logo=r)
![Status](https://img.shields.io/badge/status-in%20development-yellow)

---

## Overview

This pipeline processes LFQ proteomics data from MaxQuant output through quality control, normalization, differential expression analysis, and gene set enrichment analysis (GSEA) using GO, KEGG, and Reactome databases.

```
MaxQuant output (.xlsx)
    │
    ├── [Step 1] Preprocessing             → log2 transform, zero/NA/constant removal
    │
    ├── [Step 2] Quality Control           → heatmaps, Pearson correlation, PCA
    │
    ├── [Step 3] Outlier Removal           → visual QC, re-run after exclusion
    │
    ├── [Step 4] Normalization             → proteoDA: RLR normalization
    │                                        normalization and QC reports
    │
    ├── [Step 5] Differential Expression   → limma: MO vs SC contrast
    │                                        volcano plot, DEP heatmap, PCA
    │
    └── [Step 6] Pathway Enrichment        → GSEA: GO, KEGG, Reactome
                                             dotplot, ridgeplot, heatplot,
                                             cnetplot, enrichment map
```

---

## Repository Structure

```
proteomics-lfq-pipeline/
├── 01_qc_normalization.R    # Preprocessing, QC, normalization, limma
├── 02_dep_analysis.R        # DEP visualization: volcano, heatmap, PCA
├── 03_gsea.R                # GSEA: GO, KEGG, Reactome
├── data/                    # Place input .xlsx file here
├── results/                 # Output files
├── README.md
├── .gitignore
└── LICENSE
```

---

## Requirements

- R ≥ 4.1

### R dependencies

```r
install.packages(c("readxl", "ggplot2", "tidyr", "dplyr", "pheatmap", "factoextra", "plotly", "ggfortify", "RColorBrewer",
                   "readr", "tibble", "ggrepel", "writexl"))

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("preprocessCore", "affy", "vsn", "ComplexHeatmap", "clusterProfiler", "enrichplot",
                        "gprofiler2", "AnnotationDbi", "org.Rn.eg.db", "ReactomePA"))

# proteoDA (normalization and differential expression)
devtools::install_github("ByrumLab/proteoDA", dependencies = TRUE,
                         build_vignettes = TRUE, force = TRUE)
```

---

## Usage

### Step 1 — QC, normalization, and differential expression

Update the file path at the top of the script:

```r
FILE_PATH <- "data/proteomics_data.xlsx"
```

**Outputs:**
```
Normalization_report/               # PDF reports (normalization + QC)
limma_rlr_BH/per_contrast_results/  # MO_vs_SC.csv
figuras/pre-processing/             # Heatmaps, PCA (pre/post outlier removal)
```

### Step 2 — DEP visualization

**Outputs:**
```
analise_DEP/logFC_percentil.xlsx
analise_DEP/genes_lFC.xlsx
analise_DEP/contagem_degs_sig.xlsx
analise_DEP/01_volcano_plot.tiff
analise_DEP/02_heatmap_DEP.tiff
analise_DEP/03_pca_pvalue.tiff
```

### Step 3 — GSEA (GO, KEGG, Reactome)

**Outputs:**
```
gsea_GO/        # GO results (CSV) and 13 figures
gsea_kegg/      # KEGG results (CSV) and 13 figures
gsea_reactome/  # Reactome results (CSV) and 13 figures
```

---

## Analysis Parameters

| Parameter           | Value  | Step                       |
|---------------------|--------|----------------------------|
| Min valid values    | ≥ 1    | Preprocessing              |
| Min group proportion | 50%  | proteoDA filter            |
| Normalization method | RLR  | proteoDA normalize_data    |
| p-value threshold   | 0.05   | limma + GSEA               |
| logFC threshold     | 0      | DEP categorization         |
| adj. method (limma) | BH     | extract_DA_results         |
| GSEA ranking metric | t-statistic | bitr + gene_list      |
| Min gene set size   | 5      | GSEA (all databases)       |
| Max gene set size   | 800    | GSEA (all databases)       |

Modify based on your data.

---

## Tools and References

| Tool            | Reference |
|-----------------|-----------|
| proteoDA        | [Byrum Lab, GitHub](https://github.com/ByrumLab/proteoDA) |
| limma           | [Ritchie et al., 2015](https://doi.org/10.1093/nar/gkv007) |
| clusterProfiler | [Yu et al., 2012](https://doi.org/10.1089/omi.2011.01) |
| ReactomePA      | [Yu & He, 2016](https://doi.org/10.1039/C5MB00663E) |
| org.Rn.eg.db    | [Bioconductor](https://bioconductor.org/packages/org.Rn.eg.db) |
| pheatmap        | [Kolde, 2019](https://CRAN.R-project.org/package=pheatmap) |
| ggplot2         | [Wickham, 2016](https://doi.org/10.1007/978-3-319-24277-4) |

---

## Notes

- All scripts use relative paths from the project root. Run from the `proteomics-lfq-pipeline/` directory.
- Update `FILE_PATH` in `01_qc_normalization.R` to point to your MaxQuant Excel output.
- Outlier sample `SC-1` was removed based on visual QC (correlation heatmap and PCA). Update the sample name if applying to a different dataset.
- GSEA uses the **t-statistic** from limma as the ranking metric, which captures both effect size and significance. This is preferred over logFC alone.
- The organism database used is `org.Rn.eg.db` (rat). For other organisms, replace with the appropriate `org.Xx.eg.db` package and update the `organism` argument in `gseKEGG` and `gsePathway`.

---

## Author

**Gabrieli Bovi**
Bioinformatics | Proteomics | Multi-omics Data Analysis
🔗 [github.com/gabrielibovi-bioinfo](https://github.com/gabrielibovi-bioinfo)

---

## License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.
