# Manuscript code and demo data

This repository contains analysis code and representative demo data associated with the manuscript.

Currently included module:

- Bulk RNA-seq analysis for Fig. 4, Supplementary Fig. 4 and Supplementary Fig. 5

Additional analysis modules will be added later.

---

# Repository structure

```text
manuscript-code/
  README.md

  code/
    bulkRNAseq/
      Fig4_SFig4_SFig5/
        original/
          01_bulkRNAseq_Fig4_SFig4_SFig5_full_original_analysis_github.R

        demo/
          01_bulkRNAseq_Fig4_SFig4_SFig5_full_demo_analysis.R

  demo_data/
    bulkRNAseq_Fig4_SFig4_SFig5_demo/
      gene_counts_demo.csv
      metadata_demo.csv
      gene_annotation_demo.csv

  results/
    bulkRNAseq_Fig4_SFig4_SFig5_demo/
      Generated after running the demo script
```

---

# Bulk RNA-seq analysis for Fig. 4, Supplementary Fig. 4 and Supplementary Fig. 5

## Overview

This module contains the bulk RNA-seq downstream analysis code used for Fig. 4, Supplementary Fig. 4 and Supplementary Fig. 5.

The repository provides two versions of the analysis script:

```text
code/bulkRNAseq/Fig4_SFig4_SFig5/original/
  01_bulkRNAseq_Fig4_SFig4_SFig5_full_original_analysis_github.R
```

This script is intended for the full original analysis using the complete featureCounts gene-level count matrix, sample metadata and GENCODE GTF annotation file.

```text
code/bulkRNAseq/Fig4_SFig4_SFig5/demo/
  01_bulkRNAseq_Fig4_SFig4_SFig5_full_demo_analysis.R
```

This script uses the demo dataset provided in this repository and is intended to test whether the workflow can be executed successfully.

The demo dataset is a representative subset of the processed bulk RNA-seq data. It is not expected to reproduce the exact quantitative results reported in the manuscript.

---

## Demo data

The demo dataset is located in:

```text
demo_data/bulkRNAseq_Fig4_SFig4_SFig5_demo/
```

It contains the following files:

```text
gene_counts_demo.csv
metadata_demo.csv
gene_annotation_demo.csv
```

The demo count matrix contains 5,000 genes and 24 samples across six groups:

```text
DMSO_un
GSK_un
PBMC_un
DMSO_CD19
GSK_CD19
PBMC_CD19
```

Each group contains four biological replicates.

The demo data are provided only for testing the workflow and checking the expected input/output structure.

Because the demo dataset is a subset of the complete dataset, some differential expression results, KEGG/GO enrichment results and selected pathway plots may differ from the full manuscript analysis.

---

## Analysis included in the demo script

The demo script performs the following analyses:

1. Read demo count matrix, metadata and gene annotation
2. Generate basic QC summaries and mean ± SEM plots
3. Construct a DESeq2 object and perform VST normalization
4. Generate six-group PCA plot
5. Generate Pearson correlation heatmap based on the top 1000 highly variable genes
6. Perform pairwise DESeq2 differential expression analysis
7. Perform three-group overall DESeq2 likelihood ratio test analysis
8. Generate pairwise volcano plots
9. Generate three-group top1000 DEG heatmaps with selected pathway-symbol gene labels
10. Perform KEGG and GO enrichment analysis
11. Generate selected KEGG pathway bar plots when the selected pathways are present in the demo enrichment results
12. Save R session information

Because the demo dataset is a representative subset of the full dataset, some selected KEGG pathway plots may not be generated if the corresponding pathways are not enriched in the demo results. This is expected and does not indicate an error in the workflow.

---

## How to run the demo

Open R or RStudio and set the working directory to the repository root:

```r
setwd("path/to/manuscript-code")
```

Then run:

```r
source("code/bulkRNAseq/Fig4_SFig4_SFig5/demo/01_bulkRNAseq_Fig4_SFig4_SFig5_full_demo_analysis.R")
```

The demo results will be written to:

```text
results/bulkRNAseq_Fig4_SFig4_SFig5_demo/
```

Expected output folders include:

```text
00_basic_data/
00_QC/
01_PCA/
02_Pearson_top1000/
03_DESeq2_DEG/
04_pairwise_Volcano/
05_Three_group_DEG_heatmap/
06_KEGG_GO/
07_Selected_pathway/
sessionInfo.txt
```

---

## Expected runtime

The demo analysis was tested on a Windows 11 laptop with the following specifications:

```text
CPU: Intel Core i5-10300H @ 2.50 GHz
RAM: 64 GB
Operating system: Windows 11, 64-bit
```

The demo analysis takes approximately:

```text
15 minutes
```

Runtime may vary depending on the computer, internet connection and whether the required R/Bioconductor packages have already been installed.

The KEGG/GO enrichment steps may take longer on the first run because they require access to online annotation resources.

---

## Software requirements

The workflow was developed and tested in R.

Tested environment:

```text
R version 4.5.1
Platform: x86_64-w64-mingw32/x64
Operating system: Windows 11 x64
```

Main R packages directly used in the scripts:

```text
DESeq2 v1.48.2
ggplot2 v4.0.2
dplyr v1.1.4
tidyr v1.3.2
readr v2.1.6
pheatmap v1.0.13
ComplexHeatmap v2.24.1
circlize v0.4.17
clusterProfiler v4.16.0
org.Hs.eg.db v3.21.0
data.table v1.18.0
stringr v1.6.0
RColorBrewer v1.1-3
```

The script will attempt to install missing CRAN and Bioconductor packages automatically.

After running the script, complete R session information, including package dependencies, is saved to:

```text
results/bulkRNAseq_Fig4_SFig4_SFig5_demo/sessionInfo.txt
```

---

## Original full analysis script

The original full analysis script is provided at:

```text
code/bulkRNAseq/Fig4_SFig4_SFig5/original/
  01_bulkRNAseq_Fig4_SFig4_SFig5_full_original_analysis_github.R
```

This script is intended for reproducing the full bulk RNA-seq analysis using the complete dataset.

The original analysis script requires the following full input files:

```text
featureCounts gene-level count matrix
sample metadata file
GENCODE GTF annotation file
```

These full input files are not included in this repository.

Before running the original script, users should modify the following paths according to their local environment:

```r
raw_count_file <- "path/to/gene_counts.txt"
metadata_file  <- "path/to/metadata.csv"
gtf_file       <- "path/to/gencode.annotation.gtf.gz"

result_dir <- "path/to/output/bulkRNA_Fig4_SFig4_SFig5"
```

The original script performs the same analysis workflow as the demo script but uses the complete bulk RNA-seq dataset.

---

## Notes on reproducibility

The demo analysis is designed to verify that the code can be executed and that the expected output structure is generated.

Because the demo dataset is a subset of the full dataset:

- PCA and correlation patterns may be similar but not identical to the manuscript figures.
- DEG numbers may differ from the full analysis.
- Some KEGG/GO pathways may not be enriched.
- Some selected pathway plots may not be generated.
- Three-group DEG heatmaps and selected pathway plots may contain fewer labeled features than the full manuscript analysis.

This behavior is expected and does not indicate an error in the workflow.

---

# Planned additional modules

The following sections are placeholders for additional analysis modules that will be added later.

## Bulk RNA-seq analysis for Fig. 2

To be added.

Planned contents:

```text
code/bulkRNAseq/Fig2/original/
code/bulkRNAseq/Fig2/demo/
demo_data/bulkRNAseq_Fig2_demo/
```

---

## ATAC-seq analysis for Fig. 4, Supplementary Fig. 4 and Supplementary Fig. 5

To be added.

Planned contents:

```text
code/ATACseq/Fig4_SFig4_SFig5/original/
code/ATACseq/Fig4_SFig4_SFig5/demo/
demo_data/ATACseq_Fig4_SFig4_SFig5_demo/
```

---

## ATAC-seq analysis for Fig. 2

To be added.

Planned contents:

```text
code/ATACseq/Fig2/original/
code/ATACseq/Fig2/demo/
demo_data/ATACseq_Fig2_demo/
```

---

## scRNA-seq analysis

To be added.

Planned contents:

```text
code/scRNAseq/original/
code/scRNAseq/demo/
demo_data/scRNAseq_demo/
```

---

## Additional analysis scripts

To be added.

This section may include additional plotting scripts, integrative RNA-seq/ATAC-seq analysis scripts, gene set scoring scripts or other supporting analyses associated with the manuscript.

---

# License

This repository is currently provided for peer review and reproducibility purposes. A formal license will be added before public release.



