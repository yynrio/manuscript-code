# Manuscript code and demo data

This repository contains analysis code and representative demo data associated with the manuscript.

Currently included modules:

- Bulk RNA-seq analysis for Fig. 4, Supplementary Fig. 4 and Supplementary Fig. 5
- Bulk ATAC-seq analysis for Fig. 4, Supplementary Fig. 4 and Supplementary Fig. 5
- HSPC reference comparison analysis for Extended Data Fig. 1a,b
- TCR-seq analysis for Extended Data Fig. 3

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

    ATACseq/
      Fig4_SFig4_SFig5/
        original/
          ATACseq_Fig4_SFig4_SFig5_original_analysis_GitHub.R

        demo/
          ATACseq_Fig4_SFig4_SFig5_demo_analysis.R

    HSPC_reference_comparison/
      HSPC_reference_comparison_SFig1_a_b/
        original/
          HSPC_reference_comparison_ExtendedData_Fig1ab_full_original_analysis_github.R

        demo/
          HSPC_reference_comparison_SFig1_a_b_demo_analysis.R

    TCRseq/
      ExtendedData_Fig3/
        original/
          TCRseq_ExtendedData_Fig3_original_analysis_GitHub.py
          TCRseq_ExtendedData_Fig3_original_analysis_GitHub.ipynb

        demo/
          TCRseq_ExtendedData_Fig3_demo_analysis.py
          TCRseq_ExtendedData_Fig3_demo_analysis.ipynb
          requirements.txt

  demo_data/
    bulkRNAseq_Fig4_SFig4_SFig5_demo/
      gene_counts_demo.csv
      metadata_demo.csv
      gene_annotation_demo.csv

    ATACseq_Fig4_SFig4_SFig5_demo/
      peak_counts_samples_demo.txt
      metadata_demo.csv
      README_demo_data.txt

      IGVgene_demo/
        graph_Core_genes_list_with_transcript_id_demo.xlsx
        ncbiRefSeqSelect_demo.txt.gz

        bigwig_demo/
          DMSO-1_demo.bw
          GSK-1_demo.bw
          PBMC-1_demo.bw
          DMSO-CD19-1_demo.bw
          GSK-CD19-1_demo.bw
          PBMC-CD19-1_demo.bw

    HSPC_reference_comparison_SFig1_a_b_demo/
      HSPC_waves.txt
      GSE253355_BM_HSPC_demo.rds
      GSE162950_fetal_liver_demo.rds
      Tera_integrate_Label_demo.rds

    TCRseq_ExtendedData_Fig3_demo/
      TCR_Clones_Collated_demo.xlsx
      demo_data_summary.csv

  results/
    bulkRNAseq_Fig4_SFig4_SFig5_demo/
      Generated after running the bulk RNA-seq demo script

    ATACseq_Fig4_SFig4_SFig5_demo/
      Generated after running the ATAC-seq demo script

    HSPC_reference_comparison_SFig1_a_b_demo/
      Generated after running the HSPC reference comparison demo script

    TCRseq_ExtendedData_Fig3_demo/
      Generated after running the TCR-seq demo script
```

---

# System requirements and installation

## Hardware requirements

No non-standard hardware is required. The demo workflows were tested on a standard Windows 11 laptop.

Tested hardware:

```text
CPU: Intel Core i5-10300H @ 2.50 GHz
RAM: 64 GB
Operating system: Windows 11, 64-bit
```

## Installation time

The R scripts attempt to install missing CRAN and Bioconductor packages automatically. The TCR-seq Python demo uses packages listed in its `requirements.txt` file and can be installed with `pip`.

Typical installation time on a standard desktop computer:

```text
20-60 minutes
```

Installation time may vary depending on the computer, internet connection and whether the required R/Bioconductor or Python packages have already been installed.

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
2. Generate basic QC summaries and mean +/- SEM plots
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

## How to run the bulk RNA-seq demo

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

# Bulk ATAC-seq analysis for Fig. 4, Supplementary Fig. 4 and Supplementary Fig. 5

## Overview

This module contains the bulk ATAC-seq downstream analysis code used for Fig. 4, Supplementary Fig. 4 and Supplementary Fig. 5.

The repository provides two versions of the analysis script:

```text
code/ATACseq/Fig4_SFig4_SFig5/original/
  ATACseq_Fig4_SFig4_SFig5_original_analysis_GitHub.R
```

This script is intended for the full original analysis using the complete peak-level count matrix, sample metadata, gene-coordinate table, RefSeq Select annotation file and bigWig files.

```text
code/ATACseq/Fig4_SFig4_SFig5/demo/
  ATACseq_Fig4_SFig4_SFig5_demo_analysis.R
```

This script uses the demo dataset provided in this repository and is intended to test whether the ATAC-seq workflow can be executed successfully.

The demo dataset is a representative subset of the processed ATAC-seq data. It is not expected to reproduce the exact quantitative results reported in the manuscript.

---

## Demo data

The ATAC-seq demo dataset is located in:

```text
demo_data/ATACseq_Fig4_SFig4_SFig5_demo/
```

It contains the following files:

```text
peak_counts_samples_demo.txt
metadata_demo.csv
README_demo_data.txt

IGVgene_demo/
  graph_Core_genes_list_with_transcript_id_demo.xlsx
  ncbiRefSeqSelect_demo.txt.gz

  bigwig_demo/
    DMSO-1_demo.bw
    GSK-1_demo.bw
    PBMC-1_demo.bw
    DMSO-CD19-1_demo.bw
    GSK-CD19-1_demo.bw
    PBMC-CD19-1_demo.bw
```

The demo peak count matrix contains all samples and the top 5,000 highly variable filtered peaks across six groups:

```text
DMSO_un
GSK_un
PBMC_un
DMSO_CD19
GSK_CD19
PBMC_CD19
```

The demo data also include small regional bigWig files for IGV-like ATAC signal plotting. These bigWig files contain only selected genomic regions and are intended for testing the plotting workflow.

The demo data are provided only for testing the workflow and checking the expected input/output structure.

Because the demo dataset is a subset of the complete dataset, differential accessibility results, KEGG/GO enrichment results and selected pathway plots may differ from the full manuscript analysis.

---

## Analysis included in the demo script

The ATAC-seq demo script performs the following analyses:

1. Read demo peak-level count matrix and sample metadata
2. Clean count matrix and match metadata
3. Generate total-count QC plot
4. Filter low-count peaks
5. Generate six-group PCA plot
6. Generate Pearson correlation heatmap based on the top 1000 highly variable peaks and group-average signals
7. Annotate filtered peaks using ChIPseeker
8. Construct a DESeq2 object
9. Perform pairwise DESeq2 differential accessibility analysis
10. Generate pairwise DAR volcano plots
11. Summarize significant DAR annotation distributions
12. Perform KEGG and GO enrichment analysis
13. Generate selected GO pathway bar plots when the selected GO terms are present in the demo enrichment results
14. Generate IGV-like ATAC signal plots from small regional demo bigWig files
15. Save R session information

Because the demo dataset is a representative subset of the full dataset, some comparisons may contain fewer significant DARs, and some KEGG/GO enrichment results or selected GO pathway plots may not be generated. This is expected and does not indicate an error in the workflow.

---

## How to run the ATAC-seq demo

Open R or RStudio and set the working directory to the repository root:

```r
setwd("path/to/manuscript-code")
```

Then run:

```r
source("code/ATACseq/Fig4_SFig4_SFig5/demo/ATACseq_Fig4_SFig4_SFig5_demo_analysis.R")
```

The demo results will be written to:

```text
results/ATACseq_Fig4_SFig4_SFig5_demo/
```

Expected output folders include:

```text
00_basic_data/
00_QC/
01_PCA/
02_Pearson_top1000/
03_peak_annotation/
04_DESeq2_DAR/
05_pairwise_Volcano/
06_DAR_annotation_summary/
07_KEGG_GO/
08_Selected_pathway/
09_IGV_like_bigWig_plot/
sessionInfo.txt
```

---

## Expected runtime

The ATAC-seq demo analysis was tested on a Windows 11 laptop with the following specifications:

```text
CPU: Intel Core i5-10300H @ 2.50 GHz
RAM: 64 GB
Operating system: Windows 11, 64-bit
```

The demo analysis takes approximately:

```text
10-20 minutes
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

Main R packages directly used in the ATAC-seq scripts:

```text
DESeq2 v1.48.2
ChIPseeker v1.44.0
clusterProfiler v4.16.0
org.Hs.eg.db v3.21.0
TxDb.Hsapiens.UCSC.hg38.knownGene v3.21.0
GenomicRanges v1.60.0
GenomeInfoDb v1.44.3
Gviz v1.52.0
rtracklayer v1.68.0
ggplot2 v4.0.2
ggrepel v0.9.6
pheatmap v1.0.13
data.table v1.18.0
dplyr v1.1.4
tidyr v1.3.2
readxl v1.4.5
```

The script will attempt to install missing CRAN and Bioconductor packages automatically.

After running the script, complete R session information, including package dependencies, is saved to:

```text
results/ATACseq_Fig4_SFig4_SFig5_demo/sessionInfo.txt
```

---

## Original full analysis script

The original full ATAC-seq analysis script is provided at:

```text
code/ATACseq/Fig4_SFig4_SFig5/original/
  ATACseq_Fig4_SFig4_SFig5_original_analysis_GitHub.R
```

This script is intended for reproducing the full ATAC-seq analysis using the complete dataset.

The original analysis script requires the following full input files:

```text
peak-level count matrix
sample metadata file
bigWig files
gene-coordinate table for IGV-like plots
RefSeq Select annotation file
```

These full input files are not included in this repository.

Before running the original script, users should modify the following paths according to their local environment:

```r
count_file    <- "path/to/peak_counts_samples.txt"
metadata_file <- "path/to/metadata.csv"

bw_dir <- "path/to/bigwig"

bw_files <- c(
  "DMSO-1"       = file.path(bw_dir, "DMSO-1.bw"),
  "GSK-1"        = file.path(bw_dir, "GSK-1.bw"),
  "PBMC-1"       = file.path(bw_dir, "PBMC-1.bw"),
  "DMSO-CD19-1"  = file.path(bw_dir, "DMSO-CD19-1.bw"),
  "GSK-CD19-1"   = file.path(bw_dir, "GSK-CD19-1.bw"),
  "PBMC-CD19-1"  = file.path(bw_dir, "PBMC-CD19-1.bw")
)

igv_gene_file <- "path/to/graph_Core_genes_list_with_transcript_id.xlsx"
refseq_file   <- "path/to/ncbiRefSeqSelect.txt.gz"

result_dir <- "path/to/output/ATACseq_Fig4_SFig4_SFig5"
```

The original script performs the same analysis workflow as the demo script but uses the complete ATAC-seq dataset.

---

## Notes on reproducibility

The ATAC-seq demo analysis is designed to verify that the code can be executed and that the expected output structure is generated.

Because the demo dataset is a subset of the full dataset:

- PCA and correlation patterns may be similar but not identical to the manuscript figures.
- DAR numbers may differ from the full analysis.
- Genomic annotation distributions may differ from the full analysis.
- Some KEGG/GO pathways may not be enriched.
- Some selected GO pathway plots may not be generated.
- IGV-like signal plots are generated from small regional demo bigWig files and are intended only for testing the plotting workflow.

This behavior is expected and does not indicate an error in the workflow.

---

# HSPC reference comparison analysis for Extended Data Fig. 1a,b

## Overview

This module contains the HSPC reference comparison analysis code used for Extended Data Fig. 1a,b.

Our previous study established that transient inactivation of SP140 enables the generation of hematopoietic stem and progenitor cells (HSPCs) with robust serial transplantation capacity from hPSC-derived teratomas:

```text
Liu, X. et al. Transient SP140 inhibition unlocks hematopoietic stem cell fate from human pluripotent stem cells. Blood 147, 1584-1597 (2026).
```

In the current analysis, teratoma-derived HSPCs from this previous study were compared with fetal liver HSPC and adult bone marrow HSPC reference datasets. The comparative transcriptomic analysis showed that teratoma-derived HSPCs are transcriptionally more similar to adult bone marrow HSPCs than to fetal liver HSPCs.

The fetal liver HSPC reference was derived from:

```text
Calvanese, V. et al. Mapping human haematopoietic stem cells from haemogenic endothelium to birth. Nature 604, 534-540 (2022).
```

The adult bone marrow HSPC reference was derived from:

```text
Bandyopadhyay, S. et al. Mapping the cellular biogeography of human bone marrow niches using single-cell transcriptomics and proteomic imaging. Cell 187, 3120-3140 e3129 (2024).
```

This analysis supports the rationale that SP140 inhibition during early differentiation may confer adult-like T-cell developmental potential to hPSCs.

The repository provides two versions of the analysis script:

```text
code/HSPC_reference_comparison/HSPC_reference_comparison_SFig1_a_b/original/
  HSPC_reference_comparison_ExtendedData_Fig1ab_full_original_analysis_github.R
```

This script is intended for the full original analysis using the complete Seurat objects for adult bone marrow HSPCs, fetal liver HSPCs and teratoma-derived HSPCs.

```text
code/HSPC_reference_comparison/HSPC_reference_comparison_SFig1_a_b/demo/
  HSPC_reference_comparison_SFig1_a_b_demo_analysis.R
```

This script uses the downsampled demo Seurat objects provided in this repository and is intended to test whether the workflow can be executed successfully.

The demo dataset is a representative reduced version of the full HSPC reference comparison dataset. It is not expected to reproduce the exact quantitative results reported in the manuscript.

---

## Demo data

The HSPC reference comparison demo dataset is located in:

```text
demo_data/HSPC_reference_comparison_SFig1_a_b_demo/
```

It contains the following files:

```text
HSPC_waves.txt
GSE253355_BM_HSPC_demo.rds
GSE162950_fetal_liver_demo.rds
Tera_integrate_Label_demo.rds
```

The demo data include downsampled Seurat objects for:

```text
Adult bone marrow HSPCs
Fetal liver HSPCs
Teratoma-derived HSPCs
```

The fetal liver HSPC groups correspond to the following developmental stages:

```text
FL-HSPC-1: Carnegie Stage 17
FL-HSPC-2: week 8 of gestation
FL-HSPC-3: week 11 of gestation
FL-HSPC-4: week 15 of gestation
```

The adult bone marrow HSPC group is labelled as:

```text
BM-HSPC
```

The teratoma-derived HSPC group is labelled as:

```text
Teratoma-HSPC
```

The demo data are provided only for testing the workflow and checking the expected input/output structure.

Because the demo dataset is a subset of the complete dataset, dot plot patterns and pseudo-bulk clustering results may differ from the full manuscript analysis.

---

## Analysis included in the demo script

The HSPC reference comparison demo script performs the following analyses:

1. Read HSPC wave gene list
2. Read demo adult bone marrow HSPC Seurat object
3. Read demo fetal liver HSPC Seurat object
4. Read demo teratoma-derived HSPC Seurat object
5. Assign HSPC reference group labels
6. Generate HSPC wave gene dot plot
7. Generate pseudo-bulk expression profiles
8. Perform unsupervised hierarchical clustering
9. Save R session information

The demo workflow is designed to reproduce the input/output structure of the full analysis and to verify that the analysis code can be executed successfully.

---

## How to run the HSPC reference comparison demo

Open R or RStudio and set the working directory to the repository root:

```r
setwd("path/to/manuscript-code")
```

Then run:

```r
source("code/HSPC_reference_comparison/HSPC_reference_comparison_SFig1_a_b/demo/HSPC_reference_comparison_SFig1_a_b_demo_analysis.R")
```

The demo results will be written to:

```text
results/HSPC_reference_comparison_SFig1_a_b_demo/
```

Expected output folders include:

```text
00_basic_data/
01_ExtendedData_Fig1a_HSPC_waves_dotplot/
02_ExtendedData_Fig1b_pseudobulk_hclust/
sessionInfo.txt
```

Expected main output files include:

```text
01_ExtendedData_Fig1a_HSPC_waves_dotplot/
  ExtendedData_Fig1a_HSPC_reference_comparison_HSPC_waves_DotPlot.pdf
  ExtendedData_Fig1a_HSPC_reference_comparison_HSPC_waves_DotPlot.png

02_ExtendedData_Fig1b_pseudobulk_hclust/
  ExtendedData_Fig1b_HSPC_reference_comparison_pseudobulk_hclust_tree.pdf
  ExtendedData_Fig1b_HSPC_reference_comparison_pseudobulk_hclust_tree.png
```

---

## Expected runtime

The HSPC reference comparison demo analysis was tested on a Windows 11 laptop.

The demo analysis takes approximately:

```text
1-5 minutes
```

Runtime may vary depending on the computer and whether the required R packages have already been installed.

---

## Software requirements

The workflow was developed and tested in R.

Tested environment:

```text
R version 4.5.1
Platform: x86_64-w64-mingw32/x64
Operating system: Windows 11 x64
```

Main R packages directly used in the HSPC reference comparison script:

```text
Seurat v5.4.0
SeuratObject v5.3.0
ggplot2 v4.0.2
Matrix v1.7-4
pheatmap v1.0.13
```

The script will attempt to install missing CRAN packages automatically.

After running the script, complete R session information, including package dependencies, is saved to:

```text
results/HSPC_reference_comparison_SFig1_a_b_demo/sessionInfo.txt
```

---

## Original full analysis script

The original full HSPC reference comparison analysis script is provided at:

```text
code/HSPC_reference_comparison/HSPC_reference_comparison_SFig1_a_b/original/
  HSPC_reference_comparison_ExtendedData_Fig1ab_full_original_analysis_github.R
```

This script is intended for reproducing the full HSPC reference comparison analysis using the complete Seurat objects.

The original analysis script requires the following full input files:

```text
HSPC_waves.txt
GSE253355_Normal_Bone_Marrow_Atlas_Seurat_SB_v2.rds
GSE162950_seurat_object.rds
Tera_integrate_Label.rds
```

Alternatively, the fetal liver Seurat object can be provided as:

```text
GSE162950seurat_object.Rdata
```

These full input files are not included in this repository due to file size limitations.

Before running the original script, users should modify the following paths according to their local environment:

```r
input_dir <- "path/to/HSPC_reference_comparison_input"

result_dir <- "path/to/output/HSPC_reference_comparison_ExtendedData_Fig1ab"
```

The original script performs the same analysis workflow as the demo script but uses the complete HSPC reference comparison dataset.

---

## Notes on reproducibility

The HSPC reference comparison demo analysis is designed to verify that the code can be executed and that the expected output structure is generated.

Because the demo dataset is a reduced subset of the full dataset:

- HSPC wave gene expression patterns may be similar but not identical to the manuscript figure.
- Pseudo-bulk hierarchical clustering may differ from the full analysis.
- The demo output is intended to demonstrate workflow execution and file structure rather than reproduce the exact final manuscript figure.

This behavior is expected and does not indicate an error in the workflow.

---

## TCR-seq analysis for Extended Data Fig. 3

### Overview

This module contains the TCR-seq analysis code used for Extended Data Fig. 3.

The analysis focuses on TCR gene fragment usage, TCR repertoire diversity and V/J segment usage patterns across hPSC-derived T cells and reference thymic or peripheral T-cell populations.

The repository provides two versions of the analysis script:

```text
code/TCRseq/ExtendedData_Fig3/original/
  TCRseq_ExtendedData_Fig3_original_analysis_GitHub.py
  TCRseq_ExtendedData_Fig3_original_analysis_GitHub.ipynb
```

These scripts are intended for the full original analysis using the complete TCR clone table.

```text
code/TCRseq/ExtendedData_Fig3/demo/
  TCRseq_ExtendedData_Fig3_demo_analysis.py
  TCRseq_ExtendedData_Fig3_demo_analysis.ipynb
  requirements.txt
```

These scripts use the demo dataset provided in this repository and are intended to test whether the TCR-seq workflow can be executed successfully.

The demo dataset is a reduced version of the processed TCR clone table. It is not expected to reproduce the exact quantitative results reported in the manuscript.

---

### Demo data

The TCR-seq demo dataset is located in:

```text
demo_data/TCRseq_ExtendedData_Fig3_demo/
```

It contains the following files:

```text
TCR_Clones_Collated_demo.xlsx
demo_data_summary.csv
```

The demo workbook has the same sheet names and column structure as the full analysis workbook, but each sheet contains a reduced number of rows for testing the workflow.

The sheets used by the demo script are:

```text
PBMC Collated 2 filtered
LD-SC40 Revision
Fetal_w13_FCAImm
Fetal_w14_FCAImm
Fetal_thy_w17
Postnatal_thy_10m
Postnatal_30m_T06
Adult_thy
DMSO-T
Inhibitor-T
```

The sheets PBMC Collated 2 filtered, LD-SC40 Revision, Fetal_w13_FCAImm, Fetal_w14_FCAImm, Fetal_thy_w17, Postnatal_thy_10m, Postnatal_30m_T06 and Adult_thy contain published reference TCR datasets used as in vivo-differentiated reference populations, including primary human thymocytes and PBMC CD8SP T cells. The primary thymocyte references were derived from Park et al. (Science, 2020), and the PBMC CD8SP T-cell reference was derived from Park et al. (Communications Biology, 2023).

The sheets DMSO-T and Inhibitor-T contain the TCR-seq data generated in this study. DMSO-T corresponds to DMSO-treated hPSC-derived CD8 T cells, and Inhibitor-T corresponds to GSK761-treated hPSC-derived CD8 T cells.

The demo data are provided only for testing the workflow and checking the expected input/output structure.

Because the demo dataset is a subset of the complete dataset, TCR diversity metrics, segment usage profiles, PCA results and density plots may differ from the full manuscript analysis.

---

### Analysis included in the demo script

The TCR-seq demo script performs the following analyses:

1. Read the reduced TCR clone workbook
2. Assign sample subsets and cell-type labels
3. Remove records with missing paired CDR3-alpha or CDR3-beta information
4. Clean and standardize TCR V and J gene fragment names
5. Calculate CDR3 amino-acid and nucleotide lengths for downstream repertoire summaries
6. Calculate TCR diversity metrics, including Shannon entropy and Simpson diversity
7. Compare DMSO and inhibitor-treated CD8 T-cell repertoires using a Mann-Whitney U test
8. Summarize TCR V and J segment usage
9. Perform PCA based on combined TRBV, TRBJ, TRAV and TRAJ usage
10. Generate smoothed positional density plots for TCR segment usage
11. Generate TCR V and J segment usage heatmaps
12. Save Python and package version information

Because the demo dataset is a reduced dataset, the output figures are intended to demonstrate code execution and workflow structure rather than reproduce the full manuscript figures exactly.

---

### How to run the TCR-seq demo

The demo can be run either as a Python script or as a Jupyter notebook.

#### Option 1: run the Python script

Open a terminal and set the working directory to the repository root:

```bash
cd path/to/manuscript-code
```

Install the required Python packages:

```bash
pip install -r code/TCRseq/ExtendedData_Fig3/demo/requirements.txt
```

Then run:

```bash
python code/TCRseq/ExtendedData_Fig3/demo/TCRseq_ExtendedData_Fig3_demo_analysis.py
```

The demo results will be written to:

```text
results/TCRseq_ExtendedData_Fig3_demo/
```

#### Option 2: run the notebook

Open the whole repository folder in VS Code or Jupyter, then open and run:

```text
code/TCRseq/ExtendedData_Fig3/demo/TCRseq_ExtendedData_Fig3_demo_analysis.ipynb
```

The notebook uses the same demo data and writes output files to the same results folder:

```text
results/TCRseq_ExtendedData_Fig3_demo/
```

---

### Expected output files

Expected output files include:

```text
software_versions.csv
Fig.A - Simpson diversity.pdf
Fig.C - Va Simple.pdf
Fig.C - Vb Simple.pdf
Fig.D - Ja Simple.pdf
Fig.D - Jb Simple.pdf
Fig.E - Position density.pdf
Fig.F - PCA.pdf
```

The exact number and appearance of output figures may vary depending on the reduced demo data.

---

### Expected runtime

The TCR-seq demo analysis was tested on a Windows 11 computer.

The demo analysis takes approximately:

```text
1-5 minutes
```

Runtime may vary depending on the computer and whether the required Python packages have already been installed.

---

### Software requirements

The workflow was developed and tested in Python.

Tested environment:

```text
Python 3.13.12
Operating system: Windows 11
```

Main Python packages directly used in the script:

```text
numpy
pandas
scipy
matplotlib
seaborn
scikit-learn
openpyxl
ipython
jupyter
```

In the tested environment, the following package versions were recorded:

```text
numpy 2.4.4
pandas 3.0.2
scipy 1.17.1
matplotlib 3.10.9
seaborn 0.13.2
scikit-learn 1.8.0
openpyxl 3.1.5
ipython 9.13.0
```

The script saves the detected Python and package versions to:

```text
results/TCRseq_ExtendedData_Fig3_demo/software_versions.csv
```

---

### Original full analysis script

The original full TCR-seq analysis scripts are provided at:

```text
code/TCRseq/ExtendedData_Fig3/original/
  TCRseq_ExtendedData_Fig3_original_analysis_GitHub.py
  TCRseq_ExtendedData_Fig3_original_analysis_GitHub.ipynb
```

These scripts are intended for reproducing the full TCR gene fragment usage analysis using the complete TCR clone table.

The original analysis requires the complete input workbook:

```text
TCR Clones Collated.xlsx
```

This full input workbook is not included in this repository.

Before running the original analysis scripts, users should modify the input and output paths according to their local environment.

---

### Notes on reproducibility

The TCR-seq demo analysis is designed to verify that the code can be executed and that the expected output structure is generated.

Because the demo dataset is a reduced subset of the full dataset:

- Simpson diversity values may differ from the full analysis.
- TCR V and J segment usage heatmaps may differ from the full manuscript figures.
- PCA patterns may be similar but not identical to the full analysis.
- Positional density plots may differ from the full manuscript figures.

This behavior is expected and does not indicate an error in the workflow.

---

## Additional analysis scripts

To be added.

This section may include additional plotting scripts, integrative RNA-seq/ATAC-seq analysis scripts, gene set scoring scripts or other supporting analyses associated with the manuscript.

---

# Planned additional modules

The following sections are placeholders for additional analysis modules that will be added later.

## Bulk RNA-seq analysis for Fig. 1

To be added.

Planned contents:

```text
code/bulkRNAseq/Fig1/original/
code/bulkRNAseq/Fig1/demo/
demo_data/bulkRNAseq_Fig1_demo/
```

---

## scRNA_seq analysis for SFig. 2

To be added.

Planned contents:

```text
code/scRNA_seq/SFig. 2/original/
code/scRNA_seq/SFig. 2/demo/
demo_data/scRNA_seq_SFig. 2_demo/
```

---

# License

The analysis code in this repository is released under the MIT License. See the LICENSE file for details.

The demo data are provided to test and reproduce the example workflows described in this repository.
