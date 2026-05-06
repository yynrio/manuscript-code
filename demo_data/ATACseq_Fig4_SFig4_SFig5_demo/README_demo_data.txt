ATAC-seq demo data for Fig. 4 / Supplementary Fig. 4 / Supplementary Fig. 5

This folder contains reduced demo data for testing the ATAC-seq analysis workflow.
The demo data are intended to verify code execution and output structure, not to reproduce the exact full manuscript statistics.

Demo data contents:

1. peak_counts_samples_demo.txt
   All samples with the top 5000 highly variable filtered peaks.

2. metadata_demo.csv
   Sample metadata for the demo analysis.

3. IGVgene_demo/
   Demo gene-coordinate table, RefSeq Select annotation subset, and small regional bigWig files for IGV-like plots.

Notes:
- Full manuscript results should be generated using the complete peak count matrix and complete bigWig files.
- This README_demo_data.txt is only for the demo data folder and does not conflict with the repository-level README.md.
