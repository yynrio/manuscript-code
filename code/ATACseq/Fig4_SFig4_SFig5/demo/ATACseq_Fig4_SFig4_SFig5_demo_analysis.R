############################################################
############################################################
## ATAC-seq demo analysis
## Fig. 4 / Supplementary Fig. 4 / Supplementary Fig. 5
##
## This script performs the demo bulk ATAC-seq downstream analysis
## using reduced demo data.
##
## Demo data:
## 1. All samples are retained.
## 2. The peak count matrix contains the top 5000 highly variable
##    filtered peaks.
## 3. Small regional bigWig files are used for IGV-like plots.
##
## Main analysis steps are kept consistent with the full ATAC-seq
## analysis script.
##
## The demo dataset is intended to test the workflow and
## input/output structure. It is not expected to reproduce
## the exact quantitative results in the manuscript.
##
## Main analysis steps:
## 1. Read demo peak-level count matrix and sample metadata
## 2. Clean count matrix and match metadata
## 3. Generate total-count QC plot
## 4. Filter low-count peaks
## 5. Generate six-group PCA plot
## 6. Generate Pearson correlation heatmap based on
##    top 1000 highly variable peaks and group-average signals
## 7. Annotate filtered peaks using ChIPseeker
## 8. Construct DESeq2 object
## 9. Perform pairwise DESeq2 differential accessibility analysis
## 10. Generate pairwise DAR volcano plots
## 11. Summarize significant DAR annotation distributions
## 12. Perform KEGG and GO enrichment analysis
## 13. Generate selected GO pathway bar plots
## 14. Generate IGV-like ATAC signal plots from small regional
##     demo bigWig files


############################################################

rm(list = ls())
gc()


##############################
##  0. Load packages
##############################

packages_cran <- c(
  "ggplot2",
  "ggrepel",
  "data.table",
  "dplyr",
  "tidyr",
  "pheatmap",
  "readxl"
)

packages_bioc <- c(
  "DESeq2",
  "ChIPseeker",
  "GenomicRanges",
  "GenomeInfoDb",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "org.Hs.eg.db",
  "clusterProfiler",
  "Gviz",
  "rtracklayer"
)

for (pkg in packages_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in packages_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
  library(readxl)
  
  library(DESeq2)
  library(ChIPseeker)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(Gviz)
  library(rtracklayer)
})

options(stringsAsFactors = FALSE)
options(ucscChromosomeNames = FALSE)
options(timeout = 1200)
options(download.file.method = "libcurl")

##############################
## 0.1. Input files
##############################

work_dir <- "G:/GitHub/manuscript-code/demo_data/ATACseq_Fig4_SFig4_SFig5_demo"
setwd(work_dir)

count_file <- file.path(
  work_dir,
  "peak_counts_samples_demo.txt"
)

metadata_file <- file.path(
  work_dir,
  "metadata_demo.csv"
)

bw_dir <- file.path(
  work_dir,
  "IGVgene_demo",
  "bigwig_demo"
)

bw_files <- c(
  "DMSO-1"       = file.path(bw_dir, "DMSO-1_demo.bw"),
  "GSK-1"        = file.path(bw_dir, "GSK-1_demo.bw"),
  "PBMC-1"       = file.path(bw_dir, "PBMC-1_demo.bw"),
  "DMSO-CD19-1"  = file.path(bw_dir, "DMSO-CD19-1_demo.bw"),
  "GSK-CD19-1"   = file.path(bw_dir, "GSK-CD19-1_demo.bw"),
  "PBMC-CD19-1"  = file.path(bw_dir, "PBMC-CD19-1_demo.bw")
)

igv_gene_file <- file.path(
  work_dir,
  "IGVgene_demo",
  "graph_Core_genes_list_with_transcript_id_demo.xlsx"
)

refseq_file <- file.path(
  work_dir,
  "IGVgene_demo",
  "ncbiRefSeqSelect_demo.txt.gz"
)

##############################
## 1. Output directories
##############################

result_dir <- "G:/GitHub/manuscript-code/results/ATACseq_Fig4_SFig4_SFig5_demo"

basic_dir   <- file.path(result_dir, "00_basic_data")
qc_dir      <- file.path(result_dir, "00_QC")
pca_dir     <- file.path(result_dir, "01_PCA")
pearson_dir <- file.path(result_dir, "02_Pearson_top1000")

peak_anno_dir    <- file.path(result_dir, "03_peak_annotation")
dar_dir          <- file.path(result_dir, "04_DESeq2_DAR")
volcano_dir      <- file.path(result_dir, "05_pairwise_Volcano")
anno_summary_dir <- file.path(result_dir, "06_DAR_annotation_summary")
kegg_go_dir      <- file.path(result_dir, "07_KEGG_GO")

selected_pathway_dir <- file.path(result_dir, "08_Selected_pathway")
igv_dir <- file.path(result_dir, "09_IGV_like_bigWig_plot")

for (d in c(
  result_dir,
  basic_dir,
  qc_dir,
  pca_dir,
  pearson_dir,
  peak_anno_dir,
  dar_dir,
  volcano_dir,
  anno_summary_dir,
  kegg_go_dir,
  selected_pathway_dir,
  igv_dir
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

cat("Demo input directory:\n")
cat(work_dir, "\n\n")

cat("Demo output directory:\n")
cat(result_dir, "\n")

############################################################
## 2. Read metadata
############################################################

meta <- read.csv(
  metadata_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

need_cols <- c(
  "sample_id",
  "source",
  "stim",
  "replicate",
  "group"
)

missing_cols <- setdiff(need_cols, colnames(meta))

if (length(missing_cols) > 0) {
  stop(
    paste0(
      "metadata is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  )
}

if (anyDuplicated(meta$sample_id) > 0) {
  stop("metadata$sample_id contains duplicated sample names.")
}

meta$sample_id <- as.character(meta$sample_id)
meta$source <- as.character(meta$source)
meta$stim <- as.character(meta$stim)
meta$replicate <- factor(meta$replicate)
meta$group_raw <- as.character(meta$group)

cat("metadata dimension:\n")
print(dim(meta))

cat("metadata columns:\n")
print(colnames(meta))

############################################################
## 3. Read peak count matrix
############################################################

count_df <- data.table::fread(
  count_file,
  data.table = FALSE,
  check.names = FALSE
)

cat("Raw peak count matrix dimension:\n")
print(dim(count_df))

cat("Raw peak count matrix first columns:\n")
print(colnames(count_df)[1:min(10, ncol(count_df))])

## Remove empty column generated by export, if present.
if ("V22" %in% colnames(count_df)) {
  count_df <- count_df[, !colnames(count_df) %in% "V22", drop = FALSE]
  cat("Removed empty column V22.\n")
}

need_peak_cols <- c("chr", "start", "end")
missing_peak_cols <- setdiff(need_peak_cols, colnames(count_df))

if (length(missing_peak_cols) > 0) {
  stop(
    paste0(
      "peak count matrix is missing required peak columns: ",
      paste(missing_peak_cols, collapse = ", ")
    )
  )
}

count_df$chr <- as.character(count_df$chr)
count_df$start <- suppressWarnings(as.numeric(count_df$start))
count_df$end <- suppressWarnings(as.numeric(count_df$end))

############################################################
## 4. Extract count matrix and match metadata
############################################################

sample_cols <- which(colnames(count_df) %in% meta$sample_id)

cat("Detected sample columns:\n")
print(colnames(count_df)[sample_cols])

if (length(sample_cols) != nrow(meta)) {
  
  missing_in_count <- setdiff(meta$sample_id, colnames(count_df))
  missing_in_meta <- setdiff(colnames(count_df)[sample_cols], meta$sample_id)
  
  cat("Samples in metadata but not count matrix:\n")
  print(missing_in_count)
  
  cat("Samples in count matrix but not metadata:\n")
  print(missing_in_meta)
  
  stop("Sample columns do not match metadata rows.")
}

count_mat <- count_df[, sample_cols, drop = FALSE]
count_mat <- as.matrix(count_mat)
mode(count_mat) <- "numeric"

peak_ids <- paste(
  count_df$chr,
  count_df$start,
  count_df$end,
  sep = "_"
)

if (anyDuplicated(peak_ids) > 0) {
  stop("Duplicated peak IDs were detected. Please check chr/start/end columns.")
}

rownames(count_mat) <- peak_ids

## Reorder count matrix according to metadata.
count_mat <- count_mat[, meta$sample_id, drop = FALSE]

if (!all(colnames(count_mat) == meta$sample_id)) {
  stop("Sample order mismatch after reordering count matrix.")
}

rownames(meta) <- meta$sample_id
meta <- meta[colnames(count_mat), , drop = FALSE]

cat("count_mat dimension:\n")
print(dim(count_mat))

cat("count_mat sample names:\n")
print(colnames(count_mat))

############################################################
## 5. Set metadata factors and six-group names
############################################################

meta$source_raw <- as.character(meta$source)
meta$source_raw[meta$source_raw == "GSK761"] <- "GSK"

meta$source <- factor(
  meta$source_raw,
  levels = c("DMSO", "GSK", "PBMC")
)

meta$stim <- factor(
  as.character(meta$stim),
  levels = c("Unstim", "Stim")
)

meta$replicate <- factor(meta$replicate)

meta$group <- ifelse(
  meta$source_raw == "DMSO" & meta$stim == "Unstim", "DMSO_un",
  ifelse(
    meta$source_raw == "GSK" & meta$stim == "Unstim", "GSK_un",
    ifelse(
      meta$source_raw == "PBMC" & meta$stim == "Unstim", "PBMC_un",
      ifelse(
        meta$source_raw == "DMSO" & meta$stim == "Stim", "DMSO_CD19",
        ifelse(
          meta$source_raw == "GSK" & meta$stim == "Stim", "GSK_CD19",
          ifelse(
            meta$source_raw == "PBMC" & meta$stim == "Stim", "PBMC_CD19",
            NA
          )
        )
      )
    )
  )
)

if (any(is.na(meta$group))) {
  stop("meta$group contains NA. Please check source/stim columns.")
}

group_levels_atac <- c(
  "DMSO_un",
  "GSK_un",
  "PBMC_un",
  "DMSO_CD19",
  "GSK_CD19",
  "PBMC_CD19"
)

meta$group <- factor(
  meta$group,
  levels = group_levels_atac
)

cat("source × stim table:\n")
print(table(meta$source, meta$stim))

cat("ATAC six-group table:\n")
print(table(meta$group))

############################################################
## 6. Save basic data before filtering
############################################################

peak_info_raw <- data.frame(
  Peakid = rownames(count_mat),
  chr = count_df$chr,
  start = count_df$start,
  end = count_df$end,
  stringsAsFactors = FALSE
)

saveRDS(
  count_mat,
  file = file.path(basic_dir, "peak_counts_clean_matrix_raw.rds")
)

saveRDS(
  meta,
  file = file.path(basic_dir, "metadata_clean_ATAC.rds")
)

saveRDS(
  peak_info_raw,
  file = file.path(basic_dir, "peak_id_coordinate_map_raw.rds")
)

write.csv(
  meta,
  file = file.path(basic_dir, "metadata_clean_ATAC.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  peak_info_raw,
  file = file.path(basic_dir, "peak_id_coordinate_map_raw.csv"),
  row.names = FALSE,
  quote = FALSE
)

############################################################
## 7. QC plots: mean ± SEM bar plot with sample dots
############################################################

sample_library_size <- colSums(count_mat)
detected_peaks <- colSums(count_mat > 0)

qc_basic_atac <- data.frame(
  sample_id = colnames(count_mat),
  total_counts = sample_library_size,
  detected_peaks = detected_peaks,
  source = meta$source,
  stim = meta$stim,
  replicate = meta$replicate,
  group = meta$group,
  stringsAsFactors = FALSE
)

qc_basic_atac$group <- factor(
  qc_basic_atac$group,
  levels = group_levels_atac
)

write.csv(
  qc_basic_atac,
  file = file.path(qc_dir, "ATAC_qc_basic_sample_summary.csv"),
  row.names = FALSE,
  quote = FALSE
)

qc_summary_atac <- qc_basic_atac %>%
  group_by(group) %>%
  summarise(
    mean_total_counts = mean(total_counts),
    sem_total_counts = sd(total_counts) / sqrt(n()),
    mean_detected_peaks = mean(detected_peaks),
    sem_detected_peaks = sd(detected_peaks) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

write.csv(
  qc_summary_atac,
  file = file.path(qc_dir, "ATAC_qc_group_mean_sem_summary.csv"),
  row.names = FALSE,
  quote = FALSE
)

group_colors <- c(
  "DMSO_un" = "#7b95c6",
  "GSK_un" = "#49c2d9",
  "PBMC_un" = "#a1d8e8",
  "DMSO_CD19" = "#67a583",
  "GSK_CD19" = "#a2c986",
  "PBMC_CD19" = "#d0e2c0"
)

p_total_counts_atac <- ggplot() +
  geom_col(
    data = qc_summary_atac,
    aes(x = group, y = mean_total_counts, fill = group),
    color = "black",
    width = 0.7
  ) +
  geom_errorbar(
    data = qc_summary_atac,
    aes(
      x = group,
      ymin = mean_total_counts - sem_total_counts,
      ymax = mean_total_counts + sem_total_counts
    ),
    width = 0.2,
    linewidth = 0.5
  ) +
  geom_jitter(
    data = qc_basic_atac,
    aes(x = group, y = total_counts),
    width = 0.12,
    size = 1.8,
    color = "black"
  ) +
  scale_fill_manual(values = group_colors) +
  theme_bw(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major.x = element_blank()
  ) +
  xlab(NULL) +
  ylab("Total counts (mean ± SEM)") +
  ggtitle("ATAC-seq")

ggsave(
  filename = file.path(qc_dir, "ATAC_QC_total_counts_mean_SEM.png"),
  plot = p_total_counts_atac,
  width = 5.5,
  height = 4.5,
  dpi = 300
)

ggsave(
  filename = file.path(qc_dir, "ATAC_QC_total_counts_mean_SEM.pdf"),
  plot = p_total_counts_atac,
  width = 5.5,
  height = 4.5
)

############################################################
## 8. Filter low-count peaks
############################################################
## Keep peaks with count >= 10 in at least 3 samples.

keep <- rowSums(count_mat >= 10) >= 3
count_mat_filt <- count_mat[keep, , drop = FALSE]

cat("Peaks before filtering:", nrow(count_mat), "\n")
cat("Peaks after filtering:", nrow(count_mat_filt), "\n")
cat("Retained proportion:", round(nrow(count_mat_filt) / nrow(count_mat) * 100, 2), "%\n")

filter_summary <- data.frame(
  filter_rule = "count >= 10 in at least 3 samples",
  peaks_before_filtering = nrow(count_mat),
  peaks_after_filtering = nrow(count_mat_filt),
  peaks_removed = sum(!keep),
  retained_percent = round(nrow(count_mat_filt) / nrow(count_mat) * 100, 2),
  stringsAsFactors = FALSE
)

write.csv(
  filter_summary,
  file = file.path(qc_dir, "ATAC_peak_filtering_summary.csv"),
  row.names = FALSE,
  quote = FALSE
)

peak_info_filt <- peak_info_raw[
  match(rownames(count_mat_filt), peak_info_raw$Peakid),
  ,
  drop = FALSE
]

if (!all(peak_info_filt$Peakid == rownames(count_mat_filt))) {
  stop("Filtered peak information does not match filtered count matrix.")
}

saveRDS(
  count_mat_filt,
  file = file.path(basic_dir, "peak_counts_clean_matrix_filtered.rds")
)

saveRDS(
  peak_info_filt,
  file = file.path(basic_dir, "peak_id_coordinate_map_filtered.rds")
)

write.csv(
  count_mat_filt,
  file = file.path(basic_dir, "peak_counts_clean_matrix_filtered.csv"),
  quote = FALSE
)

write.csv(
  peak_info_filt,
  file = file.path(basic_dir, "peak_id_coordinate_map_filtered.csv"),
  row.names = FALSE,
  quote = FALSE
)

############################################################
## 9. log2 transformation
############################################################
## log2(count + 1)

log_mat <- log2(count_mat_filt + 1)

saveRDS(
  log_mat,
  file = file.path(basic_dir, "peak_counts_log2_count_plus_1_filtered.rds")
)

cat("log_mat dimension:\n")
print(dim(log_mat))

############################################################
## 10. Six-group PCA using top 5000 variable peaks
############################################################

meta_pca <- meta[colnames(log_mat), , drop = FALSE]

if (!all(rownames(meta_pca) == colnames(log_mat))) {
  stop("meta_pca order does not match log_mat columns.")
}

peak_var_pca <- apply(log_mat, 1, var, na.rm = TRUE)
peak_var_pca <- sort(peak_var_pca, decreasing = TRUE)

top_n_pca <- min(5000, length(peak_var_pca))
top_peaks_pca <- names(peak_var_pca)[1:top_n_pca]

pca_input <- log_mat[top_peaks_pca, , drop = FALSE]

cat("Number of variable peaks used for ATAC six-group PCA:", nrow(pca_input), "\n")

pca_res <- prcomp(
  t(pca_input),
  center = TRUE,
  scale. = FALSE
)

percent_var <- round(
  100 * (pca_res$sdev^2 / sum(pca_res$sdev^2))[1:3],
  2
)

pca_data <- data.frame(
  sample_id = rownames(pca_res$x),
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  PC3 = pca_res$x[, 3],
  source = meta_pca$source,
  stim = meta_pca$stim,
  group = meta_pca$group,
  stringsAsFactors = FALSE
)

pca_data$source <- factor(
  pca_data$source,
  levels = c("DMSO", "GSK", "PBMC")
)

pca_data$stim <- factor(
  pca_data$stim,
  levels = c("Unstim", "Stim")
)

pca_data$group <- factor(
  pca_data$group,
  levels = group_levels_atac
)

source_colors <- c(
  "DMSO" = "black",
  "GSK" = "red",
  "PBMC" = "blue"
)

stim_shapes <- c(
  "Unstim" = 16,
  "Stim" = 17
)

p_pca <- ggplot(
  pca_data,
  aes(PC1, PC2, color = source, shape = stim)
) +
  geom_point(size = 2.5) +
  scale_color_manual(values = source_colors) +
  scale_shape_manual(values = stim_shapes) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  theme_bw(base_size = 14) +
  ggtitle("ATAC-seq")

ggsave(
  filename = file.path(pca_dir, "ATAC_PCA_six_groups_top5000.png"),
  plot = p_pca,
  width = 6,
  height = 4,
  dpi = 300
)

ggsave(
  filename = file.path(pca_dir, "ATAC_PCA_six_groups_top5000.pdf"),
  plot = p_pca,
  width = 6,
  height = 4
)

write.csv(
  pca_data,
  file = file.path(pca_dir, "ATAC_PCA_six_groups_top5000_coordinates.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  data.frame(
    Peakid = names(peak_var_pca),
    variance = as.numeric(peak_var_pca)
  ),
  file = file.path(pca_dir, "ATAC_PCA_peak_variance_rank.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  data.frame(Peakid = top_peaks_pca),
  file = file.path(pca_dir, "ATAC_PCA_top5000_variable_peaks.csv"),
  row.names = FALSE,
  quote = FALSE
)

cat("ATAC six-group PCA finished.\n")
cat("ATAC PCA percent variance:\n")
print(percent_var)

############################################################
## 11. Top1000 highly variable peak Pearson correlation heatmap
############################################################

meta_pearson <- meta[colnames(log_mat), , drop = FALSE]

meta_pearson$group <- factor(
  meta_pearson$group,
  levels = group_levels_atac
)

meta_pearson <- meta_pearson[
  order(meta_pearson$group),
  ,
  drop = FALSE
]

log_mat_6group <- log_mat[, rownames(meta_pearson), drop = FALSE]

if (!all(colnames(log_mat_6group) == rownames(meta_pearson))) {
  stop("log_mat_6group order does not match meta_pearson rows.")
}

peak_var <- apply(log_mat_6group, 1, var, na.rm = TRUE)

top1000_peaks <- names(sort(peak_var, decreasing = TRUE))[
  1:min(1000, length(peak_var))
]

log_top1000 <- log_mat_6group[top1000_peaks, , drop = FALSE]

group_mean_mat <- sapply(group_levels_atac, function(g) {
  samples_g <- rownames(meta_pearson)[meta_pearson$group == g]
  rowMeans(log_top1000[, samples_g, drop = FALSE])
})

colnames(group_mean_mat) <- group_levels_atac

pearson_cor <- cor(
  group_mean_mat,
  method = "pearson",
  use = "pairwise.complete.obs"
)

write.csv(
  group_mean_mat,
  file.path(pearson_dir, "ATAC_top1000_group_mean_log2_matrix.csv"),
  quote = FALSE
)

write.csv(
  pearson_cor,
  file.path(pearson_dir, "ATAC_top1000_group_mean_Pearson_correlation.csv"),
  quote = FALSE
)

write.csv(
  data.frame(
    Peakid = names(sort(peak_var, decreasing = TRUE)),
    variance = as.numeric(sort(peak_var, decreasing = TRUE))
  ),
  file = file.path(pearson_dir, "ATAC_6group_peak_variance_rank.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  data.frame(Peakid = top1000_peaks),
  file = file.path(pearson_dir, "ATAC_6group_top1000_variable_peaks.csv"),
  row.names = FALSE,
  quote = FALSE
)

png(
  filename = file.path(pearson_dir, "ATAC_top1000_group_mean_Pearson_heatmap.png"),
  width = 2200,
  height = 2000,
  res = 300
)

pheatmap(
  pearson_cor,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.3f",
  border_color = "white",
  main = "ATAC-seq"
)

dev.off()

pdf(
  file = file.path(pearson_dir, "ATAC_top1000_group_mean_Pearson_heatmap.pdf"),
  width = 7,
  height = 6.5
)

pheatmap(
  pearson_cor,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  display_numbers = TRUE,
  number_format = "%.3f",
  border_color = "white",
  main = "ATAC-seq"
)

dev.off()

############################################################
## 12. Load packages for peak annotation and DAR analysis
############################################################

packages_cran_extra <- c(
  "ggplot2",
  "ggrepel",
  "dplyr",
  "tidyr"
)

packages_bioc <- c(
  "DESeq2",
  "ChIPseeker",
  "GenomicRanges",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "org.Hs.eg.db"
)

for (pkg in packages_cran_extra) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

for (pkg in packages_bioc) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(DESeq2)
  library(ChIPseeker)
  library(GenomicRanges)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
})

############################################################
## 13. Output directories for peak annotation, DAR and volcano
############################################################

peak_anno_dir <- file.path(result_dir, "03_peak_annotation")
dar_dir <- file.path(result_dir, "04_DESeq2_DAR")
volcano_dir <- file.path(result_dir, "05_pairwise_Volcano")
anno_summary_dir <- file.path(result_dir, "06_DAR_annotation_summary")

dir.create(peak_anno_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dar_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(anno_summary_dir, recursive = TRUE, showWarnings = FALSE)

dar_unstim_dir <- file.path(dar_dir, "unstimulated_three_groups")
dar_stim_dir <- file.path(dar_dir, "CD19_stimulated_three_groups")
dar_stim_unstim_dir <- file.path(dar_dir, "stim_vs_unstim")

dir.create(dar_unstim_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dar_stim_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dar_stim_unstim_dir, recursive = TRUE, showWarnings = FALSE)

dar_unstim_all_dir <- file.path(dar_unstim_dir, "DAR_all_ATAC")
dar_unstim_sig_dir <- file.path(dar_unstim_dir, "DAR_sig_ATAC")

dar_stim_all_dir <- file.path(dar_stim_dir, "DAR_all_ATAC")
dar_stim_sig_dir <- file.path(dar_stim_dir, "DAR_sig_ATAC")

dar_stim_unstim_all_dir <- file.path(dar_stim_unstim_dir, "DAR_all_ATAC")
dar_stim_unstim_sig_dir <- file.path(dar_stim_unstim_dir, "DAR_sig_ATAC")

dir.create(dar_unstim_all_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dar_unstim_sig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dar_stim_all_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dar_stim_sig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dar_stim_unstim_all_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(dar_stim_unstim_sig_dir, recursive = TRUE, showWarnings = FALSE)

volcano_unstim_dir <- file.path(volcano_dir, "unstimulated_three_groups")
volcano_stim_dir <- file.path(volcano_dir, "CD19_stimulated_three_groups")
volcano_stim_unstim_dir <- file.path(volcano_dir, "stim_vs_unstim")

dir.create(volcano_unstim_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(volcano_stim_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(volcano_stim_unstim_dir, recursive = TRUE, showWarnings = FALSE)

anno_unstim_dir <- file.path(anno_summary_dir, "unstimulated_three_groups")
anno_stim_dir <- file.path(anno_summary_dir, "CD19_stimulated_three_groups")
anno_stim_unstim_dir <- file.path(anno_summary_dir, "stim_vs_unstim")

dir.create(anno_unstim_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(anno_stim_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(anno_stim_unstim_dir, recursive = TRUE, showWarnings = FALSE)

############################################################
## 14. Helper functions for safe export
############################################################

sanitize_for_excel <- function(df) {
  
  df[] <- lapply(df, function(x) {
    
    if (is.factor(x)) {
      x <- as.character(x)
    }
    
    if (is.character(x)) {
      x <- gsub("\r\n|\r|\n", " ", x)
      x <- gsub("\t", " ", x)
      x <- trimws(x)
    }
    
    return(x)
  })
  
  return(df)
}

safe_write_csv <- function(df, file) {
  
  df <- sanitize_for_excel(df)
  
  write.csv(
    df,
    file = file,
    row.names = FALSE,
    quote = TRUE,
    na = "",
    fileEncoding = "UTF-8"
  )
}

safe_write_tsv <- function(df, file) {
  
  df <- sanitize_for_excel(df)
  
  write.table(
    df,
    file = file,
    sep = "\t",
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE,
    na = "",
    fileEncoding = "UTF-8"
  )
}

############################################################
## 15. Peak annotation categories
############################################################

anno_levels <- c(
  "Promoter (<=1kb)",
  "Promoter (1-2kb)",
  "Promoter (2-3kb)",
  "5' UTR",
  "3' UTR",
  "1st Exon",
  "Other Exon",
  "1st Intron",
  "Other Intron",
  "Downstream",
  "Distal Intergenic"
)

anno_merged_levels <- c(
  "Promoter (<=3kb)",
  "UTR",
  "Exon",
  "Intron",
  "Downstream",
  "Distal Intergenic"
)

simplify_annotation <- function(x) {
  
  x <- as.character(x)
  
  out <- dplyr::case_when(
    grepl("^Promoter \\(<=1kb\\)", x) ~ "Promoter (<=1kb)",
    grepl("^Promoter \\(1-2kb\\)", x) ~ "Promoter (1-2kb)",
    grepl("^Promoter \\(2-3kb\\)", x) ~ "Promoter (2-3kb)",
    grepl("^5' UTR", x) ~ "5' UTR",
    grepl("^3' UTR", x) ~ "3' UTR",
    grepl("Exon", x) & grepl("exon 1 of", x) ~ "1st Exon",
    grepl("Exon", x) ~ "Other Exon",
    grepl("Intron", x) & grepl("intron 1 of", x) ~ "1st Intron",
    grepl("Intron", x) ~ "Other Intron",
    grepl("^Downstream", x) ~ "Downstream",
    grepl("^Distal Intergenic", x) ~ "Distal Intergenic",
    TRUE ~ NA_character_
  )
  
  return(out)
}

merge_annotation <- function(x) {
  
  x <- as.character(x)
  
  out <- dplyr::case_when(
    x %in% c(
      "Promoter (<=1kb)",
      "Promoter (1-2kb)",
      "Promoter (2-3kb)"
    ) ~ "Promoter (<=3kb)",
    x %in% c("5' UTR", "3' UTR") ~ "UTR",
    x %in% c("1st Exon", "Other Exon") ~ "Exon",
    x %in% c("1st Intron", "Other Intron") ~ "Intron",
    x == "Downstream" ~ "Downstream",
    x == "Distal Intergenic" ~ "Distal Intergenic",
    TRUE ~ NA_character_
  )
  
  return(out)
}

############################################################
## 16. Annotate all filtered peaks using ChIPseeker
############################################################

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

peak_info_for_anno <- peak_info_filt

peak_info_for_anno$chr_for_anno <- as.character(peak_info_for_anno$chr)

peak_info_for_anno$chr_for_anno[
  !grepl("^chr", peak_info_for_anno$chr_for_anno)
] <- paste0(
  "chr",
  peak_info_for_anno$chr_for_anno[
    !grepl("^chr", peak_info_for_anno$chr_for_anno)
  ]
)

peak_info_for_anno <- peak_info_for_anno[
  !is.na(peak_info_for_anno$chr_for_anno) &
    !is.na(peak_info_for_anno$start) &
    !is.na(peak_info_for_anno$end),
  ,
  drop = FALSE
]

gr_all_peaks <- GRanges(
  seqnames = peak_info_for_anno$chr_for_anno,
  ranges = IRanges(
    start = peak_info_for_anno$start,
    end = peak_info_for_anno$end
  )
)

mcols(gr_all_peaks)$Peakid <- peak_info_for_anno$Peakid

cat("Running ChIPseeker annotation for filtered ATAC peaks...\n")

peakAnno_all <- annotatePeak(
  peak = gr_all_peaks,
  TxDb = txdb,
  tssRegion = c(-3000, 3000),
  annoDb = "org.Hs.eg.db"
)

peak_anno_df <- as.data.frame(peakAnno_all)

if (!"Peakid" %in% colnames(peak_anno_df)) {
  peak_anno_df$Peakid <- mcols(gr_all_peaks)$Peakid
}

anno_keep_cols <- c(
  "Peakid",
  "annotation",
  "geneChr",
  "geneStart",
  "geneEnd",
  "geneLength",
  "geneStrand",
  "geneId",
  "transcriptId",
  "distanceToTSS",
  "SYMBOL",
  "GENENAME"
)

anno_keep_cols <- anno_keep_cols[
  anno_keep_cols %in% colnames(peak_anno_df)
]

peak_anno_df <- peak_anno_df[, anno_keep_cols, drop = FALSE]

peak_anno_df$anno_simple <- simplify_annotation(peak_anno_df$annotation)

peak_anno_df$anno_simple <- factor(
  peak_anno_df$anno_simple,
  levels = anno_levels
)

peak_anno_df$anno_simple <- as.character(peak_anno_df$anno_simple)

peak_anno_df$anno_merged <- merge_annotation(peak_anno_df$anno_simple)

peak_anno_df$anno_merged <- factor(
  peak_anno_df$anno_merged,
  levels = anno_merged_levels
)

peak_anno_df$anno_merged <- as.character(peak_anno_df$anno_merged)

peak_annotation_all <- dplyr::left_join(
  peak_info_filt,
  peak_anno_df,
  by = "Peakid"
)

safe_write_csv(
  peak_annotation_all,
  file = file.path(
    peak_anno_dir,
    "ATAC_all_filtered_peaks_ChIPseeker_annotation.csv"
  )
)

saveRDS(
  peak_annotation_all,
  file = file.path(
    peak_anno_dir,
    "ATAC_all_filtered_peaks_ChIPseeker_annotation.rds"
  )
)

safe_write_csv(
  peak_annotation_all,
  file = file.path(
    basic_dir,
    "ATAC_all_filtered_peaks_ChIPseeker_annotation.csv"
  )
)

saveRDS(
  peak_annotation_all,
  file = file.path(
    basic_dir,
    "ATAC_all_filtered_peaks_ChIPseeker_annotation.rds"
  )
)

cat("Peak annotation finished.\n")

############################################################
## 17. Summarize annotation distribution for all filtered peaks
############################################################

summarise_annotation_distribution <- function(df,
                                              comparison,
                                              result_type,
                                              anno_col,
                                              anno_levels_use,
                                              summary_type) {
  
  if (is.null(df) || nrow(df) == 0 || !anno_col %in% colnames(df)) {
    return(data.frame())
  }
  
  df_use <- df %>%
    dplyr::filter(
      !is.na(.data[[anno_col]]),
      .data[[anno_col]] %in% anno_levels_use
    )
  
  if (nrow(df_use) == 0) {
    return(data.frame())
  }
  
  if ("direction" %in% colnames(df_use)) {
    
    out <- df_use %>%
      dplyr::mutate(
        comparison = comparison,
        result_type = result_type,
        direction = as.character(direction),
        annotation_summary_type = summary_type,
        annotation_category = .data[[anno_col]]
      ) %>%
      dplyr::count(
        comparison,
        result_type,
        direction,
        annotation_summary_type,
        annotation_category,
        name = "peak_count"
      ) %>%
      dplyr::group_by(
        comparison,
        result_type,
        direction,
        annotation_summary_type
      ) %>%
      dplyr::mutate(
        percent = peak_count / sum(peak_count) * 100
      ) %>%
      dplyr::ungroup()
    
  } else {
    
    out <- df_use %>%
      dplyr::mutate(
        comparison = comparison,
        result_type = result_type,
        annotation_summary_type = summary_type,
        annotation_category = .data[[anno_col]]
      ) %>%
      dplyr::count(
        comparison,
        result_type,
        annotation_summary_type,
        annotation_category,
        name = "peak_count"
      ) %>%
      dplyr::group_by(
        comparison,
        result_type,
        annotation_summary_type
      ) %>%
      dplyr::mutate(
        percent = peak_count / sum(peak_count) * 100
      ) %>%
      dplyr::ungroup()
  }
  
  out$annotation_category <- factor(
    out$annotation_category,
    levels = anno_levels_use
  )
  
  out <- out[
    order(
      out$comparison,
      out$result_type,
      out$annotation_summary_type,
      out$annotation_category
    ),
    ,
    drop = FALSE
  ]
  
  out$annotation_category <- as.character(out$annotation_category)
  
  return(out)
}

check_promoter_merge <- function(df,
                                 comparison,
                                 result_type) {
  
  if (is.null(df) || nrow(df) == 0) {
    return(data.frame())
  }
  
  if (!all(c("anno_simple", "anno_merged") %in% colnames(df))) {
    return(data.frame())
  }
  
  df_use <- df
  df_use$comparison <- comparison
  df_use$result_type <- result_type
  
  group_cols <- c("comparison", "result_type")
  
  if ("direction" %in% colnames(df_use)) {
    df_use$direction <- as.character(df_use$direction)
    group_cols <- c(group_cols, "direction")
  }
  
  split_fac <- do.call(
    interaction,
    c(df_use[group_cols], drop = TRUE, sep = "___")
  )
  
  split_list <- split(df_use, split_fac)
  
  out_list <- lapply(split_list, function(one) {
    
    promoter_simple_sum <- sum(
      one$anno_simple %in% c(
        "Promoter (<=1kb)",
        "Promoter (1-2kb)",
        "Promoter (2-3kb)"
      ),
      na.rm = TRUE
    )
    
    promoter_merged_count <- sum(
      one$anno_merged == "Promoter (<=3kb)",
      na.rm = TRUE
    )
    
    base <- data.frame(
      comparison = unique(one$comparison)[1],
      result_type = unique(one$result_type)[1],
      promoter_simple_sum = promoter_simple_sum,
      promoter_merged_count = promoter_merged_count,
      check_pass = promoter_simple_sum == promoter_merged_count,
      stringsAsFactors = FALSE
    )
    
    if ("direction" %in% colnames(one)) {
      base$direction <- unique(one$direction)[1]
    }
    
    return(base)
  })
  
  out <- dplyr::bind_rows(out_list)
  
  return(out)
}

all_peaks_fine_summary <- summarise_annotation_distribution(
  df = peak_annotation_all,
  comparison = "all_filtered_peaks",
  result_type = "all_filtered_peaks",
  anno_col = "anno_simple",
  anno_levels_use = anno_levels,
  summary_type = "fine_annotation"
)

all_peaks_merged_summary <- summarise_annotation_distribution(
  df = peak_annotation_all,
  comparison = "all_filtered_peaks",
  result_type = "all_filtered_peaks",
  anno_col = "anno_merged",
  anno_levels_use = anno_merged_levels,
  summary_type = "merged_annotation"
)

all_peaks_promoter_check <- check_promoter_merge(
  df = peak_annotation_all,
  comparison = "all_filtered_peaks",
  result_type = "all_filtered_peaks"
)

safe_write_csv(
  all_peaks_fine_summary,
  file = file.path(
    anno_summary_dir,
    "ATAC_all_filtered_peaks_fine_annotation_summary.csv"
  )
)

safe_write_csv(
  all_peaks_merged_summary,
  file = file.path(
    anno_summary_dir,
    "ATAC_all_filtered_peaks_merged_annotation_summary_no_Others.csv"
  )
)

safe_write_csv(
  all_peaks_promoter_check,
  file = file.path(
    anno_summary_dir,
    "ATAC_all_filtered_peaks_promoter_merge_check.csv"
  )
)

############################################################
## 18. Construct DESeq2 object
############################################################

meta_deseq <- meta

rownames(meta_deseq) <- meta_deseq$sample_id

meta_deseq <- meta_deseq[
  colnames(count_mat_filt),
  ,
  drop = FALSE
]

meta_deseq$group <- factor(
  meta_deseq$group,
  levels = group_levels_atac
)

if (!all(rownames(meta_deseq) == colnames(count_mat_filt))) {
  stop("meta_deseq order does not match count_mat_filt columns.")
}

cat("ATAC group levels:\n")
print(levels(meta_deseq$group))

cat("ATAC group table:\n")
print(table(meta_deseq$group))

dds_atac <- DESeqDataSetFromMatrix(
  countData = round(count_mat_filt),
  colData = meta_deseq,
  design = ~ group
)

cat("Original ATAC dds dimension:\n")
print(dim(dds_atac))

cat("Running DESeq2 model for ATAC DAR analysis...\n")

dds_atac <- DESeq(dds_atac)

saveRDS(
  dds_atac,
  file = file.path(basic_dir, "ATAC_dds_DESeq_done.rds")
)

writeLines(
  text = resultsNames(dds_atac),
  con = file.path(basic_dir, "ATAC_resultsNames.txt")
)

writeLines(
  text = levels(meta_deseq$group),
  con = file.path(basic_dir, "ATAC_group_levels.txt")
)

safe_write_csv(
  as.data.frame(table(meta_deseq$group)),
  file = file.path(basic_dir, "ATAC_group_table.csv")
)

cat("ATAC DESeq2 model finished.\n")
print(resultsNames(dds_atac))

############################################################
## 19. Volcano plot and BED export functions
############################################################

export_sig_bed_atac <- function(sig_df, outdir_sig, out_name) {
  
  if (nrow(sig_df) == 0) {
    return(NULL)
  }
  
  bed_df <- data.frame(
    chr = sig_df$chr,
    start = pmax(sig_df$start - 1, 0),
    end = sig_df$end,
    Peakid = sig_df$Peakid,
    log2FoldChange = sig_df$log2FoldChange,
    padj = sig_df$padj,
    direction = sig_df$direction,
    anno_simple = if ("anno_simple" %in% colnames(sig_df)) sig_df$anno_simple else NA,
    anno_merged = if ("anno_merged" %in% colnames(sig_df)) sig_df$anno_merged else NA,
    SYMBOL = if ("SYMBOL" %in% colnames(sig_df)) sig_df$SYMBOL else NA,
    stringsAsFactors = FALSE
  )
  
  safe_write_tsv(
    bed_df,
    file = file.path(outdir_sig, paste0("ATAC_DAR_", out_name, "_sig.bed"))
  )
}

plot_volcano_atac <- function(res_df,
                              comp_name,
                              outdir_plot,
                              padj_cutoff = 0.05,
                              lfc_cutoff = 1,
                              top_n_label = 0) {
  
  plot_df <- res_df
  
  plot_df <- plot_df[
    !is.na(plot_df$log2FoldChange) &
      !is.na(plot_df$padj),
    ,
    drop = FALSE
  ]
  
  if (nrow(plot_df) == 0) {
    cat("Volcano skipped:", comp_name, "has no available data.\n")
    return(NULL)
  }
  
  plot_df$change <- "NS"
  
  plot_df$change[
    plot_df$padj < padj_cutoff &
      plot_df$log2FoldChange >= lfc_cutoff
  ] <- "Up"
  
  plot_df$change[
    plot_df$padj < padj_cutoff &
      plot_df$log2FoldChange <= -lfc_cutoff
  ] <- "Down"
  
  plot_df$change <- factor(
    plot_df$change,
    levels = c("Up", "Down", "NS")
  )
  
  plot_df$neglog10padj <- -log10(plot_df$padj)
  
  finite_vals <- plot_df$neglog10padj[is.finite(plot_df$neglog10padj)]
  finite_max <- if (length(finite_vals) > 0) max(finite_vals, na.rm = TRUE) else 1
  plot_df$neglog10padj[!is.finite(plot_df$neglog10padj)] <- finite_max + 1
  
  plot_df$label <- ""
  
  sig_df <- plot_df[plot_df$change != "NS", , drop = FALSE]
  
  if (nrow(sig_df) > 0 && top_n_label > 0) {
    
    sig_df <- sig_df[
      order(sig_df$padj, -abs(sig_df$log2FoldChange)),
      ,
      drop = FALSE
    ]
    
    top_n_label <- min(top_n_label, nrow(sig_df))
    
    use_label <- if ("SYMBOL" %in% colnames(sig_df)) {
      ifelse(
        is.na(sig_df$SYMBOL) | sig_df$SYMBOL == "",
        sig_df$Peakid,
        sig_df$SYMBOL
      )
    } else {
      sig_df$Peakid
    }
    
    plot_df$label[
      match(sig_df$Peakid[1:top_n_label], plot_df$Peakid)
    ] <- use_label[1:top_n_label]
  }
  
  up_n <- sum(plot_df$change == "Up")
  down_n <- sum(plot_df$change == "Down")
  
  volcano_colors <- c(
    "Up" = "#E64B35",
    "Down" = "#4DBBD5",
    "NS" = "grey75"
  )
  
  p <- ggplot(
    plot_df,
    aes(x = log2FoldChange, y = neglog10padj)
  ) +
    geom_point(
      aes(color = change),
      size = 1.8,
      alpha = 0.8
    ) +
    scale_color_manual(values = volcano_colors) +
    geom_vline(
      xintercept = c(-lfc_cutoff, lfc_cutoff),
      linetype = "dashed",
      color = "grey40",
      linewidth = 0.5
    ) +
    geom_hline(
      yintercept = -log10(padj_cutoff),
      linetype = "dashed",
      color = "grey40",
      linewidth = 0.5
    ) +
    ggrepel::geom_text_repel(
      data = subset(plot_df, label != ""),
      aes(label = label),
      size = 3,
      max.overlaps = Inf,
      box.padding = 0.3,
      point.padding = 0.2,
      show.legend = FALSE
    ) +
    annotate(
      "text",
      x = min(plot_df$log2FoldChange, na.rm = TRUE),
      y = max(plot_df$neglog10padj, na.rm = TRUE),
      label = paste0("Down: ", down_n),
      hjust = 0,
      vjust = 1,
      size = 4.5,
      color = "#4DBBD5"
    ) +
    annotate(
      "text",
      x = max(plot_df$log2FoldChange, na.rm = TRUE),
      y = max(plot_df$neglog10padj, na.rm = TRUE),
      label = paste0("Up: ", up_n),
      hjust = 1,
      vjust = 1,
      size = 4.5,
      color = "#E64B35"
    ) +
    theme_bw(base_size = 14) +
    labs(
      title = paste0(comp_name, " ATAC volcano"),
      x = "log2FoldChange",
      y = "-log10(adjusted p-value)",
      color = NULL
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
  
  ggsave(
    filename = file.path(outdir_plot, paste0("ATAC_DAR_", comp_name, "_volcano.pdf")),
    plot = p,
    width = 7,
    height = 6
  )
  
  ggsave(
    filename = file.path(outdir_plot, paste0("ATAC_DAR_", comp_name, "_volcano.png")),
    plot = p,
    width = 7,
    height = 6,
    dpi = 300
  )
  
  return(p)
}

############################################################
## 20. Pairwise DAR function
############################################################

run_dar_pair_atac <- function(dds,
                              group1,
                              group2,
                              peak_annotation_all,
                              outdir_all,
                              outdir_sig,
                              outdir_plot,
                              out_name = NULL,
                              padj_cutoff = 0.05,
                              lfc_cutoff = 1) {
  
  if (is.null(out_name)) {
    out_name <- paste0(group1, "_vs_", group2)
  }
  
  group_levels_dds <- levels(colData(dds)$group)
  
  if (!(group1 %in% group_levels_dds)) {
    stop(paste("group1 not found in dds group levels:", group1))
  }
  
  if (!(group2 %in% group_levels_dds)) {
    stop(paste("group2 not found in dds group levels:", group2))
  }
  
  res <- results(
    dds,
    contrast = c("group", group1, group2)
  )
  
  res_df <- as.data.frame(res)
  res_df$Peakid <- rownames(res_df)
  
  keep_cols <- c(
    "Peakid",
    "baseMean",
    "log2FoldChange",
    "lfcSE",
    "stat",
    "pvalue",
    "padj"
  )
  
  keep_cols <- keep_cols[keep_cols %in% colnames(res_df)]
  res_basic <- res_df[, keep_cols, drop = FALSE]
  
  res_anno <- dplyr::left_join(
    res_basic,
    peak_annotation_all,
    by = "Peakid"
  )
  
  res_anno <- res_anno[
    order(res_anno$padj, na.last = TRUE),
    ,
    drop = FALSE
  ]
  
  safe_write_csv(
    res_anno,
    file = file.path(outdir_all, paste0("ATAC_DAR_", out_name, "_all.csv"))
  )
  
  res_sig <- subset(
    res_anno,
    !is.na(padj) &
      padj < padj_cutoff &
      abs(log2FoldChange) >= lfc_cutoff
  )
  
  if (nrow(res_sig) > 0) {
    
    res_sig$direction <- ifelse(
      res_sig$log2FoldChange >= lfc_cutoff,
      "Up",
      ifelse(res_sig$log2FoldChange <= -lfc_cutoff, "Down", NA)
    )
    
    res_sig <- res_sig[
      !is.na(res_sig$direction),
      ,
      drop = FALSE
    ]
    
    res_sig <- res_sig[
      order(abs(res_sig$log2FoldChange), decreasing = TRUE),
      ,
      drop = FALSE
    ]
    
  } else {
    
    res_sig$direction <- character(0)
  }
  
  safe_write_csv(
    res_sig,
    file = file.path(outdir_sig, paste0("ATAC_DAR_", out_name, "_sig.csv"))
  )
  
  export_sig_bed_atac(
    sig_df = res_sig,
    outdir_sig = outdir_sig,
    out_name = out_name
  )
  
  plot_volcano_atac(
    res_df = res_anno,
    comp_name = out_name,
    outdir_plot = outdir_plot,
    padj_cutoff = padj_cutoff,
    lfc_cutoff = lfc_cutoff,
    top_n_label = 0
  )
  
  summary_df <- data.frame(
    comparison = out_name,
    total_tested = length(unique(res_anno$Peakid[!is.na(res_anno$Peakid)])),
    total_sig = length(unique(res_sig$Peakid[!is.na(res_sig$Peakid)])),
    up = length(unique(res_sig$Peakid[res_sig$direction == "Up" & !is.na(res_sig$Peakid)])),
    down = length(unique(res_sig$Peakid[res_sig$direction == "Down" & !is.na(res_sig$Peakid)])),
    padj_cutoff = padj_cutoff,
    lfc_cutoff = lfc_cutoff,
    stringsAsFactors = FALSE
  )
  
  return(list(
    res = res,
    res_anno = res_anno,
    res_sig = res_sig,
    summary = summary_df
  ))
}

############################################################
## 21. Collect significant DAR annotation summary only
############################################################

collect_annotation_summary <- function(results_list,
                                       outdir,
                                       set_name,
                                       anno_levels,
                                       anno_merged_levels) {
  
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  fine_list <- list()
  merged_list <- list()
  promoter_check_list <- list()
  
  for (comp_name in names(results_list)) {
    
    res_one <- results_list[[comp_name]]
    
    sig_df <- res_one$res_sig
    
    fine_list[[comp_name]] <- summarise_annotation_distribution(
      df = sig_df,
      comparison = comp_name,
      result_type = "significant_DARs",
      anno_col = "anno_simple",
      anno_levels_use = anno_levels,
      summary_type = "fine_annotation"
    )
    
    merged_list[[comp_name]] <- summarise_annotation_distribution(
      df = sig_df,
      comparison = comp_name,
      result_type = "significant_DARs",
      anno_col = "anno_merged",
      anno_levels_use = anno_merged_levels,
      summary_type = "merged_annotation"
    )
    
    promoter_check_list[[comp_name]] <- check_promoter_merge(
      df = sig_df,
      comparison = comp_name,
      result_type = "significant_DARs"
    )
  }
  
  fine_df <- dplyr::bind_rows(fine_list)
  merged_df <- dplyr::bind_rows(merged_list)
  promoter_check_df <- dplyr::bind_rows(promoter_check_list)
  
  ## Double check: keep significant DARs only
  if (nrow(fine_df) > 0 && "result_type" %in% colnames(fine_df)) {
    fine_df <- fine_df %>%
      dplyr::filter(result_type == "significant_DARs")
  }
  
  if (nrow(merged_df) > 0 && "result_type" %in% colnames(merged_df)) {
    merged_df <- merged_df %>%
      dplyr::filter(result_type == "significant_DARs")
  }
  
  if (nrow(promoter_check_df) > 0 && "result_type" %in% colnames(promoter_check_df)) {
    promoter_check_df <- promoter_check_df %>%
      dplyr::filter(result_type == "significant_DARs")
  }
  
  safe_write_csv(
    fine_df,
    file = file.path(
      outdir,
      paste0("ATAC_", set_name, "_sig_DAR_fine_annotation_summary.csv")
    )
  )
  
  safe_write_csv(
    merged_df,
    file = file.path(
      outdir,
      paste0("ATAC_", set_name, "_sig_DAR_merged_annotation_summary.csv")
    )
  )
  
  safe_write_csv(
    promoter_check_df,
    file = file.path(
      outdir,
      paste0("ATAC_", set_name, "_sig_DAR_promoter_merge_check.csv")
    )
  )
  
  return(list(
    fine = fine_df,
    merged = merged_df,
    promoter_check = promoter_check_df
  ))
}

############################################################
## 22. Run unstimulated three-group pairwise DAR
############################################################

unstim_comparison_list_atac <- list(
  list(group1 = "GSK_un", group2 = "DMSO_un", out_name = "GSK_un_vs_DMSO_un"),
  list(group1 = "PBMC_un", group2 = "DMSO_un", out_name = "PBMC_un_vs_DMSO_un"),
  list(group1 = "GSK_un", group2 = "PBMC_un", out_name = "GSK_un_vs_PBMC_un")
)

unstim_dar_results_atac <- list()

for (comp in unstim_comparison_list_atac) {
  
  cat("Running unstim pairwise:", comp$out_name, "\n")
  
  unstim_dar_results_atac[[comp$out_name]] <- run_dar_pair_atac(
    dds = dds_atac,
    group1 = comp$group1,
    group2 = comp$group2,
    peak_annotation_all = peak_annotation_all,
    outdir_all = dar_unstim_all_dir,
    outdir_sig = dar_unstim_sig_dir,
    outdir_plot = volcano_unstim_dir,
    out_name = comp$out_name,
    padj_cutoff = 0.05,
    lfc_cutoff = 1
  )
}

unstim_summary_all_atac <- do.call(
  rbind,
  lapply(unstim_dar_results_atac, function(x) x$summary)
)

safe_write_csv(
  unstim_summary_all_atac,
  file = file.path(dar_unstim_dir, "ATAC_DAR_summary_all_unstim.csv")
)

saveRDS(
  list(
    pairwise = unstim_dar_results_atac,
    summary_all = unstim_summary_all_atac
  ),
  file = file.path(dar_unstim_dir, "ATAC_DAR_results_all_unstim.rds")
)

anno_unstim_summary <- collect_annotation_summary(
  results_list = unstim_dar_results_atac,
  outdir = anno_unstim_dir,
  set_name = "unstimulated_three_groups",
  anno_levels = anno_levels,
  anno_merged_levels = anno_merged_levels
)

############################################################
## 23. Run CD19-stimulated three-group pairwise DAR
############################################################

stim_comparison_list_atac <- list(
  list(group1 = "GSK_CD19", group2 = "DMSO_CD19", out_name = "GSK_CD19_vs_DMSO_CD19"),
  list(group1 = "PBMC_CD19", group2 = "DMSO_CD19", out_name = "PBMC_CD19_vs_DMSO_CD19"),
  list(group1 = "GSK_CD19", group2 = "PBMC_CD19", out_name = "GSK_CD19_vs_PBMC_CD19")
)

stim_dar_results_atac <- list()

for (comp in stim_comparison_list_atac) {
  
  cat("Running stim pairwise:", comp$out_name, "\n")
  
  stim_dar_results_atac[[comp$out_name]] <- run_dar_pair_atac(
    dds = dds_atac,
    group1 = comp$group1,
    group2 = comp$group2,
    peak_annotation_all = peak_annotation_all,
    outdir_all = dar_stim_all_dir,
    outdir_sig = dar_stim_sig_dir,
    outdir_plot = volcano_stim_dir,
    out_name = comp$out_name,
    padj_cutoff = 0.05,
    lfc_cutoff = 1
  )
}

stim_summary_all_atac <- do.call(
  rbind,
  lapply(stim_dar_results_atac, function(x) x$summary)
)

safe_write_csv(
  stim_summary_all_atac,
  file = file.path(dar_stim_dir, "ATAC_DAR_summary_all_stim.csv")
)

saveRDS(
  list(
    pairwise = stim_dar_results_atac,
    summary_all = stim_summary_all_atac
  ),
  file = file.path(dar_stim_dir, "ATAC_DAR_results_all_stim.rds")
)

anno_stim_summary <- collect_annotation_summary(
  results_list = stim_dar_results_atac,
  outdir = anno_stim_dir,
  set_name = "CD19_stimulated_three_groups",
  anno_levels = anno_levels,
  anno_merged_levels = anno_merged_levels
)

############################################################
## 24. Run stimulation response pairwise DAR
############################################################

stim_unstim_comparison_list_atac <- list(
  list(group1 = "DMSO_CD19", group2 = "DMSO_un", out_name = "DMSO_CD19_vs_DMSO_un"),
  list(group1 = "GSK_CD19", group2 = "GSK_un", out_name = "GSK_CD19_vs_GSK_un"),
  list(group1 = "PBMC_CD19", group2 = "PBMC_un", out_name = "PBMC_CD19_vs_PBMC_un")
)

stim_unstim_dar_results_atac <- list()

for (comp in stim_unstim_comparison_list_atac) {
  
  cat("Running stim vs unstim pairwise:", comp$out_name, "\n")
  
  stim_unstim_dar_results_atac[[comp$out_name]] <- run_dar_pair_atac(
    dds = dds_atac,
    group1 = comp$group1,
    group2 = comp$group2,
    peak_annotation_all = peak_annotation_all,
    outdir_all = dar_stim_unstim_all_dir,
    outdir_sig = dar_stim_unstim_sig_dir,
    outdir_plot = volcano_stim_unstim_dir,
    out_name = comp$out_name,
    padj_cutoff = 0.05,
    lfc_cutoff = 1
  )
}

stim_unstim_summary_all_atac <- do.call(
  rbind,
  lapply(stim_unstim_dar_results_atac, function(x) x$summary)
)

safe_write_csv(
  stim_unstim_summary_all_atac,
  file = file.path(dar_stim_unstim_dir, "ATAC_DAR_summary_all_stim_vs_unstim.csv")
)

saveRDS(
  list(
    stim_vs_unstim = stim_unstim_dar_results_atac,
    summary_all = stim_unstim_summary_all_atac
  ),
  file = file.path(dar_stim_unstim_dir, "ATAC_DAR_results_all_stim_vs_unstim.rds")
)

anno_stim_unstim_summary <- collect_annotation_summary(
  results_list = stim_unstim_dar_results_atac,
  outdir = anno_stim_unstim_dir,
  set_name = "stim_vs_unstim",
  anno_levels = anno_levels,
  anno_merged_levels = anno_merged_levels
)

############################################################
## 25. Save combined DAR and significant DAR annotation summaries
############################################################
## The final annotation distribution plots shown in Fig. 4f
## and Fig. 4j were generated in Origin software using the
## exported significant DAR merged annotation summary tables:
## 1. ATAC_unstimulated_three_groups_sig_DAR_merged_annotation_summary.csv
## 2. ATAC_CD19_stimulated_three_groups_sig_DAR_merged_annotation_summary.csv
############################################################

dar_summary_all_atac <- rbind(
  unstim_summary_all_atac,
  stim_summary_all_atac,
  stim_unstim_summary_all_atac
)

safe_write_csv(
  dar_summary_all_atac,
  file = file.path(dar_dir, "ATAC_DAR_summary_all.csv")
)

anno_sig_fine_all <- dplyr::bind_rows(
  anno_unstim_summary$fine,
  anno_stim_summary$fine,
  anno_stim_unstim_summary$fine
)

anno_sig_merged_all <- dplyr::bind_rows(
  anno_unstim_summary$merged,
  anno_stim_summary$merged,
  anno_stim_unstim_summary$merged
)

promoter_sig_check_all <- dplyr::bind_rows(
  anno_unstim_summary$promoter_check,
  anno_stim_summary$promoter_check,
  anno_stim_unstim_summary$promoter_check
)

## Double check: keep significant DARs only
if (nrow(anno_sig_fine_all) > 0 && "result_type" %in% colnames(anno_sig_fine_all)) {
  anno_sig_fine_all <- anno_sig_fine_all %>%
    dplyr::filter(result_type == "significant_DARs")
}

if (nrow(anno_sig_merged_all) > 0 && "result_type" %in% colnames(anno_sig_merged_all)) {
  anno_sig_merged_all <- anno_sig_merged_all %>%
    dplyr::filter(result_type == "significant_DARs")
}

if (nrow(promoter_sig_check_all) > 0 && "result_type" %in% colnames(promoter_sig_check_all)) {
  promoter_sig_check_all <- promoter_sig_check_all %>%
    dplyr::filter(result_type == "significant_DARs")
}

safe_write_csv(
  anno_sig_fine_all,
  file = file.path(
    anno_summary_dir,
    "ATAC_sig_DAR_fine_annotation_summary.csv"
  )
)

safe_write_csv(
  anno_sig_merged_all,
  file = file.path(
    anno_summary_dir,
    "ATAC_sig_DAR_merged_annotation_summary.csv"
  )
)

safe_write_csv(
  promoter_sig_check_all,
  file = file.path(
    anno_summary_dir,
    "ATAC_sig_DAR_promoter_merge_check.csv"
  )
)

saveRDS(
  list(
    dds_atac = dds_atac,
    meta_deseq = meta_deseq,
    peak_annotation_all = peak_annotation_all,
    unstim_pairwise = unstim_dar_results_atac,
    stim_pairwise = stim_dar_results_atac,
    stim_vs_unstim = stim_unstim_dar_results_atac,
    dar_summary_all = dar_summary_all_atac,
    sig_DAR_annotation_fine_all = anno_sig_fine_all,
    sig_DAR_annotation_merged_all = anno_sig_merged_all,
    sig_DAR_promoter_merge_check_all = promoter_sig_check_all
  ),
  file = file.path(
    dar_dir,
    "ATAC_DAR_results_all_with_sig_annotation_summary.rds"
  )
)

############################################################
## 26. KEGG and GO enrichment analysis for pairwise DARs
############################################################
## KEGG/GO enrichment is performed using geneId annotated
## from significant DAR-associated peaks.
############################################################

kegg_go_dir <- file.path(result_dir, "07_KEGG_GO")

enrich_unstim_dir <- file.path(kegg_go_dir, "unstimulated_three_groups")
enrich_stim_dir <- file.path(kegg_go_dir, "CD19_stimulated_three_groups")
enrich_stim_unstim_dir <- file.path(kegg_go_dir, "stim_vs_unstim")

dir.create(kegg_go_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(enrich_unstim_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(enrich_stim_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(enrich_stim_unstim_dir, recursive = TRUE, showWarnings = FALSE)

options(timeout = 1200)
options(download.file.method = "libcurl")

packages_bioc_enrich <- c(
  "clusterProfiler",
  "org.Hs.eg.db"
)

for (pkg in packages_bioc_enrich) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

############################################################
## 26.1 Helper functions for enrichment
############################################################

extract_entrez_from_atac_sig <- function(sig_df,
                                         direction_keep = c("Up", "Down")) {
  
  direction_keep <- match.arg(direction_keep)
  
  if (is.null(sig_df) || nrow(sig_df) == 0) {
    return(character(0))
  }
  
  if (!"direction" %in% colnames(sig_df)) {
    return(character(0))
  }
  
  if (!"geneId" %in% colnames(sig_df)) {
    return(character(0))
  }
  
  df_use <- sig_df %>%
    dplyr::filter(direction == direction_keep)
  
  if (nrow(df_use) == 0) {
    return(character(0))
  }
  
  gene_vec <- as.character(df_use$geneId)
  gene_vec <- gene_vec[!is.na(gene_vec) & gene_vec != ""]
  
  ## A peak may be assigned to multiple genes.
  gene_vec <- unlist(strsplit(gene_vec, "/"))
  gene_vec <- trimws(gene_vec)
  gene_vec <- gene_vec[gene_vec != ""]
  
  ## Remove possible decimal suffix, such as 105378267.0
  gene_vec <- sub("\\.0+$", "", gene_vec)
  
  ## Keep numeric Entrez IDs only.
  gene_vec <- gene_vec[grepl("[0-9]", gene_vec)]
  gene_vec <- gsub("[^0-9]", "", gene_vec)
  gene_vec <- gene_vec[gene_vec != ""]
  
  gene_vec <- unique(gene_vec)
  
  return(gene_vec)
}

parse_ratio_num_den <- function(x) {
  
  x <- as.character(x)
  x[is.na(x)] <- NA_character_
  
  num <- suppressWarnings(as.numeric(sub("/.*", "", x)))
  den <- suppressWarnings(as.numeric(sub(".*/", "", x)))
  
  return(list(num = num, den = den))
}

add_enrich_metrics <- function(df) {
  
  if (is.null(df) || nrow(df) == 0) {
    return(df)
  }
  
  need_cols <- c("GeneRatio", "BgRatio", "Count")
  
  for (cc in need_cols) {
    if (!cc %in% colnames(df)) {
      df[[cc]] <- NA
    }
  }
  
  gr <- parse_ratio_num_den(df$GeneRatio)
  br <- parse_ratio_num_den(df$BgRatio)
  
  gene_hits <- gr$num
  input_size <- gr$den
  bg_hits <- br$num
  bg_size <- br$den
  
  count_num <- suppressWarnings(as.numeric(df$Count))
  k_obs <- ifelse(!is.na(count_num), count_num, gene_hits)
  
  rich_factor <- ifelse(
    !is.na(k_obs) & !is.na(bg_hits) & bg_hits > 0,
    k_obs / bg_hits,
    NA_real_
  )
  
  fold_enrichment <- ifelse(
    !is.na(k_obs) & !is.na(input_size) & input_size > 0 &
      !is.na(bg_hits) & bg_hits > 0 &
      !is.na(bg_size) & bg_size > 0,
    (k_obs / input_size) / (bg_hits / bg_size),
    NA_real_
  )
  
  expected_k <- ifelse(
    !is.na(input_size) & input_size > 0 &
      !is.na(bg_hits) & bg_hits >= 0 &
      !is.na(bg_size) & bg_size > 0,
    input_size * (bg_hits / bg_size),
    NA_real_
  )
  
  z_score <- ifelse(
    !is.na(k_obs) & !is.na(expected_k) & expected_k > 0,
    (k_obs - expected_k) / sqrt(expected_k),
    NA_real_
  )
  
  df$RichFactor <- rich_factor
  df$FoldEnrichment <- fold_enrichment
  df$zScore <- z_score
  
  return(df)
}

safe_enrich_kegg_atac <- function(gene_ids) {
  
  if (length(gene_ids) == 0) {
    return(NULL)
  }
  
  kk <- tryCatch(
    enrichKEGG(
      gene = gene_ids,
      organism = "hsa",
      keyType = "ncbi-geneid",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff = 0.2
    ),
    error = function(e) {
      cat("KEGG enrichment failed:", conditionMessage(e), "\n")
      return(NULL)
    }
  )
  
  if (is.null(kk)) {
    return(NULL)
  }
  
  kk_df <- as.data.frame(kk)
  
  if (nrow(kk_df) == 0) {
    return(NULL)
  }
  
  return(kk)
}

safe_enrich_go_atac <- function(gene_ids,
                                ont_use = c("BP", "CC", "MF")) {
  
  ont_use <- match.arg(ont_use)
  
  if (length(gene_ids) == 0) {
    return(NULL)
  }
  
  ego <- tryCatch(
    enrichGO(
      gene = gene_ids,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = ont_use,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff = 0.2,
      readable = TRUE
    ),
    error = function(e) {
      cat("GO-", ont_use, " enrichment failed: ", conditionMessage(e), "\n", sep = "")
      return(NULL)
    }
  )
  
  if (is.null(ego)) {
    return(NULL)
  }
  
  ego_df <- as.data.frame(ego)
  
  if (nrow(ego_df) == 0) {
    return(NULL)
  }
  
  return(ego)
}

format_kegg_table_atac <- function(kegg_obj,
                                   comparison,
                                   direction_tag) {
  
  out <- data.frame(
    comparison = character(0),
    direction = character(0),
    ID = character(0),
    Description = character(0),
    GeneRatio = character(0),
    BgRatio = character(0),
    RichFactor = numeric(0),
    FoldEnrichment = numeric(0),
    zScore = numeric(0),
    pvalue = numeric(0),
    p.adjust = numeric(0),
    qvalue = numeric(0),
    geneID = character(0),
    Count = numeric(0),
    stringsAsFactors = FALSE
  )
  
  if (is.null(kegg_obj)) {
    return(out)
  }
  
  df <- as.data.frame(kegg_obj, stringsAsFactors = FALSE)
  
  if (nrow(df) == 0) {
    return(out)
  }
  
  df$comparison <- as.character(comparison)
  df$direction <- as.character(direction_tag)
  
  df <- add_enrich_metrics(df)
  
  std_cols <- c(
    "comparison",
    "direction",
    "ID",
    "Description",
    "GeneRatio",
    "BgRatio",
    "RichFactor",
    "FoldEnrichment",
    "zScore",
    "pvalue",
    "p.adjust",
    "qvalue",
    "geneID",
    "Count"
  )
  
  for (cc in std_cols) {
    if (!cc %in% colnames(df)) {
      if (cc %in% c(
        "RichFactor", "FoldEnrichment", "zScore",
        "pvalue", "p.adjust", "qvalue", "Count"
      )) {
        df[[cc]] <- NA_real_
      } else {
        df[[cc]] <- NA_character_
      }
    }
  }
  
  df <- df[, std_cols, drop = FALSE]
  
  df$RichFactor <- suppressWarnings(as.numeric(df$RichFactor))
  df$FoldEnrichment <- suppressWarnings(as.numeric(df$FoldEnrichment))
  df$zScore <- suppressWarnings(as.numeric(df$zScore))
  df$pvalue <- suppressWarnings(as.numeric(df$pvalue))
  df$p.adjust <- suppressWarnings(as.numeric(df$p.adjust))
  df$qvalue <- suppressWarnings(as.numeric(df$qvalue))
  df$Count <- suppressWarnings(as.numeric(df$Count))
  
  df <- df[order(df$p.adjust, df$pvalue), , drop = FALSE]
  rownames(df) <- NULL
  
  return(df)
}

format_go_table_atac <- function(go_obj,
                                 comparison,
                                 direction_tag,
                                 ont_tag) {
  
  out <- data.frame(
    comparison = character(0),
    direction = character(0),
    ONTOLOGY = character(0),
    ID = character(0),
    Description = character(0),
    GeneRatio = character(0),
    BgRatio = character(0),
    RichFactor = numeric(0),
    FoldEnrichment = numeric(0),
    zScore = numeric(0),
    pvalue = numeric(0),
    p.adjust = numeric(0),
    qvalue = numeric(0),
    geneID = character(0),
    Count = numeric(0),
    stringsAsFactors = FALSE
  )
  
  if (is.null(go_obj)) {
    return(out)
  }
  
  df <- as.data.frame(go_obj, stringsAsFactors = FALSE)
  
  if (nrow(df) == 0) {
    return(out)
  }
  
  df$comparison <- as.character(comparison)
  df$direction <- as.character(direction_tag)
  df$ONTOLOGY <- as.character(ont_tag)
  
  df <- add_enrich_metrics(df)
  
  std_cols <- c(
    "comparison",
    "direction",
    "ONTOLOGY",
    "ID",
    "Description",
    "GeneRatio",
    "BgRatio",
    "RichFactor",
    "FoldEnrichment",
    "zScore",
    "pvalue",
    "p.adjust",
    "qvalue",
    "geneID",
    "Count"
  )
  
  for (cc in std_cols) {
    if (!cc %in% colnames(df)) {
      if (cc %in% c(
        "RichFactor", "FoldEnrichment", "zScore",
        "pvalue", "p.adjust", "qvalue", "Count"
      )) {
        df[[cc]] <- NA_real_
      } else {
        df[[cc]] <- NA_character_
      }
    }
  }
  
  df <- df[, std_cols, drop = FALSE]
  
  df$RichFactor <- suppressWarnings(as.numeric(df$RichFactor))
  df$FoldEnrichment <- suppressWarnings(as.numeric(df$FoldEnrichment))
  df$zScore <- suppressWarnings(as.numeric(df$zScore))
  df$pvalue <- suppressWarnings(as.numeric(df$pvalue))
  df$p.adjust <- suppressWarnings(as.numeric(df$p.adjust))
  df$qvalue <- suppressWarnings(as.numeric(df$qvalue))
  df$Count <- suppressWarnings(as.numeric(df$Count))
  
  df <- df[
    order(df$ONTOLOGY, df$direction, df$p.adjust, df$pvalue),
    ,
    drop = FALSE
  ]
  
  rownames(df) <- NULL
  
  return(df)
}

############################################################
## 26.2 Run KEGG/GO for one pairwise comparison
############################################################

run_kegg_go_one_atac <- function(sig_df,
                                 comp_name,
                                 outdir_kegg,
                                 outdir_go) {
  
  cat("Running ATAC KEGG/GO enrichment:", comp_name, "\n")
  
  if (is.null(sig_df) || nrow(sig_df) == 0) {
    
    cat("  Skipped: empty significant DAR table.\n")
    
    empty_kegg <- format_kegg_table_atac(NULL, comp_name, "Up")
    empty_go <- format_go_table_atac(NULL, comp_name, "Up", "BP")
    
    safe_write_csv(
      empty_kegg,
      file = file.path(
        outdir_kegg,
        paste0("ATAC_", comp_name, "_KEGG_up_down_combined.csv")
      )
    )
    
    safe_write_csv(
      empty_go,
      file = file.path(
        outdir_go,
        paste0("ATAC_", comp_name, "_GO_BP_CC_MF_up_down_combined.csv")
      )
    )
    
    summary_df <- data.frame(
      comparison = comp_name,
      n_input_rows = 0,
      n_up_genes = 0,
      n_down_genes = 0,
      n_kegg_up = 0,
      n_kegg_down = 0,
      n_go_bp_up = 0,
      n_go_bp_down = 0,
      n_go_cc_up = 0,
      n_go_cc_down = 0,
      n_go_mf_up = 0,
      n_go_mf_down = 0,
      stringsAsFactors = FALSE
    )
    
    return(list(
      summary = summary_df,
      kegg = empty_kegg,
      go = empty_go
    ))
  }
  
  if (!"direction" %in% colnames(sig_df)) {
    
    cat("  Skipped: direction column not found.\n")
    
    empty_kegg <- format_kegg_table_atac(NULL, comp_name, "Up")
    empty_go <- format_go_table_atac(NULL, comp_name, "Up", "BP")
    
    summary_df <- data.frame(
      comparison = comp_name,
      n_input_rows = nrow(sig_df),
      n_up_genes = 0,
      n_down_genes = 0,
      n_kegg_up = 0,
      n_kegg_down = 0,
      n_go_bp_up = 0,
      n_go_bp_down = 0,
      n_go_cc_up = 0,
      n_go_cc_down = 0,
      n_go_mf_up = 0,
      n_go_mf_down = 0,
      stringsAsFactors = FALSE
    )
    
    return(list(
      summary = summary_df,
      kegg = empty_kegg,
      go = empty_go
    ))
  }
  
  sig_df$direction <- trimws(as.character(sig_df$direction))
  sig_df$direction[sig_df$direction %in% c("up", "UP")] <- "Up"
  sig_df$direction[sig_df$direction %in% c("down", "DOWN")] <- "Down"
  
  genes_up <- extract_entrez_from_atac_sig(
    sig_df = sig_df,
    direction_keep = "Up"
  )
  
  genes_down <- extract_entrez_from_atac_sig(
    sig_df = sig_df,
    direction_keep = "Down"
  )
  
  cat("  Input significant DAR rows:", nrow(sig_df), "\n")
  cat("  Up genes:", length(genes_up), "\n")
  cat("  Down genes:", length(genes_down), "\n")
  
  ##########################################################
  ## KEGG
  ##########################################################
  
  kk_up <- safe_enrich_kegg_atac(genes_up)
  kk_down <- safe_enrich_kegg_atac(genes_down)
  
  kegg_up_df <- format_kegg_table_atac(
    kegg_obj = kk_up,
    comparison = comp_name,
    direction_tag = "Up"
  )
  
  kegg_down_df <- format_kegg_table_atac(
    kegg_obj = kk_down,
    comparison = comp_name,
    direction_tag = "Down"
  )
  
  kegg_combined_df <- dplyr::bind_rows(
    kegg_down_df,
    kegg_up_df
  )
  
  safe_write_csv(
    kegg_combined_df,
    file = file.path(
      outdir_kegg,
      paste0("ATAC_", comp_name, "_KEGG_up_down_combined.csv")
    )
  )
  
  ##########################################################
  ## GO BP / CC / MF
  ##########################################################
  
  go_bp_up <- safe_enrich_go_atac(genes_up, ont_use = "BP")
  go_bp_down <- safe_enrich_go_atac(genes_down, ont_use = "BP")
  
  go_cc_up <- safe_enrich_go_atac(genes_up, ont_use = "CC")
  go_cc_down <- safe_enrich_go_atac(genes_down, ont_use = "CC")
  
  go_mf_up <- safe_enrich_go_atac(genes_up, ont_use = "MF")
  go_mf_down <- safe_enrich_go_atac(genes_down, ont_use = "MF")
  
  go_bp_up_df <- format_go_table_atac(go_bp_up, comp_name, "Up", "BP")
  go_bp_down_df <- format_go_table_atac(go_bp_down, comp_name, "Down", "BP")
  
  go_cc_up_df <- format_go_table_atac(go_cc_up, comp_name, "Up", "CC")
  go_cc_down_df <- format_go_table_atac(go_cc_down, comp_name, "Down", "CC")
  
  go_mf_up_df <- format_go_table_atac(go_mf_up, comp_name, "Up", "MF")
  go_mf_down_df <- format_go_table_atac(go_mf_down, comp_name, "Down", "MF")
  
  go_combined_df <- dplyr::bind_rows(
    go_bp_down_df,
    go_bp_up_df,
    go_cc_down_df,
    go_cc_up_df,
    go_mf_down_df,
    go_mf_up_df
  )
  
  safe_write_csv(
    go_combined_df,
    file = file.path(
      outdir_go,
      paste0("ATAC_", comp_name, "_GO_BP_CC_MF_up_down_combined.csv")
    )
  )
  
  ##########################################################
  ## Summary
  ##########################################################
  
  summary_df <- data.frame(
    comparison = comp_name,
    n_input_rows = nrow(sig_df),
    n_up_genes = length(genes_up),
    n_down_genes = length(genes_down),
    n_kegg_up = nrow(kegg_up_df),
    n_kegg_down = nrow(kegg_down_df),
    n_go_bp_up = nrow(go_bp_up_df),
    n_go_bp_down = nrow(go_bp_down_df),
    n_go_cc_up = nrow(go_cc_up_df),
    n_go_cc_down = nrow(go_cc_down_df),
    n_go_mf_up = nrow(go_mf_up_df),
    n_go_mf_down = nrow(go_mf_down_df),
    stringsAsFactors = FALSE
  )
  
  return(list(
    summary = summary_df,
    kegg = kegg_combined_df,
    go = go_combined_df
  ))
}

############################################################
## 26.3 Run KEGG/GO for one comparison class
############################################################

run_kegg_go_batch_atac <- function(results_list,
                                   outdir_base,
                                   batch_tag) {
  
  outdir_kegg <- file.path(outdir_base, "KEGG")
  outdir_go <- file.path(outdir_base, "GO")
  outdir_sum <- file.path(outdir_base, "summary")
  
  dir.create(outdir_kegg, recursive = TRUE, showWarnings = FALSE)
  dir.create(outdir_go, recursive = TRUE, showWarnings = FALSE)
  dir.create(outdir_sum, recursive = TRUE, showWarnings = FALSE)
  
  all_summary_list <- list()
  all_kegg_list <- list()
  all_go_list <- list()
  
  for (comp_name in names(results_list)) {
    
    sig_df <- results_list[[comp_name]]$res_sig
    
    res_one <- run_kegg_go_one_atac(
      sig_df = sig_df,
      comp_name = comp_name,
      outdir_kegg = outdir_kegg,
      outdir_go = outdir_go
    )
    
    all_summary_list[[comp_name]] <- res_one$summary
    all_kegg_list[[comp_name]] <- res_one$kegg
    all_go_list[[comp_name]] <- res_one$go
  }
  
  all_summary_df <- dplyr::bind_rows(all_summary_list)
  all_kegg_df <- dplyr::bind_rows(all_kegg_list)
  all_go_df <- dplyr::bind_rows(all_go_list)
  
  safe_write_csv(
    all_summary_df,
    file = file.path(
      outdir_sum,
      paste0("ATAC_", batch_tag, "_KEGG_GO_summary.csv")
    )
  )
  
  safe_write_csv(
    all_kegg_df,
    file = file.path(
      outdir_sum,
      paste0("ATAC_", batch_tag, "_ALL_KEGG_up_down_combined.csv")
    )
  )
  
  safe_write_csv(
    all_go_df,
    file = file.path(
      outdir_sum,
      paste0("ATAC_", batch_tag, "_ALL_GO_BP_CC_MF_up_down_combined.csv")
    )
  )
  
  saveRDS(
    list(
      summary = all_summary_df,
      kegg = all_kegg_df,
      go = all_go_df
    ),
    file = file.path(
      outdir_base,
      paste0("ATAC_", batch_tag, "_KEGG_GO_results.rds")
    )
  )
  
  return(list(
    summary = all_summary_df,
    kegg = all_kegg_df,
    go = all_go_df
  ))
}

############################################################
## 26.4 Run KEGG/GO enrichment for three comparison classes
############################################################

cat("Running KEGG/GO enrichment for unstimulated three-group pairwise DARs...\n")

enrich_unstim_results_atac <- run_kegg_go_batch_atac(
  results_list = unstim_dar_results_atac,
  outdir_base = enrich_unstim_dir,
  batch_tag = "unstimulated_three_groups"
)

cat("Running KEGG/GO enrichment for CD19-stimulated three-group pairwise DARs...\n")

enrich_stim_results_atac <- run_kegg_go_batch_atac(
  results_list = stim_dar_results_atac,
  outdir_base = enrich_stim_dir,
  batch_tag = "CD19_stimulated_three_groups"
)

cat("Running KEGG/GO enrichment for stimulation response pairwise DARs...\n")

enrich_stim_unstim_results_atac <- run_kegg_go_batch_atac(
  results_list = stim_unstim_dar_results_atac,
  outdir_base = enrich_stim_unstim_dir,
  batch_tag = "stim_vs_unstim"
)

############################################################
## 26.5 Save combined KEGG/GO enrichment summaries
############################################################

enrich_summary_all_atac <- dplyr::bind_rows(
  enrich_unstim_results_atac$summary,
  enrich_stim_results_atac$summary,
  enrich_stim_unstim_results_atac$summary
)

enrich_kegg_all_atac <- dplyr::bind_rows(
  enrich_unstim_results_atac$kegg,
  enrich_stim_results_atac$kegg,
  enrich_stim_unstim_results_atac$kegg
)

enrich_go_all_atac <- dplyr::bind_rows(
  enrich_unstim_results_atac$go,
  enrich_stim_results_atac$go,
  enrich_stim_unstim_results_atac$go
)

safe_write_csv(
  enrich_summary_all_atac,
  file = file.path(
    kegg_go_dir,
    "ATAC_ALL_pairwise_KEGG_GO_summary.csv"
  )
)

safe_write_csv(
  enrich_kegg_all_atac,
  file = file.path(
    kegg_go_dir,
    "ATAC_ALL_pairwise_KEGG_up_down_combined.csv"
  )
)

safe_write_csv(
  enrich_go_all_atac,
  file = file.path(
    kegg_go_dir,
    "ATAC_ALL_pairwise_GO_BP_CC_MF_up_down_combined.csv"
  )
)

saveRDS(
  list(
    unstimulated_three_groups = enrich_unstim_results_atac,
    CD19_stimulated_three_groups = enrich_stim_results_atac,
    stim_vs_unstim = enrich_stim_unstim_results_atac,
    summary_all = enrich_summary_all_atac,
    kegg_all = enrich_kegg_all_atac,
    go_all = enrich_go_all_atac
  ),
  file = file.path(
    kegg_go_dir,
    "ATAC_ALL_pairwise_KEGG_GO_results.rds"
  )
)

############################################################
## 27. Generate selected GO pathway bar plots
############################################################

selected_pathway_dir <- file.path(result_dir, "08_Selected_pathway")
selected_go_base_dir <- file.path(selected_pathway_dir, "GO_selected_bar")

selected_go_unstim_dir <- file.path(
  selected_go_base_dir,
  "unstimulated_three_groups"
)

selected_go_stim_dir <- file.path(
  selected_go_base_dir,
  "CD19_stimulated_three_groups"
)

dir.create(selected_go_unstim_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(selected_go_stim_dir, recursive = TRUE, showWarnings = FALSE)

padj_cutoff <- 0.05
sig_line <- -log10(padj_cutoff)

col_up <- "#4EA660"
col_down <- "#5292F7"

############################################################
## 27.1 Read GO enrichment table
############################################################

read_go_enrich_csv <- function(file) {
  
  if (!file.exists(file)) {
    message("File does not exist: ", file)
    return(NULL)
  }
  
  finfo <- file.info(file)
  
  if (is.na(finfo$size) || finfo$size <= 5) {
    message("Empty file, skipped: ", basename(file))
    return(NULL)
  }
  
  df <- tryCatch(
    data.table::fread(
      file,
      data.table = FALSE,
      check.names = FALSE
    ),
    error = function(e) {
      message("Failed to read file, skipped: ", basename(file))
      return(NULL)
    }
  )
  
  if (is.null(df) || nrow(df) == 0 || ncol(df) == 0) {
    message("No content in file, skipped: ", basename(file))
    return(NULL)
  }
  
  colnames(df) <- trimws(colnames(df))
  
  if ("ONTOLOGY" %in% colnames(df) && !"ontology" %in% colnames(df)) {
    df$ontology <- df$ONTOLOGY
  }
  
  if (!"ontology" %in% colnames(df)) {
    df$ontology <- NA_character_
  }
  
  need_cols <- c(
    "comparison",
    "direction",
    "Description",
    "p.adjust",
    "Count"
  )
  
  miss_cols <- setdiff(need_cols, colnames(df))
  
  if (length(miss_cols) > 0) {
    message(
      "Missing columns, skipped: ",
      basename(file),
      " -> ",
      paste(miss_cols, collapse = ", ")
    )
    return(NULL)
  }
  
  df$comparison <- trimws(as.character(df$comparison))
  df$direction <- trimws(as.character(df$direction))
  df$Description <- trimws(as.character(df$Description))
  df$ontology <- trimws(as.character(df$ontology))
  df$p.adjust <- suppressWarnings(as.numeric(df$p.adjust))
  df$Count <- suppressWarnings(as.numeric(df$Count))
  
  df$direction <- dplyr::case_when(
    tolower(df$direction) %in% c("up", "upregulated") ~ "Up",
    tolower(df$direction) %in% c("down", "downregulated") ~ "Down",
    TRUE ~ df$direction
  )
  
  df <- df %>%
    dplyr::filter(
      !is.na(comparison),
      !is.na(direction),
      !is.na(Description),
      !is.na(p.adjust),
      !is.na(Count)
    ) %>%
    dplyr::filter(direction %in% c("Up", "Down"))
  
  if (nrow(df) == 0) {
    message("No Up/Down GO terms, skipped: ", basename(file))
    return(NULL)
  }
  
  df$p.adjust[df$p.adjust <= 0 | is.na(df$p.adjust)] <- 1e-300
  df$neglog10_padj <- -log10(df$p.adjust)
  
  return(df)
}

############################################################
## 27.2 Extract comparison name from GO output file
############################################################

get_go_comp_name_from_file <- function(x) {
  
  x <- basename(x)
  x <- sub("\\.csv$", "", x, ignore.case = TRUE)
  x <- sub("^ATAC_", "", x)
  x <- sub("_GO_BP_CC_MF_up_down_combined$", "", x)
  
  return(x)
}

############################################################
## 27.3 Selected GO term lists
############################################################

selected_go_unstim <- list(
  
  "PBMC_un_vs_DMSO_un" = data.frame(
    Description = c(
      "lymphocyte differentiation",
      "small GTPase-mediated signal transduction",
      "alpha-beta T cell differentiation",
      "positive regulation of cell activation",
      "regulation of autophagy",
      "mitotic cell cycle phase transition",
      "regulation of CD8-positive, alpha-beta T cell activation",
      "chromosome segregation",
      "positive regulation of alpha-beta T cell activation",
      "positive regulation of cell-cell adhesion",
      "regulation of cell cycle phase transition",
      "lymphocyte proliferation",
      "alpha-beta T cell activation"
    ),
    ontology = c(
      "BP", "BP", "BP", "BP", "BP",
      "BP", "BP", "BP", "BP", "BP",
      "BP", "BP", "BP"
    ),
    stringsAsFactors = FALSE
  ),
  
  "GSK_un_vs_DMSO_un" = data.frame(
    Description = c(
      "lymphocyte differentiation",
      "small GTPase-mediated signal transduction",
      "alpha-beta T cell differentiation",
      "positive regulation of cell activation",
      "regulation of autophagy",
      "mitotic cell cycle phase transition",
      "regulation of CD8-positive, alpha-beta T cell activation",
      "chromosome segregation",
      "positive regulation of alpha-beta T cell activation",
      "positive regulation of cell-cell adhesion",
      "regulation of cell cycle phase transition",
      "lymphocyte proliferation",
      "alpha-beta T cell activation"
    ),
    ontology = c(
      "BP", "BP", "BP", "BP", "BP",
      "BP", "BP", "BP", "BP", "BP",
      "BP", "BP", "BP"
    ),
    stringsAsFactors = FALSE
  )
)

selected_go_stim <- list(
  
  "PBMC_CD19_vs_DMSO_CD19" = data.frame(
    Description = c(
      "positive regulation of cell adhesion",
      "alpha-beta T cell activation",
      "immune response-activating cell surface receptor signaling pathway",
      "regulation of T cell activation",
      "regulation of leukocyte cell-cell adhesion",
      "alpha-beta T cell activation involved in immune response",
      "alpha-beta T cell differentiation involved in immune response",
      "immune response-regulating cell surface receptor signaling pathway",
      "T-helper cell differentiation",
      "T cell differentiation involved in immune response",
      "alpha-beta T cell differentiation"
    ),
    ontology = c(
      "BP", "BP", "BP", "BP", "BP",
      "BP", "BP", "BP", "BP", "BP", "BP"
    ),
    stringsAsFactors = FALSE
  ),
  
  "GSK_CD19_vs_DMSO_CD19" = data.frame(
    Description = c(
      "positive regulation of cell adhesion",
      "alpha-beta T cell activation",
      "immune response-activating cell surface receptor signaling pathway",
      "regulation of T cell activation",
      "regulation of leukocyte cell-cell adhesion",
      "alpha-beta T cell activation involved in immune response",
      "alpha-beta T cell differentiation involved in immune response",
      "immune response-regulating cell surface receptor signaling pathway",
      "T-helper cell differentiation",
      "T cell differentiation involved in immune response",
      "alpha-beta T cell differentiation"
    ),
    ontology = c(
      "BP", "BP", "BP", "BP", "BP",
      "BP", "BP", "BP", "BP", "BP", "BP"
    ),
    stringsAsFactors = FALSE
  )
)

############################################################
## 27.4 Pick selected GO Up terms
############################################################

pick_selected_go_up <- function(df,
                                comp_name,
                                selected_go_list) {
  
  if (!comp_name %in% names(selected_go_list)) {
    message("selected_go_list does not include this comparison, skipped: ", comp_name)
    return(NULL)
  }
  
  sel_df <- selected_go_list[[comp_name]]
  
  df2 <- df %>%
    dplyr::filter(comparison == comp_name) %>%
    dplyr::filter(direction == "Up") %>%
    dplyr::filter(p.adjust < padj_cutoff) %>%
    dplyr::left_join(
      sel_df %>% dplyr::mutate(.selected = TRUE),
      by = c("Description", "ontology")
    ) %>%
    dplyr::filter(.selected %in% TRUE)
  
  if (nrow(df2) == 0) {
    return(NULL)
  }
  
  df2 <- df2 %>%
    dplyr::group_by(
      Description,
      ontology,
      direction
    ) %>%
    dplyr::slice_min(
      order_by = p.adjust,
      n = 1,
      with_ties = FALSE
    ) %>%
    dplyr::ungroup()
  
  df2 <- df2 %>%
    dplyr::mutate(
      term_show = paste0(Description, " (", ontology, ")")
    ) %>%
    dplyr::arrange(neglog10_padj, term_show)
  
  rownames(df2) <- NULL
  
  return(df2)
}

############################################################
## 27.5 Plot selected GO Up terms
############################################################

plot_selected_go_up_bar <- function(plot_df,
                                    comp_name,
                                    outdir) {
  
  if (is.null(plot_df) || nrow(plot_df) == 0) {
    message("Empty plot data, skipped: ", comp_name)
    return(NULL)
  }
  
  plot_df <- plot_df %>%
    dplyr::arrange(neglog10_padj, term_show) %>%
    dplyr::mutate(
      y = seq_len(n()),
      signed_x = neglog10_padj,
      count_label = as.character(Count)
    )
  
  xmax <- max(plot_df$neglog10_padj, na.rm = TRUE)
  xmax <- max(xmax, sig_line)
  
  label_offset <- xmax * 0.04
  name_offset <- xmax * 0.03
  
  plot_df <- plot_df %>%
    dplyr::mutate(
      count_x = signed_x + label_offset,
      path_x = -name_offset,
      path_hjust = 1
    )
  
  fig_h <- max(6, 0.45 * nrow(plot_df) + 2.5)
  
  p <- ggplot(plot_df) +
    geom_segment(
      aes(
        x = 0,
        xend = signed_x,
        y = y,
        yend = y
      ),
      linewidth = 14,
      lineend = "butt",
      color = col_up
    ) +
    geom_vline(
      xintercept = 0,
      color = "grey35",
      linewidth = 0.7
    ) +
    geom_vline(
      xintercept = sig_line,
      linetype = "dashed",
      color = "grey55",
      linewidth = 0.6
    ) +
    geom_text(
      aes(
        x = path_x,
        y = y,
        label = term_show,
        hjust = path_hjust
      ),
      size = 3.8,
      color = "black"
    ) +
    geom_text(
      aes(
        x = count_x,
        y = y,
        label = count_label
      ),
      size = 4.2,
      color = "black"
    ) +
    scale_x_continuous(
      limits = c(-(xmax * 1.35), xmax * 1.35),
      breaks = pretty(c(0, xmax), n = 6),
      labels = function(x) sprintf("%.1f", x)
    ) +
    scale_y_continuous(
      breaks = plot_df$y,
      labels = rep("", nrow(plot_df)),
      expand = expansion(mult = c(0.03, 0.03))
    ) +
    labs(
      title = paste0(comp_name, " GO enrichment (Up)"),
      x = expression(-log[10](padj)),
      y = NULL
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 15
      ),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      legend.position = "none"
    )
  
  out_prefix <- file.path(
    outdir,
    paste0(comp_name, "_selected_GO_up_bar")
  )
  
  ggsave(
    filename = paste0(out_prefix, ".pdf"),
    plot = p,
    width = 14,
    height = fig_h,
    units = "in"
  )
  
  ggsave(
    filename = paste0(out_prefix, ".png"),
    plot = p,
    width = 14,
    height = fig_h,
    units = "in",
    dpi = 300
  )
  
  safe_write_csv(
    plot_df %>%
      dplyr::select(
        comparison,
        direction,
        ontology,
        Description,
        term_show,
        Count,
        p.adjust,
        neglog10_padj,
        y
      ),
    file = paste0(out_prefix, "_plot_table.csv")
  )
  
  return(p)
}

############################################################
## 27.6 Run selected GO pathway plots
############################################################

run_selected_go_plot <- function(go_dir,
                                 selected_go_list,
                                 outdir) {
  
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  
  files_go <- list.files(
    go_dir,
    pattern = "^ATAC_.*_GO_BP_CC_MF_up_down_combined\\.csv$",
    full.names = TRUE
  )
  
  if (length(files_go) == 0) {
    warning("No GO csv files found in: ", go_dir)
    return(NULL)
  }
  
  for (ff in files_go) {
    
    comp_name <- get_go_comp_name_from_file(ff)
    
    if (!comp_name %in% names(selected_go_list)) {
      message("Not a selected comparison, skipped: ", comp_name)
      next
    }
    
    message("Processing selected GO plot: ", comp_name)
    
    go_df <- read_go_enrich_csv(ff)
    
    if (is.null(go_df) || nrow(go_df) == 0) {
      message("No GO result, skipped: ", comp_name)
      next
    }
    
    plot_df <- pick_selected_go_up(
      df = go_df,
      comp_name = comp_name,
      selected_go_list = selected_go_list
    )
    
    selected_terms <- selected_go_list[[comp_name]]
    
    if (!is.null(plot_df) && nrow(plot_df) > 0) {
      missing_terms <- setdiff(
        selected_terms$Description,
        plot_df$Description
      )
    } else {
      missing_terms <- selected_terms$Description
    }
    
    if (length(missing_terms) > 0) {
      message(
        "  Selected GO terms not found: ",
        paste(missing_terms, collapse = " ; ")
      )
    }
    
    if (is.null(plot_df) || nrow(plot_df) == 0) {
      message("No selected significant GO terms, skipped: ", comp_name)
      next
    }
    
    plot_selected_go_up_bar(
      plot_df = plot_df,
      comp_name = comp_name,
      outdir = outdir
    )
  }
  
  return(invisible(TRUE))
}

run_selected_go_plot(
  go_dir = file.path(enrich_unstim_dir, "GO"),
  selected_go_list = selected_go_unstim,
  outdir = selected_go_unstim_dir
)

run_selected_go_plot(
  go_dir = file.path(enrich_stim_dir, "GO"),
  selected_go_list = selected_go_stim,
  outdir = selected_go_stim_dir
)


############################################################
## 28. Generate IGV-like ATAC signal plots from bigWig files
############################################################

genome <- "hg38"

outdir_unstim <- file.path(igv_dir, "unstim")
outdir_stim   <- file.path(igv_dir, "stim")

dir.create(outdir_unstim, recursive = TRUE, showWarnings = FALSE)
dir.create(outdir_stim,   recursive = TRUE, showWarnings = FALSE)

cat("Checking bigWig files:\n")
print(file.exists(bw_files))

cat("\nChecking input gene Excel file:\n")
print(file.exists(igv_gene_file))

cat("\nChecking RefSeq Select file:\n")
print(file.exists(refseq_file))

if (!all(file.exists(bw_files))) {
  stop("Some bigWig files are missing. Please check bw_dir and bw_files.")
}

if (!file.exists(igv_gene_file)) {
  stop("Input gene Excel file does not exist: ", igv_gene_file)
}

if (!file.exists(refseq_file)) {
  stop("RefSeq Select file does not exist: ", refseq_file)
}

############################################################
## 28.1 Colors
############################################################

col_dmso <- "#BDBDBD"
col_gsk  <- "#4EA660"
col_pbmc <- "#5292F7"

col_refseq_select <- "#000"

############################################################
## 28.2 Read gene coordinate table
############################################################

gene_raw <- readxl::read_excel(igv_gene_file)

cat("\nExcel column names:\n")
print(colnames(gene_raw))

need_cols <- c(
  "Gene",
  "transcript_id",
  "chr",
  "graph start",
  "graph end",
  "range",
  "strand",
  "position",
  "ref_txStart",
  "ref_txEnd",
  "group"
)

miss_cols <- setdiff(need_cols, colnames(gene_raw))

if (length(miss_cols) > 0) {
  stop("Excel file is missing columns: ", paste(miss_cols, collapse = ", "))
}

gene_df <- gene_raw[, need_cols]
gene_df <- as.data.frame(gene_df)

gene_df$Gene <- as.character(gene_df$Gene)
gene_df$transcript_id <- as.character(gene_df$transcript_id)
gene_df$chr <- as.character(gene_df$chr)
gene_df$`graph start` <- as.numeric(gene_df$`graph start`)
gene_df$`graph end` <- as.numeric(gene_df$`graph end`)
gene_df$range <- as.character(gene_df$range)
gene_df$strand <- as.character(gene_df$strand)
gene_df$position <- as.character(gene_df$position)
gene_df$ref_txStart <- as.numeric(gene_df$ref_txStart)
gene_df$ref_txEnd <- as.numeric(gene_df$ref_txEnd)
gene_df$group <- as.character(gene_df$group)

gene_df$chr <- sub("^chr", "", gene_df$chr, ignore.case = TRUE)

gene_df <- gene_df[
  !is.na(gene_df$Gene) & gene_df$Gene != "" &
    !is.na(gene_df$chr) & gene_df$chr != "" &
    !is.na(gene_df$`graph start`) &
    !is.na(gene_df$`graph end`),
  ,
  drop = FALSE
]

cat("\nGene coordinate table preview:\n")
print(head(gene_df))

############################################################
## 28.3 Read RefSeq Select annotation
############################################################

ref_df <- read.delim(
  gzfile(refseq_file),
  sep = "\t",
  header = FALSE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

if (ncol(ref_df) != 16) {
  stop("RefSeq Select file should have 16 columns. Current column number: ", ncol(ref_df))
}

colnames(ref_df) <- c(
  "bin",
  "transcript_id",
  "chr",
  "strand",
  "txStart",
  "txEnd",
  "cdsStart",
  "cdsEnd",
  "exonCount",
  "exonStarts",
  "exonEnds",
  "score",
  "Gene",
  "cdsStartStat",
  "cdsEndStat",
  "exonFrames"
)

ref_df$transcript_id <- as.character(ref_df$transcript_id)
ref_df$Gene <- as.character(ref_df$Gene)
ref_df$chr <- sub("^chr", "", as.character(ref_df$chr), ignore.case = TRUE)
ref_df$strand <- as.character(ref_df$strand)
ref_df$txStart <- as.numeric(ref_df$txStart)
ref_df$txEnd <- as.numeric(ref_df$txEnd)
ref_df$exonCount <- as.numeric(ref_df$exonCount)
ref_df$exonStarts <- as.character(ref_df$exonStarts)
ref_df$exonEnds <- as.character(ref_df$exonEnds)

cat("\nRefSeq Select preview:\n")
print(head(ref_df[, c("Gene", "transcript_id", "chr", "strand", "txStart", "txEnd")]))

############################################################
## 28.4 Parse y-axis range
############################################################

parse_ylim <- function(x) {
  
  x <- gsub(" ", "", x)
  x <- gsub("–", "-", x)
  x <- gsub("—", "-", x)
  
  if (is.na(x) || x == "") {
    return(c(0, 30))
  }
  
  sp <- strsplit(x, "-", fixed = TRUE)[[1]]
  
  if (length(sp) != 2) {
    return(c(0, 30))
  }
  
  out <- suppressWarnings(as.numeric(sp))
  
  if (any(is.na(out))) {
    return(c(0, 30))
  }
  
  return(out)
}

ylim_mat <- t(sapply(gene_df$range, parse_ylim))
gene_df$ylim1 <- ylim_mat[, 1]
gene_df$ylim2 <- ylim_mat[, 2]

############################################################
## 28.5 Helper functions
############################################################

read_bw_region <- function(bw_file, roi) {
  
  bw_obj <- rtracklayer::BigWigFile(bw_file)
  
  bw_si <- tryCatch(
    GenomeInfoDb::seqinfo(bw_obj),
    error = function(e) NULL
  )
  
  if (is.null(bw_si)) {
    return(
      rtracklayer::import(
        bw_file,
        which = roi,
        format = "BigWig"
      )
    )
  }
  
  bw_seqlevels <- GenomeInfoDb::seqlevels(bw_si)
  
  roi_chr <- as.character(GenomeInfoDb::seqnames(roi))
  roi_chr_nochr <- sub("^chr", "", roi_chr, ignore.case = TRUE)
  roi_chr_chr <- paste0("chr", roi_chr_nochr)
  
  matched_chr <- ifelse(
    roi_chr_chr %in% bw_seqlevels,
    roi_chr_chr,
    ifelse(
      roi_chr_nochr %in% bw_seqlevels,
      roi_chr_nochr,
      roi_chr
    )
  )
  
  roi2 <- GenomicRanges::GRanges(
    seqnames = matched_chr,
    ranges = IRanges::ranges(roi)
  )
  
  common_seq <- intersect(
    GenomeInfoDb::seqlevels(roi2),
    GenomeInfoDb::seqlevels(bw_si)
  )
  
  if (length(common_seq) > 0) {
    GenomeInfoDb::seqlevels(
      roi2,
      pruning.mode = "coarse"
    ) <- common_seq
    
    GenomeInfoDb::seqinfo(roi2) <- bw_si[
      GenomeInfoDb::seqlevels(roi2)
    ]
  }
  
  gr <- rtracklayer::import(
    bw_obj,
    which = roi2,
    format = "BigWig"
  )
  
  return(gr)
}

fix_chr_name <- function(gr) {
  
  if (length(gr) == 0) {
    return(gr)
  }
  
  old_seq <- GenomeInfoDb::seqlevels(gr)
  new_seq <- ifelse(grepl("^chr", old_seq), old_seq, paste0("chr", old_seq))
  GenomeInfoDb::seqlevels(gr) <- new_seq
  
  return(gr)
}

parse_num_vec <- function(x) {
  
  if (is.na(x) || x == "") {
    return(numeric(0))
  }
  
  x <- trimws(as.character(x))
  x <- gsub(",$", "", x)
  x <- gsub("\\s+", "", x)
  
  if (x == "") {
    return(numeric(0))
  }
  
  out <- suppressWarnings(as.numeric(strsplit(x, ",", fixed = TRUE)[[1]]))
  out <- out[!is.na(out)]
  
  return(out)
}

make_bw_track2 <- function(gr,
                           track_name,
                           fill_col,
                           bw_ylim,
                           genome,
                           chr_now) {
  
  bw_ylim <- as.numeric(bw_ylim)
  
  if (length(bw_ylim) != 2 || any(is.na(bw_ylim))) {
    bw_ylim <- c(0, 30)
  }
  
  yticks_use <- unique(c(bw_ylim[1], bw_ylim[2]))
  
  track_now <- DataTrack(
    range = gr,
    data = if (length(gr) > 0) gr$score else numeric(0),
    genome = genome,
    chromosome = paste0("chr", chr_now),
    name = track_name,
    type = "histogram",
    ylim = bw_ylim,
    baseline = bw_ylim[1],
    col.histogram = fill_col,
    fill.histogram = fill_col,
    col.axis = "black",
    background.title = "white",
    col.title = "black",
    fontsize = 8,
    lwd = 0.3
  )
  
  displayPars(track_now) <- list(
    col.histogram = fill_col,
    fill.histogram = fill_col,
    col.axis = "black",
    cex.axis = 0.9,
    fontface = 1,
    lwd = 0.3,
    lwd.baseline = 0.6,
    baseline = bw_ylim[1],
    ylim = bw_ylim,
    background.title = "white",
    col.title = "black",
    fontsize = 8,
    yTicksAt = yticks_use,
    col.grid = NA
  )
  
  return(track_now)
}

make_refseq_track <- function(gene_symbol,
                              transcript_id_now,
                              chr_now,
                              from_now,
                              to_now,
                              genome,
                              ref_df) {
  
  sub_df <- ref_df[
    ref_df$transcript_id == transcript_id_now,
    ,
    drop = FALSE
  ]
  
  if (nrow(sub_df) == 0) {
    sub_df <- ref_df[
      ref_df$Gene == gene_symbol &
        ref_df$chr == chr_now,
      ,
      drop = FALSE
    ]
  }
  
  if (nrow(sub_df) == 0) {
    return(
      AnnotationTrack(
        start = from_now,
        end = to_now,
        chromosome = paste0("chr", chr_now),
        genome = genome,
        name = "RefSeq\nSelect",
        id = paste0(gene_symbol, "\n(not found)"),
        shape = "box",
        fill = col_refseq_select,
        col = col_refseq_select,
        background.title = "white",
        col.title = "black",
        fontsize = 7
      )
    )
  }
  
  sub_df <- sub_df[1, , drop = FALSE]
  
  exon_starts <- parse_num_vec(sub_df$exonStarts)
  exon_ends <- parse_num_vec(sub_df$exonEnds)
  
  n_exon <- min(length(exon_starts), length(exon_ends))
  
  if (n_exon == 0) {
    return(
      AnnotationTrack(
        start = sub_df$txStart + 1,
        end = sub_df$txEnd,
        chromosome = paste0("chr", chr_now),
        genome = genome,
        name = "RefSeq\nSelect",
        id = gene_symbol,
        shape = "box",
        fill = col_refseq_select,
        col = col_refseq_select,
        background.title = "white",
        col.title = "black",
        fontsize = 7
      )
    )
  }
  
  exon_starts <- exon_starts[seq_len(n_exon)] + 1
  exon_ends <- exon_ends[seq_len(n_exon)]
  
  exon_df <- data.frame(
    start = exon_starts,
    end = exon_ends,
    stringsAsFactors = FALSE
  )
  
  exon_df <- exon_df[
    !is.na(exon_df$start) &
      !is.na(exon_df$end) &
      exon_df$end >= exon_df$start - 1,
    ,
    drop = FALSE
  ]
  
  if (nrow(exon_df) == 0) {
    return(
      AnnotationTrack(
        start = sub_df$txStart + 1,
        end = sub_df$txEnd,
        chromosome = paste0("chr", chr_now),
        genome = genome,
        name = "RefSeq\nSelect",
        id = gene_symbol,
        shape = "box",
        fill = col_refseq_select,
        col = col_refseq_select,
        background.title = "white",
        col.title = "black",
        fontsize = 7
      )
    )
  }
  
  exon_gr <- GRanges(
    seqnames = paste0("chr", chr_now),
    ranges = IRanges(
      start = exon_df$start,
      end = exon_df$end
    ),
    strand = sub_df$strand
  )
  
  mcols(exon_gr)$transcript <- rep(gene_symbol, length(exon_gr))
  mcols(exon_gr)$gene <- rep(gene_symbol, length(exon_gr))
  mcols(exon_gr)$symbol <- rep(gene_symbol, length(exon_gr))
  mcols(exon_gr)$feature <- rep("exon", length(exon_gr))
  mcols(exon_gr)$exon <- seq_along(exon_gr)
  
  refTrack <- GeneRegionTrack(
    exon_gr,
    genome = genome,
    chromosome = paste0("chr", chr_now),
    start = from_now,
    end = to_now,
    name = "RefSeq\nSelect",
    transcriptAnnotation = "symbol",
    showId = TRUE,
    geneSymbol = TRUE,
    background.title = "white",
    col.title = "black",
    fontsize = 7,
    fill = col_refseq_select,
    col = col_refseq_select,
    stacking = "dense"
  )
  
  return(refTrack)
}

get_bw_set <- function(plot_group) {
  
  if (tolower(plot_group) == "unstim") {
    return(list(
      files = c(
        "DMSO-1" = bw_files["DMSO-1"],
        "GSK-1" = bw_files["GSK-1"],
        "PBMC-1" = bw_files["PBMC-1"]
      ),
      names = c("DMSO", "GSK", "PBMC"),
      cols = c(col_dmso, col_gsk, col_pbmc)
    ))
  }
  
  if (tolower(plot_group) == "stim") {
    return(list(
      files = c(
        "DMSO-CD19-1" = bw_files["DMSO-CD19-1"],
        "GSK-CD19-1" = bw_files["GSK-CD19-1"],
        "PBMC-CD19-1" = bw_files["PBMC-CD19-1"]
      ),
      names = c("DMSO_CD19", "GSK_CD19", "PBMC_CD19"),
      cols = c(col_dmso, col_gsk, col_pbmc)
    ))
  }
  
  stop("plot_group must be 'unstim' or 'stim'.")
}

############################################################
## 28.6 Build tracks and plot one gene
############################################################

build_track_list <- function(gene_now,
                             transcript_id_now,
                             chr_now,
                             from_now,
                             to_now,
                             bw_ylim,
                             plot_group) {
  
  roi <- GRanges(
    seqnames = chr_now,
    ranges = IRanges(
      start = from_now,
      end = to_now
    )
  )
  
  bw_set <- get_bw_set(plot_group)
  
  track_list <- list()
  
  for (i in seq_along(bw_set$files)) {
    
    gr_now <- read_bw_region(bw_set$files[i], roi)
    gr_now <- fix_chr_name(gr_now)
    
    track_now <- make_bw_track2(
      gr = gr_now,
      track_name = bw_set$names[i],
      fill_col = bw_set$cols[i],
      bw_ylim = bw_ylim,
      genome = genome,
      chr_now = chr_now
    )
    
    track_list[[length(track_list) + 1]] <- track_now
  }
  
  refTrack <- make_refseq_track(
    gene_symbol = gene_now,
    transcript_id_now = transcript_id_now,
    chr_now = chr_now,
    from_now = from_now,
    to_now = to_now,
    genome = genome,
    ref_df = ref_df
  )
  
  track_list[[length(track_list) + 1]] <- refTrack
  
  return(track_list)
}

plot_one_gene <- function(gene_row,
                          outdir,
                          plot_group) {
  
  gene_now <- gene_row$Gene
  transcript_id_now <- gene_row$transcript_id
  chr_now <- gene_row$chr
  from_now <- gene_row$`graph start`
  to_now <- gene_row$`graph end`
  bw_ylim <- c(gene_row$ylim1, gene_row$ylim2)
  
  if (is.na(from_now) || is.na(to_now) || to_now < from_now) {
    stop(
      paste0(
        "Coordinate error: ",
        gene_now,
        " | graph start = ",
        from_now,
        " | graph end = ",
        to_now
      )
    )
  }
  
  cat("\n=============================\n")
  cat("Plotting: ", plot_group, " | ", gene_now, " | ", transcript_id_now, "\n", sep = "")
  cat("=============================\n")
  
  track_list <- build_track_list(
    gene_now = gene_now,
    transcript_id_now = transcript_id_now,
    chr_now = chr_now,
    from_now = from_now,
    to_now = to_now,
    bw_ylim = bw_ylim,
    plot_group = plot_group
  )
  
  size_vec <- c(1, 1, 1, 0.7)
  
  pdf(
    file.path(outdir, paste0(gene_now, "_", plot_group, ".pdf")),
    width = 1.9,
    height = 2.6
  )
  
  plotTracks(
    track_list,
    from = from_now,
    to = to_now,
    chromosome = paste0("chr", chr_now),
    background.panel = "white",
    col.grid = NA,
    col.axis = "black",
    main = gene_now,
    sizes = size_vec,
    add53 = FALSE,
    add35 = FALSE
  )
  
  dev.off()
  
  png(
    file.path(outdir, paste0(gene_now, "_", plot_group, ".png")),
    width = 470,
    height = 650,
    res = 220
  )
  
  plotTracks(
    track_list,
    from = from_now,
    to = to_now,
    chromosome = paste0("chr", chr_now),
    background.panel = "white",
    col.grid = NA,
    col.axis = "black",
    main = gene_now,
    sizes = size_vec,
    add53 = FALSE,
    add35 = FALSE
  )
  
  dev.off()
  
  cat(gene_now, " plot finished.\n")
}

############################################################
## 28.7 Run IGV-like plots by group
############################################################

run_plot_by_group <- function(plot_group,
                              outdir_use) {
  
  cat("\n====================================================\n")
  cat("Processing group: ", plot_group, "\n", sep = "")
  cat("====================================================\n")
  
  gene_sub <- gene_df[
    tolower(gene_df$group) == tolower(plot_group),
    ,
    drop = FALSE
  ]
  
  if (nrow(gene_sub) == 0) {
    cat("No genes found for group = ", plot_group, ". Skipped.\n", sep = "")
    return(NULL)
  }
  
  cat("\nNumber of genes in this group: ", nrow(gene_sub), "\n", sep = "")
  print(gene_sub$Gene)
  
  for (i in seq_len(nrow(gene_sub))) {
    plot_one_gene(
      gene_row = gene_sub[i, ],
      outdir = outdir_use,
      plot_group = plot_group
    )
  }
  
  write.csv(
    gene_sub,
    file = file.path(outdir_use, paste0("gene_table_", plot_group, ".csv")),
    row.names = FALSE,
    quote = TRUE
  )
  
  cat("\nGroup ", plot_group, " finished.\n", sep = "")
}

run_plot_by_group(
  plot_group = "unstim",
  outdir_use = outdir_unstim
)

run_plot_by_group(
  plot_group = "stim",
  outdir_use = outdir_stim
)


############################################################
## Save session information
############################################################

writeLines(
  capture.output(sessionInfo()),
  con = file.path(result_dir, "sessionInfo.txt")
)

cat("Session information saved to:\n")
cat(file.path(result_dir, "sessionInfo.txt"), "\n")




















