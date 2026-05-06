############################################################
## Bulk RNA-seq full original analysis
## Fig. 4 / Supplementary Fig. 4 / Supplementary Fig. 5
##
## This script performs the full bulk RNA-seq downstream
## analysis using the complete featureCounts gene-level
## count matrix and sample metadata.
##
## Main analysis steps:
## 1. Read featureCounts gene-level count matrix
## 2. Clean sample names and match metadata
## 3. Extract gene annotation from GENCODE GTF
## 4. Generate basic QC summaries and mean ± SEM plots
## 5. Construct DESeq2 object and perform VST normalization
## 6. Generate six-group PCA plot
## 7. Generate Pearson correlation heatmap based on
##    top 1000 highly variable genes
## 8. Perform pairwise DESeq2 differential expression analysis
## 9. Perform three-group overall DESeq2 LRT analysis
##    for unstimulated and CD19-stimulated groups
## 10. Generate pairwise volcano plots
## 11. Generate three-group top1000 DEG heatmaps with
##     selected pathway-symbol gene labels
## 12. Perform KEGG and GO enrichment analysis
## 13. Generate selected KEGG pathway bar plots
##
## Note:
## The input paths should be modified according to the local
## environment before running this script.
############################################################

rm(list = ls())
gc()

############################################################
## 0. Load packages
############################################################

packages_cran <- c(
  "ggplot2",
  "dplyr",
  "tidyr",
  "pheatmap"
)

packages_bioc <- c(
  "DESeq2"
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
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(pheatmap)
})

############################################################
## 1. Input files and output directories
############################################################

## Please modify these paths according to your local environment.

raw_count_file <- "path/to/gene_counts.txt"
metadata_file  <- "path/to/metadata.csv"
gtf_file       <- "path/to/gencode.annotation.gtf.gz"

result_dir <- "path/to/output/bulkRNA_Fig4_SFig4_SFig5"

basic_dir   <- file.path(result_dir, "00_basic_data")
qc_dir      <- file.path(result_dir, "00_QC")
pca_dir     <- file.path(result_dir, "01_PCA")
pearson_dir <- file.path(result_dir, "02_Pearson_top1000")

dir.create(basic_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(qc_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pca_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(pearson_dir, recursive = TRUE, showWarnings = FALSE)

cat("Output directory:\n")
cat(result_dir, "\n")

############################################################
## 2. Read featureCounts raw count file
############################################################

raw_counts <- read.table(
  raw_count_file,
  header = TRUE,
  sep = "\t",
  comment.char = "#",
  stringsAsFactors = FALSE,
  check.names = FALSE
)

cat("raw_counts dimension:\n")
print(dim(raw_counts))

cat("raw_counts first columns:\n")
print(colnames(raw_counts)[1:min(10, ncol(raw_counts))])

############################################################
## 3. Extract count matrix
############################################################
## featureCounts first 6 columns are annotation columns:
## Geneid, Chr, Start, End, Strand, Length
## Count columns start from column 7.

count_df <- raw_counts[, c(1, 7:ncol(raw_counts))]
colnames(count_df)[1] <- "Geneid"

old_names <- colnames(count_df)[-1]

## Remove path before file name
name_step1 <- sub("^.*/", "", old_names)

## Remove STAR BAM suffix
## Keep sample suffix such as _L1 if present
new_names <- sub(
  "\\.Aligned.sortedByCoord\\.out\\.bam$",
  "",
  name_step1
)

colnames(count_df) <- c("Geneid", new_names)

rownames(count_df) <- count_df$Geneid

count_mat <- count_df[, -1, drop = FALSE]
count_mat <- as.matrix(count_mat)
mode(count_mat) <- "integer"

cat("count_mat dimension:\n")
print(dim(count_mat))

cat("count_mat sample names:\n")
print(colnames(count_mat))

############################################################
## 4. Read and check metadata
############################################################

meta <- read.csv(
  metadata_file,
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

need_cols <- c(
  "sample",
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

if (anyDuplicated(meta$sample) > 0) {
  stop("metadata$sample contains duplicated sample names.")
}

rownames(meta) <- meta$sample

missing_in_meta <- setdiff(colnames(count_mat), rownames(meta))
missing_in_count <- setdiff(rownames(meta), colnames(count_mat))

cat("Samples in count matrix but not metadata:\n")
print(missing_in_meta)

cat("Samples in metadata but not count matrix:\n")
print(missing_in_count)

if (length(missing_in_meta) > 0 || length(missing_in_count) > 0) {
  stop("Sample names do not match between count matrix and metadata.")
}

meta <- meta[colnames(count_mat), , drop = FALSE]

if (!all(rownames(meta) == colnames(count_mat))) {
  stop("Sample order mismatch after reordering metadata.")
}

############################################################
## 5. Set metadata factors
############################################################

meta$source <- factor(
  meta$source,
  levels = c("DMSO", "GSK", "PBMC")
)

meta$stim <- factor(
  meta$stim,
  levels = c("Unstim", "Stim")
)

meta$replicate <- factor(meta$replicate)

meta$group <- factor(
  meta$group,
  levels = c(
    "DMSO_un",
    "GSK_un",
    "PBMC_un",
    "DMSO_CD19",
    "GSK_CD19",
    "PBMC_CD19"
  )
)

cat("source × stim table:\n")
print(table(meta$source, meta$stim))

cat("group table:\n")
print(table(meta$group))

group_n <- table(meta$group)

if (!all(group_n == 4)) {
  warning("Not all groups contain four samples. Please check metadata.")
}

############################################################
## 6. Extract gene annotation from GTF
############################################################

extract_gtf_field <- function(attribute, key) {
  
  pattern1 <- paste0(key, " \"([^\"]+)\"")
  pattern2 <- paste0(key, " ([^;]+)")
  
  out <- rep(NA_character_, length(attribute))
  
  idx1 <- grepl(pattern1, attribute)
  out[idx1] <- sub(
    paste0(".*", pattern1, ".*"),
    "\\1",
    attribute[idx1]
  )
  
  idx2 <- is.na(out) & grepl(pattern2, attribute)
  out[idx2] <- sub(
    paste0(".*", pattern2, ".*"),
    "\\1",
    attribute[idx2]
  )
  
  out <- trimws(out)
  return(out)
}

gtf <- read.delim(
  gtf_file,
  header = FALSE,
  sep = "\t",
  comment.char = "#",
  stringsAsFactors = FALSE,
  quote = ""
)

colnames(gtf) <- c(
  "seqname",
  "source",
  "feature",
  "start",
  "end",
  "score",
  "strand",
  "frame",
  "attribute"
)

gtf_gene <- gtf[gtf$feature == "gene", , drop = FALSE]

gene_anno_gtf <- data.frame(
  Geneid_raw = extract_gtf_field(gtf_gene$attribute, "gene_id"),
  gene_name  = extract_gtf_field(gtf_gene$attribute, "gene_name"),
  gene_type  = extract_gtf_field(gtf_gene$attribute, "gene_type"),
  stringsAsFactors = FALSE
)

if (all(is.na(gene_anno_gtf$gene_type))) {
  gene_anno_gtf$gene_type <- extract_gtf_field(
    gtf_gene$attribute,
    "gene_biotype"
  )
}

gene_anno_gtf$Geneid <- sub("\\..*$", "", gene_anno_gtf$Geneid_raw)

gene_anno_gtf <- gene_anno_gtf[, c(
  "Geneid_raw",
  "Geneid",
  "gene_name",
  "gene_type"
)]

gene_anno_gtf <- unique(gene_anno_gtf)

cat("gene annotation dimension:\n")
print(dim(gene_anno_gtf))

cat("gene_type table:\n")
print(table(gene_anno_gtf$gene_type, useNA = "ifany"))

############################################################
## 7. Save basic data needed for downstream analysis
############################################################

gene_id_raw <- rownames(count_mat)
gene_id_noversion <- sub("\\..*$", "", gene_id_raw)

gene_info <- data.frame(
  Geneid_raw = gene_id_raw,
  Geneid = gene_id_noversion,
  stringsAsFactors = FALSE
)

saveRDS(
  count_mat,
  file = file.path(basic_dir, "gene_counts_clean_matrix.rds")
)

saveRDS(
  meta,
  file = file.path(basic_dir, "metadata_clean.rds")
)

write.csv(
  meta,
  file = file.path(basic_dir, "metadata_clean.csv"),
  row.names = FALSE,
  quote = FALSE
)

saveRDS(
  gene_info,
  file = file.path(basic_dir, "gene_id_version_map.rds")
)

write.csv(
  gene_info,
  file = file.path(basic_dir, "gene_id_version_map.csv"),
  row.names = FALSE,
  quote = FALSE
)

saveRDS(
  gene_anno_gtf,
  file = file.path(basic_dir, "gene_annotation_from_gtf.rds")
)

write.csv(
  gene_anno_gtf,
  file = file.path(basic_dir, "gene_annotation_from_gtf.csv"),
  row.names = FALSE,
  quote = FALSE
)

############################################################
## 8. Basic QC summary
############################################################

sample_library_size <- colSums(count_mat)
detected_genes <- colSums(count_mat > 0)

qc_basic <- data.frame(
  sample = colnames(count_mat),
  total_counts = sample_library_size,
  detected_genes = detected_genes,
  source = meta$source,
  stim = meta$stim,
  replicate = meta$replicate,
  group = meta$group,
  stringsAsFactors = FALSE
)

qc_basic$group_plot <- as.character(qc_basic$group)
qc_basic$group_plot[qc_basic$group_plot == "GSK_un"] <- "GSK761"
qc_basic$group_plot[qc_basic$group_plot == "GSK_CD19"] <- "GSK761_CD19"
qc_basic$group_plot[qc_basic$group_plot == "DMSO_un"] <- "DMSO"
qc_basic$group_plot[qc_basic$group_plot == "PBMC_un"] <- "PBMC"

qc_basic$group_plot <- factor(
  qc_basic$group_plot,
  levels = c(
    "DMSO",
    "GSK761",
    "PBMC",
    "DMSO_CD19",
    "GSK761_CD19",
    "PBMC_CD19"
  )
)

write.csv(
  qc_basic,
  file = file.path(qc_dir, "qc_basic_sample_summary.csv"),
  row.names = FALSE,
  quote = FALSE
)

############################################################
## 9. QC plots: mean ± SEM bar plot with sample dots
############################################################

qc_summary <- qc_basic %>%
  group_by(group_plot) %>%
  summarise(
    mean_total_counts = mean(total_counts),
    sem_total_counts = sd(total_counts) / sqrt(n()),
    mean_detected_genes = mean(detected_genes),
    sem_detected_genes = sd(detected_genes) / sqrt(n()),
    n = n(),
    .groups = "drop"
  )

write.csv(
  qc_summary,
  file = file.path(qc_dir, "qc_group_mean_sem_summary.csv"),
  row.names = FALSE,
  quote = FALSE
)

group_colors <- c(
  "DMSO" = "#7b95c6",
  "GSK761" = "#49c2d9",
  "PBMC" = "#a1d8e8",
  "DMSO_CD19" = "#67a583",
  "GSK761_CD19" = "#a2c986",
  "PBMC_CD19" = "#d0e2c0"
)

p_total_counts <- ggplot() +
  geom_col(
    data = qc_summary,
    aes(x = group_plot, y = mean_total_counts, fill = group_plot),
    color = "black",
    width = 0.7
  ) +
  geom_errorbar(
    data = qc_summary,
    aes(
      x = group_plot,
      ymin = mean_total_counts - sem_total_counts,
      ymax = mean_total_counts + sem_total_counts
    ),
    width = 0.2,
    linewidth = 0.5
  ) +
  geom_jitter(
    data = qc_basic,
    aes(x = group_plot, y = total_counts),
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
  ggtitle("RNA-seq")

ggsave(
  filename = file.path(qc_dir, "QC_total_counts_mean_SEM.png"),
  plot = p_total_counts,
  width = 5.5,
  height = 4.5,
  dpi = 300
)

ggsave(
  filename = file.path(qc_dir, "QC_total_counts_mean_SEM.pdf"),
  plot = p_total_counts,
  width = 5.5,
  height = 4.5
)



############################################################
## 10. Build DESeq2 object and VST
############################################################

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = meta,
  design = ~ group
)

cat("Original dds dimension:\n")
print(dim(dds))

keep <- rowSums(counts(dds) >= 10) >= 2

cat("Genes retained after filtering:\n")
print(sum(keep))

cat("Genes removed after filtering:\n")
print(sum(!keep))

dds <- dds[keep, ]

cat("Filtered dds dimension:\n")
print(dim(dds))

vsd <- vst(dds, blind = TRUE)
vsd_mat <- assay(vsd)

saveRDS(
  dds,
  file = file.path(basic_dir, "dds_prefiltered.rds")
)

saveRDS(
  vsd,
  file = file.path(basic_dir, "vsd.rds")
)

############################################################
## 11. Six-group PCA
############################################################

pca_data <- plotPCA(
  vsd,
  intgroup = c("source", "stim", "group"),
  returnData = TRUE
)

percent_var <- round(100 * attr(pca_data, "percentVar"))

pca_data$source_plot <- as.character(pca_data$source)
pca_data$source_plot[pca_data$source_plot == "GSK"] <- "GSK761"

pca_data$source_plot <- factor(
  pca_data$source_plot,
  levels = c("DMSO", "GSK761", "PBMC")
)

pca_data$stim <- factor(
  pca_data$stim,
  levels = c("Unstim", "Stim")
)

source_colors <- c(
  "DMSO" = "black",
  "GSK761" = "red",
  "PBMC" = "blue"
)

stim_shapes <- c(
  "Unstim" = 16,
  "Stim" = 17
)

p_pca <- ggplot(
  pca_data,
  aes(PC1, PC2, color = source_plot, shape = stim)
) +
  geom_point(size = 2.5) +
  scale_color_manual(values = source_colors) +
  scale_shape_manual(values = stim_shapes) +
  xlab(paste0("PC1: ", percent_var[1], "% variance")) +
  ylab(paste0("PC2: ", percent_var[2], "% variance")) +
  theme_bw(base_size = 14) +
  ggtitle("RNA-seq")

ggsave(
  filename = file.path(pca_dir, "PCA_six_groups.png"),
  plot = p_pca,
  width = 6,
  height = 4,
  dpi = 300
)

ggsave(
  filename = file.path(pca_dir, "PCA_six_groups.pdf"),
  plot = p_pca,
  width = 6,
  height = 4
)

write.csv(
  pca_data,
  file = file.path(pca_dir, "PCA_six_groups_data.csv"),
  row.names = TRUE,
  quote = FALSE
)

############################################################
## 12. Top1000 highly variable gene Pearson correlation heatmap
############################################################

gene_var <- apply(vsd_mat, 1, var)

top1000_genes <- names(sort(gene_var, decreasing = TRUE))[
  1:min(1000, length(gene_var))
]

vsd_top1000 <- vsd_mat[top1000_genes, , drop = FALSE]

group_levels <- c(
  "DMSO_un",
  "GSK_un",
  "PBMC_un",
  "DMSO_CD19",
  "GSK_CD19",
  "PBMC_CD19"
)

group_mean_mat <- sapply(group_levels, function(g) {
  samples_g <- rownames(meta)[meta$group == g]
  rowMeans(vsd_top1000[, samples_g, drop = FALSE])
})

colnames(group_mean_mat) <- c(
  "DMSO",
  "GSK761",
  "PBMC",
  "DMSO_CD19",
  "GSK761_CD19",
  "PBMC_CD19"
)

pearson_cor <- cor(
  group_mean_mat,
  method = "pearson",
  use = "pairwise.complete.obs"
)

write.csv(
  group_mean_mat,
  file.path(pearson_dir, "top1000_group_mean_expression.csv"),
  quote = FALSE
)

write.csv(
  pearson_cor,
  file.path(pearson_dir, "top1000_group_mean_Pearson_correlation.csv"),
  quote = FALSE
)

png(
  filename = file.path(pearson_dir, "top1000_group_mean_Pearson_heatmap.png"),
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
  main = "RNA-seq"
)

dev.off()

pdf(
  file = file.path(pearson_dir, "top1000_group_mean_Pearson_heatmap.pdf"),
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
  main = "RNA-seq"
)

dev.off()

############################################################
## 13. DESeq2 pairwise DEG and three-group overall DEG
############################################################

deg_dir <- file.path(result_dir, "03_DESeq2_DEG")
volcano_dir <- file.path(result_dir, "04_pairwise_Volcano")
three_group_heatmap_dir <- file.path(result_dir, "05_Three_group_DEG_heatmap")

dir.create(deg_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(volcano_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(three_group_heatmap_dir, recursive = TRUE, showWarnings = FALSE)

deg_unstim_dir <- file.path(deg_dir, "unstimulated_three_groups")
deg_stim_dir <- file.path(deg_dir, "CD19_stimulated_three_groups")
deg_stim_unstim_dir <- file.path(deg_dir, "stim_vs_unstim")

dir.create(deg_unstim_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(deg_stim_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(deg_stim_unstim_dir, recursive = TRUE, showWarnings = FALSE)

deg_unstim_all_dir <- file.path(deg_unstim_dir, "DEG_all_bulk")
deg_unstim_sig_dir <- file.path(deg_unstim_dir, "DEG_sig_bulk")

deg_stim_all_dir <- file.path(deg_stim_dir, "DEG_all_bulk")
deg_stim_sig_dir <- file.path(deg_stim_dir, "DEG_sig_bulk")

deg_stim_unstim_all_dir <- file.path(deg_stim_unstim_dir, "DEG_all_bulk")
deg_stim_unstim_sig_dir <- file.path(deg_stim_unstim_dir, "DEG_sig_bulk")

dir.create(deg_unstim_all_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(deg_unstim_sig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(deg_stim_all_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(deg_stim_sig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(deg_stim_unstim_all_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(deg_stim_unstim_sig_dir, recursive = TRUE, showWarnings = FALSE)

############################################################
## 13.1 Load additional packages
############################################################

packages_bioc_extra <- c(
  "ComplexHeatmap",
  "circlize"
)

for (pkg in packages_bioc_extra) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(grDevices)
})

############################################################
## 13.2 Run DESeq2 model
############################################################

cat("Running DESeq2 model for pairwise DEG analysis...\n")

dds <- DESeq(dds)

saveRDS(
  dds,
  file = file.path(basic_dir, "dds_DESeq_done.rds")
)

writeLines(
  text = resultsNames(dds),
  con = file.path(basic_dir, "resultsNames_bulk.txt")
)

writeLines(
  text = levels(meta$group),
  con = file.path(basic_dir, "group_levels_bulk.txt")
)

write.csv(
  as.data.frame(table(meta$group)),
  file = file.path(basic_dir, "group_table_bulk.csv"),
  row.names = FALSE,
  quote = FALSE
)

cat("DESeq2 model finished.\n")
print(resultsNames(dds))

############################################################
## 13.3 Pairwise DEG function
############################################################

run_deg_pair_bulk <- function(dds,
                              group1,
                              group2,
                              gene_anno_gtf,
                              outdir_all,
                              outdir_sig,
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
  res_df$Geneid_raw <- rownames(res_df)
  res_df$Geneid <- sub("\\..*$", "", res_df$Geneid_raw)
  
  res_anno <- merge(
    res_df,
    gene_anno_gtf,
    by = c("Geneid_raw", "Geneid"),
    all.x = TRUE,
    sort = FALSE
  )
  
  keep_cols <- c(
    "Geneid_raw", "Geneid", "gene_name", "gene_type",
    "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"
  )
  keep_cols <- keep_cols[keep_cols %in% colnames(res_anno)]
  res_anno <- res_anno[, keep_cols, drop = FALSE]
  
  res_anno <- res_anno[order(res_anno$padj, na.last = TRUE), , drop = FALSE]
  
  write.csv(
    res_anno,
    file = file.path(outdir_all, paste0("DEG_", out_name, "_all_bulk.csv")),
    row.names = FALSE,
    quote = FALSE
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
    
    res_sig <- res_sig[!is.na(res_sig$direction), , drop = FALSE]
    res_sig <- res_sig[order(abs(res_sig$log2FoldChange), decreasing = TRUE), , drop = FALSE]
  } else {
    res_sig$direction <- character(0)
  }
  
  keep_cols_sig <- c(
    "Geneid_raw", "Geneid", "gene_name", "gene_type", "direction",
    "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"
  )
  keep_cols_sig <- keep_cols_sig[keep_cols_sig %in% colnames(res_sig)]
  res_sig <- res_sig[, keep_cols_sig, drop = FALSE]
  
  write.csv(
    res_sig,
    file = file.path(outdir_sig, paste0("DEG_", out_name, "_sig_bulk.csv")),
    row.names = FALSE,
    quote = FALSE
  )
  
  summary_df <- data.frame(
    comparison = out_name,
    total_tested = length(unique(res_anno$Geneid[!is.na(res_anno$Geneid)])),
    total_sig = length(unique(res_sig$Geneid[!is.na(res_sig$Geneid)])),
    up = length(unique(res_sig$Geneid[res_sig$direction == "Up" & !is.na(res_sig$Geneid)])),
    down = length(unique(res_sig$Geneid[res_sig$direction == "Down" & !is.na(res_sig$Geneid)])),
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
## 13.4 Three-group overall DEG function using LRT
############################################################

run_three_group_deg_bulk <- function(dds,
                                     meta,
                                     groups_keep,
                                     out_prefix,
                                     gene_anno_gtf,
                                     outdir_all,
                                     outdir_sig,
                                     padj_cutoff = 0.05) {
  
  keep_samples <- rownames(meta)[meta$group %in% groups_keep]
  
  if (length(keep_samples) == 0) {
    stop("groups_keep did not match any samples.")
  }
  
  dds_sub <- dds[, keep_samples]
  meta_sub <- meta[keep_samples, , drop = FALSE]
  
  meta_sub$group <- factor(
    meta_sub$group,
    levels = groups_keep
  )
  
  colData(dds_sub)$group <- meta_sub$group
  
  design(dds_sub) <- ~ group
  
  dds_sub <- DESeq(
    dds_sub,
    test = "LRT",
    reduced = ~ 1
  )
  
  res <- results(dds_sub)
  res_df <- as.data.frame(res)
  
  res_df$Geneid_raw <- rownames(res_df)
  res_df$Geneid <- sub("\\..*$", "", res_df$Geneid_raw)
  
  res_anno <- merge(
    res_df,
    gene_anno_gtf,
    by = c("Geneid_raw", "Geneid"),
    all.x = TRUE,
    sort = FALSE
  )
  
  keep_cols <- c(
    "Geneid_raw", "Geneid", "gene_name", "gene_type",
    "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"
  )
  keep_cols <- keep_cols[keep_cols %in% colnames(res_anno)]
  res_anno <- res_anno[, keep_cols, drop = FALSE]
  
  res_anno <- res_anno[order(res_anno$padj, na.last = TRUE), , drop = FALSE]
  
  write.csv(
    res_anno,
    file = file.path(outdir_all, paste0(out_prefix, "_all_bulk.csv")),
    row.names = FALSE,
    quote = FALSE
  )
  
  res_sig <- subset(
    res_anno,
    !is.na(padj) & padj < padj_cutoff
  )
  
  if (nrow(res_sig) > 0 && "log2FoldChange" %in% colnames(res_sig)) {
    res_sig <- res_sig[order(abs(res_sig$log2FoldChange), decreasing = TRUE), , drop = FALSE]
  }
  
  write.csv(
    res_sig,
    file = file.path(outdir_sig, paste0(out_prefix, "_sig_bulk.csv")),
    row.names = FALSE,
    quote = FALSE
  )
  
  summary_df <- data.frame(
    comparison = out_prefix,
    total_tested = length(unique(res_anno$Geneid[!is.na(res_anno$Geneid)])),
    total_sig = length(unique(res_sig$Geneid[!is.na(res_sig$Geneid)])),
    up = NA,
    down = NA,
    padj_cutoff = padj_cutoff,
    lfc_cutoff = NA,
    stringsAsFactors = FALSE
  )
  
  return(list(
    dds_sub = dds_sub,
    res = res,
    res_anno = res_anno,
    res_sig = res_sig,
    summary = summary_df
  ))
}

############################################################
## 13.5 Run unstimulated three-group pairwise + overall
############################################################

unstim_comparison_list_bulk <- list(
  list(group1 = "GSK_un",  group2 = "DMSO_un", out_name = "GSK761_un_vs_DMSO_un"),
  list(group1 = "PBMC_un", group2 = "DMSO_un", out_name = "PBMC_un_vs_DMSO_un"),
  list(group1 = "GSK_un",  group2 = "PBMC_un", out_name = "GSK761_un_vs_PBMC_un")
)

unstim_deg_results_bulk <- list()

for (comp in unstim_comparison_list_bulk) {
  
  cat("Running unstim pairwise:", comp$out_name, "\n")
  
  unstim_deg_results_bulk[[comp$out_name]] <- run_deg_pair_bulk(
    dds = dds,
    group1 = comp$group1,
    group2 = comp$group2,
    gene_anno_gtf = gene_anno_gtf,
    outdir_all = deg_unstim_all_dir,
    outdir_sig = deg_unstim_sig_dir,
    out_name = comp$out_name,
    padj_cutoff = 0.05,
    lfc_cutoff = 1
  )
}

cat("Running unstim 3-group overall...\n")

res_unstim_3group_bulk <- run_three_group_deg_bulk(
  dds = dds,
  meta = meta,
  groups_keep = c("DMSO_un", "GSK_un", "PBMC_un"),
  out_prefix = "DEG_unstim_3group_overall_bulk",
  gene_anno_gtf = gene_anno_gtf,
  outdir_all = deg_unstim_all_dir,
  outdir_sig = deg_unstim_sig_dir,
  padj_cutoff = 0.05
)

unstim_summary_all_bulk <- do.call(
  rbind,
  c(
    lapply(unstim_deg_results_bulk, function(x) x$summary),
    list(res_unstim_3group_bulk$summary)
  )
)

write.csv(
  unstim_summary_all_bulk,
  file = file.path(deg_unstim_dir, "DEG_summary_all_unstim_bulk.csv"),
  row.names = FALSE,
  quote = FALSE
)

saveRDS(
  list(
    pairwise = unstim_deg_results_bulk,
    overall_3group = res_unstim_3group_bulk,
    summary_all = unstim_summary_all_bulk
  ),
  file = file.path(deg_unstim_dir, "DEG_results_all_unstim_bulk.rds")
)

############################################################
## 13.6 Run CD19-stimulated three-group pairwise + overall
############################################################

stim_comparison_list_bulk <- list(
  list(group1 = "GSK_CD19",  group2 = "DMSO_CD19", out_name = "GSK761_CD19_vs_DMSO_CD19"),
  list(group1 = "PBMC_CD19", group2 = "DMSO_CD19", out_name = "PBMC_CD19_vs_DMSO_CD19"),
  list(group1 = "GSK_CD19",  group2 = "PBMC_CD19", out_name = "GSK761_CD19_vs_PBMC_CD19")
)

stim_deg_results_bulk <- list()

for (comp in stim_comparison_list_bulk) {
  
  cat("Running stim pairwise:", comp$out_name, "\n")
  
  stim_deg_results_bulk[[comp$out_name]] <- run_deg_pair_bulk(
    dds = dds,
    group1 = comp$group1,
    group2 = comp$group2,
    gene_anno_gtf = gene_anno_gtf,
    outdir_all = deg_stim_all_dir,
    outdir_sig = deg_stim_sig_dir,
    out_name = comp$out_name,
    padj_cutoff = 0.05,
    lfc_cutoff = 1
  )
}

cat("Running stim 3-group overall...\n")

res_stim_3group_bulk <- run_three_group_deg_bulk(
  dds = dds,
  meta = meta,
  groups_keep = c("DMSO_CD19", "GSK_CD19", "PBMC_CD19"),
  out_prefix = "DEG_stim_3group_overall_bulk",
  gene_anno_gtf = gene_anno_gtf,
  outdir_all = deg_stim_all_dir,
  outdir_sig = deg_stim_sig_dir,
  padj_cutoff = 0.05
)

stim_summary_all_bulk <- do.call(
  rbind,
  c(
    lapply(stim_deg_results_bulk, function(x) x$summary),
    list(res_stim_3group_bulk$summary)
  )
)

write.csv(
  stim_summary_all_bulk,
  file = file.path(deg_stim_dir, "DEG_summary_all_stim_bulk.csv"),
  row.names = FALSE,
  quote = FALSE
)

saveRDS(
  list(
    pairwise = stim_deg_results_bulk,
    overall_3group = res_stim_3group_bulk,
    summary_all = stim_summary_all_bulk
  ),
  file = file.path(deg_stim_dir, "DEG_results_all_stim_bulk.rds")
)

############################################################
## 13.7 Run stim vs unstim pairwise comparisons
############################################################

stim_unstim_comparison_list_bulk <- list(
  list(group1 = "DMSO_CD19", group2 = "DMSO_un", out_name = "DMSO_CD19_vs_DMSO_un"),
  list(group1 = "GSK_CD19",  group2 = "GSK_un",  out_name = "GSK761_CD19_vs_GSK761_un"),
  list(group1 = "PBMC_CD19", group2 = "PBMC_un", out_name = "PBMC_CD19_vs_PBMC_un")
)

stim_unstim_deg_results_bulk <- list()

for (comp in stim_unstim_comparison_list_bulk) {
  
  cat("Running stim vs unstim:", comp$out_name, "\n")
  
  stim_unstim_deg_results_bulk[[comp$out_name]] <- run_deg_pair_bulk(
    dds = dds,
    group1 = comp$group1,
    group2 = comp$group2,
    gene_anno_gtf = gene_anno_gtf,
    outdir_all = deg_stim_unstim_all_dir,
    outdir_sig = deg_stim_unstim_sig_dir,
    out_name = comp$out_name,
    padj_cutoff = 0.05,
    lfc_cutoff = 1
  )
}

stim_unstim_summary_all_bulk <- do.call(
  rbind,
  lapply(stim_unstim_deg_results_bulk, function(x) x$summary)
)

write.csv(
  stim_unstim_summary_all_bulk,
  file = file.path(deg_stim_unstim_dir, "DEG_summary_all_stim_vs_unstim_bulk.csv"),
  row.names = FALSE,
  quote = FALSE
)

saveRDS(
  list(
    stim_vs_unstim = stim_unstim_deg_results_bulk,
    summary_all = stim_unstim_summary_all_bulk
  ),
  file = file.path(deg_stim_unstim_dir, "DEG_results_all_stim_vs_unstim_bulk.rds")
)

############################################################
## 13.8 Save combined DEG summary
############################################################

deg_summary_all_bulk <- rbind(
  unstim_summary_all_bulk,
  stim_summary_all_bulk,
  stim_unstim_summary_all_bulk
)

write.csv(
  deg_summary_all_bulk,
  file = file.path(deg_dir, "DEG_summary_all_bulk.csv"),
  row.names = FALSE,
  quote = FALSE
)

cat("DESeq2 pairwise and three-group overall DEG analysis finished.\n")

############################################################
## 14. Pairwise volcano plots
############################################################

volcano_unstim_dir <- file.path(volcano_dir, "unstimulated_three_groups")
volcano_stim_dir <- file.path(volcano_dir, "CD19_stimulated_three_groups")
volcano_stim_unstim_dir <- file.path(volcano_dir, "stim_vs_unstim")

dir.create(volcano_unstim_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(volcano_stim_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(volcano_stim_unstim_dir, recursive = TRUE, showWarnings = FALSE)

plot_volcano_bulk <- function(res_anno,
                              out_name,
                              outdir,
                              padj_cutoff = 0.05,
                              lfc_cutoff = 1) {
  
  df <- as.data.frame(res_anno)
  
  df <- df[!is.na(df$log2FoldChange) & !is.na(df$pvalue), , drop = FALSE]
  
  df$change <- "NS"
  df$change[
    !is.na(df$padj) &
      df$padj < padj_cutoff &
      df$log2FoldChange >= lfc_cutoff
  ] <- "Up"
  
  df$change[
    !is.na(df$padj) &
      df$padj < padj_cutoff &
      df$log2FoldChange <= -lfc_cutoff
  ] <- "Down"
  
  df$change <- factor(df$change, levels = c("Up", "Down", "NS"))
  
  df$neglog10padj <- -log10(df$padj)
  df$neglog10padj[is.infinite(df$neglog10padj)] <- NA
  
  n_up <- length(unique(df$Geneid[df$change == "Up" & !is.na(df$Geneid)]))
  n_down <- length(unique(df$Geneid[df$change == "Down" & !is.na(df$Geneid)]))
  
  volcano_colors <- c(
    "Up" = "#E64B35",
    "Down" = "#4DBBD5",
    "NS" = "grey75"
  )
  
  p <- ggplot(df, aes(x = log2FoldChange, y = neglog10padj)) +
    geom_point(aes(color = change), size = 1.8, alpha = 0.8) +
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
    theme_bw(base_size = 14) +
    labs(
      title = paste0(out_name, " bulk volcano"),
      x = "log2FoldChange",
      y = "-log10(adjusted p-value)",
      color = NULL
    ) +
    annotate(
      "text",
      x = min(df$log2FoldChange, na.rm = TRUE),
      y = max(df$neglog10padj, na.rm = TRUE),
      label = paste0("Down: ", n_down),
      hjust = 0,
      vjust = 1,
      size = 4.5,
      color = "#4DBBD5"
    ) +
    annotate(
      "text",
      x = max(df$log2FoldChange, na.rm = TRUE),
      y = max(df$neglog10padj, na.rm = TRUE),
      label = paste0("Up: ", n_up),
      hjust = 1,
      vjust = 1,
      size = 4.5,
      color = "#E64B35"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right"
    )
  
  ggsave(
    filename = file.path(outdir, paste0(out_name, "_volcano_bulk.pdf")),
    plot = p,
    width = 7,
    height = 6
  )
  
  ggsave(
    filename = file.path(outdir, paste0(out_name, "_volcano_bulk.png")),
    plot = p,
    width = 7,
    height = 6,
    dpi = 300
  )
  
  return(p)
}

for (nm in names(unstim_deg_results_bulk)) {
  cat("Plotting unstim volcano:", nm, "\n")
  plot_volcano_bulk(
    res_anno = unstim_deg_results_bulk[[nm]]$res_anno,
    out_name = nm,
    outdir = volcano_unstim_dir,
    padj_cutoff = 0.05,
    lfc_cutoff = 1
  )
}

for (nm in names(stim_deg_results_bulk)) {
  cat("Plotting stim volcano:", nm, "\n")
  plot_volcano_bulk(
    res_anno = stim_deg_results_bulk[[nm]]$res_anno,
    out_name = nm,
    outdir = volcano_stim_dir,
    padj_cutoff = 0.05,
    lfc_cutoff = 1
  )
}

for (nm in names(stim_unstim_deg_results_bulk)) {
  cat("Plotting stim vs unstim volcano:", nm, "\n")
  plot_volcano_bulk(
    res_anno = stim_unstim_deg_results_bulk[[nm]]$res_anno,
    out_name = nm,
    outdir = volcano_stim_unstim_dir,
    padj_cutoff = 0.05,
    lfc_cutoff = 1
  )
}

cat("Pairwise volcano plots finished. Three-group overall volcano plots were not generated.\n")

############################################################
## 15. Three-group overall top1000 DEG heatmaps
############################################################

heatmap_unstim_dir <- file.path(
  three_group_heatmap_dir,
  "unstimulated_three_group_overall"
)

heatmap_stim_dir <- file.path(
  three_group_heatmap_dir,
  "CD19_stimulated_three_group_overall"
)

dir.create(heatmap_unstim_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(heatmap_stim_dir, recursive = TRUE, showWarnings = FALSE)

############################################################
## 15.1 General function for three-group top1000 heatmap
############################################################

plot_three_group_top1000_heatmap_mark_custom_updown <- function(
    comp_res,
    meta,
    vsd_mat,
    groups_keep,
    group_label_map,
    comp_name,
    title_text,
    outdir,
    up_genes = character(0),
    down_genes = character(0),
    top_n = 1000,
    gene2cat,
    category_symbols,
    category_text_colors,
    category_levels,
    symbol_labels,
    symbol_pch,
    file_suffix,
    pdf_family = "Arial"
) {
  
  if (is.null(comp_res$res_sig) || nrow(comp_res$res_sig) == 0) {
    stop(paste0(comp_name, " res_sig is empty."))
  }
  
  res_df <- comp_res$res_sig
  res_df <- res_df[!is.na(res_df$padj), , drop = FALSE]
  
  if (nrow(res_df) == 0) {
    stop(paste0(comp_name, " has no genes after padj filtering."))
  }
  
  need_cols <- c("Geneid_raw", "Geneid", "gene_name", "log2FoldChange", "padj")
  miss_cols <- setdiff(need_cols, colnames(res_df))
  
  if (length(miss_cols) > 0) {
    stop(paste0("res_sig is missing columns: ", paste(miss_cols, collapse = ", ")))
  }
  
  res_df <- res_df[order(res_df$padj, -abs(res_df$log2FoldChange)), , drop = FALSE]
  res_df <- res_df[!duplicated(res_df$Geneid), , drop = FALSE]
  
  top_n <- min(top_n, nrow(res_df))
  top_df <- res_df[seq_len(top_n), , drop = FALSE]
  
  anno_map <- unique(top_df[, c("Geneid_raw", "Geneid", "gene_name", "log2FoldChange")])
  anno_map <- anno_map[!duplicated(anno_map$Geneid_raw), , drop = FALSE]
  
  meta_sub <- meta[meta$group %in% groups_keep, , drop = FALSE]
  meta_sub$group <- factor(meta_sub$group, levels = groups_keep)
  meta_sub <- meta_sub[order(meta_sub$group), , drop = FALSE]
  
  sample_use <- rownames(meta_sub)
  sample_use <- sample_use[sample_use %in% colnames(vsd_mat)]
  
  if (length(sample_use) == 0) {
    stop(paste0(comp_name, " no matched samples in vst matrix."))
  }
  
  meta_sub <- meta_sub[sample_use, , drop = FALSE]
  
  genes_raw_in_mat <- anno_map$Geneid_raw[anno_map$Geneid_raw %in% rownames(vsd_mat)]
  
  use_geneid_mode <- FALSE
  
  if (length(genes_raw_in_mat) < 10) {
    genes_id_in_mat <- anno_map$Geneid[anno_map$Geneid %in% rownames(vsd_mat)]
    if (length(genes_id_in_mat) > length(genes_raw_in_mat)) {
      use_geneid_mode <- TRUE
    }
  }
  
  if (!use_geneid_mode) {
    plot_mat <- vsd_mat[genes_raw_in_mat, sample_use, drop = FALSE]
    row_info <- anno_map[match(rownames(plot_mat), anno_map$Geneid_raw), , drop = FALSE]
    ord <- match(row_info$Geneid_raw, top_df$Geneid_raw)
  } else {
    genes_id_in_mat <- anno_map$Geneid[anno_map$Geneid %in% rownames(vsd_mat)]
    plot_mat <- vsd_mat[genes_id_in_mat, sample_use, drop = FALSE]
    row_info <- anno_map[match(rownames(plot_mat), anno_map$Geneid), , drop = FALSE]
    ord <- match(row_info$Geneid, top_df$Geneid)
  }
  
  if (nrow(plot_mat) == 0) {
    stop(paste0(comp_name, " no matched genes for heatmap."))
  }
  
  keep_idx <- !duplicated(rownames(plot_mat))
  plot_mat <- plot_mat[keep_idx, , drop = FALSE]
  row_info <- row_info[keep_idx, , drop = FALSE]
  ord <- ord[keep_idx]
  
  plot_mat <- plot_mat[order(ord), , drop = FALSE]
  row_info <- row_info[order(ord), , drop = FALSE]
  
  plot_mat_scaled <- t(scale(t(plot_mat)))
  plot_mat_scaled[is.na(plot_mat_scaled)] <- 0
  
  gene_name_use <- row_info$gene_name
  gene_name_use[is.na(gene_name_use)] <- ""
  
  valid_label_idx <- which(
    gene_name_use != "" &
      !grepl("^ENSG", gene_name_use, ignore.case = TRUE)
  )
  
  gene_name_valid <- gene_name_use[valid_label_idx]
  
  up_idx <- valid_label_idx[gene_name_valid %in% up_genes]
  down_idx <- valid_label_idx[gene_name_valid %in% down_genes]
  
  up_idx <- sort(unique(up_idx))
  down_idx <- sort(unique(down_idx))
  
  at_idx <- c(up_idx, down_idx)
  base_label_vec <- c(gene_name_use[up_idx], gene_name_use[down_idx])
  
  label_cat <- unname(gene2cat[base_label_vec])
  label_symbol <- unname(category_symbols[label_cat])
  label_col <- unname(category_text_colors[label_cat])
  
  label_symbol[is.na(label_symbol)] <- ""
  label_col[is.na(label_col)] <- "black"
  
  label_vec <- paste0(base_label_vec, " ", label_symbol)
  
  label_group <- rep("No", nrow(row_info))
  label_group[up_idx] <- "Custom_Up"
  label_group[down_idx] <- "Custom_Down"
  
  label_category_full <- rep(NA_character_, nrow(row_info))
  label_symbol_full <- rep(NA_character_, nrow(row_info))
  
  all_labeled_idx <- c(up_idx, down_idx)
  label_category_full[all_labeled_idx] <- label_cat
  label_symbol_full[all_labeled_idx] <- label_symbol
  
  out_table <- data.frame(
    Geneid_raw = row_info$Geneid_raw,
    Geneid = row_info$Geneid,
    gene_name = row_info$gene_name,
    log2FoldChange = row_info$log2FoldChange,
    label_group = label_group,
    label_category = label_category_full,
    label_symbol = label_symbol_full,
    plot_mat,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  
  write.csv(
    out_table,
    file = file.path(outdir, paste0(comp_name, "_top1000_expression_matrix.csv")),
    row.names = FALSE,
    quote = FALSE
  )
  
  label_table <- data.frame(
    gene_name = base_label_vec,
    label_show = label_vec,
    category = label_cat,
    symbol = label_symbol,
    text_color = label_col,
    log2FoldChange = c(row_info$log2FoldChange[up_idx], row_info$log2FoldChange[down_idx]),
    label_group = c(rep("Custom_Up", length(up_idx)), rep("Custom_Down", length(down_idx))),
    stringsAsFactors = FALSE
  )
  
  write.csv(
    label_table,
    file = file.path(outdir, paste0(comp_name, "_custom_label_genes_in_top1000.csv")),
    row.names = FALSE,
    quote = FALSE
  )
  
  genes_in_top1000 <- unique(gene_name_valid)
  
  up_not_found <- setdiff(up_genes, genes_in_top1000)
  down_not_found <- setdiff(down_genes, genes_in_top1000)
  
  not_found_table <- data.frame(
    gene_name = c(up_not_found, down_not_found),
    label_group = c(
      rep("Custom_Up_not_in_top1000", length(up_not_found)),
      rep("Custom_Down_not_in_top1000", length(down_not_found))
    ),
    stringsAsFactors = FALSE
  )
  
  write.csv(
    not_found_table,
    file = file.path(outdir, paste0(comp_name, "_custom_label_genes_not_in_top1000.csv")),
    row.names = FALSE,
    quote = FALSE
  )
  
  col_fun <- circlize::colorRamp2(
    c(-2, -1, 0, 1, 2),
    c("#3B4CC0", "#88AEEF", "white", "#F4A582", "#B40426")
  )
  
  right_anno <- NULL
  
  if (length(at_idx) > 0) {
    
    gene_mark_width <- max_text_width(
      label_vec,
      gp = gpar(fontsize = 7, fontface = "bold", fontfamily = pdf_family)
    ) + unit(4, "mm")
    
    right_anno <- rowAnnotation(
      gene_mark = anno_mark(
        at = at_idx,
        labels = label_vec,
        which = "row",
        side = "right",
        labels_gp = gpar(
          fontsize = 7,
          col = "black",
          fontface = "bold",
          fontfamily = pdf_family
        ),
        link_width = unit(4.5, "mm"),
        link_gp = gpar(lwd = 0.8)
      ),
      width = gene_mark_width
    )
  }
  
  group_for_split <- factor(meta_sub$group, levels = groups_keep)
  
  bottom_ha <- HeatmapAnnotation(
    Group = anno_block(
      gp = gpar(fill = "white", col = NA),
      labels = group_label_map[groups_keep],
      labels_gp = gpar(
        fontsize = 9,
        fontface = "bold",
        col = "black",
        fontfamily = pdf_family
      ),
      labels_rot = 0
    ),
    which = "column",
    show_annotation_name = FALSE,
    height = unit(5, "mm")
  )
  
  ht <- Heatmap(
    plot_mat_scaled,
    name = "Z-score",
    col = col_fun,
    cluster_rows = TRUE,
    cluster_columns = FALSE,
    cluster_column_slices = FALSE,
    column_split = group_for_split,
    column_gap = unit(0, "mm"),
    bottom_annotation = bottom_ha,
    show_row_names = FALSE,
    show_column_names = FALSE,
    show_column_dend = FALSE,
    column_title = title_text,
    column_title_gp = gpar(fontsize = 9, fontface = "bold", fontfamily = pdf_family),
    rect_gp = gpar(col = NA),
    right_annotation = right_anno,
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 8, fontface = "plain", fontfamily = pdf_family),
      labels_gp = gpar(fontsize = 7, fontfamily = pdf_family)
    )
  )
  
  lgd_cat <- Legend(
    title = NULL,
    labels = symbol_labels,
    type = "points",
    pch = symbol_pch,
    legend_gp = gpar(col = "black", fill = "black"),
    ncol = 1,
    by_row = TRUE,
    gap = unit(1.5, "mm"),
    labels_gp = gpar(fontsize = 7, col = "black", fontfamily = pdf_family),
    grid_width = unit(4, "mm"),
    grid_height = unit(4, "mm")
  )
  
  cairo_pdf(
    filename = file.path(outdir, paste0(comp_name, "_top1000_bulk_", file_suffix, ".pdf")),
    width = 6,
    height = 4.7,
    family = pdf_family
  )
  
  draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    annotation_legend_list = list(lgd_cat),
    merge_legends = FALSE,
    align_annotation_legend = "heatmap_top",
    padding = unit(c(8, 12, 2, 2), "mm")
  )
  
  dev.off()
  
  png(
    filename = file.path(outdir, paste0(comp_name, "_top1000_bulk_", file_suffix, ".png")),
    width = 1750,
    height = 1300,
    res = 300
  )
  
  draw(
    ht,
    heatmap_legend_side = "right",
    annotation_legend_side = "right",
    annotation_legend_list = list(lgd_cat),
    merge_legends = FALSE,
    align_annotation_legend = "heatmap_top",
    padding = unit(c(8, 12, 2, 2), "mm")
  )
  
  dev.off()
  
  cat(comp_name, " three-group overall heatmap finished.\n")
  cat("Labeled up genes:", length(up_idx), "\n")
  cat("Labeled down genes:", length(down_idx), "\n")
  cat("Up genes not in top1000:", length(up_not_found), "\n")
  cat("Down genes not in top1000:", length(down_not_found), "\n")
  
  return(list(
    raw_mat = plot_mat,
    scaled_mat = plot_mat_scaled,
    row_info = row_info,
    up_label_genes = gene_name_use[up_idx],
    down_label_genes = gene_name_use[down_idx],
    up_not_found = up_not_found,
    down_not_found = down_not_found,
    label_categories = label_cat,
    label_symbols = label_symbol
  ))
}

############################################################
## 15.2 Unstimulated three-group heatmap gene categories
############################################################

up_cell_cycle <- c(
  "CDC25B", "E2F2", "CDK1", "CCNA2", "AURKB"
)

up_dna_replication <- c(
  "MCM2", "CHEK1", "FEN1", "POLD1", "RFC4"
)

up_dna_repair <- c(
  "BRCA1", "RAD51", "FANCD2", "XRCC2", "BRIP1"
)

down_Cytokine_cytokine_receptor_interaction <- c(
  "IL21R", "TNFSF4", "IL12RB2", "EBI3"
)

down_p53_pathway <- c(
  "SESN1", "TP53I3", "GADD45A", "ZMAT3", "CCNG2"
)

down_chemokine_signaling_pathway <- c(
  "CCR4", "CCR7", "SRC", "PLCG2"
)

up_genes_unstim <- unique(c(
  up_cell_cycle,
  up_dna_replication,
  up_dna_repair
))

down_genes_unstim <- unique(c(
  down_Cytokine_cytokine_receptor_interaction,
  down_p53_pathway,
  down_chemokine_signaling_pathway
))

category_levels_unstim <- c(
  "Cell cycle",
  "DNA replication",
  "DNA repair",
  "Cytokine-cytokine receptor interaction",
  "p53 pathway",
  "Chemokine signaling pathway"
)

category_symbols_unstim <- c(
  "Cell cycle" = "●",
  "DNA replication" = "■",
  "DNA repair" = "▲",
  "Cytokine-cytokine receptor interaction" = "○",
  "p53 pathway" = "□",
  "Chemokine signaling pathway" = "△"
)

gene2cat_unstim <- c(
  setNames(rep("Cell cycle", length(up_cell_cycle)), up_cell_cycle),
  setNames(rep("DNA replication", length(up_dna_replication)), up_dna_replication),
  setNames(rep("DNA repair", length(up_dna_repair)), up_dna_repair),
  setNames(
    rep(
      "Cytokine-cytokine receptor interaction",
      length(down_Cytokine_cytokine_receptor_interaction)
    ),
    down_Cytokine_cytokine_receptor_interaction
  ),
  setNames(rep("p53 pathway", length(down_p53_pathway)), down_p53_pathway),
  setNames(
    rep("Chemokine signaling pathway", length(down_chemokine_signaling_pathway)),
    down_chemokine_signaling_pathway
  )
)

category_text_colors_unstim <- setNames(
  rep("black", length(category_levels_unstim)),
  category_levels_unstim
)

symbol_pch_unstim <- c(16, 15, 17, 1, 0, 2)

group_label_map_unstim <- c(
  "DMSO_un" = "DMSO",
  "GSK_un" = "GSK761",
  "PBMC_un" = "PBMC"
)

heatmap_unstim_3group_overall <- plot_three_group_top1000_heatmap_mark_custom_updown(
  comp_res = res_unstim_3group_bulk,
  meta = meta,
  vsd_mat = vsd_mat,
  groups_keep = c("DMSO_un", "GSK_un", "PBMC_un"),
  group_label_map = group_label_map_unstim,
  comp_name = "DEG_unstim_3group_overall_bulk",
  title_text = "Unstim bulk: 3-group overall",
  outdir = heatmap_unstim_dir,
  up_genes = up_genes_unstim,
  down_genes = down_genes_unstim,
  top_n = 1000,
  gene2cat = gene2cat_unstim,
  category_symbols = category_symbols_unstim,
  category_text_colors = category_text_colors_unstim,
  category_levels = category_levels_unstim,
  symbol_labels = category_levels_unstim,
  symbol_pch = symbol_pch_unstim,
  file_suffix = "unstim_updown",
  pdf_family = "Arial"
)

############################################################
## 15.3 CD19-stimulated three-group heatmap gene categories
############################################################

up_nfkb <- c(
  "CD40", "BTK", "SYK", "BCL2A1", "CXCL8"
)

up_protein_processing_er <- c(
  "DNAJB11", "HSP90AB1", "SEC62", "LMAN1", "DDIT3"
)

up_nucleocytoplasmic_transport <- c(
  "IPO5", "KPNA2", "TNPO1", "RAN", "NUP205"
)

down_chemokine <- c(
  "CCR1", "CCR5", "CXCR6", "JAK3", "GNAI2"
)

down_nk_cytotoxicity <- c(
  "NCR2", "NCR3", "CD244", "LCK", "FCER1G"
)

down_integrin_signaling <- c(
  "ITGAL", "ITGB2", "ITGAM", "PIK3CD", "RAC2"
)

up_genes_stim <- unique(c(
  up_nfkb,
  up_protein_processing_er,
  up_nucleocytoplasmic_transport
))

down_genes_stim <- unique(c(
  down_chemokine,
  down_nk_cytotoxicity,
  down_integrin_signaling
))

category_levels_stim <- c(
  "NF-kappa B",
  "Protein processing in endoplasmic reticulum",
  "Nucleocytoplasmic transport",
  "Chemokine signaling",
  "NK mediated cytotoxicity",
  "Integrin signaling"
)

category_symbols_stim <- c(
  "NF-kappa B" = "●",
  "Protein processing in endoplasmic reticulum" = "■",
  "Nucleocytoplasmic transport" = "▲",
  "Chemokine signaling" = "○",
  "NK mediated cytotoxicity" = "□",
  "Integrin signaling" = "△"
)

gene2cat_stim <- c(
  setNames(rep("NF-kappa B", length(up_nfkb)), up_nfkb),
  setNames(
    rep("Protein processing in endoplasmic reticulum", length(up_protein_processing_er)),
    up_protein_processing_er
  ),
  setNames(
    rep("Nucleocytoplasmic transport", length(up_nucleocytoplasmic_transport)),
    up_nucleocytoplasmic_transport
  ),
  setNames(rep("Chemokine signaling", length(down_chemokine)), down_chemokine),
  setNames(rep("NK mediated cytotoxicity", length(down_nk_cytotoxicity)), down_nk_cytotoxicity),
  setNames(rep("Integrin signaling", length(down_integrin_signaling)), down_integrin_signaling)
)

category_text_colors_stim <- setNames(
  rep("black", length(category_levels_stim)),
  category_levels_stim
)

symbol_pch_stim <- c(16, 15, 17, 1, 0, 2)

group_label_map_stim <- c(
  "DMSO_CD19" = "DMSO",
  "GSK_CD19" = "GSK761",
  "PBMC_CD19" = "PBMC"
)

heatmap_stim_3group_overall <- plot_three_group_top1000_heatmap_mark_custom_updown(
  comp_res = res_stim_3group_bulk,
  meta = meta,
  vsd_mat = vsd_mat,
  groups_keep = c("DMSO_CD19", "GSK_CD19", "PBMC_CD19"),
  group_label_map = group_label_map_stim,
  comp_name = "DEG_stim_3group_overall_bulk",
  title_text = "Stim bulk: 3-group overall",
  outdir = heatmap_stim_dir,
  up_genes = up_genes_stim,
  down_genes = down_genes_stim,
  top_n = 1000,
  gene2cat = gene2cat_stim,
  category_symbols = category_symbols_stim,
  category_text_colors = category_text_colors_stim,
  category_levels = category_levels_stim,
  symbol_labels = category_levels_stim,
  symbol_pch = symbol_pch_stim,
  file_suffix = "stim_updown",
  pdf_family = "Arial"
)

cat("Three-group overall top1000 DEG heatmaps finished.\n")

############################################################
## 16. KEGG and GO enrichment analysis
############################################################

kegg_go_dir <- file.path(result_dir, "06_KEGG_GO")
selected_pathway_dir <- file.path(result_dir, "07_Selected_pathway")

dir.create(kegg_go_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(selected_pathway_dir, recursive = TRUE, showWarnings = FALSE)

enrich_unstim_dir <- file.path(kegg_go_dir, "unstimulated_three_groups")
enrich_stim_dir <- file.path(kegg_go_dir, "CD19_stimulated_three_groups")
enrich_stim_unstim_dir <- file.path(kegg_go_dir, "stim_vs_unstim")

dir.create(enrich_unstim_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(enrich_stim_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(enrich_stim_unstim_dir, recursive = TRUE, showWarnings = FALSE)

options(timeout = 1200)
options(download.file.method = "libcurl")

packages_bioc_enrich <- c(
  "clusterProfiler",
  "org.Hs.eg.db"
)

packages_cran_enrich <- c(
  "data.table",
  "stringr"
)

for (pkg in packages_cran_enrich) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

for (pkg in packages_bioc_enrich) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg)
  }
}

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(data.table)
  library(stringr)
})

############################################################
## 16.1 Helper functions for KEGG/GO enrichment
############################################################

convert_gene_to_entrez <- function(gene_vec) {
  
  gene_vec <- unique(gene_vec)
  gene_vec <- gene_vec[!is.na(gene_vec) & gene_vec != ""]
  
  if (length(gene_vec) == 0) {
    return(data.frame())
  }
  
  gene_map <- tryCatch({
    bitr(
      gene_vec,
      fromType = "ENSEMBL",
      toType = c("ENTREZID", "SYMBOL"),
      OrgDb = org.Hs.eg.db
    )
  }, error = function(e) {
    message("ID conversion failed: ", e$message)
    return(data.frame())
  })
  
  if (is.null(gene_map) || nrow(gene_map) == 0) {
    return(data.frame())
  }
  
  gene_map <- gene_map[!duplicated(gene_map$ENSEMBL), , drop = FALSE]
  rownames(gene_map) <- NULL
  
  return(gene_map)
}

get_direction_gene_map <- function(sig_df, direction_keep) {
  
  df_sub <- sig_df %>%
    dplyr::filter(direction == direction_keep)
  
  gene_vec <- unique(df_sub$Geneid)
  gene_map <- convert_gene_to_entrez(gene_vec)
  
  return(list(
    sig_df = df_sub,
    gene_vec = gene_vec,
    gene_map = gene_map
  ))
}

tidy_enrich_result <- function(ego, comparison, direction, db_type, ontology = NA) {
  
  if (is.null(ego)) {
    return(data.frame())
  }
  
  ego_df <- as.data.frame(ego)
  
  if (is.null(ego_df) || nrow(ego_df) == 0) {
    return(data.frame())
  }
  
  ego_df$comparison <- comparison
  ego_df$direction <- direction
  ego_df$db_type <- db_type
  ego_df$ontology <- ontology
  
  keep_cols <- c(
    "comparison", "direction", "db_type", "ontology",
    "ID", "Description", "GeneRatio", "BgRatio",
    "RichFactor", "FoldEnrichment", "zScore",
    "pvalue", "p.adjust", "qvalue", "geneID", "Count"
  )
  
  keep_cols <- keep_cols[keep_cols %in% colnames(ego_df)]
  ego_df <- ego_df[, keep_cols, drop = FALSE]
  
  ego_df <- ego_df %>%
    dplyr::arrange(p.adjust, pvalue)
  
  rownames(ego_df) <- NULL
  return(ego_df)
}

run_kegg_one_direction <- function(gene_map, comparison, direction) {
  
  if (is.null(gene_map) || nrow(gene_map) == 0) {
    return(list(result_df = data.frame(), n_input = 0, n_mapped = 0))
  }
  
  gene_use <- unique(gene_map$ENTREZID)
  gene_use <- gene_use[!is.na(gene_use) & gene_use != ""]
  
  if (length(gene_use) == 0) {
    return(list(result_df = data.frame(), n_input = nrow(gene_map), n_mapped = 0))
  }
  
  ekegg <- tryCatch({
    enrichKEGG(
      gene = gene_use,
      organism = "hsa",
      keyType = "kegg",
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff = 0.2
    )
  }, error = function(e) {
    message("KEGG enrichment failed - ", comparison, " - ", direction, ": ", e$message)
    return(NULL)
  })
  
  if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
    ekegg <- tryCatch({
      setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
    }, error = function(e) {
      ekegg
    })
  }
  
  ekegg_df <- tidy_enrich_result(
    ego = ekegg,
    comparison = comparison,
    direction = direction,
    db_type = "KEGG",
    ontology = "KEGG"
  )
  
  return(list(
    result_df = ekegg_df,
    n_input = nrow(gene_map),
    n_mapped = length(gene_use)
  ))
}

run_go_one_direction <- function(gene_map, comparison, direction, ont_use) {
  
  if (is.null(gene_map) || nrow(gene_map) == 0) {
    return(list(result_df = data.frame(), n_input = 0, n_mapped = 0))
  }
  
  gene_use <- unique(gene_map$ENTREZID)
  gene_use <- gene_use[!is.na(gene_use) & gene_use != ""]
  
  if (length(gene_use) == 0) {
    return(list(result_df = data.frame(), n_input = nrow(gene_map), n_mapped = 0))
  }
  
  ego <- tryCatch({
    enrichGO(
      gene = gene_use,
      OrgDb = org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = ont_use,
      pvalueCutoff = 0.05,
      pAdjustMethod = "BH",
      qvalueCutoff = 0.2,
      readable = TRUE
    )
  }, error = function(e) {
    message("GO enrichment failed - ", comparison, " - ", direction, " - ", ont_use, ": ", e$message)
    return(NULL)
  })
  
  ego_df <- tidy_enrich_result(
    ego = ego,
    comparison = comparison,
    direction = direction,
    db_type = "GO",
    ontology = ont_use
  )
  
  return(list(
    result_df = ego_df,
    n_input = nrow(gene_map),
    n_mapped = length(gene_use)
  ))
}

make_enrich_summary_row <- function(comparison, direction, db_type, ontology,
                                    n_input_gene, n_mapped_gene, result_df) {
  
  if (is.null(result_df) || nrow(result_df) == 0) {
    return(data.frame(
      comparison = comparison,
      direction = direction,
      db_type = db_type,
      ontology = ontology,
      n_input_gene = n_input_gene,
      n_mapped_gene = n_mapped_gene,
      n_enriched_terms = 0,
      top_term = NA,
      top_p_adjust = NA,
      stringsAsFactors = FALSE
    ))
  }
  
  result_df <- result_df %>%
    dplyr::arrange(p.adjust, pvalue)
  
  data.frame(
    comparison = comparison,
    direction = direction,
    db_type = db_type,
    ontology = ontology,
    n_input_gene = n_input_gene,
    n_mapped_gene = n_mapped_gene,
    n_enriched_terms = nrow(result_df),
    top_term = result_df$Description[1],
    top_p_adjust = result_df$p.adjust[1],
    stringsAsFactors = FALSE
  )
}

############################################################
## 16.2 General enrichment runner
############################################################

run_kegg_go_for_deg_results <- function(deg_results_list, outdir_base) {
  
  outdir_kegg <- file.path(outdir_base, "KEGG")
  outdir_go <- file.path(outdir_base, "GO")
  outdir_sum <- file.path(outdir_base, "summary")
  
  dir.create(outdir_kegg, recursive = TRUE, showWarnings = FALSE)
  dir.create(outdir_go, recursive = TRUE, showWarnings = FALSE)
  dir.create(outdir_sum, recursive = TRUE, showWarnings = FALSE)
  
  all_kegg_combined_list <- list()
  all_go_combined_list <- list()
  all_summary_list <- list()
  
  for (comp_name in names(deg_results_list)) {
    
    cat("Running KEGG/GO enrichment:", comp_name, "\n")
    
    sig_df <- deg_results_list[[comp_name]]$res_sig
    
    up_info <- get_direction_gene_map(sig_df, "Up")
    down_info <- get_direction_gene_map(sig_df, "Down")
    
    ## KEGG
    kegg_up <- run_kegg_one_direction(
      gene_map = up_info$gene_map,
      comparison = comp_name,
      direction = "Up"
    )
    
    kegg_down <- run_kegg_one_direction(
      gene_map = down_info$gene_map,
      comparison = comp_name,
      direction = "Down"
    )
    
    kegg_combined_df <- dplyr::bind_rows(
      kegg_up$result_df,
      kegg_down$result_df
    )
    
    if (nrow(kegg_combined_df) > 0) {
      kegg_combined_df <- kegg_combined_df %>%
        dplyr::arrange(direction, p.adjust, dplyr::desc(FoldEnrichment))
    }
    
    write.csv(
      kegg_combined_df,
      file = file.path(outdir_kegg, paste0(comp_name, "_KEGG_up_down_combined.csv")),
      row.names = FALSE,
      quote = TRUE,
      na = ""
    )
    
    all_kegg_combined_list[[comp_name]] <- kegg_combined_df
    
    sum_kegg_up <- make_enrich_summary_row(
      comparison = comp_name,
      direction = "Up",
      db_type = "KEGG",
      ontology = "KEGG",
      n_input_gene = length(unique(up_info$gene_vec)),
      n_mapped_gene = kegg_up$n_mapped,
      result_df = kegg_up$result_df
    )
    
    sum_kegg_down <- make_enrich_summary_row(
      comparison = comp_name,
      direction = "Down",
      db_type = "KEGG",
      ontology = "KEGG",
      n_input_gene = length(unique(down_info$gene_vec)),
      n_mapped_gene = kegg_down$n_mapped,
      result_df = kegg_down$result_df
    )
    
    ## GO BP / CC / MF
    go_up_bp <- run_go_one_direction(up_info$gene_map, comp_name, "Up", "BP")
    go_up_cc <- run_go_one_direction(up_info$gene_map, comp_name, "Up", "CC")
    go_up_mf <- run_go_one_direction(up_info$gene_map, comp_name, "Up", "MF")
    
    go_down_bp <- run_go_one_direction(down_info$gene_map, comp_name, "Down", "BP")
    go_down_cc <- run_go_one_direction(down_info$gene_map, comp_name, "Down", "CC")
    go_down_mf <- run_go_one_direction(down_info$gene_map, comp_name, "Down", "MF")
    
    go_combined_df <- dplyr::bind_rows(
      go_up_bp$result_df,
      go_up_cc$result_df,
      go_up_mf$result_df,
      go_down_bp$result_df,
      go_down_cc$result_df,
      go_down_mf$result_df
    )
    
    if (nrow(go_combined_df) > 0) {
      go_combined_df <- go_combined_df %>%
        dplyr::arrange(direction, ontology, p.adjust, dplyr::desc(FoldEnrichment))
    }
    
    write.csv(
      go_combined_df,
      file = file.path(outdir_go, paste0(comp_name, "_GO_BP_CC_MF_up_down_combined.csv")),
      row.names = FALSE,
      quote = TRUE,
      na = ""
    )
    
    all_go_combined_list[[comp_name]] <- go_combined_df
    
    summary_df <- dplyr::bind_rows(
      sum_kegg_up,
      sum_kegg_down,
      make_enrich_summary_row(comp_name, "Up", "GO", "BP", length(unique(up_info$gene_vec)), go_up_bp$n_mapped, go_up_bp$result_df),
      make_enrich_summary_row(comp_name, "Up", "GO", "CC", length(unique(up_info$gene_vec)), go_up_cc$n_mapped, go_up_cc$result_df),
      make_enrich_summary_row(comp_name, "Up", "GO", "MF", length(unique(up_info$gene_vec)), go_up_mf$n_mapped, go_up_mf$result_df),
      make_enrich_summary_row(comp_name, "Down", "GO", "BP", length(unique(down_info$gene_vec)), go_down_bp$n_mapped, go_down_bp$result_df),
      make_enrich_summary_row(comp_name, "Down", "GO", "CC", length(unique(down_info$gene_vec)), go_down_cc$n_mapped, go_down_cc$result_df),
      make_enrich_summary_row(comp_name, "Down", "GO", "MF", length(unique(down_info$gene_vec)), go_down_mf$n_mapped, go_down_mf$result_df)
    )
    
    write.csv(
      summary_df,
      file = file.path(outdir_sum, paste0(comp_name, "_KEGG_GO_summary.csv")),
      row.names = FALSE,
      quote = TRUE,
      na = ""
    )
    
    all_summary_list[[comp_name]] <- summary_df
  }
  
  all_summary <- dplyr::bind_rows(all_summary_list)
  
  write.csv(
    all_summary,
    file = file.path(outdir_sum, "ALL_pairwise_KEGG_GO_summary_combined.csv"),
    row.names = FALSE,
    quote = TRUE,
    na = ""
  )
  
  saveRDS(
    list(
      deg_results = deg_results_list,
      kegg_all = all_kegg_combined_list,
      go_all = all_go_combined_list,
      summary_all = all_summary
    ),
    file = file.path(outdir_base, "pairwise_KEGG_GO_results.rds")
  )
  
  return(list(
    kegg_all = all_kegg_combined_list,
    go_all = all_go_combined_list,
    summary_all = all_summary
  ))
}

############################################################
## 16.3 Run KEGG/GO for three comparison classes
############################################################

enrich_unstim_results <- run_kegg_go_for_deg_results(
  deg_results_list = unstim_deg_results_bulk,
  outdir_base = enrich_unstim_dir
)

enrich_stim_results <- run_kegg_go_for_deg_results(
  deg_results_list = stim_deg_results_bulk,
  outdir_base = enrich_stim_dir
)

enrich_stim_unstim_results <- run_kegg_go_for_deg_results(
  deg_results_list = stim_unstim_deg_results_bulk,
  outdir_base = enrich_stim_unstim_dir
)

cat("KEGG and GO enrichment analysis finished.\n")

############################################################
## 17. Selected KEGG pathway visualization
############################################################

selected_kegg_base_dir <- file.path(selected_pathway_dir, "KEGG_selected_bar")
selected_kegg_unstim_dir <- file.path(selected_kegg_base_dir, "unstimulated_three_groups")
selected_kegg_stim_dir <- file.path(selected_kegg_base_dir, "CD19_stimulated_three_groups")
selected_kegg_stim_unstim_dir <- file.path(selected_kegg_base_dir, "stim_vs_unstim")

dir.create(selected_kegg_unstim_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(selected_kegg_stim_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(selected_kegg_stim_unstim_dir, recursive = TRUE, showWarnings = FALSE)

padj_cutoff <- 0.05
sig_line <- -log10(padj_cutoff)

col_up <- "#4EA660"
col_down <- "#5292F7"

############################################################
## 17.1 Selected pathway lists
############################################################

selected_kegg_unstim <- list(
  "GSK761_un_vs_DMSO_un" = list(
    Down = c(
      "Cytokine-cytokine receptor interaction",
      "p53 signaling pathway",
      "Chemokine signaling pathway"
    ),
    Up = c(
      "Cell cycle",
      "DNA replication",
      "Fanconi anemia pathway",
      "Homologous recombination",
      "Base excision repair"
    )
  ),
  "PBMC_un_vs_DMSO_un" = list(
    Down = c(
      "p53 signaling pathway",
      "Chemokine signaling pathway",
      "Cytokine-cytokine receptor interaction"
    ),
    Up = c(
      "Cell cycle",
      "DNA replication",
      "Fanconi anemia pathway",
      "Homologous recombination",
      "Base excision repair"
    )
  )
)

selected_kegg_stim <- list(
  "PBMC_CD19_vs_DMSO_CD19" = list(
    Down = c(
      "Chemokine signaling pathway",
      "Natural killer cell mediated cytotoxicity",
      "Integrin signaling"
    ),
    Up = c(
      "Nucleocytoplasmic transport",
      "Cell cycle",
      "Protein processing in endoplasmic reticulum",
      "DNA replication",
      "NF-kappa B signaling pathway",
      "Mismatch repair"
    )
  ),
  "GSK761_CD19_vs_DMSO_CD19" = list(
    Down = c(
      "Chemokine signaling pathway",
      "Natural killer cell mediated cytotoxicity",
      "Integrin signaling"
    ),
    Up = c(
      "Cell cycle",
      "DNA replication",
      "Nucleocytoplasmic transport",
      "Protein processing in endoplasmic reticulum",
      "Mismatch repair",
      "NF-kappa B signaling pathway"
    )
  )
)

selected_kegg_stim_unstim <- list(
  "DMSO_CD19_vs_DMSO_un" = list(
    Up = c(
      "Natural killer cell mediated cytotoxicity",
      "Integrin signaling"
    ),
    Down = c(
      "Cell cycle",
      "Cytokine-cytokine receptor interaction",
      "DNA replication"
    )
  ),
  "GSK761_CD19_vs_GSK761_un" = list(
    Up = c(
      "Biosynthesis of amino acids",
      "Cytokine-cytokine receptor interaction",
      "Efferocytosis",
      "Notch signaling pathway",
      "Protein processing in endoplasmic reticulum",
      "Wnt signaling pathway"
    ),
    Down = c(
      "Apoptosis",
      "Chemokine signaling pathway",
      "Fc gamma R-mediated phagocytosis",
      "Inflammatory mediator regulation of TRP channels",
      "Integrin signaling",
      "Natural killer cell mediated cytotoxicity"
    )
  ),
  "PBMC_CD19_vs_PBMC_un" = list(
    Up = c(
      "Biosynthesis of amino acids",
      "Cytokine-cytokine receptor interaction",
      "Efferocytosis",
      "Notch signaling pathway",
      "Protein processing in endoplasmic reticulum",
      "Wnt signaling pathway"
    ),
    Down = c(
      "Apoptosis",
      "Chemokine signaling pathway",
      "Fc gamma R-mediated phagocytosis",
      "Inflammatory mediator regulation of TRP channels",
      "Integrin signaling",
      "Natural killer cell mediated cytotoxicity"
    )
  )
)

############################################################
## 17.2 Helper functions for selected KEGG plotting
############################################################

read_kegg_csv <- function(file) {
  
  df <- fread(file, data.table = FALSE, check.names = FALSE)
  colnames(df) <- trimws(colnames(df))
  
  need_cols <- c("comparison", "direction", "Description", "p.adjust", "Count")
  miss_cols <- setdiff(need_cols, colnames(df))
  
  if (length(miss_cols) > 0) {
    stop("File is missing required columns: ", basename(file), " -> ", paste(miss_cols, collapse = ", "))
  }
  
  df$comparison <- trimws(as.character(df$comparison))
  df$direction <- trimws(as.character(df$direction))
  df$Description <- trimws(as.character(df$Description))
  df$p.adjust <- suppressWarnings(as.numeric(df$p.adjust))
  df$Count <- suppressWarnings(as.numeric(df$Count))
  
  df <- df %>%
    dplyr::filter(
      !is.na(comparison),
      !is.na(direction),
      !is.na(Description),
      Description != "",
      !is.na(p.adjust),
      !is.na(Count)
    )
  
  df$direction <- dplyr::case_when(
    tolower(df$direction) == "up" ~ "Up",
    tolower(df$direction) == "down" ~ "Down",
    TRUE ~ df$direction
  )
  
  df <- df %>%
    dplyr::filter(direction %in% c("Up", "Down"))
  
  df$p.adjust[df$p.adjust <= 0] <- 1e-300
  df$neglog10_padj <- -log10(df$p.adjust)
  
  return(df)
}

pick_selected_kegg <- function(df, comp_name, selected_list_one) {
  
  down_terms <- selected_list_one$Down
  up_terms <- selected_list_one$Up
  
  df_sub <- df %>%
    dplyr::filter(comparison == comp_name) %>%
    dplyr::filter(p.adjust < padj_cutoff)
  
  df_down <- df_sub %>%
    dplyr::filter(direction == "Down", Description %in% down_terms) %>%
    dplyr::group_by(Description) %>%
    dplyr::arrange(p.adjust, dplyr::desc(Count), .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  
  if (nrow(df_down) > 0) {
    df_down$Description <- factor(df_down$Description, levels = down_terms)
    df_down <- df_down %>% dplyr::arrange(Description)
    df_down$Description <- as.character(df_down$Description)
  }
  
  df_up <- df_sub %>%
    dplyr::filter(direction == "Up", Description %in% up_terms) %>%
    dplyr::group_by(Description) %>%
    dplyr::arrange(p.adjust, dplyr::desc(Count), .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()
  
  if (nrow(df_up) > 0) {
    df_up$Description <- factor(df_up$Description, levels = up_terms)
    df_up <- df_up %>% dplyr::arrange(Description)
    df_up$Description <- as.character(df_up$Description)
  }
  
  list(
    down = df_down,
    up = df_up
  )
}

plot_selected_kegg <- function(df_down, df_up, comp_name, outdir) {
  
  if (nrow(df_down) == 0 && nrow(df_up) == 0) {
    message("No selected pathways to plot: ", comp_name)
    return(NULL)
  }
  
  if (nrow(df_down) > 0) {
    df_down <- df_down %>%
      dplyr::arrange(dplyr::desc(neglog10_padj)) %>%
      dplyr::mutate(
        y = seq(from = 1, length.out = n()),
        signed_x = -neglog10_padj,
        direction_show = "Down"
      )
  } else {
    df_down <- data.frame()
  }
  
  if (nrow(df_up) > 0) {
    start_y <- if (nrow(df_down) > 0) nrow(df_down) + 1 else 1
    
    df_up <- df_up %>%
      dplyr::arrange(dplyr::desc(neglog10_padj)) %>%
      dplyr::mutate(
        y = seq(from = start_y, length.out = n()),
        signed_x = neglog10_padj,
        direction_show = "Up"
      )
  } else {
    df_up <- data.frame()
  }
  
  plot_df <- dplyr::bind_rows(df_down, df_up)
  if (nrow(plot_df) == 0) return(NULL)
  
  plot_df$direction_show <- factor(plot_df$direction_show, levels = c("Down", "Up"))
  plot_df$count_label <- as.character(plot_df$Count)
  
  xmax <- max(plot_df$neglog10_padj, na.rm = TRUE)
  xmax <- max(xmax, sig_line)
  
  count_offset <- xmax * 0.04
  name_offset <- xmax * 0.03
  
  plot_df <- plot_df %>%
    dplyr::mutate(
      count_x = ifelse(
        direction_show == "Down",
        signed_x - count_offset,
        signed_x + count_offset
      ),
      path_x = ifelse(
        direction_show == "Down",
        name_offset,
        -name_offset
      ),
      path_hjust = ifelse(direction_show == "Down", 0, 1)
    )
  
  fig_h <- max(6, 0.48 * nrow(plot_df) + 2.5)
  
  p <- ggplot(plot_df) +
    geom_segment(
      aes(
        x = 0,
        xend = signed_x,
        y = y,
        yend = y,
        color = direction_show
      ),
      linewidth = 14,
      lineend = "butt"
    ) +
    geom_vline(
      xintercept = 0,
      color = "grey35",
      linewidth = 0.7
    ) +
    geom_vline(
      xintercept = c(-sig_line, sig_line),
      linetype = "dashed",
      color = "grey55",
      linewidth = 0.6
    ) +
    geom_text(
      aes(
        x = path_x,
        y = y,
        label = Description,
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
    scale_color_manual(
      values = c(
        "Down" = col_down,
        "Up" = col_up
      ),
      name = NULL,
      breaks = c("Down", "Up"),
      labels = c("Down", "Up")
    ) +
    scale_x_continuous(
      limits = c(-(xmax * 1.40), xmax * 1.40),
      breaks = pretty(c(-xmax, xmax), n = 6),
      labels = function(x) sprintf("%.1f", abs(x))
    ) +
    scale_y_continuous(
      breaks = plot_df$y,
      labels = rep("", nrow(plot_df)),
      expand = expansion(mult = c(0.03, 0.03))
    ) +
    labs(
      title = paste0(comp_name, " KEGG enrichment"),
      x = expression(-log[10](padj)),
      y = NULL
    ) +
    theme_classic(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      legend.position = "right",
      legend.text = element_text(size = 11),
      legend.title = element_blank()
    )
  
  out_prefix <- file.path(outdir, paste0(comp_name, "_KEGG_selected_bar"))
  
  ggsave(
    paste0(out_prefix, ".pdf"),
    plot = p,
    width = 14,
    height = fig_h,
    units = "in"
  )
  
  ggsave(
    paste0(out_prefix, ".png"),
    plot = p,
    width = 14,
    height = fig_h,
    units = "in",
    dpi = 300
  )
  
  write.csv(
    plot_df %>%
      dplyr::select(comparison, direction, Description, Count, p.adjust, neglog10_padj, y),
    file = paste0(out_prefix, "_plot_table.csv"),
    row.names = FALSE
  )
  
  return(p)
}

run_selected_kegg_plot <- function(kegg_dir, selected_kegg_list, outdir) {
  
  files_kegg <- list.files(kegg_dir, pattern = "\\.csv$", full.names = TRUE)
  
  if (length(files_kegg) == 0) {
    warning("No KEGG csv files found in: ", kegg_dir)
    return(NULL)
  }
  
  kegg_all <- dplyr::bind_rows(lapply(files_kegg, read_kegg_csv))
  
  for (comp_name in names(selected_kegg_list)) {
    
    message("Processing selected KEGG plot: ", comp_name)
    
    picked <- pick_selected_kegg(
      df = kegg_all,
      comp_name = comp_name,
      selected_list_one = selected_kegg_list[[comp_name]]
    )
    
    miss_down <- setdiff(selected_kegg_list[[comp_name]]$Down, picked$down$Description)
    miss_up <- setdiff(selected_kegg_list[[comp_name]]$Up, picked$up$Description)
    
    if (length(miss_down) > 0) {
      message("  Down not found: ", paste(miss_down, collapse = " ; "))
    }
    
    if (length(miss_up) > 0) {
      message("  Up not found: ", paste(miss_up, collapse = " ; "))
    }
    
    plot_selected_kegg(
      df_down = picked$down,
      df_up = picked$up,
      comp_name = comp_name,
      outdir = outdir
    )
  }
}

############################################################
## 17.3 Run selected KEGG plots
############################################################

run_selected_kegg_plot(
  kegg_dir = file.path(enrich_unstim_dir, "KEGG"),
  selected_kegg_list = selected_kegg_unstim,
  outdir = selected_kegg_unstim_dir
)

run_selected_kegg_plot(
  kegg_dir = file.path(enrich_stim_dir, "KEGG"),
  selected_kegg_list = selected_kegg_stim,
  outdir = selected_kegg_stim_dir
)

run_selected_kegg_plot(
  kegg_dir = file.path(enrich_stim_unstim_dir, "KEGG"),
  selected_kegg_list = selected_kegg_stim_unstim,
  outdir = selected_kegg_stim_unstim_dir
)

cat("Selected KEGG pathway visualization finished.\n")

############################################################
## 18. Final message
############################################################

cat("\n========================================\n")
cat("Bulk RNA-seq original analysis finished successfully.\n")
cat("Output directory:\n")
cat(result_dir, "\n\n")

cat("Generated folders:\n")
cat("00_basic_data\n")
cat("00_QC\n")
cat("01_PCA\n")
cat("02_Pearson_top1000\n")
cat("03_DESeq2_DEG\n")
cat("04_pairwise_Volcano\n")
cat("05_Three_group_DEG_heatmap\n")
cat("06_KEGG_GO\n")
cat("07_Selected_pathway\n")

cat("\nMain completed analyses:\n")
cat("1. Basic count matrix processing and metadata matching\n")
cat("2. QC mean ± SEM plots\n")
cat("3. Six-group PCA\n")
cat("4. Top1000 Pearson correlation heatmap\n")
cat("5. Pairwise DESeq2 DEG analysis\n")
cat("6. Three-group overall LRT DEG analysis\n")
cat("7. Pairwise volcano plots\n")
cat("8. Three-group overall top1000 DEG heatmaps\n")
cat("9. KEGG/GO enrichment analysis\n")
cat("10. Selected KEGG pathway visualization\n")
cat("========================================\n")

############################################################
## Save session information
############################################################

writeLines(
  capture.output(sessionInfo()),
  con = file.path(result_dir, "sessionInfo.txt")
)

cat("Session information saved to:\n")
cat(file.path(result_dir, "sessionInfo.txt"), "\n")

