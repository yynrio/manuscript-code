############################################################
############################################################
## HSPC reference comparison demo analysis
## Extended Data Fig. 1a,b
##
## This script performs the demo HSPC reference comparison
## analysis using reduced demo Seurat objects.
##
## Demo data:
## 1. Adult bone marrow HSPCs are downsampled from the full
##    adult bone marrow reference Seurat object.
## 2. Fetal liver HSPCs are downsampled from the full fetal
##    liver reference Seurat object.
## 3. Teratoma-derived HSPCs are downsampled from the full
##    teratoma-derived HSPC Seurat object.
##
## Main analysis steps are kept consistent with the full
## original HSPC reference comparison script.
##
## The demo dataset is intended to test the workflow and
## input/output structure. It is not expected to reproduce
## the exact quantitative results in the manuscript.
##
## Main analysis steps:
## 1. Read HSPC wave gene list
## 2. Read demo adult bone marrow HSPC Seurat object
## 3. Read demo fetal liver HSPC Seurat object
## 4. Read demo teratoma-derived HSPC Seurat object
## 5. Assign HSPC reference group labels
## 6. Generate HSPC wave gene dot plot
## 7. Generate pseudo-bulk expression profiles
## 8. Perform unsupervised hierarchical clustering
##
## Figure outputs:
## Extended Data Fig. 1a: HSPC wave gene dot plot
## Extended Data Fig. 1b: pseudo-bulk hierarchical clustering
############################################################
############################################################

rm(list = ls())
gc()

############################################################
## 0. Load packages
############################################################

packages_cran <- c(
  "ggplot2"
)

packages_other <- c(
  "Seurat",
  "Matrix"
)

for (pkg in packages_cran) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

for (pkg in packages_other) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(Matrix)
})

options(stringsAsFactors = FALSE)

############################################################
## 0.1. Input files
############################################################

demo_data_dir <- "demo_data/HSPC_reference_comparison_SFig1_a_b_demo"

hspc_wave_file <- file.path(
  demo_data_dir,
  "HSPC_waves.txt"
)

bm_hspc_rds_file <- file.path(
  demo_data_dir,
  "GSE253355_BM_HSPC_demo.rds"
)

fetal_liver_rds_file <- file.path(
  demo_data_dir,
  "GSE162950_fetal_liver_demo.rds"
)

teratoma_hspc_rds_file <- file.path(
  demo_data_dir,
  "Tera_integrate_Label_demo.rds"
)

############################################################
## 1. Output directories
############################################################

result_dir <- "results/HSPC_reference_comparison_SFig1_a_b_demo"

basic_dir <- file.path(result_dir, "00_basic_data")
fig_a_dir <- file.path(result_dir, "01_ExtendedData_Fig1a_HSPC_waves_dotplot")
fig_b_dir <- file.path(result_dir, "02_ExtendedData_Fig1b_pseudobulk_hclust")

for (d in c(
  result_dir,
  basic_dir,
  fig_a_dir,
  fig_b_dir
)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

cat("Demo input directory:\n")
cat(demo_data_dir, "\n\n")

cat("Demo output directory:\n")
cat(result_dir, "\n\n")

############################################################
## 2. Check input files
############################################################

input_files <- c(
  hspc_wave_file,
  bm_hspc_rds_file,
  fetal_liver_rds_file,
  teratoma_hspc_rds_file
)

if (!all(file.exists(input_files))) {
  cat("Input file check:\n")
  print(data.frame(
    file = input_files,
    exists = file.exists(input_files)
  ))
  stop("Some input files are missing. Please check demo_data_dir.")
}

############################################################
## 3. Helper functions
############################################################

get_rna_data_matrix <- function(obj) {
  
  DefaultAssay(obj) <- "RNA"
  
  mat <- tryCatch(
    {
      GetAssayData(obj, assay = "RNA", layer = "data")
    },
    error = function(e) {
      GetAssayData(obj, assay = "RNA", slot = "data")
    }
  )
  
  return(mat)
}

make_pseudo_bulk <- function(obj, group_col, groups_keep) {
  
  DefaultAssay(obj) <- "RNA"
  
  expr_mat <- get_rna_data_matrix(obj)
  meta_use <- obj@meta.data
  
  if (!(group_col %in% colnames(meta_use))) {
    stop(paste0("group_col not found in metadata: ", group_col))
  }
  
  pseudo_list <- list()
  
  for (g in groups_keep) {
    
    cells_g <- rownames(meta_use)[
      which(as.character(meta_use[[group_col]]) == g)
    ]
    
    cells_g <- intersect(cells_g, colnames(expr_mat))
    
    cat("Pseudo-bulk group:", g, " cells:", length(cells_g), "\n")
    
    if (length(cells_g) == 0) {
      warning(paste0("No cells found for group: ", g))
      next
    }
    
    if (length(cells_g) == 1) {
      pseudo_list[[g]] <- as.numeric(expr_mat[, cells_g])
    } else {
      pseudo_list[[g]] <- Matrix::rowMeans(expr_mat[, cells_g, drop = FALSE])
    }
  }
  
  if (length(pseudo_list) == 0) {
    stop("No pseudo-bulk groups were generated.")
  }
  
  pseudo_df <- as.data.frame(pseudo_list)
  pseudo_df$geneID <- rownames(expr_mat)
  
  pseudo_df <- pseudo_df[, c("geneID", names(pseudo_list)), drop = FALSE]
  
  return(pseudo_df)
}

############################################################
## 4. Read HSPC wave gene list
############################################################

HSPCwaves_geneSet <- read.table(
  hspc_wave_file,
  header = FALSE,
  sep = "\t",
  stringsAsFactors = FALSE
)

hspc_wave_genes <- as.character(HSPCwaves_geneSet[, 1])
hspc_wave_genes <- trimws(hspc_wave_genes)
hspc_wave_genes <- hspc_wave_genes[hspc_wave_genes != ""]
hspc_wave_genes <- unique(hspc_wave_genes)

write.csv(
  data.frame(gene = hspc_wave_genes),
  file = file.path(basic_dir, "HSPC_waves_gene_list.csv"),
  row.names = FALSE,
  quote = FALSE
)

cat("Number of HSPC wave genes:\n")
print(length(hspc_wave_genes))

############################################################
## 5. Adult bone marrow HSPC reference
############################################################

cat("Reading adult bone marrow HSPC demo object...\n")

seurat.BM <- readRDS(bm_hspc_rds_file)

BM.HSPC.cells <- colnames(seurat.BM)[
  which(
    seurat.BM$cluster_anno_l1 == "HSPC" &
      seurat.BM$orig.ident != "H14_MACS"
  )
]

cat("BM-HSPC cells:\n")
print(length(BM.HSPC.cells))

seurat.BM_HSPC <- subset(
  seurat.BM,
  cells = BM.HSPC.cells
)

seurat.BM_HSPC$group_id <- "BM_HSPC"
seurat.BM_HSPC$group_label <- "BM-HSPC"

seurat.BM_HSPC <- RenameCells(
  seurat.BM_HSPC,
  add.cell.id = "BM"
)

rm(seurat.BM)
gc()

############################################################
## 6. Fetal liver HSPC reference
############################################################

cat("Reading fetal liver HSPC demo object...\n")

seurat.Liver <- readRDS(fetal_liver_rds_file)

Liver.CS17 <- colnames(seurat.Liver)[
  which(seurat.Liver$orig.ident == "liver-6wk-563")
]

Liver.8w <- colnames(seurat.Liver)[
  which(seurat.Liver$orig.ident == "liver-8wk-553")
]

Liver.11w <- colnames(seurat.Liver)[
  which(seurat.Liver$orig.ident == "liver-11wk-569")
]

Liver.15w <- colnames(seurat.Liver)[
  which(seurat.Liver$orig.ident == "liver-15wk-101")
]

cat("FL-HSPC-1 cells, CS17:\n")
print(length(Liver.CS17))

cat("FL-HSPC-2 cells, week 8:\n")
print(length(Liver.8w))

cat("FL-HSPC-3 cells, week 11:\n")
print(length(Liver.11w))

cat("FL-HSPC-4 cells, week 15:\n")
print(length(Liver.15w))

FL.keep.cells <- c(
  Liver.CS17,
  Liver.8w,
  Liver.11w,
  Liver.15w
)

seurat.Liver_HSPC <- subset(
  seurat.Liver,
  cells = FL.keep.cells
)

seurat.Liver_HSPC$group_id <- NA_character_
seurat.Liver_HSPC$group_label <- NA_character_

seurat.Liver_HSPC$group_id[
  seurat.Liver_HSPC$orig.ident == "liver-6wk-563"
] <- "FL_HSPC_1"

seurat.Liver_HSPC$group_label[
  seurat.Liver_HSPC$orig.ident == "liver-6wk-563"
] <- "FL-HSPC-1"

seurat.Liver_HSPC$group_id[
  seurat.Liver_HSPC$orig.ident == "liver-8wk-553"
] <- "FL_HSPC_2"

seurat.Liver_HSPC$group_label[
  seurat.Liver_HSPC$orig.ident == "liver-8wk-553"
] <- "FL-HSPC-2"

seurat.Liver_HSPC$group_id[
  seurat.Liver_HSPC$orig.ident == "liver-11wk-569"
] <- "FL_HSPC_3"

seurat.Liver_HSPC$group_label[
  seurat.Liver_HSPC$orig.ident == "liver-11wk-569"
] <- "FL-HSPC-3"

seurat.Liver_HSPC$group_id[
  seurat.Liver_HSPC$orig.ident == "liver-15wk-101"
] <- "FL_HSPC_4"

seurat.Liver_HSPC$group_label[
  seurat.Liver_HSPC$orig.ident == "liver-15wk-101"
] <- "FL-HSPC-4"

seurat.Liver_HSPC <- RenameCells(
  seurat.Liver_HSPC,
  add.cell.id = "FL"
)

rm(seurat.Liver)
gc()

############################################################
## 7. Teratoma-derived HSPC
############################################################

cat("Reading teratoma-derived HSPC demo object...\n")

seurat.Tera <- readRDS(teratoma_hspc_rds_file)

Tera.HSPC.cells <- colnames(seurat.Tera)[
  which(seurat.Tera$myLabel == "HSPC")
]

cat("Teratoma-HSPC cells:\n")
print(length(Tera.HSPC.cells))

seurat.Tera_HSPC <- subset(
  seurat.Tera,
  cells = Tera.HSPC.cells
)

seurat.Tera_HSPC$group_id <- "Teratoma_HSPC"
seurat.Tera_HSPC$group_label <- "Teratoma-HSPC"

seurat.Tera_HSPC <- RenameCells(
  seurat.Tera_HSPC,
  add.cell.id = "Teratoma"
)

rm(seurat.Tera)
gc()

############################################################
## 8. Combine HSPC objects
############################################################

cat("Merging HSPC reference objects...\n")

combined_HSPC <- merge(
  seurat.BM_HSPC,
  y = c(seurat.Liver_HSPC, seurat.Tera_HSPC)
)

group_id_levels <- c(
  "FL_HSPC_1",
  "FL_HSPC_2",
  "FL_HSPC_3",
  "FL_HSPC_4",
  "BM_HSPC",
  "Teratoma_HSPC"
)

group_label_levels_display <- c(
  "FL-HSPC-1",
  "FL-HSPC-2",
  "FL-HSPC-3",
  "FL-HSPC-4",
  "BM-HSPC",
  "Teratoma-HSPC"
)

combined_HSPC$group_label <- factor(
  combined_HSPC$group_label,
  levels = rev(group_label_levels_display)
)

combined_HSPC$group_id <- factor(
  combined_HSPC$group_id,
  levels = group_id_levels
)

cell_number_table <- as.data.frame(table(combined_HSPC$group_label))
colnames(cell_number_table) <- c("group", "cell_number")

write.csv(
  cell_number_table,
  file = file.path(basic_dir, "HSPC_reference_group_cell_numbers.csv"),
  row.names = FALSE,
  quote = FALSE
)

cat("Final HSPC reference group cell numbers:\n")
print(cell_number_table)

############################################################
## 9. Extended Data Fig. 1a: HSPC wave dot plot
############################################################

cat("Generating Extended Data Fig. 1a dot plot...\n")

DefaultAssay(combined_HSPC) <- "RNA"

genes_found <- hspc_wave_genes[
  hspc_wave_genes %in% rownames(combined_HSPC)
]

genes_missing <- hspc_wave_genes[
  !hspc_wave_genes %in% rownames(combined_HSPC)
]

write.csv(
  data.frame(gene = genes_found),
  file = file.path(basic_dir, "HSPC_waves_genes_found_in_Seurat_object.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  data.frame(gene = genes_missing),
  file = file.path(basic_dir, "HSPC_waves_genes_missing_in_Seurat_object.csv"),
  row.names = FALSE,
  quote = FALSE
)

cat("Found HSPC wave genes:\n")
print(genes_found)

cat("Missing HSPC wave genes:\n")
print(genes_missing)

p_dot <- DotPlot(
  combined_HSPC,
  features = genes_found,
  group.by = "group_label",
  dot.scale = 6,
  cols = c("#FFFFFF", "#E64B35")
) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(
      color = "black",
      face = "bold"
    ),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "top",
    plot.title = element_text(
      hjust = 0,
      face = "bold"
    )
  ) +
  labs(
    title = "HSPC waves"
  )

ggsave(
  filename = file.path(
    fig_a_dir,
    "ExtendedData_Fig1a_HSPC_reference_comparison_HSPC_waves_DotPlot.pdf"
  ),
  plot = p_dot,
  width = 11,
  height = 4.5
)

ggsave(
  filename = file.path(
    fig_a_dir,
    "ExtendedData_Fig1a_HSPC_reference_comparison_HSPC_waves_DotPlot.png"
  ),
  plot = p_dot,
  width = 11,
  height = 4.5,
  dpi = 300
)

print(p_dot)

############################################################
## 10. Pseudo-bulk expression matrix
############################################################

cat("Generating pseudo-bulk expression matrix...\n")

pseudo_BM <- make_pseudo_bulk(
  obj = seurat.BM_HSPC,
  group_col = "group_id",
  groups_keep = "BM_HSPC"
)

pseudo_FL <- make_pseudo_bulk(
  obj = seurat.Liver_HSPC,
  group_col = "group_id",
  groups_keep = c(
    "FL_HSPC_1",
    "FL_HSPC_2",
    "FL_HSPC_3",
    "FL_HSPC_4"
  )
)

pseudo_Tera <- make_pseudo_bulk(
  obj = seurat.Tera_HSPC,
  group_col = "group_id",
  groups_keep = "Teratoma_HSPC"
)

pseudo_all <- Reduce(
  function(x, y) merge(x, y, by = "geneID"),
  list(pseudo_FL, pseudo_BM, pseudo_Tera)
)

pseudo_all <- pseudo_all[
  rowSums(pseudo_all[, -1, drop = FALSE] > 0) > 1,
  ,
  drop = FALSE
]

pseudo_all_log <- log2(pseudo_all[, -1, drop = FALSE] + 1)

pseudo_col_label_map <- c(
  "FL_HSPC_1" = "FL-HSPC-1",
  "FL_HSPC_2" = "FL-HSPC-2",
  "FL_HSPC_3" = "FL-HSPC-3",
  "FL_HSPC_4" = "FL-HSPC-4",
  "BM_HSPC" = "BM-HSPC",
  "Teratoma_HSPC" = "Teratoma-HSPC"
)

colnames(pseudo_all_log) <- pseudo_col_label_map[colnames(pseudo_all_log)]

write.csv(
  pseudo_all,
  file = file.path(basic_dir, "pseudo_bulk_expression_matrix_raw_group_ids.csv"),
  row.names = FALSE,
  quote = FALSE
)

write.csv(
  pseudo_all_log,
  file = file.path(basic_dir, "pseudo_bulk_expression_matrix_log2_display_labels.csv"),
  row.names = TRUE,
  quote = FALSE
)

############################################################
## 11. Extended Data Fig. 1b: pseudo-bulk hclust
############################################################

cat("Generating Extended Data Fig. 1b pseudo-bulk hclust tree...\n")

pseudo_all_log_t <- as.data.frame(t(pseudo_all_log))

dist_pseudo <- dist(
  pseudo_all_log_t,
  method = "euclidean"
)

sampleTree_pseudo <- hclust(
  dist_pseudo,
  method = "average"
)

pdf(
  file = file.path(
    fig_b_dir,
    "ExtendedData_Fig1b_HSPC_reference_comparison_pseudobulk_hclust_tree.pdf"
  ),
  width = 5.5,
  height = 4.5
)

par(mar = c(4, 4, 3, 2))

plot(
  sampleTree_pseudo,
  main = "Pseudo-bulk HSPC reference comparison",
  sub = "",
  xlab = "",
  cex.lab = 1.1,
  cex.axis = 1.0,
  cex.main = 1.1,
  hang = -1
)

dev.off()

png(
  filename = file.path(
    fig_b_dir,
    "ExtendedData_Fig1b_HSPC_reference_comparison_pseudobulk_hclust_tree.png"
  ),
  width = 1600,
  height = 1300,
  res = 300
)

par(mar = c(4, 4, 3, 2))

plot(
  sampleTree_pseudo,
  main = "Pseudo-bulk HSPC reference comparison",
  sub = "",
  xlab = "",
  cex.lab = 1.1,
  cex.axis = 1.0,
  cex.main = 1.1,
  hang = -1
)

dev.off()

write.csv(
  data.frame(
    sample = sampleTree_pseudo$labels,
    stringsAsFactors = FALSE
  ),
  file = file.path(basic_dir, "pseudobulk_hclust_sample_labels.csv"),
  row.names = FALSE,
  quote = FALSE
)

saveRDS(
  sampleTree_pseudo,
  file = file.path(basic_dir, "pseudobulk_hclust_tree.rds")
)

############################################################
## 12. Save session information
############################################################

writeLines(
  capture.output(sessionInfo()),
  con = file.path(result_dir, "sessionInfo.txt")
)

cat("\n========================================\n")
cat("HSPC reference comparison demo analysis finished successfully.\n")
cat("Output directory:\n")
cat(result_dir, "\n\n")

cat("Generated folders:\n")
cat("00_basic_data\n")
cat("01_ExtendedData_Fig1a_HSPC_waves_dotplot\n")
cat("02_ExtendedData_Fig1b_pseudobulk_hclust\n")
cat("========================================\n")