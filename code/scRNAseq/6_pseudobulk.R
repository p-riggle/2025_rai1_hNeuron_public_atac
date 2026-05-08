# pseudobulk method

# specify 

# Load libraries
# Enter commands in R (or R studio, if installed)
library(Seurat)
library(BPCells)
# install.packages("tidyverse")
library(tidyverse)
library(DESeq2)
setwd("/nfs/turbo/umms-siwase/SHIGEKI2/rai1_hesc")

# load the R.data
geo_so = readRDS('data/scRNAseq/iwaselab/rdata/geo_so_sct_unintegrated_with_annos.rds')
unique(geo_so@meta.data$cell_type)

# check the number of cells
geo_so
ncol(geo_so)
length(Cells(geo_so))

# edit metadata

# clustMetaDat <- geo_so@meta.data
# Sample <- paste0(clustMetaDat$orig.ident, ".", clustMetaDat$genotype)
# Clust <- clustMetaDat$seurat_clusters
# aDF <- data.frame(Sample,Clust)
# clustMetaDat <- cbind(clustMetaDat, aDF)
# geo_so@meta.data <- clustMetaDat
# head(clustMetaDat)

head(geo_so@meta.data)
meta_data <- geo_so@meta.data

# aggregate count data by cluster
DefaultAssay(geo_so) <- "RNA"
cts <- AggregateExpression(geo_so, assays = "RNA", slot = "counts", return.seurat = F, group.by =c("orig.ident", "cell_type"))
cts <- cts$RNA
colnames(cts)

# split two time points
cts_d20 <- cts[, c(1:42)]

colnames(cts_d20)

cts_d60 <- cts[, c(43:73)]

colnames(cts_d60)

# save cts objects
saveRDS(cts, "data/scRNAseq/iwaselab/cts")
saveRDS(cts_d20, "data/scRNAseq/iwaselab/cts_d20")
saveRDS(cts_d60, "data/scRNAseq/iwaselab/cts_d60")

# read cts objects
cts <- readRDS("data/scRNAseq/iwaselab/cts")
cts_d20 <- readRDS("data/scRNAseq/iwaselab/cts_d20")
cts_d60 <- readRDS("data/scRNAseq/iwaselab/cts_d60")

# make samplesheet for deseq2
sample_d20 <- sapply(strsplit(colnames(cts_d20), '_'),'[[',1)
cluster_d20 <- sapply(strsplit(colnames(cts_d20), '_'),'[[',2)
sample_d20
cluster_d20
genotype_d20 <- sapply(strsplit(sample_d20, '-'), '[[',2)
genotype_cluster_d20 <- paste(genotype_d20, cluster_d20, sep='_')
genotype_cluster_d20

samplesheet_d20 <- data.frame(sample_d20, genotype_cluster_d20, cluster_d20, genotype_d20, row.names = colnames(cts_d20))
samplesheet_d20

sample_d60 <- sapply(strsplit(colnames(cts_d60), '_'),'[[',1)
cluster_d60 <- sapply(strsplit(colnames(cts_d60), '_'),'[[',2)
sample_d60
cluster_d60
genotype_d60 <- sapply(strsplit(sample_d60, '-'), '[[',2)
genotype_cluster_d60 <- paste(genotype_d60, cluster_d60, sep='_')
genotype_cluster_d60

samplesheet_d60 <- data.frame(sample_d60, genotype_cluster_d60, cluster_d60, genotype_d60, row.names = colnames(cts_d60))
samplesheet_d60

sample_both <- sapply(strsplit(colnames(cts), '_'),'[[',1)
cluster_both <- sapply(strsplit(colnames(cts), '_'),'[[',2)
sample_both
cluster_both
genotype_both <- sapply(strsplit(sample_both, '-'), '[[',2)
genotype_cluster_both <- paste(genotype_both, cluster_both, sep='_')
genotype_cluster_both

samplesheet_both <- data.frame(sample_both, genotype_cluster_both, cluster_both, genotype_both, row.names = colnames(cts))
head(samplesheet_both)


# calculation
dds_d20 <- DESeqDataSetFromMatrix(countData = cts_d20, colData = samplesheet_d20, design = ~ genotype_cluster_d20)
dds_d20 <- DESeq(dds_d20)

dds_d60 <- DESeqDataSetFromMatrix(countData = cts_d60, colData = samplesheet_d60, design = ~ genotype_cluster_d60)
dds_d60 <- DESeq(dds_d60)

dds_both <- DESeqDataSetFromMatrix(countData = cts, colData = samplesheet_both, design = ~ genotype_cluster_both)

# save the dds r object
saveRDS(dds_d20, "data/scRNAseq/iwaselab/dds_d20")
saveRDS(dds_d60, "data/scRNAseq/iwaselab/dds_d60")
saveRDS(dds_both, "data/scRNAseq/iwaselab/dds_both")

##############################################
## for local project directory, start here. ##
##############################################
library(tidyverse)
library(DESeq2)

# read dds r objects

dds_d20 <- readRDS("data/scRNAseq/iwaselab/dds_d20")
dds_d60 <- readRDS("data/scRNAseq/iwaselab/dds_d60")
dds_both <- readRDS("data/scRNAseq/iwaselab/dds_both")

# PCAplots
## transform data
vsd_dds_d20 <- vst(dds_d20, blind = TRUE)  
vsd_dds_d60 <- vst(dds_d60, blind = TRUE)  
vsd_dds_both <- vst(dds_both, blind = TRUE)

head(assay(vsd_dds_d20), 3)
head(assay(vsd_dds_d60), 3)
head(assay(vsd_dds_both), 3)

#install.packages("ggplot2")
library("ggplot2")
str(vsd_dds_d20)
head(vsd_dds_d20)

# Generate the plot object (this is what you already had)
PCA_d20 <- plotPCA(vsd_dds_d20, intgroup = c("genotype_d20", "cluster_d20"))
str(PCA_d20)

# Print the plot object to the open Cairo device
print(PCA_d20)

dev.off()

# filter the PCA data to include only clusters 0,2,3,4
head(PCA_d20$data)
PCA_d20_data <- PCA_d20$data
PCA_d20_data_filtered <- PCA_d20_data[PCA_d20_data$cluster_d20 %in% c("NPC-G1", "NPC-G2/M", "Mes-WT", "Mes-KO", "IPC"), ]
PCA_d20_data_filtered

# create a shape and color mapping for the genotype and day
# shape mapping for genotypes
shape_mapping <- c("wt" = 1, "het" = 2, "ko" = 4)

# color mapping using the assigned names

color_mapping_d20 <- c(
  "NPC-G2/M" = "hotpink",
  "Mes-WT" = "forestgreen",
  "Mes-KO"  = "yellowgreen",
  "NPC-G1" = "purple",
  "IPC"        = "gold"
)

# create a new plot with the modified aesthetics
# extract percent variance explained for PC1 and PC2
percentVar_d20 <- attr(PCA_d20_data, "percentVar")
x_lab <- paste0("PC1 (", round(percentVar_d20[1] * 100, 1), "%)")
y_lab <- paste0("PC2 (", round(percentVar_d20[2] * 100, 1), "%)")

pca_d20_2 <- ggplot(PCA_d20_data_filtered, aes(x = PC1, y = PC2, shape = genotype_d20, color = cluster_d20)) +
  geom_point(size = 4, stroke = 1) +
  scale_shape_manual(values = shape_mapping) +
  scale_color_manual(values = color_mapping_d20) +
  labs(title = "PCA Plot: Day 20", x = x_lab, y = y_lab) +
  theme(text = element_text(family = "Helvetica", size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8)) +
  theme_classic()

pca_d20_2

# save the plot as a PDF file (without using Cairo)
ggsave(filename = "results/scRNAseq/7_pseudobulk/PCA_d20_plot_2.pdf", plot = pca_d20_2, width = 800/100, height = 600/100, units = "in", dpi = 100)

# day 60 pca
PCA_d60 <- plotPCA(vsd_dds_d60, intgroup = c("genotype_d60", "cluster_d60"))
PCA_d60_data <- PCA_d60$data

# Filter the PCA data to include only clusters 1 and 5
PCA_d60_data_filtered <- PCA_d60_data[PCA_d60_data$cluster_d60 %in% c("Neuron", "gRG"), ]
PCA_d60_data_filtered

# create a shape and color mapping for the genotype and day
shape_mapping <- c("wt" = 1, "het" = 2, "ko" = 4)
# color mapping using the assigned names
color_mapping_d60 <- c(
  "Neuron"    = "dodgerblue",
  "gRG" = "orangered"
)

# create pca plot of d60_filtered
# extract percent variance explained for PC1 and PC2
percentVar_d60 <- attr(PCA_d60_data_filtered, "percentVar")
x_lab_d60 <- paste0("PC1 (", round(percentVar_d60[1] * 100, 1), "%)")
y_lab_d60 <- paste0("PC2 (", round(percentVar_d60[2] * 100, 1), "%)")

pca_d60_2 <- ggplot(PCA_d60_data_filtered, aes(x = PC1, y = PC2, shape = genotype_d60, color = cluster_d60)) +
  geom_point(size = 4, stroke = 2) +
  scale_shape_manual(values = shape_mapping) +
  scale_color_manual(values = color_mapping_d60) +
  labs(title = "PCA Plot: Day 60", x = x_lab_d60, y = y_lab_d60) +
  theme(text = element_text(family = "Helvetica", size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8)) +
  theme_classic()

pca_d60_2

# save the plot as pdf (without using Cairo)
ggsave(filename = "results/scRNAseq/7_pseudobulk/PCA_d60_plot_2.pdf", plot = pca_d60_2, width = 800/100, height = 600/100, units = "in", dpi = 100)

# make pca plot of PCA_both
vsd_dds_both$genotype_both
PCA_both <- plotPCA(vsd_dds_both, intgroup = c("genotype_both", "cluster_both"))
PCA_both

PCA_both_data <- PCA_both$data
PCA_both_data

# extract day infromation from the name column and make it as a new column 'day'
PCA_both_data$day <- ifelse(grepl("day20", PCA_both_data$name), "day20", "day60")
PCA_both_data

# create PCA plot of PCA_both_data

# create a shape and color mapping for the genotype and day
shape_mapping <- c("wt" = 1, "het" = 2, "ko" = 4)

color_mapping_both <- c(
  "NPC-G2/M" = "hotpink",
  "Neuron"          = "dodgerblue",
  "Mes-WT" = "forestgreen",
  "Mes-KO"  = "yellowgreen",
  "NPC-G1" = "purple",
  "gRG"       = "orangered",
  "IPC"        = "gold"
)

# make the shape size dependent on days
shape_size_mapping <- c("day20" = 3, "day60" = 5)

# extract percent variance explained for PC1 and PC2
percentVar_both <- attr(PCA_both_data, "percentVar")
x_lab_both <- paste0("PC1 (", round(percentVar_both[1] * 100, 1), "%)")
y_lab_both <- paste0("PC2 (", round(percentVar_both[2] * 100, 1), "%)")

pca_both <- ggplot(PCA_both_data, aes(x = PC1, y = PC2, shape = genotype_both, color = cluster_both, size = day)) +
  geom_point(stroke = 2) +
  scale_shape_manual(values = shape_mapping) +
  scale_color_manual(values = color_mapping_both) +
  scale_size_manual(values = shape_size_mapping) +
  labs(title = "PCA Plot: Day 20&60", x = x_lab_both, y = y_lab_both) +
  theme(text = element_text(family = "Helvetica", size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8)) +
  theme_classic()

pca_both

# save the PCA plot as pdf
ggsave(filename = "results/scRNAseq/7_pseudobulk/PCA_both_plot.pdf", plot = pca_both, width = 800/100, height = 600/100, units = "in", dpi = 100)

dds_d20$genotype_cluster_d20

# make het restables
rstbl_NPCG2M_wt.vs.het <- results(dds_d20, contrast = c("genotype_cluster_d20", "het_NPC-G2/M", "wt_NPC-G2/M"), alpha = 0.01)
rstbl_MesWT_wt.vs.het <- results(dds_d20, contrast = c("genotype_cluster_d20", "het_Mes-WT", "wt_Mes-WT"), alpha = 0.01)
rstbl_MesKO_wt.vs.het <- results(dds_d20, contrast = c("genotype_cluster_d20", "het_Mes-KO", "wt_Mes-KO"), alpha = 0.01)
rstbl_NPCG1_wt.vs.het <- results(dds_d20, contrast = c("genotype_cluster_d20", "het_NPC-G1", "wt_NPC-G1"), alpha = 0.01)
rstbl_IPC_wt.vs.het <- results(dds_d20, contrast = c("genotype_cluster_d20", "het_IPC", "wt_IPC"), alpha = 0.01)
rstbl_Neuron_wt.vs.het <- results(dds_d60, contrast = c("genotype_cluster_d60", "het_Neuron", "wt_Neuron"), alpha = 0.01)
rstbl_gRG_wt.vs.het <- results(dds_d60, contrast = c("genotype_cluster_d60", "het_gRG", "wt_gRG"), alpha = 0.01)

summary(rstbl_NPCG2M_wt.vs.het)
summary(rstbl_MesWT_wt.vs.het)
summary(rstbl_MesKO_wt.vs.het)
summary(rstbl_NPCG1_wt.vs.het)
summary(rstbl_IPC_wt.vs.het)
summary(rstbl_Neuron_wt.vs.het)
summary(rstbl_gRG_wt.vs.het)

# make ko restables
rstbl_NPCG2M_wt.vs.ko <- results(dds_d20, contrast = c("genotype_cluster_d20", "ko_NPC-G2/M", "wt_NPC-G2/M"), alpha = 0.01)
rstbl_MesWT_wt.vs.ko <- results(dds_d20, contrast = c("genotype_cluster_d20", "ko_Mes-WT", "wt_Mes-WT"), alpha = 0.01)
rstbl_MesKO_wt.vs.ko <- results(dds_d20, contrast = c("genotype_cluster_d20", "ko_Mes-KO", "wt_Mes-KO"), alpha = 0.01)
rstbl_NPCG1_wt.vs.ko <- results(dds_d20, contrast = c("genotype_cluster_d20", "ko_NPC-G1", "wt_NPC-G1"), alpha = 0.01)
rstbl_IPC_wt.vs.ko <- results(dds_d20, contrast = c("genotype_cluster_d20", "ko_IPC", "wt_IPC"), alpha = 0.01)
rstbl_Neuron_wt.vs.ko <- results(dds_d60, contrast = c("genotype_cluster_d60", "ko_Neuron", "wt_Neuron"), alpha = 0.01)
rstbl_gRG_wt.vs.ko <- results(dds_d60, contrast = c("genotype_cluster_d60", "ko_gRG", "wt_gRG"), alpha = 0.01)

summary(rstbl_NPCG2M_wt.vs.ko)
summary(rstbl_MesWT_wt.vs.ko)
summary(rstbl_MesKO_wt.vs.ko)
summary(rstbl_NPCG1_wt.vs.ko)
summary(rstbl_IPC_wt.vs.ko)
summary(rstbl_Neuron_wt.vs.ko)
summary(rstbl_gRG_wt.vs.ko)

# write all rstbles to the results directory
# 1. Define the parts that change in the variable/file names
cell_types <- c(
  "NPCG2M", "MesWT", "MesKO", "NPCG1", 
  "IPC", "Neuron", "gRG"
)

comparisons <- c("wt.vs.het", "wt.vs.ko")

# 2. Use nested loops to iterate through all combinations
for (comp in comparisons) {
  for (cell in cell_types) {
    
    # A. Construct the R object name (e.g., "rstbl_NPCG2M_wt.vs.het")
    object_name <- paste0("rstbl_", cell, "_", comp)
    
    # B. Construct the output file path
    # The file name is identical to the object name plus the path and extension
    file_path <- paste0("results/scRNAseq/7_pseudobulk/", object_name, ".csv")
    
    # C. Get the data object using get()
    # This assumes the object (e.g., rstbl_NPCG2M_wt.vs.het) exists in your global environment.
    data_to_write <- get(object_name)
    
    # D. Write the file
    write.csv(as.data.frame(data_to_write), file = file_path)
    
    # Print a confirmation (optional)
    cat("Saved:", object_name, "to", file_path, "\n")
  }
}


# read restbls from the results directory
# 1. Define the parts that change in the file/variable names
cell_types <- c(
  "NPCG2M", "MesWT", "MesKO", "NPCG1", 
  "IPC", "Neuron", "gRG"
)

comparisons <- c("wt.vs.het", "wt.vs.ko")

# 2. Use nested loops to iterate through all combinations
for (comp in comparisons) {
  for (cell in cell_types) {
    
    # A. Construct the R object name (e.g., "rstbl_NPCG2M_wt.vs.het")
    object_name <- paste0("rstbl_", cell, "_", comp)
    
    # B. Construct the input file path
    file_path <- paste0("results/scRNAseq/7_pseudobulk/", object_name, ".csv")
    
    # C. Read the file
    data_to_assign <- read.csv(file_path, row.names = 1)
    
    # D. Assign the data to the dynamically created variable name in the environment
    # The assign() function is necessary to use a variable name stored in a string (object_name)
    assign(x = object_name, value = data_to_assign)
    
    # Print a confirmation (optional)
    cat("Read:", file_path, "into object:", object_name, "\n")
  }
}


# save the restables of each cluster as excel spreadsheet with up- and down-regulated genes in separate tabs
# install.packages("openxlsx")
library(openxlsx)

# Create a new workbook
wb <- createWorkbook()

# List of result tables for HET vs WT using meaningful cell type names
res_list_het <- list(
  "NPC-G2/M" = rstbl_NPCG2M_wt.vs.het,
  "Neuron" = rstbl_Neuron_wt.vs.het,
  "Mes-WT" = rstbl_MesWT_wt.vs.het,
  "Mes-Mut" = rstbl_MesKO_wt.vs.het,
  "NPC-G1" = rstbl_NPCG1_wt.vs.het,
  "gRG" = rstbl_gRG_wt.vs.het,
  "IPC" = rstbl_IPC_wt.vs.het
)

# Add worksheets for each cell type with up- and down-regulated genes
for (cell_type in names(res_list_het)) {
  # Clean cell type name for worksheet (remove special characters)
  clean_name <- gsub("[^A-Za-z0-9]", "_", cell_type)
  addWorksheet(wb, paste0(clean_name, "_Up"))
  addWorksheet(wb, paste0(clean_name, "_Down"))

  res <- res_list_het[[cell_type]]

  # Extract up- and down-regulated genes
  up_genes <- subset(res, padj < 0.01 & log2FoldChange > 1)
  down_genes <- subset(res, padj < 0.01 & log2FoldChange < -1)

  # Add gene names as a column and sort by absolute log2FoldChange
  if (nrow(up_genes) > 0) {
    up_genes$gene <- rownames(up_genes)
    up_genes <- up_genes[, c("gene", setdiff(names(up_genes), "gene"))]  # Move gene column to first position
    # Sort by absolute log2FoldChange (descending - highest fold changes first)
    up_genes <- up_genes[order(abs(up_genes$log2FoldChange), decreasing = TRUE), ]
  }
  if (nrow(down_genes) > 0) {
    down_genes$gene <- rownames(down_genes)
    down_genes <- down_genes[, c("gene", setdiff(names(down_genes), "gene"))]  # Move gene column to first position
    # Sort by absolute log2FoldChange (descending - highest fold changes first)
    down_genes <- down_genes[order(abs(down_genes$log2FoldChange), decreasing = TRUE), ]
  }

  # Write data to the respective worksheets
  writeData(wb, sheet = paste0(clean_name, "_Up"), up_genes)
  writeData(wb, sheet = paste0(clean_name, "_Down"), down_genes)
}

# Save the workbook
saveWorkbook(wb, file = "results/scRNAseq/7_pseudobulk/degs_het_vs_wt.xlsx", overwrite = TRUE)

# make similar workbook for KO vs WT
library(openxlsx)

# Create a new workbook
wb_ko <- createWorkbook()

# extract and save data for KO vs WT
res_list_ko <- list(
  "NPC-G2/M" = rstbl_NPCG2M_wt.vs.ko,
  "Neuron" = rstbl_Neuron_wt.vs.ko,
  "Mes-WT" = rstbl_MesWT_wt.vs.ko,
  "Mes-Mut" = rstbl_MesKO_wt.vs.ko,
  "NPC-G1" = rstbl_NPCG1_wt.vs.ko,
  "gRG" = rstbl_gRG_wt.vs.ko,
  "IPC" = rstbl_IPC_wt.vs.ko
)

# Add worksheets for each cell type with up- and down-regulated genes
for (cell_type in names(res_list_ko)) {
  # Clean cell type name for worksheet (remove special characters)
  clean_name <- gsub("[^A-Za-z0-9]", "_", cell_type)
  addWorksheet(wb_ko, paste0(clean_name, "_Up"))
  addWorksheet(wb_ko, paste0(clean_name, "_Down"))
  res <- res_list_ko[[cell_type]]

  # Extract up- and down-regulated genes
  up_genes <- subset(res, padj < 0.01 & log2FoldChange > 1)
  down_genes <- subset(res, padj < 0.01 & log2FoldChange < -1)

  # Add gene names as a column and sort by absolute log2FoldChange
  if (nrow(up_genes) > 0) {
    up_genes$gene <- rownames(up_genes)
    up_genes <- up_genes[, c("gene", setdiff(names(up_genes), "gene"))]  # Move gene column to first position
    # Sort by absolute log2FoldChange (descending - highest fold changes first)
    up_genes <- up_genes[order(abs(up_genes$log2FoldChange), decreasing = TRUE), ]
  }
  if (nrow(down_genes) > 0) {
    down_genes$gene <- rownames(down_genes)
    down_genes <- down_genes[, c("gene", setdiff(names(down_genes), "gene"))]  # Move gene column to first position
    # Sort by absolute log2FoldChange (descending - highest fold changes first)
    down_genes <- down_genes[order(abs(down_genes$log2FoldChange), decreasing = TRUE), ]
  }

  # Write data to the respective worksheets
  writeData(wb, sheet = paste0(clean_name, "_Up"), up_genes)
  writeData(wb, sheet = paste0(clean_name, "_Down"), down_genes)
}

# Save the workbook
saveWorkbook(wb, file = "results/scRNAseq/4_pseudobulk/degs_ko_vs_wt.xlsx", overwrite = TRUE)

# Load ggplot2 for plotting

library(ggplot2)

# Create the histogram (bar chart) of DEGs of HET vs WT
deg_counts_het <- data.frame(
  Cell_Type = c("NPC-G2/M", "Mes-WT", "Mes-KO", "NPC-G1", "IPC", "Neuron", "gRG"),
  Upregulated = c(
    sum(rstbl_NPCG2M_wt.vs.het$padj < 0.01 & rstbl_NPCG2M_wt.vs.het$log2FoldChange > 1, na.rm = TRUE),
    sum(rstbl_MesWT_wt.vs.het$padj < 0.01 & rstbl_MesWT_wt.vs.het$log2FoldChange > 1, na.rm = TRUE),
    sum(rstbl_MesKO_wt.vs.het$padj < 0.01 & rstbl_MesKO_wt.vs.het$log2FoldChange > 1, na.rm = TRUE),
    sum(rstbl_NPCG1_wt.vs.het$padj < 0.01 & rstbl_NPCG1_wt.vs.het$log2FoldChange > 1, na.rm = TRUE),
    sum(rstbl_IPC_wt.vs.het$padj < 0.01 & rstbl_IPC_wt.vs.het$log2FoldChange > 1, na.rm = TRUE),
    sum(rstbl_Neuron_wt.vs.het$padj < 0.01 & rstbl_Neuron_wt.vs.het$log2FoldChange > 1, na.rm = TRUE),
    sum(rstbl_gRG_wt.vs.het$padj < 0.01 & rstbl_gRG_wt.vs.het$log2FoldChange > 1, na.rm = TRUE)
  ),
  Downregulated = c(
    sum(rstbl_NPCG2M_wt.vs.het$padj < 0.01 & rstbl_NPCG2M_wt.vs.het$log2FoldChange < -1, na.rm = TRUE),
    sum(rstbl_MesWT_wt.vs.het$padj < 0.01 & rstbl_MesWT_wt.vs.het$log2FoldChange < -1, na.rm = TRUE),
    sum(rstbl_MesKO_wt.vs.het$padj < 0.01 & rstbl_MesKO_wt.vs.het$log2FoldChange < -1, na.rm = TRUE),
    sum(rstbl_NPCG1_wt.vs.het$padj < 0.01 & rstbl_NPCG1_wt.vs.het$log2FoldChange < -1, na.rm = TRUE),
    sum(rstbl_IPC_wt.vs.het$padj < 0.01 & rstbl_IPC_wt.vs.het$log2FoldChange < -1, na.rm = TRUE),
    sum(rstbl_Neuron_wt.vs.het$padj < 0.01 & rstbl_Neuron_wt.vs.het$log2FoldChange < -1, na.rm = TRUE),
    sum(rstbl_gRG_wt.vs.het$padj < 0.01 & rstbl_gRG_wt.vs.het$log2FoldChange < -1, na.rm = TRUE)
  )
)

deg_counts_het_melted <- reshape2::melt(deg_counts_het, id.vars = "Cell_Type", variable.name = "Regulation", value.name = "Count")

# Plot the bar chart
p_het <- ggplot(deg_counts_het_melted, aes(x = Cell_Type, y = Count, fill = Regulation)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "Number of DEGs in HET vs WT by Cell Type", x = "Cell Type", y = "Number of DEGs") +
  theme_minimal() +
  scale_fill_manual(values = c("Upregulated" = "steelblue", "Downregulated" = "salmon")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p_het

#############################################################################################
### go analysis of DEGs of HET vs WT in each cell type using clusterProfiler and enrichGO ###
#############################################################################################

# install.packages(c("clusterProfiler", "enrichplot", "org.Hs.eg.db"))
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(tidyverse)


# read restbls from the results directoryplot
# 1. Define the parts that change in the file/variable names
cell_types <- c(
  "NPCG2M", "MesWT", "MesKO", "NPCG1", 
  "IPC", "Neuron", "gRG"
)

comparisons <- c("wt.vs.het", "wt.vs.ko")

# 2. Use nested loops to iterate through all combinations
for (comp in comparisons) {
  for (cell in cell_types) {
    
    # A. Construct the R object name (e.g., "rstbl_NPCG2M_wt.vs.het")
    object_name <- paste0("rstbl_", cell, "_", comp)
    
    # B. Construct the input file path
    file_path <- paste0("results/scRNAseq/7_pseudobulk/", object_name, ".csv")
    
    # C. Read the file
    data_to_assign <- read.csv(file_path, row.names = 1)
    
    # D. Assign the data to the dynamically created variable name in the environment
    # The assign() function is necessary to use a variable name stored in a string (object_name)
    assign(x = object_name, value = data_to_assign)
    
    # Print a confirmation (optional)
    cat("Read:", file_path, "into object:", object_name, "\n")
  }
}

# output directory for GO analysis
go_outdir <- "results/scRNAseq/7_pseudobulk/go_bp_het_vs_wt"
if (!dir.exists(go_outdir)) {
  dir.create(go_outdir, recursive = TRUE)
}

# use HET vs WT result tables already created above
res_list_het <- list(
  "NPC-G2/M" = rstbl_NPCG2M_wt.vs.het,
  "NPC-G1" = rstbl_NPCG1_wt.vs.het,
  "Mes-WT" = rstbl_MesWT_wt.vs.het,
  "IPC" = rstbl_IPC_wt.vs.het,
  "Neuron" = rstbl_Neuron_wt.vs.het,
  "gRG" = rstbl_gRG_wt.vs.het,
  "Mes-Mut" = rstbl_MesKO_wt.vs.het
)

# explicit cell type order for big dotplots (exclude Mes-Mut)
all_cell_types <- c("NPC-G2/M", "NPC-G1", "Mes-WT", "IPC", "Neuron", "gRG")

# keep only selected cell types in requested order
res_list_het <- res_list_het[all_cell_types]

# DEG thresholds for GO input
padj_cutoff <- 0.01
lfc_cutoff <- 1

# helper: run enrichGO for one gene vector using gene-specific universe
run_enrichgo <- function(sig_symbols, universe_symbols, ontology = "BP") {
  if (length(sig_symbols) == 0) {
    return(NULL)
  }

  gene_map <- bitr(sig_symbols,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)

  universe_map <- bitr(universe_symbols,
                       fromType = "SYMBOL",
                       toType = "ENTREZID",
                       OrgDb = org.Hs.eg.db)

  if (is.null(gene_map) || nrow(gene_map) == 0 || is.null(universe_map) || nrow(universe_map) == 0) {
    return(NULL)
  }

  ego <- enrichGO(
    gene = unique(gene_map$ENTREZID),
    universe = unique(universe_map$ENTREZID),
    OrgDb = org.Hs.eg.db,
    keyType = "ENTREZID",
    ont = ontology,
    pAdjustMethod = "BH",
    pvalueCutoff = 0.05,
    qvalueCutoff = 0.2,
    minGSSize = 10,
    maxGSSize = 299,
    readable = TRUE
  )

  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    return(NULL)
  }

  ego
}

# helper: parse GeneRatio string "x/y" into numeric value
parse_gene_ratio <- function(gene_ratio_vec) {
  sapply(gene_ratio_vec, function(x) {
    if (is.na(x) || !grepl("/", x)) return(NA_real_)
    parts <- as.numeric(strsplit(x, "/")[[1]])
    if (length(parts) != 2 || any(is.na(parts)) || parts[2] == 0) return(NA_real_)
    parts[1] / parts[2]
  })
}

# helper: one-line top gene label (top five by abs(log2FoldChange))
format_top_genes_one_line <- function(gene_str, res_df, n_genes = 5) {
  genes <- unlist(strsplit(gene_str, "/"))
  sub <- res_df[res_df$gene %in% genes & !is.na(res_df$log2FoldChange), ]
  if (nrow(sub) == 0) {
    return(paste(head(genes, n_genes), collapse = ", "))
  }
  sub <- sub[order(abs(sub$log2FoldChange), decreasing = TRUE), ]
  paste(head(sub$gene, n_genes), collapse = ", ")
}

# helper: create one regulation-specific dotplot with vertical cell-type alignment
# each plot keeps GO terms on the y-axis
make_regulation_go_dotplot <- function(go_tbl, all_cell_types, title_text, regulation = c("Up", "Down")) {
  regulation <- match.arg(regulation)

  if (is.null(go_tbl) || nrow(go_tbl) == 0) {
    return(NULL)
  }

  # Check which cell types actually have data rows BEFORE any factor conversion.
  present_types <- unique(go_tbl$cell_type)          # character at this point
  missing_types <- setdiff(all_cell_types, present_types)

  # Inject one placeholder row per missing cell type so facet_grid renders
  # an empty panel for them. Do this BEFORE building term_id factors.
  if (length(missing_types) > 0) {
    padding <- data.frame(
      cell_type        = missing_types,
      GeneRatio        = NA_character_,
      p.adjust         = NA_real_,
      Description      = NA_character_,
      top_genes        = NA_character_,
      stringsAsFactors = FALSE
    )
    go_tbl <- bind_rows(go_tbl, padding)
  }

  go_tbl <- go_tbl %>%
    mutate(
      enrich_ratio   = parse_gene_ratio(GeneRatio),
      neglog10_padj  = -log10(p.adjust),
      adj_p_for_plot = pmax(p.adjust, .Machine$double.xmin),
      cell_type      = factor(cell_type, levels = all_cell_types)
    ) %>%
    group_by(cell_type) %>%
    arrange(p.adjust, Description, .by_group = TRUE) %>%
    mutate(term_rank = row_number()) %>%
    ungroup() %>%
    mutate(
      # Placeholder rows get a stable unique term_id; real rows encode rank.
      term_id = ifelse(
        is.na(Description),
        paste(as.character(cell_type), "99", "no enrichment", sep = "__"),
        paste(as.character(cell_type), sprintf("%02d", term_rank), Description, sep = "__")
      ),
      x_col = "GO"
    )

  # Build y-axis factor levels: low p-value (rank 1) at top, high p-value at
  # bottom, placeholder rows last within each cell type.
  term_levels <- go_tbl %>%
    arrange(cell_type, desc(is.na(term_rank)), desc(term_rank)) %>%
    pull(term_id) %>%
    unique()
  go_tbl$term_id <- factor(go_tbl$term_id, levels = term_levels)


  color_low <- if (regulation == "Down") "#08519c" else "#d2b005"
  color_high <- if (regulation == "Down") "#deebf7" else "#fdeebf"

  p <- ggplot(go_tbl, aes(x = x_col, y = term_id))
  p <- p + geom_point(aes(size = enrich_ratio, color = adj_p_for_plot), alpha = 0.95)
  p <- p + geom_text(
    aes(label = top_genes),
    nudge_x = 0.08,
    hjust = 0,
    size = 2.0,
    check_overlap = TRUE
  )
  p <- p + scale_color_gradient(
    name = "Adjusted p-value",
    low = color_low,
    high = color_high,
    trans = "log10",
    labels = scales::label_scientific(digits = 2)
  )
  p <- p + scale_size_continuous(name = "Enrichment ratio", range = c(2, 8))
  p <- p + scale_x_discrete(drop = FALSE)
  p <- p + facet_grid(rows = vars(cell_type), scales = "free_y", switch = "y", drop = FALSE)
  p <- p + scale_y_discrete(labels = function(x) sub("^.*__[0-9]{2}__", "", x))
  p <- p + guides(
    size = guide_legend(order = 1),
    color = guide_colorbar(order = 2, barheight = grid::unit(3, "cm"))
  )
  p <- p + labs(
    title = title_text,
    x = NULL,
    y = "GO term"
  )
  p <- p + coord_cartesian(clip = "off")
  p <- p + theme_classic()
  p <- p + theme(
    strip.background = element_rect(fill = "grey95", color = "grey70"),
    strip.text = element_text(face = "bold", size = 8),
    strip.placement = "outside",
    axis.text.x = element_text(size = 8),
    axis.text.y = element_text(size = 7),
    legend.position = "bottom",
    legend.box = "vertical",
    panel.spacing = grid::unit(0.6, "lines"),
    plot.margin = margin(5.5, 120, 5.5, 40)
  )

  return(p)
}

go_up_all <- list()
go_down_all <- list()

# run GO for each cell type (up- and down-regulated separately)
for (cell_type in names(res_list_het)) {
  res <- as.data.frame(res_list_het[[cell_type]])
  res$gene <- rownames(res)

  # define significant DEGs
  up_symbols <- res %>%
    filter(!is.na(padj), padj < padj_cutoff, log2FoldChange > lfc_cutoff) %>%
    pull(gene)

  down_symbols <- res %>%
    filter(!is.na(padj), padj < padj_cutoff, log2FoldChange < -lfc_cutoff) %>%
    pull(gene)

  # universe: all genes tested by DESeq2 for this cell type
  universe_symbols <- res %>%
    filter(!is.na(padj)) %>%
    pull(gene)

  clean_name <- gsub("[^A-Za-z0-9]", "_", cell_type)

  # run enrichGO
  ego_up <- run_enrichgo(up_symbols, universe_symbols, ontology = "BP")
  ego_down <- run_enrichgo(down_symbols, universe_symbols, ontology = "BP")

  # write enrichment tables
  if (!is.null(ego_up)) {
    up_tbl <- as.data.frame(ego_up)
    up_tbl$cell_type <- cell_type
    up_tbl$top_genes <- sapply(up_tbl$geneID, format_top_genes_one_line, res_df = res)
    write.csv(up_tbl,
              file = file.path(go_outdir, paste0(clean_name, "_UP_enrichGO_BP.csv")),
              row.names = FALSE)
    go_up_all[[cell_type]] <- up_tbl %>%
      arrange(p.adjust) %>%
      slice_head(n = 5)
  } else {
    message("No significant GO-BP enrichment for UP genes in ", cell_type)
  }

  if (!is.null(ego_down)) {
    down_tbl <- as.data.frame(ego_down)
    down_tbl$cell_type <- cell_type
    down_tbl$top_genes <- sapply(down_tbl$geneID, format_top_genes_one_line, res_df = res)
    write.csv(down_tbl,
              file = file.path(go_outdir, paste0(clean_name, "_DOWN_enrichGO_BP.csv")),
              row.names = FALSE)
    go_down_all[[cell_type]] <- down_tbl %>%
      arrange(p.adjust) %>%
      slice_head(n = 5)
  } else {
    message("No significant GO-BP enrichment for DOWN genes in ", cell_type)
  }
}

# build separate figures for Downregulated and Upregulated GO terms
go_up_tbl <- if (length(go_up_all) > 0) bind_rows(go_up_all) else NULL
go_down_tbl <- if (length(go_down_all) > 0) bind_rows(go_down_all) else NULL

p_down <- make_regulation_go_dotplot(
  go_tbl = go_down_tbl,
  all_cell_types = all_cell_types,
  title_text = "HET vs WT DEGs: GO-BP (Downregulated)",
  regulation = "Down"
)

p_up <- make_regulation_go_dotplot(
  go_tbl = go_up_tbl,
  all_cell_types = all_cell_types,
  title_text = "HET vs WT DEGs: GO-BP (Upregulated)",
  regulation = "Up"
)

if (!is.null(p_down)) {
  n_terms_down <- nrow(go_down_tbl)
  plot_h_down <- max(11, 0.20 * n_terms_down)

  ggsave(
    filename = file.path(go_outdir, "HET_vs_WT_DOWN_GO_BP_dotplot.pdf"),
    plot = p_down,
    width = 10,
    height = plot_h_down,
    limitsize = FALSE
  )
}

if (!is.null(p_up)) {
  n_terms_up <- nrow(go_up_tbl)
  plot_h_up <- max(11, 0.20 * n_terms_up)

  ggsave(
    filename = file.path(go_outdir, "HET_vs_WT_UP_GO_BP_dotplot.pdf"),
    plot = p_up,
    width = 10,
    height = plot_h_up,
    limitsize = FALSE
  )
}

if (!is.null(p_down) && !is.null(p_up)) {
  combined_h <- max(
    max(11, 0.20 * nrow(go_down_tbl)),
    max(11, 0.20 * nrow(go_up_tbl))
  )

  grDevices::pdf(
    file = file.path(go_outdir, "HET_vs_WT_UP_DOWN_GO_BP_dotplot_combined.pdf"),
    width = 20,
    height = combined_h,
    onefile = TRUE
  )

  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))
  print(p_down, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p_up, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grDevices::dev.off()
}

message("GO enrichment analysis (HET vs WT) completed. Results saved to: ", go_outdir)


##################################################################################################################
### go analysis (cellular component) of DEGs of HET vs WT in each cell type using clusterProfiler and enrichGO ###
##################################################################################################################

# output directory for GO-CC analysis
go_cc_outdir <- "results/scRNAseq/7_pseudobulk/go_cc_het_vs_wt"
if (!dir.exists(go_cc_outdir)) {
  dir.create(go_cc_outdir, recursive = TRUE)
}

go_cc_up_all   <- list()
go_cc_down_all <- list()

# run GO-CC for each cell type (up- and down-regulated separately)
for (cell_type in names(res_list_het)) {
  res <- as.data.frame(res_list_het[[cell_type]])
  res$gene <- rownames(res)

  up_symbols <- res %>%
    filter(!is.na(padj), padj < padj_cutoff, log2FoldChange > lfc_cutoff) %>%
    pull(gene)

  down_symbols <- res %>%
    filter(!is.na(padj), padj < padj_cutoff, log2FoldChange < -lfc_cutoff) %>%
    pull(gene)

  universe_symbols <- res %>%
    filter(!is.na(padj)) %>%
    pull(gene)

  clean_name <- gsub("[^A-Za-z0-9]", "_", cell_type)

  ego_up   <- run_enrichgo(up_symbols,   universe_symbols, ontology = "CC")
  ego_down <- run_enrichgo(down_symbols, universe_symbols, ontology = "CC")

  if (!is.null(ego_up)) {
    up_tbl <- as.data.frame(ego_up)
    up_tbl$cell_type  <- cell_type
    up_tbl$top_genes  <- sapply(up_tbl$geneID, format_top_genes_one_line, res_df = res)
    write.csv(up_tbl,
              file = file.path(go_cc_outdir, paste0(clean_name, "_UP_enrichGO_CC.csv")),
              row.names = FALSE)
    go_cc_up_all[[cell_type]] <- up_tbl %>%
      arrange(p.adjust) %>%
      slice_head(n = 5)
  } else {
    message("No significant GO-CC enrichment for UP genes in ", cell_type)
  }

  if (!is.null(ego_down)) {
    down_tbl <- as.data.frame(ego_down)
    down_tbl$cell_type  <- cell_type
    down_tbl$top_genes  <- sapply(down_tbl$geneID, format_top_genes_one_line, res_df = res)
    write.csv(down_tbl,
              file = file.path(go_cc_outdir, paste0(clean_name, "_DOWN_enrichGO_CC.csv")),
              row.names = FALSE)
    go_cc_down_all[[cell_type]] <- down_tbl %>%
      arrange(p.adjust) %>%
      slice_head(n = 5)
  } else {
    message("No significant GO-CC enrichment for DOWN genes in ", cell_type)
  }
}

# build dotplots for CC
go_cc_up_tbl   <- if (length(go_cc_up_all)   > 0) bind_rows(go_cc_up_all)   else NULL
go_cc_down_tbl <- if (length(go_cc_down_all) > 0) bind_rows(go_cc_down_all) else NULL

p_cc_down <- make_regulation_go_dotplot(
  go_tbl      = go_cc_down_tbl,
  all_cell_types = all_cell_types,
  title_text  = "HET vs WT DEGs: GO-CC (Downregulated)",
  regulation  = "Down"
)

p_cc_up <- make_regulation_go_dotplot(
  go_tbl      = go_cc_up_tbl,
  all_cell_types = all_cell_types,
  title_text  = "HET vs WT DEGs: GO-CC (Upregulated)",
  regulation  = "Up"
)

if (!is.null(p_cc_down)) {
  n_terms_cc_down <- nrow(go_cc_down_tbl)
  plot_h_cc_down  <- max(11, 0.20 * n_terms_cc_down)

  ggsave(
    filename  = file.path(go_cc_outdir, "HET_vs_WT_DOWN_GO_CC_dotplot.pdf"),
    plot      = p_cc_down,
    width     = 10,
    height    = plot_h_cc_down,
    limitsize = FALSE
  )
}

if (!is.null(p_cc_up)) {
  n_terms_cc_up <- nrow(go_cc_up_tbl)
  plot_h_cc_up  <- max(11, 0.20 * n_terms_cc_up)

  ggsave(
    filename  = file.path(go_cc_outdir, "HET_vs_WT_UP_GO_CC_dotplot.pdf"),
    plot      = p_cc_up,
    width     = 10,
    height    = plot_h_cc_up,
    limitsize = FALSE
  )
}

if (!is.null(p_cc_down) && !is.null(p_cc_up)) {
  combined_h_cc <- max(
    max(11, 0.20 * nrow(go_cc_down_tbl)),
    max(11, 0.20 * nrow(go_cc_up_tbl))
  )

  grDevices::pdf(
    file    = file.path(go_cc_outdir, "HET_vs_WT_UP_DOWN_GO_CC_dotplot_combined.pdf"),
    width   = 20,
    height  = combined_h_cc,
    onefile = TRUE
  )

  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))
  print(p_cc_down, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p_cc_up,   vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grDevices::dev.off()
}

message("GO-CC enrichment analysis (HET vs WT) completed. Results saved to: ", go_cc_outdir)


##################################################################################################################
### KEGG pathway analysis of DEGs of HET vs WT in each cell type using clusterProfiler and enrichKEGG ###
##################################################################################################################

# output directory
kegg_outdir <- "results/scRNAseq/7_pseudobulk/kegg_het_vs_wt"
if (!dir.exists(kegg_outdir)) {
  dir.create(kegg_outdir, recursive = TRUE)
}

# helper: run enrichKEGG for one gene vector
run_enrichkegg <- function(sig_symbols, universe_symbols) {
  if (length(sig_symbols) == 0) {
    return(NULL)
  }

  gene_map <- bitr(sig_symbols,
                   fromType = "SYMBOL",
                   toType   = "ENTREZID",
                   OrgDb    = org.Hs.eg.db)

  universe_map <- bitr(universe_symbols,
                       fromType = "SYMBOL",
                       toType   = "ENTREZID",
                       OrgDb    = org.Hs.eg.db)

  if (is.null(gene_map) || nrow(gene_map) == 0 ||
      is.null(universe_map) || nrow(universe_map) == 0) {
    return(NULL)
  }

  ek <- enrichKEGG(
    gene          = unique(gene_map$ENTREZID),
    universe      = unique(universe_map$ENTREZID),
    organism      = "hsa",
    pAdjustMethod = "BH",
    pvalueCutoff  = 0.05,
    qvalueCutoff  = 0.2,
    minGSSize     = 10,
    maxGSSize     = 299
  )

  if (is.null(ek) || nrow(as.data.frame(ek)) == 0) {
    return(NULL)
  }

  # convert ENTREZID back to gene symbols for readable labels
  ek <- setReadable(ek, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  ek
}

kegg_up_all   <- list()
kegg_down_all <- list()

for (cell_type in names(res_list_het)) {
  res <- as.data.frame(res_list_het[[cell_type]])
  res$gene <- rownames(res)

  up_symbols <- res %>%
    filter(!is.na(padj), padj < padj_cutoff, log2FoldChange > lfc_cutoff) %>%
    pull(gene)

  down_symbols <- res %>%
    filter(!is.na(padj), padj < padj_cutoff, log2FoldChange < -lfc_cutoff) %>%
    pull(gene)

  universe_symbols <- res %>%
    filter(!is.na(padj)) %>%
    pull(gene)

  clean_name <- gsub("[^A-Za-z0-9]", "_", cell_type)

  ek_up   <- run_enrichkegg(up_symbols,   universe_symbols)
  ek_down <- run_enrichkegg(down_symbols, universe_symbols)

  if (!is.null(ek_up)) {
    up_tbl <- as.data.frame(ek_up)
    up_tbl$cell_type <- cell_type
    up_tbl$top_genes <- sapply(up_tbl$geneID, format_top_genes_one_line, res_df = res)
    write.csv(up_tbl,
              file = file.path(kegg_outdir, paste0(clean_name, "_UP_enrichKEGG.csv")),
              row.names = FALSE)
    kegg_up_all[[cell_type]] <- up_tbl %>%
      arrange(p.adjust) %>%
      slice_head(n = 5)
  } else {
    message("No significant KEGG enrichment for UP genes in ", cell_type)
  }

  if (!is.null(ek_down)) {
    down_tbl <- as.data.frame(ek_down)
    down_tbl$cell_type <- cell_type
    down_tbl$top_genes <- sapply(down_tbl$geneID, format_top_genes_one_line, res_df = res)
    write.csv(down_tbl,
              file = file.path(kegg_outdir, paste0(clean_name, "_DOWN_enrichKEGG.csv")),
              row.names = FALSE)
    kegg_down_all[[cell_type]] <- down_tbl %>%
      arrange(p.adjust) %>%
      slice_head(n = 5)
  } else {
    message("No significant KEGG enrichment for DOWN genes in ", cell_type)
  }
}

kegg_up_tbl   <- if (length(kegg_up_all)   > 0) bind_rows(kegg_up_all)   else NULL
kegg_down_tbl <- if (length(kegg_down_all) > 0) bind_rows(kegg_down_all) else NULL

p_kegg_down <- make_regulation_go_dotplot(
  go_tbl         = kegg_down_tbl,
  all_cell_types = all_cell_types,
  title_text     = "HET vs WT DEGs: KEGG (Downregulated)",
  regulation     = "Down"
)

p_kegg_up <- make_regulation_go_dotplot(
  go_tbl         = kegg_up_tbl,
  all_cell_types = all_cell_types,
  title_text     = "HET vs WT DEGs: KEGG (Upregulated)",
  regulation     = "Up"
)

if (!is.null(p_kegg_down)) {
  n_terms_kegg_down <- nrow(kegg_down_tbl)
  plot_h_kegg_down  <- max(11, 0.20 * n_terms_kegg_down)

  ggsave(
    filename  = file.path(kegg_outdir, "HET_vs_WT_DOWN_KEGG_dotplot.pdf"),
    plot      = p_kegg_down,
    width     = 10,
    height    = plot_h_kegg_down,
    limitsize = FALSE
  )
}

if (!is.null(p_kegg_up)) {
  n_terms_kegg_up <- nrow(kegg_up_tbl)
  plot_h_kegg_up  <- max(11, 0.20 * n_terms_kegg_up)

  ggsave(
    filename  = file.path(kegg_outdir, "HET_vs_WT_UP_KEGG_dotplot.pdf"),
    plot      = p_kegg_up,
    width     = 10,
    height    = plot_h_kegg_up,
    limitsize = FALSE
  )
}

if (!is.null(p_kegg_down) && !is.null(p_kegg_up)) {
  combined_h_kegg <- max(
    max(11, 0.20 * nrow(kegg_down_tbl)),
    max(11, 0.20 * nrow(kegg_up_tbl))
  )

  grDevices::pdf(
    file    = file.path(kegg_outdir, "HET_vs_WT_UP_DOWN_KEGG_dotplot_combined.pdf"),
    width   = 20,
    height  = combined_h_kegg,
    onefile = TRUE
  )

  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))
  print(p_kegg_down, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p_kegg_up,   vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grDevices::dev.off()
}

message("KEGG enrichment analysis (HET vs WT) completed. Results saved to: ", kegg_outdir)


##################################################################################################################
### GO-BP / GO-CC / KEGG analysis of DEGs of WT vs KO in selected cell types (exclude Mes-WT, IPC, NPC-G1; include Mes-Mut)
##################################################################################################################

# use WT vs KO result tables and keep only requested cell types
res_list_ko <- list(
  "NPC-G2/M" = rstbl_NPCG2M_wt.vs.ko,
  "Mes-Mut" = rstbl_MesKO_wt.vs.ko,
  "Neuron" = rstbl_Neuron_wt.vs.ko,
  "gRG" = rstbl_gRG_wt.vs.ko
)

all_cell_types_ko <- c("NPC-G2/M", "Mes-Mut", "Neuron", "gRG")
res_list_ko <- res_list_ko[all_cell_types_ko]


##################################################################################################################
### GO-BP analysis (WT vs KO)
##################################################################################################################

go_bp_ko_outdir <- "results/scRNAseq/7_pseudobulk/go_bp_wt_vs_ko"
if (!dir.exists(go_bp_ko_outdir)) {
  dir.create(go_bp_ko_outdir, recursive = TRUE)
}

go_bp_ko_up_all <- list()
go_bp_ko_down_all <- list()

for (cell_type in names(res_list_ko)) {
  res <- as.data.frame(res_list_ko[[cell_type]])
  res$gene <- rownames(res)

  up_symbols <- res %>%
    filter(!is.na(padj), padj < padj_cutoff, log2FoldChange > lfc_cutoff) %>%
    pull(gene)

  down_symbols <- res %>%
    filter(!is.na(padj), padj < padj_cutoff, log2FoldChange < -lfc_cutoff) %>%
    pull(gene)

  universe_symbols <- res %>%
    filter(!is.na(padj)) %>%
    pull(gene)

  clean_name <- gsub("[^A-Za-z0-9]", "_", cell_type)

  ego_up <- run_enrichgo(up_symbols, universe_symbols, ontology = "BP")
  ego_down <- run_enrichgo(down_symbols, universe_symbols, ontology = "BP")

  if (!is.null(ego_up)) {
    up_tbl <- as.data.frame(ego_up)
    up_tbl$cell_type <- cell_type
    up_tbl$top_genes <- sapply(up_tbl$geneID, format_top_genes_one_line, res_df = res)
    write.csv(up_tbl,
              file = file.path(go_bp_ko_outdir, paste0(clean_name, "_UP_enrichGO_BP.csv")),
              row.names = FALSE)
    go_bp_ko_up_all[[cell_type]] <- up_tbl %>%
      arrange(p.adjust) %>%
      slice_head(n = 5)
  } else {
    message("No significant GO-BP enrichment for UP genes in ", cell_type, " (WT vs KO)")
  }

  if (!is.null(ego_down)) {
    down_tbl <- as.data.frame(ego_down)
    down_tbl$cell_type <- cell_type
    down_tbl$top_genes <- sapply(down_tbl$geneID, format_top_genes_one_line, res_df = res)
    write.csv(down_tbl,
              file = file.path(go_bp_ko_outdir, paste0(clean_name, "_DOWN_enrichGO_BP.csv")),
              row.names = FALSE)
    go_bp_ko_down_all[[cell_type]] <- down_tbl %>%
      arrange(p.adjust) %>%
      slice_head(n = 5)
  } else {
    message("No significant GO-BP enrichment for DOWN genes in ", cell_type, " (WT vs KO)")
  }
}

go_bp_ko_up_tbl <- if (length(go_bp_ko_up_all) > 0) bind_rows(go_bp_ko_up_all) else NULL
go_bp_ko_down_tbl <- if (length(go_bp_ko_down_all) > 0) bind_rows(go_bp_ko_down_all) else NULL

p_bp_ko_down <- make_regulation_go_dotplot(
  go_tbl = go_bp_ko_down_tbl,
  all_cell_types = all_cell_types_ko,
  title_text = "WT vs KO DEGs: GO-BP (Downregulated)",
  regulation = "Down"
)

p_bp_ko_up <- make_regulation_go_dotplot(
  go_tbl = go_bp_ko_up_tbl,
  all_cell_types = all_cell_types_ko,
  title_text = "WT vs KO DEGs: GO-BP (Upregulated)",
  regulation = "Up"
)

if (!is.null(p_bp_ko_down)) {
  n_terms_bp_ko_down <- nrow(go_bp_ko_down_tbl)
  plot_h_bp_ko_down <- max(11, 0.20 * n_terms_bp_ko_down)

  ggsave(
    filename = file.path(go_bp_ko_outdir, "WT_vs_KO_DOWN_GO_BP_dotplot.pdf"),
    plot = p_bp_ko_down,
    width = 10,
    height = plot_h_bp_ko_down,
    limitsize = FALSE
  )
}

if (!is.null(p_bp_ko_up)) {
  n_terms_bp_ko_up <- nrow(go_bp_ko_up_tbl)
  plot_h_bp_ko_up <- max(11, 0.20 * n_terms_bp_ko_up)

  ggsave(
    filename = file.path(go_bp_ko_outdir, "WT_vs_KO_UP_GO_BP_dotplot.pdf"),
    plot = p_bp_ko_up,
    width = 10,
    height = plot_h_bp_ko_up,
    limitsize = FALSE
  )
}

if (!is.null(p_bp_ko_down) && !is.null(p_bp_ko_up)) {
  combined_h_bp_ko <- max(
    max(11, 0.20 * nrow(go_bp_ko_down_tbl)),
    max(11, 0.20 * nrow(go_bp_ko_up_tbl))
  )

  grDevices::pdf(
    file = file.path(go_bp_ko_outdir, "WT_vs_KO_UP_DOWN_GO_BP_dotplot_combined.pdf"),
    width = 20,
    height = combined_h_bp_ko,
    onefile = TRUE
  )

  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))
  print(p_bp_ko_down, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p_bp_ko_up, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grDevices::dev.off()
}

message("GO-BP enrichment analysis (WT vs KO) completed. Results saved to: ", go_bp_ko_outdir)


##################################################################################################################
### GO-CC analysis (WT vs KO)
##################################################################################################################

go_cc_ko_outdir <- "results/scRNAseq/7_pseudobulk/go_cc_wt_vs_ko"
if (!dir.exists(go_cc_ko_outdir)) {
  dir.create(go_cc_ko_outdir, recursive = TRUE)
}

go_cc_ko_up_all <- list()
go_cc_ko_down_all <- list()

for (cell_type in names(res_list_ko)) {
  res <- as.data.frame(res_list_ko[[cell_type]])
  res$gene <- rownames(res)

  up_symbols <- res %>%
    filter(!is.na(padj), padj < padj_cutoff, log2FoldChange > lfc_cutoff) %>%
    pull(gene)

  down_symbols <- res %>%
    filter(!is.na(padj), padj < padj_cutoff, log2FoldChange < -lfc_cutoff) %>%
    pull(gene)

  universe_symbols <- res %>%
    filter(!is.na(padj)) %>%
    pull(gene)

  clean_name <- gsub("[^A-Za-z0-9]", "_", cell_type)

  ego_up <- run_enrichgo(up_symbols, universe_symbols, ontology = "CC")
  ego_down <- run_enrichgo(down_symbols, universe_symbols, ontology = "CC")

  if (!is.null(ego_up)) {
    up_tbl <- as.data.frame(ego_up)
    up_tbl$cell_type <- cell_type
    up_tbl$top_genes <- sapply(up_tbl$geneID, format_top_genes_one_line, res_df = res)
    write.csv(up_tbl,
              file = file.path(go_cc_ko_outdir, paste0(clean_name, "_UP_enrichGO_CC.csv")),
              row.names = FALSE)
    go_cc_ko_up_all[[cell_type]] <- up_tbl %>%
      arrange(p.adjust) %>%
      slice_head(n = 5)
  } else {
    message("No significant GO-CC enrichment for UP genes in ", cell_type, " (WT vs KO)")
  }

  if (!is.null(ego_down)) {
    down_tbl <- as.data.frame(ego_down)
    down_tbl$cell_type <- cell_type
    down_tbl$top_genes <- sapply(down_tbl$geneID, format_top_genes_one_line, res_df = res)
    write.csv(down_tbl,
              file = file.path(go_cc_ko_outdir, paste0(clean_name, "_DOWN_enrichGO_CC.csv")),
              row.names = FALSE)
    go_cc_ko_down_all[[cell_type]] <- down_tbl %>%
      arrange(p.adjust) %>%
      slice_head(n = 5)
  } else {
    message("No significant GO-CC enrichment for DOWN genes in ", cell_type, " (WT vs KO)")
  }
}

go_cc_ko_up_tbl <- if (length(go_cc_ko_up_all) > 0) bind_rows(go_cc_ko_up_all) else NULL
go_cc_ko_down_tbl <- if (length(go_cc_ko_down_all) > 0) bind_rows(go_cc_ko_down_all) else NULL

p_cc_ko_down <- make_regulation_go_dotplot(
  go_tbl = go_cc_ko_down_tbl,
  all_cell_types = all_cell_types_ko,
  title_text = "WT vs KO DEGs: GO-CC (Downregulated)",
  regulation = "Down"
)

p_cc_ko_up <- make_regulation_go_dotplot(
  go_tbl = go_cc_ko_up_tbl,
  all_cell_types = all_cell_types_ko,
  title_text = "WT vs KO DEGs: GO-CC (Upregulated)",
  regulation = "Up"
)

if (!is.null(p_cc_ko_down)) {
  n_terms_cc_ko_down <- nrow(go_cc_ko_down_tbl)
  plot_h_cc_ko_down <- max(11, 0.20 * n_terms_cc_ko_down)

  ggsave(
    filename = file.path(go_cc_ko_outdir, "WT_vs_KO_DOWN_GO_CC_dotplot.pdf"),
    plot = p_cc_ko_down,
    width = 10,
    height = plot_h_cc_ko_down,
    limitsize = FALSE
  )
}

if (!is.null(p_cc_ko_up)) {
  n_terms_cc_ko_up <- nrow(go_cc_ko_up_tbl)
  plot_h_cc_ko_up <- max(11, 0.20 * n_terms_cc_ko_up)

  ggsave(
    filename = file.path(go_cc_ko_outdir, "WT_vs_KO_UP_GO_CC_dotplot.pdf"),
    plot = p_cc_ko_up,
    width = 10,
    height = plot_h_cc_ko_up,
    limitsize = FALSE
  )
}

if (!is.null(p_cc_ko_down) && !is.null(p_cc_ko_up)) {
  combined_h_cc_ko <- max(
    max(11, 0.20 * nrow(go_cc_ko_down_tbl)),
    max(11, 0.20 * nrow(go_cc_ko_up_tbl))
  )

  grDevices::pdf(
    file = file.path(go_cc_ko_outdir, "WT_vs_KO_UP_DOWN_GO_CC_dotplot_combined.pdf"),
    width = 20,
    height = combined_h_cc_ko,
    onefile = TRUE
  )

  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))
  print(p_cc_ko_down, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p_cc_ko_up, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grDevices::dev.off()
}

message("GO-CC enrichment analysis (WT vs KO) completed. Results saved to: ", go_cc_ko_outdir)


##################################################################################################################
### KEGG analysis (WT vs KO)
##################################################################################################################

kegg_ko_outdir <- "results/scRNAseq/7_pseudobulk/kegg_wt_vs_ko"
if (!dir.exists(kegg_ko_outdir)) {
  dir.create(kegg_ko_outdir, recursive = TRUE)
}

kegg_ko_up_all <- list()
kegg_ko_down_all <- list()

for (cell_type in names(res_list_ko)) {
  res <- as.data.frame(res_list_ko[[cell_type]])
  res$gene <- rownames(res)

  up_symbols <- res %>%
    filter(!is.na(padj), padj < padj_cutoff, log2FoldChange > lfc_cutoff) %>%
    pull(gene)

  down_symbols <- res %>%
    filter(!is.na(padj), padj < padj_cutoff, log2FoldChange < -lfc_cutoff) %>%
    pull(gene)

  universe_symbols <- res %>%
    filter(!is.na(padj)) %>%
    pull(gene)

  clean_name <- gsub("[^A-Za-z0-9]", "_", cell_type)

  ek_up <- run_enrichkegg(up_symbols, universe_symbols)
  ek_down <- run_enrichkegg(down_symbols, universe_symbols)

  if (!is.null(ek_up)) {
    up_tbl <- as.data.frame(ek_up)
    up_tbl$cell_type <- cell_type
    up_tbl$top_genes <- sapply(up_tbl$geneID, format_top_genes_one_line, res_df = res)
    write.csv(up_tbl,
              file = file.path(kegg_ko_outdir, paste0(clean_name, "_UP_enrichKEGG.csv")),
              row.names = FALSE)
    kegg_ko_up_all[[cell_type]] <- up_tbl %>%
      arrange(p.adjust) %>%
      slice_head(n = 5)
  } else {
    message("No significant KEGG enrichment for UP genes in ", cell_type, " (WT vs KO)")
  }

  if (!is.null(ek_down)) {
    down_tbl <- as.data.frame(ek_down)
    down_tbl$cell_type <- cell_type
    down_tbl$top_genes <- sapply(down_tbl$geneID, format_top_genes_one_line, res_df = res)
    write.csv(down_tbl,
              file = file.path(kegg_ko_outdir, paste0(clean_name, "_DOWN_enrichKEGG.csv")),
              row.names = FALSE)
    kegg_ko_down_all[[cell_type]] <- down_tbl %>%
      arrange(p.adjust) %>%
      slice_head(n = 5)
  } else {
    message("No significant KEGG enrichment for DOWN genes in ", cell_type, " (WT vs KO)")
  }
}

kegg_ko_up_tbl <- if (length(kegg_ko_up_all) > 0) bind_rows(kegg_ko_up_all) else NULL
kegg_ko_down_tbl <- if (length(kegg_ko_down_all) > 0) bind_rows(kegg_ko_down_all) else NULL

p_kegg_ko_down <- make_regulation_go_dotplot(
  go_tbl = kegg_ko_down_tbl,
  all_cell_types = all_cell_types_ko,
  title_text = "WT vs KO DEGs: KEGG (Downregulated)",
  regulation = "Down"
)

p_kegg_ko_up <- make_regulation_go_dotplot(
  go_tbl = kegg_ko_up_tbl,
  all_cell_types = all_cell_types_ko,
  title_text = "WT vs KO DEGs: KEGG (Upregulated)",
  regulation = "Up"
)

if (!is.null(p_kegg_ko_down)) {
  n_terms_kegg_ko_down <- nrow(kegg_ko_down_tbl)
  plot_h_kegg_ko_down <- max(11, 0.20 * n_terms_kegg_ko_down)

  ggsave(
    filename = file.path(kegg_ko_outdir, "WT_vs_KO_DOWN_KEGG_dotplot.pdf"),
    plot = p_kegg_ko_down,
    width = 10,
    height = plot_h_kegg_ko_down,
    limitsize = FALSE
  )
}

if (!is.null(p_kegg_ko_up)) {
  n_terms_kegg_ko_up <- nrow(kegg_ko_up_tbl)
  plot_h_kegg_ko_up <- max(11, 0.20 * n_terms_kegg_ko_up)

  ggsave(
    filename = file.path(kegg_ko_outdir, "WT_vs_KO_UP_KEGG_dotplot.pdf"),
    plot = p_kegg_ko_up,
    width = 10,
    height = plot_h_kegg_ko_up,
    limitsize = FALSE
  )
}

if (!is.null(p_kegg_ko_down) && !is.null(p_kegg_ko_up)) {
  combined_h_kegg_ko <- max(
    max(11, 0.20 * nrow(kegg_ko_down_tbl)),
    max(11, 0.20 * nrow(kegg_ko_up_tbl))
  )

  grDevices::pdf(
    file = file.path(kegg_ko_outdir, "WT_vs_KO_UP_DOWN_KEGG_dotplot_combined.pdf"),
    width = 20,
    height = combined_h_kegg_ko,
    onefile = TRUE
  )

  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))
  print(p_kegg_ko_down, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p_kegg_ko_up, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grDevices::dev.off()
}

message("KEGG enrichment analysis (WT vs KO) completed. Results saved to: ", kegg_ko_outdir)

#########################################################################################
#### GO and KEGG analysis of DEGs between Mes-WT of WT and Mes-mut of KO ################
#########################################################################################

library(DESeq2)


# read dds r objects

dds_d20 <- readRDS("data/scRNAseq/iwaselab/dds_d20")
dds_d60 <- readRDS("data/scRNAseq/iwaselab/dds_d60")
dds_both <- readRDS("data/scRNAseq/iwaselab/dds_both")

# direct contrast: Mes-WT (WT) vs Mes-Mut (KO)
rstbl_MesWT_wt.vs.MesMut_ko <- results(
  dds_d20,
  contrast = c("genotype_cluster_d20", "ko_Mes-KO", "wt_Mes-WT"),
  alpha = 0.01
)

res_mes_wt_vs_mut <- as.data.frame(rstbl_MesWT_wt.vs.MesMut_ko)
res_mes_wt_vs_mut$gene <- rownames(res_mes_wt_vs_mut)

mes_outdir <- "results/scRNAseq/7_pseudobulk/mes_wt_vs_mes_mut_ko"
if (!dir.exists(mes_outdir)) {
  dir.create(mes_outdir, recursive = TRUE)
}

# save full DE table
write.csv(
  res_mes_wt_vs_mut,
  file = file.path(mes_outdir, "rstbl_MesWT_wt.vs.MesMut_ko.csv"),
  row.names = FALSE
)

# define significant DEGs
mes_up_symbols <- res_mes_wt_vs_mut %>%
  filter(!is.na(padj), padj < padj_cutoff, log2FoldChange > lfc_cutoff) %>%
  pull(gene)

mes_down_symbols <- res_mes_wt_vs_mut %>%
  filter(!is.na(padj), padj < padj_cutoff, log2FoldChange < -lfc_cutoff) %>%
  pull(gene)

mes_universe_symbols <- res_mes_wt_vs_mut %>%
  filter(!is.na(padj)) %>%
  pull(gene)

mes_cell_label <- "Mes-WT(WT)-vs-Mes-Mut(KO)"
mes_cell_types <- c(mes_cell_label)


##################################################################################################################
### GO-BP (Mes-WT WT vs Mes-Mut KO)
##################################################################################################################

mes_bp_up <- run_enrichgo(mes_up_symbols, mes_universe_symbols, ontology = "BP")
mes_bp_down <- run_enrichgo(mes_down_symbols, mes_universe_symbols, ontology = "BP")

mes_bp_up_tbl <- NULL
mes_bp_down_tbl <- NULL

if (!is.null(mes_bp_up)) {
  mes_bp_up_tbl <- as.data.frame(mes_bp_up)
  mes_bp_up_tbl$cell_type <- mes_cell_label
  mes_bp_up_tbl$top_genes <- sapply(mes_bp_up_tbl$geneID, format_top_genes_one_line, res_df = res_mes_wt_vs_mut)
  write.csv(
    mes_bp_up_tbl,
    file = file.path(mes_outdir, "MesWT_wt_vs_MesMut_ko_UP_enrichGO_BP.csv"),
    row.names = FALSE
  )
  mes_bp_up_tbl <- mes_bp_up_tbl %>% arrange(p.adjust) %>% slice_head(n = 10)
} else {
  message("No significant GO-BP enrichment for UP genes in Mes-WT (WT) vs Mes-Mut (KO)")
}

if (!is.null(mes_bp_down)) {
  mes_bp_down_tbl <- as.data.frame(mes_bp_down)
  mes_bp_down_tbl$cell_type <- mes_cell_label
  mes_bp_down_tbl$top_genes <- sapply(mes_bp_down_tbl$geneID, format_top_genes_one_line, res_df = res_mes_wt_vs_mut)
  write.csv(
    mes_bp_down_tbl,
    file = file.path(mes_outdir, "MesWT_wt_vs_MesMut_ko_DOWN_enrichGO_BP.csv"),
    row.names = FALSE
  )
  mes_bp_down_tbl <- mes_bp_down_tbl %>% arrange(p.adjust) %>% slice_head(n = 10)
} else {
  message("No significant GO-BP enrichment for DOWN genes in Mes-WT (WT) vs Mes-Mut (KO)")
}

p_mes_bp_down <- make_regulation_go_dotplot(
  go_tbl = mes_bp_down_tbl,
  all_cell_types = mes_cell_types,
  title_text = "Mes-WT (WT) vs Mes-Mut (KO): GO-BP (Downregulated)",
  regulation = "Down"
)

p_mes_bp_up <- make_regulation_go_dotplot(
  go_tbl = mes_bp_up_tbl,
  all_cell_types = mes_cell_types,
  title_text = "Mes-WT (WT) vs Mes-Mut (KO): GO-BP (Upregulated)",
  regulation = "Up"
)

if (!is.null(p_mes_bp_down)) {
  ggsave(
    filename = file.path(mes_outdir, "MesWT_wt_vs_MesMut_ko_DOWN_GO_BP_dotplot.pdf"),
    plot = p_mes_bp_down,
    width = 10,
    height = 11,
    limitsize = FALSE
  )
}

if (!is.null(p_mes_bp_up)) {
  ggsave(
    filename = file.path(mes_outdir, "MesWT_wt_vs_MesMut_ko_UP_GO_BP_dotplot.pdf"),
    plot = p_mes_bp_up,
    width = 10,
    height = 11,
    limitsize = FALSE
  )
}

if (!is.null(p_mes_bp_down) && !is.null(p_mes_bp_up)) {
  grDevices::pdf(
    file = file.path(mes_outdir, "MesWT_wt_vs_MesMut_ko_UP_DOWN_GO_BP_dotplot_combined.pdf"),
    width = 20,
    height = 11,
    onefile = TRUE
  )
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))
  print(p_mes_bp_down, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p_mes_bp_up, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grDevices::dev.off()
}


##################################################################################################################
### GO-CC (Mes-WT WT vs Mes-Mut KO)
##################################################################################################################

mes_cc_up <- run_enrichgo(mes_up_symbols, mes_universe_symbols, ontology = "CC")
mes_cc_down <- run_enrichgo(mes_down_symbols, mes_universe_symbols, ontology = "CC")

mes_cc_up_tbl <- NULL
mes_cc_down_tbl <- NULL

if (!is.null(mes_cc_up)) {
  mes_cc_up_tbl <- as.data.frame(mes_cc_up)
  mes_cc_up_tbl$cell_type <- mes_cell_label
  mes_cc_up_tbl$top_genes <- sapply(mes_cc_up_tbl$geneID, format_top_genes_one_line, res_df = res_mes_wt_vs_mut)
  write.csv(
    mes_cc_up_tbl,
    file = file.path(mes_outdir, "MesWT_wt_vs_MesMut_ko_UP_enrichGO_CC.csv"),
    row.names = FALSE
  )
  mes_cc_up_tbl <- mes_cc_up_tbl %>% arrange(p.adjust) %>% slice_head(n = 10)
} else {
  message("No significant GO-CC enrichment for UP genes in Mes-WT (WT) vs Mes-Mut (KO)")
}

if (!is.null(mes_cc_down)) {
  mes_cc_down_tbl <- as.data.frame(mes_cc_down)
  mes_cc_down_tbl$cell_type <- mes_cell_label
  mes_cc_down_tbl$top_genes <- sapply(mes_cc_down_tbl$geneID, format_top_genes_one_line, res_df = res_mes_wt_vs_mut)
  write.csv(
    mes_cc_down_tbl,
    file = file.path(mes_outdir, "MesWT_wt_vs_MesMut_ko_DOWN_enrichGO_CC.csv"),
    row.names = FALSE
  )
  mes_cc_down_tbl <- mes_cc_down_tbl %>% arrange(p.adjust) %>% slice_head(n = 10)
} else {
  message("No significant GO-CC enrichment for DOWN genes in Mes-WT (WT) vs Mes-Mut (KO)")
}

p_mes_cc_down <- make_regulation_go_dotplot(
  go_tbl = mes_cc_down_tbl,
  all_cell_types = mes_cell_types,
  title_text = "Mes-WT (WT) vs Mes-Mut (KO): GO-CC (Downregulated)",
  regulation = "Down"
)

p_mes_cc_up <- make_regulation_go_dotplot(
  go_tbl = mes_cc_up_tbl,
  all_cell_types = mes_cell_types,
  title_text = "Mes-WT (WT) vs Mes-Mut (KO): GO-CC (Upregulated)",
  regulation = "Up"
)

if (!is.null(p_mes_cc_down)) {
  ggsave(
    filename = file.path(mes_outdir, "MesWT_wt_vs_MesMut_ko_DOWN_GO_CC_dotplot.pdf"),
    plot = p_mes_cc_down,
    width = 10,
    height = 11,
    limitsize = FALSE
  )
}

if (!is.null(p_mes_cc_up)) {
  ggsave(
    filename = file.path(mes_outdir, "MesWT_wt_vs_MesMut_ko_UP_GO_CC_dotplot.pdf"),
    plot = p_mes_cc_up,
    width = 10,
    height = 11,
    limitsize = FALSE
  )
}

if (!is.null(p_mes_cc_down) && !is.null(p_mes_cc_up)) {
  grDevices::pdf(
    file = file.path(mes_outdir, "MesWT_wt_vs_MesMut_ko_UP_DOWN_GO_CC_dotplot_combined.pdf"),
    width = 20,
    height = 11,
    onefile = TRUE
  )
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))
  print(p_mes_cc_down, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p_mes_cc_up, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grDevices::dev.off()
}


##################################################################################################################
### KEGG (Mes-WT WT vs Mes-Mut KO)
##################################################################################################################

mes_kegg_up <- run_enrichkegg(mes_up_symbols, mes_universe_symbols)
mes_kegg_down <- run_enrichkegg(mes_down_symbols, mes_universe_symbols)

mes_kegg_up_tbl <- NULL
mes_kegg_down_tbl <- NULL

if (!is.null(mes_kegg_up)) {
  mes_kegg_up_tbl <- as.data.frame(mes_kegg_up)
  mes_kegg_up_tbl$cell_type <- mes_cell_label
  mes_kegg_up_tbl$top_genes <- sapply(mes_kegg_up_tbl$geneID, format_top_genes_one_line, res_df = res_mes_wt_vs_mut)
  write.csv(
    mes_kegg_up_tbl,
    file = file.path(mes_outdir, "MesWT_wt_vs_MesMut_ko_UP_enrichKEGG.csv"),
    row.names = FALSE
  )
  mes_kegg_up_tbl <- mes_kegg_up_tbl %>% arrange(p.adjust) %>% slice_head(n = 10)
} else {
  message("No significant KEGG enrichment for UP genes in Mes-WT (WT) vs Mes-Mut (KO)")
}

if (!is.null(mes_kegg_down)) {
  mes_kegg_down_tbl <- as.data.frame(mes_kegg_down)
  mes_kegg_down_tbl$cell_type <- mes_cell_label
  mes_kegg_down_tbl$top_genes <- sapply(mes_kegg_down_tbl$geneID, format_top_genes_one_line, res_df = res_mes_wt_vs_mut)
  write.csv(
    mes_kegg_down_tbl,
    file = file.path(mes_outdir, "MesWT_wt_vs_MesMut_ko_DOWN_enrichKEGG.csv"),
    row.names = FALSE
  )
  mes_kegg_down_tbl <- mes_kegg_down_tbl %>% arrange(p.adjust) %>% slice_head(n = 10)
} else {
  message("No significant KEGG enrichment for DOWN genes in Mes-WT (WT) vs Mes-Mut (KO)")
}

p_mes_kegg_down <- make_regulation_go_dotplot(
  go_tbl = mes_kegg_down_tbl,
  all_cell_types = mes_cell_types,
  title_text = "Mes-WT (WT) vs Mes-Mut (KO): KEGG (Downregulated)",
  regulation = "Down"
)

p_mes_kegg_up <- make_regulation_go_dotplot(
  go_tbl = mes_kegg_up_tbl,
  all_cell_types = mes_cell_types,
  title_text = "Mes-WT (WT) vs Mes-Mut (KO): KEGG (Upregulated)",
  regulation = "Up"
)

if (!is.null(p_mes_kegg_down)) {
  ggsave(
    filename = file.path(mes_outdir, "MesWT_wt_vs_MesMut_ko_DOWN_KEGG_dotplot.pdf"),
    plot = p_mes_kegg_down,
    width = 10,
    height = 11,
    limitsize = FALSE
  )
}

if (!is.null(p_mes_kegg_up)) {
  ggsave(
    filename = file.path(mes_outdir, "MesWT_wt_vs_MesMut_ko_UP_KEGG_dotplot.pdf"),
    plot = p_mes_kegg_up,
    width = 10,
    height = 11,
    limitsize = FALSE
  )
}

if (!is.null(p_mes_kegg_down) && !is.null(p_mes_kegg_up)) {
  grDevices::pdf(
    file = file.path(mes_outdir, "MesWT_wt_vs_MesMut_ko_UP_DOWN_KEGG_dotplot_combined.pdf"),
    width = 20,
    height = 11,
    onefile = TRUE
  )
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1, 2)))
  print(p_mes_kegg_down, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(p_mes_kegg_up, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
  grDevices::dev.off()
}

message("GO and KEGG analysis for Mes-WT (WT) vs Mes-Mut (KO) completed. Results saved to: ", mes_outdir)


##################################################################################################################
### Summary Excel workbooks for GO/KEGG enrichment
### 1) WT vs HET across cell types, 2) WT vs KO across cell types, 3) Mes-WT(WT) vs Mes-Mut(KO)
### Each workbook includes tabs: GO-BP, GO-CC, KEGG
##################################################################################################################

library("openxlsx")

collect_enrichment_summary <- function(dir_path, comparison_label, category_label) {
  if (!dir.exists(dir_path)) {
    return(NULL)
  }

  csv_files <- list.files(dir_path, pattern = "_enrich.*\\.csv$", full.names = TRUE)
  if (length(csv_files) == 0) {
    return(NULL)
  }

  out_list <- vector("list", length(csv_files))

  for (i in seq_along(csv_files)) {
    f <- csv_files[i]
    tbl <- tryCatch(read.csv(f, stringsAsFactors = FALSE), error = function(e) NULL)

    if (is.null(tbl) || nrow(tbl) == 0) {
      next
    }

    file_name <- basename(f)
    regulation <- ifelse(grepl("_UP_", file_name), "Up", ifelse(grepl("_DOWN_", file_name), "Down", NA_character_))

    if (!"cell_type" %in% names(tbl)) {
      tbl$cell_type <- sub("_(UP|DOWN)_.*$", "", file_name)
    }

    tbl$comparison <- comparison_label
    tbl$ontology <- category_label
    tbl$regulation <- regulation
    tbl$source_file <- file_name

    keep_cols <- c(
      "comparison", "ontology", "regulation", "cell_type",
      "ID", "Description", "GeneRatio", "BgRatio",
      "pvalue", "p.adjust", "qvalue", "Count",
      "geneID", "top_genes", "source_file"
    )

    present_cols <- intersect(keep_cols, names(tbl))
    out_list[[i]] <- tbl[, present_cols, drop = FALSE]
  }

  out_list <- out_list[!sapply(out_list, is.null)]
  if (length(out_list) == 0) {
    return(NULL)
  }

  bind_rows(out_list) %>%
    arrange(regulation, cell_type, p.adjust)
}

write_enrichment_workbook <- function(bp_tbl, cc_tbl, kegg_tbl, out_xlsx) {
  wb <- openxlsx::createWorkbook()

  openxlsx::addWorksheet(wb, "GO-BP")
  if (!is.null(bp_tbl) && nrow(bp_tbl) > 0) {
    openxlsx::writeData(wb, sheet = "GO-BP", x = bp_tbl)
  } else {
    openxlsx::writeData(wb, sheet = "GO-BP", x = data.frame(Message = "No enrichment rows found."))
  }

  openxlsx::addWorksheet(wb, "GO-CC")
  if (!is.null(cc_tbl) && nrow(cc_tbl) > 0) {
    openxlsx::writeData(wb, sheet = "GO-CC", x = cc_tbl)
  } else {
    openxlsx::writeData(wb, sheet = "GO-CC", x = data.frame(Message = "No enrichment rows found."))
  }

  openxlsx::addWorksheet(wb, "KEGG")
  if (!is.null(kegg_tbl) && nrow(kegg_tbl) > 0) {
    openxlsx::writeData(wb, sheet = "KEGG", x = kegg_tbl)
  } else {
    openxlsx::writeData(wb, sheet = "KEGG", x = data.frame(Message = "No enrichment rows found."))
  }

  openxlsx::saveWorkbook(wb, file = out_xlsx, overwrite = TRUE)
}

subset_summary_by_pattern <- function(tbl, pattern) {
  if (is.null(tbl) || nrow(tbl) == 0) {
    return(NULL)
  }
  tbl %>% filter(grepl(pattern, source_file))
}

summary_dir <- "results/scRNAseq/7_pseudobulk"
if (!dir.exists(summary_dir)) {
  dir.create(summary_dir, recursive = TRUE)
}

# 1) WT vs HET across cell types
wt_het_bp <- collect_enrichment_summary(go_outdir, "WT_vs_HET", "GO_BP")
wt_het_cc <- collect_enrichment_summary(go_cc_outdir, "WT_vs_HET", "GO_CC")
wt_het_kegg <- collect_enrichment_summary(kegg_outdir, "WT_vs_HET", "KEGG")

wt_het_xlsx <- file.path(summary_dir, "go_kegg_summary_WT_vs_HET_across_celltypes.xlsx")
write_enrichment_workbook(wt_het_bp, wt_het_cc, wt_het_kegg, wt_het_xlsx)
message("Workbook saved: ", wt_het_xlsx)

# 2) WT vs KO across cell types
wt_ko_bp <- collect_enrichment_summary(go_bp_ko_outdir, "WT_vs_KO", "GO_BP")
wt_ko_cc <- collect_enrichment_summary(go_cc_ko_outdir, "WT_vs_KO", "GO_CC")
wt_ko_kegg <- collect_enrichment_summary(kegg_ko_outdir, "WT_vs_KO", "KEGG")

wt_ko_xlsx <- file.path(summary_dir, "go_kegg_summary_WT_vs_KO_across_celltypes.xlsx")
write_enrichment_workbook(wt_ko_bp, wt_ko_cc, wt_ko_kegg, wt_ko_xlsx)
message("Workbook saved: ", wt_ko_xlsx)

# 3) Mes-WT(WT) vs Mes-Mut(KO)
mes_all <- collect_enrichment_summary(mes_outdir, "WT_MesWT_vs_KO_MesKO", "Mes_specific")
mes_bp <- subset_summary_by_pattern(mes_all, "enrichGO_BP")
mes_cc <- subset_summary_by_pattern(mes_all, "enrichGO_CC")
mes_kegg <- subset_summary_by_pattern(mes_all, "enrichKEGG")

mes_xlsx <- file.path(summary_dir, "go_kegg_summary_WT_MesWT_vs_KO_MesKO.xlsx")
write_enrichment_workbook(mes_bp, mes_cc, mes_kegg, mes_xlsx)
message("Workbook saved: ", mes_xlsx)


