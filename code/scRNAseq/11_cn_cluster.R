# Load libraries
library(Seurat)
library(BPCells)
library(tidyverse)
setwd("/nfs/turbo/umms-siwase/SHIGEKI2/rai1_hesc")

# read cts objects, i.e. aggregated expression data across cell types
cts <- readRDS("data/scRNAseq/iwaselab/cts")
head(cts)
str(cts)

# extract only WT data - fix this
cts_wt <- cts[, grepl("WT", colnames(cts))]

# load the cluster from bulk-RNAseq analysis
hcl_het <- read.table("results/cn_rnaseq/cluster/sep/mh_het_2026.04.13.txt", sep = "\t", header = TRUE, row.names = 1)
head(hcl_het)

hcl_ko <- read.table("results/cn_rnaseq/cluster/sep/mh_ko_2026.04.13.txt", sep = "\t", header = TRUE, row.names = 1)
head(hcl_ko)

# extract cluster 1 genes from hcl_het
cluster1_genes <- unique(hcl_het$gene[hcl_het$cluster == 1])
head(cluster1_genes)

# 1. Find genes that exist in BOTH the cluster list and the cts matrix
genes_to_keep <- intersect(cluster1_genes, rownames(cts_wt))

# 2. Check how many genes were lost (optional but helpful)
message("Found ", length(genes_to_keep), " out of ", length(cluster1_genes), " genes.")

# 3. Subset safely
cts_cluster1_sparse <- cts_wt[genes_to_keep, ]

library(dplyr)
library(tidyr)
library(tibble)

cts_cluster1_long <- as.data.frame(as.matrix(cts_cluster1_sparse)) %>%
  rownames_to_column(var = "Gene") %>%
  pivot_longer(-Gene, names_to = "Sample_Full", values_to = "Expression") %>%
  # Extract just the part after the underscore
  mutate(CellType = sub(".*_", "", Sample_Full))
head(cts_cluster1_long)
str(cts_cluster1_long)


# make a violin plot of the expression of the cluster 1 genes across cell types
## make y axis in log scale
ggplot(cts_cluster1_long, aes(x = CellType, y = Expression)) +
  geom_violin() +
  theme_bw() +
  labs(title = "Expression of Cluster 1 Genes Across Cell Types", x = "Cell Type", y = "Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# make a box plot of the expression of the cluster 1 genes across cell types
ggplot(cts_cluster1_long, aes(x = CellType, y = Expression)) +
  geom_boxplot() +
  theme_bw() +
  labs(title = "Expression of Cluster 1 Genes Across Cell Types", x = "Cell Type", y = "Expression") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
