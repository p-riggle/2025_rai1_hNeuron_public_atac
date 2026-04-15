# Load libraries
library(Seurat)
library(BPCells)
library(tidyverse)
setwd("/nfs/turbo/umms-siwase/SHIGEKI2/rai1_hesc")

# read cts objects, i.e. aggregated expression data across cell types
cts <- readRDS("data/scRNAseq/iwaselab/cts")
head(cts)
str(cts)

library(Matrix)
library(tidyverse)

# 1. Convert sparse matrix to a long data frame of non-zero entries
# summary() on a sparse matrix gives i (row index), j (col index), and x (value)
cts_long <- summary(cts) %>%
  as.data.frame() %>%
  mutate(
    gene = rownames(cts)[i],
    sample_id = colnames(cts)[j]
  ) %>%
  select(gene, sample_id, count = x)

# 2. Parse the sample_id into the specific metadata columns
cts_final <- cts_long %>%
  separate(
    sample_id, 
    into = c("day_raw", "genotype", "replicate_celltype"), 
    sep = "-", 
    remove = FALSE
  ) %>%
  separate(
    replicate_celltype, 
    into = c("replicate", "cell_type"), 
    sep = "_", 
    extra = "merge"
  ) %>%
  mutate(
    day = as.numeric(gsub("day", "", day_raw))
  ) %>%
  select(gene, day, genotype, replicate, cell_type, count)

# View the result
head(cts_final)
unique(cts_final$cell_type)

# load the cluster from bulk-RNAseq analysis
hcl_het <- read.table("results/cn_rnaseq/cluster/sep/mh_het_2026.04.13.txt", sep = "\t", header = TRUE, row.names = 1)
head(hcl_het)

hcl_ko <- read.table("results/cn_rnaseq/cluster/sep/mh_ko_2026.04.13.txt", sep = "\t", header = TRUE, row.names = 1)
head(hcl_ko)

# extract cluster 1 genes from hcl_het
cluster1_genes <- unique(hcl_het$gene[hcl_het$cluster == 1])
head(cluster1_genes)

# 1. & 2. (Your previous steps are fine)
genes_to_keep <- intersect(cluster1_genes, cts_final$gene)

# 3. FIX: Use filter() or %in% to subset the rows based on the gene column
cts_cluster1 <- cts_final %>% 
  filter(gene %in% genes_to_keep)

# 4. Plotting (Adding a check to ensure data exists)
if(nrow(cts_cluster1) > 0) {
  cts_cluster1_summary <- cts_cluster1 %>%
    group_by(cell_type) %>%
    summarise(mean_count = mean(count, na.rm = TRUE), .groups = "drop")
  
  ggplot(cts_cluster1_summary, aes(x = cell_type, y = mean_count)) +
    geom_col(fill = "steelblue") +  # geom_col is shorthand for geom_bar(stat="identity")
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Mean Expression of Cluster 1 Genes Across Cell Types", 
         x = "Cell Type", 
         y = "Mean Count")
} else {
  message("Warning: cts_cluster1 is empty. Check if gene names match exactly.")
}

################################
# 1. Find genes that exist in BOTH the cluster list and the cts matrix
genes_to_keep <- intersect(cluster1_genes, cts_final$gene)

# 2. Check how many genes were lost (optional but helpful)
message("Found ", length(genes_to_keep), " out of ", length(cluster1_genes), " genes.")

# 3. Subset safely
cts_cluster1 <- cts_final[genes_to_keep, ]

library(dplyr)
library(tidyr)
library(tibble)

# 4. plot the expression of cluster 1 genes across cell types

cts_cluster1_summary <- cts_cluster1 %>%
  group_by(cell_type) %>%
  summarise(mean_count = mean(count), .groups = "drop")

ggplot(cts_cluster1_summary, aes(x = cell_type, y = mean_count)) +
  geom_bar(stat = "identity") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Mean Expression of Cluster 1 Genes Across Cell Types", x = "Cell Type", y = "Mean Count")


