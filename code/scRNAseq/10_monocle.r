# set working directory
setwd("/nfs/turbo/umms-siwase/SHIGEKI2/rai1_hesc")

# load seurat
library(Seurat)
library(BPCells)
library(tidyverse)
library(monocle3)

##### section 1 #####

# load seurat object
geo_so = readRDS('data/scRNAseq/iwaselab/rdata/geo_so_sct_unintegrated_with_annos.rds')

# Load the Matrix package if it's not already
library(Matrix) 

# 1. Get the counts matrix (same as before)
counts <- GetAssayData(geo_so, assay = "RNA", slot = "counts")

# 2. (THE FIX) Explicitly convert to a standard sparse matrix
counts_sparse <- as(counts, "sparseMatrix")

# 3. Get the cell metadata
cell_meta <- geo_so@meta.data

# 4. Get the gene metadata
gene_meta <- data.frame(gene_short_name = rownames(counts_sparse))
rownames(gene_meta) <- rownames(counts_sparse)

# 5. Create the CDS object using the converted sparse matrix -> worked!
cds <- new_cell_data_set(
  expression_data = counts_sparse,
  cell_metadata = cell_meta,
  gene_metadata = gene_meta
)

# --- Add this code AFTER creating the cds object ---

# 1. Transfer the Seurat UMAP coordinates
# Make sure the reduction name matches what's in your Seurat object (e.g., "umap")
seurat_umap <- geo_so@reductions$umap@cell.embeddings
reducedDims(cds)[["UMAP"]] <- seurat_umap

# 2. Transfer the Seurat clusters (or any other metadata)
# This is often done by default since we passed the full cell_meta,
# but it's good to confirm. Your Seurat clusters are now in the cds colData.
# For example, you can now color your Monocle plot by your Seurat clusters.

plot_cells(cds, color_cells_by = "cell_type")


# 3. (Important) Manually add the cluster/partition data
# Monocle 3's trajectory functions rely on "clusters" and "partitions"
# We can set them from our Seurat clusters to start.
cds@clusters@listData[["UMAP"]]$clusters <- geo_so$cell_type

# In Monocle 3, a "partition" is a high-level group of cells. When you run the learn_graph() function to draw a trajectory, Monocle will never draw a line connecting cells that are in different partitions.
# Normally: You would run Monocle's cluster_cells() function, which automatically analyzes your data and decides if you have separate, disconnected groups of cells (e.g., T-cells and B-cells, which aren't on the same differentiation path). It would then assign them to different partitions (Partition 1, Partition 2, etc.).
# In Our Case: We are manually loading data from Seurat and want to force Monocle to analyze all the cells together.
# By setting all partitions to 1, we are essentially telling Monocle: "Treat all of these cells as one single, connected dataset. Assume they are all part of the same biological story and try to find a single trajectory graph that fits them all."

cds@clusters@listData[["UMAP"]]$partitions <- as.factor(rep(1, ncol(cds)))

# 1. Save the object
save_monocle_objects(cds = cds, directory_path = "data/scRNAseq/iwaselab/rdata/monocle3")

##### section 2 #####

# Point to the directory you created earlier
cds <- load_monocle_objects(directory_path = "data/scRNAseq/iwaselab/rdata/monocle3")

# Constructing single-cell trajectories
## Pre-process the data
## Reduce dimensionality and visualize the results
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds)

plot_cells(cds, label_groups_by_cluster=FALSE,  color_cells_by = "cell_type")

# Cluster your cells
cds <- cluster_cells(cds)
plot_cells(cds, color_cells_by = "partition")

# Learn the trajectory graph
cds <- learn_graph(cds)

colnames(colData(cds))

plot_cells(cds,
           color_cells_by = "cell_type",
           label_groups_by_cluster = TRUE,   # <-- MUST be TRUE to show labels
           label_leaves = TRUE,
           label_branch_points = TRUE)

plot_cells(cds,
           color_cells_by = "cell_type",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=3)

# This line modifies your cds object directly
colData(cds)$genotype <- factor(colData(cds)$genotype, 
                                levels = c("wt", "het", "ko"))
# Make sure ggplot2 is loaded
library(ggplot2)

# Add the + facet_wrap() layer to your plot command.
# Be sure to replace "genotype" with your actual column name.

monocle_umap <- plot_cells(cds,
                           color_cells_by = "cell_type",
                           show_trajectory_graph = FALSE,  # <-- Add this line
                           label_groups_by_cluster = TRUE,
                           group_label_size = 3,
                           label_leaves = FALSE,
                           label_branch_points = FALSE) +
  facet_wrap(~ genotype, ncol = 3)

monocle_umap_with_t <- plot_cells(cds,
           color_cells_by = "cell_type",
           label_groups_by_cluster = TRUE,
           group_label_size = 3,
           label_leaves = FALSE,
           label_branch_points = FALSE) +
  facet_wrap(~ genotype, ncol = 3)

monocle_umap
monocle_umap_with_t

ggsave(filename = 'results/scRNAseq/5_monocle3/monocle_umap.png', 
       plot = monocle_umap, 
       width = 20, height = 8, units = 'in')

ggsave(filename = 'results/scRNAseq/5_monocle3/monocle_umap_with_t.png', 
       plot = monocle_umap_with_t, 
       width = 20, height = 8, units = 'in')

# Save the object
save_monocle_objects(cds = cds, directory_path = "data/scRNAseq/iwaselab/rdata/monocle3")

##### section 3 #####

# Load the directory you created earlier
cds <- load_monocle_objects(directory_path = "data/scRNAseq/iwaselab/rdata/monocle3")

cds <- order_cells(cds)
