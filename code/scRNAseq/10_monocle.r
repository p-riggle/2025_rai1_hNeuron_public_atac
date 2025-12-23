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

# test different number of dimensions
cds2 <- preprocess_cds(cds, num_dim = 2)
cds2 <- reduce_dimension(cds2, reduction_method = "UMAP")
cds3 <- preprocess_cds(cds, num_dim = 3)
cds3 <- reduce_dimension(cds3, reduction_method = "UMAP")
cds4 <- preprocess_cds(cds, num_dim = 4)
cds4 <- reduce_dimension(cds4, reduction_method = "UMAP")
cds5 <- preprocess_cds(cds, num_dim = 5)
cds5 <- reduce_dimension(cds5, reduction_method = "UMAP")
cds10 <- preprocess_cds(cds, num_dim = 10)
cds10 <- reduce_dimension(cds10, reduction_method = "UMAP")

# visualize the variance explained by each principal component
plot_pc_variance_explained(cds2)
plot_pc_variance_explained(cds3)
plot_pc_variance_explained(cds4)
plot_pc_variance_explained(cds5)
plot_pc_variance_explained(cds10)

## save the variance plots aligned 3 plots in one row in one pdf
library(patchwork)
p1 <- plot_pc_variance_explained(cds2) +
  ggtitle("2 PCs")
p2 <- plot_pc_variance_explained(cds3) +
  ggtitle("3 PCs")
p3 <- plot_pc_variance_explained(cds4) +
  ggtitle("4 PCs")
p4 <- plot_pc_variance_explained(cds5) +
  ggtitle("5 PCs")
p5 <- plot_pc_variance_explained(cds10) +
  ggtitle("10 PCs")

combined_plot <- p1 + p2 + p3 + p4 + p5 + plot_layout(ncol = 3) 

ggsave(
  filename = "results/scRNAseq/6_monocle3/pc_variance_explained.pdf",
  plot = combined_plot,
  width = 12,
  height = 8
)

# UMAP PLOTS *rasterization did not work well here, try later if needed*

# This grabs the first layer (the dots) and wraps it in rasterise()
# rasterize_monocle <- function(p, dpi=300) {
#   p$layers[[1]] <- ggrastr::rasterise(p$layers[[1]], dpi = dpi)
#   return(p)
# }

# Generate plots without cell type labels
umap1 <- plot_cells(cds2, label_groups_by_cluster=FALSE, color_cells_by = "cell_type", label_cell_groups = FALSE)
umap2 <- plot_cells(cds3, label_groups_by_cluster=FALSE, color_cells_by = "cell_type", label_cell_groups = FALSE)
umap3 <- plot_cells(cds4, label_groups_by_cluster=FALSE, color_cells_by = "cell_type", label_cell_groups = FALSE)
umap4 <- plot_cells(cds5, label_groups_by_cluster=FALSE, color_cells_by = "cell_type", label_cell_groups = FALSE)
umap5 <- plot_cells(cds10, label_groups_by_cluster=FALSE, color_cells_by = "cell_type", label_cell_groups = FALSE) 

# save the UMAP plots in single png file with 3 plots in one row, add the number of demention as the title for each plot
umap1 <- umap1 + ggtitle("2 PCs")
umap2 <- umap2 + ggtitle("3 PCs")
umap3 <- umap3 + ggtitle("4 PCs")
umap4 <- umap4 + ggtitle("5 PCs")
umap5 <- umap5 + ggtitle("10 PCs")

combined_umap <- umap1 + umap2 + umap3 + umap4 + umap5 + plot_layout(ncol = 3)
ggsave(
  filename = "results/scRNAseq/6_monocle3/umap_combined.png",
  plot = combined_umap,
  width = 18,
  height = 8
)

# save each plot as separate png files use loop 
umap_list <- list(umap1, umap2, umap3, umap4, umap5)
for (i in 1:length(umap_list)) {

ggsave(
  filename = paste0("results/scRNAseq/6_monocle3/umap", i, ".png"),
  plot = umap_list[[i]],
  width = 12,
  height = 8
) 
}

# umap4 and umap5 look look good. 
# Cluster your cells 
cds4 <- cluster_cells(cds4)
cds5 <- cluster_cells(cds5)
cds10 <- cluster_cells(cds10) 

# plot the 
cluster4 <- plot_cells(cds4, color_cells_by = "partition", label_cell_groups = FALSE)
cluster5 <- plot_cells(cds5, color_cells_by = "partition", label_cell_groups = FALSE)
cluster10 <- plot_cells(cds10, color_cells_by = "partition", label_cell_groups = FALSE)

# save the cluster plots in single png file with 3 plots in one row with titles indicating number of dimensions used
cluster4 <- cluster4 + ggtitle("4 PCs")
cluster5 <- cluster5 + ggtitle("5 PCs")
cluster10 <- cluster10 + ggtitle("10 PCs")

cluster4

combined_cluster <- cluster4 + cluster5 + cluster10 + plot_layout(ncol = 3)
ggsave(
  filename = "results/scRNAseq/6_monocle3/cluster_combined.png",
  plot = combined_cluster,
  width = 18,
  height = 5
)

# save each plot as separate png files using loop
cluster_list <- list(cluster4, cluster5, cluster10)
for (i in 1:length(cluster_list)) {
  ggsave(
    filename = paste0("results/scRNAseq/6_monocle3/cluster_plot", i, ".png"),
    plot = cluster_list[[i]],
    width = 12,
    height = 8
  ) 
}

# Learn the trajectory graph (5 and 10 not done yet)
cds4 <- learn_graph(cds4)
cds5 <- learn_graph(cds5)
cds10 <- learn_graph(cds10)

colnames(colData(cds10))

# save the cds4, cds5, cds10 objects
save_monocle_objects(cds = cds4, directory_path = "data/scRNAseq/iwaselab/rdata/monocle3_cds4")
save_monocle_objects(cds = cds5, directory_path = "data/scRNAseq/iwaselab/rdata/monocle3_cds5")
save_monocle_objects(cds = cds10, directory_path = "data/scRNAseq/iwaselab/rdata/monocle3_cds10")

####### section 3 #####

# load the saved cds projects
cds4 <- load_monocle_objects(directory_path = "data/scRNAseq/iwaselab/rdata/monocle3_cds4")
cds5 <- load_monocle_objects(directory_path = "data/scRNAseq/iwaselab/rdata/monocle3_cds5")
cds10 <- load_monocle_objects(directory_path = "data/scRNAseq/iwaselab/rdata/monocle3_cds10")

# define the color for each cell type
# Define this once at the top of your script
my_cell_colors <- c(
  "NPC-G2/M" = "hotpink",
  "Neuron"          = "dodgerblue",
  "Mes-WT" = "forestgreen",
  "Mes-KO"  = "yellowgreen",
  "NPC-G1" = "purple",
  "Astrocyte"       = "orangered",
  "IPC"        = "gold"
)

# visualize the trajectory graph using the defined color
# also include legends

cds10_tjt <- plot_cells(cds10,
           color_cells_by = "cell_type",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           graph_label_size=3) +
  scale_color_manual(values = my_cell_colors) +
  theme(legend.position = "right")

cds4_tjt <- plot_cells(cds4,
                              color_cells_by = "cell_type",
                              label_cell_groups=FALSE,
                              label_leaves=TRUE,
                              label_branch_points=TRUE,
                              graph_label_size=3) +
  scale_color_manual(values = my_cell_colors) +
  theme(legend.position = "right")

cds5_tjt <- plot_cells(cds5,
                              color_cells_by = "cell_type",
                              label_cell_groups=FALSE,
                              label_leaves=TRUE,
                              label_branch_points=TRUE,
                              graph_label_size=3) +
  scale_color_manual(values = my_cell_colors) +
  theme(legend.position = "right")

cds4_tjt
cds5_tjt
cds10_tjt 

dev.off()

# save the trajectory plots as png files: use loop
trajectory_list <- list(cds4_tjt, cds5_tjt, cds10_tjt)
for (i in seq_along(trajectory_list)) {
  ggsave(
    filename = paste0(
      "results/scRNAseq/6_monocle3/trajectory_cds",
      c(4, 5, 10)[i],
      ".png"
    ),
    plot = trajectory_list[[i]],
    width = 10,
    height = 8,
    units = "in"
  )
}

# save the three trajectory plots in one pdf file with 3 plots in one row
library(patchwork)
combined_trajectory <- cds4_tjt + cds5_tjt + cds10_tjt + plot_layout(ncol = 3)
ggsave(filename = 'results/scRNAseq/6_monocle3/trajectory_combined.png', 
       plot = combined_trajectory, 
       width = 30, height = 8, units = 'in')

combined_trajectory


# Turn the genotypes in all three cds (4,5,10) to an ordered factor
# Fix: Update each object individually rather than using a loop with copies

colData(cds4)$genotype <- factor(colData(cds4)$genotype, 
                                levels = c("wt", "het", "ko"))
colData(cds5)$genotype <- factor(colData(cds5)$genotype, 
                                levels = c("wt", "het", "ko"))
colData(cds10)$genotype <- factor(colData(cds10)$genotype, 
                                 levels = c("wt", "het", "ko"))
   
# Make sure ggplot2 is loaded
library(ggplot2)

# Add the + facet_wrap() layer to your plot command.
# Be sure to replace "genotype" with your actual column name.
# Create plots for all three cds objects using a loop
cds_list <- list(cds4 = cds4, cds5 = cds5, cds10 = cds10)
tjt_genotype_list <- list()

for (cds_name in names(cds_list)) {
  cds_obj <- cds_list[[cds_name]]
  
  # Create the plot
  tjt_genotype_list[[cds_name]] <- plot_cells(
    cds_obj,
    color_cells_by = "cell_type",
    label_cell_groups = FALSE,
    label_leaves = TRUE,
    label_branch_points = TRUE,
    graph_label_size = 3
  ) +
    scale_color_manual(values = my_cell_colors) +
    theme(legend.position = "right") +
    facet_wrap(~ genotype, ncol = 3)
  
  # Display the plot
  print(tjt_genotype_list[[cds_name]])
  
  # Save the plot
  ggsave(
    filename = paste0("results/scRNAseq/6_monocle3/trajectory_", cds_name, "_genotype.png"),
    plot = tjt_genotype_list[[cds_name]],
    width = 30,
    height = 8,
    units = "in"
  )
}

# similarly, make trajectory plots without trajectory for cds4, cds5, cds10 with genotype facet_wrap
umap_genotype_list <- list()
for (cds_name in names(cds_list)) {
  cds_obj <- cds_list[[cds_name]]
  
  # Create the plot
  umap_genotype_list[[cds_name]] <- plot_cells(
    cds_obj,
    color_cells_by = "cell_type",
    show_trajectory_graph = FALSE,
    label_groups_by_cluster = FALSE,
    group_label_size = 0,
    label_leaves = FALSE,
    label_branch_points = FALSE
  ) +
    scale_color_manual(values = my_cell_colors) +
    theme(legend.position = "right") +
    facet_wrap(~ genotype, ncol = 3)
  
  # Display the plot
  print(umap_genotype_list[[cds_name]])
  
  # Save the plot
  ggsave(
    filename = paste0("results/scRNAseq/6_monocle3/umap_", cds_name, "_genotype.png"),
    plot = umap_genotype_list[[cds_name]],
    width = 30,
    height = 8,
    units = "in"
  )
}

##### section 4 #####

# Load the saved cds objects
cds4 <- load_monocle_objects(directory_path = "data/scRNAseq/iwaselab/rdata/monocle3_cds4")
cds5 <- load_monocle_objects(directory_path = "data/scRNAseq/iwaselab/rdata/monocle3_cds5")
cds10 <- load_monocle_objects(directory_path = "data/scRNAseq/iwaselab/rdata/monocle3_cds10")

# Order cells in pseudotime for all three cds objects using a loop
cds_list <- list(cds4 = cds4, cds5 = cds5, cds10 = cds10)
pseudotime_plots <- list()

for (cds_name in names(cds_list)) {
  cds_obj <- cds_list[[cds_name]]
  
  # First, let's check what clusters and cell types are available
  print(paste("Processing", cds_name))
  print("Available cell types:")
  print(table(colData(cds_obj)$cell_type))
  print("Available clusters:")
  print(table(colData(cds_obj)$clusters))
  
  # Check if NPC-G2/M cells exist and which clusters they belong to
  npc_g2m_cells <- colData(cds_obj)$cell_type == "NPC-G2/M"
  if (sum(npc_g2m_cells) > 0) {
    print("NPC-G2/M cells found in clusters:")
    print(table(colData(cds_obj)$clusters[npc_g2m_cells]))

    # Select a few NPC-G2/M cells as root cells
    # (not all, just a subset for stability)
    npc_g2m_cell_names <- colnames(cds_obj)[npc_g2m_cells]
    # Take first 10 NPC-G2/M cells as root cells
    root_cells <- head(npc_g2m_cell_names, 10)
    print(paste("Selected", length(root_cells), "root cells"))

    # Order cells in pseudotime using the selected root cells
    cds_obj <- order_cells(cds_obj, root_cells = root_cells)
  } else {
    print("No NPC-G2/M cells found. Using interactive root selection.")
    # If no NPC-G2/M cells, use interactive selection
    cds_obj <- order_cells(cds_obj)
  }

  # Update the original object in the list
  cds_list[[cds_name]] <- cds_obj

  # Plot cells colored by pseudotime
  pseudotime_plots[[cds_name]] <- plot_cells(
    cds_obj,
    color_cells_by = "pseudotime",
    label_cell_groups = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE,
    graph_label_size = 0
  ) +
    scale_color_viridis_c() +
    theme(legend.position = "right") +
    ggtitle(paste("Pseudotime", cds_name))

  # Display the plot
  print(pseudotime_plots[[cds_name]])

  # Save the pseudotime plot
  ggsave(
    filename = paste0("results/scRNAseq/6_monocle3/pseudotime_",
                      cds_name, ".png"),
    plot = pseudotime_plots[[cds_name]],
    width = 10,
    height = 8,
    units = "in"
  )
}

# save the three plots in one pdf file with 3 plots in one row
library(patchwork)
combined_pseudotime <- pseudotime_plots[["cds4"]] + pseudotime_plots[["cds5"]] + pseudotime_plots[["cds10"]] + plot_layout(ncol = 3)
ggsave(filename = 'results/scRNAseq/6_monocle3/pseudotime_combined.png', 
       plot = combined_pseudotime, 
       width = 30, height = 8, units = 'in')
combined_pseudotime


# Update the original objects with ordered cells
cds4 <- cds_list[["cds4"]]
cds5 <- cds_list[["cds5"]]
cds10 <- cds_list[["cds10"]]

cds4

