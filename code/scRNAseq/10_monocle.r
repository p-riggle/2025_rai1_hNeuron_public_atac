# set working directory
setwd("/nfs/turbo/umms-siwase/SHIGEKI2/rai1_hesc")

# load seurat
library(Seurat)
library(BPCells)
library(tidyverse)
library(monocle3)

preview_plot <- function(plot_obj) {
  if (interactive() && identical(Sys.getenv("MONOCLE_PREVIEW_PLOTS"), "1")) {
    print(plot_obj)
  }

  invisible(plot_obj)
}

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

preview_plot(plot_cells(cds, color_cells_by = "cell_type"))


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
preview_plot(plot_pc_variance_explained(cds2))
preview_plot(plot_pc_variance_explained(cds3))
preview_plot(plot_pc_variance_explained(cds4))
preview_plot(plot_pc_variance_explained(cds5))
preview_plot(plot_pc_variance_explained(cds10))

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

# Define this once at the top of your script
my_cell_colors <- c(
  "NPC-G2/M" = "hotpink",
  "Neuron"          = "dodgerblue",
  "Mes-WT" = "forestgreen",
  "Mes-KO"  = "yellowgreen",
  "NPC-G1" = "purple",
  "gRG"       = "orangered",
  "IPC"        = "gold"
)

# Generate UMAP plots with consistent colors using loop
cds_objects <- list(cds2 = cds2, cds3 = cds3, cds4 = cds4, cds5 = cds5, cds10 = cds10)
pc_counts <- c(2, 3, 4, 5, 10)
umap_list <- list()

for (i in seq_along(cds_objects)) {
  cds_name <- names(cds_objects)[i]
  pc_count <- pc_counts[i]
  
  umap_list[[i]] <- plot_cells(
    cds_objects[[i]],
    label_groups_by_cluster = FALSE,
    color_cells_by = "cell_type",
    label_cell_groups = FALSE
  ) +
    scale_color_manual(values = my_cell_colors) +
    ggtitle(paste(pc_count, "PCs"))
}

combined_umap <- umap_list[[1]] + umap_list[[2]] + umap_list[[3]] + umap_list[[4]] + umap_list[[5]] + plot_layout(ncol = 3)
ggsave(
  filename = "results/scRNAseq/6_monocle3/umap_combined.png",
  plot = combined_umap,
  width = 18,
  height = 8
)

# save each plot as separate png files using loop
for (i in seq_along(umap_list)) {
  ggsave(
    filename = paste0("results/scRNAseq/6_monocle3/umap", i, ".png"),
    plot = umap_list[[i]],
    width = 12,
    height = 8
  ) 
}

# Generate UMAP plots with consistent colors using loop
cds_objects <- list(cds2 = cds2, cds3 = cds3, cds4 = cds4, cds5 = cds5, cds10 = cds10)
pc_counts <- c(2, 3, 4, 5, 10)
umap_list <- list()

for (i in seq_along(cds_objects)) {
  cds_name <- names(cds_objects)[i]
  pc_count <- pc_counts[i]
  
  umap_list[[i]] <- plot_cells(
    cds_objects[[i]],
    label_groups_by_cluster = FALSE,
    color_cells_by = "cell_type",
    label_cell_groups = FALSE
  ) +
    scale_color_manual(values = my_cell_colors) +
    ggtitle(paste(pc_count, "PCs"))
}

combined_umap <- umap_list[[1]] + umap_list[[2]] + umap_list[[3]] + umap_list[[4]] + umap_list[[5]] + plot_layout(ncol = 3)
ggsave(
  filename = "results/scRNAseq/6_monocle3/umap_combined.png",
  plot = combined_umap,
  width = 18,
  height = 8
)

# save each plot as separate png files using loop
for (i in seq_along(umap_list)) {
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

preview_plot(cluster4)

combined_cluster <- cluster4 + cluster5 + cluster10 + plot_layout(ncol = 3)
ggsave(
  filename = "results/scRNAseq/6_monocle3/cluster_combined.png",
  plot = combined_cluster,
  width = 18,
  height = 5
)

# save each plot as separate png files using loop
cluster_list <- list(cluster4, cluster5, cluster10)
for (i in seq_along(cluster_list)) {
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
  "gRG"       = "orangered",
  "IPC"        = "gold"
)

# make the umap plots with the new color scheme

# Generate UMAP plots with consistent colors using loop
cds_objects <- list(cds4 = cds4, cds5 = cds5, cds10 = cds10)
pc_counts <- c(4, 5, 10)
umap_list <- list()

for (i in seq_along(cds_objects)) {
  cds_name <- names(cds_objects)[i]
  pc_count <- pc_counts[i]
  
  umap_list[[i]] <- plot_cells(
    cds_objects[[i]],
    label_groups_by_cluster = FALSE,
    color_cells_by = "cell_type",
    label_cell_groups = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE,
    graph_label_size = 0
  ) +
    scale_color_manual(values = my_cell_colors) +
    ggtitle(paste(pc_count, "PCs"))
}

combined_umap <- umap_list[[1]] + umap_list[[2]] + umap_list[[3]] + plot_layout(ncol = 3)
ggsave(
  filename = "results/scRNAseq/6_monocle3/umap_combined.png",
  plot = combined_umap,
  width = 18,
  height = 5
)

# save each plot as separate png files using loop
for (i in seq_along(umap_list)) {
  ggsave(
    filename = paste0("results/scRNAseq/6_monocle3/umap_pc", pc_counts[i], ".png"),
    plot = umap_list[[i]],
    width = 12,
    height = 8
  )
}

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

preview_plot(cds4_tjt)
preview_plot(cds5_tjt)
preview_plot(cds10_tjt)

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

preview_plot(combined_trajectory)


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

# Convert genotype to ordered factor for proper facet ordering
colData(cds4)$genotype <- factor(colData(cds4)$genotype, 
                                 levels = c("wt", "het", "ko"))
colData(cds5)$genotype <- factor(colData(cds5)$genotype, 
                                 levels = c("wt", "het", "ko"))
colData(cds10)$genotype <- factor(colData(cds10)$genotype, 
                                  levels = c("wt", "het", "ko"))

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
  # Filter for WT NPC-G2/M cells specifically
  npc_g2m_cells <- colData(cds_obj)$cell_type == "NPC-G2/M" & colData(cds_obj)$genotype == "wt"
  if (sum(npc_g2m_cells) > 0) {
    print("WT NPC-G2/M cells found in clusters:")
    print(table(colData(cds_obj)$clusters[npc_g2m_cells]))

    # Select a few WT NPC-G2/M cells as root cells from specific UMAP region
    # (not all, just a subset for stability from the highlighted region)
    npc_g2m_cell_names <- colnames(cds_obj)[npc_g2m_cells]
    
    # Get UMAP coordinates for filtering to specific region
    umap_coords <- reducedDims(cds_obj)[["UMAP"]]
    npc_g2m_coords <- umap_coords[npc_g2m_cells, ]
    
    # Filter for cells in the upper right region (highlighted in yellow)
    # Adjust these coordinates based on your specific region of interest
    region_filter <- npc_g2m_coords[, 1] > 5 & npc_g2m_coords[, 2] > 5
    regional_npc_g2m_cells <- npc_g2m_cell_names[region_filter]
    
    print(paste("Found", length(regional_npc_g2m_cells), "WT NPC-G2/M cells in target region"))

    # Better approach: randomly sample 10 cells from the regional subset
    set.seed(42)  # For reproducible results
    if (length(regional_npc_g2m_cells) >= 10) {
      root_cells <- sample(regional_npc_g2m_cells, 10)
    } else if (length(regional_npc_g2m_cells) > 0) {
      root_cells <- regional_npc_g2m_cells  # Use all regional cells if less than 10
    } else {
      # Fall back to any WT NPC-G2/M cells if none in region
      root_cells <- sample(npc_g2m_cell_names, min(10, length(npc_g2m_cell_names)))
      print("No cells found in target region, using random WT NPC-G2/M cells")
    }
    print(paste("Selected", length(root_cells), "WT root cells from target region"))

    # Order cells in pseudotime using the selected root cells
    cds_obj <- order_cells(cds_obj, root_cells = root_cells)
  } else {
    print("No WT NPC-G2/M cells found. Using interactive root selection.")
    # If no WT NPC-G2/M cells, use interactive selection
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
    graph_label_size = 0,
    show_trajectory_graph = FALSE
  ) +
    scale_color_viridis_c(name = "Pseudotime") +
    theme(legend.position = "right") +
    ggtitle(paste("Pseudotime", cds_name))

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
combined_pseudotime <- pseudotime_plots[["cds4"]] +
  pseudotime_plots[["cds5"]] +
  pseudotime_plots[["cds10"]] +
  plot_layout(ncol = 3)
ggsave(
  filename = "results/scRNAseq/6_monocle3/pseudotime_combined.png",
  plot = combined_pseudotime,
  width = 30,
  height = 8,
  units = "in"
)
preview_plot(combined_pseudotime)

# Create pseudotime plots with genotype faceting for all three cds objects using loop
pseudotime_genotype_plots <- list()

for (cds_name in names(cds_list)) {
  # Create pseudotime plot with genotype faceting
  pseudotime_genotype_plots[[cds_name]] <- plot_cells(
    cds_list[[cds_name]],
    color_cells_by = "pseudotime",
    label_cell_groups = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE,
    graph_label_size = 0,
    show_trajectory_graph = FALSE
  ) +
    scale_color_viridis_c(name = "Pseudotime") +
    theme(legend.position = "right") +
    ggtitle(paste("Pseudotime", cds_name)) +
    facet_wrap(~ genotype, ncol = 3)
  
  # Save the plot
  ggsave(
    filename = paste0("results/scRNAseq/6_monocle3/pseudotime_", cds_name, "_genotype.png"),
    plot = pseudotime_genotype_plots[[cds_name]],
    width = 30,
    height = 8,
    units = "in"
  )
}

# Update the original objects with ordered cells
cds4 <- cds_list[["cds4"]]
cds5 <- cds_list[["cds5"]]
cds10 <- cds_list[["cds10"]]

# save the updated cds objects with pseudotime
save_monocle_objects(cds = cds4, directory_path = "data/scRNAseq/iwaselab/rdata/monocle3_cds4_pseudotime")
save_monocle_objects(cds = cds5, directory_path = "data/scRNAseq/iwaselab/rdata/monocle3_cds5_pseudotime")
save_monocle_objects(cds = cds10, directory_path = "data/scRNAseq/iwaselab/rdata/monocle3_cds10_pseudotime")

##### section 5 #####

# Load the saved cds objects with pseudotime
cds4 <- load_monocle_objects(directory_path = "data/scRNAseq/iwaselab/rdata/monocle3_cds4_pseudotime")
cds5 <- load_monocle_objects(directory_path = "data/scRNAseq/iwaselab/rdata/monocle3_cds5_pseudotime")
cds10 <- load_monocle_objects(directory_path = "data/scRNAseq/iwaselab/rdata/monocle3_cds10_pseudotime")

# make stack plots for binned pseudotime for cds5. Three stacked bar plots represent three genotypes (wt, het, ko). 
library(ggplot2)
library(dplyr)
# Function to create stacked bar plot for a given cds object
create_stacked_bar_plot <- function(cds_obj, cds_name) {
  # Extract pseudotime and genotype information
  cell_data <- data.frame(
    pseudotime = pseudotime(cds_obj),
    genotype = colData(cds_obj)$genotype
  )
  
  # Remove cells with NA or infinite pseudotime
  cell_data <- cell_data[!is.na(cell_data$pseudotime) & is.finite(cell_data$pseudotime), ]
  
  # Check if we have enough data and variation for binning
  if (nrow(cell_data) == 0) {
    stop(paste("No valid pseudotime data found for", cds_name))
  }
  
  # Check if pseudotime has variation
  if (diff(range(cell_data$pseudotime)) == 0) {
    # All pseudotime values are the same, create a single bin
    cell_data$pseudotime_bin <- factor("Single Value")
    print(paste("Warning:", cds_name, "has no pseudotime variation, using single bin"))
  } else {
    # Create pseudotime bins with error handling
    n_unique <- length(unique(cell_data$pseudotime))
    n_bins <- min(10, max(2, n_unique))  # Use fewer bins if not enough unique values
    
    # Create bins and get the actual ranges
    cell_data$pseudotime_bin <- cut(cell_data$pseudotime, 
                                   breaks = n_bins,
                                   include.lowest = TRUE)
    
    # Convert factor levels to show numerical ranges
    levels(cell_data$pseudotime_bin) <- gsub("\\(|\\]", "", levels(cell_data$pseudotime_bin))
    levels(cell_data$pseudotime_bin) <- gsub(",", "-", levels(cell_data$pseudotime_bin))
    
    # Round the numbers for cleaner display
    new_labels <- levels(cell_data$pseudotime_bin)
    for (i in seq_along(new_labels)) {
      parts <- strsplit(new_labels[i], "-")[[1]]
      if (length(parts) == 2) {
        new_labels[i] <- paste(round(as.numeric(parts[1]), 1), 
                              round(as.numeric(parts[2]), 1), 
                              sep = "-")
      }
    }
    levels(cell_data$pseudotime_bin) <- new_labels
  }
                                 
  # Count cells in each pseudotime bin for each genotype
  count_data <- cell_data %>%
    group_by(genotype, pseudotime_bin) %>%
    summarise(cell_count = n(), .groups = "drop")
  
  # Create stacked bar plot
  plot <- ggplot(count_data, aes(x = pseudotime_bin, y = cell_count, fill = genotype)) +
    geom_bar(stat = "identity") +
    labs(title = paste("Stacked Bar Plot of Pseudotime Bins -",
                      cds_name),
          x = "Pseudotime Bins",
          y = "Cell Count") +
    theme_minimal() +
    scale_fill_manual(values = c("wt" = "#7b7b7b", "het" = "#f1d174", "ko" = "#bd6565")) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  return(plot)
}

# Create and save stacked bar plots for each cds object
cds_list <- list(cds4 = cds4, cds5 = cds5, cds10 = cds10)
for (cds_name in names(cds_list)) {
  plot <- create_stacked_bar_plot(cds_list[[cds_name]], cds_name)
  ggsave(
    filename = paste0("results/scRNAseq/6_monocle3/stacked_bar_pseudotime_", cds_name, ".pdf"),
    plot = plot,
    width = 10,
    height = 6,
    units = "in"
  )
}

# save the three stacked bar plots in one pdf file with 3 plots in one row
library(patchwork)
stacked_bar_plots <- list()
for (cds_name in names(cds_list)) {
  plot <- create_stacked_bar_plot(cds_list[[cds_name]], cds_name)
  stacked_bar_plots[[cds_name]] <- plot
}
combined_stacked_bar <- stacked_bar_plots[["cds4"]] +
  stacked_bar_plots[["cds5"]] +
  stacked_bar_plots[["cds10"]] +
  plot_layout(ncol = 3)
ggsave(
  filename = "results/scRNAseq/6_monocle3/stacked_bar_pseudotime_combined.pdf",
  plot = combined_stacked_bar,
  width = 30,
  height = 8,
  units = "in"
)
preview_plot(combined_stacked_bar)

# in the above plots, x-axis is the bined pseudotime, and each stack is genotype. Now I am making stack plots where x-axis is genotype, and each stack is binned pseudotime. 

# Function to create stacked bar plot for a given cds object with genotype on x-axis
create_stacked_bar_plot_genotype_x <- function(cds_obj, cds_name) {
  # Extract pseudotime and genotype information
  cell_data <- data.frame(
    pseudotime = pseudotime(cds_obj),
    genotype = colData(cds_obj)$genotype
  )
  # Remove cells with NA or infinite pseudotime
  cell_data <- cell_data[!is.na(cell_data$pseudotime) & is.finite(cell_data$pseudotime), ]
  # Check if we have enough data and variation for binning
  if (nrow(cell_data) == 0) {
    stop(paste("No valid pseudotime data found for", cds_name))
  }
  # Check if pseudotime has variation
  if (diff(range(cell_data$pseudotime)) == 0) {
    # All pseudotime values are the same, create a single bin
    cell_data$pseudotime_bin <- factor("Single Value")
    print(paste("Warning:", cds_name, "has no pseudotime variation, using
  single bin"))
  } else {
    # Create pseudotime bins with error handling
    n_unique <- length(unique(cell_data$pseudotime))
    n_bins <- min(10, max(2, n_unique))  # Use fewer bins if not enough unique values
    # Create bins and get the actual ranges
    cell_data$pseudotime_bin <- cut(cell_data$pseudotime, 
                                   breaks = n_bins,
                                   include.lowest = TRUE)
    # Convert factor levels to show numerical ranges
    levels(cell_data$pseudotime_bin) <- gsub("\\(|\\]", "", levels(cell_data$pseudotime_bin))
    levels(cell_data$pseudotime_bin) <- gsub(",", "-", levels(cell_data$pseudotime_bin))
    # Round the numbers for cleaner display
    new_labels <- levels(cell_data$pseudotime_bin)
    for (i in seq_along(new_labels)) {
      parts <- strsplit(new_labels[i], "-")[[1]]
      if (length(parts) == 2) {
        new_labels[i] <- paste(round(as.numeric(parts[1]), 1), 
                              round(as.numeric(parts[2]), 1), 
                              sep = "-")
      }
    }
    levels(cell_data$pseudotime_bin) <- new_labels
  }
  # Count cells in each pseudotime bin for each genotype
  count_data <- cell_data %>%
    group_by(genotype, pseudotime_bin) %>%
    summarise(cell_count = n(), .groups = "drop") %>%
    group_by(genotype) %>%
    mutate(percentage = cell_count / sum(cell_count) * 100) %>%
    ungroup()
  # Create stacked bar plot with genotype on x-axis
  plot <- ggplot(count_data, aes(x = genotype, y = percentage, fill = pseudotime_bin)) +
    geom_bar(stat = "identity") +
    labs(title = paste("Stacked Bar Plot of Genotype -",
                      cds_name),
          x = "Genotype",
          y = "Percentage of Cells") +
    theme_minimal() +
    scale_fill_viridis_d(name = "Pseudotime Bins") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  return(plot)
}
# Create and save stacked bar plots for each cds object with genotype on x-axis
cds_list <- list(cds4 = cds4, cds5 = cds5, cds10 = cds10)
for (cds_name in names(cds_list)) {
  plot <- create_stacked_bar_plot_genotype_x(cds_list[[cds_name]], cds_name)
  ggsave(
    filename = paste0("results/scRNAseq/6_monocle3/st_pt_geno_", cds_name, ".pdf"),
    plot = plot,
    width = 10,
    height = 6,
    units = "in"
  )
}
# save the three stacked bar plots in one pdf file with 3 plots in one row
library(patchwork)
stacked_bar_genotype_x_plots <- list()
for (cds_name in names(cds_list)) {
  plot <- create_stacked_bar_plot_genotype_x(cds_list[[cds_name]], cds_name)
  stacked_bar_genotype_x_plots[[cds_name]] <- plot
}
combined_stacked_bar_genotype_x <- stacked_bar_genotype_x_plots[["cds4"]] +
  stacked_bar_genotype_x_plots[["cds5"]] +
  stacked_bar_genotype_x_plots[["cds10"]] +
  plot_layout(ncol = 3)
ggsave(
  filename = "results/scRNAseq/6_monocle3/st_pt_geno_combined.pdf",
  plot = combined_stacked_bar_genotype_x,
  width = 30,
  height = 8,
  units = "in"
)
preview_plot(combined_stacked_bar_genotype_x)

###### section 6 #########

# set working directory
setwd("/nfs/turbo/umms-siwase/SHIGEKI2/rai1_hesc")

# load seurat
library(Seurat)
library(BPCells)
library(tidyverse)
library(monocle3)

# Load the saved cds objects with pseudotime
cds5 <- load_monocle_objects(directory_path = "data/scRNAseq/iwaselab/rdata/monocle3_cds5_pseudotime")

# next subset the cds5 into the 6 cell types. Then assess whether the distribution of pseudotime in each cell types differ by genotypes. 

# Convert the cds5 to a data frame for easier manipulation
cds5_df <- as.data.frame(colData(cds5))
cds5_df$pseudotime <- pseudotime(cds5)
# Check the structure of the data frame
str(cds5_df)

# save the cds5_df as rdata file
save(cds5_df, file = "data/scRNAseq/iwaselab/rdata/cds5_pseudotime_dataframe.RData")

# read the cds5_df from rdata file
load("data/scRNAseq/iwaselab/rdata/cds5_pseudotime_dataframe.RData")

# Create violin plots of pseudotime by genotype for each cell type
# Using base R statistical tests instead of ggpubr (which requires newer ggplot2)

# Create the violin subdirectory if it doesn't exist
if (!dir.exists("results/scRNAseq/6_monocle3/violin")) {
  dir.create("results/scRNAseq/6_monocle3/violin", recursive = TRUE)
}

cell_types <- unique(cds5_df$cell_type)
for (cell_type in cell_types) {
  # Subset data for the current cell type
  subset_df <- cds5_df %>% filter(cell_type == !!cell_type)
  
  # Create a filesystem-safe filename
  clean_cell_type <- gsub("[/\\-]", "_", cell_type)
  
  # Calculate statistical comparisons using base R
  comparisons <- list(c("wt", "het"), c("wt", "ko"), c("het", "ko"))
  p_values <- c()
  
  for (comp in comparisons) {
    group1_data <- subset_df$pseudotime[subset_df$genotype == comp[1]]
    group2_data <- subset_df$pseudotime[subset_df$genotype == comp[2]]
    
    # Remove NA values
    group1_data <- group1_data[!is.na(group1_data)]
    group2_data <- group2_data[!is.na(group2_data)]
    
    # Perform Wilcoxon test if both groups have data
    if (length(group1_data) > 0 & length(group2_data) > 0) {
      test_result <- wilcox.test(group1_data, group2_data)
      p_val <- test_result$p.value
      p_values <- c(p_values, ifelse(p_val < 0.001, "***", 
                                    ifelse(p_val < 0.01, "**", 
                                          ifelse(p_val < 0.05, "*", "ns"))))
    } else {
      p_values <- c(p_values, "ns")
    }
  }
  
  # Create violin plot
  p <- ggplot(subset_df, aes(x = genotype, y = pseudotime)) +
    geom_violin(trim = FALSE, fill = "lightblue") +
    geom_boxplot(width = 0.1, position = position_dodge(0.9), outlier.shape = NA) +
    labs(title = paste("Pseudotime Distribution in", cell_type),
         x = "Genotype",
         y = "Pseudotime") +
    theme_classic()
  
  # Add manual statistical annotations
  y_max <- max(subset_df$pseudotime, na.rm = TRUE)
  y_range <- max(subset_df$pseudotime, na.rm = TRUE) - min(subset_df$pseudotime, na.rm = TRUE)
  
  # Add comparison annotations manually
  p <- p + annotate("text", x = 1.5, y = y_max + y_range * 0.1, 
                   label = paste("wt vs het:", p_values[1]), size = 3) +
           annotate("text", x = 2, y = y_max + y_range * 0.15, 
                   label = paste("wt vs ko:", p_values[2]), size = 3) +
           annotate("text", x = 2.5, y = y_max + y_range * 0.1, 
                   label = paste("het vs ko:", p_values[3]), size = 3)
  
  # Save the plot with cleaned filename
  ggsave(
    filename = paste0("results/scRNAseq/6_monocle3/violin/pseudotime_violin_", clean_cell_type, ".pdf"),
    plot = p,
    width = 8,
    height = 6
  )
}

# combine all the violin plots into one big violoin plots. For color use the scheme below. And save it as one pdf file. 
# Define this once at the top of your script

my_cell_colors <- c(
  "NPC-G2/M" = "hotpink",
  "Neuron"          = "dodgerblue",
  "Mes-WT" = "forestgreen",
  "Mes-KO"  = "yellowgreen",
  "NPC-G1" = "purple",
  "gRG"       = "orangered",
  "IPC"        = "gold"
)

# Create one big violin plot grouped by genotype within each cell type
# Set the order of cell types
cds5_df$cell_type <- factor(cds5_df$cell_type, 
                           levels = c("NPC-G2/M", "NPC-G1", "Mes-WT", "Mes-KO", "IPC", "Neuron", "gRG"))

big_violin_plot <- ggplot(cds5_df, aes(x = interaction(genotype, cell_type), y = pseudotime, fill = cell_type)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.02, outlier.shape = NA, alpha = 0) +
  scale_fill_manual(values = my_cell_colors) +
  labs(title = "Pseudotime Distribution by Genotype Within Each Cell Type",
       x = "Genotype . Cell Type", 
       y = "Pseudotime") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the big violin plot
ggsave("results/scRNAseq/6_monocle3/violin/big_pseudotime_violin_plot.pdf", 
       plot = big_violin_plot, width = 20, height = 5)


###################################################################

# Track the expression of neuroprogenitor marker genes (e.g. PAX6, SOX2) and neuronal marker genes (e.g. MAP2, RBFOX3) along pseudotime in day20 cell types (NPC-G2/M, NPC-G1, gRG, IPC, Neuron).
# We do this for each genotype separately.

# set working directory
setwd("/nfs/turbo/umms-siwase/SHIGEKI2/rai1_hesc")

# load seurat
library(Seurat)
library(BPCells)
library(tidyverse)
library(monocle3)

# Define this once at the top of your script
my_cell_colors <- c(
  "NPC-G2/M" = "hotpink",
  "Neuron" = "dodgerblue",
  "Mes-WT" = "forestgreen",
  "Mes-KO"  = "yellowgreen",
  "NPC-G1" = "purple",
  "gRG" = "orangered",
  "IPC"        = "gold"
)

# Load the saved cds objects with pseudotime
cds5 <- load_monocle_objects(directory_path = "data/scRNAseq/iwaselab/rdata/monocle3_cds5_pseudotime")

# print the list of metadata columns in cds5
print("Metadata columns in cds5:")
print(colnames(colData(cds5)))

# check what values are in the day column
print(unique(colData(cds5)$day))
# [1] day20 day60
# Levels: day20 day60

# subset to day20 before marker analysis
day20_cells <- as.character(colData(cds5)$day) == "day20"
if (!any(day20_cells)) {
  stop("No day20 cells found in cds5.")
}
cds5_day20 <- cds5[, day20_cells]

# set marker genes 
marker_genes_requested <- c("PAX6", "SOX2", "MAP2", "RBFOX3")
marker_output_dir <- "results/scRNAseq/6_monocle3/marker_expression_pseudotime"

dir.create(marker_output_dir, recursive = TRUE, showWarnings = FALSE)

marker_cells <- colData(cds5_day20)$genotype %in% c("wt", "het", "ko") &
  !is.na(pseudotime(cds5_day20)) &
  is.finite(pseudotime(cds5_day20))

if (!any(marker_cells)) {
  stop("No valid cells found for marker pseudotime plots in day20.")
}

gene_map <- tibble(requested_gene = marker_genes_requested) %>%
  mutate(
    gene_row = dplyr::coalesce(
      rownames(cds5_day20)[match(requested_gene, rowData(cds5_day20)$gene_short_name)],
      rownames(cds5_day20)[match(requested_gene, rownames(cds5_day20))]
    )
  )

marker_cds <- cds5_day20[gene_map$gene_row, colnames(cds5_day20)[marker_cells]]
rowData(marker_cds)$gene_short_name <- gene_map$requested_gene[
  match(rownames(marker_cds), gene_map$gene_row)
]

colData(marker_cds)$genotype <- factor(
  as.character(colData(marker_cds)$genotype),
  levels = c("wt", "het", "ko"),
  labels = c("WT", "Het", "KO")
)

if (all(is.na(colData(marker_cds)$genotype))) {
  stop("No cells with WT/Het/KO genotype labels were found after filtering.")
}

combined_marker_plot <- plot_genes_in_pseudotime(
  marker_cds,
  color_cells_by = "cell_type",
  min_expr = 0.5
) +
  scale_color_manual(values = my_cell_colors, drop = FALSE) +
  facet_grid(gene_short_name ~ genotype) +
  scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
  labs(
    title = "Marker expression along pseudotime",
    subtitle = "All day20 cell types",
    x = "pseudotime",
    y = "expression"
  ) +
  theme_classic(base_size = 11) +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = file.path(marker_output_dir, "cds5_marker_expression_pseudotime_WT_Het_KO_by_gene.pdf"),
  plot = combined_marker_plot,
  width = 18,
  height = 12
)

ggsave(
  filename = file.path(marker_output_dir, "cds5_marker_expression_pseudotime_WT_Het_KO_by_gene.png"),
  plot = combined_marker_plot,
  width = 18,
  height = 12,
  dpi = 300
)

preview_plot(combined_marker_plot)
  