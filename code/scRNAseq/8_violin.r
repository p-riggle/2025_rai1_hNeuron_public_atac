# violin plots of gene of interest between clusters
# Load libraries
# Enter commands in R (or R studio, if installed)
library(Seurat)
library(BPCells)
library(tidyverse)

# load the R.data
geo_so = readRDS('data/scRNAseq/iwaselab/rdata/geo_so_sct_unintegrated_with_markers_0.1.rds')
geo_so

# plot RAI1 level

VlnPlot(
  object = geo_so,
  features = "RAI1",
  group.by = "seurat_clusters", # Use the default clustering results
  pt.size = 0, # Setting pt.size to 0 removes all individual data points
  layer = "scale.data"
) +
  # Use ggplot2 functions to customize the plot's appearance
  ggtitle("RAI1 Expression by Cluster") +
  theme(plot.title = element_text(hjust = 0.5)) + # Center the plot title
  NoLegend() + # Remove the legend as the clusters are already labeled on the x-axis
  
  # This line adds a log10 scale to the y-axis, which will expand the
  # view of the lower expression values and compress the higher ones.
  # We add 1 to the expression values to avoid taking the log of zero,
  # which is mathematically undefined.
  scale_y_log10()


viol_rai1 <- VlnPlot(
  object = geo_so,
  features = "RAI1",
  group.by = "seurat_clusters", # Use the default clustering results
  pt.size = 0, # Setting pt.size to 0 removes all individual data points
  layer = "scale.data"
) +
  # Use ggplot2 functions to customize the plot's appearance
  ggtitle("RAI1 Expression by Cluster") +
  theme(plot.title = element_text(hjust = 0.5)) + # Center the plot title
  NoLegend() 

ggsave(filename = "results/scRNAseq/4_pseudobulk/viol_rai1.pdf", plot = viol_rai1, width = 5, height = 4, units = "in")

# plot BDNF beween genotypes
VlnPlot(
  object = geo_so,
  features = "BDNF",
  group.by = "seurat_clusters", # Use the default clustering results
  split.by = "genotype", # This is the new argument to split the plots by genotype
  pt.size = 0, # Setting pt.size to 0 removes all individual data points
  layer = "scale.data", # Plot the continuous, scaled data
  combine = TRUE # Return a list of plots instead of a single combined plot
) +
  # Use ggplot2 functions to customize the plot's appearance
  ggtitle("BDNF Expression by Cluster") +
  theme(plot.title = element_text(hjust = 0.5)) + # Center the plot title
  NoLegend() + 
  scale_y_log10()

# 1. Reorder the genotypes
# We use the factor() function to set the specific order of the levels.
# This ensures that 'wt', 'het', and 'ko' appear in the desired order on the plot.
geo_so$genotype <- factor(geo_so$genotype, levels = c("wt", "het", "ko"))

# Now, generate the violin plot with the desired modifications
VlnPlot(
  object = geo_so,
  features = "BDNF",
  group.by = "seurat_clusters", # Use the default clustering results
  split.by = "genotype", # This is the new argument to split the plots by genotype
  pt.size = 0, # Setting pt.size to 0 removes all individual data points
  layer = "scale.data", # Plot the continuous, scaled data
  combine = TRUE # Return a single combined plot
) +
  # Use ggplot2 functions to customize the plot's appearance
  ggtitle("BDNF Expression by Cluster") +
  theme(plot.title = element_text(hjust = 0.5)) + # Center the plot title
  # 2. Show the legend by removing NoLegend()
  
  # 3. Change the colors of the violin plots
  # We use scale_fill_manual to assign specific colors to each genotype.
  # The colors are matched to the levels we set earlier: wt, het, and ko.
  scale_fill_manual(values = c("wt" = "darkblue", "het" = "blue", "ko" = "lightblue")) +
  
  # Add the log scale to the y-axis
  scale_y_log10()

# Now, generate the violin plot with the desired modifications
viol_bdnf <- VlnPlot(
  object = geo_so,
  features = "BDNF",
  group.by = "seurat_clusters", # Use the default clustering results
  split.by = "genotype", # This is the new argument to split the plots by genotype
  pt.size = 0, # Setting pt.size to 0 removes all individual data points
  layer = "scale.data", # Plot the continuous, scaled data
  combine = TRUE # Return a single combined plot
) +
  # Use ggplot2 functions to customize the plot's appearance
  ggtitle("BDNF Expression by Cluster") +
  theme(plot.title = element_text(hjust = 0.5)) + # Center the plot title
  # 2. Show the legend by removing NoLegend()
  
  # 3. Change the colors of the violin plots
  # We use scale_fill_manual to assign specific colors to each genotype.
  # The colors are matched to the levels we set earlier: wt, het, and ko.
  scale_fill_manual(values = c("wt" = "darkblue", "het" = "blue", "ko" = "lightblue")) +
  
  # Add the log scale to the y-axis
  scale_y_log10() +
  
  # Add geom_boxplot to overlay the boxplots on top of the violin plots
  geom_boxplot(width = 0.5, outlier.shape = NA, position = position_dodge(0.9))

ggsave(filename = "results/scRNAseq/4_pseudobulk/viol_bdnf.pdf", plot = viol_bdnf, width = 6, height = 4, units = "in")

