# violin plots of gene of interest between clusters
# Load libraries
# Enter commands in R (or R studio, if installed)
library(Seurat)
library(BPCells)
library(tidyverse)

# load the R.data
geo_so = readRDS('data/scRNAseq/iwaselab/rdata/geo_so_sct_unintegrated_with_annos.rds')
unique(geo_so$cell_type)
head(geo_so@meta.data)
geo_so

# plot RAI1 level
# 1. Reorder the cell_types
# We use the factor() function to set the specific order of the levels.

# define colors
my_cell_colors <- c(
  "NPC-G2/M" = "hotpink",
  "Neuron"          = "dodgerblue",
  "Mes-WT" = "forestgreen",
  "Mes-KO"  = "yellowgreen",
  "NPC-G1" = "purple",
  "Astrocyte"       = "orangered",
  "IPC"        = "gold"
)

geo_so$cell_type <- factor(geo_so$cell_type, levels = c("NPC-G2/M", "NPC-G1", "IPC", "Neuron", "Mes-WT", "Mes-KO", "Astrocyte"))
head(geo_so@meta.data)

viol_rai1 <- VlnPlot(
  object = geo_so,
  features = "RAI1",
  group.by = "cell_type", # Use the default clustering results
  pt.size = 0.1, alpha = 0.1 # Setting pt.size to 0 removes all individual data points
) +
  # Use ggplot2 functions to customize the plot's appearance
  ggtitle("RAI1 Expression by Cluster") +
  theme(plot.title = element_text(hjust = 0.5)) + # Center the plot title
  NoLegend() 

viol_rai1

# make violin plot using "SCT" instead of "RNA" data
create_gene_vlnplot_log <- function(geo_so, gene_name) {
  # Generate the violin plot using the input gene name
  vln_plot <- VlnPlot(
    object = geo_so,
    features = gene_name,         # Use the function's gene_name argument
    group.by = "cell_type",
    pt.size = 0.1, alpha = 0.005,
    slot = "scale.data",  # Use slot instead of layer
    combine = TRUE, 
    col = my_cell_colors,
  ) +
    # Customize the plot using ggplot2 layers
    ggtitle(paste(gene_name, "Expression by Cluster")) + # Create a dynamic title
    theme(plot.title = element_text(hjust = 0.5)) 
    # stat_summary(
    #   fun = median,
    #   geom = "crossbar",
    #   width = 0.5,
    #   col = "black",
    #   position = position_dodge(0.9)
    # )
  
  # Return the final plot object
  return(vln_plot)
}
RAI1_violin_log <- create_gene_vlnplot_log(geo_so, "RAI1")
RAI1_violin_log

# save the plot as pdf
ggsave("results/scRNAseq/4_celltype/RAI1_violin_log.pdf", plot = RAI1_violin_log)

# expression is low. lets make dot plot
# 1. Subset the Seurat object to keep only the 'wt' genotype cells
wt_so <- subset(geo_so, subset = genotype == "wt")

# 2. Create the DotPlot for the RAI1 gene, grouped by cell type
DotPlot(wt_so, features = "RAI1", group.by = "cell_type") + RotatedAxis()




# plot BDNF beween genotypes
# 1. Reorder the genotypes
# We use the factor() function to set the specific order of the levels.
# This ensures that 'wt', 'het', and 'ko' appear in the desired order on the plot.
geo_so$genotype <- factor(geo_so$genotype, levels = c("wt", "het", "ko"))

# Now, generate the violin plot with the desired modifications
# This step should be done once before using the function
# with linear scale 
create_gene_vlnplot <- function(geo_so, gene_name) {
  # Generate the violin plot using the input gene name
  vln_plot <- VlnPlot(
    object = geo_so,
    features = gene_name,         # Use the function's gene_name argument
    group.by = "cell_type",
    split.by = "genotype",
    pt.size = 0.1, alpha = 0.1,
    slot = "data",  # Use slot instead of layer
    combine = TRUE
  ) +
    # Customize the plot using ggplot2 layers
    ggtitle(paste(gene_name, "Expression by Cluster")) + # Create a dynamic title
    theme(plot.title = element_text(hjust = 0.5)) +
    scale_fill_manual(values = c("wt" = "darkgreen", "het" = "lightgreen", "ko" = "aquamarine")) +
    stat_summary(
      fun = median,
      geom = "crossbar",
      width = 0.5,
      color = "black",
      position = position_dodge(0.9)
    )
  
  # Return the final plot object
  return(vln_plot)
}

BDNF_violin <- create_gene_vlnplot(geo_so, "BDNF")
BDNF_violin

# expression low, lets make a dot plot
# Create the dot plot for BDNF, grouping by cell type and splitting by genotype
# Corrected command with the 'cols' parameter specified
# Add the 'cols' parameter to define a color gradient for expression
# Use a three-color gradient to satisfy the error message
# Make sure ggplot2 is loaded (usually loaded with Seurat)
# Make sure ggplot2 is loaded
library(ggplot2)

neuron_deg <- c("GRM5", "GRIN2D", "GRIA4")

DotPlot(geo_so, 
        features = neuron_deg, group.by = "cell_type", split.by = "genotype", cols = c("green", "green", "green")) + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        labs(title = '') + 
        coord_flip()

SOST_violin_log <- create_gene_vlnplot_log(geo_so, "SOST")
SOST_violin <- create_gene_vlnplot(geo_so, "SOST")

SOST_violin_log
SOST_violin

adcy2_violin_log <- create_gene_vlnplot_log(geo_so, "ADCY2")
adcy2_violin <- create_gene_vlnplot(geo_so, "ADCY2")

adcy2_violin_log
adcy2_violin

grm5_violin <- create_gene_vlnplot_log(geo_so, "GRM5")
grm5_violin_log <- create_gene_vlnplot(geo_so, "GRM5")

grm5_violin
grm5_violin_log 


gria4_violin <- create_gene_vlnplot_log(geo_so, "GRIA4")
gria4_violin_log <- create_gene_vlnplot(geo_so, "GRIA4")

gria4_violin
gria4_violin_log 

# Make sure to replace "YourGeneName" with the actual gene (e.g., "Adcy2")
# Add the 'cols' parameter to define the color scale
# The code to prepare your ordered labels remains the same
# geo_so$celltype_genotype <- factor(geo_so$celltype_genotype, levels = existing_labels)

DotPlot(
  geo_so,
  features = "ADCY2",
  group.by = "celltype_genotype", # Use our ordered column
  cols = c("lightgrey", "red")
) + 
  # This line rotates the x-axis labels to prevent them from overlapping
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  xlab("Cell Type and Genotype") # Optional: Rename the x-axis label
