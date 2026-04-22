
# assign cell type
# Load libraries
library(Seurat)
library(BPCells)
library(tidyverse)
setwd("/nfs/turbo/umms-siwase/SHIGEKI2/rai1_hesc")

# load the R.data
geo_so = readRDS('data/scRNAseq/iwaselab/rdata/geo_so_sct_unintegrated_with_markers_0.1.rds')
geo_so

# make csv file for cell type - clus# Create a data frame with your cluster assignments
cluster_data <- data.frame(
  cluster = c(0, 1, 2, 3, 4, 5, 6),
  cell_type = c("NPC-G2/M", "Neuron", "Mes-WT", "Mes-KO", "NPC-G1", "gRG", "IPC")
)

# Write the data frame to a CSV file
write.csv(cluster_data, file = "results/scRNAseq/4_celltype/cluster_annotations.csv", row.names = FALSE)

# Load the annotations from prepared file
celltype_annos = read_csv('results/scRNAseq/4_celltype/cluster_annotations.csv') %>%
  mutate(cluster=factor(cluster))
head(celltype_annos)

# Merge cell types in but as a new table to slide into @meta.data
head(geo_so@meta.data)
copy_metadata = geo_so@meta.data
new_metadata = copy_metadata %>% 
  left_join(celltype_annos, by = c('seurat_clusters' = 'cluster'))

#  We are implicitly relying on the same row order!
rownames(new_metadata) = rownames(geo_so@meta.data)

# Replace the meta.data
geo_so@meta.data = new_metadata 

head(geo_so@meta.data)
unique(geo_so$cell_type)

# save the r data
saveRDS(geo_so, file = 'data/scRNAseq/iwaselab/rdata/geo_so_sct_unintegrated_with_annos.rds')

##################################################################
#### start here ##################################################
##################################################################

# load the r data
geo_so = readRDS('data/scRNAseq/iwaselab/rdata/geo_so_sct_unintegrated_with_annos.rds')

# # Visualize UMAP clusters: the assigned celltype labels (use sans font for better readability)

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

library(ggplot2)

pre_integration_umap_plot_clusters = 
  DimPlot(geo_so, 
          group.by = 'cell_type',
          label = TRUE,
          reduction = 'umap.unintegrated.sct.rpca_0.1',
          raster = TRUE,            # Keep dots rasterized
          raster.dpi = c(300, 300), 
          cols = my_cell_colors
  ) + 
  # 1. Explicitly add the Title and Axis Labels
  ggtitle("Unintegrated UMAP Clusters") + 
  labs(x = "UMAP_1", y = "UMAP_2") +
  
  # 2. Force the axis lines and text to be visible
  theme(
    # Set font to Arial, but fallback to sans if it fails
    text = element_text(family = "sans", color = "black"), 
    
    # Ensure titles and labels are black and visible
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10, color = "black"),
    
    # Explicitly draw the X and Y axis lines
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black")
  ) +
  NoLegend()

pre_integration_umap_plot_clusters

ggsave(filename = 'results/scRNAseq/2_pca/6_umap_unintegrated_sct_clusters_0.1_celltype.pdf',
       plot = pre_integration_umap_plot_clusters,
       width = 6, height = 6, units = 'in')

##############################################################
### make marker expression panel (cell type identification) ##
##############################################################

# Create lists of markers that show specific expression in each cluster

markers = list()
markers[['NPC']] = c('PAX6', 'SOX2', 'HES5') 
markers[['Mitotic']] = c('MKI67')
markers[['Mesenchimal']] = c('FN1', 'COL1A1','COL1A2','SNAI2')
markers[['Neurons']] = c('DCX', 'ELAVL2', 'ELAVL3', 'GAD1')
markers[['IPC']] = c('EOMES','NHLH1')
markers[['gRG']] = c('NFIA', 'NFIB', 'SOX9', 'S100B')

markers

# make dotplot for the markers
##########################################################################

library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)
library(patchwork)

# 1. Define your cluster-to-label mapping
cluster_data <- data.frame(
  cluster = c(0, 1, 2, 3, 4, 5, 6),
  cell_type_label = c("NPC-G2/M (0)", "Neuron (1)", "Mes-WT (2)", "Mes-KO (3)", "NPC-G1 (4)", "gRG (5)", "IPC (6)")
)

# 2. Re-define the Data Function (Fixed the rename issue)
get_dotplot_data <- function(object, assay, features, group.by) {
  DefaultAssay(object) <- assay
  clean_feats <- unique(unname(unlist(features)))
  
  df <- FetchData(object, vars = c(group.by, clean_feats))
  
  df_long <- df %>%
    rownames_to_column("Cell") %>%
    pivot_longer(cols = all_of(clean_feats), names_to = "Gene", values_to = "Expression")
  
  # Fix the rename issue by using a simpler approach
  names(df_long)[names(df_long) == group.by] <- "Group"
  
  df_summary <- df_long %>%
    group_by(Group, Gene) %>%
    summarise(
      Pct_Expressed = sum(Expression > 0) / n() * 100,
      Avg_Expression = mean(Expression),
      .groups = 'drop'
    ) %>%
    group_by(Gene) %>%
    mutate(Avg_Scaled = scale(Avg_Expression))
  
  return(df_summary)
}

# 3. Generate Data & Apply Labels
features_list <- markers
features_list

# --- SCT Data ---
sct_data <- get_dotplot_data(geo_so, assay = "SCT", features = features_list, group.by = "seurat_clusters")

# Merge with labels
# We convert 'Group' to numeric to match your cluster_data table
sct_data <- sct_data %>%
  mutate(Group_Numeric = as.numeric(as.character(Group))) %>% 
  left_join(cluster_data, by = c("Group_Numeric" = "cluster"))

# --- RNA Data ---
rna_data <- get_dotplot_data(geo_so, assay = "RNA", features = features_list, group.by = "seurat_clusters")

# Merge with labels
rna_data <- rna_data %>%
  mutate(Group_Numeric = as.numeric(as.character(Group))) %>% 
  left_join(cluster_data, by = c("Group_Numeric" = "cluster"))


# 4. Define Plotting Function (Updated to use 'cell_type_label')
make_labeled_plot <- function(data, title) {
  
  # Set order of genes
  data$Gene <- factor(data$Gene, levels = unique(unname(unlist(markers))))
  
  # Set order of Cell Types (Optional: arrange 0 to 6)
  # We use the order from your cluster_data dataframe to keep them sorted 0-6
  data$cell_type_label <- factor(data$cell_type_label, levels = cluster_data$cell_type_label)
  
  ggplot(data, aes(x = Gene, y = cell_type_label)) +
    geom_point(aes(size = Pct_Expressed, color = Avg_Scaled)) +
    scale_color_gradient(low = "lightgrey", high = "blue") +
    theme_bw() +
    labs(
      title = title,
      x = "Features", y = "Cell Type",
      color = "Avg Exp\n(Scaled)", size = "% Expressed"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color="black"),
      panel.grid.major = element_line(color = "grey90")
    )
}

# 5. Create and Combine
p1 <- make_labeled_plot(sct_data, "SCT Assay")
p2 <- make_labeled_plot(rna_data, "RNA Assay")

# 6. Save
combined_plot <- p1 + p2
ggsave("results/scRNAseq/4_celltype/marker_dotplot_labeled.pdf", plot = combined_plot, width = 14, height = 3.5)

# Display just to confirm
combined_plot

############################################################
### Make the violin plots for marker genes ################
############################################################

# Violin plots for marker genes in SCT assay
marker_genes <- unname(unlist(markers))
violin_plot <- VlnPlot(geo_so, features = marker_genes, assay = "SCT", group.by = "cell_type", pt.size = 0, col = my_cell_colors) +
  theme(
    text = element_text(family = "sans", color = "black"),
    plot.title = element_text(hjust = 0.5, size = 14),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10, color = "black")
  )
violin_plot
# Save the violin plot
ggsave("results/scRNAseq/4_celltype/marker_violin_plots.pdf", plot = violin_plot, width = 12, height = 8)   

############################################################
#### Calculate mimimal and maximal number genes detected in each cell (nFeature_RNA) ################ 

# load the R.data
geo_so = readRDS('data/scRNAseq/iwaselab/rdata/geo_so_sct_unintegrated_with_annos.rds')
head(geo_so@meta.data)
geo_so

# calculate the number of genes expressed in each cell (nFeature_RNA) with more than 1 read

# 1. Access the raw RNA counts matrix
# We use the 'RNA' assay and 'counts' layer for the truest representation of reads
counts_matrix <- Seurat::GetAssayData(geo_so, assay = "RNA", layer = "counts")

# 2. Calculate the number of genes with > 0 reads for each cell
# rowSums is faster, but since we want genes per cell, we operate on columns
genes_per_cell <- Matrix::colSums(counts_matrix > 0)

# 3. Get the minimum and maximum
min_genes <- min(genes_per_cell)
max_genes <- max(genes_per_cell)

# Display results
cat("Min genes detected:", min_genes, "\n")
cat("Max genes detected:", max_genes, "\n")

#############################################################
### Analyze cell type composition between genotypes #########
#############################################################

# get cell type proportions for each sample
celltype_proportions <- prop.table(table(geo_so$orig.ident, geo_so$cell_type), margin = 1)
celltype_proportions_df <- as.data.frame(celltype_proportions)
head(celltype_proportions_df)
colnames(celltype_proportions_df) <- c("Sample", "CellType", "Proportion")

# add genotype information by splitting on underscore and extracting genotype
celltype_proportions_df$Genotype <- sapply(strsplit(as.character(celltype_proportions_df$Sample), "_"), "[[", 2)

# add timepoint information as well
celltype_proportions_df$Timepoint <- sapply(strsplit(as.character(celltype_proportions_df$Sample), "_"), "[[", 1)

# reorder genotype factor levels: WT, Het, KO
celltype_proportions_df$Genotype <- factor(celltype_proportions_df$Genotype, levels = c("wt", "het", "ko"))

# save the dataframe
write.csv(celltype_proportions_df, file = "results/scRNAseq/4_celltype/celltype_proportions.csv", row.names = FALSE)

# read the dataframe and re-apply factor ordering
celltype_proportions_df <- read.csv("results/scRNAseq/4_celltype/celltype_proportions.csv")

# calculate the percentatage of each cell type for the three genotypes
celltype_proportions_df$Percentage <- celltype_proportions_df$Proportion * 100

# take the average of the two replicates 
celltype_proportions_avg <- celltype_proportions_df %>%
  group_by(Genotype, CellType) %>%
  summarize(Percentage = mean(Percentage))

celltype_proportions_avg

# scale that percentage such that when all cell types are summed for a genotype, it equals 100%
celltype_proportions_avg <- celltype_proportions_avg %>%
  group_by(Genotype) %>%
  mutate(Percentage = Percentage / sum(Percentage) * 100) 

# re-apply the genotype factor ordering after reading from CSV
celltype_proportions_avg$Genotype <- factor(celltype_proportions_avg$Genotype, levels = c("wt", "het", "ko"))

# make stacked bar plot of cell type proportions (%) by genotype
barplot_celltype_proportion <- ggplot(celltype_proportions_avg, 
                                      aes(x = Genotype, y = Percentage, fill = CellType)) +
  geom_bar(stat = "identity") +
  # This is where you apply your custom color vector
  scale_fill_manual(values = my_cell_colors) + 
  theme_minimal() +
  labs(title = "Cell Type Proportions by Genotype", x = "Genotype", y = "Percentage (%)") +
  scale_x_discrete(labels = c("wt" = "WT", "het" = "Het", "ko" = "KO")) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = c(0, 0))

# Display the plot
barplot_celltype_proportion

# save the plot as a PDF file
ggsave(filename = "results/scRNAseq/4_celltype/barplot_celltype_proportion.pdf",
       plot = barplot_celltype_proportion,
       width = 8, height = 6, units = "in")

#################################################################
### cell type compostion analysis: two day points separately ####
#################################################################

library(dplyr)
library(ggplot2)

# 1. Take the average percentage, grouping by both Genotype AND Timepoint
celltype_proportions_avg <- celltype_proportions_df %>%
  group_by(Timepoint, Genotype, CellType) %>%
  summarize(Percentage = mean(Percentage), .groups = 'drop')

# 2. Scale within each (Timepoint + Genotype) combination to ensure they sum to 100%
celltype_proportions_avg <- celltype_proportions_avg %>%
  group_by(Timepoint, Genotype) %>%
  mutate(Percentage = Percentage / sum(Percentage) * 100) %>%
  ungroup()

# 3. Ensure Genotype factor ordering is preserved
celltype_proportions_avg$Genotype <- factor(celltype_proportions_avg$Genotype, levels = c("wt", "het", "ko"))

# 4. Create the plot using facet_wrap to separate Day 20 and Day 60
barplot_celltype_proportion <- ggplot(celltype_proportions_avg, 
                                      aes(x = Genotype, y = Percentage, fill = CellType)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = my_cell_colors) + 
  facet_wrap(~Timepoint) +  # This creates two side-by-side plots
  theme_minimal() +
  labs(title = "Cell Type Proportions by Genotype and Timepoint", 
       x = "Genotype", 
       y = "Percentage (%)") +
  scale_x_discrete(labels = c("wt" = "WT", "het" = "Het", "ko" = "KO")) +
  scale_y_continuous(limits = c(0, 100.1), breaks = seq(0, 100, 20), expand = c(0, 0)) +
  theme(panel.spacing = unit(2, "lines")) # Adds a bit of space between the two plots

# Display
print(barplot_celltype_proportion)

# 5. Save the faceted plot
ggsave(filename = "results/scRNAseq/4_celltype/barplot_celltype_proportion_day_split.pdf",
       plot = barplot_celltype_proportion,
       width = 10, height = 6, units = "in")

############################################################
#### make replicate dimention plots with new color ########
############################################################

pre_integration_umap_plot_replicate_split = 
  DimPlot(geo_so,
          group.by = 'cell_type',
          split.by = 'replicate',
          label = FALSE,
          cols = my_cell_colors,
          raster = TRUE,           # Keep dots rasterized
          raster.dpi = c(300, 300),
          reduction = 'umap.unintegrated.sct.rpca_0.1') + 
  # 1. Explicitly add the Title and Axis Labels
  labs(x = "UMAP_1", y = "UMAP_2") +
  
  # 2. Force the axis lines and text to be visible
  theme(
    text = element_text(family = "sans", color = "black"), 
    axis.line = element_line(color = "black", size = 0.5),
    axis.ticks = element_line(color = "black")
  ) +
  NoLegend() +
  # 3. Force vertical alignment
  facet_wrap(~replicate, ncol = 1) 

pre_integration_umap_plot_replicate_split


# save the plot as pdf
# save the plot as a PDF file
ggsave(filename = "results/scRNAseq/4_celltype/umap_replicates.pdf",
       plot = pre_integration_umap_plot_replicate_split,
       width = 4, height = 8, units = "in")

############################################################
#### make demention plots of genotypes and day with assigned cell types 
############################################################

# Define the desired order of your combined factor levels
# Day 20 plots on the top row (WT, HET, KO)
# Day 60 plots on the bottom row (WT, HET, KO)

# Create a new combined column
geo_so$day_genotype <- paste(geo_so$day, geo_so$genotype, sep = "_")

head(geo_so@meta.data)

desired_order <- c(
  "day20_wt", "day20_het", "day20_ko", # Top row
  "day60_wt", "day60_het", "day60_ko"  # Bottom row
)

# Convert your combined metadata column to a factor with the specified levels
geo_so$day_genotype <- factor(geo_so$day_genotype, levels = desired_order)

# Verify the levels (optional, but good for debugging)
levels(geo_so$day_genotype)


pre_integration_umap_plot_day_genotype_split =
  DimPlot(
    geo_so,
    group.by = "cell_type",
    split.by = "day_genotype",
    label = TRUE,
    reduction = "umap.unintegrated.sct.rpca_0.1",
    cols = my_cell_colors, 
    raster = TRUE,            # Keep dots rasterized
    raster.dpi = c(300, 300),
    ncol = 3
  ) + 
  # 1. Explicitly add the Title and Axis Labels
  labs(x = "UMAP_1", y = "UMAP_2") +
  # 2. Force the axis lines and text to be visible
  theme(
    text = element_text(family = "sans", color = "black")
  ) +
  NoLegend() 

pre_integration_umap_plot_day_genotype_split

ggsave(filename = 'results/scRNAseq/4_celltype/umap_genotype_day_split_0.1.pdf',
       plot = pre_integration_umap_plot_day_genotype_split,
       width = 12, height = 8, units = 'in')


