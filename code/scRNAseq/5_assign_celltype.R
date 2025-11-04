
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
  cell_type = c("NPC-G2/M", "Neuron", "Mes-WT", "Mes-Mut", "NPC-G1", "Astrocyte", "IPC")
)

# Write the data frame to a CSV file
write.csv(cluster_data, file = "results/scRNAseq/5_celltype/cluster_annotations.csv", row.names = FALSE)

# Load the annotations from prepared file
celltype_annos = read_csv('results/scRNAseq/5_celltype/cluster_annotations.csv') %>%
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

# load the r data
geo_so = readRDS('data/scRNAseq/iwaselab/rdata/geo_so_sct_unintegrated_with_annos.rds')

# Visualize UMAP clusters: the assigned celltype labels (use Arial font)
pre_integration_umap_plot_clusters = 
  DimPlot(geo_so, 
          group.by = 'cell_type',
          label = TRUE,
          reduction = 'umap.unintegrated.sct.rpca_0.1') + 
  NoLegend() +
  theme(
    text = element_text(family = "Arial"),
    plot.title = element_text(family = "Arial"),
    axis.title = element_text(family = "Arial"),
    axis.text = element_text(family = "Arial"),
    legend.text = element_text(family = "Arial")
  )

pre_integration_umap_plot_clusters

ggsave(filename = 'results/scRNAseq/2_pca/5_umap_unintegrated_sct_clusters_0.1.png',
       plot = pre_integration_umap_plot_clusters,
       width = 6, height = 6, units = 'in')


#############################################################
### Analyze cell type composition between genotypes #########
#############################################################

# load the R.data
geo_so = readRDS('data/scRNAseq/iwaselab/rdata/geo_so_sct_unintegrated_with_annos.rds')
head(geo_so@meta.data)

# get cell type proportions for each sample
celltype_proportions <- prop.table(table(geo_so$orig.ident, geo_so$cell_type), margin = 1)
celltype_proportions_df <- as.data.frame(celltype_proportions)
colnames(celltype_proportions_df) <- c("Sample", "CellType", "Proportion")

# add genotype information by splitting on underscore and extracting genotype
celltype_proportions_df$Genotype <- sapply(strsplit(as.character(celltype_proportions_df$Sample), "_"), "[[", 2)

# add timepoint information as well
celltype_proportions_df$Timepoint <- sapply(strsplit(as.character(celltype_proportions_df$Sample), "_"), "[[", 1)

# reorder genotype factor levels: WT, Het, KO
celltype_proportions_df$Genotype <- factor(celltype_proportions_df$Genotype, levels = c("wt", "het", "ko"))

# save the dataframe
write.csv(celltype_proportions_df, file = "results/scRNAseq/5_celltype/celltype_proportions.csv", row.names = FALSE)

# read the dataframe and re-apply factor ordering
celltype_proportions_df <- read.csv("results/scRNAseq/5_celltype/celltype_proportions.csv")

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
barplot_celltype_proportion <- ggplot(celltype_proportions_avg, aes(x = Genotype, y = Percentage, fill = CellType)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Cell Type Proportions by Genotype", x = "Genotype", y = "Percentage (%)") +
  scale_x_discrete(labels = c("wt" = "WT", "het" = "Het", "ko" = "KO")) +  # Clean labels for display
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 20), expand = c(0, 0))  # Force y-axis to 0-100% with no expansion

barplot_celltype_proportion

# save the plot as a PDF file
ggsave(filename = "results/scRNAseq/5_celltype/barplot_celltype_proportion.pdf",
       plot = barplot_celltype_proportion,
       width = 8, height = 6, units = "in")

