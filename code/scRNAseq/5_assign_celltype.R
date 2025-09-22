
# assign cell type
# Load libraries
library(Seurat)
library(BPCells)
library(tidyverse)

# load the R.data
geo_so = readRDS('data/scRNAseq/iwaselab/rdata/geo_so_sct_unintegrated_with_markers_0.1.rds')
geo_so

# make csv file for cell type - clus# Create a data frame with your cluster assignments
cluster_data <- data.frame(
  cluster = c(0, 1, 2, 3, 4, 5, 6),
  cell_type = c("NPC", "Neuron", "T2", "T3", "NPC-G1", "OPC", "IPC")
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

# Visualize UMAP clusters: the assigned celltype labels
pre_integration_umap_plot_clusters = 
  DimPlot(geo_so, 
          group.by = 'cell_type.x',
          label = TRUE,
          reduction = 'umap.unintegrated.sct.rpca_0.1') + 
  NoLegend()
pre_integration_umap_plot_clusters

ggsave(filename = 'results/scRNAseq/2_pca/5_umap_unintegrated_sct_clusters_0.1.png',
       plot = pre_integration_umap_plot_clusters,
       width = 6, height = 6, units = 'in')
