
# Load libraries
library(Seurat)
library(BPCells)
library(tidyverse)
# =========================================================================
# PCA and Integration
# =========================================================================
geo_so = readRDS('data/scRNAseq/iwaselab/rdata/geo_so_sct_normalized.rds')

# Build a PCA and add it to the Seurat object
geo_so = RunPCA(geo_so, reduction.name = 'unintegrated.sct.pca')
geo_so

pcs = 13 #rebecca used 19

# =========================================================================
# Clustering and Projection
# =========================================================================

# Create KNN graph for PCs with `FindNeighbors()`
geo_so = FindNeighbors(geo_so, dims = 1:pcs, reduction = 'unintegrated.sct.pca')

# Then generate clusters
geo_so = FindClusters(geo_so,
                      resolution = 0.1,
                      cluster.name = 'unintegrated.sct.rpca.clusters_0.1')

head(geo_so@meta.data)

# -------------------------------------------------------------------------
# Create UMAP reduction
geo_so = RunUMAP(geo_so,
                 dims = 1:pcs,
                 reduction = 'unintegrated.sct.pca',
                 reduction.name = 'umap.unintegrated.sct.rpca_0.1')


# Note a third reduction has been added: `umap.unintegrated.sct.rpca`
geo_so 

# save R data with clustering and dementional reduction

saveRDS(geo_so, file = 'data/scRNAseq/iwaselab/rdata/geo_so_sct_clustered_0.1.rds')

geo_so = readRDS('data/scRNAseq/iwaselab/rdata/geo_so_sct_clustered_0.1.rds')


# -------------------------------------------------------------------------
# Visualize UMAP clusters: ID labels
pre_integration_umap_plot_clusters = 
  DimPlot(geo_so, 
          group.by = 'seurat_clusters',
          label = TRUE,
          reduction = 'umap.unintegrated.sct.rpca_0.1') + 
  NoLegend()
pre_integration_umap_plot_clusters

ggsave(filename = 'results/scRNAseq/2_pca/5_umap_unintegrated_sct_clusters_0.1.png',
       plot = pre_integration_umap_plot_clusters,
       width = 6, height = 6, units = 'in')
# -------------------------------------------------------------------------
# Visualize UMAP clusters: clusters with labels, split by condition
pre_integration_umap_plot_split_clusters = 
  DimPlot(geo_so,
          group.by = 'seurat_clusters',
          split.by = 'day',
          label = TRUE,
          reduction = 'umap.unintegrated.sct.rpca_0.1') + 
  NoLegend()
pre_integration_umap_plot_split_clusters

ggsave(filename = 'results/scRNAseq/2_pca/6_umap_unintegrated_sct_split_clusters_0.1.png',
       plot = pre_integration_umap_plot_split_clusters, 
       width = 14, height = 6, units = 'in')
dev.off()  
# Visualize UMAP clusters: genotype labels
pre_integration_umap_plot_genotype = 
  DimPlot(geo_so,
          group.by = 'genotype',
          label = FALSE,
          reduction = 'umap.unintegrated.sct.rpca_0.1')
pre_integration_umap_plot_genotype

ggsave(filename = 'results/scRNAseq/2_pca/8_umap_unintegrated_sct_genotype_0.1.png',
       plot = pre_integration_umap_plot_genotype,
       width = 8, height = 6, units = 'in')

# Visualize UMAP clusters: split by genotype labels
pre_integration_umap_plot_genotype_split = 
  DimPlot(geo_so,
          group.by = 'seurat_clusters',
          split.by = 'genotype',
          label = FALSE,
          reduction = 'umap.unintegrated.sct.rpca_0.1')
pre_integration_umap_plot_genotype_split

ggsave(filename = 'results/scRNAseq/2_pca/8_umap_unintegrated_sct_genotype_split_0.1.png',
       plot = pre_integration_umap_plot_genotype_split,
       width = 12, height = 6, units = 'in')

dev.off()
# Visualize UMAP clusters: split by genotype and day labels

# Assuming 'day' and 'genotype' are already in geo_so@meta.data
# If not, make sure they are added first.
# Example: geo_so$day <- some_vector_of_day_values

# Create a new combined column
geo_so$day_genotype <- paste(geo_so$day, geo_so$genotype, sep = "_")

head(geo_so@meta.data)

pre_integration_umap_plot_genotype_day_split = 
  DimPlot(geo_so,
          group.by = 'seurat_clusters',
          split.by = 'day_genotype',
          label = FALSE,
          reduction = 'umap.unintegrated.sct.rpca_0.1')
pre_integration_umap_plot_genotype_day_split

ggsave(filename = 'results/scRNAseq/2_pca/8_umap_unintegrated_sct_genotype_split_0.1.png',
       plot = pre_integration_umap_plot_genotype_split,
       width = 12, height = 6, units = 'in')

# Define the desired order of your combined factor levels
# Day 20 plots on the top row (WT, HET, KO)
# Day 60 plots on the bottom row (WT, HET, KO)
desired_order <- c(
  "day20_wt", "day20_het", "day20_ko", # Top row
  "day60_wt", "day60_het", "day60_ko"  # Bottom row
)

# Convert your combined metadata column to a factor with the specified levels
geo_so$day_genotype <- factor(geo_so$day_genotype, levels = desired_order)

# Verify the levels (optional, but good for debugging)
levels(geo_so$day_genotype)

# -------------------------------------------------------------------------
# Visualize UMAP clusters with the custom split order and layout
# -------------------------------------------------------------------------
pre_integration_umap_plot_day_genotype_split =
  DimPlot(
    geo_so,
    group.by = "seurat_clusters",
    split.by = "day_genotype",
    label = FALSE,
    reduction = "umap.unintegrated.sct.rpca_0.1",
    ncol = 3
  )

pre_integration_umap_plot_day_genotype_split

ggsave(filename = 'results/scRNAseq/2_pca/8_umap_unintegrated_sct_genotype_day_split_0.1.png',
       plot = pre_integration_umap_plot_day_genotype_split,
       width = 12, height = 8, units = 'in')

# separate by replicate to see if any outliers generated the new cluster in KO
unique(geo_so@meta.data$replicate)

pre_integration_umap_plot_replicate_split = 
  DimPlot(geo_so,
          group.by = 'seurat_clusters',
          split.by = 'replicate',
          label = FALSE,
          reduction = 'umap.unintegrated.sct.rpca_0.1')
pre_integration_umap_plot_replicate_split

ggsave(filename = 'results/scRNAseq/2_pca/8_umap_unintegrated_sct_replicate_split_0.1.png',
       plot = pre_integration_umap_plot_replicate_split,
       width = 12, height = 6, units = 'in')

# =========================================================================
# Marker identification and visualization
# =========================================================================

geo_so = readRDS('data/scRNAseq/iwaselab/rdata/geo_so_sct_clustered_0.1.rds')

geo_so

# Prep for cluster comparisons to find empirical markers
colnames(geo_so@meta.data)
geo_so = SetIdent(geo_so, value = 'unintegrated.sct.rpca.clusters_0.1')
geo_so = PrepSCTFindMarkers(geo_so)

# Run comparisons for each cluster to generate markers
geo_markers = FindAllMarkers(geo_so, only.pos = TRUE, min.pct = 0.05)

# Write out full cluster marker results to file
write_csv(geo_markers, file = 'results/scRNAseq/3_markers/marker_genes_0.1res.csv')

# write out full cluster marker results into a excel file with each cluster in a separate sheet
library(openxlsx)
# Create a new workbook
wb <- createWorkbook()
# Get unique cluster identities
clusters <- unique(geo_markers$cluster)
# Loop through each cluster and add a worksheet
for (clust in clusters) {
  # Filter markers for the current cluster
  cluster_markers <- geo_markers %>% filter(cluster == clust)
  # Add a new worksheet named after the cluster
  addWorksheet(wb, paste0("Cluster_", clust))
  # Write the cluster markers to the worksheet
  writeData(wb, sheet = paste0("Cluster_", clust), cluster_markers)
}
# Save the workbook to a file
saveWorkbook(wb, file = 'results/scRNAseq/3_markers/marker_genes_0.1res.xlsx', overwrite = TRUE)
# -------------------------------------------------------------------------

# Take a look at the first few rows of the result
head(geo_markers)
# -------------------------------------------------------------------------
# Identify marker genes for each cluster
# Create table of top 5 markers per cluster (using default ranking)
top_5 = geo_markers %>% filter(p_val_adj < 0.01) %>% group_by(cluster) %>% slice_head(n = 5)

# Look at results
head(top_5, n = 10)# -------------------------------------------------------------------------
# Visualize top marker genes as dot plot
top_5_sct_dot_plot = DotPlot(geo_so, features = unique(top_5$gene)) + 
  theme(panel.background = element_rect(fill = "white"), axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  labs(title = 'Top 5 Cluster Genes by FDR and avg_log2FC') + coord_flip()
top_5_sct_dot_plot
# -------------------------------------------------------------------------
# Save dot plot of top marker genes
ggsave(filename = 'results/scRNAseq/3_markers/markers_top_5_sct_dot_plot_0.1.png', 
       plot = top_5_sct_dot_plot, 
       width = 8, height = 18, units = 'in') 

# -------------------------------------------------------------------------
# Save Seurat object and gene marker data
saveRDS(geo_so, file = 'data/scRNAseq/iwaselab/rdata/geo_so_sct_unintegrated_with_markers_0.1.rds')
saveRDS(geo_markers, file = 'data/scRNAseq/iwaselab/rdata/geo_markers_0.1.rds')


# =========================================================================
# Cell Type Annotation
# =========================================================================

# Load scCATCH
library(scCATCH)

# check that cell identities are set to expected resolution 
all(Idents(geo_so) == geo_so$unintegrated.sct.rpca.clusters_0.1)

# -------------------------------------------------------------------------
# Annotate clusters using scCATCH

# create scCATCH object, using count data
geo_catch = createscCATCH(data = geo_so@assays$SCT@counts, cluster = as.character(Idents(geo_so)))

# add marker genes to use for predictions
geo_catch@markergene = geo_markers

# Access the cellmatch dataset
data(cellmatch)

# Get unique tissue types from the cellmatch dataset
available_tissues <- unique(cellmatch$tissue)

# Print the list of available tissue types
print(available_tissues)

# specify tissues/cell-types from the scCATCH reference
tissues_of_interest = c('Brain', 'Embryonic brain', 'Midbrain', 'Embryonic prefrontal cortex', 'Hippocampus', 'Cerebellum', 'Neural tube')
species_tissue_mask = cellmatch$species == 'Human' & cellmatch$tissue %in% tissues_of_interest
geo_catch@marker = cellmatch[species_tissue_mask, ]

# run scCATCH to generate predictions
geo_catch = findcelltype(geo_catch)

# look at the predictions
sc_pred_tbl = geo_catch@celltype %>% 
  select(cluster, cell_type, celltype_score) %>%
  rename(cell_type_prediction = cell_type)
sc_pred_tbl

write_csv(sc_pred_tbl, file = 'results/scRNAseq/3_markers/scCATCH_cluster_predictions_0.1.csv')

# -------------------------------------------------------------------------
# Check mitochondrial gene expression
percent_mito_plot = FeaturePlot(geo_so, features='percent.mt')
percent_mito_plot

ggsave(filename = 'results/scRNAseq/3_markers/mitoGenes_0.1.png', 
       plot = percent_mito_plot, 
       width = 8, height = 8, units = 'in') 

# check RAI1 expression

rai1_plot = FeaturePlot(geo_so, features = 'RAI1')
rai1_plot

ggsave(filename = 'results/scRNAseq/3_markers/rai1_0.1.png', 
       plot = rai1_plot, 
       width = 8, height = 8, units = 'in') 
