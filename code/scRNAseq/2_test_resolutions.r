
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


# Note a third reduction has been added: `umap.integrated.sct.rpca`
geo_so 

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
pre_integration_umap_plot_genotype_split = 
  DimPlot(geo_so,
          group.by = 'seurat_clusters',
          split.by = 'genotype', 'day',
          label = FALSE,
          reduction = 'umap.unintegrated.sct.rpca_0.1')
pre_integration_umap_plot_genotype_day_split

ggsave(filename = 'results/scRNAseq/2_pca/8_umap_unintegrated_sct_genotype_split_0.1.png',
       plot = pre_integration_umap_plot_genotype_split,
       width = 12, height = 6, units = 'in')

