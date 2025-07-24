
# Load libraries
library(Seurat)
library(BPCells)
library(tidyverse)

options(future.globals.maxSize = 1e9)

# -------------------------------------------------------------------------
# To load data from 10X Cell Ranger

# Collect the input directories
#   Each sample dir contains barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz.
#   Naming the sample_dirs vector makes Seurat name the
#     samples in the corresponding manner, which is nice for us.
sample_dirs = list(
  day20_wt_rep1  = "data/scRNAseq/10x_analysis_11534-BZ/Sample_11534-BZ-Pool01/per_sample_outs/11534-BZ-2/count/sample_filtered_feature_bc_matrix", 
  day20_wt_rep2  = "data/scRNAseq/10x_analysis_11534-BZ/Sample_11534-BZ-Pool01/per_sample_outs/11534-BZ-4/count/sample_filtered_feature_bc_matrix",
  day20_ko_rep1  = "data/scRNAseq/10x_analysis_11534-BZ/Sample_11534-BZ-Pool01/per_sample_outs/11534-BZ-6/count/sample_filtered_feature_bc_matrix",
  day20_ko_rep2 = "data/scRNAseq/10x_analysis_11534-BZ/Sample_11534-BZ-Pool01/per_sample_outs/11534-BZ-7/count/sample_filtered_feature_bc_matrix",
  day20_het_rep1 = "data/scRNAseq/10x_analysis_11534-BZ/Sample_11534-BZ-Pool01/per_sample_outs/11534-BZ-10/count/sample_filtered_feature_bc_matrix",
  day20_het_rep2 = "data/scRNAseq/10x_analysis_11534-BZ/Sample_11534-BZ-Pool01/per_sample_outs/11534-BZ-12/count/sample_filtered_feature_bc_matrix",
  day60_wt_rep1  = "data/scRNAseq/10x_analysis_11534-BZ/Sample_11534-BZ-Pool01/per_sample_outs/11534-BZ-14/count/sample_filtered_feature_bc_matrix",
  day60_wt_rep2  = "data/scRNAseq/10x_analysis_11534-BZ/Sample_11534-BZ-Pool01/per_sample_outs/11534-BZ-16/count/sample_filtered_feature_bc_matrix",
  day60_ko_rep1  = "data/scRNAseq/10x_analysis_11534-BZ/Sample_11534-BZ-Pool01/per_sample_outs/11534-BZ-18/count/sample_filtered_feature_bc_matrix",
  day60_ko_rep2 = "data/scRNAseq/10x_analysis_11534-BZ/Sample_11534-BZ-Pool01/per_sample_outs/11534-BZ-20/count/sample_filtered_feature_bc_matrix",
  day60_het_rep1 = "data/scRNAseq/10x_analysis_11534-BZ/Sample_11534-BZ-Pool01/per_sample_outs/11534-BZ-22/count/sample_filtered_feature_bc_matrix",
  day60_het_rep2 = "data/scRNAseq/10x_analysis_11534-BZ/Sample_11534-BZ-Pool01/per_sample_outs/11534-BZ-24/count/sample_filtered_feature_bc_matrix"
)

# Create the expression matrix from sample dirs
#   Read10X needs a *vector* instead of a *list*, so we use *unlist* to convert
geo_mat = Read10X(data.dir = unlist(sample_dirs))

# -------------------------------------------------------------------------
# Build BPCells input dir
# To use BPCells (to save some memory), you can transform 
# the expression matrix data structure above into BPCells files.
# BPCells uses these files for processing, but you typically never look at 
# their contents directly
write_matrix_dir(mat = geo_mat, dir = 'data/scRNAseq/iwaselab', overwrite = TRUE)

# -------------------------------------------------------------------------
# Cleanup
# Since we'll now be reading in from BPCells files, we will remove geo_mat
# from the environment and then prompt RStudio to run a "garbage collection"
# to free up unused memory
rm(geo_mat)
gc()

# -------------------------------------------------------------------------
# Create expression matrix and Seurat object from BPCells files
geo_mat = open_matrix_dir(dir = 'data/scRNAseq/iwaselab')
# -------------------------------------------------------------------------
# Create seurat object
geo_so = CreateSeuratObject(counts = geo_mat, min.cells = 1, min.features = 50)
geo_so

# -------------------------------------------------------------------------

head(geo_so@meta.data)

# add day, genotype, and replicate in meta.data

# Assuming 'geo_so' is your Seurat object

# Extract rownames
cell_names <- rownames(geo_so@meta.data)

# Initialize empty vectors to store extracted information
days_formatted <- c()
genotypes <- c()
replicates_formatted <- c()

# Loop through each cell name to extract information
for (name in cell_names) {
  # Split the string by '_'
  parts <- strsplit(name, "_")[[1]]
  
  # Extract 'day' and format it as "dayXX"
  # Assuming the day is always in the first part and starts with "day"
  if (length(parts) >= 1 && grepl("^day", parts[1])) {
    day_val_raw <- parts[1] # e.g., "day20"
    days_formatted <- c(days_formatted, day_val_raw)
  } else {
    days_formatted <- c(days_formatted, NA) # Or handle as per your data
  }
  
  
  # Extract 'genotype' - usually the second part
  if (length(parts) >= 2) {
    genotype_val <- parts[2]
    genotypes <- c(genotypes, genotype_val)
  } else {
    genotypes <- c(genotypes, NA)
  }
  
  
  # Extract 'replicate' and format it as "repX"
  # Assuming replicate is in the third part and might be "rep1" or just "1"
  if (length(parts) >= 3) {
    replicate_part <- parts[3] # e.g., "rep1" or "1"
    # Check if it already has "rep", if not, add it
    if (!grepl("^rep", replicate_part)) {
      replicate_formatted_val <- paste0("rep", gsub("[^0-9]", "", replicate_part))
    } else {
      replicate_formatted_val <- replicate_part
    }
    replicates_formatted <- c(replicates_formatted, replicate_formatted_val)
  } else {
    replicates_formatted <- c(replicates_formatted, NA)
  }
}

# Add the new columns to the meta.data
geo_so@meta.data$day <- days_formatted
geo_so@meta.data$genotype <- genotypes
geo_so@meta.data$replicate <- replicates_formatted

# You can then check the updated meta.data
head(geo_so@meta.data)

# Initialize an empty vector to store the new orig.ident values
new_orig_idents <- c()

# Loop through each cell name to extract the desired parts
for (name in cell_names) {
  # Split the string by '_'
  parts <- strsplit(name, "_")[[1]]
  
  # We want to reconstruct "day[day_num]_[genotype]_[replicate_num]"
  # based on the example "day20_wt_rep1_AAACAAGCAAGCTAATACTTTAGG-1"
  # This typically corresponds to the first three parts of the split string.
  
  # Check if there are at least 3 parts to avoid errors
  if (length(parts) >= 3) {
    reconstructed_id <- paste(parts[1], parts[2], parts[3], sep = "_")
    # Clean up the replicate part if it still contains non-numeric characters
    reconstructed_id <- gsub("rep([0-9]+)", "rep\\1", reconstructed_id)
    new_orig_idents <- c(new_orig_idents, reconstructed_id)
  } else {
    # Handle cases where the rowname might not have enough parts
    # For now, we'll assign NA or keep the original, depending on your preference
    # Let's assign an NA for simplicity, or you could keep the existing orig.ident value
    new_orig_idents <- c(new_orig_idents, NA) # Or geo_so@meta.data[name, "orig.ident"]
  }
}

# Assign the new values to the orig.ident column
geo_so@meta.data$day <- days
geo_so@meta.data$genotype <- genotypes
geo_so@meta.data$replicate <- replicates
geo_so@meta.data$orig.ident <- new_orig_idents

# --- NEW ADDITION: Convert columns to factors ---
geo_so@meta.data$day <- factor(geo_so@meta.data$day)
geo_so@meta.data$genotype <- factor(geo_so@meta.data$genotype)
geo_so@meta.data$replicate <- factor(geo_so@meta.data$replicate)
# --------------------------------------------------

# You can then check the updated meta.data and their types
head(geo_so@meta.data)
str(geo_so@meta.data) # Use str() to see the data types of columns

# You can then check the updated meta.data
head(geo_so@meta.data)
unique(geo_so@meta.data$day)
unique(geo_so@meta.data$genotype)
unique(geo_so@meta.data$replicate)

# Save the Seurat object
saveRDS(geo_so, file = 'data/scRNAseq/iwaselab/rdata/geo_so_unfiltered.rds')

#-----------------------------

# Load the Seurat object
geo_so = readRDS('data/scRNAseq/iwaselab/rdata/geo_so_unfiltered.rds')
head(geo_so@meta.data)

# -------------------------------------------------------------------------
cell_counts_pre_tbl = geo_so@meta.data %>% count(orig.ident, name = 'prefilter_cells')
cell_counts_pre_tbl

# -------------------------------------------------------------------------
# Review feature violin plots
# Make a violin plot of nFeature_RNA grouped by orig.ident
VlnPlot(
  object = geo_so,
  features = "nFeature_RNA",
  group.by = "orig.ident",
  layer = "count",
  pt.size = 0 # Optional: Adjust point size for individual cells, set to 0 to hide them
) +
  ggtitle("nFeature_RNA by Sample") + # Optional: Add a title
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Optional: Rotate x-axis labels for readability

ggsave(filename = 'results/scRNAseq/1_qc/1_nFeature_violin.png',
       width = 12, height = 6, units = 'in')

# Make a violin plot of nCount_RNA grouped by orig.ident
VlnPlot(
  object = geo_so,
  features = "nCount_RNA",
  group.by = "orig.ident",
  layer = "count",
  pt.size = 0 # Optional: Adjust point size for individual cells, set to 0 to hide them
) +
  ggtitle("nCount_RNA by Sample") + # Optional: Add a title
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Optional: Rotate x-axis labels for readability

ggsave(filename = 'results/scRNAseq/1_qc/1_nCount_violin.png',
       width = 12, height = 6, units = 'in')
# -------------------------------------------------------------------------
# Consider mitochondrial transcripts
# We use "MT" because this is human, depending on the organism, this might need to be changed
geo_so$percent.mt = PercentageFeatureSet(geo_so, pattern = '^MT-')

# Use summary() to quickly check the range of values
summary(geo_so$percent.mt)

# Make a violin plot of percent.mt grouped by orig.ident
VlnPlot(
  object = geo_so,
  features = "percent.mt",
  group.by = "orig.ident",
  layer = "count",
  pt.size = 0 # Optional: Adjust point size for individual cells, set to 0 to hide them
) +
  ggtitle("Percent Mitochondria by Sample") + # Optional: Add a title
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) # Optional: Rotate x-axis labels for readability

ggsave(filename = 'results/scRNAseq/1_qc/1_mt_violin.png',
       width = 12, height = 6, units = 'in')

# -------------------------------------------------------------------------
# Consider excluding suspect cells
subset(geo_so, subset = nFeature_RNA > 500 & percent.mt < 15)@meta.data %>% 
  count(orig.ident, name = 'postfilter_cells')


# -------------------------------------------------------------------------
# Examine remaining cell counts
cell_counts_post_tbl = geo_so@meta.data %>% count(orig.ident, name = 'postfilter_cells')
cell_counts_post_tbl

# -------------------------------------------------------------------------
# Show cell counts before and after filtering
cell_counts_tbl = cell_counts_pre_tbl %>% left_join(cell_counts_post_tbl, by = 'orig.ident')
cell_counts_tbl

# -------------------------------------------------------------------------
write_csv(cell_counts_tbl, file = 'results/scRNAseq/1_qc/1_cell_filtering_counts.csv')

# -------------------------------------------------------------------------
# Save the current Seurat object
saveRDS(geo_so, file = 'data/scRNAseq/iwaselab/rdata/geo_so_filtered.rds')

# =========================================================================
# Normalization
# =========================================================================

# -------------------------------------------------------------------------
# readin the so 
geo_so = readRDS('data/scRNAseq/iwaselab/rdata/geo_so_filtered.rds')

# Separate sample data into layers
geo_so[['RNA']] = split(geo_so[['RNA']], f = geo_so$orig.ident)
geo_so

# -------------------------------------------------------------------------
# Normalize the data with SCTransform
geo_so = SCTransform(geo_so)
geo_so

# -------------------------------------------------------------------------
# Save Seurat object
saveRDS(geo_so, file = 'data/scRNAseq/iwaselab/rdata/geo_so_sct_normalized.rds')

# =========================================================================
# PCA and Integration
# =========================================================================

geo_so = readRDS('data/scRNAseq/iwaselab/rdata/geo_so_sct_normalized.rds')

# Build a PCA and add it to the Seurat object
geo_so = RunPCA(geo_so, reduction.name = 'unintegrated.sct.pca')
geo_so

# -------------------------------------------------------------------------
# look at first PC alone first
heatmap_1 <- DimHeatmap(geo_so, dims=1, cells=500, balanced=TRUE, fast = FALSE, reduction = 'unintegrated.sct.pca') # doesn't include which PC(s) are included in plot
heatmap_1 + 
  labs(title="PC1: top 30 genes x top 500 cells") + 
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) 

# save to file
ggsave(filename = 'results/scRNAseq/2_pca/1_top30.png',
       width = 6, height = 4.5, units = 'in')
#  -------------------------------------------------------
# plot cell by gene heatmaps for larger set of dimensions
DimHeatmap(geo_so, dims=1:18, cells=500, balanced=TRUE, combine = FALSE, reduction = 'unintegrated.sct.pca')

# need to use png() to write to file because this isn't a ggplot
png(filename = 'results/scRNAseq/2_pca/2_qc_pca_heatmap-pc1-18.png', width = 12, height = 30, units = 'in', res = 300)
DimHeatmap(geo_so, dims=1:18, cells=500, balanced=TRUE, reduction = 'unintegrated.sct.pca')
dev.off()
# -------------------------------------------------------------------------
# Visualize how many PCs to include using an elbow plot
ElbowPlot(geo_so, ndims = 50, reduction = 'unintegrated.sct.pca')
ggsave(filename = 'results/scRNAseq/2_pca/3_qc_sct_elbow_plot.png',
       width = 8, height = 8, units = 'in')
# -------------------------------------------------------------------------
# Estimate the number of PCs to use for clustering with a function
check_pcs = function(so, reduction) {
  # quantitative check for number of PCs to include
  pct = so@reductions[[reduction]]@stdev / sum(so@reductions[[reduction]]@stdev) * 100
  cum = cumsum(pct)
  co1 = which(cum > 90 & pct < 5)[1]
  co2 = sort(which((pct[1:length(pct)-1] - pct[2:length(pct)]) > .1), decreasing = T)[1] + 1
  pcs = min(co1, co2) 
  
  return(pcs)
}

# Apply function to our data
pcs = check_pcs(geo_so, 'unintegrated.sct.pca')
pcs

pcs = 13 #rebecca used 19
  # -------------------------------------------------------------------------
str(geo_so@meta.data)
  # Visualize PCA (PC1 and PC2, unintegrated)
  # first label by sample
  DimPlot(geo_so, reduction = 'unintegrated.sct.pca', group.by = 'genotype') 

  ggsave(filename = 'results/scRNAseq/2_pca/4_qc_pca_plot_unintegrated_sct_genotype.png',
         width = 7, height = 6, units = 'in')
  
  # then label by day, which looks more useful
  DimPlot(geo_so, reduction = 'unintegrated.sct.pca', group.by = 'day')
  
  # save plot labeled by day
  ggsave(filename = 'results/scRNAseq/2_pca/4_qc_pca_plot_unintegrated_sct_day.png',
         width = 7, height = 6, units = 'in')
  # =========================================================================
  # Clustering and Projection
  # =========================================================================
  
  # Create KNN graph for PCs with `FindNeighbors()`
  geo_so = FindNeighbors(geo_so, dims = 1:pcs, reduction = 'unintegrated.sct.pca')
  
  # Then generate clusters
  geo_so = FindClusters(geo_so,
                        resolution = 0.4,
                        cluster.name = 'unintegrated.sct.rpca.clusters')
  
  # look at meta.data to see cluster labels
  head(geo_so@meta.data)

  # -------------------------------------------------------------------------
  # Create UMAP reduction
  geo_so = RunUMAP(geo_so,
                   dims = 1:pcs,
                   reduction = 'unintegrated.sct.pca',
                   reduction.name = 'umap.unintegrated.sct.rpca')
  
  
  # Note a third reduction has been added: `umap.integrated.sct.rpca`
  geo_so 
  
  # -------------------------------------------------------------------------
  # Visualize UMAP clusters: ID labels
  pre_integration_umap_plot_clusters = 
    DimPlot(geo_so, 
            group.by = 'seurat_clusters',
            label = TRUE,
            reduction = 'umap.unintegrated.sct.rpca') + 
    NoLegend()
  pre_integration_umap_plot_clusters
  
  ggsave(filename = 'results/scRNAseq/2_pca/5_umap_unintegrated_sct_clusters.png',
         plot = pre_integration_umap_plot_clusters,
         width = 6, height = 6, units = 'in')
  # -------------------------------------------------------------------------
  # Visualize UMAP clusters: clusters with labels, split by condition
  pre_integration_umap_plot_split_clusters = 
    DimPlot(geo_so,
            group.by = 'seurat_clusters',
            split.by = 'day',
            label = TRUE,
            reduction = 'umap.unintegrated.sct.rpca') + 
    NoLegend()
  pre_integration_umap_plot_split_clusters
  
  ggsave(filename = 'results/scRNAseq/2_pca/6_umap_unintegrated_sct_split_clusters.png',
         plot = pre_integration_umap_plot_split_clusters, 
         width = 14, height = 6, units = 'in')
dev.off()  
# -------------------------------------------------------------------------
# Check mitochondrial gene expression
percent_mito_plot = FeaturePlot(geo_so, features='percent.mt')

percent_mito_plot
# save to file
ggsave(filename = 'results/scRNAseq/2_pca/7_percent_umap_mito_plot.png', 
       plot = percent_mito_plot, 
       width = 6, height = 6, units = 'in')

# -------------------------------------------------------------------------
# Visualize UMAP clusters: genotype labels
pre_integration_umap_plot_genotype = 
  DimPlot(geo_so,
          group.by = 'genotype',
          label = FALSE,
          reduction = 'umap.unintegrated.sct.rpca')
pre_integration_umap_plot_genotype

ggsave(filename = 'results/scRNAseq/2_pca/8_umap_unintegrated_sct_genotype.png',
       plot = pre_integration_umap_plot_genotype,
       width = 8, height = 6, units = 'in')

# Visualize UMAP clusters: split by genotype labels
pre_integration_umap_plot_genotype_split = 
  DimPlot(geo_so,
          group.by = 'seurat_clusters',
          split.by = 'genotype',
          label = FALSE,
          reduction = 'umap.unintegrated.sct.rpca')
pre_integration_umap_plot_genotype_split

ggsave(filename = 'results/scRNAseq/2_pca/8_umap_unintegrated_sct_genotype_split.png',
       plot = pre_integration_umap_plot_genotype_split,
       width = 12, height = 6, units = 'in')

# -------------------------------------------------------------------------
# Save Seurat object
saveRDS(geo_so, file = 'data/scRNAseq/iwaselab/rdata/geo_so_sct_clustered.rds')


