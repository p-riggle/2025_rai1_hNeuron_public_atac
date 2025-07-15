
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



