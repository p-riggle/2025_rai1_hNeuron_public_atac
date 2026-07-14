# the dataset files contained in objects "combined" and "combined_FPKM" were created in code/Fig7_huang/0_combinefiles_ENCODE_WT.R script.
# this code reads in those files, makes a deseq object, does a PCA plot, and extracts the top 300 PC1 genes for plotting in a lineplot.
# this is the beginning of the pipeline following the established method from iNeuron data analysis, now applied to ENCODE and Huang datasets.

# load libraries
library(DESeq2)
library(ggplot2)
library(tidyr)
library(readr)

# read in the combined count data
combined <- read.table("data/huang/combined_ENCODE_data.tsv")

# make the experimental design file
design <- data.frame(row.names = colnames(combined), conditions = c( "count_p10", "count_p14", "count_p25", "count_p36", 
                                                                     "count_p60", "count_p10_rep2", "count_p14_rep2", "count_p25_rep2", "count_p36_rep2", "count_p60_rep2"))

# create a new column in design to designate the timepoint
design$timepoint <- gsub("^count_p([0-9]+).*", "\\1", design$conditions)
design 

# make conditions a factor (not necessary, but then deseq wont give you the warning message)
design$conditions <- as.factor(design$conditions)

# create deseq object
dds <- DESeqDataSetFromMatrix(combined, colData=design, design = ~ timepoint) 


# variance stabilizing transformation
vsd <- vst(dds, blind = TRUE)
head(assay(vsd), 3)

# get PCA data
PCA <- plotPCA(vsd, intgroup = c("conditions", "timepoint")) # using top 500 features by variance
PCA_data <- plotPCA(vsd, intgroup = c("conditions", "timepoint"), returnData = TRUE) 
percentVar <- round(100 * attr(PCA_data, "percentVar"))

# plot PCA
# create the color mapping
color_mapping <- c("10" = "dark grey", "14" = "blue", "25" = "purple", "36" = "red", "60" = "green")

pca_plot <- ggplot(PCA_data, aes(x = PC1, y = PC2, color = timepoint)) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  geom_point(size = 4, shape = 1, stroke = 1.5) +
  scale_color_manual(values = color_mapping) +
  ggtitle("PCA Plot of RNA-seq Data (ENCODE cortex WT)") +
  theme(text = element_text(family = "Arial", size = 12),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8)) +
  theme_classic()
pca_plot

# save plot
ggsave(filename = sprintf("result/Fig7_huang/ENCODE_WT_PCA_Fig7A.pdf", Sys.Date()), plot = pca_plot, width = 6, height = 4.5, units = "in")


# now extract the PC1 genes, and plot as lineplot
# extract gene loadings to PC1

# 1. Get variance-stabilized assay matrix
vsd_mat <- assay(vsd) 

# 2. Filter out genes with zero variance across samples
gene_vars <- apply(vsd_mat, 1, var)
vsd_mat_filtered <- vsd_mat[gene_vars > 0, ]

# 3. Run PCA (samples as rows, genes as columns) to give the gene loadings
pca_res <- prcomp(t(vsd_mat_filtered), scale. = TRUE)

# 4. Extract loadings for PC1
loadings_pc1 <- pca_res$rotation[, 1]

# plot PC1 eigenvalue loadings as a histogram
hist(loadings_pc1, breaks = 50, main = "Histogram of PC1 Loadings", xlab = "PC1 Loadings")

# 5. Get 300 genes with highest absolute loading for PC1
top300_pc1_genes <- names(sort(abs(loadings_pc1), decreasing = TRUE)[1:300])

# plot top300 pc1 gene eigenvalues as a histogram
hist(loadings_pc1[top300_pc1_genes], breaks = 50, main = "Histogram of Top 300 PC1 Loadings", xlab = "Top 300 PC1 Loadings")

head(top300_pc1_genes)

# save the top 300 PC1 genes to a file
write.table(
  top300_pc1_genes, 
  file = "result/huang/january_huang/encode_top300_PC1_genes.txt", 
  sep = "\t", quote = FALSE, row.names = FALSE
)

# now, plot the top300 PC1 genes from the ENCODE data in a lineplot to visualize. Plot the FPKM values, zscore normalized. 
# read in combined_FPKM
combined_FPKM <- read.delim("data/huang/combined_FPKM_encode_2reps.tsv", sep = "\t", header = TRUE)
# subset the encode FPKM data for the PC1 genes 

filtered_combined <- combined_FPKM[combined_FPKM$ensembl_id %in% top300_pc1_genes, ] 
# zscore normalize
filtered_combined_z <- filtered_combined

filtered_combined_z[, c("FPKM_p10", "FPKM_p14", "FPKM_p25", "FPKM_p36", "FPKM_p60")] <- t(apply(
  filtered_combined_z[, c("FPKM_p10", "FPKM_p14", "FPKM_p25", "FPKM_p36", "FPKM_p60")], 1,
  function(x) (x - mean(x)) / sd(x)
))

# save filtered_combined_z for use in clustering_lineplot_encode.R script
write_tsv(filtered_combined_z, file = "result/huang/january_huang/encode_pc1_top300_genes_FPKM_zscore.tsv")

# plot 
# ### RESHAPE TO LONG FORMAT FOR PLOTTING
long_filtered <- pivot_longer(
  filtered_combined_z,
  cols = starts_with("FPKM"),
  names_to = "Timepoint",
  values_to = "Zscore"
)
long_filtered$Timepoint <- factor(
  long_filtered$Timepoint,
  levels = c("FPKM_p10", "FPKM_p14", "FPKM_p25", "FPKM_p36", "FPKM_p60")
)

### PLOT LINEPLOT to assess expression pattern across time

lineplot_pc1 <- ggplot(long_filtered, aes(x = Timepoint, y = Zscore, group = ensembl_id)) +
  geom_line(color = "gray") +
  geom_point(size=0.5, color = "black") +
  theme_bw() +
  labs(
    x = "Timepoint",
    y = "Z-score",
    title = "PC1 genes, WT, Z-score normalized)"
  )
lineplot_pc1

ggsave(filename = sprintf("result/huang/january_huang/encode_pc1300genes_lineplot.pdf", Sys.Date()), plot = lineplot_pc1, width = 6, height = 4.5, units = "in")


#--------------------------------------------
