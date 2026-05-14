# Fig 7, paper
# create wildtype and rai1 ko PCA plot 
# for Kasper/Iwase bulk RNAseq data to generate PCA plot

# load all libraries
library(tidyverse)
library(DESeq2)
library(dplyr)
library(clusterProfiler)
library(org.Mm.eg.db)
library(tibble)

#import gene count file and metadata
count <- read.delim("data/kasper/count/combined_counts_with_symbols.txt", header=TRUE, stringsAsFactors=FALSE)
# preparing the count file by remove duplicates, removing some columns and convert to rounded integer
count$gene_symbol <- NULL
count$gene_id <- NULL 
count_matrix <- as.data.frame(count[!duplicated(count$ensembl_gene_id), ])

rownames(count_matrix) <- NULL

count_df <- count_matrix %>%
  column_to_rownames("ensembl_gene_id")

count_df[] <- lapply(count_df, as.numeric)
count_df <- round(count_df)
count_df <- as.matrix(count_df)
storage.mode(count_df) <- "integer"

#load metadata
metadata <- read.delim("data/kasper/count/count_metadata.txt", header=TRUE, stringsAsFactors=FALSE)
metadata <- metadata %>% column_to_rownames("Sample")

# Ensure samples in metadata and counts are in the same order
all(colnames(count_df) %in% rownames(metadata)) #TRUE
head(count_df)
view(count_df)

# make the Condition and Timepoint columns as factors
metadata$Condition <- as.factor(metadata$Condition)
metadata$Timepoint <- as.factor(metadata$Timepoint)

# remove the KO data
metadata <- metadata[metadata$Condition %in% c("WT"), ]
count_df <- count_df[, colnames(count_df) %in% rownames(metadata)]
# check dimensions
dim(count_df) #  57010   12
dim(metadata) # 12  2

# make deseq2 dataset and normalize
dds <- DESeqDataSetFromMatrix(
  countData = count_df,
  colData = metadata,
  design = ~ Timepoint
)

# Variance Stabilizing Transformation (VST)
vsd <- vst(dds, blind=FALSE)

# Get PCA
pcaData <- plotPCA(vsd, intgroup=c("Condition","Timepoint"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
# save PCAData

# Plot
# make plot points shape and color to match publication figs

pca_plot <- ggplot(pcaData, aes(PC1, PC2, color=Timepoint, shape=Condition)) +
  geom_point(size=6, stroke=0.75) +
  scale_shape_manual(values = c("WT" = 1)) +
  scale_color_manual(values = c("0" = "grey", "8" = "#377EB8", "12" = "purple", "23" = "red")) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  labs(shape="Genotype", color="Day") +
  theme_minimal()   

pca_plot
# save plot
ggsave("result/kasper/january_kasper/PCA_WT_kasper.pdf", plot = pca_plot, width = 6, height = 5)


# end of script