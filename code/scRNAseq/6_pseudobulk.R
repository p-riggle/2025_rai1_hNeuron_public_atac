# pseudobulk method

# specify 

# Load libraries
# Enter commands in R (or R studio, if installed)
library(Seurat)
library(BPCells)
library(tidyverse)

# load the R.data
geo_so = readRDS('data/scRNAseq/iwaselab/rdata/geo_so_sct_unintegrated_with_markers_0.1.rds')

# define the output directory
dirout <- ("/nfs/turbo/umms-siwase/SHIGEKI2/rai1_hesc/results/scRNAseq/4_pseudobulk/")

# edit metadata

# clustMetaDat <- geo_so@meta.data
# Sample <- paste0(clustMetaDat$orig.ident, ".", clustMetaDat$genotype)
# Clust <- clustMetaDat$seurat_clusters
# aDF <- data.frame(Sample,Clust)
# clustMetaDat <- cbind(clustMetaDat, aDF)
# geo_so@meta.data <- clustMetaDat
# head(clustMetaDat)

head(geo_so@meta.data)
meta_data <- geo_so@meta.data

# aggregate count data by cluster
DefaultAssay(geo_so) <- "RNA"
cts <- AggregateExpression(geo_so, assays = "RNA", slot = "counts", return.seurat = F, group.by =c("orig.ident", "seurat_clusters"))
cts <- cts$RNA
colnames(cts)

# split two time points

cts_d20 <- cts[, c(1:42)]

colnames(cts_d20)

cts_d60 <- cts[, c(43:73)]

colnames(cts_d60)

# # make DEseq metadata
# sample <- sapply(strsplit(colnames(cts), "\\."), "[[",1)
# Gtype <- sapply(strsplit(colnames(cts), "\\."), "[[",2)
# Clust <- paste0("Clust_",sapply(strsplit( Gtype , "_"), "[[",2))
# samplesheet <- data.frame(sample,Gtype,Clust,  Group=colnames(cts), row.names=colnames(cts) )
# head(samplesheet)
# samplesheet

sample_d20 <- sapply(strsplit(colnames(cts_d20), '_'),'[[',1)
cluster_d20 <- sapply(strsplit(colnames(cts_d20), '_'),'[[',2)
sample_d20
cluster_d20
genotype_d20 <- sapply(strsplit(sample_d20, '-'), '[[',2)
genotype_cluster_d20 <- paste(genotype_d20, cluster_d20, sep='_')
genotype_cluster_d20

samplesheet_d20 <- data.frame(sample_d20, genotype_cluster_d20, cluster_d20, genotype_d20, row.names = colnames(cts_d20))
samplesheet_d20

sample_d60 <- sapply(strsplit(colnames(cts_d60), '_'),'[[',1)
cluster_d60 <- sapply(strsplit(colnames(cts_d60), '_'),'[[',2)
sample_d60
cluster_d60
genotype_d60 <- sapply(strsplit(sample_d60, '-'), '[[',2)
genotype_cluster_d60 <- paste(genotype_d60, cluster_d60, sep='_')
genotype_cluster_d60

samplesheet_d60 <- data.frame(sample_d60, genotype_cluster_d60, cluster_d60, genotype_d60, row.names = colnames(cts_d60))
samplesheet_d60


library(DESeq2)
# calculation
dds_d20 <- DESeqDataSetFromMatrix(countData = cts_d20, colData = samplesheet_d20, design = ~ genotype_cluster_d20)
dds_d20 <- DESeq(dds_d20)

dds_d60 <- DESeqDataSetFromMatrix(countData = cts_d60, colData = samplesheet_d60, design = ~ genotype_cluster_d60)
dds_d60 <- DESeq(dds_d60)

# save the dds r object
saveRDS(dds_d20, "data/scRNAseq/iwaselab/dds_d20")
saveRDS(dds_d60, "data/scRNAseq/iwaselab/dds_d60")

# read dds r objects
dds_d20 <- readRDS("data/scRNAseq/iwaselab/dds_d20")
dds_d60 <- readRDS("data/scRNAseq/iwaselab/dds_d60")

# PCAplots
## transform data
vsd_dds_d20 <- vst(dds_d20, blind = TRUE)  
vsd_dds_d60 <- vst(dds_d60, blind = TRUE)  

head(assay(vsd_dds_d20), 3)
head(assay(vsd_dds_d60), 3)

#install.packages("ggplot2")
library("ggplot2")
str(vsd_dds_d20)

PCA_d20 <- plotPCA(vsd_dds_d20, intgroup = c("genotype_d20", "cluster_d20"))
PCA_data_d20 <- plotPCA(vsd_dds_d20, intgroup = c("genotype_d20", "cluster_d20"), returnData = TRUE)
PCA_d20

# create a shape and color mapping for the genotype and day
shape_mapping <- c("wt" = 1, "het" = 2, "ko" = 4)
color_mapping_d20 <- c("0" = "red", "2" = "green", "3" = "darkgreen", "4" = "blue")

# create a new plot with the modified aesthetics
pca_d20 <- ggplot(PCA_data_d20, aes(x = PC1, y = PC2, shape = genotype_d20, color = cluster_d20)) +
  geom_point(size = 5, stroke = 2) +
  scale_shape_manual(values = shape_mapping) +
  scale_color_manual(values = color_mapping_d20) +
  ggtitle("PCA Plot: Day 20") +
  theme(text = element_text(family = "Helvetica", size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8)) +
  theme_classic()

pca_d20

# day 60 pca
PCA_d60 <- plotPCA(vsd_dds_d60, intgroup = c("genotype_d60", "cluster_d60"))
PCA_data_d60 <- plotPCA(vsd_dds_d60, intgroup = c("genotype_d60", "cluster_d60"), returnData = TRUE)
PCA_d60

# Filter the PCA data to include only clusters 1 and 5
PCA_data_d60_filtered <- PCA_data_d60[PCA_data_d60$cluster_d60 %in% c("1", "5"), ]

# create a shape and color mapping for the genotype and day
shape_mapping <- c("wt" = 1, "het" = 2, "ko" = 4)
color_mapping_d60 <- c("1" = "gold", "5" = "purple")

# create a new plot with the modified aesthetics
pca_d60 <- ggplot(PCA_data_d60_filtered, aes(x = PC1, y = PC2, shape = genotype_d60, color = cluster_d60)) +
  geom_point(size = 5, stroke = 2) +
  scale_shape_manual(values = shape_mapping) +
  scale_color_manual(values = color_mapping_d60) +
  ggtitle("PCA Plot: Day 60") +
  theme(text = element_text(family = "Helvetica", size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8)) +
  theme_classic()

pca_d60

ggsave(filename = "results/scRNAseq/4_pseudobulk/pca2_d20.pdf", plot = pca_d20, width = 5, height = 4, units = "in")
ggsave(filename = "results/scRNAseq/4_pseudobulk/pca2_d60.pdf", plot = pca_d60, width = 5, height = 4, units = "in")

# make het restables
rstbl_c0_wt.vs.het <- results(dds_d20, contrast = c("genotype_cluster_d20", "het_0", "wt_0"), alpha = 0.01)
rstbl_c2_wt.vs.het <- results(dds_d20, contrast = c("genotype_cluster_d20", "het_2", "wt_2"), alpha = 0.01)
rstbl_c3_wt.vs.het <- results(dds_d20, contrast = c("genotype_cluster_d20", "het_3", "wt_3"), alpha = 0.01)
rstbl_c4_wt.vs.het <- results(dds_d20, contrast = c("genotype_cluster_d20", "het_4", "wt_4"), alpha = 0.01)
rstbl_c6_wt.vs.het <- results(dds_d20, contrast = c("genotype_cluster_d20", "het_6", "wt_6"), alpha = 0.01)

rstbl_c1_wt.vs.het <- results(dds_d60, contrast = c("genotype_cluster_d60", "het_1", "wt_1"), alpha = 0.01)
rstbl_c5_wt.vs.het <- results(dds_d60, contrast = c("genotype_cluster_d60", "het_5", "wt_5"), alpha = 0.01)

summary(rstbl_c0_wt.vs.het)
summary(rstbl_c2_wt.vs.het)
summary(rstbl_c3_wt.vs.het)
summary(rstbl_c4_wt.vs.het)
summary(rstbl_c6_wt.vs.het)

summary(rstbl_c1_wt.vs.het)
summary(rstbl_c5_wt.vs.het)

# make ko restables
rstbl_c0_wt.vs.ko <- results(dds_d20, contrast = c("genotype_cluster_d20", "ko_0", "wt_0"), alpha = 0.01)
rstbl_c2_wt.vs.ko <- results(dds_d20, contrast = c("genotype_cluster_d20", "ko_2", "wt_2"), alpha = 0.01)
rstbl_c3_wt.vs.ko <- results(dds_d20, contrast = c("genotype_cluster_d20", "ko_3", "wt_3"), alpha = 0.01)
rstbl_c4_wt.vs.ko <- results(dds_d20, contrast = c("genotype_cluster_d20", "ko_4", "wt_4"), alpha = 0.01)
rstbl_c6_wt.vs.ko <- results(dds_d20, contrast = c("genotype_cluster_d20", "ko_6", "wt_6"), alpha = 0.01)

rstbl_c1_wt.vs.ko <- results(dds_d60, contrast = c("genotype_cluster_d60", "ko_1", "wt_1"), alpha = 0.01)
rstbl_c5_wt.vs.ko <- results(dds_d60, contrast = c("genotype_cluster_d60", "ko_5", "wt_5"), alpha = 0.01)

summary(rstbl_c0_wt.vs.ko)
summary(rstbl_c2_wt.vs.ko)
summary(rstbl_c3_wt.vs.ko)
summary(rstbl_c4_wt.vs.ko)
summary(rstbl_c6_wt.vs.ko)

summary(rstbl_c1_wt.vs.ko)
summary(rstbl_c5_wt.vs.ko)

# write all rstbles to the results directory
write.csv(as.data.frame(rstbl_c0_wt.vs.het), file = "results/scRNAseq/4_pseudobulk/rstbl_c0_wt.vs.het.csv")
write.csv(as.data.frame(rstbl_c2_wt.vs.het), file = "results/scRNAseq/4_pseudobulk/rstbl_c2_wt.vs.het.csv")
write.csv(as.data.frame(rstbl_c3_wt.vs.het), file = "results/scRNAseq/4_pseudobulk/rstbl_c3_wt.vs.het.csv")
write.csv(as.data.frame(rstbl_c4_wt.vs.het), file = "results/scRNAseq/4_pseudobulk/rstbl_c4_wt.vs.het.csv")
write.csv(as.data.frame(rstbl_c6_wt.vs.het), file = "results/scRNAseq/4_pseudobulk/rstbl_c6_wt.vs.het.csv")
write.csv(as.data.frame(rstbl_c1_wt.vs.het), file = "results/scRNAseq/4_pseudobulk/rstbl_c1_wt.vs.het.csv")
write.csv(as.data.frame(rstbl_c5_wt.vs.het), file = "results/scRNAseq/4_pseudobulk/rstbl_c5_wt.vs.het.csv")

write.csv(as.data.frame(rstbl_c0_wt.vs.ko), file = "results/scRNAseq/4_pseudobulk/rstbl_c0_wt.vs.ko.csv")
write.csv(as.data.frame(rstbl_c2_wt.vs.ko), file = "results/scRNAseq/4_pseudobulk/rstbl_c2_wt.vs.ko.csv")
write.csv(as.data.frame(rstbl_c3_wt.vs.ko), file = "results/scRNAseq/4_pseudobulk/rstbl_c3_wt.vs.ko.csv")
write.csv(as.data.frame(rstbl_c4_wt.vs.ko), file = "results/scRNAseq/4_pseudobulk/rstbl_c4_wt.vs.ko.csv")
write.csv(as.data.frame(rstbl_c6_wt.vs.ko), file = "results/scRNAseq/4_pseudobulk/rstbl_c6_wt.vs.ko.csv")
write.csv(as.data.frame(rstbl_c1_wt.vs.ko), file = "results/scRNAseq/4_pseudobulk/rstbl_c1_wt.vs.ko.csv")
write.csv(as.data.frame(rstbl_c5_wt.vs.ko), file = "results/scRNAseq/4_pseudobulk/rstbl_c5_wt.vs.ko.csv")

# read restbls
rstbl_c0_wt.vs.het <- read.csv("results/scRNAseq/4_pseudobulk/rstbl_c0_wt.vs.het.csv", row.names = 1)
rstbl_c2_wt.vs.het <- read.csv("results/scRNAseq/4_pseudobulk/rstbl_c2_wt.vs.het.csv", row.names = 1)
rstbl_c3_wt.vs.het <- read.csv("results/scRNAseq/4_pseudobulk/rstbl_c3_wt.vs.het.csv", row.names = 1)
rstbl_c4_wt.vs.het <- read.csv("results/scRNAseq/4_pseudobulk/rstbl_c4_wt.vs.het.csv", row.names = 1)
rstbl_c6_wt.vs.het <- read.csv("results/scRNAseq/4_pseudobulk/rstbl_c6_wt.vs.het.csv", row.names = 1)
rstbl_c1_wt.vs.het <- read.csv("results/scRNAseq/4_pseudobulk/rstbl_c1_wt.vs.het.csv", row.names = 1)
rstbl_c5_wt.vs.het <- read.csv("results/scRNAseq/4_pseudobulk/rstbl_c5_wt.vs.het.csv", row.names = 1)

rstbl_c0_wt.vs.ko <- read.csv("results/scRNAseq/4_pseudobulk/rstbl_c0_wt.vs.ko.csv", row.names = 1)
rstbl_c1_wt.vs.ko <- read.csv("results/scRNAseq/4_pseudobulk/rstbl_c1_wt.vs.ko.csv", row.names = 1)
rstbl_c2_wt.vs.ko <- read.csv("results/scRNAseq/4_pseudobulk/rstbl_c2_wt.vs.ko.csv", row.names = 1)
rstbl_c3_wt.vs.ko <- read.csv("results/scRNAseq/4_pseudobulk/rstbl_c3_wt.vs.ko.csv", row.names = 1)
rstbl_c4_wt.vs.ko <- read.csv("results/scRNAseq/4_pseudobulk/rstbl_c4_wt.vs.ko.csv", row.names = 1)
rstbl_c5_wt.vs.ko <- read.csv("results/scRNAseq/4_pseudobulk/rstbl_c5_wt.vs.ko.csv", row.names = 1)
rstbl_c6_wt.vs.ko <- read.csv("results/scRNAseq/4_pseudobulk/rstbl_c6_wt.vs.ko.csv", row.names = 1)


# make deg number bar plot
# make a list of your results tables
res_list_het <- list(
  c0 = rstbl_c0_wt.vs.het,
  c1 = rstbl_c1_wt.vs.het,
  c2 = rstbl_c2_wt.vs.het,
  c3 = rstbl_c3_wt.vs.het,
  c4 = rstbl_c4_wt.vs.het,
  c5 = rstbl_c5_wt.vs.het,
  c6 = rstbl_c6_wt.vs.het
)


# Extract and count up- and down-regulated genes for each cluster
counts_df_het <- data.frame(
  cluster = character(),
  direction = character(),
  count = integer()
)

for (cluster_name in names(res_list_het)) {
  res <- res_list_het[[cluster_name]]

  # Remove rows with NA p-values
  res <- na.omit(as.data.frame(res))

  # Count up-regulated genes (padj < 0.05 and log2FoldChange > 0)
  up_genes <- sum(res$padj < 0.01 & res$log2FoldChange > 1)

  # Count down-regulated genes (padj < 0.05 and log2FoldChange < 0)
  down_genes <- sum(res$padj < 0.01 & res$log2FoldChange < -1)

  # Append to the dataframe
  counts_df_het <- rbind(counts_df_het,
                     data.frame(cluster = cluster_name, direction = "up", count = up_genes),
                     data.frame(cluster = cluster_name, direction = "down", count = down_genes)
  )
}

counts_df_het

# Load ggplot2 for plotting

library(ggplot2)

# Create the histogram (bar chart)
## het
degN_barplot_het <- ggplot(counts_df_het, aes(x = cluster, y = count, fill = direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Number of Up- and Down-regulated Genes per Cluster (HET vs. WT)",
    x = "Cluster",
    y = "Number of Genes",
    fill = "Regulation Direction"
  ) +
  scale_fill_manual(values = c("up" = "gold", "down" = "blue")) +
  geom_text(aes(label = count), position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  theme_classic()

degN_barplot_het

dev.off()

ggsave(filename = paste0(DirOutQC, "degN_bar_het.pdf"), plot = degN_barplot_het, width = 6, height = 4, units = "in")

## homo
# make deg number bar plot
# make a list of your results tables
res_list_ko <- list(
  c0 = rstbl_c0_wt.vs.ko,
  c1 = rstbl_c1_wt.vs.ko,
  c2 = rstbl_c2_wt.vs.ko,
  c3 = rstbl_c3_wt.vs.ko,
  c4 = rstbl_c4_wt.vs.ko,
  c5 = rstbl_c5_wt.vs.ko,
  c6 = rstbl_c6_wt.vs.ko
)


# Extract and count up- and down-regulated genes for each cluster
counts_df_ko <- data.frame(
  cluster = character(),
  direction = character(),
  count = integer()
)

for (cluster_name in names(res_list_het)) {
  res <- res_list_ko[[cluster_name]]
  
  # Remove rows with NA p-values
  res <- na.omit(as.data.frame(res))
  
  # Count up-regulated genes (padj < 0.05 and log2FoldChange > 0)
  up_genes <- sum(res$padj < 0.01 & res$log2FoldChange > 1)
  
  # Count down-regulated genes (padj < 0.05 and log2FoldChange < 0)
  down_genes <- sum(res$padj < 0.01 & res$log2FoldChange < -1)
  
  # Append to the dataframe
  counts_df_ko <- rbind(counts_df_ko,
                         data.frame(cluster = cluster_name, direction = "up", count = up_genes),
                         data.frame(cluster = cluster_name, direction = "down", count = down_genes)
  )
}

counts_df_ko

# Create the histogram (bar chart)

degN_barplot_ko <- ggplot(counts_df_ko, aes(x = cluster, y = count, fill = direction)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(
    title = "Number of Up- and Down-regulated Genes per Cluster (KO vs. WT)",
    x = "Cluster",
    y = "Number of Genes",
    fill = "Regulation Direction"
  ) +
  scale_fill_manual(values = c("up" = "gold", "down" = "blue")) +
  geom_text(aes(label = count), position = position_dodge(width = 0.9), vjust = -0.5, size = 3) +
  theme_classic()

degN_barplot_ko

ggsave(filename = "results/scRNAseq/4_pseudobulk/degN_bar_ko.pdf", plot = degN_barplot_ko, width = 6, height = 4, units = "in")



