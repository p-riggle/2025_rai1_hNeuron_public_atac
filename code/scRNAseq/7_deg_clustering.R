# perform clustering to group genes based on differential expression between days and genotypes.
# ref https://2-bitbio.com/2017/04/clustering-rnaseq-data-making-heatmaps.html


library(DESeq2)

# read dds r objects
dds_d20 <- readRDS("data/scRNAseq/iwaselab/dds_d20")
dds_d60 <- readRDS("data/scRNAseq/iwaselab/dds_d60")

# transform data
vsd_dds_d20 <- vst(dds_d20, blind = TRUE)  
vsd_dds_d60 <- vst(dds_d60, blind = TRUE)  

# save transformed data as data frame
vsd_dds_d20_df <- as.data.frame(assay(vsd_dds_d20))
vsd_dds_d60_df <- as.data.frame(assay(vsd_dds_d60))

head(vsd_dds_d20_df)
head(vsd_dds_d60_df)

# calculate average of two replicates
# Calculate average of two replicates for all d20 conditions (from _0 to _6)
## day20
for (i in 0:6) {
    for (cond in c("wt", "het", "ko")) {
        colname <- paste0("day20-", cond, "_", i)
        rep1 <- paste0("day20-", cond, "-rep1_", i)
        rep2 <- paste0("day20-", cond, "-rep2_", i)
        vsd_dds_d20_df[[colname]] <- rowMeans(vsd_dds_d20_df[, c(rep1, rep2)], na.rm = FALSE)
    }
}

## day60
for (i in 0:6) {
    for (cond in c("wt", "het", "ko")) {
        colname <- paste0("day60-", cond, "_", i)
        rep1 <- paste0("day60-", cond, "-rep1_", i)
        rep2 <- paste0("day60-", cond, "-rep2_", i)
        # Only calculate if both replicates exist in the data frame
        if (all(c(rep1, rep2) %in% colnames(vsd_dds_d60_df))) {
            vsd_dds_d60_df[[colname]] <- rowMeans(vsd_dds_d60_df[, c(rep1, rep2)], na.rm = FALSE)
        }
    }
}

# remove the replicate columns
vsd_dds_d20_df <- vsd_dds_d20_df[, !grepl("-rep[12]", colnames(vsd_dds_d20_df))]
vsd_dds_d60_df <- vsd_dds_d60_df[, !grepl("-rep[12]", colnames(vsd_dds_d60_df))]
head(vsd_dds_d20_df)
head(vsd_dds_d60_df)

# read restbls
## het
rstbl_c0_wt.vs.het <- read.csv("results/scRNAseq/4_pseudobulk/rstbl_c0_wt.vs.het.csv", row.names = 1)
rstbl_c2_wt.vs.het <- read.csv("results/scRNAseq/4_pseudobulk/rstbl_c2_wt.vs.het.csv", row.names = 1)
rstbl_c3_wt.vs.het <- read.csv("results/scRNAseq/4_pseudobulk/rstbl_c3_wt.vs.het.csv", row.names = 1)
rstbl_c4_wt.vs.het <- read.csv("results/scRNAseq/4_pseudobulk/rstbl_c4_wt.vs.het.csv", row.names = 1)
rstbl_c6_wt.vs.het <- read.csv("results/scRNAseq/4_pseudobulk/rstbl_c6_wt.vs.het.csv", row.names = 1)
rstbl_c1_wt.vs.het <- read.csv("results/scRNAseq/4_pseudobulk/rstbl_c1_wt.vs.het.csv", row.names = 1)
rstbl_c5_wt.vs.het <- read.csv("results/scRNAseq/4_pseudobulk/rstbl_c5_wt.vs.het.csv", row.names = 1)

## ko
rstbl_c0_wt.vs.ko <- read.csv("results/scRNAseq/4_pseudobulk/rstbl_c0_wt.vs.ko.csv", row.names = 1)
rstbl_c2_wt.vs.ko <- read.csv("results/scRNAseq/4_pseudobulk/rstbl_c2_wt.vs.ko.csv", row.names = 1)
rstbl_c3_wt.vs.ko <- read.csv("results/scRNAseq/4_pseudobulk/rstbl_c3_wt.vs.ko.csv", row.names = 1)
rstbl_c4_wt.vs.ko <- read.csv("results/scRNAseq/4_pseudobulk/rstbl_c4_wt.vs.ko.csv", row.names = 1)
rstbl_c6_wt.vs.ko <- read.csv("results/scRNAseq/4_pseudobulk/rstbl_c6_wt.vs.ko.csv", row.names = 1)
rstbl_c1_wt.vs.ko <- read.csv("results/scRNAseq/4_pseudobulk/rstbl_c1_wt.vs.ko.csv", row.names = 1)
rstbl_c5_wt.vs.ko <- read.csv("results/scRNAseq/4_pseudobulk/rstbl_c5_wt.vs.ko.csv", row.names = 1)

# for now, only analyze het data
# combine het restbls into a single data frame for easier handling
het_rstbls <- cbind.data.frame(
  c0 = rstbl_c0_wt.vs.het,
  c1 = rstbl_c1_wt.vs.het,
  c2 = rstbl_c2_wt.vs.het,
  c3 = rstbl_c3_wt.vs.het,
  c4 = rstbl_c4_wt.vs.het,
  c5 = rstbl_c5_wt.vs.het,
  c6 = rstbl_c6_wt.vs.het
)

head(het_rstbls)

# remove non-needed columns: lfcSE, stat, pvalue, baseMean
cols_to_remove <- c("lfcSE", "stat", "pvalue", "baseMean")
het_rstbls <- het_rstbls[, !grepl(paste0("(", paste(cols_to_remove, collapse="|"), ")$"), colnames(het_rstbls))]
head(het_rstbls)

## change log2FoldChange to log2fc
colnames(het_rstbls) <- gsub("log2FoldChange", "log2fc", colnames(het_rstbls))
head(het_rstbls)

# select degs with a stringent cutoff (padj < 0.01, |log2fc| > 1)
deg_c0 <- subset(het_rstbls, c0.padj < 0.01 & abs(c0.log2fc) > 1)
deg_c1 <- subset(het_rstbls, c1.padj < 0.01 & abs(c1.log2fc) > 1)
deg_c2 <- subset(het_rstbls, c2.padj < 0.01 & abs(c2.log2fc) > 1)
deg_c3 <- subset(het_rstbls, c3.padj < 0.01 & abs(c3.log2fc) > 1)
deg_c4 <- subset(het_rstbls, c4.padj < 0.01 & abs(c4.log2fc) > 1)
deg_c5 <- subset(het_rstbls, c5.padj < 0.01 & abs(c5.log2fc) > 1)
deg_c6 <- subset(het_rstbls, c6.padj < 0.01 & abs(c6.log2fc) > 1)

# combine rownames of selected degs into a list
## day20
deg_d20_het <- c(rownames(deg_c0), rownames(deg_c2),
  rownames(deg_c3), rownames(deg_c4),rownames(deg_c6)
)

head(deg_d20_het)
length(deg_d20_het)

# remove duplicates
deg_d20_het <- unique(deg_d20_het)
length(deg_d20_het)

## day60
deg_d60_het <- c(rownames(deg_c1), rownames(deg_c5))
head(deg_d60_het)
length(deg_d60_het)

deg_d60_het <- unique(deg_d60_het)
length(deg_d60_het)

# Subset vsd_dds_d20 only with het DEGs
## day20
head(vsd_dds_d20_df)
deg_vsd_d20_het <- vsd_dds_d20_df[rownames(vsd_dds_d20_df) %in% deg_d20_het, ]
nrow(deg_vsd_d20_het)
head(deg_vsd_d20_het)

## remove _1 and _5 columns
deg_vsd_d20_het <- deg_vsd_d20_het[, !grepl("_1$|_5$", colnames(deg_vsd_d20_het))]
head(deg_vsd_d20_het)

## day60
head(vsd_dds_d60_df)
deg_vsd_d60_het <- vsd_dds_d60_df[rownames(vsd_dds_d60_df) %in% deg_d60_het, ]
nrow(deg_vsd_d60_het)
head(deg_vsd_d60_het)

## remove _0, _2, _3, _4, and _6 columns
deg_vsd_d60_het <- deg_vsd_d60_het[, !grepl("_0$|_2$|_3$|_4$|_6$", colnames(deg_vsd_d60_het))]
head(deg_vsd_d60_het)

# convert data.frame to matrix
deg_vsd_d20_het <- as.matrix(deg_vsd_d20_het)
deg_vsd_d60_het <- as.matrix(deg_vsd_d60_het)

head(deg_vsd_d20_het)

# scale data to z scores
deg_d20_het_z <- t(scale(t(deg_vsd_d20_het)))
deg_d60_het_z <- t(scale(t(deg_vsd_d60_het)))

# check if na values are present
sum(is.na(deg_d20_het_z))
sum(is.na(deg_d60_het_z))


# how many cluters are good for K-means clustering?
# SSE method
## day20
head(deg_d20_het_z)
wss_d20_het <- (nrow(deg_d20_het_z)-1)*sum(apply(deg_d20_het_z, 2, var))
for (i in 1:20) wss_d20_het[i] <- sum(kmeans(deg_d20_het_z, centers=i)$withinss)
str(wss_d20_het)
plot(1:20, wss_d20_het, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")
# ~5 cluster seems appropriate
dev.off()

## day60
head(deg_d60_het_z)
wss_d60_het <- (nrow(deg_d60_het_z)-1)*sum(apply(deg_d60_het_z, 2, var))
for (i in 1:20) wss_d60_het[i] <- sum(kmeans(deg_d60_het_z, centers=i)$withinss)
str(wss_d60_het)
plot(1:20, wss_d60_het, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")
# ~5 cluster seems appropriate
dev.off()

# hierchical clustering
hp_d20_het <- hclust(as.dist(1-cor(t(deg_d20_het_z), method="pearson")), method="complete") # Cluster rows by Pearson correlation.
str(hp_d20_het)

hp_d60_het <- hclust(as.dist(1-cor(t(deg_d60_het_z), method="pearson")), method="complete") # Cluster rows by Pearson correlation.
str(hp_d60_het)

# save the cluster object
save(hp_d20_het, file = "results/scRNAseq/4_pseudobulk/hp_d20_het.RData")
save(hp_d60_het, file = "results/scRNAseq/4_pseudobulk/hp_d60_het.RData")

# read the cluster objects
load("results/scRNAseq/4_pseudobulk/hp_d20_het.RData")
load("results/scRNAseq/4_pseudobulk/hp_d60_het.RData")

# visualize the dendrograms
tree_deg_d20_het <- as.dendrogram(hp_d20_het, method = "average")
tree_deg_d60_het <- as.dendrogram(hp_d60_het, method = "average")

plot(tree_deg_d20_het, main = "Hierarchical clustering of HET-DEGs (Day 20)", leaflab = "none", ylab="Height")
plot(tree_deg_d60_het, main = "Hierarchical clustering of HET-DEGs (Day 60)", leaflab = "none", ylab="Height")
dev.off()

# cut the dendrogram to generate 5 (d20) or 4 (d60) clusters
hclust_d20_1.82 <- cutree(hp_d20_het, h=1.82)
hclust_d60_1.6 <- cutree(hp_d60_het, h=1.6)

unique(hclust_d20_1.82)
unique(hclust_d60_1.6)

# plot tree again and add the cluster information
library('dendextend')
library(gridExtra)
library(ggplot2)

# make line plot for each cluster
# install.packages("tidyverse")
library(tidyverse)
str(hclust_d20_1.82)
head(hclust_d20_1.82)
class(hclust_d20_1.82)

# convert the Z-scores matric to data.frame
deg_d20_het_z_df <- as.data.frame(deg_d20_het_z)
deg_d60_het_z_df <- as.data.frame(deg_d60_het_z)
head(deg_d20_het_z_df)
deg_d20_het_z_df$gene <- row.names(deg_d20_het_z_df) 
deg_d60_het_z_df$gene <- row.names(deg_d60_het_z_df)
head(deg_d20_het_z_df)

# assign cluster labels to each gene using cutree
hp_d20_het_clusters <- cutree(hp_d20_het, h = 1.82)
hp_d60_het_clusters <- cutree(hp_d60_het, h = 1.6)
head(hp_d20_het_clusters)

# convert cluster assignments to data frames for joining
hp_d20_het_df <- data.frame(gene = names(hp_d20_het_clusters), cluster = hp_d20_het_clusters)
hp_d60_het_df <- data.frame(gene = names(hp_d60_het_clusters), cluster = hp_d60_het_clusters)

# add cluster names to each gene
h_d20_het <- deg_d20_het_z_df %>% left_join(hp_d20_het_df, by = "gene")
h_d60_het <- deg_d60_het_z_df %>% left_join(hp_d60_het_df, by = "gene")

# save the cluster information
write.table(h_d20_het, file = "results/scRNAseq/4_pseudobulk/h_d20_het.txt", sep = "\t", row.names = FALSE)
write.table(h_d60_het, file = "results/scRNAseq/4_pseudobulk/h_d60_het.txt", sep = "\t", row.names = FALSE)

# read the cluster information file
h_d20_het <- read.table("results/scRNAseq/4_pseudobulk/h_d20_het.txt", header = TRUE, sep = "\t")
h_d60_het <- read.table("results/scRNAseq/4_pseudobulk/h_d60_het.txt", header = TRUE, sep = "\t")

head(h_d20_het)
head(h_d60_het)

# make sample sheet
head(h_d20_het)
sample_d20 <- colnames(h_d20_het)[1:15]
celltype <- c(rep(c(0, 2, 3, 4, 6), each = 3))
genotype <- c(rep(c("wt", "het", "ko"), 5))
traits_d20 <- as.data.frame(cbind(sample_d20, celltype, genotype))
traits_d20

head(h_d60_het)
sample_d60 <- colnames(h_d60_het)[1:6]
celltype <- c(rep(c(1, 5), each = 3))
genotype <- c(rep(c("wt", "het", "ko"), 2))
traits_d60 <- as.data.frame(cbind(sample_d60, celltype, genotype))
traits_d60

# make cluster matrix with long format for heatmap
mh_d20 <- h_d20_het  %>%
    # convert to long format
    pivot_longer(cols = colnames(h_d20_het)[1:15], names_to = "sample_d20", values_to = "z") %>%
    # join with sample info table
    full_join(traits_d20, by = "sample_d20") %>%
    # for each gene, genotype, and day
    group_by(gene, celltype, genotype)

head(mh_d20)

mh_d60 <- h_d60_het  %>%
    # convert to long format
    pivot_longer(cols = colnames(h_d60_het)[1:6], names_to = "sample_d60", values_to = "z") %>%
    # join with sample info table
    full_join(traits_d60, by = "sample_d60") %>%
    # for each gene, genotype, and day
    group_by(gene, celltype, genotype)

head(mh_d60)

# write cluster dataframe for pathway analysis
write.table(mh_d20, file = "results/scRNAseq/4_pseudobulk/mh_d20.txt", sep = "\t")
write.table(mh_d60, file = "results/scRNAseq/4_pseudobulk/mh_d60.txt", sep = "\t")

# Plot gene expression of each cluster ####
color1 <- c("#757582", "#bc9f0dff", "#0062ffff")

  # Convert the genotype column to factor with the specified order
  genotype_order <- c("wt", "het", "ko")
  mh_d20$genotype <- factor(mh_d20$genotype, levels = genotype_order)

  # line plot
lp_h_d20 <- mh_d20 %>% 
    ggplot(aes(celltype, z, group = genotype)) +
    geom_line(stat = "summary", fun = "median", size = 1.0, aes(col = genotype), alpha = 0.3) +
    geom_point(stat = "summary", fun = "median", size = 4, aes(col = genotype), shape = 16) +
    scale_color_manual(values = color1) +
    scale_x_discrete(limits = c("0", "2", "3", "4", "6")) +
    xlab("Cell type") +
    ylab("Expression") +
    facet_grid(cols = vars(cluster)) +
    theme(
        panel.background = element_rect(fill = NA, color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))

lp_h_d20

# day60 lineplot

  # Convert the genotype column to factor with the specified order
  genotype_order <- c("wt", "het", "ko")
  mh_d60$genotype <- factor(mh_d60$genotype, levels = genotype_order)

lp_h_d60 <- mh_d60 %>% 
    ggplot(aes(celltype, z, group = genotype)) +
    geom_line(stat = "summary", fun = "median", size = 1.0, aes(col = genotype), alpha = 0.3) +
    geom_point(stat = "summary", fun = "median", size = 4, aes(col = genotype), shape = 16) +
    scale_color_manual(values = color1) +
    scale_x_discrete(limits = c("1", "5")) +
    xlab("Cell type") +
    ylab("Expression") +
    facet_grid(cols = vars(cluster)) +
    theme(
        panel.background = element_rect(fill = NA, color = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))

lp_h_d60

# save lineplot as pdf d20 and d60
ggsave("results/scRNAseq/4_pseudobulk/lp_h_d20.pdf", lp_h_d20, width = 8, height = 6)
ggsave("results/scRNAseq/4_pseudobulk/lp_h_d60.pdf", lp_h_d60, width = 8, height = 6)

# know the number of genes in each cluster
table(h_d20_het$cluster)
table(h_d60_het$cluster)

# heatmap
## define color for heatmap
library("colorspace")
hcl_palettes(plot = TRUE)
br <- diverging_hcl(10, palette = "Blue-Red3")

library(ComplexHeatmap)
library(gridtext)
library(dendextend)

str(hp_d20_het)
head(deg_d20_het_z)
head(hp_d20_het_clusters)

hmap_day20 <- Heatmap(
    deg_d20_het_z,
    name = "Z score",
    row_names_gp = gpar(fontsize = 6.5),
    row_split = factor(hp_d20_het_clusters),
    col = br,
    show_row_names = FALSE,
    cluster_columns = FALSE,
    # You can adjust the column_split if you want to split columns by celltype/genotype
    # column_split = rep(c("0", "2", "3", "4", "6"), each = 3),
    column_names_side = "top",
    column_names_rot = 0,
    width = unit(15, "cm")
)

hmap_day20 

hmap_day60 <- Heatmap(
    deg_d60_het_z,
    name = "Z score",
    row_names_gp = gpar(fontsize = 6.5),
    row_split = factor(hp_d60_het_clusters),
    col = br,
    show_row_names = FALSE,
    cluster_columns = FALSE,
    # You can adjust the column_split if you want to split columns by celltype/genotype
    # column_split = rep(c("0", "2", "3", "4", "6"), each = 3),
    column_names_side = "top",
    column_names_rot = 0,
    width = unit(15, "cm")
)

hmap_day60
