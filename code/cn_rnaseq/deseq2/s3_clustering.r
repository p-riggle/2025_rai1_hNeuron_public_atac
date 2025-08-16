# perform clustering to group genes based on differential expression between days and genotypes.

# ref https://2-bitbio.com/2017/04/clustering-rnaseq-data-making-heatmaps.html

# import tpm file

tpm <- read.delim("data/cn_rnaseq/blueberry_output/combine_counts/gene_TPM.txt", header = TRUE, sep = "\t")
head(tpm)
colnames(tpm) <- gsub("Rai1.|.TPM|d", "", colnames(tpm))
colnames(tpm)[1] <- "gene_id"
# remove het.d34.rep1 and het.d60.rep4
tpm[, c("het_34_rep1", "het_60_rep4")] <- NULL
colnames(tpm)

#calcurate tpm average between replicates
tpm$het_0 <- rowMeans(tpm[2:5])
tpm$het_20 <- rowMeans(tpm[6:9])
tpm$het_34 <- rowMeans(tpm[10:12])
tpm$het_60 <- rowMeans(tpm[13:15])
tpm$ko_0 <- rowMeans(tpm[16:19])
tpm$ko_20 <- rowMeans(tpm[20:23])
tpm$ko_34 <- rowMeans(tpm[24:27])
tpm$ko_60 <- rowMeans(tpm[28:31])
tpm$wt_0 <- rowMeans(tpm[32:35])
tpm$wt_20 <- rowMeans(tpm[36:39])
tpm$wt_34 <- rowMeans(tpm[40:43])
tpm$wt_60 <- rowMeans(tpm[44:47])

colnames(tpm)
# remove the replicate data
tpm[, c(2:47)] <- NULL
# reorder column
tpm <- tpm[, c(1, 10:13, 2:5, 6:9)]
head(tpm)

# import restables
rstbl_d0_het <- read.table("results/deseq2/rstbl_d0_wt.vs.het.txt", header = TRUE)
rstbl_d20_het <- read.table("results/deseq2/rstbl_d20_wt.vs.het.txt", header = TRUE)
rstbl_d34_het <- read.table("results/deseq2/rstbl_d34_wt.vs.het.txt", header = TRUE)
rstbl_d60_het <- read.table("results/deseq2/rstbl_d60_wt.vs.het.txt", header = TRUE)

rstbl_d0_ko <- read.table("results/deseq2/rstbl_d0_wt.vs.ko.txt", header = TRUE)
rstbl_d20_ko <- read.table("results/deseq2/rstbl_d20_wt.vs.ko.txt", header = TRUE)
rstbl_d34_ko <- read.table("results/deseq2/rstbl_d34_wt.vs.ko.txt", header = TRUE)
rstbl_d60_ko <- read.table("results/deseq2/rstbl_d60_wt.vs.ko.txt", header = TRUE)

#combine the 8 restables
rstbls <- cbind.data.frame(list(rstbl_d0_het, rstbl_d0_ko, rstbl_d20_het, rstbl_d20_ko, rstbl_d34_het, rstbl_d34_ko, rstbl_d60_het, rstbl_d60_ko))

head(rstbls)
colnames(rstbls)
# remove columns non needed
rstbls[, c(4, 7, 11, 14, 18, 21, 25, 28, 32, 35, 39, 42, 38, 49, 45)] <- NULL
colnames(rstbls) <- c("bmean_d0_het", "log2fc_d0_het", "lfcse_d0_het", "pvalue_d0_het", "padj_d0_het", "bmean_d20_het", "log2fc_d20_het", "lfcse_d20_het", "pvalue_d20_het", "padj_d20_het",
"bmean_d34_het", "log2fc_d34_het", "lfcse_d34_het", "pvalue_d34_het", "padj_d34_het",
"bmean_d60_het", "log2fc_d60_het", "lfcse_d60_het", "pvalue_d60_het", "padj_d60_het",
"bmean_d0_ko", "log2fc_d0_ko", "lfcse_d0_ko", "pvalue_d0_ko", "padj_d0_ko", "bmean_d20_ko", "log2fc_d20_ko", "lfcse_d20_ko", "pvalue_d20_ko", "padj_d20_ko",
"bmean_d34_ko", "log2fc_d34_ko", "lfcse_d34_ko", "pvalue_d34_ko", "padj_d34_ko",
"bmean_d60_ko", "log2fc_d60_ko", "lfcse_d60_ko", "pvalue_d60_ko", "padj_d60_ko", "entrezid")

colnames(rstbls)

#Select all DEGs
degs <- subset(rstbls, padj_d0_het < 0.05 | padj_d20_het < 0.05 | padj_d34_het < 0.05 | padj_d60_het < 0.05 | padj_d0_ko < 0.05 | padj_d20_ko < 0.05 | padj_d34_ko < 0.05 | padj_d60_ko < 0.05)
nrow(rstbls)
nrow(degs)
prcnt_degs <- nrow(degs)/nrow(rstbls)*100
prcnt_degs

# too many degs now impose fc cutoffs
degs <- subset(degs, abs(log2fc_d0_het) > 4 | abs(log2fc_d20_het) > 4  | abs(log2fc_d34_het) > 4 | abs(log2fc_d60_het) > 4 
| abs(log2fc_d0_ko) > 4 | abs(log2fc_d20_ko) > 4 | abs(log2fc_d34_ko) > 4 | abs(log2fc_d60_ko) > 4)
nrow(degs)
prcnt_degs2 <- nrow(degs)/nrow(rstbls)*100
prcnt_degs2

# Subset tpm to keep only the rows with symbols in degs
deg_tpm <- tpm[tpm$gene_id %in% rownames(degs), ]
row.names(deg_tpm) <- deg_tpm$gene_id
deg_tpm$gene_id <- NULL
nrow(deg_tpm)
head(deg_tpm)

#reorder columns by genotype -> days
deg_tpm <- deg_tpm[, c(1, 5, 9, 2, 6, 10, 3, 7, 11, 4, 8, 12)]
colnames(deg_tpm)
colnames(deg_tpm) <- c("WT_d0", "HET_d0", "KO_d0", "WT_d20", "HET_d20", "KO_d20", "WT_d34", "HET_d34", "KO_d34", "WT_d60", "HET_d60", "KO_d60")
# convert dataframe into matrix
deg_tpm <- as.matrix(deg_tpm)

# scale data to z scores 
deg_z <- t(scale(t(deg_tpm)))
deg_z <- deg_z[complete.cases(deg_z), ]

#How many clusters are good for K-mean clustering? ####
#I will try to address this according based on the tutorial presented at https://scienceparkstudygroup.github.io/rna-seq-lesson/08-cluster-analysis/index.html#4-looking-at-the-individual-clusters

# SSE method
head(deg_z)
wss2 <- (nrow(deg_z)-1)*sum(apply(deg_z, 2, var))
for (i in 1:20) wss2[i] <- sum(kmeans(deg_z, centers=i)$withinss)
str(wss2)
plot(1:20, wss2, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")
# ~7 clusters seem appropriate

# Average silhouette width method
library(cluster)
sil <- rep(0, 20)
sil

# repeat k-means for 1:20 and extract silhouette:
for(i in 2:20){
  k1to20_1 <- kmeans(deg_z, centers = i, nstart = 25, iter.max = 20)
  ss1 <- silhouette(k1to20_1$cluster, dist(deg_z))
  sil[i] <- mean(ss1[, 3])
}

# Plot the  average silhouette width
plot(1:20, sil, type = "b", pch = 19, xlab = "Number of clusters k", ylab="Average silhouette width")
abline(v = which.max(sil), lty = 2)
cat("Average silhouette width optimal number of clusters:", which.max(sil), "\n")
# it says 2 clusters but that is too small

# hierchical clustering
hr <- hclust(as.dist(1-cor(t(deg_z), method="pearson")), method="complete") # Cluster rows by Pearson correlation.
str(hr)
hc <- hclust(as.dist(1-cor(deg_z, method="spearman")), method="complete") # Clusters columns by Spearman correlation.

# visualize the tree
# make a tree
tree_deg <- as.dendrogram(hr, method="average")

plot(tree_deg,
     leaflab = "none",
     main = "DEG Clustering",
     ylab = "Height")

# cuttrees
hclust1.6 <- cutree(hr, h=1.6)
hclust1.7 <- cutree(hr, h=1.7)
hclust1.8 <- cutree(hr, h=1.8)

# dinamic tree cutting ####
# install.packages('dendextend')
# install.packages("dynamicTreeCut")
library('dynamicTreeCut')

library('dendextend')
clusDyn <- cutreeDynamic(hr, distM = as.matrix(as.dist(1-cor(t(deg_z)))), method = "hybrid")

head(clusDyn)
# dynamic cuts generate too many clusters

# k=7 clustering ####
set.seed(20)
kclust7 <- kmeans(deg_z, centers=7, nstart = 1000, iter.max = 20)
kclusters <- kclust7$cluster
head(kclust7)
head(kclusters)


# plot tree again andadd the cluster vectors 
pdf(sprintf("results/cluster/%s_tree_deg.pdf", Sys.Date()), width = 10, height = 10)

plot(tree_deg,
     leaflab = "none",
     main = "DEG Clustering",
     ylab = "Height")

bars <- cbind(hclust1.6, hclust1.7, hclust1.8, clusDyn, kclusters)

colored_bars(bars, tree_deg, sort_by_lables_order = TRUE, y_shift = -0.1, rowLabels = c("h=1.6", "h=1.7", "h=1.8", "Dynamic", "k=7"), cex.rowLables = 0.7)

#this will add lines showing the cut heights
abline(h=1.6, lty = 2, col="grey")
abline(h=1.7, lty = 2, col="grey")
abline(h=1.8, lty = 2, col="grey")

dev.off()

# Employ tava method to plot median and all genes of each cluster ####
# cut tree and assign cluster# to genes
library(tidyverse)

head(hclust1.6)

hclust1.6 <- cutree(hr, h = 1.6) %>% 
  # turn the named vector into a tibble
  enframe() %>% 
  # rename some of the columns
  rename(gene = name, cluster = value)
hclust1.7 <- cutree(hr, h = 1.7) %>% 
  # turn the named vector into a tibble
  enframe() %>% 
  # rename some of the columns
  rename(gene = name, cluster = value)
hclust1.8 <- cutree(hr, h = 1.8) %>% 
  # turn the named vector into a tibble
  enframe() %>% 
  # rename some of the columns
  rename(gene = name, cluster = value)
kclust7 <- kclust7$cluster %>%
  # turn the named vector into a tibble
  enframe() 
  # rename some of the columns 
  colnames(kclust7)
  head(kclust7)
  colnames(kclust7) <- c('gene', 'cluster')

#covert the Zsocre matrix to data.frame
deg_z_df <- as.data.frame(deg_z)
deg_z_df$gene <- row.names(deg_z_df)
str(deg_z_df)
head(deg_z_df)
str(hclust1.6)
head(hclust1.6)

#add culster names to each gene
c1.6 <- deg_z_df %>%
     inner_join(hclust1.6, by = "gene")
c1.7 <- deg_z_df %>%
     inner_join(hclust1.7, by = "gene")
c1.8 <- deg_z_df %>%
     inner_join(hclust1.8, by = "gene")
k7 <- deg_z_df %>%
     inner_join(kclust7, by = "gene")
head(c1.6)
str(c1.6)


#making the traits file
sample <- colnames(c1.6[c(1:12)])
day <- c(rep(c(0, 20, 34, 60), each = 3))
day
genotype <- c(rep(c("WT", "HET", "KO"), 4))
genotype
traits <- as.data.frame(cbind(sample, day, genotype))
traits
class(traits$day)
traits$day <- as.numeric(as.character(traits$day))
str(traits)

#make all Cluster matrix 
make_cluster_matrix <- function(data) {
  data %>%
    # convert to long format
    pivot_longer(cols = WT_d0:KO_d60, names_to = "sample", values_to = "z") %>%
    # join with sample info table
    full_join(traits, by = ("sample")) %>%
    # for each gene, genotype, and day
    group_by(gene, day, genotype)
}

cm1.6 <- make_cluster_matrix(c1.6)
head(cm1.6)
cm1.7 <- make_cluster_matrix(c1.7)
cm1.8 <- make_cluster_matrix(c1.8)
cmk7 <- make_cluster_matrix(k7)
head(cmk7)
cmk7$cluster

# write cluster dataframe for pathway analysis
write.table(cm1.8, file = "results/cluster/clh.txt", sep = "\t")
write.table(cmk7, file = "results/cluster/clk.txt", sep = "\t")

# Plot gene expression of each cluster ####
library("colorspace")
my_colors2 <- sequential_hcl(3, palette = "TealGrn") 
my_colors2


cluster_lplot <- function(data) {
  # Convert the genotype column to factor with the specified order
  genotype_order <- c("WT", "HET", "KO")
  data$genotype <- factor(data$genotype, levels = genotype_order)
  # line plot
  data %>% 
  ggplot(aes(day, z, group = genotype)) +
  geom_line(stat = "summary", fun = "median", size = 1.0, aes(col = genotype)) +
  scale_color_manual(values = my_colors2) +
  scale_x_continuous(minor_breaks = NULL, breaks= c(0, 20, 34, 60)) +
  xlab("Time") +
  ylab("Expression") +
  facet_grid(cols = vars(cluster)) +
  theme(
    panel.background = element_rect(fill = NA, color = NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(colour = "black"))
}

lp1.6 <- cluster_lplot(cm1.6) 
lp1.7 <- cluster_lplot(cm1.7)
lp1.8 <- cluster_lplot(cm1.8)
lpk7 <- cluster_lplot(cmk7)
lp1.8
lpk7

# make cluster line plots
library(gridExtra)
pdf(sprintf("results/cluster/%s_cluster_lineplot.pdf", Sys.Date()), width = 12, height = 9)
grid.arrange(
  arrangeGrob(lp1.6, ncol = 1, top = "Hiearchical: 1.6"), 
  arrangeGrob(lp1.7, ncol = 1, top = "Hiearchical:1.7"),
  arrangeGrob(lp1.8, ncol = 1, top = "Hiearchical:1.8"), 
  arrangeGrob(lpk7, ncol = 1, top = "K-means: 7"),
  ncol = 1, heights = c(1, 1, 1, 1)
)
dev.off()

pdf(sprintf("results/cluster/%s_lineplot2.pdf", Sys.Date()), width = 12, height = 9)
grid.arrange(
  arrangeGrob(lp1.8, ncol = 1, top = "Hiearchical:1.8"), 
  arrangeGrob(lpk7, ncol = 1, top = "K-means: 7"),
  ncol = 1, heights = c(1, 1)
)
dev.off()

# know the number of genes in each cluster
table(c1.8$cluster)
table(k7$cluster)

# define color for heatmap
# install.packages("colorspace")
library("colorspace")
hcl_palettes(plot = TRUE)
br <- diverging_hcl(10, palette = "Blue-Red3") 


# try ComplexHeatmap 
# ref: https://www.datanovia.com/en/lessons/heatmap-in-r-static-and-interactive-visualization/
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("ComplexHeatmap")

library(ComplexHeatmap)

df <- scale(deg_z)

# To specify a custom colors, you must use the the colorRamp2() function [circlize package], as follow:
library(circlize)
library("RColorBrewer")

# We can also customize the appearance of dendograms using the function color_branches() [dendextend package]:
# install.packages("gridtext")
library(gridtext)
library(dendextend)
str(hr)
head(df)
split <- data.frame(cutree(hr, h=1.8))
head(split)
colnames(df) = rep(c("WT", "HET", "KO"), 4)

pdf(sprintf("results/cluster/%s_degheatmap_hc.pdf", Sys.Date()), width = 9, height = 10)

hm1 <- Heatmap(df, name = "Z score (Hiearchical)", 
        row_names_gp = gpar(fontsize = 6.5),
        row_split = split,
        col = br,
        show_row_names = FALSE,
        row_title = "C%s",
        cluster_columns = FALSE,
        column_split = rep(c("d0", "d20", "d34", "d60"), each = 3),
        column_names_side = "top",
        column_names_rot = 0, 
        width = unit(15, "cm"))
hm1 

dev.off() 

# split by k-mean clustering
library("cluster")

pdf(sprintf("results/cluster/%s_degheatmap_km.pdf", Sys.Date()), width = 9, height = 10)

hm2 <- Heatmap(df, name = "Z score (K-means)", 
        row_names_gp = gpar(fontsize = 6.5),
        row_split = k7$cluster,
        col = br,
        show_row_names = FALSE,
        row_title = "C%s",
        cluster_columns = FALSE,
        column_split = rep(c("d0", "d20", "d34", "d60"), each = 3),
        column_names_side = "top",
        column_names_rot = 0,
        width = unit(15, "cm")
)

hm2

dev.off() 


