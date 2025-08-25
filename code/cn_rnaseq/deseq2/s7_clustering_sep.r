# perform clusteirng on the DE genes HET and KO separtely. 

# import tpm file
tpm <- read.delim("data/blueberry_output/combine_counts/gene_TPM.txt", header = TRUE, sep = "\t")
colnames(tpm) <- gsub("Rai1.|.TPM|d", "", colnames(tpm))
colnames(tpm)[1] <- "gene_id"

# remove het.d34.rep1 and het.d60.rep4
tpm[, c("het_34_rep1", "het_60_rep4")] <- NULL
colnames(tpm)

# calrulate tpm average between replicates
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

# select degs with a stringent cutoff (padj < 0.01, |log2fc| > 2)
d0_het_deg3 <- subset(rstbl_d0_het, padj < 0.01 & abs(log2FoldChange) > 2)
d0_ko_deg3 <- subset(rstbl_d0_ko, padj < 0.01 & abs(log2FoldChange) > 2)
d20_het_deg3 <- subset(rstbl_d20_het, padj < 0.01 & abs(log2FoldChange) > 2)
d20_ko_deg3 <- subset(rstbl_d20_ko, padj < 0.01 & abs(log2FoldChange) > 2)
d34_het_deg3 <- subset(rstbl_d34_het, padj < 0.01 & abs(log2FoldChange) > 2)
d34_ko_deg3 <- subset(rstbl_d34_ko, padj < 0.01 & abs(log2FoldChange) > 2)
d60_het_deg3 <- subset(rstbl_d60_het, padj < 0.01 & abs(log2FoldChange) > 2)
d60_ko_deg3 <- subset(rstbl_d60_ko, padj < 0.01 & abs(log2FoldChange) > 2)

# combine rownames of selected degs into a list
het_deg3 <- c(rownames(d0_het_deg3), rownames(d20_het_deg3), rownames(d34_het_deg3), rownames(d60_het_deg3))
ko_deg3 <- c(rownames(d0_ko_deg3), rownames(d20_ko_deg3), rownames(d34_ko_deg3), rownames(d60_ko_deg3))

# remove the duplicated genes, which change at multiple time points 
het_deg3 <- unique(het_deg3)
het_deg3
ko_deg3 <- unique(ko_deg3)
ko_deg3

# #Select het DEGs
# degs_het <- subset(rstbls, padj_d0_het < 0.01 | padj_d20_het < 0.01 | padj_d34_het < 0.01 | padj_d60_het < 0.01)
# nrow(rstbls)
# nrow(degs_het)
# prcnt_degs_het <- nrow(degs_het)/nrow(rstbls)*100
# prcnt_degs_het

# # further select with log2fc > 2
# degs_het2 <- subset(degs_het, abs(log2fc_d0_het) > 2 | abs(log2fc_d20_het) > 2  | abs(log2fc_d34_het) > 2 | abs(log2fc_d60_het) > 2)
# nrow(degs_het2)
# prcnt_degs_het2 <- nrow(degs_het2)/nrow(rstbls)*100
# prcnt_degs_het2

# #Select KO DEGs
# degs_ko <- subset(rstbls, padj_d0_ko < 0.01 | padj_d20_ko < 0.01 | padj_d34_ko < 0.01 | padj_d60_ko < 0.01)
# nrow(degs_ko)
# prcnt_degs_ko <- nrow(degs_ko)/nrow(rstbls)*100
# prcnt_degs_ko

# # further select with log2fc > 3
# degs_ko2 <- subset(degs_ko, abs(log2fc_d0_ko) > 2 | abs(log2fc_d20_ko) > 2  | abs(log2fc_d34_ko) > 2 | abs(log2fc_d60_ko) > 2)
# nrow(degs_ko2)
# prcnt_degs_ko2 <- nrow(degs_ko2)/nrow(rstbls)*100
# prcnt_degs_ko2

# Subset tpm table only with het DEGs
deg_het_tpm <- tpm[tpm$gene_id %in% het_deg3, ]
row.names(deg_het_tpm) <- deg_het_tpm$gene_id
deg_het_tpm$gene_id <- NULL
nrow(deg_het_tpm)
head(deg_het_tpm)

# Subset tpm table only with ko DEGs
deg_ko_tpm <- tpm[tpm$gene_id %in% ko_deg3, ]
row.names(deg_ko_tpm) <- deg_ko_tpm$gene_id
deg_ko_tpm$gene_id <- NULL
nrow(deg_ko_tpm)
head(deg_ko_tpm)

# reorder columns by genoetype -> days
deg_het_tpm <- deg_het_tpm[, c(1, 5, 2, 6, 3, 7, 4, 8)]
colnames(deg_het_tpm) <- c("WT_d0", "HET_d0", "WT_d20", "HET_d20", "WT_d34", "HET_d34", "WT_d60", "HET_d60")

deg_ko_tpm <- deg_ko_tpm[, c(1, 9, 2, 10, 3, 11, 4, 12)]
colnames(deg_ko_tpm) <- c("WT_d0", "KO_d0", "WT_d20", "KO_d20", "WT_d34", "KO_d34", "WT_d60", "KO_d60")

# convert datafrmae to matrix 
deg_het_tpm <- as.matrix(deg_het_tpm)
deg_ko_tpm <- as.matrix(deg_ko_tpm)

# scale data to z scores
deg_het_z <- t(scale(t(deg_het_tpm)))
deg_ko_z <- t(scale(t(deg_ko_tpm)))

# check if na values are present
sum(is.na(deg_het_z))
sum(is.na(deg_ko_z))

# remove na values
deg_ko_z <- na.omit(deg_ko_z)

# how many cluters are good for K-means clustering?
# SSE method
head(deg_het_z)
wss_het <- (nrow(deg_het_z)-1)*sum(apply(deg_het_z, 2, var))
for (i in 1:20) wss_het[i] <- sum(kmeans(deg_het_z, centers=i)$withinss)
str(wss_het)
plot(1:20, wss_het, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")
# ~5 cluster seems appropriate 
dev.off()

head(deg_ko_z)
wss_ko <- (nrow(deg_ko_z)-1)*sum(apply(deg_ko_z, 2, var))
for (i in 1:20) wss_ko[i] <- sum(kmeans(deg_ko_z, centers=i)$withinss)
str(wss_ko)
plot(1:20, wss_ko, type="b", xlab="Number of Clusters", ylab="Within groups sum of squares")
# ~5 cluster seems appropriate 

# hierchical clustering
hr_het <- hclust(as.dist(1-cor(t(deg_het_z), method="pearson")), method="complete") # Cluster rows by Pearson correlation.
str(hr_het)
hc_het <- hclust(as.dist(1-cor(deg_het_z, method="spearman")), method="complete") # Clusters columns by Spearman correlation.

hr_ko <- hclust(as.dist(1-cor(t(deg_ko_z), method="pearson")), method="complete") # Cluster rows by Pearson correlation.
str(hr_ko)
hc_ko <- hclust(as.dist(1-cor(deg_ko_z, method="spearman")), method="complete") # Clusters columns by Spearman correlation.

# visualize the dendrograms
tree_deg_het <- as.dendrogram(hr_het, method = "average")
tree_deg_ko <- as.dendrogram(hr_ko, method = "average")

# save the dendrograms as a pdf 
pdf(sprintf("results/cluster/sep/%s-dendrograms_sep.pdf", Sys.Date()), width=10, height=5)
par(mfrow=c(1,2))
plot(tree_deg_het, main = "Hierarchical clustering of HET-DEGs", leaflab = "none", ylab="Height")
plot(tree_deg_ko, main = "Hierarchical clustering of KO-DEGs", leaflab = "none", ylab="Height")
dev.off()

# cut the dendrogram to generate 6 clusters
hclust_het_1.75 <- cutree(hr_het, h=1.75)
hclust_ko_1.8 <- cutree(hr_ko, h=1.8)

unique(hclust_het_1.75)
unique(hclust_ko_1.8)

# K-means clustering K = 6 or 5
set.seed(20)
kclust_het <- kmeans(deg_het_z, centers = 6, nstart=1000, iter.max=20)
kclust_ko <- kmeans(deg_ko_z, centers = 5, nstart=1000, iter.max=20)
k_het <- kclust_het$cluster
k_ko <- kclust_ko$cluster

head(kclust_het)

# plot tre again and add the cluster information
library('dendextend')
library(gridExtra)
library(ggplot2)

pdf(sprintf("results/cluster/sep/%s-dendrograms_sep.pdf", Sys.Date()), width=10, height=5)
par(mfrow=c(1,2))

c1 <- plot(tree_deg_het, main = "Hierarchical clustering of HET-DEGs", leaflab = "none", ylab="Height")
c2 <- plot(tree_deg_ko, main = "Hierarchical clustering of KO-DEGs", leaflab = "none", ylab="Height")
barko <- cbind(hclust_ko_1.8, k_ko)
colored_bars(barko, tree_deg_ko, sort_by_lables_order = TRUE, y_shift = -0.1, rowLabels = c("h=1.8", "k=5"), cex.rowLables = 0.7)

dev.off()

# make line plot for each cluster
install.packages("tidyverse")
library(tidyverse)
str(hclust_het_1.75)
head(hclust_het_1.75)
class(hclust_het_1.75)
str(hr_het)


hclust_het <- cutree(hr_het, h = 1.75) %>% enframe() %>% rename(gene = name, cluster = value)
hclust_ko <- cutree(hr_ko, h = 1.8) %>% enframe() %>% rename(gene = name, cluster = value)
k_het <- k_het %>% enframe() %>% rename(gene = name, cluster = value)
k_ko <- k_ko %>% enframe() %>% rename(gene = name, cluster = value)
hclust_het
str(k_het)
str(k_ko)

# convert the Z-scores matric to data.frame
deg_het_z_df <- as.data.frame(deg_het_z)
deg_ko_z_df <- as.data.frame(deg_ko_z)
deg_het_z_df$gene <- row.names(deg_het_z_df) 
deg_ko_z_df$gene <- row.names(deg_ko_z_df)

str(deg_het_z_df)
head(deg_het_z_df)

# add culster names to each gene
h_het <- deg_het_z_df %>% left_join(hclust_het, by = "gene")
k_het <- deg_het_z_df %>% left_join(k_het, by = "gene")
h_ko <- deg_ko_z_df %>% left_join(hclust_ko, by = "gene")
k_ko <- deg_ko_z_df %>% left_join(k_ko, by = "gene")

str(h_het)
str(k_het)
str(h_ko)
str(k_ko)

# make the traits file 
sample_het <- colnames(deg_het_z_df)[1:8]
day <- c(rep(c(0, 20, 34, 60), each = 2))
genotype_het <- c(rep(c("WT", "HET"), 4))
traits_het <- as.data.frame(cbind(sample_het, day, genotype_het))
traits_het

sample_ko <- colnames(deg_ko_z_df)[1:8]
genotype_ko <- c(rep(c("WT", "KO"), 4))
traits_ko <- as.data.frame(cbind(sample_ko, day, genotype_ko))
traits_ko

class(traits_het$day)
traits_het$day <- as.numeric(as.character(traits_het$day))
traits_ko$day <- as.numeric(as.character(traits_ko$day))

# make cluster matrix 
make_cluster_matrix_het <- function(data) {
  data %>%
    # convert to long format
    pivot_longer(cols = WT_d0:HET_d60, names_to = "sample_het", values_to = "z") %>%
    # join with sample info table
    full_join(traits_het, by = ("sample_het")) %>%
    # for each gene, genotype, and day
    group_by(gene, day, genotype_het)
}

mh_het <- make_cluster_matrix_het(h_het) 
mk_het <- make_cluster_matrix_het(k_het)

make_cluster_matrix_ko <- function(data) {
  data %>%
    # convert to long format
    pivot_longer(cols = WT_d0:KO_d60, names_to = "sample_ko", values_to = "z") %>%
    # join with sample info table
    full_join(traits_ko, by = ("sample_ko")) %>%
    # for each gene, genotype, and day
    group_by(gene, day, genotype_ko)
}

mh_ko <- make_cluster_matrix_ko(h_ko)
mk_ko <- make_cluster_matrix_ko(k_ko)

# check cluster matrix
mh_het
mk_het
mh_ko
mk_ko

# write cluster dataframe for pathway analysis
write.table(mh_het, file = "results/cluster/sep/mh_het.txt", sep = "\t")
write.table(mk_het, file = "results/cluster/sep/mk_het.txt", sep = "\t")
write.table(mh_ko, file = "results/cluster/sep/mh_ko.txt", sep = "\t")
write.table(mk_ko, file = "results/cluster/sep/mk_ko.txt", sep = "\t")

# Plot gene expression of each cluster ####
color1 <- c("#757582", "#5876D1") 

cluster_lplot1 <- function(data) {
  # Convert the genotype column to factor with the specified order
  genotype_order <- c("WT", "HET")
  data$genotype_het <- factor(data$genotype_het, levels = genotype_order)
  # line plot
  data %>% 
  ggplot(aes(day, z, group = genotype_het)) +
  geom_line(stat = "summary", fun = "median", size = 1.0, aes(col = genotype_het)) +
  scale_color_manual(values = color1) +
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

lp_h_het <- cluster_lplot1(mh_het)
lp_k_het <- cluster_lplot1(mk_het)

cluster_lplot2 <- function(data) {
  # Convert the genotype column to factor with the specified order
  genotype_order1 <- c("WT", "KO")
  data$genotype_ko <- factor(data$genotype_ko, levels = genotype_order1)
  # line plot
  data %>% 
  ggplot(aes(day, z, group = genotype_ko)) +
  geom_line(stat = "summary", fun = "median", size = 1.0, aes(col = genotype_ko)) +
  scale_color_manual(values = color1) +
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

lp_h_ko <- cluster_lplot2(mh_ko)
lp_k_ko <- cluster_lplot2(mk_ko)

# arrange lp_h_het and lp_k_het in two rows and save them as one pdf.
# add the title "Hiearchical clustering" to lp_h_het and "K-means clustering" to lp_k_het

library(gridExtra)
pdf(sprintf("results/cluster/sep/%s_cluster_lineplot_het.pdf", Sys.Date()), width = 12, height = 9)

lp_h_het <- lp_h_het + ggtitle("Hierarchical clustering")
lp_k_het <- lp_k_het + ggtitle("K-means clustering")

grid.arrange(lp_h_het, lp_k_het, ncol = 1)
dev.off()

# save ko line plots
pdf(sprintf("results/cluster/sep/%s_cluster_lineplot_ko.pdf", Sys.Date()), width = 12, height = 9)

lp_h_ko <- lp_h_ko + ggtitle("Hierarchical clustering")
lp_k_ko <- lp_k_ko + ggtitle("K-means clustering")

grid.arrange(lp_h_ko, lp_k_ko, ncol = 1)
dev.off()
dev.off()
dev.off()

# know the number of genes in each cluster
table(mh_het$cluster)
table(mk_het$cluster)
table(mh_ko$cluster)
table(mk_ko$cluster)

# heatmap
## define color for heatmap
library("colorspace")
hcl_palettes(plot = TRUE)
br <- diverging_hcl(10, palette = "Blue-Red3")

library(ComplexHeatmap)
library(gridtext)
library(dendextend)

str(hr_het)
head(deg_het_z)
head(hclust_het_1.75)

library(ComplexHeatmap)

drawheatmap <- function (matrix, split) {Heatmap(matrix, name = "Z score", 
        row_names_gp = gpar(fontsize = 6.5),
        row_split = split,
        col = br,
        show_row_names = FALSE,
        row_title = "C%s",
        cluster_columns = FALSE,
        column_split = rep(c("d0", "d20", "d34", "d60"), each = 2),
        column_names_side = "top",
        column_names_rot = 0, 
        width = unit(15, "cm"))}

# K-measns cluster separation (we may need to fix this part for k-means heatmap to work)
k_het1 <- kclust_het$cluster
k_ko1 <- kclust_ko$cluster

str(hclust_het_1.75)
str(k_het1)

hm_het_hi <- drawheatmap(deg_het_z, hclust_het_1.75)
hm_ko_hi <- drawheatmap(deg_ko_z, hclust_ko_1.8)
hm_het_k <- drawheatmap(deg_het_z, k_het1)
hm_ko_k <- drawheatmap(deg_ko_z, k_ko1)

# save heatmap
pdf(sprintf("results/cluster/sep/%s_degheatmap_het_hi.pdf", Sys.Date()), width = 9, height = 10)
hm_het_hi
dev.off()

pdf(sprintf("results/cluster/sep/%s_degheatmap_ko_hi.pdf", Sys.Date()), width = 9, height = 10)
hm_ko_hi
dev.off()

pdf(sprintf("results/cluster/sep/%s_degheatmap_het_k.pdf", Sys.Date()), width = 9, height = 10)
hm_het_k
dev.off()

pdf(sprintf("results/cluster/sep/%s_degheatmap_ko_k.pdf", Sys.Date()), width = 9, height = 10)
hm_ko_k
dev.off()
