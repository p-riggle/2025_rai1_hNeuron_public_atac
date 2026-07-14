# this code reads in top300 pc1 genes from ENCODE data (identified in 1_PCA_ENCODE.R),
# performs kmeans clustering, and plots lineplots for each cluster (Fig 7B)
# makes boxplot comparing Rai1 ko and wt from Huang et al data (Fig 7D)

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(reshape2)

# import zscore of FPKM of PC1 genes from 1_PCA_ENCODE_WT.R script
pc1_genes <- read.table("result/Fig7_huang/encode_pc1_top300_genes_FPKM_zscore.tsv", header = TRUE) # data comes from 1_deseq2_ENCODE_jan script 

head(pc1_genes)

row.names(pc1_genes) <- pc1_genes$ensembl_id
pc1_genes$ensembl_id <- NULL
nrow(pc1_genes) # 300 as expected
head(pc1_genes)

# convert datafrmae to matrix 
pc1_genes <- as.matrix(pc1_genes)

# check if na values are present
sum(is.na(pc1_genes)) # none are present

# K-means clustering K = 2
set.seed(20)
kclust <- kmeans(pc1_genes, centers = 2, nstart=1000, iter.max=20)
k <- kclust$cluster

#visualize the kmeans centers to interpret:
centers_df <- as.data.frame(kclust$centers)
centers_df$cluster <- factor(1:2)

centers_long <- centers_df %>% pivot_longer(-cluster, names_to="timepoint", values_to="Zscore")


kmeans_centers_plot <- ggplot(centers_long, aes(x=timepoint, y=Zscore, color=cluster, group=cluster)) +
  geom_line(size=1.5) +
  geom_point(size=3) +
  theme_minimal() +
  labs(title="K-means Cluster Centers (Gene Expression over Time)", x="Timepoint", y="Z-score")

kmeans_centers_plot

# how many genes in each cluster?
table(k) # cluster 1 has 249 genes, cluster 2 has 51 genes

head(kclust)

# using the k means clustering only, make the clusters for the ENCODE wildtype and  plot as lineplot

k <- enframe(k)

k <- dplyr::rename(k, gene = name, cluster = value)

# make pc1 genes dataframe
pc1_genes_df <- as.data.frame(pc1_genes)

pc1_genes_df2 <- pc1_genes_df %>%
  tibble::rownames_to_column(var = "ensembl_id")

k <- dplyr::rename(k, ensembl_id = gene)
# add cluster names to each gene
pc1_genes_df2 <- pc1_genes_df2 %>% left_join(k, by = "ensembl_id")
str(k)

# make traits table
sample <- colnames(pc1_genes_df2)[2:6]
day <- as.numeric(sub("FPKM_p", "", sample)) 
genotype <- rep("WT", length(sample))

traits <- data.frame(
  sample = sample,
  day = day,
  genotype = genotype
)

traits 

# make cluster matrix
make_cluster_matrix <- function(data) {
  data %>%
    # convert to long format
    pivot_longer(cols = FPKM_p10:FPKM_p60, names_to = "sample", values_to = "z") %>%
    # join with sample info table
    full_join(traits, by = ("sample")) %>%
    # for each gene, genotype, and day
    group_by(ensembl_id, day, genotype)
}

mk <- make_cluster_matrix(pc1_genes_df2)

mk

# write cluster dataframe for pathway analysis
write.table(mk, file = "result/Fig7_huang/mk_ENCODE_WTctx_pc1300genesclusterdataframe.txt", sep = "\t")

# plot gene expression of each cluster!
color1 <- c("#7b7b7b") # grey for wildtype

cluster_lplot <- function(data) {
  # Convert the genotype column to factor with the specified order
  genotype_order <- c("WT")
  data$genotype <- factor(data$genotype, levels = genotype_order)
  # line plot
  data %>% 
    ggplot(aes(day, z, group = genotype)) +
    geom_line(stat = "summary", fun = "median", linewidth = 1.0, aes(col = genotype)) +
    scale_color_manual(values = color1) +
    scale_x_continuous(minor_breaks = NULL, breaks= c(10, 14, 25, 36, 60)) +
    xlab("Time") +
    ylab("Expression") +
    facet_grid(cols = vars(cluster)) +
    theme(
      panel.background = element_rect(fill = NA, color = NA),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"))
}

lp_k <- cluster_lplot(mk)

lp_k

# add grey shading around the line for each cluster indicating the expression of all genes in the cluster
lp_k_genes <- lp_k + geom_line(data = mk, aes(x = day, y = z, group = ensembl_id), color = "grey", alpha = 0.1)
lp_k_genes

# save lineplot
ggsave("result/Fig7_huang/ENCODE_cluster_lineplot_WTctxPC1300genes_Fig7B.pdf", plot = lp_k_genes, width = 12, height = 9, dpi = 300)

# ---------------------------------------------------------------------------------
## plot the cluster 1 genes in WT and KO RAI1 cortex (huang data) as a boxplot
# 1. subset the cluster 1 and 2 genes into the Huang data
huang_ctx <- read.csv("data/Fig7_huang/cleaned_huang_data.csv")

cluster1_genes_zscore <- mk[mk$cluster == 1, ]

# subset out genes in cluster 1 that are in the Huang dataset (3wk cortex)
subset_huang_ctx <- huang_ctx[huang_ctx$gene_id %in% cluster1_genes_zscore$ensembl_id, ]

# get rid of the non-cortex columns 4-7
subset_huang_ctx <- subset_huang_ctx[, -c(4:7)]

# make boxplot of WT and KO 3 week ctx FPKM values
# zscore 
subset_huang_ctx$WT_ctx_z <- scale(subset_huang_ctx$WT_ctx)
subset_huang_ctx$KO_ctx_z <- scale(subset_huang_ctx$KO_ctx)

head(subset_huang_ctx)

# get rid of the FPKM columns 
subset_huang_ctx_z <- subset_huang_ctx[, -c(2:3)]

# plot the cluster 1 pc1 genes in the Huang Ko and WT 3 wk cortex

filtered_plot_z <- melt(subset_huang_ctx_z,
                        measure.vars=c("WT_ctx_z", "KO_ctx_z"),
                        variable.name = "Genotype",
                        value.name = "Zscore"
)

filtered_plot_z$Genotype <- factor(filtered_plot_z$Genotype,
                                   levels=c("WT_ctx_z", "KO_ctx_z")
)

filtered_plot_z$Genotype <- sub("_.*", "", filtered_plot_z$Genotype)
filtered_plot_z$Genotype <- factor(filtered_plot_z$Genotype, levels = c("WT", "KO"))

# make boxplot of cluster 1 pc1 genes in Huang 3 wk ctx
boxplot2 <- ggplot(filtered_plot_z, aes(x = Genotype, y = Zscore, fill = Genotype)) +
  geom_boxplot(trim = FALSE, color = "#333333") +
  scale_fill_manual(values = c(WT="#7b7b7b", KO="#bd6565")) +
  # geom_jitter(width = 0.2, size = 1, alpha = 0.7) + # Shows data points
  theme_minimal() +
  labs(title = "Zscore of p21 RAI1 cKO mouse ctx",
       x = "Genotype",
       y = "Z-Score") +
  theme(text = element_text(size = 14))

boxplot2
ggsave("result/Fig7_huang/boxplotZscore_ENCODEcluster1pc1300genes_Huang3wkctx_Fig7D.pdf", plot = boxplot2, width = 10, height = 6, dpi = 300)

# mann whitney u test of difference between WT and KO for cluster 1 pc1 genes in Huang 3 wk ctx
wilcox.test(Zscore ~ Genotype, data = filtered_plot_z)
# W = 26007, p-value = 0.639


# _________________________________________________________________
## plot the cluster 2 pc1 genes in WT and KO RAI1 cortex (huang data) as a boxplot
# 1. subset the cluster 2 genes into the Huang data
huang_ctx <- read.csv("data/huang/cleaned_huang_data.csv")

cluster2_genes_zscore <- mk[mk$cluster == 2, ]

# subset out genes in cluster 1 that are in the Huang dataset (3wk cortex)
subset_huang_ctx <- huang_ctx[huang_ctx$gene_id %in% cluster2_genes_zscore$ensembl_id, ]

# get rid of the non-cortex columns 4-7
subset_huang_ctx <- subset_huang_ctx[, -c(4:7)]

head(subset_huang_ctx)
# zscore 
subset_huang_ctx$WT_ctx_z <- scale(subset_huang_ctx$WT_ctx)
subset_huang_ctx$KO_ctx_z <- scale(subset_huang_ctx$KO_ctx)

head(subset_huang_ctx)

# get rid of the FPKM  and log2 columns 4-5
subset_huang_ctx_z <- subset_huang_ctx[, -c(2:3)]

# plot the boxlot of cluster 2 pc1 genes in the Huang Ko and WT 3 wk cortex

filtered_plot_z <- melt(subset_huang_ctx_z,
                        measure.vars=c("WT_ctx_z", "KO_ctx_z"),
                        variable.name = "Genotype",
                        value.name = "Zscore"
)

filtered_plot_z$Genotype <- factor(filtered_plot_z$Genotype,
                                   levels=c("WT_ctx_z", "KO_ctx_z")
)

filtered_plot_z$Genotype <- sub("_.*", "", filtered_plot_z$Genotype)
filtered_plot_z$Genotype <- factor(filtered_plot_z$Genotype, levels = c("WT", "KO"))

# make boxplot
boxplot3 <- ggplot(filtered_plot_z, aes(x = Genotype, y = Zscore, fill = Genotype)) +
  geom_boxplot(trim = FALSE, color = "#333333") +
  scale_fill_manual(values = c(WT="#7b7b7b", KO="#bd6565")) +
  #geom_jitter(width = 0.2, size = 1, alpha = 0.7) + # Shows data points
  theme_minimal() +
  labs(title = "Zscore of p21 RAI1 cKO mouse ctx",
       x = "Genotype",
       y = "Z-Score") +
  theme(text = element_text(size = 14))

boxplot3 

ggsave("result/Fig7_huang/boxplotZscore_ENCODEcluster2pc1300genes_Huang3wkctx_Fig7D.pdf", plot = boxplot3, width = 10, height = 6, dpi = 300)

# mann whitney u test for difference between WT and KO for cluster 2 pc1 genes in Huang 3 wk ctx
wilcox.test(Zscore ~ Genotype, data = filtered_plot_z)
# W = 945, p-value = 0.8636

# end of script