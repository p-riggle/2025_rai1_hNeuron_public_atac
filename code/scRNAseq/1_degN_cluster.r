# load DEseq2
library(DESeq2)

# test import the csv file 
csv20 <- read.table("results/scRNAseq/pseudobulk/d20/Day20_rlogNorm.csv", header = TRUE, sep = ",")
head(csv20)

# import the deseq2 results table (d20)
rstb_d20_c0_het <- read.table("results/scRNAseq/pseudobulk/d20/Clust.0_Day20_HET_vs_WT.txt", header = TRUE, sep = "\t")
rstb_d20_c0_homo <- read.table("results/scRNAseq/pseudobulk/d20/Clust.0_Day20_KO_vs_WT.txt", header = TRUE, sep = "\t")
rstb_d20_c1_het <- read.table("results/scRNAseq/pseudobulk/d20/Clust.1_Day20_HET_vs_WT.txt", header = TRUE, sep = "\t")
rstb_d20_c1_homo <- read.table("results/scRNAseq/pseudobulk/d20/Clust.1_Day20_KO_vs_WT.txt", header = TRUE, sep = "\t")
rstb_d20_c2_het <- read.table("results/scRNAseq/pseudobulk/d20/Clust.2_Day20_HET_vs_WT.txt", header = TRUE, sep = "\t")
rstb_d20_c2_homo <- read.table("results/scRNAseq/pseudobulk/d20/Clust.2_Day20_KO_vs_WT.txt", header = TRUE, sep = "\t")
rstb_d20_c3_het <- read.table("results/scRNAseq/pseudobulk/d20/Clust.3_Day20_HET_vs_WT.txt", header = TRUE, sep = "\t")
rstb_d20_c3_homo <- read.table("results/scRNAseq/pseudobulk/d20/Clust.3_Day20_KO_vs_WT.txt", header = TRUE, sep = "\t")
rstb_d20_c4_het <- read.table("results/scRNAseq/pseudobulk/d20/Clust.4_Day20_HET_vs_WT.txt", header = TRUE, sep = "\t")
rstb_d20_c4_homo <- read.table("results/scRNAseq/pseudobulk/d20/Clust.4_Day20_KO_vs_WT.txt", header = TRUE, sep = "\t")

# import the deseq2 results table (d60)
rstb_d60_c0_het <- read.table("results/scRNAseq/pseudobulk/d60/Clust.0_Day60_HET_vs_WT.txt", header = TRUE, sep = "\t")
rstb_d60_c0_homo <- read.table("results/scRNAseq/pseudobulk/d60/Clust.0_Day60_KO_vs_WT.txt", header = TRUE, sep = "\t")
rstb_d60_c1_het <- read.table("results/scRNAseq/pseudobulk/d60/Clust.1_Day60_HET_vs_WT.txt", header = TRUE, sep = "\t")
rstb_d60_c1_homo <- read.table("results/scRNAseq/pseudobulk/d60/Clust.1_Day60_KO_vs_WT.txt", header = TRUE, sep = "\t")
rstb_d60_c2_het <- read.table("results/scRNAseq/pseudobulk/d60/Clust.2_Day60_HET_vs_WT.txt", header = TRUE, sep = "\t")
rstb_d60_c2_homo <- read.table("results/scRNAseq/pseudobulk/d60/Clust.2_Day60_KO_vs_WT.txt", header = TRUE, sep = "\t")
rstb_d60_c3_het <- read.table("results/scRNAseq/pseudobulk/d60/Clust.3_Day60_HET_vs_WT.txt", header = TRUE, sep = "\t")
rstb_d60_c3_homo <- read.table("results/scRNAseq/pseudobulk/d60/Clust.3_Day60_KO_vs_WT.txt", header = TRUE, sep = "\t")
rstb_d60_c4_het <- read.table("results/scRNAseq/pseudobulk/d60/Clust.4_Day60_HET_vs_WT.txt", header = TRUE, sep = "\t")
rstb_d60_c4_homo <- read.table("results/scRNAseq/pseudobulk/d60/Clust.4_Day60_KO_vs_WT.txt", header = TRUE, sep = "\t")
rstb_d60_c5_het <- read.table("results/scRNAseq/pseudobulk/d60/Clust.5_Day60_HET_vs_WT.txt", header = TRUE, sep = "\t")
rstb_d60_c5_homo <- read.table("results/scRNAseq/pseudobulk/d60/Clust.5_Day60_KO_vs_WT.txt", header = TRUE, sep = "\t")
rstb_d60_c6_het <- read.table("results/scRNAseq/pseudobulk/d60/Clust.6_Day60_HET_vs_WT.txt", header = TRUE, sep = "\t")
rstb_d60_c6_homo <- read.table("results/scRNAseq/pseudobulk/d60/Clust.6_Day60_KO_vs_WT.txt", header = TRUE, sep = "\t")
rstb_d60_c7_het <- read.table("results/scRNAseq/pseudobulk/d60/Clust.7_Day60_HET_vs_WT.txt", header = TRUE, sep = "\t")
rstb_d60_c7_homo <- read.table("results/scRNAseq/pseudobulk/d60/Clust.7_Day60_KO_vs_WT.txt", header = TRUE, sep = "\t")

colnames(deg_d20_c0_het)

# from the results tables extract the differentially expressed genes (padj < 0.01, abs(log2FoldChange)>2) and make tables 
deg_d20_c0_het <- subset(rstb_d20_c0_het, padj < 0.01, abs(log2FoldChange)>2)
deg_d20_c0_homo <- subset(rstb_d20_c0_homo, padj < 0.01, abs(log2FoldChange)>2)
deg_d20_c1_het <- subset(rstb_d20_c1_het, padj < 0.01, abs(log2FoldChange)>2)
deg_d20_c1_homo <- subset(rstb_d20_c1_homo, padj < 0.01, abs(log2FoldChange)>2)
deg_d20_c2_het <- subset(rstb_d20_c2_het, padj < 0.01, abs(log2FoldChange)>2)
deg_d20_c2_homo <- subset(rstb_d20_c2_homo, padj < 0.01, abs(log2FoldChange)>2)
deg_d20_c3_het <- subset(rstb_d20_c3_het, padj < 0.01, abs(log2FoldChange)>2)
deg_d20_c3_homo <- subset(rstb_d20_c3_homo, padj < 0.01, abs(log2FoldChange)>2)
deg_d20_c4_het <- subset(rstb_d20_c4_het, padj < 0.01, abs(log2FoldChange)>2)
deg_d20_c4_homo <- subset(rstb_d20_c4_homo, padj < 0.01, abs(log2FoldChange)>2)

deg_d60_c0_het <- subset(rstb_d60_c0_het, padj < 0.01, abs(log2FoldChange)>2)
deg_d60_c0_homo <- subset(rstb_d60_c0_homo, padj < 0.01, abs(log2FoldChange)>2)
deg_d60_c1_het <- subset(rstb_d60_c1_het, padj < 0.01, abs(log2FoldChange)>2)
deg_d60_c1_homo <- subset(rstb_d60_c1_homo, padj < 0.01, abs(log2FoldChange)>2)
deg_d60_c2_het <- subset(rstb_d60_c2_het, padj < 0.01, abs(log2FoldChange)>2)
deg_d60_c2_homo <- subset(rstb_d60_c2_homo, padj < 0.01, abs(log2FoldChange)>2)
deg_d60_c3_het <- subset(rstb_d60_c3_het, padj < 0.01, abs(log2FoldChange)>2)
deg_d60_c3_homo <- subset(rstb_d60_c3_homo, padj < 0.01, abs(log2FoldChange)>2)
deg_d60_c4_het <- subset(rstb_d60_c4_het, padj < 0.01, abs(log2FoldChange)>2)
deg_d60_c4_homo <- subset(rstb_d60_c4_homo, padj < 0.01, abs(log2FoldChange)>2)
deg_d60_c5_het <- subset(rstb_d60_c5_het, padj < 0.01, abs(log2FoldChange)>2)
deg_d60_c5_homo <- subset(rstb_d60_c5_homo, padj < 0.01, abs(log2FoldChange)>2)
deg_d60_c6_het <- subset(rstb_d60_c6_het, padj < 0.01, abs(log2FoldChange)>2)
deg_d60_c6_homo <- subset(rstb_d60_c6_homo, padj < 0.01, abs(log2FoldChange)>2)
deg_d60_c7_het <- subset(rstb_d60_c7_het, padj < 0.01, abs(log2FoldChange)>2)
deg_d60_c7_homo <- subset(rstb_d60_c7_homo, padj < 0.01, abs(log2FoldChange)>2)

# count the row numbers of deg tables and make a dataframe summering the deg numbers for each samples. denote day:d20 or d60, cluster:c0-c7, genotype: het or homo. 

het_degN_d20 <- c(nrow(deg_d20_c0_het), nrow(deg_d20_c1_het), nrow(deg_d20_c2_het), nrow(deg_d20_c3_het), nrow(deg_d20_c4_het))
ko_degN_d20 <- c(nrow(deg_d20_c0_homo), nrow(deg_d20_c1_homo), nrow(deg_d20_c2_homo), nrow(deg_d20_c3_homo), nrow(deg_d20_c4_homo))

het_degN_d60 <- c(nrow(deg_d60_c0_het), nrow(deg_d60_c1_het), nrow(deg_d60_c2_het), nrow(deg_d60_c3_het), nrow(deg_d60_c4_het), nrow(deg_d60_c5_het), nrow(deg_d60_c6_het), nrow(deg_d60_c7_het))
ko_degN_d60 <- c(nrow(deg_d60_c0_homo), nrow(deg_d60_c1_homo), nrow(deg_d60_c2_homo), nrow(deg_d60_c3_homo), nrow(deg_d60_c4_homo), nrow(deg_d60_c5_homo), nrow(deg_d60_c6_homo), nrow(deg_d60_c7_homo))

het_degN_d20
ko_degN_d20
het_degN_d60
ko_degN_d60

# make a table of deg numbers of d20 datasets
deg_d20 <- rbind(het_degN_d20, ko_degN_d20)
colnames(deg_d20) <- c("Cluster0", "Cluster1", "Cluster2", "Cluster3", "Cluster4") 
rownames(deg_d20) <- c("RAI1-Het", "RAI1-Homo") 
deg_d20

# make a table of deg number of d60 datasets
deg_d60 <- rbind(het_degN_d60, ko_degN_d60)
colnames(deg_d60) <- c("Cluster0", "Cluster1", "Cluster2", "Cluster3", "Cluster4", "Cluster5", "Cluster6", "Cluster7")
rownames(deg_d60) <- c("RAI1-Het", "RAI1-Homo")
deg_d60

# rename the clusters 
colnames(deg_d20) <- c("NPC", "Mesenchymal", "NPC-G2M", "NPC-S", "Neuron")
colnames(deg_d60) <- c("Neuron", "Forebrain Glia", "MGP", "NPC", "Astrocyte", "NPC-CC", "Apoptotic", "Immune")

deg_d20
deg_d60

# make a bar plots of deg_d20 and deg_d60 separately. het bars as light blue, homo bars as deep blue
library(ggplot2)

# Bar plot for d20
deg_d20_df <- as.data.frame(t(deg_d20))
deg_d20_df$Cluster <- rownames(deg_d20_df)
deg_d20_df <- reshape2::melt(deg_d20_df, id.vars = "Cluster")

deg_d20_df$Cluster <- factor(deg_d20_df$Cluster, levels = c("NPC-G2M", "NPC-S", "NPC", "Neuron", "Mesenchymal"))

d20_barpot <- ggplot(deg_d20_df, aes(x = Cluster, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("RAI1-Het" = "lightblue", "RAI1-Homo" = "blue")) +
    labs(title = "DEG Counts for d20", x = "Cluster", y = "DEG Count") +
    theme_minimal()

d20_barpot

# Bar plot for d60
deg_d60_df <- as.data.frame(t(deg_d60))
deg_d60_df$Cluster <- rownames(deg_d60_df)
deg_d60_df <- reshape2::melt(deg_d60_df, id.vars = "Cluster")

deg_d60_df$Cluster <- factor(deg_d60_df$Cluster, levels = c("NPC-CC", "NPC", "Neuron", "Forebrain Glia", "Astrocyte", "MGP","Immune", "Apoptotic"))

d60_barpot <- ggplot(deg_d60_df, aes(x = Cluster, y = value, fill = variable)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("RAI1-Het" = "lightblue", "RAI1-Homo" = "blue")) +
    labs(title = "DEG Counts for d60", x = "Cluster", y = "DEG Count") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

d60_barpot 

# align the d20_barplot and d60_barplot and save as one PDF file.
library(gridExtra)

pdf("results/scRNAseq/pseudobulk/deg_barplots.pdf", width = 10, height = 8)
grid.arrange(d20_barpot, d60_barpot, ncol = 1)
dev.off()


