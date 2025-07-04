# load DEseq2
library(DESeq2)
library(ggplot2)
library(reshape2)

# ...existing code for loading data...

# Separate DEGs into upregulated and downregulated
deg_d20_c0_het_up <- subset(rstb_d20_c0_het, padj < 0.01 & log2FoldChange > 2)
deg_d20_c0_het_down <- subset(rstb_d20_c0_het, padj < 0.01 & log2FoldChange < -2)
deg_d20_c0_homo_up <- subset(rstb_d20_c0_homo, padj < 0.01 & log2FoldChange > 2)
deg_d20_c0_homo_down <- subset(rstb_d20_c0_homo, padj < 0.01 & log2FoldChange < -2)

deg_d20_c1_het_up <- subset(rstb_d20_c1_het, padj < 0.01 & log2FoldChange > 2)
deg_d20_c1_het_down <- subset(rstb_d20_c1_het, padj < 0.01 & log2FoldChange < -2)
deg_d20_c1_homo_up <- subset(rstb_d20_c1_homo, padj < 0.01 & log2FoldChange > 2)
deg_d20_c1_homo_down <- subset(rstb_d20_c1_homo, padj < 0.01 & log2FoldChange < -2)

deg_d20_c2_het_up <- subset(rstb_d20_c2_het, padj < 0.01 & log2FoldChange > 2)
deg_d20_c2_het_down <- subset(rstb_d20_c2_het, padj < 0.01 & log2FoldChange < -2)
deg_d20_c2_homo_up <- subset(rstb_d20_c2_homo, padj < 0.01 & log2FoldChange > 2)
deg_d20_c2_homo_down <- subset(rstb_d20_c2_homo, padj < 0.01 & log2FoldChange < -2)

deg_d20_c3_het_up <- subset(rstb_d20_c3_het, padj < 0.01 & log2FoldChange > 2)
deg_d20_c3_het_down <- subset(rstb_d20_c3_het, padj < 0.01 & log2FoldChange < -2)
deg_d20_c3_homo_up <- subset(rstb_d20_c3_homo, padj < 0.01 & log2FoldChange > 2)
deg_d20_c3_homo_down <- subset(rstb_d20_c3_homo, padj < 0.01 & log2FoldChange < -2)

deg_d20_c4_het_up <- subset(rstb_d20_c4_het, padj < 0.01 & log2FoldChange > 2)
deg_d20_c4_het_down <- subset(rstb_d20_c4_het, padj < 0.01 & log2FoldChange < -2)
deg_d20_c4_homo_up <- subset(rstb_d20_c4_homo, padj < 0.01 & log2FoldChange > 2)
deg_d20_c4_homo_down <- subset(rstb_d20_c4_homo, padj < 0.01 & log2FoldChange < -2)

deg_d60_c0_het_up <- subset(rstb_d60_c0_het, padj < 0.01 & log2FoldChange > 2)
deg_d60_c0_het_down <- subset(rstb_d60_c0_het, padj < 0.01 & log2FoldChange < -2)
deg_d60_c0_homo_up <- subset(rstb_d60_c0_homo, padj < 0.01 & log2FoldChange > 2)
deg_d60_c0_homo_down <- subset(rstb_d60_c0_homo, padj < 0.01 & log2FoldChange < -2)

deg_d60_c1_het_up <- subset(rstb_d60_c1_het, padj < 0.01 & log2FoldChange > 2)
deg_d60_c1_het_down <- subset(rstb_d60_c1_het, padj < 0.01 & log2FoldChange < -2)
deg_d60_c1_homo_up <- subset(rstb_d60_c1_homo, padj < 0.01 & log2FoldChange > 2)
deg_d60_c1_homo_down <- subset(rstb_d60_c1_homo, padj < 0.01 & log2FoldChange < -2)

deg_d60_c2_het_up <- subset(rstb_d60_c2_het, padj < 0.01 & log2FoldChange > 2)
deg_d60_c2_het_down <- subset(rstb_d60_c2_het, padj < 0.01 & log2FoldChange < -2)
deg_d60_c2_homo_up <- subset(rstb_d60_c2_homo, padj < 0.01 & log2FoldChange > 2)
deg_d60_c2_homo_down <- subset(rstb_d60_c2_homo, padj < 0.01 & log2FoldChange < -2)

deg_d60_c3_het_up <- subset(rstb_d60_c3_het, padj < 0.01 & log2FoldChange > 2)
deg_d60_c3_het_down <- subset(rstb_d60_c3_het, padj < 0.01 & log2FoldChange < -2)
deg_d60_c3_homo_up <- subset(rstb_d60_c3_homo, padj < 0.01 & log2FoldChange > 2)
deg_d60_c3_homo_down <- subset(rstb_d60_c3_homo, padj < 0.01 & log2FoldChange < -2)

deg_d60_c4_het_up <- subset(rstb_d60_c4_het, padj < 0.01 & log2FoldChange > 2)
deg_d60_c4_het_down <- subset(rstb_d60_c4_het, padj < 0.01 & log2FoldChange < -2)
deg_d60_c4_homo_up <- subset(rstb_d60_c4_homo, padj < 0.01 & log2FoldChange > 2)
deg_d60_c4_homo_down <- subset(rstb_d60_c4_homo, padj < 0.01 & log2FoldChange < -2)

deg_d60_c5_het_up <- subset(rstb_d60_c5_het, padj < 0.01 & log2FoldChange > 2)
deg_d60_c5_het_down <- subset(rstb_d60_c5_het, padj < 0.01 & log2FoldChange < -2)
deg_d60_c5_homo_up <- subset(rstb_d60_c5_homo, padj < 0.01 & log2FoldChange > 2)
deg_d60_c5_homo_down <- subset(rstb_d60_c5_homo, padj < 0.01 & log2FoldChange < -2)

deg_d60_c6_het_up <- subset(rstb_d60_c6_het, padj < 0.01 & log2FoldChange > 2)
deg_d60_c6_het_down <- subset(rstb_d60_c6_het, padj < 0.01 & log2FoldChange < -2)
deg_d60_c6_homo_up <- subset(rstb_d60_c6_homo, padj < 0.01 & log2FoldChange > 2)
deg_d60_c6_homo_down <- subset(rstb_d60_c6_homo, padj < 0.01 & log2FoldChange < -2)

deg_d60_c7_het_up <- subset(rstb_d60_c7_het, padj < 0.01 & log2FoldChange > 2)
deg_d60_c7_het_down <- subset(rstb_d60_c7_het, padj < 0.01 & log2FoldChange < -2)
deg_d60_c7_homo_up <- subset(rstb_d60_c7_homo, padj < 0.01 & log2FoldChange > 2)
deg_d60_c7_homo_down <- subset(rstb_d60_c7_homo, padj < 0.01 & log2FoldChange < -2)

# Count the number of upregulated and downregulated genes
het_degN_d20_up <- c(nrow(deg_d20_c0_het_up), nrow(deg_d20_c1_het_up), nrow(deg_d20_c2_het_up), nrow(deg_d20_c3_het_up), nrow(deg_d20_c4_het_up))
het_degN_d20_down <- c(nrow(deg_d20_c0_het_down), nrow(deg_d20_c1_het_down), nrow(deg_d20_c2_het_down), nrow(deg_d20_c3_het_down), nrow(deg_d20_c4_het_down))
ko_degN_d20_up <- c(nrow(deg_d20_c0_homo_up), nrow(deg_d20_c1_homo_up), nrow(deg_d20_c2_homo_up), nrow(deg_d20_c3_homo_up), nrow(deg_d20_c4_homo_up))
ko_degN_d20_down <- c(nrow(deg_d20_c0_homo_down), nrow(deg_d20_c1_homo_down), nrow(deg_d20_c2_homo_down), nrow(deg_d20_c3_homo_down), nrow(deg_d20_c4_homo_down))

het_degN_d60_up <- c(nrow(deg_d60_c0_het_up), nrow(deg_d60_c1_het_up), nrow(deg_d60_c2_het_up), nrow(deg_d60_c3_het_up), nrow(deg_d60_c4_het_up), nrow(deg_d60_c5_het_up), nrow(deg_d60_c6_het_up), nrow(deg_d60_c7_het_up))
het_degN_d60_down <- c(nrow(deg_d60_c0_het_down), nrow(deg_d60_c1_het_down), nrow(deg_d60_c2_het_down), nrow(deg_d60_c3_het_down), nrow(deg_d60_c4_het_down), nrow(deg_d60_c5_het_down), nrow(deg_d60_c6_het_down), nrow(deg_d60_c7_het_down))
ko_degN_d60_up <- c(nrow(deg_d60_c0_homo_up), nrow(deg_d60_c1_homo_up), nrow(deg_d60_c2_homo_up), nrow(deg_d60_c3_homo_up), nrow(deg_d60_c4_homo_up), nrow(deg_d60_c5_homo_up), nrow(deg_d60_c6_homo_up), nrow(deg_d60_c7_homo_up))
ko_degN_d60_down <- c(nrow(deg_d60_c0_homo_down), nrow(deg_d60_c1_homo_down), nrow(deg_d60_c2_homo_down), nrow(deg_d60_c3_homo_down), nrow(deg_d60_c4_homo_down), nrow(deg_d60_c5_homo_down), nrow(deg_d60_c6_homo_down), nrow(deg_d60_c7_homo_down))

# Create data frames for plotting
deg_d20_up_down <- data.frame(
    Cluster = rep(c("NPC", "Mesenchymal", "NPC-G2M", "NPC-S", "Neuron"), each = 2),
    Type = rep(c("Upregulated", "Downregulated"), times = 5),
    RAI1_Het = c(het_degN_d20_up, het_degN_d20_down),
    RAI1_Homo = c(ko_degN_d20_up, ko_degN_d20_down)
)

deg_d60_up_down <- data.frame(
    Cluster = rep(c("Neuron", "Forebrain Glia", "MGP", "NPC", "Astrocyte", "NPC-CC", "Apoptotic", "Immune"), each = 2),
    Type = rep(c("Upregulated", "Downregulated"), times = 8),
    RAI1_Het = c(het_degN_d60_up, het_degN_d60_down),
    RAI1_Homo = c(ko_degN_d60_up, ko_degN_d60_down)
)

# Melt data frames for ggplot
deg_d20_up_down_melt <- melt(deg_d20_up_down, id.vars = c("Cluster", "Type"))
deg_d60_up_down_melt <- melt(deg_d60_up_down, id.vars = c("Cluster", "Type"))

deg_d20_up_down_melt

# Plot for d20
ggplot(deg_d20_up_down_melt, aes(x = Cluster, y = value, fill = interaction(Type, variable))) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c("Upregulated.RAI1_Het" = "lightblue", "Downregulated.RAI1_Het" = "blue", "Upregulated.RAI1_Homo" = "lightcoral", "Downregulated.RAI1_Homo" = "red")) +
    labs(title = "DEG Counts for d20", x = "Cluster", y = "DEG Count") +
    theme_minimal()

# Plot for d60
ggplot(deg_d60_up_down_melt, aes(x = Cluster, y = value, fill = interaction(Type, variable))) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c("Upregulated.RAI1_Het" = "lightblue", "Downregulated.RAI1_Het" = "blue", "Upregulated.RAI1_Homo" = "lightcoral", "Downregulated.RAI1_Homo" = "red")) +
    labs(title = "DEG Counts for d60", x = "Cluster", y = "DEG Count") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
