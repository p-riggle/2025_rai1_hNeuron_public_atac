# ===============================
# Load required libraries
# ===============================
setwd("/Users/satabdi/Downloads")
library(ggplot2)
library(tidyr)

# Desired cell order
cell_order <- c(
  "NPC-G2/M",
  "NPC-G1",
  "Mes-WT",
  "Mes-KO",
  "IPC",
  "Neuron",
  "gRG"
)

# ===============================
# Read CSVs
# ===============================

# HET vs WT
rstbl_NPCG2M_wt.vs.het <- read.csv("rstbl_NPCG2M_wt.vs.het.csv")
rstbl_NPCG1_wt.vs.het  <- read.csv("rstbl_NPCG1_wt.vs.het.csv")
rstbl_MesWT_wt.vs.het  <- read.csv("rstbl_MesWT_wt.vs.het.csv")
rstbl_MesKO_wt.vs.het  <- read.csv("rstbl_MesKO_wt.vs.het.csv")
rstbl_IPC_wt.vs.het    <- read.csv("rstbl_IPC_wt.vs.het.csv")
rstbl_Neuron_wt.vs.het <- read.csv("rstbl_Neuron_wt.vs.het.csv")
rstbl_Astro_wt.vs.het  <- read.csv("rstbl_Astro_wt.vs.het.csv")

# KO vs WT
rstbl_NPCG2M_wt.vs.ko <- read.csv("rstbl_NPCG2M_wt.vs.ko.csv")
rstbl_NPCG1_wt.vs.ko  <- read.csv("rstbl_NPCG1_wt.vs.ko.csv")
rstbl_MesWT_wt.vs.ko  <- read.csv("rstbl_MesWT_wt.vs.ko.csv")
rstbl_MesKO_wt.vs.ko  <- read.csv("rstbl_MesKO_wt.vs.ko.csv")
rstbl_IPC_wt.vs.ko    <- read.csv("rstbl_IPC_wt.vs.ko.csv")
rstbl_Neuron_wt.vs.ko <- read.csv("rstbl_Neuron_wt.vs.ko.csv")
rstbl_Astro_wt.vs.ko  <- read.csv("rstbl_Astro_wt.vs.ko.csv")

# ===============================
# Count DEGs HET vs WT
# ===============================
deg_counts_het <- data.frame(
  Cell_Type = cell_order,
  
  Upregulated = c(
    sum(rstbl_NPCG2M_wt.vs.het$padj < 0.01 & rstbl_NPCG2M_wt.vs.het$log2FoldChange > 1, na.rm=TRUE),
    sum(rstbl_NPCG1_wt.vs.het$padj < 0.01 & rstbl_NPCG1_wt.vs.het$log2FoldChange > 1, na.rm=TRUE),
    sum(rstbl_MesWT_wt.vs.het$padj < 0.01 & rstbl_MesWT_wt.vs.het$log2FoldChange > 1, na.rm=TRUE),
    sum(rstbl_MesKO_wt.vs.het$padj < 0.01 & rstbl_MesKO_wt.vs.het$log2FoldChange > 1, na.rm=TRUE),
    sum(rstbl_IPC_wt.vs.het$padj < 0.01 & rstbl_IPC_wt.vs.het$log2FoldChange > 1, na.rm=TRUE),
    sum(rstbl_Neuron_wt.vs.het$padj < 0.01 & rstbl_Neuron_wt.vs.het$log2FoldChange > 1, na.rm=TRUE),
    sum(rstbl_Astro_wt.vs.het$padj < 0.01 & rstbl_Astro_wt.vs.het$log2FoldChange > 1, na.rm=TRUE)
  ),
  
  Downregulated = c(
    sum(rstbl_NPCG2M_wt.vs.het$padj < 0.01 & rstbl_NPCG2M_wt.vs.het$log2FoldChange < -1, na.rm=TRUE),
    sum(rstbl_NPCG1_wt.vs.het$padj < 0.01 & rstbl_NPCG1_wt.vs.het$log2FoldChange < -1, na.rm=TRUE),
    sum(rstbl_MesWT_wt.vs.het$padj < 0.01 & rstbl_MesWT_wt.vs.het$log2FoldChange < -1, na.rm=TRUE),
    sum(rstbl_MesKO_wt.vs.het$padj < 0.01 & rstbl_MesKO_wt.vs.het$log2FoldChange < -1, na.rm=TRUE),
    sum(rstbl_IPC_wt.vs.het$padj < 0.01 & rstbl_IPC_wt.vs.het$log2FoldChange < -1, na.rm=TRUE),
    sum(rstbl_Neuron_wt.vs.het$padj < 0.01 & rstbl_Neuron_wt.vs.het$log2FoldChange < -1, na.rm=TRUE),
    sum(rstbl_Astro_wt.vs.het$padj < 0.01 & rstbl_Astro_wt.vs.het$log2FoldChange < -1, na.rm=TRUE)
  )
)

deg_long_het <- pivot_longer(
  deg_counts_het,
  cols=c("Upregulated","Downregulated"),
  names_to="Regulation",
  values_to="Count"
)

deg_long_het$Cell_Type <- factor(deg_long_het$Cell_Type, levels=cell_order)

# ===============================
# Plot HET vs WT
# ===============================
p_het <- ggplot(deg_long_het, aes(x=Cell_Type, y=Count, fill=Regulation)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_manual(values=c("Upregulated"="yellow","Downregulated"="blue")) +
  labs(title="Number of DEGs (HET vs WT)", x="Cell Type", y="Number of DEGs") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1))

p_het

ggsave("DEG_barplot_HET_vs_WT.pdf", p_het, width=8, height=5)

# ===============================
# Count DEGs KO vs WT
# ===============================
deg_counts_ko <- data.frame(
  Cell_Type = cell_order,
  
  Upregulated = c(
    sum(rstbl_NPCG2M_wt.vs.ko$padj < 0.01 & rstbl_NPCG2M_wt.vs.ko$log2FoldChange > 1, na.rm=TRUE),
    sum(rstbl_NPCG1_wt.vs.ko$padj < 0.01 & rstbl_NPCG1_wt.vs.ko$log2FoldChange > 1, na.rm=TRUE),
    sum(rstbl_MesWT_wt.vs.ko$padj < 0.01 & rstbl_MesWT_wt.vs.ko$log2FoldChange > 1, na.rm=TRUE),
    sum(rstbl_MesKO_wt.vs.ko$padj < 0.01 & rstbl_MesKO_wt.vs.ko$log2FoldChange > 1, na.rm=TRUE),
    sum(rstbl_IPC_wt.vs.ko$padj < 0.01 & rstbl_IPC_wt.vs.ko$log2FoldChange > 1, na.rm=TRUE),
    sum(rstbl_Neuron_wt.vs.ko$padj < 0.01 & rstbl_Neuron_wt.vs.ko$log2FoldChange > 1, na.rm=TRUE),
    sum(rstbl_Astro_wt.vs.ko$padj < 0.01 & rstbl_Astro_wt.vs.ko$log2FoldChange > 1, na.rm=TRUE)
  ),
  
  Downregulated = c(
    sum(rstbl_NPCG2M_wt.vs.ko$padj < 0.01 & rstbl_NPCG2M_wt.vs.ko$log2FoldChange < -1, na.rm=TRUE),
    sum(rstbl_NPCG1_wt.vs.ko$padj < 0.01 & rstbl_NPCG1_wt.vs.ko$log2FoldChange < -1, na.rm=TRUE),
    sum(rstbl_MesWT_wt.vs.ko$padj < 0.01 & rstbl_MesWT_wt.vs.ko$log2FoldChange < -1, na.rm=TRUE),
    sum(rstbl_MesKO_wt.vs.ko$padj < 0.01 & rstbl_MesKO_wt.vs.ko$log2FoldChange < -1, na.rm=TRUE),
    sum(rstbl_IPC_wt.vs.ko$padj < 0.01 & rstbl_IPC_wt.vs.ko$log2FoldChange < -1, na.rm=TRUE),
    sum(rstbl_Neuron_wt.vs.ko$padj < 0.01 & rstbl_Neuron_wt.vs.ko$log2FoldChange < -1, na.rm=TRUE),
    sum(rstbl_Astro_wt.vs.ko$padj < 0.01 & rstbl_Astro_wt.vs.ko$log2FoldChange < -1, na.rm=TRUE)
  )
)

deg_long_ko <- pivot_longer(
  deg_counts_ko,
  cols=c("Upregulated","Downregulated"),
  names_to="Regulation",
  values_to="Count"
)

deg_long_ko$Cell_Type <- factor(deg_long_ko$Cell_Type, levels=cell_order)

# ===============================
# Plot KO vs WT
# ===============================
p_ko <- ggplot(deg_long_ko, aes(x=Cell_Type, y=Count, fill=Regulation)) +
  geom_bar(stat="identity", position="dodge") +
  scale_fill_manual(values=c("Upregulated"="yellow","Downregulated"="blue")) +
  labs(title="Number of DEGs (KO vs WT)", x="Cell Type", y="Number of DEGs") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=45, hjust=1))

p_ko

ggsave("DEG_barplot_KO_vs_WT.pdf", p_ko, width=8, height=5)
