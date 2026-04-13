#Q: How many genes go up or down at 1 hr and 4 hrs? ####
# redo the analysis with het and homo separately. 
# this script generates PCA plot and bar plot of DEgene numbers. 

## intall DEseq2
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")

## Load DEseq2 ####
library(DESeq2)

#Load forestmngr to do rounding of numbers
library(dplyr)

#import gene count file
count <- read.delim("data/blueberry_output/combine_counts/gene_expected_count.txt")
colnames(count) <- gsub(".expected_count","",colnames(count))
colnames(count)

#set the Geneid as rownames
rownames(count) <- count$gene_id
count$gene_id <- NULL

#round numbers and make gene_id as rownames
count<- round(count)
View(count)

#Reorder the columns
colnames(count)
count_het <- count[, c(33:48, 1:16)] 
count_ko <- count[, c(33:48, 17:32)] 
colnames(count_het)
colnames(count_ko)

#remove the two het samples
count_het$Rai1.het_d34_rep1 <- NULL
count_het$Rai1.het_d60_rep4 <- NULL
colnames(count_het)

#generate experimental design file so called “se”. 
design_het <- data.frame(row.names = colnames(count_het), conditions = c(
"Rai1.wt_d0", "Rai1.wt_d0", "Rai1.wt_d0", "Rai1.wt_d0", "Rai1.wt_d20", "Rai1.wt_d20", "Rai1.wt_d20", "Rai1.wt_d20",
"Rai1.wt_d34", "Rai1.wt_d34", "Rai1.wt_d34", "Rai1.wt_d34",
"Rai1.wt_d60", "Rai1.wt_d60", "Rai1.wt_d60", "Rai1.wt_d60",
"Rai1.het_d0", "Rai1.het_d0", "Rai1.het_d0", "Rai1.het_d0",   "Rai1.het_d20", "Rai1.het_d20", "Rai1.het_d20", "Rai1.het_d20",
"Rai1.het_d34", "Rai1.het_d34", "Rai1.het_d34",
"Rai1.het_d60", "Rai1.het_d60", "Rai1.het_d60"))

design_ko <- data.frame(row.names = colnames(count_ko), conditions = c(
"Rai1.wt_d0", "Rai1.wt_d0", "Rai1.wt_d0", "Rai1.wt_d0",   "Rai1.wt_d20", "Rai1.wt_d20", "Rai1.wt_d20", "Rai1.wt_d20",
"Rai1.wt_d34", "Rai1.wt_d34", "Rai1.wt_d34", "Rai1.wt_d34",
"Rai1.wt_d60", "Rai1.wt_d60", "Rai1.wt_d60", "Rai1.wt_d60",
"Rai1.ko_d0", "Rai1.ko_d0", "Rai1.ko_d0", "Rai1.ko_d0",   "Rai1.ko_d20", "Rai1.ko_d20", "Rai1.ko_d20", "Rai1.ko_d20",
"Rai1.ko_d34", "Rai1.ko_d34", "Rai1.ko_d34", "Rai1.ko_d34",
"Rai1.ko_d60", "Rai1.ko_d60", "Rai1.ko_d60", "Rai1.ko_d60"))

#check the design
design_het
design_ko

# create a new column in the design data frame to indicate the genotype
design_het$genotype <- gsub("Rai1.([a-z]+)_d.*", "\\1", design_het$conditions)
design_ko$genotype <- gsub("Rai1.([a-z]+)_d.*", "\\1", design_ko$conditions)

design_het$day <- gsub("Rai1.[a-z]+_d([0-9]+)", "\\1", design_het$conditions)
design_ko$day <- gsub("Rai1.[a-z]+_d([0-9]+)", "\\1", design_ko$conditions)
head(design_het)
head(design_ko)

#Run DEseq2
dds_het <- DESeqDataSetFromMatrix(count_het, colData=design_het, design = ~ conditions) 
dds_ko <- DESeqDataSetFromMatrix(count_ko, colData=design_ko, design = ~ conditions)

dds_het_run <- DESeq(dds_het, betaPrior = TRUE)
dds_ko_run <- DESeq(dds_ko, betaPrior = TRUE)

# PCA plots
vsd_dds_het <- vst(dds_het, blind = TRUE)
head(assay(vsd_dds_het), 3)
vsd_dds_ko <- vst(dds_ko, blind = TRUE)
head(assay(vsd_dds_ko), 3)
str(vsd_dds_het)

library("ggplot2")

PCA_het <- plotPCA(vsd_dds_het, intgroup = c("genotype", "day"))
PCA_data_het <- plotPCA(vsd_dds_het, intgroup = c("genotype", "day"), returnData = TRUE)

PCA_ko <- plotPCA(vsd_dds_ko, intgroup = c("genotype", "day"))
PCA_data_ko <- plotPCA(vsd_dds_ko, intgroup = c("genotype", "day"), returnData = TRUE)

# create a shape and color mapping for the genotype and day
shape_mapping <- c("wt" = 1, "het" = 2, "ko" = 4)
color_mapping <- c("0" = "light blue", "20" = "cyan", "34" = "pink", "60" = "red")

# plot PCA
pca_het <- ggplot(PCA_data_het, aes(x = PC1, y = PC2, shape = genotype, color = day)) +
  geom_point(size = 4) +
  scale_shape_manual(values = shape_mapping) +
  scale_color_manual(values = color_mapping) +
  ggtitle("PCA Plot of RNA-seq Data (WT vs. Het)") +
  theme(text = element_text(family = "Helvetica", size = 10),
  axis.title = element_text(size = 10),
  axis.text = element_text(size = 8)) +
  theme_classic()
pca_het

pca_ko <- ggplot(PCA_data_ko, aes(x = PC1, y = PC2, shape = genotype, color = day)) +
  geom_point(size = 4) +
  scale_shape_manual(values = shape_mapping) +
  scale_color_manual(values = color_mapping) +
  ggtitle("PCA Plot of RNA-seq Data (WT vs. KO)") +
  theme(text = element_text(family = "Helvetica", size = 10),
  axis.title = element_text(size = 10),
  axis.text = element_text(size = 8)) +
  theme_classic()
pca_ko

# save the PCA plots
ggsave(filename = sprintf("results/deseq2/sep/%s_pca_het.pdf", Sys.Date()), plot = pca_het, width = 6, height = 4.5, units = "in")
ggsave(filename = sprintf("results/deseq2/sep/%s_pca_ko.pdf", Sys.Date()), plot = pca_ko, width = 6, height = 4.5, units = "in")

sprintf("results/cluster/%s_dp_kegg.pdf", Sys.Date()

# read the restable \
rstbl_d0_wt.vs.het <- read.delim("results/deseq2/rstbl_d0_wt.vs.het.txt")
rstbl_d0_wt.vs.ko <- read.delim("results/deseq2/rstbl_d0_wt.vs.ko.txt")
rstbl_d20_wt.vs.het <- read.delim("results/deseq2/rstbl_d20_wt.vs.het.txt")
rstbl_d20_wt.vs.ko <- read.delim("results/deseq2/rstbl_d20_wt.vs.ko.txt")
rstbl_d34_wt.vs.het <- read.delim("results/deseq2/rstbl_d34_wt.vs.het.txt")
rstbl_d34_wt.vs.ko <- read.delim("results/deseq2/rstbl_d34_wt.vs.ko.txt")
rstbl_d60_wt.vs.het <- read.delim("results/deseq2/rstbl_d60_wt.vs.het.txt")
rstbl_d60_wt.vs.ko <- read.delim("results/deseq2/rstbl_d60_wt.vs.ko.txt")

# select degs 
up_d0_wt.vs.het <- subset(rstbl_d0_wt.vs.het, padj < 0.05 & log2FoldChange > 0)
up_d0_wt.vs.ko <- subset(rstbl_d0_wt.vs.ko, padj < 0.05 & log2FoldChange > 0)
up_d20_wt.vs.het <- subset(rstbl_d20_wt.vs.het, padj < 0.05 & log2FoldChange > 0)
up_d20_wt.vs.ko <- subset(rstbl_d20_wt.vs.ko, padj < 0.05 & log2FoldChange > 0)
up_d34_wt.vs.het <- subset(rstbl_d34_wt.vs.het, padj < 0.05 & log2FoldChange > 0)
up_d34_wt.vs.ko <- subset(rstbl_d34_wt.vs.ko, padj < 0.05 & log2FoldChange > 0)
up_d60_wt.vs.het <- subset(rstbl_d60_wt.vs.het, padj < 0.05 & log2FoldChange > 0)
up_d60_wt.vs.ko <- subset(rstbl_d60_wt.vs.ko, padj < 0.05 & log2FoldChange > 0)

down_d0_wt.vs.het <- subset(rstbl_d0_wt.vs.het, padj < 0.05 & log2FoldChange < 0)
down_d0_wt.vs.ko <- subset(rstbl_d0_wt.vs.ko, padj < 0.05 & log2FoldChange < 0)
down_d20_wt.vs.het <- subset(rstbl_d20_wt.vs.het, padj < 0.05 & log2FoldChange < 0)
down_d20_wt.vs.ko <- subset(rstbl_d20_wt.vs.ko, padj < 0.05 & log2FoldChange < 0)
down_d34_wt.vs.het <- subset(rstbl_d34_wt.vs.het, padj < 0.05 & log2FoldChange < 0)
down_d34_wt.vs.ko <- subset(rstbl_d34_wt.vs.ko, padj < 0.05 & log2FoldChange < 0)
down_d60_wt.vs.het <- subset(rstbl_d60_wt.vs.het, padj < 0.05 & log2FoldChange < 0)
down_d60_wt.vs.ko <- subset(rstbl_d60_wt.vs.ko, padj < 0.05 & log2FoldChange < 0)

# plot deg numbers
## making data frame for the deg number plot
condition <- c("up_d0_wt.vs.het", "up_d0_wt.vs.ko", "up_d20_wt.vs.het", "up_d20_wt.vs.ko", "up_d34_wt.vs.het", "up_d34_wt.vs.ko", "up_d60_wt.vs.het", "up_d60_wt.vs.ko", "down_d0_wt.vs.het", "down_d0_wt.vs.ko", "down_d20_wt.vs.het", 
"down_d20_wt.vs.ko", "down_d34_wt.vs.het", "down_d34_wt.vs.ko", 
"down_d60_wt.vs.het", "down_d60_wt.vs.ko")

direction <- c("Up", "Up", "Up", "Up", "Up", "Up", "Up", "Up", "Down", "Down", "Down", "Down", "Down", "Down", "Down", "Down")

time <- c("d0", "d0", "d20", "d20", "d34", "d34", "d60", "d60", "d0", "d0", "d20", "d20", "d34", "d34", "d60", "d60")

genotype <- c("Heterozygous", "Homozygous", "Heterozygous", "Homozygous", "Heterozygous", "Homozygous", "Heterozygous", "Homozygous", "Heterozygous", "Homozygous", "Heterozygous", "Homozygous", "Heterozygous", "Homozygous", "Heterozygous", "Homozygous")

gene_n <- c(nrow(up_d0_wt.vs.het), nrow(up_d0_wt.vs.ko),
nrow(up_d20_wt.vs.het), nrow(up_d20_wt.vs.ko),
nrow(up_d34_wt.vs.het), nrow(up_d34_wt.vs.ko),
nrow(up_d60_wt.vs.het), nrow(up_d60_wt.vs.ko),
nrow(down_d0_wt.vs.het), nrow(down_d0_wt.vs.ko),
nrow(down_d20_wt.vs.het), nrow(down_d20_wt.vs.ko),
nrow(down_d34_wt.vs.het), nrow(down_d34_wt.vs.ko),
nrow(down_d60_wt.vs.het), nrow(down_d60_wt.vs.ko))

degene_n <- data.frame(row.names = condition, direction = direction, genotype = genotype, day = time, gene_n = gene_n)
degene_n

# plot deg numbers
library(ggplot2)

# Filter the dataframe to include only "wt" and "het" genotypes
filtered_df_het <- degene_n[degene_n$genotype %in% c("Heterozygous"), ]

# Create the bar plot
bar_het <- ggplot(filtered_df_het, aes(x = day, y = gene_n, fill = direction)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ genotype) +
  scale_y_continuous(labels = function(x) format(x, big.mark = ",")) +
  labs(x = "Days in culture", y = "Number of genes", fill = "Direction") +
  theme(text = element_text(family = "Helvetica", size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

bar_het

filtered_df_ko <- degene_n[degene_n$genotype %in% c("Homozygous"), ]

bar_ko <- ggplot(filtered_df_ko, aes(x = day, y = gene_n, fill = direction)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ genotype) +
  scale_y_continuous(labels = function(x) format(x, big.mark = ",")) +
  labs(x = "Days in culture", y = "Number of genes", fill = "Direction") +
  theme(text = element_text(family = "Helvetica", size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

bar_ko

# save the Bar plots
ggsave(filename = sprintf("results/deseq2/sep/%s_deg_n_het.pdf", Sys.Date()), plot = bar_het, width = 4, height = 4.5, units = "in")
ggsave(filename = sprintf("results/deseq2/sep/%s_deg_n_ko.pdf", Sys.Date()), plot = bar_ko, width = 4, height = 4.5, units = "in")

# select degs with more stringent cutoff padj < 0.01)
up_d0_wt.vs.het_2 <- subset(rstbl_d0_wt.vs.het, padj < 0.01 & log2FoldChange > 0)
up_d0_wt.vs.ko_2 <- subset(rstbl_d0_wt.vs.ko, padj < 0.01 & log2FoldChange > 0)
up_d20_wt.vs.het_2 <- subset(rstbl_d20_wt.vs.het, padj < 0.01 & log2FoldChange > 0)
up_d20_wt.vs.ko_2 <- subset(rstbl_d20_wt.vs.ko, padj < 0.01 & log2FoldChange > 0)
up_d34_wt.vs.het_2 <- subset(rstbl_d34_wt.vs.het, padj < 0.01 & log2FoldChange > 0)
up_d34_wt.vs.ko_2 <- subset(rstbl_d34_wt.vs.ko, padj < 0.01 & log2FoldChange > 0)
up_d60_wt.vs.het_2 <- subset(rstbl_d60_wt.vs.het, padj < 0.01 & log2FoldChange > 0)
up_d60_wt.vs.ko_2 <- subset(rstbl_d60_wt.vs.ko, padj < 0.01 & log2FoldChange > 0)

down_d0_wt.vs.het_2 <- subset(rstbl_d0_wt.vs.het, padj < 0.01 & log2FoldChange < 0)
down_d0_wt.vs.ko_2 <- subset(rstbl_d0_wt.vs.ko, padj < 0.01 & log2FoldChange < 0)
down_d20_wt.vs.het_2 <- subset(rstbl_d20_wt.vs.het, padj < 0.01 & log2FoldChange < 0)
down_d20_wt.vs.ko_2 <- subset(rstbl_d20_wt.vs.ko, padj < 0.01 & log2FoldChange < 0)
down_d34_wt.vs.het_2 <- subset(rstbl_d34_wt.vs.het, padj < 0.01 & log2FoldChange < 0)
down_d34_wt.vs.ko_2 <- subset(rstbl_d34_wt.vs.ko, padj < 0.01 & log2FoldChange < 0)
down_d60_wt.vs.het_2 <- subset(rstbl_d60_wt.vs.het, padj < 0.01 & log2FoldChange < 0)
down_d60_wt.vs.ko_2 <- subset(rstbl_d60_wt.vs.ko, padj < 0.01 & log2FoldChange < 0)

# generate data frame for the deg number plot with more stringent cutoff
gene_n2 <- c(nrow(up_d0_wt.vs.het_2), nrow(up_d0_wt.vs.ko_2),
nrow(up_d20_wt.vs.het_2), nrow(up_d20_wt.vs.ko_2),
nrow(up_d34_wt.vs.het_2), nrow(up_d34_wt.vs.ko_2),
nrow(up_d60_wt.vs.het_2), nrow(up_d60_wt.vs.ko_2),
nrow(down_d0_wt.vs.het_2), nrow(down_d0_wt.vs.ko_2),
nrow(down_d20_wt.vs.het_2), nrow(down_d20_wt.vs.ko_2),
nrow(down_d34_wt.vs.het_2), nrow(down_d34_wt.vs.ko_2),
nrow(down_d60_wt.vs.het_2), nrow(down_d60_wt.vs.ko_2))

degene_n2 <- data.frame(row.names = condition, direction = direction, genotype = genotype, day = time, gene_n = gene_n2)
degene_n2

# plot deg numbers
library(ggplot2)

# Filter the dataframe to include only "wt" and "het" genotypes
filtered_df_het2 <- degene_n2[degene_n2$genotype %in% c("Heterozygous"), ]

# Create the bar plot of the deg numbers with more stringent cutoff
bar_het2 <- ggplot(filtered_df_het2, aes(x = day, y = gene_n, fill = direction)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ genotype) +
  scale_y_continuous(labels = function(x) format(x, big.mark = ",")) +
  labs(x = "Days in culture", y = "Number of genes", fill = "Direction") +
  theme(text = element_text(family = "Helvetica", size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

bar_het2

filtered_df_ko2 <- degene_n2[degene_n2$genotype %in% c("Homozygous"), ]

bar_ko2 <- ggplot(filtered_df_ko2, aes(x = day, y = gene_n, fill = direction)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ genotype) +
  scale_y_continuous(labels = function(x) format(x, big.mark = ",")) +
  labs(x = "Days in culture", y = "Number of genes", fill = "Direction") +
  theme(text = element_text(family = "Helvetica", size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

bar_ko2

# save the Bar plots
ggsave(filename = sprintf("results/deseq2/sep/%s_deg_n_het_2.pdf", Sys.Date()), plot = bar_het2, width = 4, height = 4.5, units = "in")
ggsave(filename = sprintf("results/deseq2/sep/%s_deg_n_ko_2.pdf", Sys.Date()), plot = bar_ko2, width = 4, height = 4.5, units = "in")

# further impose fold change cutoff (|log2fc| > 2)
# select degs with more stringent cutoff (|log2fc| > 2)
up_d0_wt.vs.het_3 <- subset(rstbl_d0_wt.vs.het, padj < 0.01 & log2FoldChange > 2)
up_d0_wt.vs.ko_3 <- subset(rstbl_d0_wt.vs.ko, padj < 0.01 & log2FoldChange > 2)
up_d20_wt.vs.het_3 <- subset(rstbl_d20_wt.vs.het, padj < 0.01 & log2FoldChange > 2)
up_d20_wt.vs.ko_3 <- subset(rstbl_d20_wt.vs.ko, padj < 0.01 & log2FoldChange > 2)
up_d34_wt.vs.het_3 <- subset(rstbl_d34_wt.vs.het, padj < 0.01 & log2FoldChange > 2)
up_d34_wt.vs.ko_3 <- subset(rstbl_d34_wt.vs.ko, padj < 0.01 & log2FoldChange > 2)
up_d60_wt.vs.het_3 <- subset(rstbl_d60_wt.vs.het, padj < 0.01 & log2FoldChange > 2)
up_d60_wt.vs.ko_3 <- subset(rstbl_d60_wt.vs.ko, padj < 0.01 & log2FoldChange > 2)

down_d0_wt.vs.het_3 <- subset(rstbl_d0_wt.vs.het, padj < 0.01 & log2FoldChange < -2)
down_d0_wt.vs.ko_3 <- subset(rstbl_d0_wt.vs.ko, padj < 0.01 & log2FoldChange < -2)
down_d20_wt.vs.het_3 <- subset(rstbl_d20_wt.vs.het, padj < 0.01 & log2FoldChange < -2)
down_d20_wt.vs.ko_3 <- subset(rstbl_d20_wt.vs.ko, padj < 0.01 & log2FoldChange < -2)
down_d34_wt.vs.het_3 <- subset(rstbl_d34_wt.vs.het, padj < 0.01 & log2FoldChange < -2)
down_d34_wt.vs.ko_3 <- subset(rstbl_d34_wt.vs.ko, padj < 0.01 & log2FoldChange < -2)
down_d60_wt.vs.het_3 <- subset(rstbl_d60_wt.vs.het, padj < 0.01 & log2FoldChange < -2)
down_d60_wt.vs.ko_3 <- subset(rstbl_d60_wt.vs.ko, padj < 0.01 & log2FoldChange < -2)


# generate data frame for the deg number plot with more stringent cutoff
gene_n3 <- c(nrow(up_d0_wt.vs.het_3), nrow(up_d0_wt.vs.ko_3),
nrow(up_d20_wt.vs.het_3), nrow(up_d20_wt.vs.ko_3),
nrow(up_d34_wt.vs.het_3), nrow(up_d34_wt.vs.ko_3),
nrow(up_d60_wt.vs.het_3), nrow(up_d60_wt.vs.ko_3),
nrow(down_d0_wt.vs.het_3), nrow(down_d0_wt.vs.ko_3),
nrow(down_d20_wt.vs.het_3), nrow(down_d20_wt.vs.ko_3),
nrow(down_d34_wt.vs.het_3), nrow(down_d34_wt.vs.ko_3),
nrow(down_d60_wt.vs.het_3), nrow(down_d60_wt.vs.ko_3))

degene_n3 <- data.frame(row.names = condition, direction = direction, genotype = genotype, day = time, gene_n = gene_n3)
degene_n3

# plot deg numbers
library(ggplot2)

# Filter the dataframe to include only "wt" and "het" genotypes
filtered_df_het3 <- degene_n3[degene_n3$genotype %in% c("Heterozygous"), ]

# Create the bar plot of the deg numbers with more stringent cutoff
bar_het3 <- ggplot(filtered_df_het3, aes(x = day, y = gene_n, fill = direction)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ genotype) +
  scale_y_continuous(labels = function(x) format(x, big.mark = ",")) +
  labs(x = "Days in culture", y = "Number of genes", fill = "Direction") +
  theme(text = element_text(family = "Helvetica", size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

bar_het3

filtered_df_ko3 <- degene_n3[degene_n3$genotype %in% c("Homozygous"), ]

bar_ko3 <- ggplot(filtered_df_ko3, aes(x = day, y = gene_n, fill = direction)) +
  geom_bar(stat = "identity") +
  facet_grid(. ~ genotype) +
  scale_y_continuous(labels = function(x) format(x, big.mark = ",")) +
  labs(x = "Days in culture", y = "Number of genes", fill = "Direction") +
  theme(text = element_text(family = "Helvetica", size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8))

bar_ko3

# save the Bar plots
ggsave(filename = sprintf("results/deseq2/sep/%s_deg_n_het_3.pdf", Sys.Date()), plot = bar_het3, width = 4, height = 4.5, units = "in")
ggsave(filename = sprintf("results/deseq2/sep/%s_deg_n_ko_3.pdf", Sys.Date()), plot = bar_ko3, width = 4, height = 4.5, units = "in")
