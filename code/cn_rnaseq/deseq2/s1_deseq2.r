#Q: How many genes go up or down at 1 hr and 4 hrs? ####
# redo the analysis without het.d34.rep1 and het.d60.rep4
# this script generates PCA plot and bar plot of DEgene numbers. 

## Load DEseq2 ####
if (!require("BiocManager", quietly = TRUE))
   install.packages("BiocManager")

BiocManager::install("DESeq2")
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
count <- count[, c(33:48, 1:16, 17:32)] 
colnames(count)

#remove the two samples
count$Rai1.het_d34_rep1 <- NULL
count$Rai1.het_d60_rep4 <- NULL
colnames(count)

#generate experimental design file so called “se”. 
design <- data.frame(row.names = colnames(count), conditions = c(
"Rai1.wt_d0", "Rai1.wt_d0", "Rai1.wt_d0", "Rai1.wt_d0",   "Rai1.wt_d20", "Rai1.wt_d20", "Rai1.wt_d20", "Rai1.wt_d20",
"Rai1.wt_d34", "Rai1.wt_d34", "Rai1.wt_d34", "Rai1.wt_d34",
"Rai1.wt_d60", "Rai1.wt_d60", "Rai1.wt_d60", "Rai1.wt_d60",
"Rai1.het_d0", "Rai1.het_d0", "Rai1.het_d0", "Rai1.het_d0",   "Rai1.het_d20", "Rai1.het_d20", "Rai1.het_d20", "Rai1.het_d20",
"Rai1.het_d34", "Rai1.het_d34", "Rai1.het_d34",
"Rai1.het_d60", "Rai1.het_d60", "Rai1.het_d60",
"Rai1.ko_d0", "Rai1.ko_d0", "Rai1.ko_d0", "Rai1.ko_d0",   "Rai1.ko_d20", "Rai1.ko_d20", "Rai1.ko_d20", "Rai1.ko_d20",
"Rai1.ko_d34", "Rai1.ko_d34", "Rai1.ko_d34", "Rai1.ko_d34",
"Rai1.ko_d60", "Rai1.ko_d60", "Rai1.ko_d60", "Rai1.ko_d60"))

#check the design
design

#Run DEseq2 
# create a new column in the design data frame to indicate the genotype
design$genotype <- gsub("Rai1.([a-z]+)_d.*", "\\1", design$conditions)
# create a new column in the design data frame to indicate the day
design$day <- gsub("Rai1.[a-z]+_d([0-9]+)", "\\1", design$conditions)
head(design)

#Run DEseq2 
dds <- DESeqDataSetFromMatrix(count, colData=design, design = ~ conditions)
dds_run <- DESeq(dds, betaPrior = TRUE)
resultsNames(dds_run)

##PCAplots ####
vsd_dds <- vst(dds, blind = TRUE)  
head(assay(vsd_dds), 3)
#install.packages("ggplot2")
library("ggplot2")
str(vsd_dds)

PCA <- plotPCA(vsd_dds, intgroup = c("genotype", "day"))
PCA_data <- plotPCA(vsd_dds, intgroup = c("genotype", "day"), returnData = TRUE)

# create a shape and color mapping for the genotype and day
shape_mapping <- c("wt" = 1, "het" = 2, "ko" = 4)
color_mapping <- c("0" = "light blue", "20" = "cyan", "34" = "pink", "60" = "red")

# create a new plot with the modified aesthetics
pca2 <- ggplot(PCA_data, aes(x = PC1, y = PC2, shape = genotype, color = day)) +
  geom_point(size = 4) +
  scale_shape_manual(values = shape_mapping) +
  scale_color_manual(values = color_mapping) +
  ggtitle("PCA Plot of RNA-seq Data") +
  theme(text = element_text(family = "Helvetica", size = 10),
  axis.title = element_text(size = 10),
  axis.text = element_text(size = 8)) +
  theme_classic()

pca2

ggsave(filename = "results/deseq2/pca2.pdf", plot = pca2, width = 6, height = 4.5, units = "in")

#Identifying Bic-induced genes in DMSO ####
rstbl_d0_wt.vs.het <- results(dds_run, contrast = c("conditions", "Rai1.het_d0", "Rai1.wt_d0"), alpha = 0.05)
rstbl_d0_wt.vs.ko <- results(dds_run, contrast = c("conditions", "Rai1.ko_d0", "Rai1.wt_d0"), alpha = 0.05)
rstbl_d20_wt.vs.het <- results(dds_run, contrast = c("conditions", "Rai1.het_d20", "Rai1.wt_d20"), alpha = 0.05)
rstbl_d20_wt.vs.ko <- results(dds_run, contrast = c("conditions", "Rai1.ko_d20", "Rai1.wt_d20"), alpha = 0.05)
rstbl_d34_wt.vs.het <- results(dds_run, contrast = c("conditions", "Rai1.het_d34", "Rai1.wt_d34"), alpha = 0.05)
rstbl_d34_wt.vs.ko <- results(dds_run, contrast = c("conditions", "Rai1.ko_d34", "Rai1.wt_d34"), alpha = 0.05)
rstbl_d60_wt.vs.het <- results(dds_run, contrast = c("conditions", "Rai1.het_d60", "Rai1.wt_d60"), alpha = 0.05)
rstbl_d60_wt.vs.ko <- results(dds_run, contrast = c("conditions", "Rai1.ko_d60", "Rai1.wt_d60"), alpha = 0.05)

summary(rstbl_d0_wt.vs.het)
summary(rstbl_d0_wt.vs.ko)
summary(rstbl_d20_wt.vs.het)
summary(rstbl_d20_wt.vs.ko)
summary(rstbl_d34_wt.vs.het)
summary(rstbl_d34_wt.vs.ko)
summary(rstbl_d60_wt.vs.het)
summary(rstbl_d60_wt.vs.ko)

#Add gene symbol and entrezID and print result tables
library(AnnotationDbi)

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("org.Hs.eg.db")

library(org.Hs.eg.db)
columns(org.Hs.eg.db) 

#print the keytypes 

## Write restables for Bic-regulated genes

add_symbol <- function(restable, filename) {
  restable$entrezid <- mapIds(org.Hs.eg.db, keys=row.names(restable), column="ENTREZID", keytype="SYMBOL", multiVals="first")
  write.table(restable, file = filename, sep = "\t")
}

restable_list <- list(rstbl_d0_wt.vs.het, rstbl_d0_wt.vs.ko, rstbl_d20_wt.vs.het, rstbl_d20_wt.vs.ko, rstbl_d34_wt.vs.het, rstbl_d34_wt.vs.ko, rstbl_d60_wt.vs.het, rstbl_d60_wt.vs.ko)

restable_names <- c("rstbl_d0_wt.vs.het", "rstbl_d0_wt.vs.ko", "rstbl_d20_wt.vs.het", "rstbl_d20_wt.vs.ko", "rstbl_d34_wt.vs.het", "rstbl_d34_wt.vs.ko", "rstbl_d60_wt.vs.het", "rstbl_d60_wt.vs.ko")

for (i in 1:length(restable_list)) {
  restable <- restable_list[[i]]
  filename <- paste0("results/deseq2/", restable_names[i], ".txt")
  add_symbol(restable, filename)
}

# read the restable
rstbl_d0_wt.vs.het <- read.delim("results/deseq2/rstbl_d0_wt.vs.het.txt")
rstbl_d0_wt.vs.ko <- read.delim("results/deseq2/rstbl_d0_wt.vs.ko.txt")
rstbl_d20_wt.vs.het <- read.delim("results/deseq2/rstbl_d20_wt.vs.het.txt")
rstbl_d20_wt.vs.ko <- read.delim("results/deseq2/rstbl_d20_wt.vs.ko.txt")
rstbl_d34_wt.vs.het <- read.delim("results/deseq2/rstbl_d34_wt.vs.het.txt")
rstbl_d34_wt.vs.ko <- read.delim("results/deseq2/rstbl_d34_wt.vs.ko.txt")
rstbl_d60_wt.vs.het <- read.delim("results/deseq2/rstbl_d60_wt.vs.het.txt")
rstbl_d60_wt.vs.ko <- read.delim("results/deseq2/rstbl_d60_wt.vs.ko.txt")

head(rstbl_d0_wt.vs.het)

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

# plot DEgene numbers
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

library(ggplot2)
bar1 <- ggplot(degene_n, aes(x = day, y = gene_n, fill = direction)) +
  geom_bar(stat = "identity") +
  facet_grid(.~ genotype) +
  scale_y_continuous(labels = function(x) format(x, big.mark = ",")) +
  labs(x = "Days in culture", y = "Number of genes", fill = "Direction") +
  theme(text = element_text(family = "Helvetica", size = 10),
  axis.title = element_text(size = 10),
  axis.text = element_text(size = 8)) 

bar1

ggsave(filename = "results/deseq2/degene_n_bar.pdf", plot = bar1, width = 6, height = 4, units = "in")

#################################################
##################################################

##load gplot and make the heatmap
install.packages("gplots")
library(gplots)
heatmap.2(FC_all_Bic_genes_HM_DMX, Rowv=TRUE, scale="none", trace="none", col=my_palette, xlab="none", 
          ylab="phylum", margins=c(10,5), key.ylab = NA, keysize = 0.5, density.info = "none", cexRow = 0.2)

heatmap.2(FC_all_Bic_genes_HM_DMX, col=my_palette, xlab="none", 
          keysize = 0.5, density.info = "none", cexRow = 0.2)

#test if Bic-regulated genes are behaving differently MM401 culture compared to DMSO culture
##select Bic-regulated genes in DMSO at the All Restable
#Select All Bic-regulated genes
Bic_UP_1hr_DMSO <- subset(All_Bic_restable, padj_1hr_vs_0h_DMSO < 0.05 & log2FC_1hr_vs_0h_DMSO > 0)
Bic_DOWN_1hr_DMSO <- subset(All_Bic_restable, padj_1hr_vs_0h_DMSO < 0.05 & log2FC_1hr_vs_0h_DMSO < 0)
Bic_UP_4hr_DMSO <- subset(All_Bic_restable, padj_4hr_vs_0h_DMSO < 0.05 & log2FC_4hr_vs_0h_DMSO > 0)                        
Bic_DOWN_4hr_DMSO <- subset(All_Bic_restable, padj_4hr_vs_0h_DMSO < 0.05 & log2FC_4hr_vs_0h_DMSO < 0)                        

r4r
