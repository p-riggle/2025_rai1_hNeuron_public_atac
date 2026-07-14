#Q: How many genes go up or down at d0, d3 and d7? ####
# this script generates PCA plot and bar plot of DEgene numbers foe ineuron RNA-seq Data set 

## Load DEseq2 ####
# installed.packages("DESeq2")
# install.packages("BiocManager")
# install.packages("dplyr")

library(DESeq2)

#Load forestmngr to do rounding of numbers
library(dplyr)

#import gene count file
count <- read.delim("data/ineuron/gene_expected_count.annot.txt")
colnames(count) <- gsub(".expected_count","",colnames(count))
colnames(count)


#set the geneid as rownames
rownames(count) <- count$gene_id
count$gene_id <- NULL
head(count)

# Keep Only Numeric Columns Before Rounding.
count <- count[, sapply(count, is.numeric)]
count <- round(count)

#This ensures that round() is applied only to numeric data. Check Data Types Before Proceeding.

str(count)
View(count)


# Check the updated data
head(count)


#Reorder the columns
colnames(count)
count <- count[, c(1,12,23,31,5:8,18:21,32:35,9:13,22:26,36,2:4,14:17,27:30)] 
colnames(count)

#Remove Duplicates
count <- count[, -c(20, 23)]
colnames(count)

#generate experimental design file so called “se”. 
design <- data.frame(row.names = colnames(count), conditions = c(
  "Rai1.wt_d0", "Rai1.wt_d0", "Rai1.wt_d0", "Rai1.wt_d0",
  "Rai1.wt_d3", "Rai1.wt_d3", "Rai1.wt_d3", "Rai1.wt_d3",
  "Rai1.wt_d7", "Rai1.wt_d7", "Rai1.wt_d7", "Rai1.wt_d7",
  "Rai1.het_d0", "Rai1.het_d0", "Rai1.het_d0", "Rai1.het_d0",   
  "Rai1.het_d3", "Rai1.het_d3", "Rai1.het_d3", "Rai1.het_d3",
  "Rai1.het_d7", "Rai1.het_d7", "Rai1.het_d7", "Rai1.het_d7",
  "Rai1.ko_d0", "Rai1.ko_d0", "Rai1.ko_d0", "Rai1.ko_d0",   
  "Rai1.ko_d3", "Rai1.ko_d3", "Rai1.ko_d3", "Rai1.ko_d3",
  "Rai1.ko_d7", "Rai1.ko_d7", "Rai1.ko_d7", "Rai1.ko_d7")
  
)

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
color_mapping <- c("0" = "lightblue", "3" = "cyan",  "7" = "red")

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

ggsave(filename = "results/ineuron/deseq2/pca2.pdf", plot = pca2, width = 6, height = 4.5, units = "in")

#Identifying Bic-induced genes in DMSO ####

rstbl_d0_wt.vs.het <- results(dds_run, contrast = c("conditions", "Rai1.het_d0", "Rai1.wt_d0"), alpha = 0.05)
rstbl_d0_wt.vs.ko <- results(dds_run, contrast = c("conditions", "Rai1.ko_d0", "Rai1.wt_d0"), alpha = 0.05)
rstbl_d3_wt.vs.het <- results(dds_run, contrast = c("conditions", "Rai1.het_d3", "Rai1.wt_d3"), alpha = 0.05)
rstbl_d3_wt.vs.ko <- results(dds_run, contrast = c("conditions", "Rai1.ko_d3", "Rai1.wt_d3"), alpha = 0.05)
rstbl_d7_wt.vs.het <- results(dds_run, contrast = c("conditions", "Rai1.het_d7", "Rai1.wt_d7"), alpha = 0.05)
rstbl_d7_wt.vs.ko <- results(dds_run, contrast = c("conditions", "Rai1.ko_d7", "Rai1.wt_d7"), alpha = 0.05)

summary(rstbl_d0_wt.vs.het)
summary(rstbl_d0_wt.vs.ko)
summary(rstbl_d3_wt.vs.het)
summary(rstbl_d3_wt.vs.ko)
summary(rstbl_d7_wt.vs.het)
summary(rstbl_d7_wt.vs.ko)


#Add ENSEMBL and entrezID and print result tables
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("AnnotationDbi")
library(AnnotationDbi)


library(org.Hs.eg.db)
columns(org.Hs.eg.db) 

#print the keytypes 

## Write restables for Bic-regulated genes

add_symbol <- function(restable, filename) {
  restable$entrezid <- mapIds(org.Hs.eg.db, keys=row.names(restable), column="ENTREZID", keytype="ENSEMBL", multiVals="first")
  write.table(restable, file = filename, sep = "\t")
}
restable_list <- list(rstbl_d0_wt.vs.het, rstbl_d0_wt.vs.ko, rstbl_d3_wt.vs.het, rstbl_d3_wt.vs.ko, rstbl_d7_wt.vs.het, rstbl_d7_wt.vs.ko)

restable_names <- c("rstbl_d0_wt.vs.het", "rstbl_d0_wt.vs.ko", "rstbl_d3_wt.vs.het", "rstbl_d3_wt.vs.ko", "rstbl_d7_wt.vs.het", "rstbl_d7_wt.vs.ko")

#Verify that restable_list is populated:
length(restable_list)

library(AnnotationDbi)
library(org.Hs.eg.db)
#Run the loop:
for (i in 1:length(restable_list)) {
  restable <- restable_list[[i]]
  filename <- paste0("results/ineuron/deseq2/", restable_names[i], ".txt")
  add_symbol(restable, filename)
}



# read the restable
rstbl_d0_wt.vs.het <- read.delim("results/ineuron/deseq2/rstbl_d0_wt.vs.het.txt")
rstbl_d0_wt.vs.ko <- read.delim("results/ineuron/deseq2/rstbl_d0_wt.vs.ko.txt")
rstbl_d3_wt.vs.het <- read.delim("results/ineuron/deseq2/rstbl_d3_wt.vs.het.txt")
rstbl_d3_wt.vs.ko <- read.delim("results/ineuron/deseq2/rstbl_d3_wt.vs.ko.txt")
rstbl_d7_wt.vs.het <- read.delim("results/ineuron/deseq2/rstbl_d7_wt.vs.het.txt")
rstbl_d7_wt.vs.ko <- read.delim("results/ineuron/deseq2/rstbl_d7_wt.vs.ko.txt")


head(rstbl_d0_wt.vs.het)

up_d0_wt.vs.het <- subset(rstbl_d0_wt.vs.het, padj < 0.05 & log2FoldChange > 0)
up_d0_wt.vs.ko <- subset(rstbl_d0_wt.vs.ko, padj < 0.05 & log2FoldChange > 0)
up_d3_wt.vs.het <- subset(rstbl_d3_wt.vs.het, padj < 0.05 & log2FoldChange > 0)
up_d3_wt.vs.ko <- subset(rstbl_d3_wt.vs.ko, padj < 0.05 & log2FoldChange > 0)
up_d7_wt.vs.het <- subset(rstbl_d7_wt.vs.het, padj < 0.05 & log2FoldChange > 0)
up_d7_wt.vs.ko <- subset(rstbl_d7_wt.vs.ko, padj < 0.05 & log2FoldChange > 0)


down_d0_wt.vs.het <- subset(rstbl_d0_wt.vs.het, padj < 0.05 & log2FoldChange < 0)
down_d0_wt.vs.ko <- subset(rstbl_d0_wt.vs.ko, padj < 0.05 & log2FoldChange < 0)
down_d3_wt.vs.het <- subset(rstbl_d3_wt.vs.het, padj < 0.05 & log2FoldChange < 0)
down_d3_wt.vs.ko <- subset(rstbl_d3_wt.vs.ko, padj < 0.05 & log2FoldChange < 0)
down_d7_wt.vs.het <- subset(rstbl_d7_wt.vs.het, padj < 0.05 & log2FoldChange < 0)
down_d7_wt.vs.ko <- subset(rstbl_d7_wt.vs.ko, padj < 0.05 & log2FoldChange < 0)


# plot DEgene numbers
condition <- c("up_d0_wt.vs.het", "up_d0_wt.vs.ko", "up_d3_wt.vs.het", "up_d3_wt.vs.ko", "up_d7_wt.vs.het", "up_d7_wt.vs.ko", "down_d0_wt.vs.het", "down_d0_wt.vs.ko", "down_d3_wt.vs.het", 
               "down_d3_wt.vs.ko", "down_d7_wt.vs.het", "down_d7_wt.vs.ko" )

direction <- c("Up", "Up", "Up", "Up", "Up", "Up", "Down", "Down", "Down", "Down", "Down", "Down")

time <- c("d0", "d0", "d3", "d3", "d7", "d7", "d0", "d0", "d3", "d3", "d7", "d7")

genotype <- c("Heterozygous", "Homozygous", "Heterozygous", "Homozygous", "Heterozygous", "Homozygous", "Heterozygous", "Homozygous", "Heterozygous", "Homozygous", "Heterozygous", "Homozygous")

gene_n <- c(nrow(up_d0_wt.vs.het), nrow(up_d0_wt.vs.ko),
            nrow(up_d3_wt.vs.het), nrow(up_d3_wt.vs.ko),
            nrow(up_d7_wt.vs.het), nrow(up_d7_wt.vs.ko),
            nrow(down_d0_wt.vs.het), nrow(down_d0_wt.vs.ko),
            nrow(down_d3_wt.vs.het), nrow(down_d3_wt.vs.ko),
            nrow(down_d7_wt.vs.het), nrow(down_d7_wt.vs.ko))


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

ggsave(filename = "results/ineuron/deseq2/degene_n_bar.pdf", plot = bar1, width = 6, height = 4, units = "in")

# save the deg list as a excel file with each tab wt.vs.het and wt.vs.ko for each day
# install.packages("openxlsx")
library(openxlsx)
wb <- createWorkbook()
addWorksheet(wb, "d0_wt.vs.het")
writeData(wb, "d0_wt.vs.het", rstbl_d0_wt.vs.het)
addWorksheet(wb, "d0_wt.vs.ko")
writeData(wb, "d0_wt.vs.ko", rstbl_d0_wt.vs.ko)
addWorksheet(wb, "d3_wt.vs.het")
writeData(wb, "d3_wt.vs.het", rstbl_d3_wt.vs.het)
addWorksheet(wb, "d3_wt.vs.ko")
writeData(wb, "d3_wt.vs.ko", rstbl_d3_wt.vs.ko)
addWorksheet(wb, "d7_wt.vs.het")
writeData(wb, "d7_wt.vs.het", rstbl_d7_wt.vs.het)
addWorksheet(wb, "d7_wt.vs.ko")
writeData(wb, "d7_wt.vs.ko", rstbl_d7_wt.vs.ko)
saveWorkbook(wb, "results/ineuron/deseq2/degene_list.xlsx", overwrite = TRUE)


