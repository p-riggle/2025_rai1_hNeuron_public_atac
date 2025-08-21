#Q: How many genes go up or down at d0, d3 and d7? ####
# this script generates PCA plot and bar plot of DEgene numbers. 

## Load DEseq2 ####
#BiocManager::install("DESeq2")
#install.packages("BiocManager")
#install.packages("dplyr")


#BiocManager::install("DESeq2")
library(DESeq2)

#Load forestmngr to do rounding of numbers
library(dplyr)

#import gene count file
count <- read.delim("/Users/satabdi/Downloads/rscripts/ineurons_gene_expected_count.annot (1).txt")
colnames(count) <- gsub(".expected_count","",colnames(count))
colnames(count)

#set the Geneid as rownames
rownames(count) <- count$gene_id
count$gene_id <- NULL


# Keep Only Numeric Columns Before Rounding. Instead of manually dropping each column, select only numeric ones:
count <- count[, sapply(count, is.numeric)]
count <- round(count)

#This ensures that round() is applied only to numeric data. Check Data Types Before Proceeding.
#To confirm which columns are non-numeric:
str(count)

# Remove ".expected_count" from the column names
colnames(count) <- gsub(".expected_count", "", colnames(count))

# Set the gene_id as rownames
rownames(count) <- count$gene_id
count$gene_id <- NULL
  
# Subset only the numeric columns (from column 1 onwards, based on your data structure)
numeric_columns <- count[, 1:ncol(count)]

# Round the numeric values
 numeric_columns <- round(numeric_columns)
 
# Replace the original numeric columns with the rounded values
count[, 1:ncol(count)] <- numeric_columns


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

ggsave(filename = "/Users/satabdi/Downloads/rscripts/Results/pca2.pdf", plot = pca2, width = 6, height = 4.5, units = "in")

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


#Add gene symbol and entrezID and print result tables
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("AnnotationDbi")
library(AnnotationDbi)


library(org.Hs.eg.db)
columns(org.Hs.eg.db) 

#print the keytypes 

## Write restables for Bic-regulated genes

add_symbol <- function(restable, filename) {
  restable$entrezid <- mapIds(org.Hs.eg.db, keys=row.names(restable), column="ENTREZID", keytype="SYMBOL", multiVals="first")
  write.table(restable, file = filename, sep = "\t")
}

restable_list <- list(rstbl_d0_wt.vs.het, rstbl_d0_wt.vs.ko, rstbl_d3_wt.vs.het, rstbl_d3_wt.vs.ko, rstbl_d7_wt.vs.het, rstbl_d7_wt.vs.ko)

restable_names <- c("rstbl_d0_wt.vs.het", "rstbl_d0_wt.vs.ko", "rstbl_d3_wt.vs.het", "rstbl_d3_wt.vs.ko", "rstbl_d7_wt.vs.het", "rstbl_d7_wt.vs.ko")

#Verify that restable_list is populated:
length(restable_list)

#3. Modify the add_symbol function:
add_symbol <- function(restable, filename) {
  # Check that SYMBOL column exists
  if ("SYMBOL" %in% colnames(restable)) {
    restable$entrezid <- mapIds(org.Hs.eg.db, keys = restable$SYMBOL, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  } else {
    warning("SYMBOL column not found in the data frame")
  }
  write.table(restable, file = filename, sep = "\t", row.names = FALSE)
}

#Run the loop:
for (i in 1:length(restable_list)) {
  restable <- restable_list[[i]]
  filename <- paste0("/Users/satabdi/Downloads/rscripts/Results/", restable_names[i], ".txt")
  add_symbol(restable, filename)
}

# read the restable
rstbl_d0_wt.vs.het <- read.delim("/Users/satabdi/Downloads/rscripts/Results/rstbl_d0_wt.vs.het.txt")
rstbl_d0_wt.vs.ko <- read.delim("/Users/satabdi/Downloads/rscripts/Results/rstbl_d0_wt.vs.ko.txt")
rstbl_d3_wt.vs.het <- read.delim("/Users/satabdi/Downloads/rscripts/Results/rstbl_d3_wt.vs.het.txt")
rstbl_d3_wt.vs.ko <- read.delim("/Users/satabdi/Downloads/rscripts/Results/rstbl_d3_wt.vs.ko.txt")
rstbl_d7_wt.vs.het <- read.delim("/Users/satabdi/Downloads/rscripts/Results/rstbl_d7_wt.vs.het.txt")
rstbl_d7_wt.vs.ko <- read.delim("/Users/satabdi/Downloads/rscripts/Results/rstbl_d7_wt.vs.ko.txt")


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

ggsave(filename = "/Users/satabdi/Downloads/rscripts/Results/degene_n_bar.pdf", plot = bar1, width = 6, height = 4, units = "in")

#################################################
##################################################

##load gplot and make the heatmap
install.packages("gplots")
library(gplots)
heatmap.2(FC_all_Bic_genes_HM_DMX, Rowv=TRUE, scale="none", trace="none", col=my_palette, xlab="none", 
          ylab="phylum", margins=c(10,5), key.ylab = NA, keysize = 0.5, density.info = "none", cexRow = 0.2)

heatmap.2(FC_all_Bic_genes_HM_DMX, col=my_palette, xlab="none", 
          keysize = 0.5, density.info = "none", cexRow = 0.2)



