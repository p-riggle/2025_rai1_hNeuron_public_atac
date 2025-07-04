# differential gene expression analysis of Kaster et al. 2023 RNAseq data

library(DESeq2)
library(dplyr)

# load data 
count2 <- read.delim("data/blueberry_kaster2023rnaseq/combine_counts/gene_expected_count.txt")
head(count2)
colnames(count2) <- gsub(".expected_count","",colnames(count2))
colnames(count2)

#set the geneid as rownames
rownames(count2) <- count2$gene_id
count2$gene_id <- NULL

# round numbers
count2 <- round(count2)
colnames(count2)

# reorder columns
count2 <- count2[, c(8,10,12,7,9,11,2,4,6,1,3,5)]
colnames(count2)

# generate experimental design2 file
design2 <- data.frame(row.names = colnames(count2), genotype = c("wt", "wt", "wt", "wt", "wt", "wt", "ko", "ko", "ko", "ko", "ko", "ko"))

design2$day <- c("8", "8", "8", "12", "12", "12", "8", "8", "8", "12", "12", "12")

design2$condition <- c("wt_d8", "wt_d8", "wt_d8", "wt_d12", "wt_d12", "wt_d12", "ko_d8", "ko_d8", "ko_d8", "ko_d12", "ko_d12", "ko_d12")

design2 

dds2 <- DESeqDataSetFromMatrix(countData = count2, colData = design2, design = ~ condition)
dds2_run <- DESeq(dds2, betaPrior = TRUE)

# pca plots 
vsd_dds2 <- vst(dds2, blind = TRUE)
head(assay(vsd_dds2), 3)

library(ggplot2)
str(vsd_dds2)

pca <- plotPCA(vsd_dds2, intgroup = c("condition"))
pca_data <- plotPCA(vsd_dds2, intgroup = c("genotype", "day"), returnData = TRUE)
pca
pca_data

# creat a shape and color mapping
color_mapping2 <- c("8" = "blue", "12" = "red")
shape_mapping2 <- c("wt" = 1, "ko" = 4)

# a new pca plot with modified aestherics
pca3 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = day, shape = genotype)) +
  geom_point(size = 4) +
  scale_color_manual(values = color_mapping2) +
  scale_shape_manual(values = shape_mapping2) +
  ggtitle("PCA Plot of Kaster et al 2023 RNA-seq Data") +
  theme(text = element_text(family = "Helvetica", size = 10),
  axis.title = element_text(size = 10),
  axis.text = element_text(size = 8)) +
  theme_classic()

pca3

# save the pca plot
ggsave("results/deseq2/kaster_2023/pca_plot_kaster2023.pdf", pca3, width = 5, height = 5, dpi = 300)

# differential gene expression analysis
rstbl_d8_wt.vs.ko <- results(dds2_run, contrast = c("condition", "wt_d8", "ko_d8"))
rstbl_d12_wt.vs.ko <- results(dds2_run, contrast = c("condition", "wt_d12", "ko_d12"))

rstbl_d8 <- as.data.frame(rstbl_d8_wt.vs.ko)
rstbl_d12 <- as.data.frame(rstbl_d12_wt.vs.ko)

# add symbol and save result table
library(AnnotationDbi)
library(org.Mm.eg.db)
columns(org.Mm.eg.db)

# remove version number
head(rstbl_d8)
rownames(rstbl_d8)<- gsub("\\..*", "", row.names(rstbl_d8))
rownames(rstbl_d12) <- gsub("\\..*", "", row.names(rstbl_d12))

rstbl_d8$symbol <- mapIds(org.Mm.eg.db, keys = row.names(rstbl_d8), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

rstbl_d12$symbol <- mapIds(org.Mm.eg.db, keys = row.names(rstbl_d12), column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")

# save the result table
write.table(rstbl_d8, file = "results/deseq2/kaster_2023/d8_wt.vs.ko.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
write.table(rstbl_d12, file = "results/deseq2/kaster_2023/d12_wt.vs.ko.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

# count DEGs
nrow(subset(rstbl_d8, padj < 0.05))
nrow(subset(rstbl_d12, padj < 0.05))
nrow(subset(rstbl_d8, padj < 0.1))
nrow(subset(rstbl_d12, padj < 0.1))
nrow(subset(rstbl_d8, padj < 0.2))
nrow(subset(rstbl_d12, padj < 0.2))

# make a barplot of deg numbers obtained with different p value cutoffs
library(ggplot2)    
library(reshape2)

# create a data frame
deg_num <- data.frame(pvalue = c(0.05, 0.1, 0.2), d8 = c(4, 4, 5), d12 = c(106, 219, 695))

# reshape the data frame
deg_num_melt <- melt(deg_num, id.vars = "pvalue", variable.name = "day", value.name =   
"deg_num")

deg_num_melt
deg_num_melt$pvalue <- as.factor(deg_num_melt$pvalue) 

# make a barplot
# Define a color palette with less vivid colors
my_colors <- c("#3d58f1", "#2d049e")

# Make a bar plot with the updated color palette
deg_num_plot <- ggplot(deg_num_melt, aes(x = pvalue, y = deg_num, fill = day)) +
  geom_bar(stat = "identity", position = "dodge") +
  scale_fill_manual(values = my_colors) +
  ggtitle("Number of DEGs Obtained with Different P Value Cutoffs") +
  theme(text = element_text(family = "Helvetica", size = 10),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 8)) +
  theme_classic()

deg_num_plot

# save the plot
ggsave("results/deseq2/kaster_2023/deg_num_plot.pdf", deg_num_plot, width = 5, height = 5, dpi = 300)

deg_num_plot

