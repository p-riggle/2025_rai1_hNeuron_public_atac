# differential gene expression analysis of Huang et al. 2016 RNAseq data

library(DESeq2)
library(dplyr)

# load data 
count1 <- read.delim("data/blueberry_huang2016rnaseq/combine_counts/gene_expected_count.txt")
head(count1)
colnames(count1) <- gsub(".expected_count","",colnames(count1))
colnames(count1)

#set the geneid as rownames
rownames(count1) <- count1$gene_id
count1$gene_id <- NULL

# round numbers
count1 <- round(count1)
View(count1)

# reorder columns
count1 <- count1[, c(3:4, 1:2)]
colnames(count1)

# generate experimental design1 file
design1 <- data.frame(row.names = colnames(count1), condition = c("wt", "wt", "rai1-cko", "rai1-cko"))
design1

# run deseq2
design1$genotype <- c("wt", "wt", "rai1-cko", "rai1-cko")
design1

dds1 <- DESeqDataSetFromMatrix(countData = count1, colData = design1, design = ~ genotype)
dds1_run <- DESeq(dds1, betaPrior = TRUE)
resultsNames(dds1_run)

# pca plots
vsd_dds1 <- vst(dds1, blind = TRUE)
head(assay(vsd_dds1), 3)

library(ggplot2)
str(vsd_dds1)

pca <- plotPCA(vsd_dds1, intgroup = c("genotype"))
pca_data <- plotPCA(vsd_dds1, intgroup = c("genotype"), returnData = TRUE)
pca
pca_data

# creat a shae and color mapping
color_mapping1 <- c("wt" = "black", "rai1-cko" = "red")

# a new pca plot with modified aestherics 
pca2 <- ggplot(pca_data, aes(x = PC1, y = PC2, color = genotype)) +
  geom_point(size = 4) +
  scale_color_manual(values = color_mapping1) +
  ggtitle("PCA Plot of Huan et al 2016 RNA-seq Data") +
  theme(text = element_text(family = "Helvetica", size = 10),
  axis.title = element_text(size = 10),
  axis.text = element_text(size = 8)) +
  theme_classic()
pca2

ggsave(filename = "results/deseq2/huang_2016/pca2.pdf", plot = pca2, width = 6, height = 4.5, units = "in")

# identify differentially expressed genes
rstbl_wt.vs.cko <- results(dds1_run, contrast = c("genotype", "rai1-cko", "wt"), alpha = 0.05)
summary(rstbl_wt.vs.cko)

rstbl1 <- as.data.frame(rstbl_wt.vs.cko)

# add gene symbol and entrezID, print result tables
library(AnnotationDbi)

# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# BiocManager::install("org.Mm.eg.db")

library(org.Mm.eg.db)
columns(org.Mm.eg.db)

# remove version number from rownames
rownames(rstbl1) <- gsub("\\..*", "", rownames(rstbl1))
class(rownames(rstbl1))
head(rstbl1)

rstbl1$symbol <- mapIds(org.Mm.eg.db, keys=row.names(rstbl1), column="SYMBOL", keytype="ENSEMBL", multiVals="first")
head(rstbl1)

deg2 <- subset(rstbl1, padj < 0.2)
nrow(deg2)

# write table of degs 
write.table(rstbl1, file = "results/deseq2/huang_2016/huang_2016_rstbl.txt", sep = "\t")
write.table(deg2, file = "results/deseq2/huang_2016/huang_2016_deg_0.2.txt", sep = "\t")


