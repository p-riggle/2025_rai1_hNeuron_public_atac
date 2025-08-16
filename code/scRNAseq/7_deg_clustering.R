# perform clustering to group genes based on differential expression between days and genotypes.
# ref https://2-bitbio.com/2017/04/clustering-rnaseq-data-making-heatmaps.html


library(DESeq2)

# import rlog transformed count data
rld <- readRDS('data/scRNAseq/iwaselab/rld')
head (rld)
rld <- as.data.frame(assay(rld))
head(rld)






