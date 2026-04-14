# what biological pathways are enriched in each cluster?
# this code generates the GO cellular component plot in Fig7C

# install cluterProfiler
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("clusterProfiler")
library(clusterProfiler)

BiocManager::install("org.Mm.eg.db")

library(org.Mm.eg.db)
library(AnnotationDbi)

library(ggplot2)
library(dplyr) 

#import the kmeans cluster dataframe.
clk <- read.delim("result/Fig7_huang/mk_ENCODE_WTctx_pc1300genesclusterdataframe.txt")

str(clk)

# extract unique genes from each cluster.

clk1 <- clk %>%
  dplyr::filter(cluster == 1) %>%
  dplyr::pull(ensembl_id) %>%
  unique()

clk2 <- clk %>%
  dplyr::filter(cluster == 2) %>%
  dplyr::pull(ensembl_id) %>%
  unique()

# combine the lists into a list, and name the clusters
clk <- list(Cluster1 = clk1, Cluster2 = clk2)

str(clk)

# convert ENSEMBL to entrez ids.

clk_entrez <- lapply(clk, function(x) {
  entrez <- mapIds(
    org.Mm.eg.db, 
    keys = unique(x),
    column = "ENTREZID",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  entrez <- unlist(entrez)
  entrez <- entrez[!is.na(entrez)]
  return(entrez)
})

str(clk_entrez)

# run clusterProfiler GO.

clkgo <- compareCluster(geneCluster = clk_entrez, fun = enrichGO, OrgDb='org.Mm.eg.db', ont= "CC", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2)
clkgo <- setReadable(clkgo, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
head(clkgo)

dpgo1 <- dotplot(clkgo, showCategory = 5, title = "Kmeans clustering")
dpgo1

# save the GO plot
ggsave(sprintf("result/Fig7_huang/%s_encodepc1top300genes_dp_kmeans_GO_Fig7C.pdf", Sys.Date()), plot = dpgo1, width = 5, height = 5)

# save clkgo result as csv
write.csv(as.data.frame(clkgo), file = sprintf("result/Fig7_huang/%s_huangpc1top300genes_kmeans_GO_resultsMarch.csv", Sys.Date()))

# end of script
