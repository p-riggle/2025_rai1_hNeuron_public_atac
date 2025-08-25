# what biological pathways are enriched in each cluster?
library(clusterProfiler)

# import the cluster matrix
mh_d20 <- read.delim("results/scRNAseq/4_pseudobulk/mh_d20.txt")
head(mh_d20)

mh_d60 <- read.delim("results/scRNAseq/4_pseudobulk/mh_d60.txt")
head(mh_d60)

# extract unique genes from each cluster.
# install.packages("dplyr")
library(dplyr)

d20_list <- list()
for (i in 1:6) {
    d20_list[[i]] <- mh_d20[mh_d20$cluster == i,] %>% select(gene) %>% unlist() %>% unique()
}
names(d20_list) <- paste0("mh_d20_c", 1:6)
list2env(d20_list, envir = .GlobalEnv)

d60_list <- list()
for (i in 1:4) {
  d60_list[[i]] <- mh_d60[mh_d60$cluster == i,] %>% select(gene) %>% unlist() %>% unique()
}
names(d60_list) <- paste0("mh_d60_c", 1:4)
list2env(d60_list, envir = .GlobalEnv)

str(d20_list)

# combine the lists into a list.
day20_cgenes <- list(mh_d20_c1, mh_d20_c2, mh_d20_c3, mh_d20_c4, mh_d20_c5, mh_d20_c6)
day60_cgenes <- list(mh_d60_c1, mh_d60_c2, mh_d60_c3, mh_d60_c4)

# name the clusters
names(day20_cgenes) <- c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Cluster5", "Cluster6")
names(day60_cgenes) <- c("Cluster1", "Cluster2", "Cluster3", "Cluster4")

str(day20_cgenes)
str(day60_cgenes)

# convert symbols to entrez ids
library(org.Hs.eg.db)

day20_cgenes <- lapply(day20_cgenes, function(x) {
    entrez <- mapIds(org.Hs.eg.db, x, "ENTREZID", "SYMBOL")
    entrez <- unlist(entrez)
    entrez <- entrez[!is.na(entrez)]
    return(entrez)
})

day60_cgenes <- lapply(day60_cgenes, function(x) {
    entrez <- mapIds(org.Hs.eg.db, x, "ENTREZID", "SYMBOL")
    entrez <- unlist(entrez)
    entrez <- entrez[!is.na(entrez)]
    return(entrez)
})

str(day20_cgenes)
str(day60_cgenes)

# run clusterProfiler KEGG.
kegg <- function(cmatrix) {compareCluster(geneCluster = cmatrix, fun = enrichKEGG) %>% setReadable(OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
}

kegg_d20 <- kegg(day20_cgenes)
kegg_d60 <- kegg(day60_cgenes)

dp_kegg_d20 <- dotplot(kegg_d20, showCategory = 10, title = "KEGG enrichment of RAI1-Het: Hiearchical Clusters")
dp_kegg_d20
