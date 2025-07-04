# what biological pathways are enriched in each cluster?
library(clusterProfiler)

# import the cluster matrix
het_hi <- read.delim("results/cluster/sep/mh_het.txt")
ko_hi <- read.delim("results/cluster/sep/mh_ko.txt")
het_k <- read.delim("results/cluster/sep/mk_het.txt")
ko_k <- read.delim("results/cluster/sep/mk_ko.txt")

str(het_hi)
str(het_k)

# extract unique genes from each cluster.
install.packages("dplyr")
library(dplyr)

het_hi_list <- list()
for (i in 1:6) {
  het_hi_list[[i]] <- het_hi[het_hi$cluster == i,] %>% select(gene) %>% unlist()
}
names(het_hi_list) <- paste0("het_hi", 1:6)
list2env(het_hi_list, envir = .GlobalEnv)

het_k_list <- list()
for (i in 1:6) {
  het_k_list[[i]] <- het_k[het_k$cluster == i,] %>% select(gene) %>% unlist()
}
names(het_k_list) <- paste0("het_k", 1:6)
list2env(het_k_list, envir = .GlobalEnv)

ko_hi_list <- list()
for (i in 1:5) {
  ko_hi_list[[i]] <- ko_hi[ko_hi$cluster == i,] %>% select(gene) %>% unlist()
}
names(ko_hi_list) <- paste0("ko_hi", 1:5)
list2env(ko_hi_list, envir = .GlobalEnv)

ko_k_list <- list()
for (i in 1:5) {
  ko_k_list[[i]] <- ko_k[ko_k$cluster == i,] %>% select(gene) %>% unlist()
}
names(ko_k_list) <- paste0("ko_k", 1:5)
list2env(ko_k_list, envir = .GlobalEnv)

# combine the lists into a list.
het_hi_cgenes <- list(het_hi1, het_hi2, het_hi3, het_hi4, het_hi5, het_hi6)
het_k_cgenes <- list(het_k1, het_k2, het_k3, het_k4, het_k5, het_k6)
ko_hi_cgenes <- list(ko_hi1, ko_hi2, ko_hi3, ko_hi4, ko_hi5)
ko_k_cgenes <- list(ko_k1, ko_k2, ko_k3, ko_k4, ko_k5)

# name the clusters
names(het_hi_cgenes) <- c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Cluster5", "Cluster6")
names(het_k_cgenes) <- c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Cluster5", "Cluster6")
names(ko_hi_cgenes) <- c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Cluster5")
names(ko_k_cgenes) <- c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Cluster5")

het_hi_cgenes

# convert symbols to entrez ids
library(org.Hs.eg.db)

het_hi_cgenes <- lapply(het_hi_cgenes, function(x) {
    entrez <- mapIds(org.Hs.eg.db, x, "ENTREZID", "SYMBOL")
    entrez <- unlist(entrez)
    entrez <- entrez[!is.na(entrez)]
    return(entrez)
})

het_k_cgenes <- lapply(het_k_cgenes, function(x) {
    entrez <- mapIds(org.Hs.eg.db, x, "ENTREZID", "SYMBOL")
    entrez <- unlist(entrez)
    entrez <- entrez[!is.na(entrez)]
    return(entrez)
})

ko_hi_cgenes <- lapply(ko_hi_cgenes, function(x) {
    entrez <- mapIds(org.Hs.eg.db, x, "ENTREZID", "SYMBOL")
    entrez <- unlist(entrez)
    entrez <- entrez[!is.na(entrez)]
    return(entrez)
})

ko_k_cgenes <- lapply(ko_k_cgenes, function(x) {
    entrez <- mapIds(org.Hs.eg.db, x, "ENTREZID", "SYMBOL")
    entrez <- unlist(entrez)
    entrez <- entrez[!is.na(entrez)]
    return(entrez)
})

# run clusterProfiler KEGG.
kegg <- function(cmatrix) {compareCluster(geneCluster = cmatrix, fun = enrichKEGG) %>% setReadable(OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
}

kegg_het_hi <- kegg(het_hi_cgenes)
kegg_het_k <- kegg(het_k_cgenes)
kegg_ko_hi <- kegg(ko_hi_cgenes)
kegg_ko_k <- kegg(ko_k_cgenes)

# make dotplots
dp_het_hi <- dotplot(kegg_het_hi, showCategory = 10, title = "KEGG enrichment of RAI1-Het: Hiearchical Clusters")
dp_het_k <- dotplot(kegg_het_k, showCategory = 10, title = "KEGG enrichment of RAI1-Het: K-mean Clusters")

dp_ko_hi <- dotplot(kegg_ko_hi, showCategory = 10, title = "KEGG enrichment of RAI1-KO: Hiearchical Clusters")
dp_ko_k <- dotplot(kegg_ko_k, showCategory = 10, title = "KEGG enrichment of RAI1-KO: K-mean clusters")
dp_ko_k

# save the two dotplots arranged side by side
library(gridExtra)
pdf(sprintf("results/cluster/sep/%s_kegg_het.pdf", Sys.Date()), width = 20, height = 15)
grid.arrange(dp_het_hi, dp_het_k, ncol = 2)
dev.off()

pdf(sprintf("results/cluster/sep/%s_kegg_ko.pdf", Sys.Date()), width = 20, height = 15)
grid.arrange(dp_ko_hi, dp_ko_k, ncol = 2)
dev.off()

# run clusterProfiler GO. 
go <- function(cmatrix) {compareCluster(geneCluster = cmatrix, fun = enrichGO, OrgDb = org.Hs.eg.db) %>% setReadable(OrgDb = org.Hs.eg.db, keyType = "ENTREZID")}

go_het_hi <- go(het_hi_cgenes)
go_het_k <- go(het_k_cgenes)
go_ko_hi <- go(ko_hi_cgenes)
go_ko_k <- go(ko_k_cgenes)

# make dotplots
dpgo_het_hi <- dotplot(go_het_hi, showCategory = 10, title = "GO enrichment of RAI1-Het: Hiearchical Clusters")
dpgo_het_k <- dotplot(go_het_k, showCategory = 10, title = "GO enrichment of RAI1-Het: K-mean Clusters")
dpgo_ko_hi <- dotplot(go_ko_hi, showCategory = 10, title = "GO enrichment of RAI1-KO: Hiearchical Clusters")
dpgo_ko_k <- dotplot(go_ko_k, showCategory = 10, title = "GO enrichment of RAI1-KO: K-mean clusters")

# save the two dotplots arranged side by side
pdf(sprintf("results/cluster/sep/%s_go_het.pdf", Sys.Date()), width = 20, height = 15)
grid.arrange(dpgo_het_hi, dpgo_het_k, ncol = 2)
dev.off()

pdf(sprintf("results/cluster/sep/%s_go_ko.pdf", Sys.Date()), width = 20, height = 15)
grid.arrange(dpgo_ko_hi, dpgo_ko_k, ncol = 2)
dev.off()

# run clusterProfiler Go of cellular components.
gocc <- function(cmatrix) {compareCluster(geneCluster = cmatrix, fun = enrichGO, OrgDb = org.Hs.eg.db, ont = "CC") %>% setReadable(OrgDb = org.Hs.eg.db, keyType = "ENTREZID")}

gocc_het_hi <- gocc(het_hi_cgenes)
gocc_het_k <- gocc(het_k_cgenes)
gocc_ko_hi <- gocc(ko_hi_cgenes)
gocc_ko_k <- gocc(ko_k_cgenes)

# make dotplots
dpgocc_het_hi <- dotplot(gocc_het_hi, showCategory = 10, title = "GO-CC enrichment of RAI1-Het: Hiearchical Clusters")
dpgocc_het_k <- dotplot(gocc_het_k, showCategory = 10, title = "GO-CC enrichment of RAI1-Het: K-mean Clusters")
dpgocc_ko_hi <- dotplot(gocc_ko_hi, showCategory = 10, title = "GO-CC enrichment of RAI1-KO: Hiearchical Clusters")
dpgocc_ko_k <- dotplot(gocc_ko_k, showCategory = 10, title = "GO-CC enrichment of RAI1-KO: K-mean clusters")

# save the two dotplots arranged side by side
pdf(sprintf("results/cluster/sep/%s_gocc_het.pdf", Sys.Date()), width = 20, height = 15)
grid.arrange(dpgocc_het_hi, dpgocc_het_k, ncol = 2)
dev.off()

pdf(sprintf("results/cluster/sep/%s_gocc_ko.pdf", Sys.Date()), width = 20, height = 15)
grid.arrange(dpgocc_ko_hi, dpgocc_ko_k, ncol = 2)
dev.off()



