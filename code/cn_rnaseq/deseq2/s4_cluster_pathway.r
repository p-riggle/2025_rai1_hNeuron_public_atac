# what biological pathways are enriched in each cluster?

# install cluterProfiler
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# BiocManager::install("clusterProfiler")
library(clusterProfiler)
browseVignettes("clusterProfiler")

data(gcSample)
str(gcSample)

#import the hiearchical cluster dataframe.
clh <- read.delim("results/cluster/clh.txt")
#import the k-mean cluster dataframe.
clk <- read.delim("results/cluster/clk.txt")
str(clh)
str(clk)

# extract unique genes from each cluster.
library(dplyr)
cl1 <- clh[clh$cluster == 1,] %>% select(gene) %>% unlist()
cl2 <- clh[clh$cluster == 2,] %>% select(gene) %>% unlist()
cl3 <- clh[clh$cluster == 3,] %>% select(gene) %>% unlist()
cl4 <- clh[clh$cluster == 4,] %>% select(gene) %>% unlist()
cl5 <- clh[clh$cluster == 5,] %>% select(gene) %>% unlist()
cl6 <- clh[clh$cluster == 6,] %>% select(gene) %>% unlist()
cl7 <- clh[clh$cluster == 7,] %>% select(gene) %>% unlist()

clk1 <- clk[clk$cluster == 1,] %>% select(gene) %>% unlist()
clk2 <- clk[clk$cluster == 2,] %>% select(gene) %>% unlist()
clk3 <- clk[clk$cluster == 3,] %>% select(gene) %>% unlist()
clk4 <- clk[clk$cluster == 4,] %>% select(gene) %>% unlist()
clk5 <- clk[clk$cluster == 5,] %>% select(gene) %>% unlist()
clk6 <- clk[clk$cluster == 6,] %>% select(gene) %>% unlist()
clk7 <- clk[clk$cluster == 7,] %>% select(gene) %>% unlist()

# combine the lists into a list.
cl <- list(cl1, cl2, cl3, cl4, cl5, cl6, cl7)
clk <- list(clk1, clk2, clk3, clk4, clk5, clk6, clk7)
str(cl)
str(clk)

# name the clusters. 
names(cl) <- c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Cluster5", "Cluster6", "Cluster7")
names(clk) <- c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Cluster5", "Cluster6", "Cluster7")

# comvert symbols to entrez ids.
## first install org.hs.eg.db
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

# convert symbols to entrez ids.
cl <- lapply(cl, function(x) {
    entrez <- mapIds(org.Hs.eg.db, x, "ENTREZID", "SYMBOL")
    entrez <- unlist(entrez)
    entrez <- entrez[!is.na(entrez)]
    return(entrez)
})

clk <- lapply(clk, function(x) {
    entrez <- mapIds(org.Hs.eg.db, x, "ENTREZID", "SYMBOL")
    entrez <- unlist(entrez)
    entrez <- entrez[!is.na(entrez)]
    return(entrez)
})
str(cl)
str(clk)

# run clusterProfiler KEGG.
clhk <- compareCluster(geneCluster = cl, fun = enrichKEGG)
clhk <- setReadable(clhk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
head(clhk)

clkk <- compareCluster(geneCluster = clk, fun = enrichKEGG)
clkk <- setReadable(clkk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
head(clkk)

# make dotplot.

dp1 <- dotplot(clhk, showCategory = 10, title = "Hierarchical clustering")
dp2 <- dotplot(clkk, showCategory = 10, title = "K-means clustering")

# save the two dotplots arranged side by side.
pdf(sprintf("results/cluster/%s_dp_kegg.pdf", Sys.Date()), width = 20, height = 15)
grid.arrange(dp1, dp2, ncol = 2)
dev.off()

# make cnetplot.
str(clhk)
cn1 <- cnetplot(clhk, showCategory = 10, title = "Hierarchical clustering")
cn2 <- cnetplot(clkk, showCategory = 10, title - "K-means clustering")

pdf(sprintf("results/cluster/%s_cnet_h_kegg.pdf", Sys.Date()), width = 10, height = 10)
cn1
dev.off()

pdf(sprintf("results/cluster/%s_cnet_k_kegg.pdf", Sys.Date()), width = 10, height = 10)
cn2
dev.off()

# run clusterProfiler GO.
clhgo <- compareCluster(geneCluster = cl, fun = enrichGO, OrgDb='org.Hs.eg.db')
clhgo <- setReadable(clhk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
head(clhgo)

clkgo <- compareCluster(geneCluster = clk, fun = enrichGO, OrgDb='org.Hs.eg.db')
clkgo <- setReadable(clkk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
head(clkgo)

dpgo1 <- dotplot(clhgo, showCategory = 10, title = "Hierarchical clustering")
dpgo1
