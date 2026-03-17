#For KO samples 

# what biological pathways are enriched in each cluster?

# install cluterProfiler
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# BiocManager::install("clusterProfiler")
library(clusterProfiler)


#import the hiearchical cluster dataframe.
clh <- read.delim("/Users/satabdi/Downloads/mh_ko.txt")
#import the k-mean cluster dataframe.
clk <- read.delim("/Users/satabdi/Downloads/mk_ko.txt")
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


clk1 <- clk[clk$cluster == 1,] %>% select(gene) %>% unlist()
clk2 <- clk[clk$cluster == 2,] %>% select(gene) %>% unlist()
clk3 <- clk[clk$cluster == 3,] %>% select(gene) %>% unlist()
clk4 <- clk[clk$cluster == 4,] %>% select(gene) %>% unlist()
clk5 <- clk[clk$cluster == 5,] %>% select(gene) %>% unlist()
clk6 <- clk[clk$cluster == 6,] %>% select(gene) %>% unlist()
clk7 <- clk[clk$cluster == 7,] %>% select(gene) %>% unlist()

# combine the lists into a list.
cl <- list(cl1, cl2, cl3, cl4, cl5, cl6)
clk <- list(clk1, clk2, clk3, clk4, clk5, clk6, clk7)
str(cl)
str(clk)

# name the clusters. 
names(cl) <- c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Cluster5", "Cluster6")
names(clk) <- c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Cluster5", "Cluster6", "Cluster7")

# comvert symbols to entrez ids.
## first install org.hs.eg.db
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# BiocManager::install("org.Hs.eg.db")
library(org.Hs.eg.db)

# convert symbols to entrez ids.
cl <- lapply(cl, function(x) {
  
  # remove NA or empty entries
  x <- x[!is.na(x)]
  
  # skip if cluster has no genes
  if(length(x) == 0){
    return(NULL)
  }
  
  entrez <- mapIds(org.Hs.eg.db,
                   keys = x,
                   column = "ENTREZID",
                   keytype = "SYMBOL",
                   multiVals = "first")
  
  entrez <- unlist(entrez)
  entrez <- entrez[!is.na(entrez)]
  
  return(entrez)
  
})

clk <- lapply(clk, function(x) {
  
  # remove NA values
  x <- x[!is.na(x)]
  
  # skip empty clusters
  if(length(x) == 0){
    return(NULL)
  }
  
  entrez <- mapIds(org.Hs.eg.db,
                   keys = x,
                   column = "ENTREZID",
                   keytype = "SYMBOL",
                   multiVals = "first")
  
  entrez <- unlist(entrez)
  entrez <- entrez[!is.na(entrez)]
  
  return(entrez)
  
})
str(cl)
str(clk)

# run clusterProfiler KEGG.
clhgo <- compareCluster(
  geneCluster = cl,
  fun = enrichGO,
  OrgDb = org.Hs.eg.db,
  ont = "CC",
  keyType = "ENTREZID"
)
clhgo <- setReadable(clhgo, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
head(clhgo)

clkk_cc <- compareCluster(
  geneCluster = clk,
  fun = enrichGO,
  OrgDb = org.Hs.eg.db,
  ont = "CC",
  keyType = "ENTREZID"
)

clkk_cc <- setReadable(clkk_cc, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

head(clkk_cc)

# make dotplot.

dp1 <- dotplot(clhgo, showCategory = 10, title = "Hierarchical clustering")
dp2 <- dotplot(clkk_cc, showCategory = 10, title = "K-means clustering")

# save the two dotplots arranged side by side.
install.packages("gridExtra")
library(gridExtra)
pdf(sprintf("/Users/satabdi/Downloads/%s_dp_ko_kegg.pdf", Sys.Date()), width = 20, height = 15)
grid.arrange(dp1, dp2, ncol = 2)
dev.off()


clhgo_df <- as.data.frame(clhgo)
clkkgo_df <- as.data.frame(clkk_cc)


top4_genes_clh <- clhgo_df %>%
  arrange(p.adjust) %>% 
  separate_rows(geneID, sep = "/") %>%
  group_by(Cluster, Description) %>%
  slice_head(n = 4) %>%
  select(Cluster,
         Cellular_Component = Description,
         Gene = geneID,
         p.adjust)

top4_genes_clh

top4_genes_clk <- clkkgo_df %>%
  arrange(p.adjust) %>%
  separate_rows(geneID, sep = "/") %>%
  group_by(Cluster, Description) %>%
  slice_head(n = 4) %>%
  select(Cluster,
         Cellular_Component = Description,
         Gene = geneID,
         p.adjust)

top4_genes_clk

write.csv(top4_genes_clh,
          "/Users/satabdi/Downloads/top4_genes_CC_clh_KO.csv",
          row.names = FALSE)

write.csv(top4_genes_clk,
          "/Users/satabdi/Downloads/top4_genes_CC_clk_KO.csv",
          row.names = FALSE)

top4_CC_clh <- clhgo_df %>%
  group_by(Cluster) %>%
  arrange(p.adjust, .by_group = TRUE) %>%
  slice_head(n = 4) %>%
  select(Cluster,
         Cellular_Component = Description,
         p.adjust,
         GeneRatio,
         Count)

top4_CC_clh

top4_CC_clk <- clkkgo_df %>%
  group_by(Cluster) %>%
  arrange(p.adjust, .by_group = TRUE) %>%
  slice_head(n = 4) %>%
  select(Cluster,
         Cellular_Component = Description,
         p.adjust,
         GeneRatio,
         Count)

top4_CC_clk

write.csv(top4_CC_clh,
          "/Users/satabdi/Downloads/top4_CC_per_cluster_clh_KO.csv",
          row.names = FALSE)

write.csv(top4_CC_clk,
          "/Users/satabdi/Downloads/top4_CC_per_cluster_clk_KO.csv",
          row.names = FALSE)