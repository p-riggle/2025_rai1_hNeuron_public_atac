# ==================================================
# For KO samples
# GO Cellular Component enrichment
# ==================================================

# install cluterProfiler
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("clusterProfiler")

library(clusterProfiler)
library(dplyr)
library(tidyr)
library(org.Hs.eg.db)
library(gridExtra)

# ==================================================
# Import cluster data
# ==================================================

clh <- read.delim("/Users/satabdi/Downloads/mh_ko.txt")
clk <- read.delim("/Users/satabdi/Downloads/mk_ko.txt")
str(clh)
str(clk)

# ==================================================
# Extract unique genes from each cluster
# ==================================================

cl1 <- clh[clh$cluster == 1,] %>% dplyr::select(gene) %>% unlist()
cl2 <- clh[clh$cluster == 2,] %>% dplyr::select(gene) %>% unlist()
cl3 <- clh[clh$cluster == 3,] %>% dplyr::select(gene) %>% unlist()
cl4 <- clh[clh$cluster == 4,] %>% dplyr::select(gene) %>% unlist()
cl5 <- clh[clh$cluster == 5,] %>% dplyr::select(gene) %>% unlist()
cl6 <- clh[clh$cluster == 6,] %>% dplyr::select(gene) %>% unlist()

clk1 <- clk[clk$cluster == 1,] %>% dplyr::select(gene) %>% unlist()
clk2 <- clk[clk$cluster == 2,] %>% dplyr::select(gene) %>% unlist()
clk3 <- clk[clk$cluster == 3,] %>% dplyr::select(gene) %>% unlist()
clk4 <- clk[clk$cluster == 4,] %>% dplyr::select(gene) %>% unlist()
clk5 <- clk[clk$cluster == 5,] %>% dplyr::select(gene) %>% unlist()
clk6 <- clk[clk$cluster == 6,] %>% dplyr::select(gene) %>% unlist()
clk7 <- clk[clk$cluster == 7,] %>% dplyr::select(gene) %>% unlist()

cl <- list(cl1, cl2, cl3, cl4, cl5, cl6)
clk <- list(clk1, clk2, clk3, clk4, clk5, clk6, clk7)

names(cl)  <- c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Cluster5", "Cluster6")
names(clk) <- c("Cluster1", "Cluster2", "Cluster3", "Cluster4", "Cluster5", "Cluster6", "Cluster7")

# ==================================================
# Convert gene symbols -> Entrez IDs
# ==================================================

cl <- lapply(cl, function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(character(0))
  
  entrez <- mapIds(org.Hs.eg.db,
                   keys = x,
                   column = "ENTREZID",
                   keytype = "SYMBOL",
                   multiVals = "first")
  
  entrez <- unlist(entrez)
  entrez <- entrez[!is.na(entrez)]
  unique(entrez)
})

clk <- lapply(clk, function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(character(0))
  
  entrez <- mapIds(org.Hs.eg.db,
                   keys = x,
                   column = "ENTREZID",
                   keytype = "SYMBOL",
                   multiVals = "first")
  
  entrez <- unlist(entrez)
  entrez <- entrez[!is.na(entrez)]
  unique(entrez)
})

str(cl)
str(clk)

# ==================================================
# Run GO Cellular Component enrichment
# ==================================================

clhgo <- compareCluster(
  geneCluster   = cl,
  fun           = enrichGO,
  OrgDb         = org.Hs.eg.db,
  ont           = "CC",
  keyType       = "ENTREZID",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  minGSSize     = 5
)
clhgo <- setReadable(clhgo, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
head(clhgo)

clkk_cc <- compareCluster(
  geneCluster   = clk,
  fun           = enrichGO,
  OrgDb         = org.Hs.eg.db,
  ont           = "CC",
  keyType       = "ENTREZID",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  minGSSize     = 5
)
clkk_cc <- setReadable(clkk_cc, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
head(clkk_cc)

# ==================================================
# Dotplots
# ==================================================

dp1 <- dotplot(clhgo, showCategory = 10, title = "Hierarchical clustering")
dp2 <- dotplot(clkk_cc, showCategory = 10, title = "K-means clustering")

pdf(sprintf("/Users/satabdi/Downloads/%s_dp_ko_kegg.pdf", Sys.Date()), width = 20, height = 15)
grid.arrange(dp1, dp2, ncol = 2)
dev.off()

# ==================================================
# Convert results to data frames
# ==================================================

clhgo_df  <- as.data.frame(clhgo)  %>% filter(!is.na(p.adjust))
clkkgo_df <- as.data.frame(clkk_cc) %>% filter(!is.na(p.adjust))

# ==================================================
# Top 4 Cellular Components per Cluster (p.adjust < 0.05)
# ==================================================

top4_CC_clh <- clhgo_df %>%
  filter(p.adjust < 0.05) %>%
  group_by(Cluster) %>%
  arrange(p.adjust, .by_group = TRUE) %>%
  slice_head(n = 4) %>%
  ungroup() %>%
  dplyr::select(Cluster,
                Cellular_Component = Description,
                p.adjust,
                GeneRatio,
                Count)

top4_CC_clh

top4_CC_clk <- clkkgo_df %>%
  filter(p.adjust < 0.05) %>%
  group_by(Cluster) %>%
  arrange(p.adjust, .by_group = TRUE) %>%
  slice_head(n = 4) %>%
  ungroup() %>%
  dplyr::select(Cluster,
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

# ==================================================
# Top 4 Genes per Top 4 CC per Cluster (p.adjust < 0.05)
# ==================================================

top4_genes_clh <- clhgo_df %>%
  filter(p.adjust < 0.05) %>%
  group_by(Cluster) %>%
  arrange(p.adjust, .by_group = TRUE) %>%
  slice_head(n = 4) %>%
  ungroup() %>%
  separate_rows(geneID, sep = "/") %>%
  group_by(Cluster, Description) %>%
  slice_head(n = 4) %>%
  ungroup() %>%
  dplyr::select(Cluster,
                Cellular_Component = Description,
                Gene = geneID,
                p.adjust)

top4_genes_clh

top4_genes_clk <- clkkgo_df %>%
  filter(p.adjust < 0.05) %>%
  group_by(Cluster) %>%
  arrange(p.adjust, .by_group = TRUE) %>%
  slice_head(n = 4) %>%
  ungroup() %>%
  separate_rows(geneID, sep = "/") %>%
  group_by(Cluster, Description) %>%
  slice_head(n = 4) %>%
  ungroup() %>%
  dplyr::select(Cluster,
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