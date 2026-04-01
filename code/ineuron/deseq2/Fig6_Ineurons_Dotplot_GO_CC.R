# ==================================================
# 0️⃣ Load Libraries
# ==================================================

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(gridExtra)
library(dplyr)
library(tidyr)
library(stringr)

# ==================================================
# 1️⃣ Read Data
# ==================================================

clh_raw <- read.delim("/Users/satabdi/Downloads/cluster_clh.txt")
clk_raw <- read.delim("/Users/satabdi/Downloads/cluster_clk.txt")

# ==================================================
# 2️⃣ Split Clusters
# ==================================================

cl_list <- clh_raw %>%
  group_split(cluster) %>%
  setNames(paste0("Cluster", unique(clh_raw$cluster))) %>%
  lapply(function(df) df$gene)

clk_list <- clk_raw %>%
  group_split(cluster) %>%
  setNames(paste0("Cluster", unique(clk_raw$cluster))) %>%
  lapply(function(df) df$gene)

# ==================================================
# 3️⃣ Map SYMBOL → ENTREZ
# ==================================================

map_genes_to_entrez <- function(gene_list) {
  if (length(gene_list) == 0) return(character(0))
  
  entrez <- mapIds(
    org.Hs.eg.db,
    keys = gene_list,
    column = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
  )
  
  entrez <- unlist(entrez)
  entrez <- entrez[!is.na(entrez)]
  unique(entrez)
}

cl_list  <- lapply(cl_list, map_genes_to_entrez)
clk_list <- lapply(clk_list, map_genes_to_entrez)

# Remove clusters with very few genes
cl_list  <- cl_list[sapply(cl_list, length) > 5]
clk_list <- clk_list[sapply(clk_list, length) > 5]

# ==================================================
# 4️⃣ GO Cellular Component Enrichment
# ==================================================

clhgo_cc <- compareCluster(
  geneCluster   = cl_list,
  fun           = enrichGO,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  minGSSize     = 5
)

clhgo_cc <- setReadable(clhgo_cc, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

clkgo_cc <- compareCluster(
  geneCluster   = clk_list,
  fun           = enrichGO,
  OrgDb         = org.Hs.eg.db,
  keyType       = "ENTREZID",
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  minGSSize     = 5
)

clkgo_cc <- setReadable(clkgo_cc, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# ==================================================
# 5️⃣ Convert Results to Data Frame
# ==================================================

clh_cc_df <- as.data.frame(clhgo_cc) %>%
  filter(!is.na(p.adjust))

clk_cc_df <- as.data.frame(clkgo_cc) %>%
  filter(!is.na(p.adjust))


# ==================================================
# 6️⃣ Get Top 4 Cellular Components per Cluster
# ==================================================
get_top4_CC <- function(df) {
  df %>%
    dplyr::group_by(Cluster) %>%
    dplyr::arrange(p.adjust, .by_group = TRUE) %>%
    dplyr::slice_head(n = 4) %>%
    dplyr::ungroup() %>%
    dplyr::select(Cluster,
                  Description,
                  p.adjust,
                  GeneRatio,
                  BgRatio,
                  Count)
}

clh_top4_CC <- get_top4_CC(clh_cc_df)
clk_top4_CC <- get_top4_CC(clk_cc_df)

# ==================================================
# 7️⃣ Get Top 4 Genes per Top 4 CC per Cluster
# ==================================================

get_top4_genes_per_CC <- function(df) {
  
  df %>%
    
    dplyr::group_by(Cluster) %>%
    dplyr::arrange(p.adjust, .by_group = TRUE) %>%
    dplyr::slice_head(n = 4) %>%
    
    
    tidyr::separate_rows(geneID, sep = "/") %>%
    
    
    dplyr::group_by(Cluster, Description) %>%
    dplyr::slice_head(n = 4) %>%
    
    dplyr::ungroup() %>%
    
    
    dplyr::select(Cluster,
                  Description,
                  geneID,
                  p.adjust)
}


clh_top4_genes <- get_top4_genes_per_CC(clh_cc_df)
clk_top4_genes <- get_top4_genes_per_CC(clk_cc_df)
# ==================================================
# 8️⃣ Dotplots
# ==================================================

dpgo_cc1 <- dotplot(
  clhgo_cc,
  showCategory = 10,
  font.size = 12,
  title = "Hierarchical Clustering – GO Cellular Component"
)

dpgo_cc2 <- dotplot(
  clkgo_cc,
  showCategory = 10,
  font.size = 12,
  title = "K-means Clustering – GO Cellular Component"
)

# ==================================================
# 9️⃣ Save Figure
# ==================================================

pdf(sprintf(
  "/Users/satabdi/Downloads/%s_GO_CC_Dotplots.pdf",
  Sys.Date()
), width = 20, height = 15)

grid.arrange(dpgo_cc1, dpgo_cc2, ncol = 2)
dev.off()

# ==================================================
# 🔟 Save Results Tables
# ==================================================

write.csv(clh_top4_CC,
          "/Users/satabdi/Downloads/CLH_Top_CellComponent_perCluster.csv",
          row.names = FALSE)

write.csv(clk_top4_CC,
          "/Users/satabdi/Downloads/CLK_Top_CellComponent_perCluster.csv",
          row.names = FALSE)

write.csv(clh_top4_genes,
          "/Users/satabdi/Downloads/CLH_Top4_Genes_per_CC.csv",
          row.names = FALSE)

write.csv(clk_top4_genes,
          "/Users/satabdi/Downloads/CLK_Top4_Genes_per_CC.csv",
          row.names = FALSE)

