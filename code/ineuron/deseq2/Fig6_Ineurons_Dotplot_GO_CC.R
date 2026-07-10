# ==================================================
#  GO Cellular Component enrichment 
# ==================================================

# ==================================================
#  Load Libraries
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
# Read Data
# ==================================================

clh_raw <- read.delim("/Users/satabdi/Downloads/cluster_clh.txt")
clk_raw <- read.delim("/Users/satabdi/Downloads/cluster_clk.txt")

# check
head(clh_raw$gene)

# ==================================================
#  Split Clusters
# ==================================================

cl_list <- clh_raw %>%
  group_split(cluster) %>%
  setNames(paste0("Cluster", sort(unique(clh_raw$cluster)))) %>%
  lapply(function(df) unique(df$gene))

clk_list <- clk_raw %>%
  group_split(cluster) %>%
  setNames(paste0("Cluster", sort(unique(clk_raw$cluster)))) %>%
  lapply(function(df) unique(df$gene))

# confirm before mapping
str(cl_list)

# ==================================================
# Map ENSEMBL ID to ENTREZ ID
# ==================================================

map_genes_to_entrez <- function(gene_list) {
  if (length(gene_list) == 0) return(character(0))
  
  entrez <- mapIds(
    org.Hs.eg.db,
    keys = gene_list,
    column = "ENTREZID",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
  
  entrez <- unlist(entrez)
  entrez <- entrez[!is.na(entrez)]
  unique(entrez)
}

cl_list  <- lapply(cl_list, map_genes_to_entrez)
clk_list <- lapply(clk_list, map_genes_to_entrez)

# ==================================================
#  GO Cellular Component Enrichment
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

# Keep ENTREZ-based data frames for joining/ranking
clh_cc_df_entrez <- as.data.frame(clhgo_cc) %>% filter(!is.na(p.adjust))
clk_cc_df_entrez <- as.data.frame(clkgo_cc) %>% filter(!is.na(p.adjust))

# ==================================================
# Load DEG results table for per-gene ranking
# ==================================================

degs_table <- read.delim("/Users/satabdi/Downloads/combined_restables.txt", row.names = 1)
degs_table$ensembl_id <- rownames(degs_table)

# identify all padj columns 
padj_cols <- grep("^padj_", colnames(degs_table), value = TRUE)

# ==================================================
# Convert a vector of Entrez IDs to Ensembl IDs
# ==================================================

entrez_to_ensembl <- function(entrez_ids) {
  mapIds(org.Hs.eg.db, keys = entrez_ids, column = "ENSEMBL",
         keytype = "ENTREZID", multiVals = "first")
}

# ==================================================
# Convert a vector of Ensembl IDs to gene symbols
# ==================================================

ensembl_to_symbol <- function(ensembl_ids) {
  mapIds(org.Hs.eg.db, keys = ensembl_ids, column = "SYMBOL",
         keytype = "ENSEMBL", multiVals = "first")
}


# ==================================================
# Get Top 4 Cellular Components per Cluster (p.adjust < 0.05)
# ==================================================

get_top4_CC <- function(df) {
  df %>%
    dplyr::filter(p.adjust < 0.05) %>%              
    dplyr::group_by(Cluster) %>%
    dplyr::arrange(p.adjust, .by_group = TRUE) %>%
    dplyr::slice_head(n = 4) %>%
    dplyr::ungroup() %>%
    dplyr::select(Cluster, Description, p.adjust, GeneRatio, BgRatio, Count)
}

clh_top4_CC <- get_top4_CC(clh_cc_df_entrez)
clk_top4_CC <- get_top4_CC(clk_cc_df_entrez)

# ==================================================
# Get Top 4 Genes per Top 4 CC per Cluster (p.adjust < 0.05)
# ==================================================

get_top4_genes_per_CC_ranked <- function(df, degs_table, padj_cols) {
  df %>%
    dplyr::filter(p.adjust < 0.05) %>%             
    dplyr::group_by(Cluster) %>%
    dplyr::arrange(p.adjust, .by_group = TRUE) %>%
    dplyr::slice_head(n = 4) %>%
    dplyr::ungroup() %>%
    tidyr::separate_rows(geneID, sep = "/") %>%
    dplyr::rename(entrez_id = geneID) %>%
    dplyr::mutate(ensembl_id = entrez_to_ensembl(entrez_id)) %>%
    dplyr::left_join(
      degs_table %>% dplyr::select(ensembl_id, dplyr::all_of(padj_cols)),
      by = "ensembl_id"
    ) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(min_padj = suppressWarnings(min(dplyr::c_across(dplyr::all_of(padj_cols)), na.rm = TRUE))) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Cluster, Description) %>%
    dplyr::arrange(min_padj, .by_group = TRUE) %>%
    dplyr::slice_head(n = 4) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(symbol = ensembl_to_symbol(ensembl_id)) %>%
    dplyr::select(Cluster, Description, ensembl_id, symbol, min_padj, p.adjust)
}

clh_top4_genes <- get_top4_genes_per_CC_ranked(clh_cc_df_entrez, degs_table, padj_cols)
clk_top4_genes <- get_top4_genes_per_CC_ranked(clk_cc_df_entrez, degs_table, padj_cols)
# ==================================================
# Attach gene symbols to the CC summary tables 
# ==================================================

add_symbol_list <- function(cc_df, full_entrez_df) {
  cc_df %>%
    dplyr::left_join(
      full_entrez_df %>%
        dplyr::select(Cluster, Description, geneID) %>%
        dplyr::rowwise() %>%
        dplyr::mutate(
          symbol_list = paste(
            ensembl_to_symbol(entrez_to_ensembl(str_split(geneID, "/")[[1]])),
            collapse = "/"
          )
        ) %>%
        dplyr::ungroup() %>%
        dplyr::select(Cluster, Description, symbol_list),
      by = c("Cluster", "Description")
    )
}

clh_top4_CC <- add_symbol_list(clh_top4_CC, clh_cc_df_entrez)
clk_top4_CC <- add_symbol_list(clk_top4_CC, clk_cc_df_entrez)

# ==================================================
# Convert to symbol-readable versions for the dotplots
# ==================================================

clhgo_cc_readable <- setReadable(clhgo_cc, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
clkgo_cc_readable <- setReadable(clkgo_cc, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")

# ==================================================
# Dotplots
# ==================================================

dpgo_cc1 <- dotplot(
  clhgo_cc_readable,
  showCategory = 10,
  font.size = 12,
  title = "Hierarchical Clustering \u2013 GO Cellular Component"
)

dpgo_cc2 <- dotplot(
  clkgo_cc_readable,
  showCategory = 10,
  font.size = 12,
  title = "K-means Clustering \u2013 GO Cellular Component"
)

# ==================================================
#  Save Figure
# ==================================================

pdf(sprintf(
  "/Users/satabdi/Downloads/%s_GO_CC_Dotplots.pdf",
  Sys.Date()
), width = 20, height = 15)

grid.arrange(dpgo_cc1, dpgo_cc2, ncol = 2)
dev.off()

# ==================================================
# Save Results Tables (Ensembl ID retained + symbol attached)
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

