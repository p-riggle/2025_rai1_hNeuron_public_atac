# ==============================
# Load libraries
# ==============================

library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# ==============================
# Load datasets
# ==============================

astro_het <- read.csv("/Users/satabdi/Downloads/rstbl_Astro_wt.vs.het.csv")
astro_ko  <- read.csv("/Users/satabdi/Downloads/rstbl_Astro_wt.vs.ko.csv")

neuron_het <- read.csv("/Users/satabdi/Downloads/rstbl_Neuron_wt.vs.het.csv")
neuron_ko  <- read.csv("/Users/satabdi/Downloads/rstbl_Neuron_wt.vs.ko.csv")

# ==============================
# Fix gene column
# ==============================

astro_het  <- astro_het  %>% dplyr::rename(gene = X)
astro_ko   <- astro_ko   %>% dplyr::rename(gene = X)

neuron_het <- neuron_het %>% dplyr::rename(gene = X)
neuron_ko  <- neuron_ko  %>% dplyr::rename(gene = X)


# ==============================
# Filter significant genes
# ==============================

astro_deg_het <- astro_het %>% filter(padj < 0.05 & abs(log2FoldChange) > 1)
astro_deg_ko  <- astro_ko  %>% filter(padj < 0.05 & abs(log2FoldChange) > 1)

neuron_deg_het <- neuron_het %>% filter(padj < 0.05 & abs(log2FoldChange) > 1)
neuron_deg_ko  <- neuron_ko  %>% filter(padj < 0.05 & abs(log2FoldChange) > 1)

# ==============================
# Separate UP and DOWN genes
# ==============================

astro_up_het   <- astro_deg_het %>% filter(log2FoldChange > 1) %>% pull(gene)
astro_down_het <- astro_deg_het %>% filter(log2FoldChange < -1) %>% pull(gene)

astro_up_ko   <- astro_deg_ko %>% filter(log2FoldChange > 1) %>% pull(gene)
astro_down_ko <- astro_deg_ko %>% filter(log2FoldChange < -1) %>% pull(gene)

neuron_up_het   <- neuron_deg_het %>% filter(log2FoldChange > 1) %>% pull(gene)
neuron_down_het <- neuron_deg_het %>% filter(log2FoldChange < -1) %>% pull(gene)

neuron_up_ko   <- neuron_deg_ko %>% filter(log2FoldChange > 1) %>% pull(gene)
neuron_down_ko <- neuron_deg_ko %>% filter(log2FoldChange < -1) %>% pull(gene)

# ==============================
# Convert gene IDs
# ==============================

convert_ids <- function(genes) {
  bitr(genes,
       fromType = "SYMBOL",
       toType   = "ENTREZID",
       OrgDb    = org.Hs.eg.db) %>%
    pull(ENTREZID)
}

astro_up_het_ids   <- convert_ids(astro_up_het)
astro_down_het_ids <- convert_ids(astro_down_het)

astro_up_ko_ids   <- convert_ids(astro_up_ko)
astro_down_ko_ids <- convert_ids(astro_down_ko)

neuron_up_het_ids   <- convert_ids(neuron_up_het)
neuron_down_het_ids <- convert_ids(neuron_down_het)

neuron_up_ko_ids   <- convert_ids(neuron_up_ko)
neuron_down_ko_ids <- convert_ids(neuron_down_ko)

# ==============================
# GO Cellular Component enrichment
# ==============================

run_CC <- function(ids) {
  enrichGO(
    gene          = ids,
    OrgDb         = org.Hs.eg.db,
    ont           = "CC",
    pAdjustMethod = "BH",
    qvalueCutoff  = 0.05,
    readable      = TRUE
  )
}

astro_up_het_CC   <- run_CC(astro_up_het_ids)
astro_down_het_CC <- run_CC(astro_down_het_ids)

astro_up_ko_CC   <- run_CC(astro_up_ko_ids)
astro_down_ko_CC <- run_CC(astro_down_ko_ids)

neuron_up_het_CC   <- run_CC(neuron_up_het_ids)
neuron_down_het_CC <- run_CC(neuron_down_het_ids)

neuron_up_ko_CC   <- run_CC(neuron_up_ko_ids)
neuron_down_ko_CC <- run_CC(neuron_down_ko_ids)

# ==============================
# Plot results
# ==============================

pdf("/Users/satabdi/Downloads/GO_CC_Astro_Neuron.pdf", width=10, height=8)

dotplot(astro_up_het_CC, showCategory=15) +
  ggtitle("Astrocyte HET – Upregulated Cellular Components")

dotplot(astro_down_het_CC, showCategory=15) +
  ggtitle("Astrocyte HET – Downregulated Cellular Components")

dotplot(astro_up_ko_CC, showCategory=15) +
  ggtitle("Astrocyte KO – Upregulated Cellular Components")

dotplot(astro_down_ko_CC, showCategory=15) +
  ggtitle("Astrocyte KO – Downregulated Cellular Components")

dotplot(neuron_up_het_CC, showCategory=15) +
  ggtitle("Neuron HET – Upregulated Cellular Components")

dotplot(neuron_down_het_CC, showCategory=15) +
  ggtitle("Neuron HET – Downregulated Cellular Components")

dotplot(neuron_up_ko_CC, showCategory=15) +
  ggtitle("Neuron KO – Upregulated Cellular Components")

dotplot(neuron_down_ko_CC, showCategory=15) +
  ggtitle("Neuron KO – Downregulated Cellular Components")

dev.off()

#Extracting top 4 genes per CC based on padj
get_top_genes_per_CC <- function(enrich_result, deg_table, top_n = 4) {
  
  cc_table <- enrich_result@result
  
  gene_list <- cc_table %>%
    dplyr::select(ID, Description, geneID)
  
  gene_expanded <- gene_list %>%
    tidyr::separate_rows(geneID, sep = "/")
  
  gene_expanded <- gene_expanded %>%
    dplyr::rename(gene = geneID)
  
  merged <- gene_expanded %>%
    dplyr::left_join(deg_table, by = "gene")
  
  top_genes <- merged %>%
    dplyr::group_by(Description) %>%
    dplyr::arrange(padj) %>%
    dplyr::slice_head(n = top_n) %>%
    dplyr::ungroup()
  
  return(top_genes)
}


astro_up_het_top <- get_top_genes_per_CC(astro_up_het_CC, astro_deg_het)
astro_up_ko_top  <- get_top_genes_per_CC(astro_up_ko_CC, astro_deg_ko)

neuron_up_het_top <- get_top_genes_per_CC(neuron_up_het_CC, neuron_deg_het)
neuron_up_ko_top  <- get_top_genes_per_CC(neuron_up_ko_CC, neuron_deg_ko)

astro_down_het_top <- get_top_genes_per_CC(astro_down_het_CC, astro_deg_het)
astro_down_ko_top  <- get_top_genes_per_CC(astro_down_ko_CC, astro_deg_ko)

neuron_down_het_top <- get_top_genes_per_CC(neuron_down_het_CC, neuron_deg_het)
neuron_down_ko_top  <- get_top_genes_per_CC(neuron_down_ko_CC, neuron_deg_ko)

#Save the results
write.csv(astro_up_het_top,
          "/Users/satabdi/Downloads/Astro_HET_up_top4_genes_per_CC.csv",
          row.names = FALSE)

write.csv(astro_down_het_top,
          "/Users/satabdi/Downloads/Astro_HET_down_top4_genes_per_CC.csv",
          row.names = FALSE)

write.csv(astro_up_ko_top,
          "/Users/satabdi/Downloads/Astro_KO_up_top4_genes_per_CC.csv",
          row.names = FALSE)

write.csv(astro_down_ko_top,
          "/Users/satabdi/Downloads/Astro_KO_down_top4_genes_per_CC.csv",
          row.names = FALSE)

write.csv(neuron_up_het_top,
          "/Users/satabdi/Downloads/Neuron_HET_up_top4_genes_per_CC.csv",
          row.names = FALSE)

write.csv(neuron_down_het_top,
          "/Users/satabdi/Downloads/Neuron_HET_down_top4_genes_per_CC.csv",
          row.names = FALSE)

write.csv(neuron_up_ko_top,
          "/Users/satabdi/Downloads/Neuron_KO_up_top4_genes_per_CC.csv",
          row.names = FALSE)

write.csv(neuron_down_ko_top,
          "/Users/satabdi/Downloads/Neuron_KO_down_top4_genes_per_CC.csv",
          row.names = FALSE)

# ==========================================
# Function to extract Top 4 CC terms
# ==========================================

get_top_CC <- function(enrich_result){
  
  enrich_result@result %>%
    arrange(p.adjust) %>%
    slice_head(n = 4)
  
}

# ==========================================
# Extract Top 4 Cellular Components
# ==========================================

astro_het_up_topCC   <- get_top_CC(astro_up_het_CC)
astro_het_down_topCC <- get_top_CC(astro_down_het_CC)

astro_ko_up_topCC   <- get_top_CC(astro_up_ko_CC)
astro_ko_down_topCC <- get_top_CC(astro_down_ko_CC)

neuron_het_up_topCC   <- get_top_CC(neuron_up_het_CC)
neuron_het_down_topCC <- get_top_CC(neuron_down_het_CC)

neuron_ko_up_topCC   <- get_top_CC(neuron_up_ko_CC)
neuron_ko_down_topCC <- get_top_CC(neuron_down_ko_CC)

# ==========================================
# View results
# ==========================================

astro_het_up_topCC
astro_het_down_topCC

astro_ko_up_topCC
astro_ko_down_topCC

neuron_het_up_topCC
neuron_het_down_topCC

neuron_ko_up_topCC
neuron_ko_down_topCC

# ==========================================
# Save results
# ==========================================

write.csv(astro_het_up_topCC, "/Users/satabdi/Downloads/Astro_HET_Up_Top4_CC.csv", row.names=FALSE)
write.csv(astro_het_down_topCC, "/Users/satabdi/Downloads/Astro_HET_Down_Top4_CC.csv", row.names=FALSE)

write.csv(astro_ko_up_topCC, "/Users/satabdi/Downloads/Astro_KO_Up_Top4_CC.csv", row.names=FALSE)
write.csv(astro_ko_down_topCC, "/Users/satabdi/Downloads/Astro_KO_Down_Top4_CC.csv", row.names=FALSE)

write.csv(neuron_het_up_topCC, "/Users/satabdi/Downloads/Neuron_HET_Up_Top4_CC.csv", row.names=FALSE)
write.csv(neuron_het_down_topCC, "/Users/satabdi/Downloads/Neuron_HET_Down_Top4_CC.csv", row.names=FALSE)

write.csv(neuron_ko_up_topCC, "/Users/satabdi/Downloads/Neuron_KO_Up_Top4_CC.csv", row.names=FALSE)
write.csv(neuron_ko_down_topCC, "/Users/satabdi/Downloads/Neuron_KO_Down_Top4_CC.csv", row.names=FALSE)

