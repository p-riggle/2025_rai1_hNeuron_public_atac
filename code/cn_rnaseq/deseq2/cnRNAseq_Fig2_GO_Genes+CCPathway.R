# ==================================================
# 1️⃣ Load Libraries
# ==================================================
# Install BiocManager if not installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# Install Bioconductor packages
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("AnnotationDbi")

# Install CRAN packages
install.packages("dplyr")
install.packages("tidyr")



library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(tidyr)
library(AnnotationDbi)

# ==================================================
# 2️⃣ Import Cluster Data
# ==================================================

clh <- read.delim("/Users/satabdi/Downloads/mh_het.txt")

# ==================================================
# 3️⃣ Extract Genes for Each Cluster
# ==================================================

clusters <- split(clh$gene, clh$cluster)

# Name clusters
names(clusters) <- paste0("Cluster", names(clusters))

# ==================================================
# 4️⃣ Convert Gene Symbols → Entrez IDs
# ==================================================

clusters_entrez <- lapply(clusters, function(x) {
  
  entrez <- mapIds(org.Hs.eg.db,
                   keys = x,
                   column = "ENTREZID",
                   keytype = "SYMBOL",
                   multiVals = "first")
  
  entrez <- na.omit(entrez)
  return(entrez)
  
})

# ==================================================
# 5️⃣ Run GO Enrichment (Cellular Component)
# ==================================================

go_cc <- compareCluster(
  geneCluster = clusters_entrez,
  fun = enrichGO,
  OrgDb = org.Hs.eg.db,
  ont = "CC",
  keyType = "ENTREZID",
  pAdjustMethod = "BH"
)

# Convert Entrez → Gene Symbols
go_cc <- setReadable(go_cc,
                     OrgDb = org.Hs.eg.db,
                     keyType = "ENTREZID")

# ==================================================
# 6️⃣ Convert Results to Dataframe
# ==================================================

cc_results <- as.data.frame(go_cc)

# ==================================================
# 7️⃣ Extract Top 4 Significant Genes per Cellular Component
# ==================================================

top4_genes_cc <- cc_results %>%
  arrange(p.adjust) %>%
  separate_rows(geneID, sep = "/") %>%
  group_by(Cluster, Description) %>%
  slice_head(n = 4) %>%
  ungroup() %>%
  select(Cluster, Cellular_Component = Description, Gene = geneID, p.adjust)

# ==================================================
# 8️⃣ Save Results
# ==================================================

write.csv(top4_genes_cc,
          "/Users/satabdi/Downloads/top4_genes_per_CC.csv",
          row.names = FALSE)

# ==================================================
# 9️⃣ Optional: Visualize GO Results
# ==================================================

dotplot(go_cc, showCategory = 10,
        title = "GO Cellular Component Enrichment")