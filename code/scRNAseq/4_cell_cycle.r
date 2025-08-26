
# https://satijalab.org/seurat/articles/cell_cycle_vignette.html

# Load libraries
library(Seurat)
library(BPCells)
library(tidyverse)

# load the R.data
geo_so = readRDS('data/scRNAseq/iwaselab/rdata/geo_so_sct_unintegrated_with_markers_0.1.rds')

# human CC genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# cellcycle scoring
sctSeuFilt <- CellCycleScoring(geo_so, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

cell_cycle_scores <- sctSeuFilt@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(
    avg_S_Score = mean(S.Score),
    avg_G2M_Score = mean(G2M.Score)
  )

install.packages("myPlottingPackage")
library(myPlottingPackage)

main_violin_s.score <- VlnPlot(sctSeuFilt, features = "S.Score", group.by = "seurat_clusters", pt.size = 0)
main_violin_s.score

ggsave("results/scRNAseq/4_cellcycle/main_violin_s.score.pdf", main_violin_s.score, width = 8, height = 6)

main_violin_g2m.score <- VlnPlot(sctSeuFilt, features = "G2M.Score", group.by = "seurat_clusters", pt.size = 0)
main_violin_g2m.score

ggsave("results/scRNAseq/4_cellcycle/main_violin_g2m.score.pdf", main_violin_g2m.score, width = 8, height = 6)

main_dim <- DimPlot(sctSeuFilt, group.by = "Phase")
main_dim

sctSeuFilt@meta.data
sctSeuFilt$genotype <- factor(sctSeuFilt$genotype, levels = c("wt", "het", "ko"))
main_dim_genotype <- DimPlot(sctSeuFilt, group.by = "Phase", split.by = "genotype")
main_dim_genotype

ggsave("results/scRNAseq/4_cellcycle/main_dim.pdf", main_dim, width = 5, height = 5)
ggsave("results/scRNAseq/4_cellcycle/main_dim_genotype.pdf", main_dim_genotype, width = 12, height = 5)

