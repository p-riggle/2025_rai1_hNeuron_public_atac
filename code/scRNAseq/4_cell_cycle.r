# set working directory
setwd("/nfs/turbo/umms-siwase/SHIGEKI2/rai1_hesc")

# https://satijalab.org/seurat/articles/cell_cycle_vignette.html

# Load libraries
library(Seurat)
library(BPCells)
library(tidyverse)

# load the R.data
geo_so = readRDS('data/scRNAseq/iwaselab/rdata/geo_so_sct_unintegrated_with_annos.rds')

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

#########################################################
### umap plot split by gentype ###########################
#########################################################

main_dim_genotype <- DimPlot(
    sctSeuFilt, 
    group.by = "Phase", 
    split.by = "genotype", 
    label = FALSE,
    raster = TRUE,            # Keep dots rasterized
    raster.dpi = c(300, 300)
    ) +
  labs(x = "UMAP_1", y = "UMAP_2") +
  theme(
    text = element_text(family = "sans", color = "black")
  ) 

main_dim_genotype

ggsave("results/scRNAseq/5_cellcycle/main_dim.pdf", main_dim, width = 5, height = 5)
ggsave("results/scRNAseq/5_cellcycle/main_dim_genotype.pdf", main_dim_genotype, width = 12, height = 5)

####################################################
#### make stack bar plot for cell cycle phases ########
#######################################################

# make stack bar plots to compare the genotypes for both day 20 and day 60
cell_cycle_df <- sctSeuFilt@meta.data %>%
  select(genotype, Phase, day) %>%
  count(genotype, Phase, day, name = "cell_count") %>%
  group_by(genotype, day) %>%
  mutate(percentage = cell_count / sum(cell_count) * 100) %>%
  ungroup()

# Create stacked bar plot with cell counts
stack_bar_count <- ggplot(cell_cycle_df, aes(x = genotype, y = cell_count, fill = Phase)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ day, scales = "free_y") +
  scale_fill_manual(values = c("G1" = "#ff6b6b", "S" = "#4dabf7", "G2M" = "#51cf66")) +
  labs(title = "Cell Cycle Phase Distribution by Genotype",
       x = "Genotype", 
       y = "Cell Count",
       fill = "Cell Cycle Phase") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

# Create stacked bar plot with percentages
stack_bar_percent <- ggplot(cell_cycle_df, aes(x = genotype, y = percentage, fill = Phase)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ day) +
  scale_fill_manual(values = c("G1" = "#ff6b6b", "S" = "#4dabf7", "G2M" = "#51cf66")) +
  labs(title = "Cell Cycle Phase Distribution by Genotype (Percentage)",
       x = "Genotype", 
       y = "Percentage of Cells",
       fill = "Cell Cycle Phase") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

# Save the plots
ggsave("results/scRNAseq/5_cellcycle/stack_bar_count.pdf", 
       plot = stack_bar_count, width = 10, height = 6)
ggsave("results/scRNAseq/5_cellcycle/stack_bar_percent.pdf", 
       plot = stack_bar_percent, width = 10, height = 6)

# Also create a single plot combining both days
stack_bar_combined <- ggplot(cell_cycle_df, aes(x = interaction(genotype, day), y = percentage, fill = Phase)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("G1" = "#ff6b6b", "S" = "#4dabf7", "G2M" = "#51cf66")) +
  labs(title = "Cell Cycle Phase Distribution Across Genotypes and days",
       x = "Genotype . day", 
       y = "Percentage of Cells",
       fill = "Cell Cycle Phase") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("results/scRNAseq/5_cellcycle/stack_bar_combined.pdf", 
       plot = stack_bar_combined, width = 12, height = 6) 

