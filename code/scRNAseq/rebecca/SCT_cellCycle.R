
# Take top ORA results from each cluster / subcluster.
# For given cluster, intersect seurat markers with the genes from those top ORA results for dot plots 
# 1- Seurat Markers in HypoMarkers
# 2- 50% of the cells in the cluster express the gene
# 3- remove duplicates ( gene is significant for only one cluster)
# 4- top ranked positive fold change


# module load singularity
#  singularity exec -B /run/user/114119227 /nfs/mm-isilon/bioinfcore/Common/singularity/single_cell_0.12.0.sif R
# https://nbisweden.github.io/excelerate-scRNAseq/session-qc/Quality_control.html
# https://hbctraining.github.io/scRNA-seq_online/lessons/seurat_subclustering.html

source("./scripts/SCfunctions.R")

# ------

DirOutQC <- "./outputs/pdfs/SCT_cellCycle/"

# human CC genes
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

# for mouse
#s.genes <- c("Mcm5", "Pcna", "Tyms", "Fen1", "Mcm2", "Mcm4", 
#           "Rrm1", "Ung", "Gins2", "Mcm6", "Cdca7", "Dtl", 
#           "Prim1", "Uhrf1", "Mlf1ip", "Hells", "Rfc2", "Rpa2", 
#           "Nasp", "Rad51ap1", "Gmnn", "Wdr76", "Slbp", "Ccne2", 
#           "Ubr7", "Pold3", "Msh2", "Atad2", "Rad51", "Rrm2", 
#           "Cdc45", "Cdc6", "Exo1", "Tipin", "Dscc1", "Blm", 
#           "Casp8ap2", "Usp1", "Clspn", "Pola1", "Chaf1b", "Brip1", 
#           "E2f8")
#g2m.genes <- c("Hmgb2", "Cdk1", "Nusap1", "Ube2c", "Birc5", "Tpx2", "Top2a", 
#           "Ndc80", "Cks2", "Nuf2", "Cks1b", "Mki67", "Tmpo", "Cenpf", 
#           "Tacc3", "Fam64a", "Smc4", "Ccbn2", "Ckap2l", "Ckap2", "Aurkb", 
#           "Bub1", "Kif11", "Anp32e", "Tubb4b", "Gtse1", "Kif20b", "Hjurp", 
#           "Cdca3", "Hn1", "Cdc20", "Ttk", "Cdc25c", "Kif2c", "Rangap1", 
#           "Ncapd2", "Dlgap5", "Cdca2", "Cdca8", "Ect2", "Kif23", "Hmmer", 
#           "Aurka", "Pscr1", "Anln", "Lbr", "Ckap5", "Cenpe", "Ctcf", 
#	"Nek2", "G2e3", "Gas2l3", "Cbx5", "Cenpa")

# 


load(file="./outputs/SeuObjs_harmony_Mt20_lambda/integrated20.RData")
load(file="./outputs/SeuObjs_harmony_Mt20_lambda/integrated60.RData")

m20 <- integrated20@meta.data
m60 <- integrated60@meta.data



#  <<<<<<<<<<<<

#DirOutQC <- "./outputs/QCplots_harmony_Mt20_lambda/CellCycleViolins/"
DirOutQCBAK <- DirOutQC



DirOutQC <- paste0(DirOutQCBAK,"CC_D20/")
sctSeuFilt <- integrated20

sctSeuFilt <- CellCycleScoring(sctSeuFilt, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


cell_cycle_scores <- sctSeuFilt@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(
    avg_S_Score = mean(S.Score),
    avg_G2M_Score = mean(G2M.Score)
  )
plotList <- list()
plotList[["Main_Violin_s.score"]] <- VlnPlot(sctSeuFilt, features = "S.Score", group.by = "seurat_clusters")
plotList[["Main_Violin_g2m.score"]] <- VlnPlot(sctSeuFilt, features = "G2M.Score", group.by = "seurat_clusters")
plotList[["Main_Dim_CCphase"]] <- DimPlot(sctSeuFilt, group.by = "Phase")
doPlots(plotList)
splitSeu <- SplitObject(sctSeuFilt, split.by = "orig.ident")
plotList <- list()
for(ii in 1:length(splitSeu)) {
    aNam <- paste0("Main_Dim__",  names(splitSeu)[ii])
    plotList[[aNam]]  <- DimPlot(splitSeu[[names(splitSeu)[ii]]], group.by = "Phase")
    aNam <- paste0("Main_Violin_s.score_",  names(splitSeu)[ii])
    plotList[[aNam]]  <- VlnPlot(splitSeu[[names(splitSeu)[ii]]], features = "S.Score",    group.by = "seurat_clusters")
    aNam <- paste0("Main_Violin_g2m.score_",  names(splitSeu)[ii])
    plotList[[aNam]]  <- VlnPlot(splitSeu[[names(splitSeu)[ii]]],  features = "G2M.Score",  group.by = "seurat_clusters")

}
doPlots(plotList)













