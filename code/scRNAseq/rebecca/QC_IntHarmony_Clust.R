
# module load singularity
#  singularity exec -B /run/user/114119227 /nfs/mm-isilon/bioinfcore/Common/singularity/single_cell_0.12.0.sif R
# https://nbisweden.github.io/excelerate-scRNAseq/session-qc/Quality_control.html

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(gridExtra)  
library(sctransform)
library(future)
library(scCATCH)

source("./scripts/SCfunctions.R")
options(future.globals.maxSize = 10000 * 1024^2) # avoid RAM limits in future package for SCT


# /nfs/mm-isilon/bioinfcore/ActiveProjects/Iwase_siwase_SC1_damki_11534-BZ/inputs/11534-BZ/10x_analysis_11534-BZ/Sample_11534-BZ-Pool01/per_sample_outs/
# 
# 11534-BZ-2/count/sample_filtered_feature_bc_matrix

Markers <- c("RAI1","NANOG","POU5F1","SOX2","GD3S","SMAD1","PAX6","FOXG1","OTX1","OTX2","EMX1","DIXDC1","NKX2-1","HOXB4","ISL1","FGF2","NOG","TBR1","BCL11B","SATB2","POU3F2","CUX1")
# not found in 20 day  NANOG, POU5F1, GD3S 

# ------
Dir <- "./inputs/11534-BZ/10x_analysis_11534-BZ/Sample_11534-BZ-Pool01/per_sample_outs/"
DirOutQC <- "./outputs/pdfs/QC_IntHarmony_Clust/"   # location for putting plots
SampID <- c("11534-BZ-2","11534-BZ-4","11534-BZ-6","11534-BZ-7","11534-BZ-10","11534-BZ-12","11534-BZ-14","11534-BZ-16","11534-BZ-18","11534-BZ-20","11534-BZ-22","11534-BZ-24")
SampName <- c("BC001","BC002","BC003","BC004","BC005","BC006","BC007","BC008","BC009","BC010","BC011","BC012")
Group <- c("Tmp1","Tmp1","TMP1",  "Tmp2","Tmp2","TMP2", "Tmp3","Tmp3","TMP3",  "Tmp4","Tmp4","TMP4")
Genotype <- c("WT","WT","KO","KO","HET","HET","WT","WT","KO","KO","HET","HET")
Day <- c("D20","D20","D20","D20","D20","D20","D60","D60","D60","D60","D60","D60")


# load Seurat Objects if restarting : (these are already filtered, so proceed to line 150
load(file="/nfs/mm-isilon/bioinfcore/ActiveProjects/Iwase_siwase_SC1_damki_11534-BZ/outputs/SeuObjs_harmony_Mt20_lambda/Seu_Day20.RData")
load(file="/nfs/mm-isilon/bioinfcore/ActiveProjects/Iwase_siwase_SC1_damki_11534-BZ/outputs/SeuObjs_harmony_Mt20_lambda/Seu_Day60.RData")

m20 <- data.filt.20@meta.data
m60 <- data.filt.60@meta.data

#An object of class Seurat 
#34536 features across 70472 samples within 2 assays 
#Active assay: SCT (16454 features, 3000 variable features)
# 3 layers present: counts, data, scale.data
# 1 other assay present: RNA
# 3 dimensional reductions calculated: pca, umap, tsne

table(m20$Group)
#HET.D20  KO.D20  WT.D20 
#  30615   10495   29362 
table(m20$orig.ident)
#BC001 BC002 BC003 BC004 BC005 BC006 
#17862 11500  5826  4669  8941 21674 

table(m60$Group)
#HET.D60  KO.D60  WT.D60 
#   6658    9222   11874 
table(m60$orig.ident)
#BC007 BC008 BC009 BC010 BC011 BC012 
# 5515  6359  4733  4489  2968  3690 

########

SeuListNoFilt <- list()
for(ii in 1:length(SampID)) {
     dat <- Read10X(data.dir = paste0(Dir, SampID[ii], "/count/sample_filtered_feature_bc_matrix/") )
     cat(SampID[ii], nrow(dat), ncol(dat), "\n")  
     mySeu <- CreateSeuratObject(counts = dat, project=SampName[ii])
     mySeu[["percent.mt"]] <- PercentageFeatureSet(mySeu, pattern = "^MT-")
#     mySeu[["percent.ribo"]] <- PercentageFeatureSet(mySeu, pattern = "^RP[SL]")
     mySeu[["Day"]] <- rep(Day[ii], ncol(mySeu))
     mySeu[["Genotype"]] <- rep(Genotype[ii], ncol(mySeu))
     mySeu[["Group"]] <- rep( paste(Genotype[ii], Day[ii], sep="."), ncol(mySeu))
#     head(mySeu@meta.data)
     aSamp <- SampName[[ii]]
     SeuListNoFilt[[aSamp]] <- mySeu
}

# 11534-BZ-2 18082 17935
# 11534-BZ-4 18082 11550
# 11534-BZ-6 18082 5860
# 11534-BZ-7 18082 4699
# 11534-BZ-10 18082 9038
# 11534-BZ-12 18082 21988
# 11534-BZ-14 18082 5542
# 11534-BZ-16 18082 6412
# 11534-BZ-18 18082 5066
# 11534-BZ-20 18082 5107
# 11534-BZ-22 18082 3132
# 11534-BZ-24 18082 3862


alldata <- merge( SeuListNoFilt[[SampName[[1]]]], 
c( SeuListNoFilt[[SampName[[2]]]], SeuListNoFilt[[SampName[[3]]]], SeuListNoFilt[[SampName[[4]]]], SeuListNoFilt[[SampName[[5]]]], SeuListNoFilt[[SampName[[6]]]],
SeuListNoFilt[[SampName[[7]]]], SeuListNoFilt[[SampName[[8]]]], SeuListNoFilt[[SampName[[9]]]], SeuListNoFilt[[SampName[[10]]]], 
SeuListNoFilt[[SampName[[11]]]], SeuListNoFilt[[SampName[[12]]]]   ),
add.cell.ids=SampName)

#day20Mer <-  merge( SeuListNoFilt[[SampName[[1]]]],
#c( SeuListNoFilt[[SampName[[2]]]], SeuListNoFilt[[SampName[[3]]]], SeuListNoFilt[[SampName[[4]]]], SeuListNoFilt[[SampName[[5]]]], SeuListNoFilt[[SampName[[6]]]] ), add.cell.ids=SampName[1:6])
#
#day60Mer <-  merge( SeuListNoFilt[[SampName[[7]]]],
#c( SeuListNoFilt[[SampName[[8]]]], SeuListNoFilt[[SampName[[9]]]], SeuListNoFilt[[SampName[[10]]]], SeuListNoFilt[[SampName[[11]]]], SeuListNoFilt[[SampName[[12]]]] ), add.cell.ids=SampName[7:12])
#


#> alldata
#An object of class Seurat 
#19059 features across 23955 samples within 1 assay 
#Active assay: RNA (19059 features, 0 variable features)
# 5 layers present: counts.MSM.1f, counts.MSM.2f, counts.C57.1f, counts.CBA.1f, counts.CBA.2f

#> table(Idents(alldata))
#MSM.1f MSM.2f C57.1f CBA.1f CBA.2f 
#  4576   3837   6936   7362   1244 

head(alldata@meta.data)
mm1 <- alldata@meta.data

plot_list <- list()
plot_list[["Vln_nFeat"]] <- VlnPlot(alldata, features = "nFeature_RNA", pt.size = 0.01, alpha=0, add=FALSE, raster=FALSE) + NoLegend()
plot_list[["Vln_nCount"]] <- VlnPlot(alldata, features = "nCount_RNA", pt.size = 0.01,alpha=0,  add=FALSE, raster=FALSE) + NoLegend()
plot_list[["Vln_percMito"]] <- VlnPlot(alldata, features = "percent.mt", pt.size =0.01, alpha=0,add=FALSE, raster=FALSE) + NoLegend()
doPlots(plot_list)

plot_list <- list()
plot_list[["Scatt_Count_Feat"]] <- FeatureScatter(alldata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot_list[["Scatt_Feat_Mito"]] <- FeatureScatter(alldata, feature1 = "nFeature_RNA", feature2 = "percent.mt")
doPlots(plot_list)

plot_list <- list()
plot_list[["Gp.Scatt_Count_Feat"]] <- FeatureScatter(alldata, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by="Group")
plot_list[["Gp.Scatt_Feat_Mito"]] <- FeatureScatter(alldata, feature1 = "nFeature_RNA", feature2 = "percent.mt", group.by="Group" )
doPlots(plot_list)


# based on summaries, only filter on Mt
# dim(alldata)   [1]  18082 100191

selMt <-   WhichCells(alldata, expression = percent.mt <= 20)

toKeep <- selMt
length(toKeep)  # just Mt: 98226 (5%), 100132 (20%) 
data.filt <- subset(alldata, cells = toKeep)
mm2 <- data.filt@meta.data
data.filt.20 <- subset(data.filt, cells=rownames(mm2)[which(mm2$Day=="D20")] )  #  18082 70472 (5%), 71037 (20%)
data.filt.60 <- subset(data.filt, cells=rownames(mm2)[which(mm2$Day=="D60")] )  # 18082 27754 (5%) , 29095 (20%)

plot_list <- list()
plot_list[["Vln_nFeat_20_mtFilt"]] <- VlnPlot(data.filt.20, features = "nFeature_RNA", pt.size = 0.01, alpha=0, add=FALSE, raster=FALSE) + NoLegend()
plot_list[["Vln_nCount_20_mtFilt"]] <- VlnPlot(data.filt.20, features = "nCount_RNA", pt.size = 0.01,alpha=0,  add=FALSE, raster=FALSE) + NoLegend()
plot_list[["Vln_percMito_20_mtFilt"]] <- VlnPlot(data.filt.20, features = "percent.mt", pt.size =0.01, alpha=0,add=FALSE, raster=FALSE) + NoLegend()
plot_list[["Gp.Scatt_Count_Feat_20_mtFilt"]] <- FeatureScatter(data.filt.20, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by="Group")
plot_list[["Gp.Scatt_Feat_Mito_20_mtFilt"]] <- FeatureScatter(data.filt.20 , feature1 = "nFeature_RNA", feature2 = "percent.mt", group.by="Group" )
doPlots(plot_list)

plot_list <- list()
plot_list[["Vln_nFeat_60_mtFilt"]] <- VlnPlot(data.filt.60, features = "nFeature_RNA", pt.size = 0.01, alpha=0, add=FALSE, raster=FALSE) + NoLegend()
plot_list[["Vln_nCount_60_mtFilt"]] <- VlnPlot(data.filt.60, features = "nCount_RNA", pt.size = 0.01,alpha=0,  add=FALSE, raster=FALSE) + NoLegend()
plot_list[["Vln_percMito_60_mtFilt"]] <- VlnPlot(data.filt.60, features = "percent.mt", pt.size =0.01, alpha=0,add=FALSE, raster=FALSE) + NoLegend()
plot_list[["Gp.Scatt_Count_Feat_60_mtFilt"]] <- FeatureScatter(data.filt.60, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by="Group")
plot_list[["Gp.Scatt_Feat_Mito_60_mtFilt"]] <- FeatureScatter(data.filt.60 , feature1 = "nFeature_RNA", feature2 = "percent.mt", group.by="Group" )
doPlots(plot_list)



# pre-integration day 20 and day 60 separately norm SCT
Seu <- data.filt.60
m0 <- Seu@meta.data
Seu <- SCTransform(Seu)
Seu <- RunPCA(Seu, npcs = 30)
myPCs <- check_pcs(Seu)  # 19 for both
Seu <- RunUMAP(Seu, dims = 1:myPCs)
Seu <- RunTSNE(Seu, dims = 1:myPCs)
#  load("./outputs/SeuObjs_harmony_Mt20_lambda/Seu_Day60.RData")  # Seu<- data.filt.60  #

plot_list <- list()
  plot_list[["UMAP_D60"]]   <- DimPlot(Seu, reduction = "umap", group.by = c("Group"))
  plot_list[["TSNE_D60"]]   <- DimPlot(Seu, reduction = "tsne", group.by = c("Group"))
  plot_list[["PCA_D60"]]   <- DimPlot(Seu, reduction = "pca", group.by = c("Group"))
doPlots(plot_list)
data.filt.60 <- Seu
#save(data.filt.60, file=paste0(DirOutQC, "/Seu_Day60.RData"))


Seu <- data.filt.20
m0 <- Seu@meta.data
Seu <- SCTransform(Seu)
Seu <- RunPCA(Seu, npcs = 30)
myPCs <- check_pcs(Seu)  # 19 for both
Seu <- RunUMAP(Seu, dims = 1:myPCs)
Seu <- RunTSNE(Seu, dims = 1:myPCs)
#  load("./outputs/SeuObjs_harmony_Mt20_lambda/Seu_Day20.RData")  # data.filt.20  # 

plot_list <- list()
  plot_list[["UMAP_D20"]]   <- DimPlot(Seu, reduction = "umap", group.by = c("Group"))
  plot_list[["TSNE_D20"]]   <- DimPlot(Seu, reduction = "tsne", group.by = c("Group"))
  plot_list[["PCA_D20"]]   <- DimPlot(Seu, reduction = "pca", group.by = c("Group"))
doPlots(plot_list)
data.filt.20 <- Seu
#save(data.filt.20, file=paste0(DirOutQC, "/Seu_Day20.RData"))

# --------
# HERE ---
Seu <- data.filt.60
m0 <- Seu@meta.data
seuWT <- subset(Seu, cells=rownames(m0)[which(m0$Genotype=="WT")] )  
seuKO <- subset(Seu, cells=rownames(m0)[which(m0$Genotype=="KO")] )  
seuHET <- subset(Seu, cells=rownames(m0)[which(m0$Genotype=="HET")] )  
plot_list <- list()
  plot_list[["UMAP_WT60_byIdent"]]   <- DimPlot(seuWT, reduction = "umap", group.by = c("orig.ident"))
  plot_list[["UMAP_KO60_byIdent"]]   <- DimPlot(seuKO, reduction = "umap", group.by = c("orig.ident"))
  plot_list[["UMAP_HET60_byIdent"]]   <- DimPlot(seuHET, reduction = "umap", group.by = c("orig.ident"))
doPlots(plot_list)

Seu <- data.filt.20
m0 <- Seu@meta.data
seuWT <- subset(Seu, cells=rownames(m0)[which(m0$Genotype=="WT")] )  
seuKO <- subset(Seu, cells=rownames(m0)[which(m0$Genotype=="KO")] )
seuHET <- subset(Seu, cells=rownames(m0)[which(m0$Genotype=="HET")] )
plot_list <- list()
  plot_list[["UMAP_WT20_byIdent"]]   <- DimPlot(seuWT, reduction = "umap", group.by = c("orig.ident"))
  plot_list[["UMAP_KO20_byIdent"]]   <- DimPlot(seuKO, reduction = "umap", group.by = c("orig.ident"))
  plot_list[["UMAP_HET20_byIdent"]]   <- DimPlot(seuHET, reduction = "umap", group.by = c("orig.ident"))
doPlots(plot_list)






# -----

integrated20 <- IntegrateLayers(object = data.filt.20 , method = HarmonyIntegration, npcs=myPCs, lambda=0.5,
                       orig.reduction = "pca", new.reduction = 'harmony', assay = "SCT")
integrated20 <- FindNeighbors(integrated20, reduction = "harmony", dims = 1:myPCs)
integrated20 <- FindClusters(integrated20, resolution = 0.1, cluster.name = "harmony_clusters")  # 5 clusters

integrated60 <- IntegrateLayers(object = data.filt.60 , method = HarmonyIntegration, npcs=myPCs, lambda=0.5,
                       orig.reduction = "pca", new.reduction = 'harmony', assay = "SCT")
integrated60 <- FindNeighbors(integrated60, reduction = "harmony", dims = 1:myPCs)
integrated60 <- FindClusters(integrated60, resolution = 0.1, cluster.name = "harmony_clusters")  # 8 clusters



# -----

load("./inputs/scCatchDat/cellmatch.rda")
iKeep <-  grep("brain", cellmatch$tissue, ignore.case = T, value = F)
cellmatchBrain <- cellmatch[iKeep, ]
brainTiss <- unique(grep("brain", cellmatch$tissue, ignore.case = T, value = T))
# There are 700 potential marker genes in CellMatch database for Human on Brain, Fetal brain, Midbrain, Embryonic brain.

nnamm <- "integrated20"
adb <- "_BrainTis"
obj <- integrated20
mm20 <- obj@meta.data
normMat <- LayerData(obj, assay = "SCT", layer = "data")
normMat <- rev_gene(data = normMat, data_type = "data", species = "Human", geneinfo = geneinfo)
clust <- mm20$seurat_clusters
SCobj <- createscCATCH(data = normMat, cluster = as.vector(clust))
SCobjHu <- findmarkergene(object = SCobj, species = "Human", marker = cellmatch, tissue = brainTiss)
SCobjHu <- findcelltype(object = SCobjHu)
aDF <- cbind(Sample=rep(nnamm, nrow(SCobjHu@celltype)), SCobjHu@celltype)
write.csv(aDF, file=paste0(DirOutQC,"ScCatch_", nnamm,adb, ".csv"))

nnamm <- "integrated60"
adb <- "_BrainTis"
obj <- integrated60
mm60 <- obj@meta.data
normMat <- LayerData(obj, assay = "SCT", layer = "data")
normMat <- rev_gene(data = normMat, data_type = "data", species = "Human", geneinfo = geneinfo)
clust <- mm60$seurat_clusters
SCobj <- createscCATCH(data = normMat, cluster = as.vector(clust))
SCobjHu <- findmarkergene(object = SCobj, species = "Human", marker = cellmatch, tissue = brainTiss)
SCobjHu <- findcelltype(object = SCobjHu)
aDF <- cbind(Sample=rep(nnamm, nrow(SCobjHu@celltype)), SCobjHu@celltype)
write.csv(aDF, file=paste0(DirOutQC,"ScCatch_", nnamm,adb, ".csv"))


# -----


integrated20 <- RunPCA(integrated20, reduction = "harmony", dims = 1:myPCs, reduction.name = "pca.harmony")
integrated20 <- RunUMAP(integrated20, reduction = "harmony", dims = 1:myPCs, reduction.name = "umap.harmony")
integrated20 <- RunTSNE(integrated20, reduction = "harmony", dims = 1:myPCs, reduction.name = "tsne.harmony")
plot_list <- list()
plot_list[["umap.20.harmony_Group"]] <- DimPlot(integrated20, reduction = "umap.harmony", group.by = "Group",shuffle = TRUE) #+ scale_color_viridis_d()
plot_list[["tsne.20.harmony_Group"]] <- DimPlot(integrated20, reduction = "tsne.harmony", group.by = "Group",shuffle = TRUE) #+ scale_color_viridis_d()
plot_list[["pca.20.harmony_Group"]] <- DimPlot(integrated20, reduction = "pca.harmony", group.by = "Group",shuffle = TRUE) #+ scale_color_viridis_d()
plot_list[["umap.20.harmony_Clust01"]] <- DimPlot(integrated20, reduction = "umap.harmony", group.by = "seurat_clusters",shuffle = TRUE) #+ scale_color_viridis_d()
plot_list[["tsne.20.harmony_Clust01"]] <- DimPlot(integrated20, reduction = "tsne.harmony", group.by = "seurat_clusters",shuffle = TRUE) #+ scale_color_viridis_d()
plot_list[["pca.20.harmony_Clust01"]] <- DimPlot(integrated20, reduction = "pca.harmony", group.by = "seurat_clusters",shuffle = TRUE) #+ scale_color_viridis_d()
doPlots(plot_list)



integrated60 <- RunPCA(integrated60, reduction = "harmony", dims = 1:myPCs, reduction.name = "pca.harmony")
integrated60 <- RunUMAP(integrated60, reduction = "harmony", dims = 1:myPCs, reduction.name = "umap.harmony")
integrated60 <- RunTSNE(integrated60, reduction = "harmony", dims = 1:myPCs, reduction.name = "tsne.harmony")
plot_list <- list()
plot_list[["umap.60.harmony_Group"]] <- DimPlot(integrated60, reduction = "umap.harmony", group.by = "Group",shuffle = TRUE) #+ scale_color_viridis_d()
plot_list[["tsne.60.harmony_Group"]] <- DimPlot(integrated60, reduction = "tsne.harmony", group.by = "Group",shuffle = TRUE) #+ scale_color_viridis_d()
plot_list[["pca.60.harmony_Group"]] <- DimPlot(integrated60, reduction = "pca.harmony", group.by = "Group",shuffle = TRUE) #+ scale_color_viridis_d()
plot_list[["umap.60.harmony_Clust01"]] <- DimPlot(integrated60, reduction = "umap.harmony", group.by = "seurat_clusters",shuffle = TRUE) #+ scale_color_viridis_d()
plot_list[["tsne.60.harmony_Clust01"]] <- DimPlot(integrated60, reduction = "tsne.harmony", group.by = "seurat_clusters",shuffle = TRUE) #+ scale_color_viridis_d()
plot_list[["pca.60.harmony_Clust01"]] <- DimPlot(integrated60, reduction = "pca.harmony", group.by = "seurat_clusters",shuffle = TRUE) #+ scale_color_viridis_d()
doPlots(plot_list)

# ----------

col <- c("blue", "red", "green","orange","orchid", "lightblue","brown", "grey", "darkgreen", "violet", "gold3","palevioletred3",  "darkblue" )
p <- DotPlot(integrated20, features = Markers, cols=col, dot.scale = 8, split.by = "seurat_clusters") + RotatedAxis()
tiff(file=paste0(DirOutQC, "dotPlot_exprMarkers_20.tiff"), width = 16, height = 8, units = 'in', res = 100)
print(p); dev.off()

p <- DotPlot(integrated60, features = Markers, cols=col, dot.scale = 8, split.by = "seurat_clusters") + RotatedAxis()
tiff(file=paste0(DirOutQC, "dotPlot_exprMarkers_60.tiff"), width = 16, height = 8, units = 'in', res = 100)
print(p); dev.off()

#   pre-integration
SCTexpr20 <- integrated20[["SCT"]]$counts  
genesum20 <- apply(SCTexpr20, 1, sum)
toRM <- which(genesum20 < 5) ; if(length(toRM) > 0) {genesum20 <- genesum20[-toRM] }
Markers20 <- intersect(names(genesum20), Markers)
key <- "_UMAP_20day"
plot_list <- list()
for(ii in 1:length(Markers20)) {
   aMar <- Markers20[ii]
   nnam <- paste0(aMar,key)
   plot_list[[nnam]] <- FeaturePlot(integrated20, features=aMar , max.cutoff = 3, cols = c("grey",   "red"), reduction = "umap")
}
doPlots(plot_list)

SCTexpr60 <- integrated60[["SCT"]]$counts  
genesum60 <- apply(SCTexpr60, 1, sum)
toRM <- which(genesum60 < 5) ; if(length(toRM) > 0) {genesum60 <- genesum60[-toRM] }
Markers60 <- intersect(names(genesum60), Markers)
key <- "_UMAP_60day"
plot_list <- list()
for(ii in 1:length(Markers60)) {
   aMar <- Markers60[ii]
   nnam <- paste0(aMar,key)
   plot_list[[nnam]] <- FeaturePlot(integrated60, features=aMar , max.cutoff = 3, cols = c("grey",   "red"), reduction = "umap")
}
doPlots(plot_list)

#   post-integration
SCTexpr20 <- integrated20[["SCT"]]$counts  
genesum20 <- apply(SCTexpr20, 1, sum)
toRM <- which(genesum20 < 5) ; if(length(toRM) > 0) {genesum20 <- genesum20[-toRM] }
Markers20 <- intersect(names(genesum20), Markers)
key <- "_UMAP_harmony_20day"
plot_list <- list()
for(ii in 1:length(Markers20)) {
   aMar <- Markers20[ii]
   nnam <- paste0(aMar,key)
   plot_list[[nnam]] <- FeaturePlot(integrated20, features=aMar , max.cutoff = 3, cols = c("grey",   "red"), reduction = "umap.harmony")
}
doPlots(plot_list)

SCTexpr60 <- integrated60[["SCT"]]$counts  
genesum60 <- apply(SCTexpr60, 1, sum)
toRM <- which(genesum60 < 5) ; if(length(toRM) > 0) {genesum60 <- genesum60[-toRM] }
Markers60 <- intersect(names(genesum60), Markers)
key <- "_UMAP_harmony_60day"
plot_list <- list()
for(ii in 1:length(Markers60)) {
   aMar <- Markers60[ii]
   nnam <- paste0(aMar,key)
   plot_list[[nnam]] <- FeaturePlot(integrated60, features=aMar , max.cutoff = 3, cols = c("grey",   "red"), reduction = "umap.harmony")
}
doPlots(plot_list)


####### 


save(integrated20, file=paste0(DirOutQC, "integrated20.RData")) ; save(integrated60, file=paste0(DirOutQC, "integrated60.RData"))

load(file=paste0("./outputs/SeuObjs_harmony_Mt20_lambda/", "integrated20.RData")) 
load(file=paste0("./outputs/SeuObjs_harmony_Mt20_lambda/"  , "integrated60.RData"))

####### 

integrated20 <- PrepSCTFindMarkers(integrated20)
markers20 <- FindAllMarkers(integrated20, only.pos = FALSE)
write.csv(markers20, file=paste0(DirOutQC,"FindAllMarkers_ALL_20.csv"))

integrated60 <- PrepSCTFindMarkers(integrated60)
markers60 <- FindAllMarkers(integrated60, only.pos = FALSE)
write.csv(markers60, file=paste0(DirOutQC,"FindAllMarkers_ALL_60.csv"))

#######

mm20 <- integrated20@meta.data
gpClust <- paste0(mm20$Group, " cluster.",    mm20$seurat_clusters )
nn20 <- as.data.frame(table(gpClust))
write.csv(nn20, file=paste0(DirOutQC, "cluster_Gp_day20.csv"))

mm60 <- integrated60@meta.data
gpClust <- paste0(mm60$Group, " cluster.",    mm60$seurat_clusters )
nn60 <- as.data.frame(table(gpClust))
write.csv(nn60, file=paste0(DirOutQC, "cluster_Gp_day60.csv"))


#########




####### 



# The RC (Relative Counts) normalization method does not correct for sequencing depth. Instead, it transforms the raw counts into relative abundances, which represent the proportion of each gene's expression relative to the total expression in that cell.
# RC is CPM








