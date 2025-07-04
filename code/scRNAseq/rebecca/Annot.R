
# module load singularity
#  singularity exec -B /run/user/114119227 /nfs/mm-isilon/bioinfcore/Common/singularity/single_cell_0.12.0.sif R
# https://nbisweden.github.io/excelerate-scRNAseq/session-qc/Quality_control.html


source("./scripts/SCfunctions.R")
options(future.globals.maxSize = 10000 * 1024^2) # avoid RAM limits in future package for SCT

DirOutQC <- "./outputs/pdfs/Annot/"
DirSeuMarkers <- "./outputs/harmony_Mt20_lambda/ClusterAnnot/"

# These are filtered, and Harmony integrated: "integrated20" "integrated60" 
load(file="/nfs/mm-isilon/bioinfcore/ActiveProjects/Iwase_siwase_SC1_damki_11534-BZ/outputs/SeuObjs_harmony_Mt20_lambda/integrated20.RData")
load(file="/nfs/mm-isilon/bioinfcore/ActiveProjects/Iwase_siwase_SC1_damki_11534-BZ/outputs/SeuObjs_harmony_Mt20_lambda/integrated60.RData")

m20 <- integrated20@meta.data
m60 <- integrated60@meta.data

cellType <- rep("Neurons", nrow(m60))
cellType[which(m60$seurat_clusters == "1")] <- "OPC"
cellType[which(m60$seurat_clusters == "2")] <- "HRGL"
cellType[which(m60$seurat_clusters == "3")] <- "NPC"
cellType[which(m60$seurat_clusters == "4")] <- "Astrocytes"
cellType[which(m60$seurat_clusters == "5")] <- "NPC-CC"
cellType[which(m60$seurat_clusters == "6")] <- "Apop"
cellType[which(m60$seurat_clusters == "7")] <- "Immune"
integrated60@meta.data <- cbind(integrated60@meta.data, cellType)

cellType <- rep("NPC-early", nrow(m20))
cellType[which(m20$seurat_clusters == "1")] <- "OPC"
cellType[which(m20$seurat_clusters == "2")] <- "NPC-G2M"
cellType[which(m20$seurat_clusters == "3")] <- "NPC-S"
cellType[which(m20$seurat_clusters == "4")] <- "Neurons"
integrated20@meta.data <- cbind(integrated20@meta.data, cellType)


SeuWT20 <- subset(integrated20, subset = Group == "WT.D20")
SeuKO20 <- subset(integrated20, subset = Group == "KO.D20")
SeuHet20 <- subset(integrated20, subset = Group == "HET.D20")
SeuWT60 <- subset(integrated60, subset = Group == "WT.D60")
SeuKO60 <- subset(integrated60, subset = Group == "KO.D60")
SeuHet60 <- subset(integrated60, subset = Group == "HET.D60")

plotList <- list()
plotList[["UMAP_20annot"]] <- DimPlot(integrated20,  reduction = "umap.harmony", group.by = c("cellType"), label=T)
plotList[["UMAP_60annot"]] <- DimPlot(integrated60,  reduction = "umap.harmony", group.by = c("cellType"), label=T)
plotList[["UMAP_WT20"]]   <- DimPlot(SeuWT20, reduction = "umap.harmony", group.by = c("cellType"), label=T)
plotList[["UMAP_KO20"]]   <- DimPlot(SeuKO20, reduction = "umap.harmony", group.by = c("cellType"), label=T)
plotList[["UMAP_HET20"]]   <- DimPlot(SeuHet20, reduction = "umap.harmony", group.by = c("cellType"), label=T)
plotList[["UMAP_WT60"]]   <- DimPlot(SeuWT60, reduction = "umap.harmony", group.by = c("cellType"), label=T)
plotList[["UMAP_KO60"]]   <- DimPlot(SeuKO60, reduction = "umap.harmony", group.by = c("cellType"), label=T)
plotList[["UMAP_HET60"]]   <- DimPlot(SeuHet60, reduction = "umap.harmony", group.by = c("cellType"), label=T)
doPlots(plotList)

save(integrated20, file=paste0(DirOutQC, "integrated20.RData"))
save(integrated60, file=paste0(DirOutQC, "integrated60.RData"))


table(SeuWT20@meta.data$harmony_clusters)
table(SeuHet20@meta.data$harmony_clusters)
table(SeuKO20@meta.data$harmony_clusters)
#WT  0     1     2     3     4 
# 10507 11248  3932  2831   962 
#Het 0     1     2     3     4 
# 12301  6788  6995  4329   590 
#KO  0    1    2    3    4 
# 2592 5330 1813  804   15 

table(SeuWT60@meta.data$harmony_clusters)
table(SeuHet60@meta.data$harmony_clusters)
table(SeuKO20@meta.data$harmony_clusters)
#WT 0    1    2    3    4    5    6    7 
# 6302 2393 1227  980  707  308   21   16 
#Het 0    1    2    3    4    5    6    7 
# 4294 1272  813  276  180  108   36    3 
#KO  0    1    2    3    4    5    6    7 
# 5733 2685  612  621  141  257   80   30 




seuMarkers20 <- read.csv(file=paste0(DirSeuMarkers, "FindAllMarkers_ALL_20.csv"), row.names=1)
seuMarkers60 <- read.csv(file=paste0(DirSeuMarkers, "FindAllMarkers_ALL_60.csv"),row.names=1)
seuMarkPos20 <- seuMarkers20[which(seuMarkers20$avg_log2FC > 0), ]
seuMarkPos60 <- seuMarkers60[which(seuMarkers60$avg_log2FC > 0), ]

SCTexpr <- integrated20[["SCT"]]$counts
genes <- rownames(SCTexpr)
background20 <- data.frame(genes, row.names=genes)

SCTexpr <- integrated20[["SCT"]]$counts
background20 <- rownames(SCTexpr)
SCTexpr <- integrated60[["SCT"]]$counts
background60 <- rownames(SCTexpr)

seuMarkPos20sig <- MarkerFilt1(seuMarkers20)
seuMarkPos60sig <- MarkerFilt1(seuMarkers60)

###########
SCTexpr20 <- integrated20[["SCT"]]$counts
allExpGenes20 <- rownames(SCTexpr20)
SCTexpr60 <- integrated60[["SCT"]]$counts
allExpGenes60 <- rownames(SCTexpr60)

###########
ShiMarkers <- ShiLiveseyMarkers()
# ZhongMajorList, ZhongPosList, uniqZhongPosList, ZhongPosSubClustList, uniqZhongPosSubClustList
load("./inputs/ZhongMarkers/ZhongLists.RData")
MarkerList <- ZhongMajorList
key <- "ZhongMajor"

qGenes <- unlist(MarkerList)
plotList <- list()
for(ii in 1:length(qGenes)) {
    aG <- qGenes[ii]
    sets_with_gene <- names(MarkerList)[sapply(MarkerList, function(x) aG %in% x)]
    tit <- paste0("FeatPlot20_",key,"_", paste(sets_with_gene, collapse=".") , "__",aG)
    if(aG %in% allExpGenes20) {
         p <- FeaturePlot(integrated20, features = aG, pt.size = 3, cols = c("lightgray", "blue"), reduction="umap.harmony")
         p <- p + theme(legend.position = "right") + theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),   axis.text.y = element_text(size = 12))
         p <- p + scale_color_gradientn(colors = c("lightgrey", "lightblue", "blue", "darkblue"))
         plotList[[tit]] <- p
    }
}
doPlots(plotList)

plotList <- list()
for(ii in 1:length(qGenes)) {
    aG <- qGenes[ii]
    sets_with_gene <- names(MarkerList)[sapply(MarkerList, function(x) aG %in% x)] 
    tit <- paste0("FeatPlot60_", key,"_",   paste(sets_with_gene, collapse=".") , "__",aG)
    if(aG %in% allExpGenes60) {
         p <- FeaturePlot(integrated60, features = aG, pt.size = 3, cols = c("lightgray", "blue"), reduction="umap.harmony")
         p <- p + theme(legend.position = "right") + theme(axis.text.x = element_text(size = 12, angle = 45, hjust = 1),   axis.text.y = element_text(size = 12))
         p <- p + scale_color_gradientn(colors = c("lightgrey", "lightblue", "blue", "darkblue"))
         plotList[[tit]] <- p
    }   
}
doPlots(plotList)


############
# move all of these feature plots to ./outputs/FeaturePlots/ZhongMajor/ FeatPlot20 and FeatPlot60
#There are 1000's so cross them with Seurat FindAllMarkers  seuMarkPos20  seuMarkPos60
# length(unlist(ZhongMajorList)) = 1023
seuMarkPos20Zhong <- MarkerFilt1(seuMarkPos20[which(seuMarkPos20$gene %in% unlist(ZhongMajorList)), ])   #  1266 before sigFilt 675 after sigFilt (MarkerFilt1)
seuMarkPos60Zhong <- MarkerFilt1(seuMarkPos60[which(seuMarkPos60$gene %in% unlist(ZhongMajorList)), ] )  #  2318 before sigFilt 1720  after sigFilt (MarkerFilt1)

duped <- seuMarkPos20Zhong$gene[which(duplicated(seuMarkPos20Zhong$gene))]  # 144
seuMarkPos20ZhongNodup <- seuMarkPos20Zhong[-which( seuMarkPos20Zhong$gene %in% duped), ] # 396 Could mean there are too many clusters

duped <- seuMarkPos60Zhong$gene[which(duplicated(seuMarkPos60Zhong$gene))]  # 1016
seuMarkPos60ZhongNodup <- seuMarkPos60Zhong[-which( seuMarkPos60Zhong$gene %in% duped), ] # 206  Could mean there are too many clusters

featFil20 <- list.files(path="./outputs/FeaturePlots/ZhongMajor/FeatPlot20")
V1 <- paste0("__",  seuMarkPos20ZhongNodup$gene , ".png")
Fil2Move <- featFil20[sapply(featFil20, function(v2_item) any(sapply(V1, grepl, x = v2_item)))]  # 396
Fil2Move <- paste0("./outputs/FeaturePlots/ZhongMajor/FeatPlot20/",Fil2Move)
new_paths <- paste0("./outputs/FeaturePlots/ZhongMajor/FeatPlot20/SeuMarkerClusters/", basename(Fil2Move))
file.rename(Fil2Move, new_paths)

featFil60 <- list.files(path="./outputs/FeaturePlots/ZhongMajor/FeatPlot60")
V1 <- paste0("__",  seuMarkPos60ZhongNodup$gene , ".png")
Fil2Move <- featFil60[sapply(featFil60, function(v2_item) any(sapply(V1, grepl, x = v2_item)))]  # 396
Fil2Move <- paste0("./outputs/FeaturePlots/ZhongMajor/FeatPlot60/",Fil2Move)
new_paths <- paste0("./outputs/FeaturePlots/ZhongMajor/FeatPlot60/SeuMarkerClusters/", basename(Fil2Move))
file.rename(Fil2Move, new_paths)


seuMarkPos20ZhongNodup$CellType <- sapply(seuMarkPos20ZhongNodup$gene, get_clusters_for_gene, ZhongMajorList = ZhongMajorList)
seuMarkPos60ZhongNodup$CellType <- sapply(seuMarkPos60ZhongNodup$gene, get_clusters_for_gene, ZhongMajorList = ZhongMajorList)

seuMarkPos20Zhong$CellType <- sapply(seuMarkPos20Zhong$gene, get_clusters_for_gene, ZhongMajorList = ZhongMajorList)
seuMarkPos60Zhong$CellType <- sapply(seuMarkPos60Zhong$gene, get_clusters_for_gene, ZhongMajorList = ZhongMajorList)


aDir <- "./outputs/QCplots_harmony_Mt20_lambda/ClusterAnnot/SeuratxZhong/"
write.csv(seuMarkPos20ZhongNodup, file=paste0(aDir, "seuMarkPos20ZhongNodup.csv"))
write.csv(seuMarkPos60ZhongNodup, file=paste0(aDir, "seuMarkPos60ZhongNodup.csv"))
write.csv(seuMarkPos20Zhong, file=paste0(aDir, "seuMarkPos20Zhong.csv"))
write.csv(seuMarkPos60Zhong, file=paste0(aDir, "seuMarkPos60Zhong.csv"))

#######

#genes in clusters tables
 sapply(SeuMarkerRes2List(seuMarkPos20), function(x) length(x))
#   0    1    2    3    4 
# 3787 4693 5509 6244 4331 


sapply(SeuMarkerRes2List(seuMarkPos60), function(x) length(x))
#   0    1    2    3    4    5    6    7 
#4110 5806 2223 4412 3943 3876 2299 2410 

#######  From Lin
# Over-representation
oraDir <- "./outputs/QCplots_harmony_Mt20_lambda/ClusterAnnot/SeuratxZhong/ORA/"

# seuMarkPos20
# background20
# seuMarkPos20Zhong

# ZhongMajorList
oraRes <- doRTTora( SeuMarkerRes2List(seuMarkPos20sig) , ZhongMajorList,  background20)      
write.csv(oraRes, file=paste0(oraDir, "ORA20_ZhongMajorSig.csv"))
oraRes <- doRTTora( SeuMarkerRes2List(seuMarkPos60sig) , ZhongMajorList,  background60)     
write.csv(oraRes, file=paste0(oraDir, "ORA60_ZhongMajorSig.csv"))

# ZhongPosList
oraRes <- doRTTora( SeuMarkerRes2List(seuMarkPos20sig) , ZhongPosList,  background20)      
write.csv(oraRes, file=paste0(oraDir, "ORA20_ZhongPosSig.csv"))
oraRes <- doRTTora( SeuMarkerRes2List(seuMarkPos60sig) , ZhongPosList,  background60)     
write.csv(oraRes, file=paste0(oraDir, "ORA60_ZhongPosSig.csv"))

# ShiMarkers - All crap results
oraRes <- doRTTora( SeuMarkerRes2List(seuMarkPos20sig) , ShiMarkers,  background20)
write.csv(oraRes, file=paste0(oraDir, "ORA20_ShiMarkers.csv"))
oraRes <- doRTTora( SeuMarkerRes2List(seuMarkPos60sig) , ShiMarkers,  background60)
write.csv(oraRes, file=paste0(oraDir, "ORA60_ShiMarkers.csv"))

# ZhongPosSubClustList
oraRes <- doRTTora( SeuMarkerRes2List(seuMarkPos20sig) , ZhongPosSubClustList,  background20)
write.csv(oraRes, file=paste0(oraDir, "ORA20_ZhongPosSubClustList.csv"))
oraRes <- doRTTora( SeuMarkerRes2List(seuMarkPos60sig) , ZhongPosSubClustList,  background60)
write.csv(oraRes, file=paste0(oraDir, "ORA60_ZhongPosSubClustList.csv"))


# Markers are mixed between clusters. Is the clustering not correct? 
# Try to cluster WT.60 alone without integration. do the markers fall out?






