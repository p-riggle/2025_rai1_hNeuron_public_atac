
library("pheatmap")
source("./scripts/SCfunctions.R")



Fin <- "./outputs/harmony_Mt20_lambda/"
DirOutQC <- "./outputs/pdfs/MarkerHeatmap/"

myMarkers <- list()
myMarkers[["early.stage.NPC"]] <- c("SOX2", "PAX6", "SIX3", "OTX2", "EMX2", "TBR1", "EOMES", "NEUROD1", "NEUROD2")
myMarkers[["early.cortical.neurons"]] <- c("FOXP2", "TBR1", "NEUROD1", "NEUROD2", "SLC17A6", "SLC32A1")
myMarkers[["interneuronsVIP"]] <- c("VIP", "SST", "NOS1", "GAD1", "GAD2", "SLC32A1", "DLX1", "DLX2")
myMarkers[["intermediate.progenitors"]] <- c("TBR2", "EOMES", "NEUROD1", "NEUROD2")
myMarkers[["developmental.cortical"]] <- c("EMX2", "OTX2", "PAX6", "SOX2", "TBR1", "SIX3", "FEZF2", "DLX1", "DLX2")
myMarkers[["glutamatergic.ExcitNeu"]] <- c("SLC17A6", "FOXP2", "TBR1", "NEUROD1", "NEUROD2", "FEZF2", "SLC17A7", "SLC17A6")
myMarkers[["GABAergic"]] <- c("GAD1", "GAD2", "SLC32A1", "GABRA1", "GABRA2", "GABRB3", "GABRG2", "SST", "DLX1", "DLX2", "FEZF2","VIP","NOS1", "PVALB", "SST")
myMarkers[["Glutamate.receptors"]] <- c("GRIN1","GRIN2A","GRIN2B","GRIN2C","GRIN2D","GRM1","GRM2","GRM3","GRM4","GRM5","SLC1A1","SLC1A2","SLC1A3","SLC1A6","SLC17A7","GLUL","GAD1","GAD2")
myMarkers[["neuronal.differentiation"]] <- c("SLC17A7", "SLC17A6", "TBR1", "FOXP2", "GAD1", "GAD2", "DLX1", "DLX2", "NEUROD1", "NEUROD2", "TBR2")
myMarkers[["neurogenesis.maturation"]] <- c("NEUROD1", "NEUROD2", "DLX1", "DLX2", "TBR1", "FOXP2")
myMarkers[["non.coding"]] <- c("TFAP2C", "NR2F1")
myMarkers[["RG.oRG"]] <- c("TRB2","SOX2","NES", "PAX6", "EMX2", "EOMES", "VLDLR", "FAM107A", "FOXM1", "TNC", "HOPX", "PROM1")
myMarkers[["astrocytes"]] <- c("GFAP","AQP4","SLCO1C1")


other <- c("NOG","PAX6","SIX3","OTX2","EMX2","FEZF2","GRIN2B","EOMES","BCL11B","HOPX","TNC","VLDLR","OLIG1","PDGFRA","S100B","OTX2","EMX2","FOXG1","FAM107A","DCX","EOMES","NEUROD1","ASCL1","NEUROG2","TBR1","EMX2","RELN","KIF5C","SYT1","SYP","KCNQ2","CAMK2A","DLG4","POU3F2","GAD1","RBFOX1","EOMES","TBR1","DCX","RELN","SYT1","SYP","GRIN2B","KCNQ2","CAMK2A","DLG4","GRM5","GAD1","FAM107A","TNC","FEZF2","OTX2","EMX2","FLNA","CDH2","PDGFRA","SLCO1C1","FOXG1","CSPG4","PAX6","SFRP1","FLNA","PDGFRA","S100A10","S100A11","PDGFRA","SOX2","HOPX","NOTCH1","VLDLR","GFAP","AQP4","POU3F2","S100B","ALDH1L1","VIM","SOX9","FN1","BCL2","CUX1","MKI67", "PROM1","TLX","VIM","NES","CD44","RCOR1","BRN2","SATB2","LX1", "NKX2.1,","HOXB4","ISL1","FOXG1",
"GRIN1","GRIN2A","GRIN2B","GRIN2C","GRIN2D","GRM1","GRM2","GRM3","GRM4","GRM5","SLC1A1","SLC1A2","SLC1A3","SLC1A6","SLC17A7","GLUL","GAD1","GAD2","VIP","NOS1","GAD1","GAD2","SLC32A1","GABRA1","GABRA2","GABRB3","GABRG2","VIP","NOS1","PVALB","SST","SOX1", "SOX2", "SOX3", "SOX4", "SOX5", "SOX6", "SOX7", "SOX8",  "SOX9", "SOX10", "SOX11", "SOX12", "SOX13", "SOX14", "SOX15", 
 "SOX17", "SOX18", "SOX19", "SOX20", "SOX21", "SOX22", "SOX23",  "SOX24", "SOX25", "SOX26", "SOX27", "SOX28","SFRP2","PTTG1","GLI3", "VGLUT2" )

allMarkers <- unique(unlist(myMarkers), other)

#####
xx <- read.csv(file=paste0(Fin,"PseudoBulk_Diffex_60/",  "Day60_rlogNorm.csv"), row.names=1)
x1 <- xx[which(rownames(xx) %in% allMarkers), ]
toRm <- which(apply(x1, 1, max) < 1)
if(length(toRm) > 0) {  x1 <- x1[-toRm, ] }
#colnames(x1) <- sub("BC00","BC", colnames(x1))

s1 <- sapply(strsplit(colnames(x1), split = "\\."), "[[", 1)
s2 <- sapply(strsplit(colnames(x1), split = "\\."), "[[", 2)
cc <- paste(s2, s1, sep="_")
colnames(x1) <- cc

plotlist <- list()
plotlist[["Day60_Marker_HeatMap"]] <-  pheatmap(x1, scale="row", cutree_cols=10  ,fontsize_row = 6, fontsize_col = 8)
doPlots(plotlist)


#####
xx <- read.csv(file=paste0(Fin,"PseudoBulk_Diffex_20/",  "Day20_rlogNorm.csv"), row.names=1)
x1 <- xx[which(rownames(xx) %in% allMarkers), ]
toRm <- which(apply(x1, 1, max) < 1)
if(length(toRm) > 0) {  x1 <- x1[-toRm, ] }
#colnames(x1) <- sub("BC00","BC", colnames(x1))

s1 <- sapply(strsplit(colnames(x1), split = "\\."), "[[", 1)
s2 <- sapply(strsplit(colnames(x1), split = "\\."), "[[", 2)
cc <- paste(s2, s1, sep="_")
colnames(x1) <- cc

plotlist <- list()
plotlist[["Day20_Marker_HeatMap"]] <-  pheatmap(x1, scale="row", cutree_cols=10, fontsize_row = 6, fontsize_col = 8)
doPlots(plotlist)


