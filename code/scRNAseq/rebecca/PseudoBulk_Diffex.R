

source("./scripts/SCfunctions.R")

# DirOutQC <- "./outputs/QCplots_harmony_Mt20_lambda/PseudoBulk_Diffex_20/"  
# load(file="./outputs/SeuObjs_harmony_Mt20_lambda/integrated20.RData")
# seurat.integrated <- integrated20
# key <- "Day20"

 DirOutQC <- "./outputs/QCplots_harmony_Mt20_lambda/PseudoBulk_Diffex_60/"
 load(file="./outputs/SeuObjs_harmony_Mt20_lambda/integrated60.RData")  
 seurat.integrated <- integrated60
 key <- "Day60"

clustMetaDat <- seurat.integrated@meta.data
Sample <- paste0(clustMetaDat$orig.ident, ".", clustMetaDat$Genotype)
Clust <- clustMetaDat$seurat_clusters
aDF <- data.frame(Sample,Clust)
clustMetaDat <- cbind(clustMetaDat, aDF)
seurat.integrated@meta.data <- clustMetaDat

DefaultAssay(seurat.integrated) <- "RNA"
cts <- AggregateExpression(seurat.integrated, assays = "RNA", slot = "counts", return.seurat = F, group.by =c("Sample", "Clust"))
cts <- cts$RNA

sample <- sapply(strsplit(colnames(cts), "\\."), "[[",1)
Gtype <- sapply(strsplit(colnames(cts), "\\."), "[[",2)
Clust <- paste0("Clust_",sapply(strsplit( Gtype , "_"), "[[",2))
samplesheet <- data.frame(sample,Gtype,Clust,  Group=colnames(cts), row.names=colnames(cts) )

dds <- DESeqDataSetFromMatrix(countData = cts, colData = samplesheet, design = ~ Gtype)
dds <- DESeq(dds)
rld <- rlog(dds, blind=FALSE)
rld_df = as.data.frame(assay(rld))
write.csv(rld_df, file=paste0(DirOutQC,  "Day20_rlogNorm.csv") )

# HET vs KO
plotList <- list()
     res_HET0_vs_KO0 <- as.data.frame(results(dds, contrast = c("Gtype","HET_0", "KO_0"))) ; res_HET0_vs_KO0 <- cbind(gene=rownames(res_HET0_vs_KO0), res_HET0_vs_KO0)
     write.table(res_HET0_vs_KO0[order(res_HET0_vs_KO0$pvalue),], file=paste(DirOutQC, "Clust.0_",key,"_HET_vs_KO.txt", sep=""), quote=F,sep="\t",row.names = F)
     plotList[[paste("Clust.0_HET_vs_KO_Volcano")]] <- doVolcano(res_HET0_vs_KO0, "Clust.0_HET_vs_KO", 2, 0.01)

     res_HET1_vs_KO1 <- as.data.frame(results(dds, contrast = c("Gtype","HET_1", "KO_1"))) ; res_HET1_vs_KO1 <- cbind(gene=rownames(res_HET1_vs_KO1), res_HET1_vs_KO1)
     write.table(res_HET1_vs_KO1[order(res_HET1_vs_KO1$pvalue),], file=paste(DirOutQC, "Clust.1_",key,"_HET_vs_KO.txt", sep=""), quote=F,sep="\t",row.names = F)
     plotList[[paste("Clust.1_HET_vs_KO_Volcano")]] <- doVolcano(res_HET1_vs_KO1, "Clust.1_HET_vs_KO", 2, 0.01)

     res_HET2_vs_KO2 <- as.data.frame(results(dds, contrast = c("Gtype","HET_2", "KO_2"))) ; res_HET2_vs_KO2 <- cbind(gene=rownames(res_HET2_vs_KO2), res_HET2_vs_KO2)
     write.table(res_HET2_vs_KO2[order(res_HET2_vs_KO2$pvalue),], file=paste(DirOutQC, "Clust.2_",key,"_HET_vs_KO.txt", sep=""), quote=F,sep="\t",row.names = F)
     plotList[[paste("Clust.2_HET_vs_KO_Volcano")]] <- doVolcano(res_HET2_vs_KO2, "Clust.2_HET_vs_KO", 2, 0.01)

     res_HET3_vs_KO3 <- as.data.frame(results(dds, contrast = c("Gtype","HET_3", "KO_3"))) ; res_HET3_vs_KO3 <- cbind(gene=rownames(res_HET3_vs_KO3), res_HET3_vs_KO3)
     write.table(res_HET3_vs_KO3[order(res_HET3_vs_KO3$pvalue),], file=paste(DirOutQC, "Clust.3_",key,"_HET_vs_KO.txt", sep=""), quote=F,sep="\t",row.names = F)
     plotList[[paste("Clust.3_HET_vs_KO_Volcano")]] <- doVolcano(res_HET3_vs_KO3, "Clust.3_HET_vs_KO", 2, 0.01)
     
     res_HET4_vs_KO4 <- as.data.frame(results(dds, contrast = c("Gtype","HET_4", "KO_4"))) ; res_HET4_vs_KO4 <- cbind(gene=rownames(res_HET4_vs_KO4), res_HET4_vs_KO4)
     write.table(res_HET4_vs_KO4[order(res_HET4_vs_KO4$pvalue),], file=paste(DirOutQC, "Clust.4_",key,"_HET_vs_KO.txt", sep=""), quote=F,sep="\t",row.names = F)
     plotList[[paste("Clust.4_HET_vs_KO_Volcano")]] <- doVolcano(res_HET4_vs_KO4, "Clust.4_HET_vs_KO", 2, 0.01)
doPlots(plotList)

if(key == "Day60") {
plotList <- list()

     res_HET5_vs_KO5 <- as.data.frame(results(dds, contrast = c("Gtype","HET_5", "KO_5"))) ; res_HET5_vs_KO5 <- cbind(gene=rownames(res_HET5_vs_KO5), res_HET5_vs_KO5)
     write.table(res_HET5_vs_KO5[order(res_HET5_vs_KO5$pvalue),], file=paste(DirOutQC, "Clust.5_",key,"_HET_vs_KO.txt", sep=""), quote=F,sep="\t",row.names = F)
     plotList[[paste("Clust.5_HET_vs_KO_Volcano")]] <- doVolcano(res_HET5_vs_KO5, "Clust.5_HET_vs_KO", 2, 0.01)

     res_HET6_vs_KO6 <- as.data.frame(results(dds, contrast = c("Gtype","HET_6", "KO_6"))) ; res_HET6_vs_KO6 <- cbind(gene=rownames(res_HET6_vs_KO6), res_HET6_vs_KO6)
     write.table(res_HET6_vs_KO6[order(res_HET6_vs_KO6$pvalue),], file=paste(DirOutQC, "Clust.6_",key,"_HET_vs_KO.txt", sep=""), quote=F,sep="\t",row.names = F)
     plotList[[paste("Clust.6_HET_vs_KO_Volcano")]] <- doVolcano(res_HET6_vs_KO6, "Clust.6_HET_vs_KO", 2, 0.01)

     res_HET7_vs_KO7 <- as.data.frame(results(dds, contrast = c("Gtype","HET_7", "KO_7"))) ; res_HET7_vs_KO7 <- cbind(gene=rownames(res_HET7_vs_KO7), res_HET7_vs_KO7)
     write.table(res_HET7_vs_KO7[order(res_HET7_vs_KO7$pvalue),], file=paste(DirOutQC, "Clust.7_",key,"_HET_vs_KO.txt", sep=""), quote=F,sep="\t",row.names = F)
     plotList[[paste("Clust.7_HET_vs_KO_Volcano")]] <- doVolcano(res_HET7_vs_KO7, "Clust.7_HET_vs_KO", 2, 0.01)

     doPlots(plotList)
}


# HET vs WT
plotList <- list()
     res_HET0_vs_WT0 <- as.data.frame(results(dds, contrast = c("Gtype","HET_0", "WT_0"))) ; res_HET0_vs_WT0 <- cbind(gene=rownames(res_HET0_vs_WT0), res_HET0_vs_WT0)
     write.table(res_HET0_vs_WT0[order(res_HET0_vs_WT0$pvalue),], file=paste(DirOutQC, "Clust.0_",key,"_HET_vs_WT.txt", sep=""), quote=F,sep="\t",row.names = F)
     plotList[[paste("Clust.0_HET_vs_WT_Volcano")]] <- doVolcano(res_HET0_vs_WT0, "Clust.0_HET_vs_WT", 2, 0.01)

     res_HET1_vs_WT1 <- as.data.frame(results(dds, contrast = c("Gtype","HET_1", "WT_1"))) ; res_HET1_vs_WT1 <- cbind(gene=rownames(res_HET1_vs_WT1), res_HET1_vs_WT1)
     write.table(res_HET1_vs_WT1[order(res_HET1_vs_WT1$pvalue),], file=paste(DirOutQC, "Clust.1_",key,"_HET_vs_WT.txt", sep=""), quote=F,sep="\t",row.names = F)
     plotList[[paste("Clust.1_HET_vs_WT_Volcano")]] <- doVolcano(res_HET1_vs_WT1, "Clust.1_HET_vs_WT", 2, 0.01)

     res_HET2_vs_WT2 <- as.data.frame(results(dds, contrast = c("Gtype","HET_2", "WT_2"))) ; res_HET2_vs_WT2 <- cbind(gene=rownames(res_HET2_vs_WT2), res_HET2_vs_WT2)
     write.table(res_HET2_vs_WT2[order(res_HET2_vs_WT2$pvalue),], file=paste(DirOutQC, "Clust.2_",key,"_HET_vs_WT.txt", sep=""), quote=F,sep="\t",row.names = F)
     plotList[[paste("Clust.2_HET_vs_WT_Volcano")]] <- doVolcano(res_HET2_vs_WT2, "Clust.2_HET_vs_WT", 2, 0.01)

     res_HET3_vs_WT3 <- as.data.frame(results(dds, contrast = c("Gtype","HET_3", "WT_3"))) ; res_HET3_vs_WT3 <- cbind(gene=rownames(res_HET3_vs_WT3), res_HET3_vs_WT3)
     write.table(res_HET3_vs_WT3[order(res_HET3_vs_WT3$pvalue),], file=paste(DirOutQC, "Clust.3_",key,"_HET_vs_WT.txt", sep=""), quote=F,sep="\t",row.names = F)
     plotList[[paste("Clust.3_HET_vs_WT_Volcano")]] <- doVolcano(res_HET3_vs_WT3, "Clust.3_HET_vs_WT", 2, 0.01)
     
     res_HET4_vs_WT4 <- as.data.frame(results(dds, contrast = c("Gtype","HET_4", "WT_4"))) ; res_HET4_vs_WT4 <- cbind(gene=rownames(res_HET4_vs_WT4), res_HET4_vs_WT4)
     write.table(res_HET4_vs_WT4[order(res_HET4_vs_WT4$pvalue),], file=paste(DirOutQC, "Clust.4_",key,"_HET_vs_WT.txt", sep=""), quote=F,sep="\t",row.names = F)
     plotList[[paste("Clust.4_HET_vs_WT_Volcano")]] <- doVolcano(res_HET4_vs_WT4, "Clust.4_HET_vs_WT", 2, 0.01)
doPlots(plotList)

if(key == "Day60") {
plotList <- list()

     res_HET5_vs_WT5 <- as.data.frame(results(dds, contrast = c("Gtype","HET_5", "WT_5"))) ; res_HET5_vs_WT5 <- cbind(gene=rownames(res_HET5_vs_WT5), res_HET5_vs_WT5)
     write.table(res_HET5_vs_WT5[order(res_HET5_vs_WT5$pvalue),], file=paste(DirOutQC, "Clust.5_",key,"_HET_vs_WT.txt", sep=""), quote=F,sep="\t",row.names = F)
     plotList[[paste("Clust.5_HET_vs_WT_Volcano")]] <- doVolcano(res_HET5_vs_WT5, "Clust.5_HET_vs_WT", 2, 0.01)

     res_HET6_vs_WT6 <- as.data.frame(results(dds, contrast = c("Gtype","HET_6", "WT_6"))) ; res_HET6_vs_WT6 <- cbind(gene=rownames(res_HET6_vs_WT6), res_HET6_vs_WT6)
     write.table(res_HET6_vs_WT6[order(res_HET6_vs_WT6$pvalue),], file=paste(DirOutQC, "Clust.6_",key,"_HET_vs_WT.txt", sep=""), quote=F,sep="\t",row.names = F)
     plotList[[paste("Clust.6_HET_vs_WT_Volcano")]] <- doVolcano(res_HET6_vs_WT6, "Clust.6_HET_vs_WT", 2, 0.01)

     res_HET7_vs_WT7 <- as.data.frame(results(dds, contrast = c("Gtype","HET_7", "WT_7"))) ; res_HET7_vs_WT7 <- cbind(gene=rownames(res_HET7_vs_WT7), res_HET7_vs_WT7)
     write.table(res_HET7_vs_WT7[order(res_HET7_vs_WT7$pvalue),], file=paste(DirOutQC, "Clust.7_",key,"_HET_vs_WT.txt", sep=""), quote=F,sep="\t",row.names = F)
     plotList[[paste("Clust.7_HET_vs_WT_Volcano")]] <- doVolcano(res_HET7_vs_WT7, "Clust.7_HET_vs_WT", 2, 0.01)

     doPlots(plotList)
}

# KO vs WT
plotList <- list()
     res_KO0_vs_WT0 <- as.data.frame(results(dds, contrast = c("Gtype","KO_0", "WT_0"))) ; res_KO0_vs_WT0 <- cbind(gene=rownames(res_KO0_vs_WT0), res_KO0_vs_WT0)
     write.table(res_KO0_vs_WT0[order(res_KO0_vs_WT0$pvalue),], file=paste(DirOutQC, "Clust.0_",key,"_KO_vs_WT.txt", sep=""), quote=F,sep="\t",row.names = F)
     plotList[[paste("Clust.0_KO_vs_WT_Volcano")]] <- doVolcano(res_KO0_vs_WT0, "Clust.0_KO_vs_WT", 2, 0.01)

     res_KO1_vs_WT1 <- as.data.frame(results(dds, contrast = c("Gtype","KO_1", "WT_1"))) ; res_KO1_vs_WT1 <- cbind(gene=rownames(res_KO1_vs_WT1), res_KO1_vs_WT1)
     write.table(res_KO1_vs_WT1[order(res_KO1_vs_WT1$pvalue),], file=paste(DirOutQC, "Clust.1_",key,"_KO_vs_WT.txt", sep=""), quote=F,sep="\t",row.names = F)
     plotList[[paste("Clust.1_KO_vs_WT_Volcano")]] <- doVolcano(res_KO1_vs_WT1, "Clust.1_KO_vs_WT", 2, 0.01)

     res_KO2_vs_WT2 <- as.data.frame(results(dds, contrast = c("Gtype","KO_2", "WT_2"))) ; res_KO2_vs_WT2 <- cbind(gene=rownames(res_KO2_vs_WT2), res_KO2_vs_WT2)
     write.table(res_KO2_vs_WT2[order(res_KO2_vs_WT2$pvalue),], file=paste(DirOutQC, "Clust.2_",key,"_KO_vs_WT.txt", sep=""), quote=F,sep="\t",row.names = F)
     plotList[[paste("Clust.2_KO_vs_WT_Volcano")]] <- doVolcano(res_KO2_vs_WT2, "Clust.2_KO_vs_WT", 2, 0.01)

     res_KO3_vs_WT3 <- as.data.frame(results(dds, contrast = c("Gtype","KO_3", "WT_3"))) ; res_KO3_vs_WT3 <- cbind(gene=rownames(res_KO3_vs_WT3), res_KO3_vs_WT3)
     write.table(res_KO3_vs_WT3[order(res_KO3_vs_WT3$pvalue),], file=paste(DirOutQC, "Clust.3_",key,"_KO_vs_WT.txt", sep=""), quote=F,sep="\t",row.names = F)
     plotList[[paste("Clust.3_KO_vs_WT_Volcano")]] <- doVolcano(res_KO3_vs_WT3, "Clust.3_KO_vs_WT", 2, 0.01)
     
     res_KO4_vs_WT4 <- as.data.frame(results(dds, contrast = c("Gtype","KO_4", "WT_4"))) ; res_KO4_vs_WT4 <- cbind(gene=rownames(res_KO4_vs_WT4), res_KO4_vs_WT4)
     write.table(res_KO4_vs_WT4[order(res_KO4_vs_WT4$pvalue),], file=paste(DirOutQC, "Clust.4_",key,"_KO_vs_WT.txt", sep=""), quote=F,sep="\t",row.names = F)
     plotList[[paste("Clust.4_KO_vs_WT_Volcano")]] <- doVolcano(res_KO4_vs_WT4, "Clust.4_KO_vs_WT", 2, 0.01)
doPlots(plotList)

if(key == "Day60") {
plotList <- list()

     res_KO5_vs_WT5 <- as.data.frame(results(dds, contrast = c("Gtype","KO_5", "WT_5"))) ; res_KO5_vs_WT5 <- cbind(gene=rownames(res_KO5_vs_WT5), res_KO5_vs_WT5)
     write.table(res_KO5_vs_WT5[order(res_KO5_vs_WT5$pvalue),], file=paste(DirOutQC, "Clust.5_",key,"_KO_vs_WT.txt", sep=""), quote=F,sep="\t",row.names = F)
     plotList[[paste("Clust.5_KO_vs_WT_Volcano")]] <- doVolcano(res_KO5_vs_WT5, "Clust.5_KO_vs_WT", 2, 0.01)

     res_KO6_vs_WT6 <- as.data.frame(results(dds, contrast = c("Gtype","KO_6", "WT_6"))) ; res_KO6_vs_WT6 <- cbind(gene=rownames(res_KO6_vs_WT6), res_KO6_vs_WT6)
     write.table(res_KO6_vs_WT6[order(res_KO6_vs_WT6$pvalue),], file=paste(DirOutQC, "Clust.6_",key,"_KO_vs_WT.txt", sep=""), quote=F,sep="\t",row.names = F)
     plotList[[paste("Clust.6_KO_vs_WT_Volcano")]] <- doVolcano(res_KO6_vs_WT6, "Clust.6_KO_vs_WT", 2, 0.01)

     res_KO7_vs_WT7 <- as.data.frame(results(dds, contrast = c("Gtype","KO_7", "WT_7"))) ; res_KO7_vs_WT7 <- cbind(gene=rownames(res_KO7_vs_WT7), res_KO7_vs_WT7)
     write.table(res_KO7_vs_WT7[order(res_KO7_vs_WT7$pvalue),], file=paste(DirOutQC, "Clust.7_",key,"_KO_vs_WT.txt", sep=""), quote=F,sep="\t",row.names = F)
     plotList[[paste("Clust.7_KO_vs_WT_Volcano")]] <- doVolcano(res_KO7_vs_WT7, "Clust.7_KO_vs_WT", 2, 0.01)

     doPlots(plotList)
}




# HET vs KO # HET vs WT     KO vs WT




