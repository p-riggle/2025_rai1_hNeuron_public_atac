
uniqList <- function(gene_list) {
   gene_counts <- table(unlist(gene_list))
   unique_genes <- names(gene_counts[gene_counts == 1])
   filtered_gene_list <- lapply(gene_list, function(gene_set) {
     gene_set[gene_set %in% unique_genes]
   })
   return(filtered_gene_list)
}

inDir <- "./inputs/ZhongMarkers/"

ZhongMajorList <- list()
Major <- read.csv(file=paste0(inDir, "major.csv"), row.names=1)
Major$cluster <- sub(" ",".", Major$cluster)
MajorPos <- Major[which(Major$avg_diff > 0),]
uniqClust <- unique(MajorPos$cluster)
for(ii in 1:length(uniqClust)) {
    clust <- uniqClust[ii]
    ZhongMajorList[[clust]] <- rownames(MajorPos)[which(MajorPos$cluster == clust)]
}
sapply(ZhongMajorList, function(x) length(x))
#       NPCs Excitatory.neurons       Interneurons                OPC 
#        142                  8                 18                 91 
# Astrocytes          Microglia 
#        326                438 
uniqZhongMajorList <- uniqList(ZhongMajorList)
sapply(uniqZhongMajorList, function(x) length(x))  # there were no overlaps

fins <- c("Astrocytes.csv","Excitatory.csv","ExcitatoryStage.csv","ExcitatoryWEEK.csv","Interneurons.csv","MicroGlia.csv","NPCs.csv","OPCs.csv")
names <- c("Astrocytes","Excitatory","ExcitatoryStage","ExcitatoryWEEK","Interneurons","MicroGlia","NPCs","OPCs")

ZhongPosList <- list()
for(ii in 1:length(fins)) {
    clust <- names[ii]
    aDF <- read.csv(file=paste0(inDir, fins[ii]))
    aDF <- aDF[which(aDF$avg_diff > 0),]
    ZhongPosList[[clust]] <- unique(aDF$gene)     # unique(rownames(aDF))
}
sapply(ZhongPosList, function(x) length(x))
# Astrocytes      Excitatory ExcitatoryStage  ExcitatoryWEEK    Interneurons 
#        431              61              22             156             120 
#  MicroGlia            NPCs            OPCs 
#        750             389             421 

uniqZhongPosList <- uniqList(ZhongPosList)
sapply(uniqZhongPosList, function(x) length(x))
# Astrocytes      Excitatory ExcitatoryStage  ExcitatoryWEEK    Interneurons 
#        286               5               1              89              54 
#  MicroGlia            NPCs            OPCs 
#        478             209             247 



ZhongPosSubClustList <- list()
uniqZhongPosSubClustList <- list()
for(ii in 1:length(fins)) {
    clust <- names[ii]
    aDF <- read.csv(file=paste0(inDir, fins[ii]))   #, row.names=1)
    aDF <- aDF[which(aDF$avg_diff > 0),]
    uniqSubClust <- unique(aDF$cluster)
    for(jj in 1:length(uniqSubClust)) {
        clust1 <- paste0(clust,".",uniqSubClust[jj])
        aDF1 <- aDF[which(aDF$cluster == uniqSubClust[jj]), ]
	ZhongPosSubClustList[[clust1]] <- unique(aDF1$gene)
    }
    uniqZhongPosSubClustList[[clust1]] <- uniqList(ZhongPosSubClustList)
}

sapply(ZhongPosSubClustList,  function(x) length(x))
sapply(uniqZhongPosSubClustList,  function(x) length(x))




save(ZhongMajorList, ZhongPosList, uniqZhongPosList, ZhongPosSubClustList, uniqZhongPosSubClustList, file="./inputs/ZhongLists.RData")

#ZhongPosSubClustList
#uniqZhongPosSubClustList

#ZhongPosList
#uniqZhongPosList

#ZhongMajorList
