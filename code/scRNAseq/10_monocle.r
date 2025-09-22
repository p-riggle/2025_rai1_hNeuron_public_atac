# load seurat
library(Seurat)

install.packages("devtools")

BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'ggrastr'))

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("beachmat")

remotes::install_github("bnprks/BPCells/r")


