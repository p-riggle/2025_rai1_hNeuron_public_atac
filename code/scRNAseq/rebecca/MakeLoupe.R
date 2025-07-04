
# singularity exec -B /run/user/114119227 /nfs/mm-isilon/bioinfcore/Common/singularity/single_cell_0.12.0.sif R
#  singularity exec -B /run/user/114119227 /nfs/mm-isilon/bioinfcore/Common/singularity/single_cell_0.11.2.sif R
# loupeR_for_single_cell_env_0.12.sif
# singularity exec -B /run/user/114119227 /nfs/mm-isilon/bioinfcore/Common/singularity/loupeR_for_single_cell_env_0.12.sif R
 
# singularity exec -B /run/user/114119227 /nfs/mm-isilon/bioinfcore/Common/singularity/loupeR_for_single_cell_env_0.12.sif R
#library("loupeR")
  


# https://www.youtube.com/watch?v=zKZ03mA3aig
# conda create -n Seurat5.1  r-seurat
#  latest Core env single_cell_env_0.11.2.sif has   Seurat_5.0.1 

#source("./scripts/SCfunctions.R")

library(Seurat)
library(loupeR)

DirOutQC <- "./outputs/QCplots/"

load(file=paste0(DirOutQC,"integrated20.RData"))
create_loupe_from_seurat(integrated20, output_dir=DirOutQC, output_name="Day.20")

load(file=paste0(DirOutQC,"integrated60.RData"))
create_loupe_from_seurat(integrated60, output_dir=DirOutQC, output_name="Day.60")






