# explore what RData files from Rebbeca has and find gene count tables to do DEG cluster analyses
# module load R
# R

load("data/scRNAseq/rebecca/SeuObjs_harmony_Mt20_lambda/integrated20.RData")

head(integrated20$counts)
