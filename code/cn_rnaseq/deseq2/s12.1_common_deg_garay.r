# check the overlap with Garay et al 2020
# import garay_deg (from Cell Reports 2020 paper, supplementary table 4, vehicle tab)
restbl_garay <- read.table("results/common_deg/garay_restbl_veh.txt", sep = "\t", header = TRUE)

# import restbls from all three datasets
restbl_huang <- read.table("results/deseq2/huang_2016/huang_2016_rstbl.txt", sep = "\t", header = TRUE, row.names = 1)
restbl_kaster_d8 <- read.table("results/deseq2/kaster_2023/d8_wt.vs.ko.txt", sep = "\t", header = TRUE, row.names = 1)
restbl_kaster_d12 <- read.table("results/deseq2/kaster_2023/d12_wt.vs.ko.txt", sep = "\t", header = TRUE, row.names = 1)
restbl_zhou_d20 <- read.table("results/deseq2/rstbl_d20_wt.vs.ko.txt", sep = "\t", header = TRUE, row.names = 1)
restbl_zhou_d34 <- read.table("results/deseq2/rstbl_d34_wt.vs.ko.txt", sep = "\t", header = TRUE, row.names = 1)
restbl_zhou_d60 <- read.table("results/deseq2/rstbl_d60_wt.vs.ko.txt", sep = "\t", header = TRUE, row.names = 1)

# add symbol column to zho data
restbl_zhou_d20$symbol <- rownames(restbl_zhou_d20)
restbl_zhou_d34$symbol <- rownames(restbl_zhou_d34)
restbl_zhou_d60$symbol <- rownames(restbl_zhou_d60)

# select degs from each rstbl 
updeg_huang <- subset(restbl_huang, padj < 0.2 & log2FoldChange > 0)
downdeg_huang <- subset(restbl_huang, padj < 0.2 & log2FoldChange < 0)
updeg_kaster_d8 <- subset(restbl_kaster_d8, padj < 0.2 & log2FoldChange > 0)
updeg_kaster_d12 <- subset(restbl_kaster_d12, padj < 0.2 & log2FoldChange > 0)
downdeg_kaster_d8 <- subset(restbl_kaster_d8, padj < 0.2 & log2FoldChange < 0)
downdeg_kaster_d12 <- subset(restbl_kaster_d12, padj < 0.2 & log2FoldChange < 0)
updeg_zhou_d20 <- subset(restbl_zhou_d20, padj < 0.05 & log2FoldChange > 1)
updeg_zhou_d34 <- subset(restbl_zhou_d34, padj < 0.05 & log2FoldChange > 1)
updeg_zhou_d60 <- subset(restbl_zhou_d60, padj < 0.05 & log2FoldChange > 1)
downdeg_zhou_d20 <- subset(restbl_zhou_d20, padj < 0.05 & log2FoldChange < -1)
downdeg_zhou_d34 <- subset(restbl_zhou_d34, padj < 0.05 & log2FoldChange < -1)
downdeg_zhou_d60 <- subset(restbl_zhou_d60, padj < 0.05 & log2FoldChange < -1)
updeg_garay <- subset(restbl_garay, padj < 0.2 & log2FoldChange > 0)
downdeg_garay <- subset(restbl_garay, padj < 0.2 & log2FoldChange < 0)
nrow(updeg_garay)
nrow(downdeg_garay)

# select degs with relaxed p value cutoff

updeg_zhou_d20_2 <- subset(restbl_zhou_d20, padj < 0.05 & log2FoldChange > 0)
updeg_zhou_d34_2 <- subset(restbl_zhou_d34, padj < 0.05 & log2FoldChange > 0)
updeg_zhou_d60_2 <- subset(restbl_zhou_d60, padj < 0.05 & log2FoldChange > 0)
downdeg_zhou_d20_2 <- subset(restbl_zhou_d20, padj < 0.05 & log2FoldChange < 0)
downdeg_zhou_d34_2 <- subset(restbl_zhou_d34, padj < 0.05 & log2FoldChange < 0)
downdeg_zhou_d60_2 <- subset(restbl_zhou_d60, padj < 0.05 & log2FoldChange < 0)

# combine multiple timepoints to one dataframe
updeg_kaster <- rbind(updeg_kaster_d8, updeg_kaster_d12)
downdeg_kaster <- rbind(downdeg_kaster_d8, downdeg_kaster_d12)
updeg_zhou <- rbind(updeg_zhou_d20, updeg_zhou_d34, updeg_zhou_d60) 
downdeg_zhou <- rbind(downdeg_zhou_d20, downdeg_zhou_d34, downdeg_zhou_d60)
downdeg_zhou_2 <- rbind(downdeg_zhou_d20_2, downdeg_zhou_d34_2, downdeg_zhou_d60_2)
updeg_zhou_2 <- rbind(updeg_zhou_d20_2, updeg_zhou_d34_2, updeg_zhou_d60_2)

# remove duplicated genes
updeg_kaster <- updeg_kaster[!duplicated(updeg_kaster$symbol), ]
downdeg_kaster <- downdeg_kaster[!duplicated(downdeg_kaster$symbol), ]
updeg_zhou <- updeg_zhou[!duplicated(updeg_zhou$symbol), ]
downdeg_zhou <- downdeg_zhou[!duplicated(downdeg_zhou$symbol), ]
updeg_zhou_2 <- updeg_zhou_2[!duplicated(updeg_zhou_2$symbol), ]
downdeg_zhou_2 <- downdeg_zhou_2[!duplicated(downdeg_zhou_2$symbol), ]

# remove NA
updeg_huang <- subset(updeg_huang, !is.na(symbol))
downdeg_huang <- subset(downdeg_huang, !is.na(symbol))
updeg_kaster <- subset(updeg_kaster, !is.na(symbol))
downdeg_kaster <- subset(downdeg_kaster, !is.na(symbol))
updeg_zhou <- subset(updeg_zhou, !is.na(symbol))
downdeg_zhou <- subset(downdeg_zhou, !is.na(symbol))
updeg_zhou_2 <- subset(updeg_zhou_2, !is.na(symbol))
downdeg_zhou_2 <- subset(downdeg_zhou_2, !is.na(symbol))

# make venn diagrams
# upregulated genes huang, kaster, and garay
library("BioVenn")

updeg_huang_symbol <- updeg_huang$symbol
updeg_kaster_symbol <- updeg_kaster$symbol
updeg_garay_symbol <- updeg_garay$symbol
updeg_garay_symbol 

venn_up_hkg <- draw.venn(updeg_huang_symbol, updeg_kaster_symbol, updeg_garay_symbol, 
title = "Common DEGs (up)", t_f = "Helvetica",
subtitle = NULL,
x_c = "#000096", xtitle = "Huang et al. 2016", xt_f = "Helvetica",
y_c = "yellow", ytitle = "Kaster et al. 2023", yt_f = "Helvetica",
z_c = "red", ztitle = "Garay et al. 2020", zt_f = "Helvetica",
nrtype = "abs", nr_f = "Helvetica",
output = "pdf", filename = sprintf("results/common_deg/%s-common_updegs_hkg_overlap.pdf", Sys.Date()))

# downregulated genes huang, kaster, and garay
downdeg_huang_symbol <- downdeg_huang$symbol
downdeg_kaster_symbol <- downdeg_kaster$symbol
downdeg_garay_symbol <- downdeg_garay$symbol
downdeg_garay_symbol

venn_down_hkg <- draw.venn(downdeg_huang_symbol, downdeg_kaster_symbol, downdeg_garay_symbol,
title = "Common DEGs (down)", t_f = "Helvetica",
subtitle = NULL,
x_c = "#000096", xtitle = "Huang et al. 2016", xt_f = "Helvetica",
y_c = "yellow", ytitle = "Kaster et al. 2023", yt_f = "Helvetica",
z_c = "red", ztitle = "Garay et al. 2020", zt_f = "Helvetica",
nrtype = "abs", nr_f = "Helvetica",
output = "pdf", filename = sprintf("results/common_deg/%s-common_downdegs_hkg_overlap.pdf", Sys.Date()))

# upregulated genes huang, zhou, and garay
# convert zhou data HUMAN gene symbols to mouse ones
library(orthogene)
updeg_zhou_m <- orthogene::convert_orthologs(gene_df = updeg_zhou,
                                        gene_input = "rownames", 
                                        gene_output = "rownames", 
                                        input_species = "human",
                                        output_species = "mouse",
                                        non121_strategy = "drop_both_species")

downdeg_zhou_m <- orthogene::convert_orthologs(gene_df = downdeg_zhou,
                                        gene_input = "rownames", 
                                        gene_output = "rownames", 
                                        input_species = "human",
                                        output_species = "mouse",
                                        non121_strategy = "drop_both_species")

updeg_zhou_symbol <- rownames(updeg_zhou_m)
updeg_zhou_symbol

venn_up_hgz <- draw.venn(updeg_huang_symbol, updeg_garay_symbol, updeg_zhou_symbol,
title = "Common DEGs (up)", t_f = "Helvetica",
subtitle = NULL,
x_c = "#000096", xtitle = "Huang et al. 2016", xt_f = "Helvetica",
y_c = "yellow", ytitle = "Garay et al. 2020", yt_f = "Helvetica",
z_c = "red", ztitle = "Zhou et al. 2023", zt_f = "Helvetica",
nrtype = "abs", nr_f = "Helvetica",
output = "pdf", filename = sprintf("results/common_deg/%s-common_updegs_hgz_overlap.pdf", Sys.Date()))

# downregulated genes huang, zhou, and garay
downdeg_zhou_symbol <- rownames(downdeg_zhou_m)
downdeg_zhou_symbol

venn_down_hgz <- draw.venn(downdeg_huang_symbol, downdeg_garay_symbol, downdeg_zhou_symbol,
title = "Common DEGs (down)", t_f = "Helvetica",
subtitle = NULL,
x_c = "#000096", xtitle = "Huang et al. 2016", xt_f = "Helvetica",
y_c = "yellow", ytitle = "Garay et al. 2020", yt_f = "Helvetica",
z_c = "red", ztitle = "Zhou et al. 2023", zt_f = "Helvetica",
nrtype = "abs", nr_f = "Helvetica",
output = "pdf", filename = sprintf("results/common_deg/%s-common_downdegs_hgz_overlap.pdf", Sys.Date()))


