# select common degs from the three datasets

# import restbls from all three datasets
restbl_huang <- read.table("results/deseq2/huang_2016/huang_2016_rstbl.txt", sep = "\t", header = TRUE, row.names = 1)
restbl_kaster_d8 <- read.table("results/deseq2/kaster_2023/d8_wt.vs.ko.txt", sep = "\t", header = TRUE, row.names = 1)
restbl_kaster_d12 <- read.table("results/deseq2/kaster_2023/d12_wt.vs.ko.txt", sep = "\t", header = TRUE, row.names = 1)
restbl_zhou_d20 <- read.table("results/deseq2/rstbl_d20_wt.vs.ko.txt", sep = "\t", header = TRUE, row.names = 1)
restbl_zhou_d34 <- read.table("results/deseq2/rstbl_d34_wt.vs.ko.txt", sep = "\t", header = TRUE, row.names = 1)
restbl_zhou_d60 <- read.table("results/deseq2/rstbl_d60_wt.vs.ko.txt", sep = "\t", header = TRUE, row.names = 1)

# add symbol column to zho data
restbl_zhou_d20$symbol <- rownames(restbl_zhou_d20)
head(restbl_zhou_d20)
restbl_zhou_d34$symbol <- rownames(restbl_zhou_d34)
restbl_zhou_d60$symbol <- rownames(restbl_zhou_d60)

# select degs from each rstbl 
head(restbl_huang)
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
head(restbl_zhou_d20)
head(downdeg_zhou_d20)
downdeg_zhou_d34 <- subset(restbl_zhou_d34, padj < 0.05 & log2FoldChange < -1)
downdeg_zhou_d60 <- subset(restbl_zhou_d60, padj < 0.05 & log2FoldChange < -1)


# combine multiple timepoints to one dataframe
updeg_kaster <- rbind(updeg_kaster_d8, updeg_kaster_d12)
downdeg_kaster <- rbind(downdeg_kaster_d8, downdeg_kaster_d12)
updeg_zhou <- rbind(updeg_zhou_d20, updeg_zhou_d34, updeg_zhou_d60) 
downdeg_zhou <- rbind(downdeg_zhou_d20, downdeg_zhou_d34, downdeg_zhou_d60)

# remove duplicated genes
updeg_kaster <- updeg_kaster[!duplicated(updeg_kaster$symbol), ]
downdeg_kaster <- downdeg_kaster[!duplicated(downdeg_kaster$symbol), ]
updeg_zhou <- updeg_zhou[!duplicated(updeg_zhou$symbol), ]
downdeg_zhou <- downdeg_zhou[!duplicated(downdeg_zhou$symbol), ]

head(updeg_zhou)
head(downdeg_zhou)

# remove NA
updeg_huang <- subset(updeg_huang, !is.na(symbol))
downdeg_huang <- subset(downdeg_huang, !is.na(symbol))
updeg_kaster <- subset(updeg_kaster, !is.na(symbol))
downdeg_kaster <- subset(downdeg_kaster, !is.na(symbol))
updeg_zhou <- subset(updeg_zhou, !is.na(symbol))
downdeg_zhou <- subset(downdeg_zhou, !is.na(symbol))

# make a bar plot of deg numbers
n_kaster <- nrow(updeg_kaster) + nrow(downdeg_kaster)
n_kaster
n_huang <- nrow(updeg_huang) + nrow(downdeg_huang)
n_huang
n_zhou <- nrow(updeg_zhou) + nrow(downdeg_zhou)
n_zhou

# make a barplot of deg numbers obtained with different p value cutoffs
library(ggplot2)    
library(reshape2)

# create a data frame
deg_num <- data.frame(dataset = c("Kaster et al.", "Huang et al.", "Zhou et al,"), down = c(nrow(downdeg_kaster), nrow(downdeg_huang), nrow(downdeg_zhou)), up = c(nrow(updeg_kaster), nrow(updeg_huang), nrow(updeg_zhou)))

deg_num

# Reshape the data into long format
df_long <- reshape2::melt(deg_num, id.vars = "dataset", variable.name = "gene", value.name = "count")

df_long

# Create the stacked bar plot
barplot_degnum <- ggplot(df_long, aes(x = dataset, y = count, fill = gene)) +
  geom_col() +
  labs(x = "Dataset", y = "Number of DEGs", title = "Number of DEGs") +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e"), labels = c("Down", "Up")) +
  theme_bw()

barplot_degnum

# save the plot
ggsave(filename = "results/common_deg/deg_num_barplot.pdf", plot = barplot_degnum, width = 6, height = 4.5, units = "in")

# creat a stacked barplot with relaxed p value cutoff of zhou data

# select degs with relaxed p value cutoff

updeg_zhou_d20_2 <- subset(restbl_zhou_d20, padj < 0.05 & log2FoldChange > 0)
updeg_zhou_d34_2 <- subset(restbl_zhou_d34, padj < 0.05 & log2FoldChange > 0)
updeg_zhou_d60_2 <- subset(restbl_zhou_d60, padj < 0.05 & log2FoldChange > 0)
downdeg_zhou_d20_2 <- subset(restbl_zhou_d20, padj < 0.05 & log2FoldChange < 0)
downdeg_zhou_d34_2 <- subset(restbl_zhou_d34, padj < 0.05 & log2FoldChange < 0)
downdeg_zhou_d60_2 <- subset(restbl_zhou_d60, padj < 0.05 & log2FoldChange < 0)

downdeg_zhou_2 <- rbind(downdeg_zhou_d20_2, downdeg_zhou_d34_2, downdeg_zhou_d60_2)
updeg_zhou_2 <- rbind(updeg_zhou_d20_2, updeg_zhou_d34_2, updeg_zhou_d60_2)

updeg_zhou_2 <- updeg_zhou_2[!duplicated(updeg_zhou_2$symbol), ]
downdeg_zhou_2 <- downdeg_zhou_2[!duplicated(downdeg_zhou_2$symbol), ]

# remove na in zhou data symbol column
updeg_zhou_2 <- subset(updeg_zhou_2, !is.na(symbol))
downdeg_zhou_2 <- subset(downdeg_zhou_2, !is.na(symbol))

# create a data frane for deg number bar plot
deg_num2 <- data.frame(dataset = c("Kaster et al.", "Huang et al.", "Zhou et al,"), down = c(nrow(downdeg_kaster), nrow(downdeg_huang), nrow(downdeg_zhou_2)), up = c(nrow(updeg_kaster), nrow(updeg_huang), nrow(updeg_zhou_2)))

deg_num2

# Reshape the data into long format
df_long2 <- reshape2::melt(deg_num2, id.vars = "dataset", variable.name = "gene", value.name = "count")

df_long2

# Create the stacked bar plot
barplot_degnum2 <- ggplot(df_long2, aes(x = dataset, y = count, fill = gene)) +
  geom_col() +
  labs(x = "Dataset", y = "Number of DEGs", title = "Number of DEGs (Zhou et al. |log2fc| > 0)") +
  scale_fill_manual(values = c("#1f77b4", "#ff7f0e"), labels = c("Down", "Up")) +
  theme_bw()

barplot_degnum2

# save the plot
ggsave(filename = "results/common_deg/deg_num_barplot2.pdf", plot = barplot_degnum2, width = 6, height = 4.5, units = "in")

# save the plot
ggsave(filename = "results/common_deg/deg_num_barplot.pdf", plot = barplot_degnum, width = 6, height = 4.5, units = "in")

# convert human gene symbols to mouse ones in zhou data
# https://neurogenomics.github.io/orthogene/reference/convert_orthologs.html

# if (!requireNamespace("BiocManager", quietly = TRUE))
#      install.packages("BiocManager")
# # orthogene is only available on Bioconductor>=3.14
# if(BiocManager::version()<"3.14") BiocManager::install(version = "3.14")
# BiocManager::install("orthogene")

library(orthogene)

# practice run
data("exp_mouse")
method <- "homologene"

gene_df <- orthogene::convert_orthologs(gene_df = exp_mouse,
                                        gene_input = "rownames", 
                                        gene_output = "rownames", 
                                        input_species = "mouse",
                                        output_species = "human",
                                        non121_strategy = "drop_both_species",
                                        method = method) 

head(exp_mouse) 
head(gene_df)

# convert zhou data HUMAN gene symbols to mouse ones
nrow(updeg_zhou)
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


# make venn diagram
# install.packages("BioVenn", dependencies = TRUE)
library("BioVenn")

updeg_huang_symbol <- updeg_huang$symbol
updeg_kaster_symbol <- updeg_kaster$symbol
updeg_zhou_symbol <- rownames(updeg_zhou_m)
updeg_zhou_symbol

venn_up <- draw.venn(updeg_huang_symbol, updeg_kaster_symbol, updeg_zhou_symbol,
title = "Common DEGs (up)", t_f = "Helvetica",
subtitle = NULL,
x_c = "#000096", xtitle = "Huang et al. 2016", xt_f = "Helvetica",
y_c = "yellow", ytitle = "Kaster et al. 2023", yt_f = "Helvetica",
z_c = "red", ztitle = "Zhou et al. 2023", zt_f = "Helvetica",
nrtype = "abs", nr_f = "Helvetica",
output = "pdf", filename = sprintf("results/common_deg/%s-common_updegs_overlap.pdf", Sys.Date()))

downdeg_huang_symbol <- downdeg_huang$symbol
downdeg_kaster_symbol <- downdeg_kaster$symbol
downdeg_zhou_symbol <- rownames(downdeg_zhou_m)

huang_zhou_up <- venn_up$xz_only
kaster_zhou_up <- venn_up$yz_only
common_up <- c(huang_zhou_up, kaster_zhou_up)
common_up

venn_down <- draw.venn(downdeg_huang_symbol, downdeg_kaster_symbol, downdeg_zhou_symbol,
title = "Common DEGs (down)", t_f = "Helvetica",
subtitle = NULL,
x_c = "#000096", xtitle = "Huang et al. 2016", xt_f = "Helvetica",
y_c = "yellow", ytitle = "Kaster et al. 2023", yt_f = "Helvetica",
z_c = "red", ztitle = "Zhou et al. 2023", zt_f = "Helvetica",
nrtype = "abs", nr_f = "Helvetica",
output = "pdf", filename = sprintf("results/common_deg/%s-common_downdegs_overlap.pdf", Sys.Date()))

huang_zhou_down <- venn_down$xz_only
kaster_zhou_down <- venn_down$yz_only
common_down <- c(huang_zhou_down, kaster_zhou_down)
common_down

# make venn with relaxed p value cutoff with zhou data.

# convert human ortholog to mouse ones
updeg_zhou_2m <- orthogene::convert_orthologs(gene_df = updeg_zhou_2,
                                        gene_input = "rownames", 
                                        gene_output = "rownames", 
                                        input_species = "human",
                                        output_species = "mouse",
                                        non121_strategy = "drop_both_species")
downdeg_zhou_2m <- orthogene::convert_orthologs(gene_df = downdeg_zhou_2,
                                        gene_input = "rownames", 
                                        gene_output = "rownames", 
                                        input_species = "human",
                                        output_species = "mouse",
                                        non121_strategy = "drop_both_species")

# make venn diagram
updeg_zhou_symbol_2 <- rownames(updeg_zhou_2m)
downdeg_zhou_symbol_2 <- rownames(downdeg_zhou_2m)

venn_up_2 <- draw.venn(updeg_huang_symbol, updeg_kaster_symbol, updeg_zhou_symbol_2,
title = "Common upDEGs Zhou log2FC > 0", t_f = "Helvetica",
subtitle = NULL,
x_c = "#000096", xtitle = "Huang et al. 2016", xt_f = "Helvetica",
y_c = "yellow", ytitle = "Kaster et al. 2023", yt_f = "Helvetica",
z_c = "red", ztitle = "Zhou et al. 2023", zt_f = "Helvetica",
nrtype = "abs", nr_f = "Helvetica",
output = "pdf", filename = sprintf("results/common_deg/%s-common_updegs_overlap2.pdf", Sys.Date()))

venn_up_2$xyz

# extract overlap genes
huang_zhou_up_2 <- venn_up_2$xz_only
kaster_zhou_up_2 <- venn_up_2$yz_only
common_up_2 <- c(huang_zhou_up_2, kaster_zhou_up_2)
huang_zhou_up_2
kaster_zhou_up_2
common_up_2

venn_down_2 <- draw.venn(downdeg_huang_symbol, downdeg_kaster_symbol, downdeg_zhou_symbol_2,
title = "Common downDEGs Zhou log2FC < 0", t_f = "Helvetica",
subtitle = NULL,
x_c = "#000096", xtitle = "Huang et al. 2016", xt_f = "Helvetica",
y_c = "yellow", ytitle = "Kaster et al. 2023", yt_f = "Helvetica",
z_c = "red", ztitle = "Zhou et al. 2023", zt_f = "Helvetica",
nrtype = "abs", nr_f = "Helvetica",
output = "pdf", filename = sprintf("results/common_deg/%s-common_downdegs_overlap2.pdf", Sys.Date()))

# extract overlap genes
huang_zhou_down_2 <- venn_down_2$xz_only
kaster_zhou_down_2 <- venn_down_2$yz_only
common_down_2 <- c(huang_zhou_down_2, kaster_zhou_down_2)

common_down_2

# pathway analysis of common degs
#combine the common deg lists into a list
class(huang_zhou_up)
class(kaster_zhou_up)
up <- list(kaster_zhou_up, huang_zhou_up)
down <- list(kaster_zhou_down, huang_zhou_down)

# name the lists
names(up) <- c("Kaster et al. 2023", "Huang et al. 2016")
names(down) <- c("Kaster et al. 2023", "Huang et al. 2016")

# convert symbols to entrez ids
library(org.Mm.eg.db)

# run clusterprofiler GO enrichment analysis
library(clusterProfiler)

str(up)
str(down)

gobp_up <- compareCluster(geneCluster = up, fun = enrichGO, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP")

gobp_down <- compareCluster(geneCluster = down, fun = enrichGO, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "BP")

# make dotplots
dplot_bp_up <- dotplot(gobp_up, showCategory = 10, title = "GO-BP enrichment of common upDEGs")
dplot_bp_up

dplot_bp_down <- dotplot(gobp_down, showCategory = 10, title = "GO-BP enrichment of common downDEGs")
dplot_bp_down

# save the two dotplots arranged side by side
library(gridExtra)
pdf(sprintf("results/common_deg/%s_gobp_dotplot.pdf", Sys.Date()), width = 15, height = 10)
grid.arrange(dplot_bp_up, dplot_bp_down, ncol = 2)
dev.off()

# run go cellular component enrichment analysis
gocc_up <- compareCluster(geneCluster = up, fun = enrichGO, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "CC")

gocc_down <- compareCluster(geneCluster = down, fun = enrichGO, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "CC")

# make dotplots
dplot_cc_up <- dotplot(gocc_up, showCategory = 10, title = "GO enrichment of common upDEGs")
dplot_cc_down <- dotplot(gocc_down, showCategory = 10, title = "GO enrichment of common downDEGs")

# save the two dotplots arranged side by side
pdf(sprintf("results/common_deg/%s_gocc_dotplot.pdf", Sys.Date()), width = 15, height = 10)
grid.arrange(dplot_cc_up, dplot_cc_down, ncol = 2)
dev.off()

# run go molecular function enrichment analysis
gomf_up <- compareCluster(geneCluster = up, fun = enrichGO, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "MF")
gomf_down <- compareCluster(geneCluster = down, fun = enrichGO, OrgDb = org.Mm.eg.db, keyType = "SYMBOL", ont = "MF")

# make dotplots
dplot_mf_up <- dotplot(gomf_up, showCategory = 10, title = "GO enrichment of common upDEGs")
dplot_mf_down <- dotplot(gomf_down, showCategory = 10, title = "GO enrichment of common downDEGs")

# save the two dotplots arranged side by side
pdf(sprintf("results/common_deg/%s_gomf_dotplot.pdf", Sys.Date()), width = 15, height = 10)
grid.arrange(dplot_mf_up, dplot_mf_down, ncol = 2)
dev.off()

# run clusterprofiler KEGG enrichment analysis
# first, convert mouse symbols to entrez ids
library(org.Mm.eg.db)
up <- lapply(up, function(x) {
    entrez <- mapIds(org.Mm.eg.db, x, "ENTREZID", "SYMBOL")
    entrez <- unlist(entrez)
    entrez <- entrez[!is.na(entrez)]
    return(entrez)
})

down <- lapply(down, function(x) {
    entrez <- mapIds(org.Mm.eg.db, x, "ENTREZID", "SYMBOL")
    entrez <- unlist(entrez)
    entrez <- entrez[!is.na(entrez)]
    return(entrez)
})

kegg_up <- compareCluster(geneCluster = up, fun = enrichKEGG) 
kegg_down <- compareCluster(geneCluster = down, fun = enrichKEGG)

# result: no enrichemnt found

# Wikipathway analysis
get_wp_organisms()
wp_up <- compareCluster(geneCluster = up, fun = enrichWP, organism = "Mus musculus", qvalueCutoff = 0.05)
wp_down <- compareCluster(geneCluster = down, fun = enrichWP, organism = "Mus musculus", qvalueCutoff = 0.05)

# make dotplots
dplot_wp_up <- dotplot(wp_up, showCategory = 10, title = "Wikipathway enrichment of common upDEGs")
dplot_wp_down <- dotplot(wp_down, showCategory = 10, title = "Wikipathway enrichment of common downDEGs")

# save the two dotplots arranged side by side
pdf(sprintf("results/common_deg/%s_wp_dotplot.pdf", Sys.Date()), width = 15, height = 6)
grid.arrange(dplot_wp_up, dplot_wp_down, ncol = 2)
dev.off()

# disease enrichemnt analysis
library(DOSE)
up_do <- compareCluster(geneCluster = up, fun = enrichDO)
down_do <- compareCluster(geneCluster = down, fun = enrichDO)

## no enrichment found 
# Reactome enrichment analysis
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# BiocManager::install("ReactomePA")

library(ReactomePA)
up_reactomePA <- compareCluster(geneCluster = up, fun = enrichPathway, organism = "mouse", pvalueCutoff = 0.05)
down_reactomePA <- compareCluster(geneCluster = down, fun = enrichPathway, organism = "mouse", pvalueCutoff = 0.05)

# make dotplots
dplot_reactomePA_up <- dotplot(up_reactomePA, showCategory = 10, title = "Reactome enrichment of common upDEGs")
dplot_reactomePA_down <- dotplot(down_reactomePA, showCategory = 10, title = "Reactome enrichment of common downDEGs")

# save the two dotplots arranged side by side
pdf(sprintf("results/common_deg/%s_reactomePA_dotplot.pdf", Sys.Date()), width = 15, height = 6)
grid.arrange(dplot_reactomePA_up, dplot_reactomePA_down, ncol = 2)
dev.off()

# MsigDB C3: regulatory target gene sets enrichment analysis
# install.packages("msigdbr")
library(msigdbr)

c3 <- msigdbr(species = "Mus musculus", category = "C3") %>% dplyr::select(gs_name, entrez_gene)
head(c3)

c3_up <- compareCluster(geneCluster = up, fun = enricher, TERM2GENE = c3)
c3_down <- compareCluster(geneCluster = down, fun = enricher, TERM2GENE = c3)

# make dotplots
dplot_c3_up <- dotplot(c3_up, showCategory = 10, title = "MsigDB C3 enrichment of common upDEGs")
dplot_c3_down <- dotplot(c3_down, showCategory = 10, title = "MsigDB C3 enrichment of common downDEGs")

# save the two dotplots arranged side by side
pdf(sprintf("results/common_deg/%s_c3_dotplot.pdf", Sys.Date()), width = 15, height = 10)
grid.arrange(dplot_c3_up, dplot_c3_down, ncol = 2)
dev.off()

# C3 enrichment analysis using loosened cutoff of Zhou DEGs

up2 <- list(kaster_zhou_up_2, huang_zhou_up_2)
down2 <- list(kaster_zhou_down_2, huang_zhou_down_2)

# name the lists
names(up2) <- c("Kaster et al. 2023", "Huang et al. 2016")
names(down2) <- c("Kaster et al. 2023", "Huang et al. 2016")

# convert symbols to entrezids
up2 <- lapply(up2, function(x) {
    entrez <- mapIds(org.Mm.eg.db, x, "ENTREZID", "SYMBOL")
    entrez <- unlist(entrez)
    entrez <- entrez[!is.na(entrez)]
    return(entrez)
})

down2 <- lapply(down2, function(x) {
    entrez <- mapIds(org.Mm.eg.db, x, "ENTREZID", "SYMBOL")
    entrez <- unlist(entrez)
    entrez <- entrez[!is.na(entrez)]
    return(entrez)
})

c3_up2 <- compareCluster(geneCluster = up2, fun = enricher, TERM2GENE = c3)
c3_down2 <- compareCluster(geneCluster = down2, fun = enricher, TERM2GENE = c3)

# make dotplots
dplot_c3_up2 <- dotplot(c3_up2, showCategory = 10, title = "MsigDB C3 enrichment of common upDEGs with Zhou DEG logFC > 0")
dplot_bp_down2 <- dotplot(c3_down2, showCategory = 10, title = "MsigDB C3 enrichment of common downDEGs with Zhou DEG logFC > 0")

# save the two dotplots arranged side by side
pdf(sprintf("results/common_deg/%s_c3_dotplot_2.pdf", Sys.Date()), width = 15, height = 10)
grid.arrange(dplot_c3_up2, dplot_bp_down2, ncol = 2)
dev.off()
