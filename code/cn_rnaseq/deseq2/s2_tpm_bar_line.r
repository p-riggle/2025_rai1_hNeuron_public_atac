# Q: Are differentiation marker genes expressed in an expacted ways?

tpm <- read.delim("data/cn_rnaseq/blueberry_output/combine_counts/gene_TPM.txt", header = TRUE, sep = "\t")
head(tpm)
colnames(tpm) <- gsub("Rai1.|.TPM|d", "", colnames(tpm))
colnames(tpm)[1] <- "gene_id"
# remove het.d34.rep1 and het.d60.rep4
tpm[, c("het_34_rep1", "het_60_rep4")] <- NULL
colnames(tpm)

# Load necessary libraries
# install.packages("tidyverse")
library(tidyverse)

# # Plot TPM values for RAI1 gene

library(dplyr)  # for data manipulation
library(tidyr)  # for data reshaping
library(ggplot2)  # for plotting

library(tidyverse)

# make function for line plotting

plot_gene_expression <- function(gene_name) {
  
  # extract columns for specified gene
  gene_data <- tpm %>%
    filter(gene_id == gene_name) %>%
    select(starts_with("wt_"), starts_with("het_"), starts_with("ko_"))
  
  # reshape data into long format for ggplot2
  gene_data_long <- gene_data %>%
    pivot_longer(cols = everything(),
                 names_to = c("genotype", "day", "replicate"),
                 names_sep = "_",
                 values_to = "tpm")
  
  # Convert the genotype column to factor with the specified order
  genotype_order <- c("wt", "het", "ko")
  gene_data_long$genotype <- factor(gene_data_long$genotype, levels = genotype_order)
  
  # Define color scheme with darker, more visible colors
  library(colorspace)
  my_colors <- c("#0c17bbff", "#722facff", "#d354f6ff")  # Dark purple (WT), medium purple (HET), darker pink/red (KO) 
  
  # calculate mean and standard error for each genotype and day
  gene_summary <- gene_data_long %>%
    group_by(genotype, day) %>%
    summarize(mean_tpm = mean(tpm),
              sem_tpm = sd(tpm) / sqrt(length(tpm))) %>% 
    ungroup() %>%
    mutate(upper = mean_tpm + sem_tpm, lower = mean_tpm - sem_tpm)
  
  # plot mean tpm with error bars and individual points
  ggplot(gene_summary, aes(x = day, y = mean_tpm, color = genotype, group = genotype)) +
    geom_line(size = 1.5) +
    geom_point(size = 2, alpha = 0.8) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, size = 1) +
    labs(x = "Day", y = "TPM", title = gene_name) +
    scale_color_manual(values = my_colors) +
    theme_classic()
}


plot_gene_expression_nolegend <- function(gene_name) {
  
  # extract columns for specified gene
  gene_data <- tpm %>%
    filter(gene_id == gene_name) %>%
    select(starts_with("wt_"), starts_with("het_"), starts_with("ko_"))
  
  # reshape data into long format for ggplot2
  gene_data_long <- gene_data %>%
    pivot_longer(cols = everything(),
                 names_to = c("genotype", "day", "replicate"),
                 names_sep = "_",
                 values_to = "tpm")
  
  # Convert the genotype column to factor with the specified order
  genotype_order <- c("wt", "het", "ko")
  gene_data_long$genotype <- factor(gene_data_long$genotype, levels = genotype_order)
  
  # Define color scheme with darker, more visible colors
  my_colors <- c("#2E004D", "#ce69e3ff", "#1788e5ff")  # Dark purple (WT), medium purple (HET), darker pink/red (KO)
  
  # calculate mean and standard error for each genotype and day
  gene_summary <- gene_data_long %>%
    group_by(genotype, day) %>%
    summarize(mean_tpm = mean(tpm),
              sem_tpm = sd(tpm) / sqrt(length(tpm))) %>% 
    ungroup() %>%
    mutate(upper = mean_tpm + sem_tpm, lower = mean_tpm - sem_tpm)
  
  # plot mean tpm with error bars and individual points
  ggplot(gene_summary, aes(x = day, y = mean_tpm, color = genotype, group = genotype)) +
    geom_line(size = 1.5) +
    geom_point(size = 2, alpha = 0.8) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, size = 1) +
    labs(x = "Day", y = "TPM", title = gene_name) +
    scale_color_manual(values = my_colors) +
    theme_classic() +
    guides(color = FALSE)
}


# make lineplots
# install.packages("ggpubr")
library("ggpubr")
library("httpgd")

# rai1
rai1 <- plot_gene_expression("RAI1") 
rai1
pdf("results/cn_rnaseq/deseq2/lp_rai1.pdf", width = 4.6, height = 4)
print(rai1)
dev.off()

rai1

# Stem cell marker
pou5f1 <- plot_gene_expression("POU5F1")
nanog <- plot_gene_expression("NANOG")

lp_stemcell <- ggarrange(pou5f1, nanog, labels = c("A", "B"), ncol = 2, nrow = 1)
lp_stemcell 
pdf(sprintf("results/deseq2/lp_stemcells.pdf"), width = 9, height = 4)
print(lp_stemcell)
dev.off()

# Progenitor cell markers
pax6 <- plot_gene_expression_nolegend("PAX6")
otx1 <- plot_gene_expression_nolegend("OTX1")
emx1 <- plot_gene_expression_nolegend("EMX1") 
foxg1 <- plot_gene_expression_nolegend("FOXG1") #also expressed in post-mitotic neurons in medial dorsal latereal pallia

# Deeper layer markers
eomes <- plot_gene_expression_nolegend("EOMES") #TBR2, 
tbr1 <- plot_gene_expression_nolegend("TBR1") 
bcl11b <- plot_gene_expression_nolegend("BCL11B") #CTIP2, layer IV. -| TBR1 at layer 5

# Upper layer markers
satb2 <- plot_gene_expression_nolegend("SATB2")
pou3f2 <- plot_gene_expression_nolegend("POU3F2") #BRN2
cux1 <- plot_gene_expression("CUX1") #BRN2

# arrage the neuronal markers into one figure
library(gridExtra)

pdf(sprintf("results/deseq2/lp_neuron.pdf"), width = 12, height = 9)
grid.arrange(
  arrangeGrob(pax6, otx1, emx1, foxg1, ncol = 4, top = "Neural Progenitor"), 
  arrangeGrob(eomes, tbr1, bcl11b, ncol = 4, top = "Deeper Layer Neurons"),
  arrangeGrob(satb2, pou3f2, cux1, ncol = 4, top = "Upper Layer Neurons"), 
  ncol = 1, heights = c(1, 1, 1)
)
dev.off()

# Astrocyte markers
gfap <- plot_gene_expression_nolegend("GFAP")
aldh1l1 <- plot_gene_expression_nolegend("ALDH1L1")
s100b <- plot_gene_expression_nolegend("S100B")
sox9 <- plot_gene_expression_nolegend("SOX9")

# Oligodendrocyte markers
olig2 <- plot_gene_expression_nolegend("OLIG2")
mbp <- plot_gene_expression_nolegend("MBP")
sox10 <- plot_gene_expression("SOX10")

# Microglia markers
sall1 <- plot_gene_expression_nolegend("SALL1")
hexb <- plot_gene_expression_nolegend("HEXB")
tmem119 <- plot_gene_expression_nolegend("TMEM119")
gpr34 <- plot_gene_expression_nolegend("GPR34")

# arrange the plots using grid.arrange
library(gridExtra)
library(ggplot2)
pdf(sprintf("results/deseq2/lp_glia.pdf"), width = 12, height = 9)
grid.arrange(
  arrangeGrob(gfap, aldh1l1, s100b, sox9, ncol = 4, top = "Astrocyte"), 
  arrangeGrob(olig2, mbp, sox10, ncol = 4, top = "Oligodendrocyte"),
  arrangeGrob(sall1, hexb, tmem119, gpr34, ncol = 4, top = "Microglia"), 
  ncol = 1, heights = c(1, 1, 1)
)
dev.off()

# circadian genes
clock <- plot_gene_expression_nolegend("CLOCK")
bmal1 <- plot_gene_expression("ARNTL")
per1 <- plot_gene_expression_nolegend("PER1")
per2 <- plot_gene_expression_nolegend("PER2")
cry1 <- plot_gene_expression_nolegend("CRY1")
cry2 <- plot_gene_expression_nolegend("CRY2")
npas2 <- plot_gene_expression_nolegend("NPAS2")
reverba <- plot_gene_expression_nolegend("NR1D1")
reverbb <- plot_gene_expression_nolegend("NR1D2")

pdf(sprintf("results/deseq2/lp_clockgenes.pdf"), width = 12, height = 9)
grid.arrange(
  arrangeGrob(clock, bmal1, ncol = 4, top = "Core Clock Genes"), 
  arrangeGrob(per1, per2, cry1, cry2, ncol = 4),
  arrangeGrob(npas2, reverba, reverbb, ncol = 4),  
  heights = c(1, 1, 1)
)
dev.off()

# cilia markers
arl13b <- plot_gene_expression("ARL13B")
arl13b

adcy3 <- plot_gene_expression("ADCY3")
adcy3

# make combined plot for cilia marker

pdf(sprintf("results/deseq2/lp_ciliagenes.pdf"), width = 8, height = 4)
grid.arrange(
  arrangeGrob(arl13b, adcy3, ncol = 2, top = "Cilia markers")
)
dev.off()

# BDNF
bdnf <- plot_gene_expression("BDNF")
bdnf

pdf(sprintf("results/cn_rnaseq/deseq2/lp_bdnf.pdf"), width = 4.6, height = 4)
print(bdnf)
dev.off()

# for KDM5C KARE grant
kdm5c <- plot_gene_expression_nolegend("KDM5C")
kdm5c
pdf(sprintf("results/deseq2/lp_kdm5c.pdf"), width = 12, height = 3)
grid.arrange(
  arrangeGrob(kdm5c, nanog, pax6, ncol=3, top = "KDM5C expression in differentiating hESCs"),
  arrangeGrob(tbr1, satb2, ncol=3, top = "NANOG"),
  heights = c(1, 1)
)
dev.off()

# plot kdm5d expression in differentiating hESCs
kdm5d <- plot_gene_expression_nolegend("KDM5D")
kdm5d

# change the color of the line plot
kdm5d + scale_colour_manual(values = c("black", "red", "blue", "green", "orange", "purple", "brown", "pink", "grey"))









