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

# make function for line plotting
lineplot.het <- function(gene_name) {
  
  # extract columns for specified gene
  gene_data <- tpm %>%
    filter(gene_id == gene_name) %>%
    select(starts_with("wt_"), starts_with("het_"))
  
  # reshape data into long format for ggplot2
  gene_data_long <- gene_data %>%
    pivot_longer(cols = everything(),
                 names_to = c("genotype", "day", "replicate"),
                 names_sep = "_",
                 values_to = "tpm")
  
  # Convert the genotype column to factor with the specified order
  genotype_order <- c("wt", "het")
  gene_data_long$genotype <- factor(gene_data_long$genotype, levels = genotype_order)
  
  # Define color scheme
  library(colorspace)
  color_het <- c("#757582", "#5876D1") 
  
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
    scale_color_manual(values = color_het) +
    theme_classic()
}

rai1_het <- lineplot.het("RAI1")
rai1_het

# Stem cell marker
library("ggpubr")
pou5f1_het <- lineplot.het("POU5F1")
nanog_het <- lineplot.het("NANOG")

# arrange the stem cell markers into one figure
lp_stemcell_het <- ggarrange(pou5f1_het, nanog_het, ncol = 2, nrow = 1)
lp_stemcell_het
pdf(sprintf("results/deseq2/sep/%s_het_lp_stemcell.pdf", Sys.Date()), width = 9, height = 4)
print(lp_stemcell_het)
dev.off()

# Progenitor cell markers
pax6_het <- lineplot.het("PAX6")
otx1_het <- lineplot.het("OTX1")
emx1_het <- lineplot.het("EMX1") 
foxg1_het <- lineplot.het("FOXG1") #also expressed in post-mitotic neurons in medial dorsal latereal pallia

# Deeper layer markers
eomes_het <- lineplot.het("EOMES") #TBR2, 
tbr1_het <- lineplot.het("TBR1") 
bcl11b_het <- lineplot.het("BCL11B") #CTIP2, layer IV. -| TBR1 at layer 5

# Upper layer markers
satb2_het <- lineplot.het("SATB2")
pou3f2_het <- lineplot.het("POU3F2") #BRN2
cux1_het <- lineplot.het("CUX1") #BRN2

# arrage the neuronal markers into one figure
library(gridExtra)

pdf(sprintf("results/deseq2/sep/%s_het_lp_neuron.pdf", Sys.Date()), width = 14, height = 9)
grid.arrange(
  arrangeGrob(pax6_het, otx1_het, emx1_het, foxg1_het, ncol = 4, top = "Neural Progenitor"), 
  arrangeGrob(eomes_het, tbr1_het, bcl11b_het, ncol = 4, top = "Deeper Layer Neurons"),
  arrangeGrob(satb2_het, pou3f2_het, cux1_het, ncol = 4, top = "Upper Layer Neurons"), 
  ncol = 1, heights = c(1, 1, 1)
)
dev.off()

# Astrocyte markers
gfap_het <- lineplot.het("GFAP")
aldh1l1_het <- lineplot.het("ALDH1L1")
s100b_het <- lineplot.het("S100B")
sox9_het <- lineplot.het("SOX9")

# Oligodendrocyte markers
olig2_het<- lineplot.het("OLIG2")
mbp_het <- lineplot.het("MBP")
sox10_het <- lineplot.het("SOX10")

# Microglia markers
sall1_het <- lineplot.het("SALL1")
hexb_het <- lineplot.het("HEXB")
tmem119_het <- lineplot.het("TMEM119")
gpr34_het <- lineplot.het("GPR34")

# arrange the glia plots using grid.arrange
pdf(sprintf("results/deseq2/sep/%s_het_lp_glia.pdf", Sys.Date()), width = 14, height = 9)
grid.arrange(
  arrangeGrob(gfap_het, aldh1l1_het, s100b_het, sox9_het, ncol = 4, top = "Astrocyte"), 
  arrangeGrob(olig2_het, mbp_het, sox10_het, ncol = 4, top = "Oligodendrocyte"),
  arrangeGrob(sall1_het, hexb_het, tmem119_het, gpr34_het, ncol = 4, top = "Microglia"), 
  ncol = 1, heights = c(1, 1, 1)
)
dev.off()


# circadian genes
clock_het <- lineplot.het("CLOCK")
bmal1_het <- lineplot.het("ARNTL")
per1_het <- lineplot.het("PER1")
per2_het <- lineplot.het("PER2")
cry1_het <- lineplot.het("CRY1")
cry2_het <- lineplot.het("CRY2")
npas2_het <- lineplot.het("NPAS2")
reverba_het <- lineplot.het("NR1D1")
reverbb_het <- lineplot.het("NR1D2")

# arrange the circadian genes into one figure
pdf(sprintf("results/deseq2/sep/%s_het_lp_circadian.pdf", Sys.Date()), width = 14, height = 9)
grid.arrange(
  arrangeGrob(clock_het, bmal1_het, ncol = 4, top = "Core Clock Genes"), 
  arrangeGrob(per1_het, per2_het, cry1_het, cry2_het, ncol = 4),
  arrangeGrob(npas2_het, reverba_het, reverbb_het, ncol = 4),  
  heights = c(1, 1, 1)
)
dev.off()

# KO plots
# make function for line plotting
lineplot.ko <- function(gene_name) {
  # extract columns for specified gene
  gene_data <- tpm %>%
    filter(gene_id == gene_name) %>%
    select(starts_with("wt_"), starts_with("ko_"))
  
  # reshape data into long format for ggplot2
  gene_data_long <- gene_data %>%
    pivot_longer(cols = everything(),
                 names_to = c("genotype", "day", "replicate"),
                 names_sep = "_",
                 values_to = "tpm")
  
  # Convert the genotype column to factor with the specified order
  genotype_order <- c("wt", "ko")
  gene_data_long$genotype <- factor(gene_data_long$genotype, levels = genotype_order)
  
  # Define color scheme
  library(colorspace)
  color_ko <- c("#757582", "#5876D1") 
  
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
    scale_color_manual(values = color_ko) +
    theme_classic()
}

rai1_ko <- lineplot.ko("RAI1")
rai1_ko

# Stem cell marker
library("ggpubr")
pou5f1_ko <- lineplot.ko("POU5F1")
nanog_ko <- lineplot.ko("NANOG")

# arrange the stem cell markers into one figure
lp_stemcell_ko <- ggarrange(pou5f1_ko, nanog_ko, ncol = 2, nrow = 1)
lp_stemcell_ko
pdf(sprintf("results/deseq2/sep/%s_ko_lp_stemcell.pdf", Sys.Date()), width = 9, height = 4)
print(lp_stemcell_ko)
dev.off()

# Progenitor cell markers
pax6_ko <- lineplot.ko("PAX6")
otx1_ko <- lineplot.ko("OTX1")
emx1_ko <- lineplot.ko("EMX1") 
foxg1_ko <- lineplot.ko("FOXG1") #also expressed in post-mitotic neurons in medial dorsal latereal pallia

# Deeper layer markers
eomes_ko <- lineplot.ko("EOMES") #TBR2, 
tbr1_ko <- lineplot.ko("TBR1") 
bcl11b_ko <- lineplot.ko("BCL11B") #CTIP2, layer IV. -| TBR1 at layer 5

# Upper layer markers
satb2_ko <- lineplot.ko("SATB2")
pou3f2_ko <- lineplot.ko("POU3F2") #BRN2
cux1_ko <- lineplot.ko("CUX1") #BRN2

# arrage the neuronal markers into one figure
library(gridExtra)

pdf(sprintf("results/deseq2/sep/%s_ko_lp_neuron.pdf", Sys.Date()), width = 14, height = 9)
grid.arrange(
  arrangeGrob(pax6_ko, otx1_ko, emx1_ko, foxg1_ko, ncol = 4, top = "Neural Progenitor"), 
  arrangeGrob(eomes_ko, tbr1_ko, bcl11b_ko, ncol = 4, top = "Deeper Layer Neurons"),
  arrangeGrob(satb2_ko, pou3f2_ko, cux1_ko, ncol = 4, top = "Upper Layer Neurons"), 
  ncol = 1, heights = c(1, 1, 1)
)
dev.off()

# Astrocyte markers
gfap_ko <- lineplot.ko("GFAP")
aldh1l1_ko <- lineplot.ko("ALDH1L1")
s100b_ko <- lineplot.ko("S100B")
sox9_ko <- lineplot.ko("SOX9")

# Oligodendrocyte markers
olig2_ko<- lineplot.ko("OLIG2")
mbp_ko <- lineplot.ko("MBP")
sox10_ko <- lineplot.ko("SOX10")

# Microglia markers
sall1_ko <- lineplot.ko("SALL1")
hexb_ko <- lineplot.ko("HEXB")
tmem119_ko <- lineplot.ko("TMEM119")
gpr34_ko <- lineplot.ko("GPR34")

# arrange the glia plots using grid.arrange
pdf(sprintf("results/deseq2/sep/%s_ko_lp_glia.pdf", Sys.Date()), width = 14, height = 9)
grid.arrange(
  arrangeGrob(gfap_ko, aldh1l1_ko, s100b_ko, sox9_ko, ncol = 4, top = "Astrocyte"), 
  arrangeGrob(olig2_ko, mbp_ko, sox10_ko, ncol = 4, top = "Oligodendrocyte"),
  arrangeGrob(sall1_ko, hexb_ko, tmem119_ko, gpr34_ko, ncol = 4, top = "Microglia"), 
  ncol = 1, heights = c(1, 1, 1)
)
dev.off()


# circadian genes
clock_ko <- lineplot.ko("CLOCK")
bmal1_ko <- lineplot.ko("ARNTL")
per1_ko <- lineplot.ko("PER1")
per2_ko <- lineplot.ko("PER2")
cry1_ko <- lineplot.ko("CRY1")
cry2_ko <- lineplot.ko("CRY2")
npas2_ko <- lineplot.ko("NPAS2")
reverba_ko <- lineplot.ko("NR1D1")
reverbb_ko <- lineplot.ko("NR1D2")

# arrange the circadian genes into one figure
pdf(sprintf("results/deseq2/sep/%s_ko_lp_circadian.pdf", Sys.Date()), width = 14, height = 9)
grid.arrange(
  arrangeGrob(clock_ko, bmal1_ko, ncol = 4, top = "Core Clock Genes"), 
  arrangeGrob(per1_ko, per2_ko, cry1_ko, cry2_ko, ncol = 4),
  arrangeGrob(npas2_ko, reverba_ko, reverbb_ko, ncol = 4),  
  heights = c(1, 1, 1)
)
dev.off()


th_ko <- lineplot.ko("TH")
th_ko
