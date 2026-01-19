# Q: Are differentiation marker genes expressed in an expected ways?

# Q: Are differentiation marker genes expressed in an expected ways?

# Read the data
tpm <- read.delim("/Users/satabdi/Downloads/rscripts/ineurons_gene_expected_count.annot (1).txt", header = TRUE, sep = "\t")
head(tpm)

# Clean column names
colnames(tpm) <- gsub("Rai1.|.TPM|d", "", colnames(tpm))

# Reorder the Column Names
colnames(tpm)[3] <- "Gene_id"
colnames(tpm)[1] <- "gene_id"
colnames(tpm)[2] <- "entrezgene_id"
colnames(tpm)[4] <- "description"
tpm <- tpm[c(1:4,5,16,27,35,9:12,22:25,36:39,13:17,26:30,40,6:8,18:21,31:34)]
colnames(tpm)



#Remove extra or duplicated Columns
tpm <- tpm[, -c(24, 27)]
colnames(tpm)

# Renaming the column Names
new_colnames <- c(
  "wt_d0_s1", "wt_d0_s2", "wt_d0_s3", "wt_d0_s4", "wt_d3_s1", "wt_d3_s2", 
  "wt_d3_s3", "wt_d3_s4", "wt_d7_s1", "wt_d7_s2", "wt_d7_s3", "wt_d7_s4", 
  "het_d0_s1", "het_d0_s2", "het_d0_s3", "het_d0_s4", "het_d3_s1", "het_d3_s2", 
  "het_d3_s3", "het_d3_s4" , "het_d7_s1", "het_d7_s2", "het_d7_s3", "het_d7_s4",
  "ko_d0_s1", "ko_d0_s2", "ko_d0_s3", "ko_d0_s4", "ko_d3_s1", "ko_d3_s2", 
  "ko_d3_s3", "ko_d3_s4" , "ko_d7_s1", "ko_d7_s2", "ko_d7_s3", "ko_d7_s4"
)


# Rename the columns
colnames(tpm)[5:40] <- new_colnames

# Load necessary libraries
#install.packages("tidyverse")
library(tidyverse)

# Make function for line plotting
plot_gene_expression <- function(gene_name) {
  
  # Extract columns for specified gene
  gene_data <- tpm %>%
    dplyr::filter(Gene_id == gene_name) %>%
    dplyr::select(starts_with("wt_"), starts_with("het_"), starts_with("ko_"))
  
  # Reshape data into long format for ggplot2
  gene_data_long <- gene_data %>%
    pivot_longer(cols = everything(),
                 names_to = c("genotype", "day", "replicate"),
                 names_sep = "_",
                 values_to = "tpm")
  
  # Convert the genotype column to factor with the specified order
  genotype_order <- c("wt", "het", "ko")
  gene_data_long$genotype <- factor(gene_data_long$genotype, levels = genotype_order)
  
  # Define color scheme
  library(colorspace)
  my_colors <- sequential_hcl(3, palette = "TealGrn") 
  
  # Calculate mean and standard error for each genotype and day
  gene_summary <- gene_data_long %>%
    group_by(genotype, day) %>%
    summarize(mean_tpm = mean(tpm),
              sem_tpm = sd(tpm) / sqrt(length(tpm)),
              .groups = "drop") %>% 
    ungroup() %>%
    mutate(upper = mean_tpm + sem_tpm, lower = mean_tpm - sem_tpm)
  
  # Plot mean tpm with error bars and individual points
  ggplot(gene_summary, aes(x = day, y = mean_tpm, color = genotype, group = genotype)) +
    geom_line(linewidth = 1.5) +  # Use linewidth instead of size
    geom_point(size = 2, alpha = 0.8) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, size = 0.5) +
    labs(x = "Day", y = "TPM", title = gene_name) +
    scale_color_manual(values = my_colors) +
    theme_classic() 
}

plot_gene_expression_nolegend <- function(gene_name) {
  
  # Extract columns for specified gene
  gene_data <- tpm %>%
    dplyr::filter(Gene_id == gene_name) %>%
    dplyr::select(starts_with("wt_"), starts_with("het_"), starts_with("ko_"))
  
  # Reshape data into long format for ggplot2
  gene_data_long <- gene_data %>%
    pivot_longer(cols = everything(),
                 names_to = c("genotype", "day", "replicate"),
                 names_sep = "_",
                 values_to = "tpm")
  
  # Convert the genotype column to factor with the specified order
  genotype_order <- c("wt", "het", "ko")
  gene_data_long$genotype <- factor(gene_data_long$genotype, levels = genotype_order)
  
  # Define color scheme
  library(colorspace)
  my_colors <- sequential_hcl(3, palette = "TealGrn") 
  
  # Calculate mean and standard error for each genotype and day
  gene_summary <- gene_data_long %>%
    group_by(genotype, day) %>%
    summarize(mean_tpm = mean(tpm),
              sem_tpm = sd(tpm) / sqrt(length(tpm)),
              .groups = "drop") %>% 
    ungroup() %>%
    mutate(upper = mean_tpm + sem_tpm, lower = mean_tpm - sem_tpm)
  
  # Plot mean tpm with error bars and individual points
  ggplot(gene_summary, aes(x = day, y = mean_tpm, color = genotype, group = genotype)) +
    geom_line(linewidth = 1.5) +  # Use linewidth instead of size
    geom_point(size = 2, alpha = 0.8) +
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, size = 0.5) +
    labs(x = "Day", y = "TPM", title = gene_name) +
    scale_color_manual(values = my_colors) +
    theme_classic() + 
    guides(color = "none")
}



# Make lineplots
#install.packages("ggpubr")
library("ggpubr")

#install.packages("colorspace")
library(colorspace)


# RAI1
rai1 <- plot_gene_expression("RAI1") 
rai1
pdf("/Users/satabdimohanty/Desktop/U-M/Iwase Lab/Project/ineurons_RNA_seq/Results/lineplot_rai1.pdf", width = 4.6, height = 4)
print(rai1)
dev.off()

#BDNF
bdnf <- plot_gene_expression("BDNF") 
bdnf
pdf("/Users/satabdimohanty/Desktop/U-M/Iwase Lab/Project/ineurons_RNA_seq/Results/lineplot_bdnf.pdf", width = 4.6, height = 4)
print(bdnf)
dev.off()

#TH
th <- plot_gene_expression_nolegend("TH") 
th
pdf("/Users/satabdimohanty/Desktop/U-M/Iwase Lab/Project/ineurons_RNA_seq/Results/lineplot_TH.pdf", width = 4.6, height = 4)
print(th)
dev.off()

# Stem cell marker
pou5f1 <- plot_gene_expression("POU5F1")
nanog <- plot_gene_expression("NANOG")

lp_stemcell <- ggarrange(pou5f1, nanog, lable = c("A", "B"), ncol = 2, nrow = 1)
lp_stemcell 
pdf(sprintf("r/Users/satabdimohanty/Desktop/U-M/Iwase Lab/Project/ineurons_RNA_seq/Results/lineplot_stemcells.pdf", Sys.Date()), width = 9, height = 4)
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
cux1 <- plot_gene_expression_nolegend("CUX1") #BRN2

# arrage the neuronal markers into one figure
library(gridExtra)

pdf(sprintf("/Users/satabdimohanty/Desktop/U-M/Iwase Lab/Project/ineurons_RNA_seq/Results/%s_lp_neuron.pdf", Sys.Date()), width = 12, height = 9)
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
sox10 <- plot_gene_expression_nolegend("SOX10")

# Microglia markers
sall1 <- plot_gene_expression_nolegend("SALL1")
hexb <- plot_gene_expression_nolegend("HEXB")
tmem119 <- plot_gene_expression_nolegend("TMEM119")
gpr34 <- plot_gene_expression_nolegend("GPR34")

# arrange the plots using grid.arrange
library(gridExtra)
library(ggplot2)
pdf(sprintf("/Users/satabdimohanty/Desktop/U-M/Iwase Lab/Project/ineurons_RNA_seq/Results/%s_lp_glia.pdf", Sys.Date()), width = 12, height = 9)
grid.arrange(
  arrangeGrob(gfap, aldh1l1, s100b, sox9, ncol = 4, top = "Astrocyte"), 
  arrangeGrob(olig2, mbp, sox10, ncol = 4, top = "Oligodendrocyte"),
  arrangeGrob(sall1, hexb, tmem119, gpr34, ncol = 4, top = "Microglia"), 
  ncol = 1, heights = c(1, 1, 1)
)
dev.off()

# circadian genes
clock <- plot_gene_expression_nolegend("CLOCK")
bmal1 <- plot_gene_expression_nolegend("BMAL1")
per1 <- plot_gene_expression_nolegend("PER1")
per2 <- plot_gene_expression_nolegend("PER2")
cry1 <- plot_gene_expression_nolegend("CRY1")
cry2 <- plot_gene_expression_nolegend("CRY2")
npas2 <- plot_gene_expression_nolegend("NPAS2")
reverba <- plot_gene_expression_nolegend("NR1D1")
reverbb <- plot_gene_expression_nolegend("NR1D2")

pdf(sprintf("/Users/satabdimohanty/Desktop/U-M/Iwase Lab/Project/ineurons_RNA_seq/Results/%s_lp_clockgenes.pdf", Sys.Date()), width = 12, height = 9)
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

pdf(sprintf("/Users/satabdimohanty/Desktop/U-M/Iwase Lab/Project/ineurons_RNA_seq/Results/%s_lp_ciliagenes.pdf", Sys.Date()), width = 8, height = 4)
grid.arrange(
  arrangeGrob(arl13b, adcy3, ncol = 2, top = "Cilia markers")
)
dev.off()



