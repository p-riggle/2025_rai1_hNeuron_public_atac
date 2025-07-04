
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
