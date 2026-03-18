# ==========================================================
# 1️⃣ Load libraries
# ==========================================================

library(dplyr)
library(tidyr)
library(ggplot2)

# ==========================================================
# 2️⃣ Gene list
# ==========================================================

genes_of_interest <- c(
  # Alpha cluster
  "PCDHA1","PCDHA2","PCDHA3","PCDHA4","PCDHA5","PCDHA6","PCDHA7","PCDHA8","PCDHA9","PCDHA10",
  "PCDHA11","PCDHA12","PCDHA13","PCDHA14","PCDHAC1","PCDHAC2",
  
  # Beta cluster
  "PCDHB1","PCDHB2","PCDHB3","PCDHB4","PCDHB5","PCDHB6","PCDHB7","PCDHB8","PCDHB9","PCDHB10",
  "PCDHB11","PCDHB12","PCDHB13","PCDHB14","PCDHB15","PCDHB16","PCDHB17","PCDHB18","PCDHB19",
  
  # Gamma A
  "PCDHGA1","PCDHGA2","PCDHGA3","PCDHGA4","PCDHGA5","PCDHGA6","PCDHGA7","PCDHGA8","PCDHGA9","PCDHGA10",
  "PCDHGA11","PCDHGA12",
  
  # Gamma B
  "PCDHGB1","PCDHGB2","PCDHGB3","PCDHGB4","PCDHGB5","PCDHGB6","PCDHGB7",
  
  # Gamma C
  "PCDHGC3","PCDHGC4","PCDHGC5",
  
  # Non-clustered
  "PCDH1","PCDH7","PCDH8","PCDH9","PCDH10","PCDH11X","PCDH11Y",
  "PCDH12","PCDH15","PCDH17","PCDH18","PCDH19","PCDH20"
)
# ==========================================================
# 3️⃣ Load datasets
# ==========================================================

astro_het  <- read.csv("/Users/satabdi/Downloads/rstbl_Astro_wt.vs.het.csv")
astro_ko   <- read.csv("/Users/satabdi/Downloads/rstbl_Astro_wt.vs.ko.csv")

ipc_het    <- read.csv("/Users/satabdi/Downloads/rstbl_IPC_wt.vs.het.csv")
ipc_ko     <- read.csv("/Users/satabdi/Downloads/rstbl_IPC_wt.vs.ko.csv")

mesko_het  <- read.csv("/Users/satabdi/Downloads/rstbl_MesKO_wt.vs.het.csv")
mesko_ko   <- read.csv("/Users/satabdi/Downloads/rstbl_MesKO_wt.vs.ko.csv")

meswt_het  <- read.csv("/Users/satabdi/Downloads/rstbl_MesWT_wt.vs.het.csv")
meswt_ko   <- read.csv("/Users/satabdi/Downloads/rstbl_MesWT_wt.vs.ko.csv")

neuron_het <- read.csv("/Users/satabdi/Downloads/rstbl_Neuron_wt.vs.het.csv")
neuron_ko  <- read.csv("/Users/satabdi/Downloads/rstbl_Neuron_wt.vs.ko.csv")

npcg1_het  <- read.csv("/Users/satabdi/Downloads/rstbl_NPCG1_wt.vs.het.csv")
npcg1_ko   <- read.csv("/Users/satabdi/Downloads/rstbl_NPCG1_wt.vs.ko.csv")

npcg2m_het <- read.csv("/Users/satabdi/Downloads/rstbl_NPCG2M_wt.vs.het.csv")
npcg2m_ko  <- read.csv("/Users/satabdi/Downloads/rstbl_NPCG2M_wt.vs.ko.csv")

# ==========================================================
# 4️⃣ Function to extract genes
# ==========================================================

extract_data <- function(df, celltype, comparison){
  
  if(!"gene" %in% colnames(df)){
    colnames(df)[1] <- "gene"
  }
  
  df %>%
    filter(gene %in% genes_of_interest) %>%
    select(gene, log2FoldChange, padj) %>%
    mutate(CellType = celltype,
           Comparison = comparison)
}

# ==========================================================
# 5️⃣ Combine all datasets
# ==========================================================

combined <- bind_rows(
  
  extract_data(astro_het,"Astro","WT_vs_HET"),
  extract_data(astro_ko,"Astro","WT_vs_KO"),
  
  extract_data(ipc_het,"IPC","WT_vs_HET"),
  extract_data(ipc_ko,"IPC","WT_vs_KO"),
  
  extract_data(meswt_het,"MesWT","WT_vs_HET"),
  extract_data(meswt_ko,"MesWT","WT_vs_KO"),
  
  extract_data(mesko_het,"MesKO","WT_vs_HET"),
  extract_data(mesko_ko,"MesKO","WT_vs_KO"),
  
  extract_data(neuron_het,"Neuron","WT_vs_HET"),
  extract_data(neuron_ko,"Neuron","WT_vs_KO"),
  
  extract_data(npcg1_het,"NPCG1","WT_vs_HET"),
  extract_data(npcg1_ko,"NPCG1","WT_vs_KO"),
  
  extract_data(npcg2m_het,"NPCG2M","WT_vs_HET"),
  extract_data(npcg2m_ko,"NPCG2M","WT_vs_KO")
)

combined <- combined %>%
  mutate(Group = paste(CellType, Comparison, sep = "_"))

# ==========================================================
# 6️⃣ Fix padj and create significance score
# ==========================================================

dot_data <- combined %>%
  mutate(
    padj = ifelse(is.na(padj),1,padj),
    padj = ifelse(padj==0,1e-300,padj),
    neglog10_padj = -log10(padj)
  )

# ==========================================================
# 7️⃣ Order group
# ==========================================================

desired_order <- c(
  "NPCG2M_WT_vs_HET","NPCG2M_WT_vs_KO",
  "NPCG1_WT_vs_HET","NPCG1_WT_vs_KO",
  "MesWT_WT_vs_HET","MesWT_WT_vs_KO",
  "MesKO_WT_vs_HET","MesKO_WT_vs_KO",
  "IPC_WT_vs_HET","IPC_WT_vs_KO",
  "Neuron_WT_vs_HET","Neuron_WT_vs_KO",
  "Astro_WT_vs_HET","Astro_WT_vs_KO"
)

dot_data$Group <- factor(dot_data$Group, levels = desired_order)
dot_data$gene  <- factor(dot_data$gene, levels = rev(genes_of_interest))

# ==========================================================
# 8️⃣ Generate Dot Plot
# ==========================================================

dot_plot <- ggplot(dot_data,
                   aes(x = Group,
                       y = gene)) +
  
  geom_point(aes(fill = log2FoldChange,
                 size = neglog10_padj),
             shape = 21,
             color = "black",
             stroke = 0.5) +
  
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    name = "log2FC"
  ) +
  
  scale_size(
    range = c(1.5,8),
    name = "-log10(padj)"
  ) +
  
  theme_minimal(base_size = 12) +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7),
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_blank()
  ) +
  
  labs(
    x = "Cell Type ",
    y = "Gene",
    title = "Protocadherin Gene Expression Changes"
  )

dot_plot
# ==========================================================
# 9️⃣ Save plot
# ==========================================================

ggsave("/Users/satabdi/Downloads/PCDH_gene_dotplot.pdf",
       dot_plot,
       width = 20,
       height = 12)