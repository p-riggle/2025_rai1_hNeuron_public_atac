# Load libraries
library(tidyverse)

# -----------------------------
# 1. Read data
# -----------------------------
df <- read.csv("/Users/satabdi/Downloads/DEseq_rlogNorm copy.csv",   check.names = FALSE)
colnames(expr)[1] <- "gene"
expr$gene <- toupper(expr$gene)

neuronal_maturation_genes <- c(
  # NaV channels
  "SCN1A", "SCN2A", "SCN4B",
  
  # Kv channels
  "KCNA2", "KCNK1",
  
  # CaV channels
  "CACNA1B", "CACNA1D", "CACNG3",
  
  # Transporters
  "SLC12A2", "SLC12A5",
  
  # ATPases
  "ATP2B2", "ATP2B4",
  
  # CaMK / synaptic signaling
  "CAMK2B",
  
  # Presynaptic / SNARE
  "SYN1", "SV2B", "STX1A", "SNAP25", "VAMP2",
  
  # Synaptic organization
  "NRXN2", "SLC17A7", "PCLO",
  
  # Post-synaptic
  "DLG4", "NLGN1", "SHANK1", "SHANK3",
  
  # Neurotransmitter receptors
  "GRIA1", "GRIK4", "GRM1", "GRIN1", "GRIN2A", "GRIN2B",
  "GABRA1", "HTR2A"
)
# -----------------------------
# 2. Gene list
# -----------------------------


library(tidyverse)


# -----------------------------
# 2. Subset maturation genes
# -----------------------------
df_sub <- df %>%
  filter(gene %in% neuronal_maturation_genes)

# -----------------------------
# 3. Z-score per gene
# -----------------------------
df_z <- df_sub %>%
  column_to_rownames("gene") %>%
  t() %>%
  scale() %>%
  t() %>%
  as.data.frame()

# -----------------------------
# 4. Convert to long format
# -----------------------------
df_long <- df_z %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "zscore")

# -----------------------------
# 5. Extract genotype + cell type
# -----------------------------
df_long <- df_long %>%
  mutate(
    genotype = case_when(
      str_detect(sample, "-wt-") ~ "WT",
      str_detect(sample, "-ko-") ~ "KO",
      str_detect(sample, "-het-") ~ "HET"
    ),
    cell_type = case_when(
      str_detect(sample, "_1$") ~ "Neuron",
      str_detect(sample, "_5$") ~ "gRG"
    )
  )

df_long$genotype <- factor(df_long$genotype, levels = c("WT", "HET", "KO"))
# -----------------------------
# 6. Average per sample
# -----------------------------
df_avg <- df_long %>%
  group_by(sample, genotype, cell_type) %>%
  summarise(mean_z = mean(zscore), .groups = "drop")

# -----------------------------
# 7. Plot
# -----------------------------
p <- ggplot(df_avg, aes(x = genotype, y = mean_z, fill = genotype)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2) +
  facet_wrap(~cell_type) +
  theme_classic() +
  labs(
    title = "Neuronal Maturation Signature (Z-score)",
    y = "Mean Z-score",
    x = "Genotype"
  )

ggsave(
  filename = "/Users/satabdi/Downloads/neuronal_maturation_zscore_boxplot.pdf",
  plot = p,   
  width = 6,
  height = 4
)




