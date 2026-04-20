# ==================================================
# Load libraries
# ==================================================
library(tidyverse)

# ==================================================
# 1. Read data
# ==================================================
df <- read.csv("/Users/satabdi/Downloads/DEseq_rlogNorm copy.csv",
               check.names = FALSE)

colnames(df)[1] <- "gene"
df$gene <- toupper(df$gene)

neuronal_maturation_genes <- c(
  # NaV channels
  "SCN1A", "SCN2A", "SCN3A", "SCN3B", "SCN4B",
  
  # Kv channels
  "KCNA2", "KCNB1", "KCNC4", "KCNK1",
  
  # CaV channels
  "CACNA1B", "CACNA1C", "CACNA1D", "CACNA1E", "CACNG3",
  
  # Transporters
  "SLC12A2", "SLC12A5", "SLC17A7",
  
  # ATPases / ion pumps
  "ATP1A1", "ATP1A3", "ATP2B1", "ATP2B2", "ATP2B4",
  
  # CaMK / signaling
  "CAMK2A", "CAMK2B",
  
  # Presynaptic / vesicle / SNARE
  "SYN1", "SYN2", "SYP",
  "SV2A", "SV2B",
  "STX1A", "STX1B",
  "SNAP25", "SNAP91",
  "SYT1",
  "VAMP2", "VTI1B",
  "RAB3A", "RAB3C",
  
  # Synaptic organization / adhesion
  "NRXN1", "NRXN2",
  "SYNGAP1", "PCLO",
  "NLGN1",
  
  # Post-synaptic scaffolding
  "DLG2", "DLG4",
  "HOMER1",
  "SHANK1", "SHANK2", "SHANK3",
  "SHISA9",
  
  # Glutamate receptors
  "GRIA1", "GRIA2",
  "GRIK4",
  "GRM1", "GRM5", "GRM7",
  "GRIN1", "GRIN2A", "GRIN2B",
  
  # GABA receptors
  "GABRA1", "GABRA2", "GABRA3", "GABRA5",
  "GABRG2", "GABRG3",
  
  # Serotonin receptors
  "HTR2A", "HTR5A"
)

genes_use <- neuronal_maturation_genes %>%
  str_split("\n") %>%
  unlist() %>%
  str_trim() %>%
  .[. != ""] %>%
  toupper()
# -----------------------------
# 2. Gene list
# -----------------------------


library(tidyverse)


# -----------------------------
# 2. Subset maturation genes
# -----------------------------
df_sub <- df %>%
  filter(gene %in% genes_use)


# ==================================================
# 4. Long format + keep ONLY gRG (_5)
# ==================================================
df_long <- df_sub %>%
  pivot_longer(-gene, names_to = "sample", values_to = "expr") %>%
  filter(str_detect(sample, "_1$")) %>%   # 🔥 KEEP ONLY gRG
  mutate(
    genotype = case_when(
      str_detect(sample, "-wt-") ~ "WT",
      str_detect(sample, "-ko-") ~ "KO",
      str_detect(sample, "-het-") ~ "HET"
    )
  )

df_long$genotype <- factor(df_long$genotype, levels = c("WT", "HET", "KO"))

# ==================================================
# 5. Z-score per gene (ONLY within gRG)
# ==================================================
df_z <- df_long %>%
  group_by(gene) %>%
  mutate(zscore = scale(expr)[,1]) %>%
  ungroup() %>%
  filter(is.finite(zscore))

# ==================================================
# 6. Sample-level mean
# ==================================================
df_avg <- df_z %>%
  group_by(sample, genotype) %>%
  summarise(mean_z = mean(zscore), .groups = "drop")

# ==================================================
# 7. Plot (cloud + boxplot + points)
# ==================================================
p <- ggplot() +
  
  # 🔥 Background gene cloud
  geom_jitter(
    data = df_z,
    aes(x = genotype, y = zscore),
    width = 0.15,
    size = 0.6,
    alpha = 0.15,
    color = "black"
  ) +
  
  # Boxplot (sample means)
  geom_boxplot(
    data = df_avg,
    aes(x = genotype, y = mean_z, fill = genotype),
    width = 0.5,
    alpha = 0.7,
    outlier.shape = NA
  ) +
  
  # Sample points
  geom_jitter(
    data = df_avg,
    aes(x = genotype, y = mean_z),
    width = 0.08,
    size = 2,
    color = "black"
  ) +
  
  scale_fill_manual(values = c(
    "WT" = "#4DBBD5",
    "HET" = "#E64B35",
    "KO" = "#00A087"
  )) +
  
  theme_classic() +
  labs(
    title = "gRG Maturation Signature (Z-score)",
    y = "Z-score",
    x = "Genotype"
  )

# ==================================================
# 8. Save
# ==================================================
ggsave(
  "/Users/satabdi/Downloads/neuronal_maturation_boxplot.pdf",
  p,
  width = 5,
  height = 4
)
