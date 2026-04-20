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

# ==================================================
# 2. Gene list
# ==================================================
astro_gene_text <- "ENO1
SLC2A5
NECAP2
PADI2
ALDH4A1
TMCO4
PLA2G5
GRHL3
CLIC4
HPCAL4
ELOVL1
TESK2
CCDC163P
PRDX1
MOB3C
TTC39A
PODN
CYP2J2
AK4
PTGFR
GBP1
GBP2
TGFBR3
TMEM56
AHCYL1
CHI3L2
ATP1A1
PDE4DIP
SELENBP1
LINGO4
S100A10
IL6R
PBXIP1
GBA
UHMK1
SELL
PRDX6
KIAA0040
BRINP2
NPL
CSRP1
LGR6
ADORA1
CHI3L1
NFASC
ITPKB
C1orf198
PFKP
AKR1C1
ITIH5
FAM171A1
SLC39A12
NEBL
ENKUR
GPR158
RASSF4
VSTM4
OGDHL
C10orf105
SPOCK2
C10orf11
KCNMA1
ANXA11
HECTD2-AS1
PGAM1
HOGA1
PPAPDC1A
KNDC1
TPP1
AMPD3
MRVI1
ZBED5-AS1
GALNT18
NAV2
FIBIN
CD44
PAMR1
FAM180B
APLNR
SLC43A1
CYB561A3
SLC22A6
PYGM
SLC29A2
PC
PRCP
TMPRSS5
ZBTB16
SCN2B
SORL1
GAPDH
C1S
C1R
A2M
EMP1
SLCO1A2
SPX
LMNTD1
CNTN1
ENDOU
SLC48A1
PFKM
SPRYD3
ZNF385A
METTL7B
NDUFA4L2
B4GALNT1
MYBPC1
DAO
RPH3A
OAS3
SDS
SDSL
TESC
PEBP1
HSPB8
HIP1R
GPC5-AS1
ABHD4
MYH7
EGLN3
HIF1A
HSPA2
RGS6
PAPLN
SLC24A4
SERPINA3
RYR3
C15orf52
MAP1A
GATM
TNFAIP8L3
CGNL1
ITGA11
PKM
ADAMTSL3
NMB
BAIAP3
GPRC5B
ACSM5
SULT1A1
MVP
MT1E
MT1M
MT1F
MT1G
TRADD
LCAT
HPR
CLEC18B
CHST6
DBNDD1
ABR
CXCL16
PLD2
FGF11
RASD1
MFAP4
SLC47A2
PIPOX
PPP1R1B
PTRF
IFI35
CCDC103
DCAKD
PNPO
COPZ2
NXPH3
PDK2
44808
TEX2
RGS9
SLC39A11
SMIM5
ITGB4
WBP2
C1QTNF1
TWSG1
SLC14A1
MRO
PRTN3
CAPS
C3
ANGPTL4
ICAM1
CYP4F3
CYP4F12
CYP4F11
TMEM38A
NWD1
CRLF1
KLHL26
LINC00906
LINC01532
SLC7A10
PLA2G4C
GYS1
CAMKMT
TTC7A
CALM2
GPR75
SPTBN1
EFEMP1
CCDC85A
TEX261
DCTN1
VAMP5
KCNIP3
FAHD2CP
FAHD2B
PSD4
CYBRD1
TTC30A
IGFBP5
ATG9A
PTPRN
SCG2
ACKR3
NEU4
SIRPA
PYGB
LBP
SDC4
SYS1
EYA2
TIAM1
C21orf62
B3GALT5-AS1
MX1
AIFM3
SEC14L2
LIMK2
MYH9
ACO2
LMCD1
HDAC11
DAZL
DLEC1
GORASP1
FAM198A
CLEC3B
LARS2
LTF
TCTA
MAPKAPK3
ABHD14B
ADAMTS9
GBE1
COL8A1
LINC01279
CCDC80
SIDT1
MYLK
ALDH1L1-AS2
CLSTN2
SLC7A14
TNIK
ECE2
DGKG
CPLX1
WFS1
PPP2R2C
CPZ
GABRG1
LINC01088
DDIT4L
SETD7
EDNRA
NR3C2
MFAP3L
PDZD2
PLCXD3
F2RL2
LINC01338
LOC101927078
GRAMD3
ALDH7A1
ABLIM3
ARHGEF37
CD74
GFPT2
NRN1
GMPR
NRSN1
HLA-E
HLA-DRA
HLA-DRB5
HLA-DRB6
HLA-DRB1
HLA-DQA2
HLA-DMB
HLA-DPA1
ENPP4
ENPP5
PAQR8
COL21A1
BCKDHB
TPD52L1
PDE7B
HIVEP2
SOD2
LINC00473
ITGB8
OSBPL3
AQP1
PGAM2
RAMP3
TNS3
PHKG1
CALN1
AGFG2
SERPINE1
NAMPT
NRCAM
IRF5
AHCYL2
STRIP2
PLXNA4
ATP6V0E2-AS1
TMEM176B
TMEM176A
FAM167A
DLC1
HR
PEBP4
DPYSL2
ADRA1A
SLC20A2
C8orf34-AS1
C8orf34
TPD52
CPQ
MATN2
RNF19A
BAALC-AS1
GSDMC
ADCY8
WISP1
ST3GAL1
LOC101927822
FAM135B
KANK1
SH3GL2
UNC13B
KLF9
GDA
ALDH1A1
TMOD1
ALAD
TLR4
BRINP1
AK1
PHYHD1
ADAMTSL2
ARSF
PDK3
PRRG1
MAOB
KCND1
HEPH
AR
SYTL4
SLC25A53
ARHGAP36
IGSF1
BGN
MPP1"

genes_use <- astro_gene_text %>%
  str_split("\n") %>%
  unlist() %>%
  str_trim() %>%
  .[. != ""] %>%
  toupper()

# ==================================================
# 3. Subset genes
# ==================================================
df_sub <- df %>%
  filter(gene %in% genes_use)

# ==================================================
# 4. Keep ONLY gRG (_5)
# ==================================================
df_long <- df_sub %>%
  pivot_longer(-gene, names_to = "sample", values_to = "expr") %>%
  filter(str_detect(sample, "_5$")) %>%   
  mutate(
    genotype = case_when(
      str_detect(sample, "-wt-") ~ "WT",
      str_detect(sample, "-ko-") ~ "KO",
      str_detect(sample, "-het-") ~ "HET"
    )
  )

df_long$genotype <- factor(df_long$genotype, levels = c("WT", "HET", "KO"))

# ==================================================
# 5. Z-score per gene ( within gRG)
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
  
  # Background gene cloud
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
  "/Users/satabdi/Downloads/gRG_maturation_boxplot_late.pdf",
  p,
  width = 5,
  height = 4
)