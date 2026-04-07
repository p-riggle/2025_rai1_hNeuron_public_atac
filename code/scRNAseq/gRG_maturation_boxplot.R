library(tidyverse)

# -----------------------------
# 1. Read data
# -----------------------------
df <- read.csv("/Users/satabdi/Downloads/DEseq_rlogNorm copy.csv",
               check.names = FALSE)


df$gene <- toupper(df$gene)

# -----------------------------
# 2. gRG maturation-gene set
# -----------------------------
astro_gene_text <-  "Myorg
ALDH1L1
AQP4
GFAP
KCNJ10
GJA1
SPARCL1
MLC1
Gabrg1
Slitrk2
Clasp2
Gramd3
Ctsa
Derl1
Hagh
P4ha1
Pgm1
Fam13c
Selenbp1
Plcd4
Smpd2
Cyp2j6
Dpm1
Dbx2
Wscd1
Cbr3
Tmem14a
H2-K1
Snhg20
Cpeb4
B4gat1
Cnp
Rnd2
Stxbp3
Rora
Gaa
Enpp5
Sorl1
Pdzph1
Aig1
Pcsk1n
Sirt2
Chpt1
Phactr3
Hepacam
Rgs20
Dbndd2
Slc6a1
Slc27a1
Ahcyl1
Plaat3
Shisa4
Aldh2
Pnisr
mt-Nd4
mt-Co1
mt-Cytb
mt-Nd2
mt-Co3
mt-Co2
mt-Atp6
Gstm5
Chchd10
Apoe
Clu
Slc1a2
Aldoc
Ldhb
Enho
Slc7a10
Gstm1
S100a1
Ntrk2
Luzp2
Atp1b2
Nkain4
Gpm6a
Fjx1
Ndrg2
Abhd3
Hacd2
Cadm2
Chst1
Chst2
Vegfa
Rnf13
Serinc1
Entpd2
Sdc4
Mgll
Car2
Fzd2
Timp4
Rsrp1
Lap3
Mir100hg
Ezr
Tlcd1
Ctsd
Il18
Mt2
Gnao1
Arhgap5
Fxyd1
Sat1
Rorb
mt-Nd4l
mt-Nd5
Dtna
Tmem47
Macf1
Camk2n1
Gria2
Itm2c
Ugp2
Csrp1
Fbxo2
Lsamp
Acsbg1
Prxl2a
Bicd1
Rnpc3
3110056K07Rik
Mir124-2hg
4632427E13Rik
Ppp2r3a
Soga3
Fuom
Lysmd2
Pnpla2
Fahd2a
Upp1
Porcn
Syngr1
Abcb9
Tmem176a
Vmp1
Tm7sf3
Cyp2d22
Hsd17b11
Tprkb
Hopx
Gpld1
Ogt
Tmcc1
Ifi27
Clmn
Gdf10
Lgr4
Nupr1
Cryab
Adk
Naaa
Mertk
Omg
Phkg1
Thrsp
Sqstm1
Hsd11b1
Pfkp
Hapln1
Fgfr3
Soat1
Tpp1
Ephx2
Ctsf
Paqr7
Hspa2
Phyh
Leng8
Syne1
Gpam
Scara3
Tle4
Gprc5b
Gm6145
Dst
Daam2
AI413582
4930402H24Rik
Ctso
Tex264
Ppp1r3g
Osbpl1a
Pfkm
Nptn
Add3
Slc41a1
Ngef
Smpd1
Dnajb9
Sirpa
Slc25a33
Mbp
Fbxo44
Tnfsf12
Ctbs
Ndp
Sec14l2
Fuca2
Mink1
4930488L21Rik
Fezf2
6430590A07Rik
Chil1
Lypd6
Tspan15
Tef
Prdm16os
Snx32
Mcur1
G0s2
Mmgt2
Laptm4b
4933431E20Rik
Mfn1
Agl
Negr1
Maob
Fam184a
Hexb
Gldc
Ephx1
Ankrd29
Gm44257
Serpinb1a
Nr1d2
Kalrn
Nkiras1
Mttp
Acat3
Fah
Gm12326
Elmo2
1700047M11Rik
Nat8
Nphp3
Tppp3
Pygm
Ccdc190
Gas6
Lmcd1
Slc35f6
Tmbim1
Usp53
Izumo4
Wnk2
B930095G15Rik
Pdk2
Snx14
Fgf1
Abr
Aifm3
Gm2115
Smim3
Pfkfb4
Cideb
Plp1
Sdsl
Ass1
Rcan2
Chst10
Fam228b
Cyp4f13
Gk
Fmn2
Usp54
Arhgef19
Sik3
Abhd14b
Nat8f5
Mast4
Emc9
Phlpp1
Sgip1
Abcc5
Tnik
Por
Eps15
Ddhd1
Pitpnc1
Fnbp1
Ttc14
Pttg1
Trpm3
Slc20a1
Gpr146
Lgi4
Cyp4f14
Hlf
Hivep3
Cyp2j9
1190005I06Rik
Pik3ip1
Caskin1
Slc48a1
Per3
Rarres2
Fam20a
Npc1
Kcna2
Cry2
Cyp7b1
Pcdh9
Adora2b
Igsf11
Ppp1r1b
Il33
Tuba4a
Neat1
Washc2
Gpr162
Crlf3
Gpr37
Kirrel2
Cgrrf1
Efcab14
Dhrs3
Mfap3l
St3gal4
Rpe65
Pde8b
Col11a2
Aldh1a1
Glrb
Sod3
Efemp1
Cyp4f15
Nr1d1
Kcnn2
Wasf3
Eps8
Mt1
Atp1a2
Plpp3
Glul
Prdx6
Btbd17
Cxcl14
Htra1
F
Scg3
Pla2g7
Prnp
S1pr1
Ntsr2
Acsl3
F3
Gpr37l1
Cldn10
Pmp22
Tagln3
Ank2
Tlr3
Oaf
Wls
Lcat
Grin2c
2900052N01Rik
mt-Atp8
Gabrb1
Dbp
Prodh
Agpat5
Appl2
Gpc5
Rida
Eva1a
Lix1
Cntfr
Cystm1
Slc39a12
Cadm1
Gjb6
Prex2
Fam107a
Malat1
mt-Nd3
S100b
Id4
Nat8f1
Dhrs7
Grina
Cmtm5
Kcnk1
Aplp1
Acsl6"

# -----------------------------
# 3. Convert to gene vector
# -----------------------------
gRG_maturation_genes <- astro_gene_text %>%
  str_split("\n") %>%
  unlist() %>%
  str_trim() %>%
  .[. != ""] %>%
  toupper()

# -----------------------------
# 4. Keep only genes present in dataset
# -----------------------------
gRG_maturation_genes <- intersect(df$gene, gRG_maturation_genes)

cat("Genes used:", length(gRG_maturation_genes), "\n")

# -----------------------------
# 5. Subset data
# -----------------------------
df_sub <- df %>%
  filter(gene %in% gRG_maturation_genes)

# -----------------------------
# 6. Z-score per gene
# -----------------------------
df_z <- df_sub %>%
  column_to_rownames("gene") %>%
  t() %>%
  scale() %>%
  t() %>%
  as.data.frame()

# -----------------------------
# 7. Convert to long format
# -----------------------------
df_long <- df_z %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "zscore")

# -----------------------------
# 8. Extract metadata
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
# 9. Compute signature score (median = robust)
# -----------------------------
df_avg <- df_long %>%
  group_by(sample, genotype, cell_type) %>%
  summarise(mean_z = median(zscore), .groups = "drop")

# -----------------------------
# 10. Plot (boxplot)
# -----------------------------
p <- ggplot(df_avg, aes(x = genotype, y = mean_z, fill = genotype)) +
  geom_boxplot(alpha = 0.7) +
  geom_jitter(width = 0.15, size = 2) +
  facet_wrap(~cell_type) +
  theme_classic() +
  labs(
    title = "gRG Maturation Signature (Z-score)",
    y = "Mean Z-score",
    x = "Genotype"
  )


print(p)

# -----------------------------
# 11. Save plot
# -----------------------------
ggsave(
  "/Users/satabdi/Downloads/gRG_maturation_boxplot_zscore.pdf",
  plot = p,
  width = 6,
  height = 4
)

