# Load libraries
library(Seurat)
library(BPCells)
library(tidyverse)

# examine the expression of known markers

geo_so = readRDS('data/scRNAseq/iwaselab/rdata/geo_so_sct_unintegrated_with_markers_0.1.rds')
geo_markers = readRDS('data/scRNAseq/iwaselab/rdata/geo_markers_0.1.rds')


# -------------------------------------------------------------------------
# Create lists of markers from Zhong et al. 2018 Nature : 8-26wk embryo prefrontal cortex

zhong_markers = list()
zhong_markers[['Microglia']] = c('PTPRC', 'P2RY12') 
zhong_markers[['NPC']] = c('PAX6','SFRP1')
zhong_markers[['OPC']] = c('OLIG1', 'PDGFRA','COL20A1','PMP2')
zhong_markers[['Excitatory neurons']] = c('NEUROD2', 'RBFOX1')
zhong_markers[['Interneurons']] = c('GAD1', 'PDE4DIP')
zhong_markers[['Astrocytes']] = c( 'GFAP','AQP4','SLCO1C1') 

zhong_markers

# -------------------------------------------------------------------------
# Plot zhong_marker expression in each cluster
zhong_markers_plot = DotPlot(geo_so, features = zhong_markers, assay = 'SCT') +
  theme(panel.background = element_rect(fill = "white"),
        text=element_text(size=10),
        axis.text.x = element_text(angle = 45, vjust = 0.5))

zhong_markers_plot = DotPlot(geo_so, features = zhong_markers, assay = 'SCT') +
  theme_bw() +
  theme(text=element_text(size=10), axis.text.x = element_text(angle = 45, vjust = 0.5))

zhong_markers_plot 

ggsave(filename = 'results/scRNAseq/3_markers/zhong_markers_0.1.png', 
       plot = zhong_markers_plot, 
       width = 12, height = 5, units = 'in') 

# plot more Zhong markers
zhong_opc = c('BCAN', 'OLIG1', 'SCRG1', 'DBI', 'TMSB10', 'PDGFRA', 'S100B', 'NKAIN4', 'PTPRZ1', 'APOD', 'NOVA1', 'PTN', 'TUBB3', 'CD24', 'TSPAN7', 'BCHE', 'COL20A1', 'GPM6B', 'LUZP2', 'SOX4', 'STMN2', 'RTN1', 'RAB31', 'RAMP1', 'B2M', 'PMP2', 'BCL11A', 'CD9', 'EPN2', 'ARL4A', 'ITGB8', 'SERPINE2', 'LIMA1', 'ITM2C', 'MLLT11', 'COL9A1', 'CNP', 'OLIG2', 'SLITRK2', 'GATM', 'TTYH1', 'CD63', 'LMF1', 'PSAT1', 'SCD5', 'TRIB2', 'CST3', 'SEMA5A', 'NEUROD6', 'DUSP6', 'MT3', 'NKX2-2', 'PON2', 'ETV5', 'SOX10', 'LRRC4C', 'SYT1', 'DCX', 'CSPG5', 'CSRP2', 'CNTN1', 'IGFBPL1', 'PCDH15', 'ELAVL4', 'CA10', 'EDNRB', 'OCIAD2', 'SEZ6L', 'NXPH1', 'RFTN2', 'SIRT2', 'PDE4B', 'SPRY4', 'ANXA5', 'SPRY1', 'PPAP2B', 'RND3', 'SCN1A', 'TSC22D4', 'ETV1', 'PLEKHB1', 'PLLP', 'ARC', 'GRID2', 'GNG7', 'SPARC', 'TMEM161B-AS1', 'FIBIN', 'RHOC', 'NSG1', 'OMG', 'KANK1', 'LINC01551', 'NEUROD2', 'SHTN1', 'ALCAM', 'C8orf46', 'SLC35F1', 'SOX8', 'LHFPL3', 'BRINP3', 'C2orf80', 'CADM2', 'SPATA6', 'EPB41L2', 'LSAMP', 'SOX6', 'GPR17', 'ENC1', 'MMP2', 'PLAT', 'ABRACL', 'ANGPTL2', 'CD82')

zhong_opc_plot = DotPlot(geo_so, features = zhong_opc, assay = 'SCT') +
  labs(title = "Zhong_opc_markers") +
  theme_bw() +
  theme(text=element_text(size=10), axis.text.x = element_text(angle = 90, vjust = 0.5))
zhong_opc_plot

zhong_astro = c('CST3', 'CLU', 'ATP1A2', 'PTN', 'SLC1A3', 'BCAN', 'HOPX', 'PON2', 'DBI', 'TNC', 'ATP1B2', 'GPM6B', 'PTPRZ1', 'S100A13', 'TTYH1', 'SPARCL1', 'PPAP2B', 'MT3', 'C1orf61', 'CXCL14', 'FABP7', 'PSAT1', 'MLC1', 'SLC1A2', 'NCAN', 'ADGRG1', 'ALDOC', 'CD9', 'TSC22D4', 'RAB31', 'STMN1', 'B2M', 'MT2A', 'RAMP1', 'EDNRB', 'SLCO1C1', 'ITM2C', 'PEA15', 'S100A16', 'TMSB10', 'GATM', 'PTMA', 'GLUL', 'SPON1', 'F3', 'SOX4', 'LGALS3', 'FAM107A', 'MFGE8', 'TUBB', 'EEPD1', 'CD63', 'GLUD1', 'APOE', 'NKAIN4', 'CNN3', 'BTBD17', 'S1PR1', 'TFPI', 'CSPG5', 'OLFM2', 'PLTP', 'SEZ6L', 'GRM3', 'SPRY1', 'MLLT11', 'PMP22', 'LRRC17', 'HES1', 'ID1', 'SLC4A4', 'AGT', 'SCRG1', 'GFAP', 'LITAF', 'TSPAN7', 'TSPAN3', 'HSPB1', 'ZFP36L2', 'HES5', 'GDPD2', 'LRIG1', 'RFX4', 'TUBB3', 'QKI', 'MSMO1', 'SERPINE2', 'COMT', 'DNER', 'HEY1', 'CYP51A1', 'SLC1A4', 'WSCD1', 'VIM', 'HN1', 'ANXA5', 'PAQR8', 'SLITRK2', 'HTRA1', 'ID4', 'CD82', 'FABP5', 'SDC3', 'AK4', 'MT1E', 'CD24', 'DTNA', 'HADHB', 'LUZP2', 'OAT', 'SIRPA', 'BMP7', 'IDI1', 'ADD3', 'CLDN10', 'LINC00943', 'SOX9', 'DDAH1', 'AIF1L', 'NPNT', 'CCDC80', 'FERMT2', 'CDC42EP4', 'SLC6A1', 'SLC9A3R1', 'BEX1', 'DIO2', 'SLC15A2', 'ZFP36L1', 'TMSB15A', 'DOK5', 'FLJ46906', 'IGFBPL1', 'RTN1', 'RGS20', 'SEMA6A', 'ADCYAP1R1', 'DDR1', 'SNX3', 'AQP4', 'PHGDH', 'ETV5', 'LYN', 'PHYHIPL', 'SCD', 'S100B', 'BCHE', 'ELOVL5', 'NDRG2', 'GNG7', 'KCNE5', 'MYO10', 'ETV1', 'ID2', 'LGALS3BP', 'LRRC3B', 'SOX11', 'C8orf4', 'PRCP', 'SFXN5', 'HLA-A', 'HSD17B12', 'LRP4', 'TIMP3', 'OCIAD2', 'RHOC', 'ANOS1', 'ELOVL2', 'GNG5', 'HAPLN3', 'TRIB2', 'ACSBG1', 'SLC3A2', 'DHCR7', 'PLEKHB1', 'SLC38A3', 'YIF1A', 'PPP1R16A', 'SOX2', 'TM7SF2', 'LIX1', 'RCN1', 'TMEM132B', 'IFI27L1', 'PTTG1IP', 'BCL11A', 'CBR1', 'NLGN3', 'CDC42EP1', 'TAGLN3', 'LGALS1', 'FBLN2', 'LRRN3', 'LRRTM3', 'ABHD3', 'CPNE5', 'SPRY2', 'HMGCS1', 'RGMA', 'C3orf70', 'JUN', 'UCHL1', 'GJA1', 'MEIS1', 'ALCAM', 'TRIM47', 'UBL3', 'PHACTR2', 'SLC25A18', 'SYPL1', 'HNMT', 'ALDH6A1', 'DARS', 'METTL7A', 'ADORA2B', 'CBS', 'EFHD2', 'GULP1', 'PREX1', 'DCXR', 'GNG11', 'RASSF4', 'ABCD3', 'GNPTG', 'LINC01351', 'RAB34', 'MGLL', 'SNX5', 'ITPR2', 'PDLIM5', 'DKK3', 'PDPN', 'TRIL', 'ZFYVE21', 'P4HA1', 'RIT2', 'TECPR2', 'SERPINB6', 'STON2', 'CSRP2', 'SEMA6D', 'GPRC5B', 'PRSS35', 'VEPH1', 'EFNA1', 'PDLIM3', 'HAPLN1', 'SLC13A5', 'RAP1GAP', 'GNG12', 'TIMP1', 'HDAC2', 'ELAVL4', 'GRB14', 'NTRK2', 'VCAM1', 'ATP13A4', 'ENHO', 'KIAA1161', 'ARHGEF26', 'ARNTL', 'HLA-C', 'DAG1', 'AASS', 'ADGRB1', 'NEUROD6', 'FADS2', 'FXYD1', 'MRC2', 'NT5E', 'TCF7L2', 'FAM84B', 'GPR37L1', 'ALPL', 'APBB2', 'C2orf72', 'DCN', 'LEPROT', 'NOTCH3', 'TNFRSF19', 'BCAP29', 'PLAGL1', 'PLIN3', 'B3GALTL', 'MASP1', 'PRODH', 'NSG1', 'BDH2', 'CAMK2G', 'IL33', 'KLHDC8A', 'ABAT', 'IFITM3', 'PDGFRB', 'H2AFZ', 'PLCD1', 'ROBO3', 'TBC1D10A', 'BBS2', 'C1orf122', 'CA12', 'LAMB2', 'FAT1', 'TSPAN6', 'ELN', 'LIPA', 'SYT1', 'NOTCH2', 'CRYL1', 'GATSL3', 'EFEMP2', 'LGI1', 'GPR137B', 'PLEC', 'SLC39A12', 'TEX264', 'KCNJ10', 'LINC00998', 'TOR1AIP1', 'BMPR1B', 'DBX2', 'HLA-B', 'RBFOX2', 'SELENBP1', 'SHISA4', 'SOX8', 'SPATA6', 'ACSL6', 'G6PC3', 'GPX3', 'HES4', 'IGFBP4', 'RGCC', 'PBXIP1', 'EMC2', 'LRRC10B', 'TTYH2', 'BRINP3', 'FAM181B', 'ITGAV', 'MAGT1', 'MMD2', 'P3H4', 'RAMP3', 'LIFR', 'PRR5', 'SASH1', 'SOAT1', 'THBS3', 'TP53I3', 'CTSL', 'PNPLA3', 'PPP1R17')
zhong_astro_plot = DotPlot(geo_so, features = zhong_astro, assay = 'SCT') +
  labs(title = "Zhong_astrocyte_markers") +
  theme_bw() +
  theme(text=element_text(size=10), axis.text.x = element_text(angle = 90, vjust = 0.5))
zhong_astro_plot


zhong_npc = c('HMGN2', 'HMGB2', 'STMN2', 'PAX6', 'H2AFZ', 'SMC4', 'TUBA1B', 'RTN1', 'VIM', 'CKS2', 'HSPB1', 'NUSAP1', 'UBE2C', 'PTTG1', 'TOP2A', 'GMNN', 'MAD2L1', 'SFRP1', 'PHGDH', 'CKS1B', 'GNG5', 'CENPF', 'ANP32E', 'KIAA0101', 'SNRPB', 'HES1', 'AURKB', 'CDK1', 'MCM7', 'PCNA', 'BIRC5', 'BTG3', 'FAM64A', 'DCX', 'UBE2T', 'EMX2', 'CENPU', 'ZWINT', 'MCM2', 'PBK', 'LDHA', 'RPA3', 'TMEM98', 'H2AFX', 'CCNB2', 'CLIC1', 'H2AFV', 'RRM2', 'UHRF1', 'GINS2', 'RCN1', 'ATAD2', 'PRC1', 'MIS18BP1', 'CDC20', 'MFAP2', 'RPS27L', 'MKI67', 'TPX2', 'PCSK1N', 'MCM3', 'MDK', 'SYT1', 'PAICS', 'MEF2C', 'C19orf48', 'NUF2', 'RRM1', 'PLIN2', 'RNASEH2A', 'KIF22', 'SHTN1', 'ID4', 'KPNA2', 'GNG2', 'C21orf58', 'GLI3', 'HMGB3', 'TMPO', 'SFRP2', 'SMC2', 'RAB3A', 'DUT', 'GTSE1', 'KIF15', 'CENPN', 'CKAP2', 'DTYMK', 'NDC80', 'CNTLN', 'CREB5', 'CCNB1', 'CENPM', 'KIF11', 'PPAP2B', 'MYT1L', 'MCM4', 'DLGAP5', 'MCM5', 'WDR34', 'CDON', 'EMX1', 'UCP2', 'NRXN1', 'CENPW', 'RUNX1T1', 'CDCA5', 'KIF2C', 'RBFOX2', 'OTX1', 'CCNA2', 'SPC25', 'CLSPN', 'FEN1', 'CENPK', 'CKAP2L', 'FBXO5', 'ITGB3BP', 'NTM', 'ASPM', 'BUB3', 'TAGLN2', 'TUBB4B', 'MAPT', 'MIAT', 'HAT1', 'LMO3', 'CDCA4', 'HELLS', 'TACC3', 'SNAP25', 'SCG3', 'BUB1', 'DEPDC1B', 'HJURP', 'SPAG5', 'ZFHX4', 'CHL1', 'ZEB2', 'CDCA7', 'DACH1', 'KNSTRN', 'NSG1', 'CDKN3', 'SOX3', 'MELK', 'NOTCH1', 'SLA', 'GAP43', 'SGOL2', 'FAM60A', 'ORC6', 'PRR11', 'MEG3', 'KIF23', 'LMO4', 'DHFR', 'HIST1H4C', 'MCM6', 'CDCA3', 'EFNB1', 'SGOL1', 'SOX9', 'APOLD1', 'IGF2BP1', 'TGIF1', 'EOMES', 'OIP5')
zhong_npc_plot = DotPlot(geo_so, features = zhong_npc, assay = 'SCT') +
  labs(title = "Zhong_npc_markers") +
  theme_bw() +
  theme(text=element_text(size=10), axis.text.x = element_text(angle = 90, vjust = 0.5))
zhong_npc_plot

zhong_ex_neuro = c('NEUROD6', 'CSRP2', 'IGFBPL1', 'SLA', 'NEUROD2', 'SOX2-OT', 'PPP2R2B', 'B2M', 'PLS3', 'RBFOX1', 'MAF', 'NEUROD1', 'ARX')
zhong_ex_neuro_plot = DotPlot(geo_so, features = zhong_ex_neuro, assay = 'SCT') +
  labs(title = "Zhong_ex_neuro_markers") +
  theme_bw() +
  theme(text=element_text(size=10), axis.text.x = element_text(angle = 90, vjust = 0.5))
zhong_ex_neuro_plot

zhong_interneuron = c('PLS3', 'ERBB4', 'DLX2', 'PTPRZ1', 'DLX1', 'GAD1', 'PDE4DIP', 'FDFT1', 'MEG3', 'IGFBPL1', 'MEF2C', 'LINC01551', 'SOX2-OT', 'DLX5', 'MAF', 'NEUROD6', 'CXCR4', 'INA', 'GAD2', 'LHX6', 'PHLDA1', 'VIM', 'PPP1R17', 'LDHA', 'PTTG1', 'CELF4', 'HMGB2', 'ARX', 'NRXN3', 'CKS2', 'GNG5', 'UBE2C', 'UBE2T', 'MDK', 'NEUROD1')
zhong_interneuron_plot = DotPlot(geo_so, features = zhong_interneuron, assay = 'SCT') +
  labs(title = "Zhong_interneuron_markers") +
  theme_bw() +
  theme(text=element_text(size=10), axis.text.x = element_text(angle = 90, vjust = 0.5))
zhong_interneuron_plot

library(cowplot)
zhong_all_marker <- plot_grid(zhong_opc_plot, zhong_astro_plot, zhong_npc_plot, zhong_ex_neuro_plot, zhong_interneuron_plot, ncol = 1)

ggsave(filename = 'results/scRNAseq/3_markers/zhong_all_markers_0.1.png', 
       plot = zhong_all_marker, 
       width = 12, height = 24, units = 'in') 

## known marker genes based on Zhong Nature paper, Yao CellStemCell paper, and Louis

library(Seurat)
library(patchwork)
library(ggplot2)

# marker gene lists
StemCell_marker = c('SOX2')
NPC_marker = c('PAX6','RBFOX1', 'DCX')
oRG_marker = c('HOPX', 'FAM107A','TNC','PMP2')
IPC_marker = c('EOMES', 'NHLH1', 'NHLH2', 'TP53I11', 'PPP1R17')
Preplate_marker = c('EOMES', 'TBR1')
Telencephalic_marker = c('LHX2', 'FOXG1', 'FEZF2')
GanglionicEminences_marker = c('DLX1', 'ASCL1', 'GAD1')
mid_hindBrain_marker = c('EN2','PAX7','TFAP2B', 'TH','PTX3', 'POU3F2', 'PAX3', 'ROBO3', 'NTRK2')
Mesenchymal_marker = c('DCN')
Dividing_cell_marker = c('MKI67')

# Create a list of plots for each cell type
stemCellPlot = FeaturePlot(geo_so, features = StemCell_marker, label = TRUE) 
NPCplot = FeaturePlot(geo_so, features = NPC_marker, label = TRUE)
oRGplot = FeaturePlot(geo_so, features = oRG_marker, label = TRUE)
IPCplot = FeaturePlot(geo_so, features = IPC_marker, label = TRUE)
PreplatePlot = FeaturePlot(geo_so, features = Preplate_marker, label = TRUE)
TelencephalicPlot = FeaturePlot(geo_so, features = Telencephalic_marker, label = TRUE)
GanglionicEminencesPlot = FeaturePlot(geo_so, features = GanglionicEminences_marker, label = TRUE)
mid_hindBrainPlot = FeaturePlot(geo_so, features = mid_hindBrain_marker, label = TRUE)
MesenchymalPlot = FeaturePlot(geo_so, features = Mesenchymal_marker, label = TRUE)
Dividing_cellPlot = FeaturePlot(geo_so, features = Dividing_cell_marker, label = TRUE)

stemCellPlot
NPCplot
PreplatePlot
TelencephalicPlot
MesenchymalPlot
Dividing_cellPlot

# save the list of plots in 

ggsave(filename = 'results/scRNAseq/3_markers/1_stem_cell_0.1.png', 
       plot = stemCellPlot,
       width = 4, height = 4, units = 'in')

library(Seurat)
library(patchwork)

# Create a list of individual FeaturePlots for each gene
NPC_individual_plots = list()
for (gene in NPC_marker) {
  NPC_individual_plots[[gene]] = FeaturePlot(geo_so, features = gene, label = TRUE) +
    ggtitle(gene) +
    theme(plot.title = element_text(hjust = 0.5))
}

# Use patchwork to combine all three plots into a single row
NPC_combined_plot = wrap_plots(NPC_individual_plots, ncol = 3) +
  plot_annotation(title = "NPC Markers") +
  theme(plot.title = element_text(size = 14, face = "bold"))

# Save the combined plot to a file
ggsave(
  filename = 'results/scRNAseq/3_markers/2_NPC_0.1.png',
  plot = NPC_combined_plot,
  width = 12, height = 4, units = 'in',
  dpi = 300
)

# oRG markers
# Create a list of individual FeaturePlots for each gene
oRG_individual_plots = list()
for (gene in oRG_marker) {
  oRG_individual_plots[[gene]] = FeaturePlot(geo_so, features = gene, label = TRUE) +
    ggtitle(gene) +
    theme(plot.title = element_text(hjust = 0.5))
}

# Use patchwork to combine all three plots into a single row
oRG_combined_plot = wrap_plots(oRG_individual_plots, ncol = 4) +
  plot_annotation(title = "oRG Markers") +
  theme(plot.title = element_text(size = 14, face = "bold"))

oRG_combined_plot

# Save the combined plot to a file
ggsave(
  filename = 'results/scRNAseq/3_markers/3_oRG_0.1.png',
  plot = oRG_combined_plot,
  width = 16, height = 4, units = 'in',
  dpi = 300
)

# IPC marker plot
# Create a list of individual FeaturePlots for each gene
IPC_individual_plots = list()
for (gene in IPC_marker) {
  IPC_individual_plots[[gene]] = FeaturePlot(geo_so, features = gene, label = TRUE) +
    ggtitle(gene) +
    theme(plot.title = element_text(hjust = 0.5))
}

# Use patchwork to combine all three plots into a single row
IPC_combined_plot = wrap_plots(IPC_individual_plots, ncol = 5) +
  plot_annotation(title = "IPC Markers") +
  theme(plot.title = element_text(size = 14, face = "bold"))

IPC_combined_plot

# Save the combined plot to a file
ggsave(
  filename = 'results/scRNAseq/3_markers/4_IPC_0.1.png',
  plot = IPC_combined_plot,
  width = 20, height = 4, units = 'in',
  dpi = 300
)

# save preplate plot
ggsave(filename = 'results/scRNAseq/3_markers/5_preplate_0.1.png', 
       plot = PreplatePlot,
       width = 8, height = 4, units = 'in')

# Telencephalic markers
# Create a list of individual FeaturePlots for each gene
tel_individual_plots = list()
for (gene in Telencephalic_marker) {
  tel_individual_plots[[gene]] = FeaturePlot(geo_so, features = gene, label = TRUE) +
    ggtitle(gene) +
    theme(plot.title = element_text(hjust = 0.5))
}

# Use patchwork to combine all three plots into a single row
tel_combined_plot = wrap_plots(tel_individual_plots, ncol = 3) +
  plot_annotation(title = "Telencephalic Markers") +
  theme(plot.title = element_text(size = 14, face = "bold"))

# Save the combined plot to a file
ggsave(
  filename = 'results/scRNAseq/3_markers/6_Telencephalic_0.1.png',
  plot = tel_combined_plot,
  width = 12, height = 4, units = 'in',
  dpi = 300
)

# GanglionicEminences_marker plots
gan_individual_plots = list()
for (gene in GanglionicEminences_marker) {
  gan_individual_plots[[gene]] = FeaturePlot(geo_so, features = gene, label = TRUE) +
    ggtitle(gene) +
    theme(plot.title = element_text(hjust = 0.5))
}

# Use patchwork to combine all three plots into a single row
gan_combined_plot = wrap_plots(gan_individual_plots, ncol = 3) +
  plot_annotation(title = "Ganglionic Eminence Markers") +
  theme(plot.title = element_text(size = 14, face = "bold"))

# Save the combined plot to a file
ggsave(
  filename = 'results/scRNAseq/3_markers/7_GanglionicEminence_0.1.png',
  plot = gan_combined_plot,
  width = 12, height = 4, units = 'in',
  dpi = 300
)

# mid/hind brain markers
# oRG markers
# Create a list of individual FeaturePlots for each gene
oRG_individual_plots = list()
for (gene in oRG_marker) {
  oRG_individual_plots[[gene]] = FeaturePlot(geo_so, features = gene, label = TRUE) +
    ggtitle(gene) +
    theme(plot.title = element_text(hjust = 0.5))
}

# Use patchwork to combine all three plots into a single row
oRG_combined_plot = wrap_plots(oRG_individual_plots, ncol = 4) +
  plot_annotation(title = "oRG Markers") +
  theme(plot.title = element_text(size = 14, face = "bold"))

oRG_combined_plot

# Save the combined plot to a file
ggsave(
  filename = 'results/scRNAseq/3_markers/3_oRG_0.1.png',
  plot = oRG_combined_plot,
  width = 16, height = 4, units = 'in',
  dpi = 300
)

# mid/hind brain marker plot
# Create a list of individual FeaturePlots for each gene
mid_hind_individual_plots = list()
for (gene in mid_hindBrain_marker) {
  mid_hind_individual_plots[[gene]] = FeaturePlot(geo_so, features = gene, label = TRUE) +
    ggtitle(gene) +
    theme(plot.title = element_text(hjust = 0.5))
}

# Use patchwork to combine all three plots into a single row
mid_hind_combined_plot = wrap_plots(mid_hind_individual_plots, ncol = 5) +
  plot_annotation(title = "Mid/Hind Brain Markers") +
  theme(plot.title = element_text(size = 14, face = "bold"))

mid_hind_combined_plot

# Save the combined plot to a file
ggsave(
  filename = 'results/scRNAseq/3_markers/8_mid_hind_0.1.png',
  plot = mid_hind_combined_plot,
  width = 20, height = 8, units = 'in',
  dpi = 300
)

# save mescenchymalplot
ggsave(
  filename = 'results/scRNAseq/3_markers/9_mesenchymal_0.1.png',
  plot = MesenchymalPlot,
  width = 4, height = 4, units = 'in',
  dpi = 300
)

# save dividing cell plot
ggsave(
  filename = 'results/scRNAseq/3_markers/10_dividing_cell_0.1.png',
  plot = Dividing_cellPlot,
  width = 4, height = 4, units = 'in',
  dpi = 300
)

