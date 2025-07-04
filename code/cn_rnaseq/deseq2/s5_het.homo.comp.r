# compare gene expression changes between RAI1-HET and RAI1-KO. Same set of genes changed?
 # 1. Venn
 # 2. Scatter

 # read restable
d0_wt.vs.het <- read.delim("results/deseq2/rstbl_d0_wt.vs.het.txt")
d0_wt.vs.ko <- read.delim("results/deseq2/rstbl_d0_wt.vs.ko.txt")
d20_wt.vs.het <- read.delim("results/deseq2/rstbl_d20_wt.vs.het.txt")
d20_wt.vs.ko <- read.delim("results/deseq2/rstbl_d20_wt.vs.ko.txt")
d34_wt.vs.het <- read.delim("results/deseq2/rstbl_d34_wt.vs.het.txt")
d34_wt.vs.ko <- read.delim("results/deseq2/rstbl_d34_wt.vs.ko.txt")
d60_wt.vs.het <- read.delim("results/deseq2/rstbl_d60_wt.vs.het.txt")
d60_wt.vs.ko <- read.delim("results/deseq2/rstbl_d60_wt.vs.ko.txt")

# # 1. Venn
library(BioVenn)
list_z <- NULL

 # create gene lists for venn diagram. list the rownames of d0_wt.vs.het and d0_wt.vs.ko their padj is lower than 0.05 
 d0_het_degenes <- rownames(subset(d0_wt.vs.het, padj < 0.05))
 d0_ko_degenes <- rownames(subset(d0_wt.vs.ko, padj < 0.05))
 d20_het_degenes <- rownames(subset(d20_wt.vs.het, padj < 0.05))
 d20_ko_degenes <- rownames(subset(d20_wt.vs.ko, padj < 0.05))
 d34_het_degenes <- rownames(subset(d34_wt.vs.het, padj < 0.05))
 d34_ko_degenes <- rownames(subset(d34_wt.vs.ko, padj < 0.05))
 d60_het_degenes <- rownames(subset(d60_wt.vs.het, padj < 0.05))
 d60_ko_degenes <- rownames(subset(d60_wt.vs.ko, padj < 0.05))

d0.het.vs.ko <- draw.venn(d0_het_degenes, d0_ko_degenes, list_z, 
title = "DEG Overlap: Day 0", t_f = "Helvetica", 
subtitle = NULL, 
x_c = "#000096", xtitle = "RAI1-HET", xt_f = "Helvetica", 
y_c = "yellow", ytitle = "RAI1-KO", yt_f = "Helvetica", 
nrtype = "abs", nr_f = "Helvetica", 
output = "pdf", 
filename = sprintf("results/het.vs.ko/%s-d0_het.vs.ko-biovenn.pdf", Sys.Date()))

d20.het.vs.ko <- draw.venn(d20_het_degenes, d20_ko_degenes, list_z, 
title = "DEG Overlap: Day 20", t_f = "Helvetica", 
subtitle = NULL, 
x_c = "#000096", xtitle = "RAI1-HET", xt_f = "Helvetica", 
y_c = "yellow", ytitle = "RAI1-KO", yt_f = "Helvetica", 
nrtype = "abs", nr_f = "Helvetica", 
output = "pdf", 
filename = sprintf("results/het.vs.ko/%s-d20_het.vs.ko-biovenn.pdf", Sys.Date()))

d34.het.vs.ko <- draw.venn(d34_het_degenes, d34_ko_degenes, list_z, 
title = "DEG Overlap: Day 34", t_f = "Helvetica", 
subtitle = NULL, 
x_c = "#000096", xtitle = "RAI1-HET", xt_f = "Helvetica", 
y_c = "yellow", ytitle = "RAI1-KO", yt_f = "Helvetica", 
nrtype = "abs", nr_f = "Helvetica", 
output = "pdf", 
filename = sprintf("results/het.vs.ko/%s-d34_het.vs.ko-biovenn.pdf", Sys.Date()))

d60.het.vs.ko <- draw.venn(d60_het_degenes, d60_ko_degenes, list_z, 
title = "DEG Overlap: Day 60", t_f = "Helvetica", 
subtitle = NULL, 
x_c = "#000096", xtitle = "RAI1-HET", xt_f = "Helvetica", 
y_c = "yellow", ytitle = "RAI1-KO", yt_f = "Helvetica", 
nrtype = "abs", nr_f = "Helvetica", 
output = "pdf", 
filename = sprintf("results/het.vs.ko/%s-d60_het.vs.ko-biovenn.pdf", Sys.Date()))
 
# combine the four Venn diagrams into one pdf file --> gave up on this, because it is not possible to combine pdf files with Rv
testvenn <- draw.venn(d0_het_degenes, d0_ko_degenes, list_z, 
title = "DEG Overlap: Day 0", t_f = "Helvetica", 
subtitle = NULL, 
x_c = "#000096", xtitle = "RAI1-HET", xt_f = "Helvetica", 
y_c = "yellow", ytitle = "RAI1-KO", yt_f = "Helvetica", 
nrtype = "abs", nr_f = "Helvetica", 
output = "svg", 
filename = sprintf("results/het.vs.ko/%s-d0_het.vs.ko-biovenn.svg", Sys.Date()))

# install.packages("grImport2")
install.packages("rsvg", dependencies = TRUE)
library("rsvg")
library(grImport2)
library(gridExtra)
library(ggplot2)

testv1 <- rsvg_svg("results/het.vs.ko/2017-11-29-d0_het.vs.ko-biovenn.svg")
testv <- readPicture("results/het.vs.ko/2017-11-29-d0_het.vs.ko-biovenn.svg.xml")

# chi square test: practice
data_frame <- read.csv("https://goo.gl/j6lRXD")
table(data_frame$treatment, data_frame$improvement)
chisq.test(data_frame$treatment, data_frame$improvement, correct = FALSE)

# chi square test: real data
# create a data frame with the DEG counts for each time point
genes <- rownames(d0_wt.vs.het)
d0_het <- d0_wt.vs.het$padj < 0.05
d0_ko <- d0_wt.vs.ko$padj < 0.05
d20_het <- d20_wt.vs.het$padj < 0.05
d20_ko <- d20_wt.vs.ko$padj < 0.05
d34_het <- d34_wt.vs.het$padj < 0.05
d34_ko <- d34_wt.vs.ko$padj < 0.05
d60_het <- d60_wt.vs.het$padj < 0.05
d60_ko <- d60_wt.vs.ko$padj < 0.05

df.chi <- data.frame(genes, d0_het, d0_ko, d20_het, d20_ko, d34_het, d34_ko, d60_het, d60_ko)
head(df.chi)
# remove NA values
df.chi <- df.chi[complete.cases(df.chi), ]
d0chi <- chisq.test(df.chi$d0_het, df.d0$d0_ko, correct = FALSE)
d20chi <- chisq.test(df.chi$d20_het, df.chi$d20_ko, correct = FALSE)
d34chi <- chisq.test(df.chi$d34_het, df.chi$d34_ko, correct = FALSE)
d60chi <- chisq.test(df.chi$d60_het, df.chi$d60_ko, correct = FALSE)

library(ggplot2)
# install.packages("ggstance")
library(ggstance)

# create the data frame
X <- c(d0chi$statistic, d20chi$statistic, d34chi$statistic, d60chi$statistic)
p-value <- c(d0chi$p.value, d20chi$p.value, d34chi$p.value, d60chi$p.value)
day <- c("0", "20", "34", "60")
df.chi <- data.frame(day, X, p-value)

# Load required packages
library(ggplot2)
# install.packages("gridExtra")
library(gridExtra)

# create a table of chi square test results and save it as a pdf file
chi.table <- gridExtra::tableGrob(df.chi, rows = NULL, theme = ttheme_default(base_size = 12))
print(chi.table)

pdf(sprintf("results/het.vs.ko/%s-chisq.pdf", Sys.Date()), width = 8, height = 6)
grid.arrange(chi.table, ncol = 1)
dev.off()

# make scatter plots of log2 fold changes of wt.vs.het and wt.vs.ko

# prepare the dataframe with all restables combined
# add sample names to colnames
colnames(d0_wt.vs.het) <- paste0("d0_wt.vs.het_", colnames(d0_wt.vs.het))
colnames(d0_wt.vs.ko) <- paste0("d0_wt.vs.ko_", colnames(d0_wt.vs.ko))
colnames(d20_wt.vs.het) <- paste0("d20_wt.vs.het_", colnames(d20_wt.vs.het))
colnames(d20_wt.vs.ko) <- paste0("d20_wt.vs.ko_", colnames(d20_wt.vs.ko))
colnames(d34_wt.vs.het) <- paste0("d34_wt.vs.het_", colnames(d34_wt.vs.het))
colnames(d34_wt.vs.ko) <- paste0("d34_wt.vs.ko_", colnames(d34_wt.vs.ko))
colnames(d60_wt.vs.het) <- paste0("d60_wt.vs.het_", colnames(d60_wt.vs.het))
colnames(d60_wt.vs.ko) <- paste0("d60_wt.vs.ko_", colnames(d60_wt.vs.ko))

 # combine all restables into one data frame using row.names  
library(dplyr)
all_rstbl <- cbind(d0_wt.vs.het, d0_wt.vs.ko, d20_wt.vs.het, d20_wt.vs.ko, d34_wt.vs.het, d34_wt.vs.ko, d60_wt.vs.het, d60_wt.vs.ko)
colnames(all_rstbl) 
head(all_rstbl)

# make a scatter plot of log2 fold changes of wt.vs.het and wt.vs.ko

# select DEGS for each time point
d0deg <- subset(all_rstbl, d0_wt.vs.het_padj < 0.05 | d0_wt.vs.ko_padj < 0.05)
d20deg <- subset(all_rstbl, d20_wt.vs.het_padj < 0.05 | d20_wt.vs.ko_padj < 0.05)
d34deg <- subset(all_rstbl, d34_wt.vs.het_padj < 0.05 | d34_wt.vs.ko_padj < 0.05)
d60deg <- subset(all_rstbl, d60_wt.vs.het_padj < 0.05 | d60_wt.vs.ko_padj < 0.05)

library(ggplot2)
library(hexbin)

d0_scatter <- ggplot(d0deg, aes(x = d0_wt.vs.het_log2FoldChange, y = d0_wt.vs.ko_log2FoldChange)) + 
  stat_binhex(bins = 60) + 
  geom_abline(intercept = 0, slope = 1, color = "gray", linetype = "dashed") + 
  geom_smooth(method = "lm", se = TRUE, color = "blue") + 
  scale_x_continuous(limits = c(-10, 10)) + 
  scale_y_continuous(limits = c(-10, 10)) +
  theme_bw() + 
  theme(text = element_text(family = "Helvetica")) +
  xlab("log2 fold change (wt vs het)") + 
  ylab("log2 fold change (wt vs ko)") + 
  ggtitle("d0")

d0_scatter

d20_scatter <- ggplot(d20deg, aes(x = d20_wt.vs.het_log2FoldChange, y = d20_wt.vs.ko_log2FoldChange)) + 
  stat_binhex(bins = 60) + 
  geom_abline(intercept = 0, slope = 1, color = "gray", linetype = "dashed") + 
  geom_smooth(method = "lm", se = TRUE, color = "blue") + 
  scale_x_continuous(limits = c(-10, 10)) + 
  scale_y_continuous(limits = c(-10, 10)) +
  theme_bw() + 
  theme(text = element_text(family = "Helvetica")) +
  xlab("log2 fold change (wt vs het)") + 
  ylab("log2 fold change (wt vs ko)") + 
  ggtitle("d20")

d20_scatter

d34_scatter <- ggplot(d34deg, aes(x = d34_wt.vs.het_log2FoldChange, y = d34_wt.vs.ko_log2FoldChange)) + 
  stat_binhex(bins = 60) + 
  geom_abline(intercept = 0, slope = 1, color = "gray", linetype = "dashed") + 
  geom_smooth(method = "lm", se = TRUE, color = "blue") + 
  scale_x_continuous(limits = c(-10, 10)) + 
  scale_y_continuous(limits = c(-10, 10)) +
  theme_bw() + 
  theme(text = element_text(family = "Helvetica")) +
  xlab("log2 fold change (wt vs het)") + 
  ylab("log2 fold change (wt vs ko)") + 
  ggtitle("d34")

d34_scatter

d60_scatter <- ggplot(d60deg, aes(x = d60_wt.vs.het_log2FoldChange, y = d60_wt.vs.ko_log2FoldChange)) + 
  stat_binhex(bins = 60) + 
  geom_abline(intercept = 0, slope = 1, color = "gray", linetype = "dashed") + 
  geom_smooth(method = "lm", se = TRUE, color = "blue") + 
  scale_x_continuous(limits = c(-10, 10)) + 
  scale_y_continuous(limits = c(-10, 10)) +
  theme_bw() + 
  theme(text = element_text(family = "Helvetica")) +
  xlab("log2 fold change (wt vs het)") + 
  ylab("log2 fold change (wt vs ko)") + 
  ggtitle("d60")

d60_scatter

# arrange the four plots in one pdf file
library(gridExtra)
pdf(sprintf("results/het.vs.ko/%s-log2fc_scatter.pdf", Sys.Date()), width = 10, height = 8)
grid.arrange(d0_scatter, d20_scatter, d34_scatter, d60_scatter, ncol = 2)
dev.off()

