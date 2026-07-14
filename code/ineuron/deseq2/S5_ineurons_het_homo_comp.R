# compare gene expression changes between RAI1-HET and RAI1-KO. Same set of genes changed?
# 1. Venn
# 2. Scatter

# read restable
d0_wt.vs.het <- read.delim("/Users/satabdi/Downloads/rstbl_d0_wt.vs.het.txt")
d0_wt.vs.ko <- read.delim("/Users/satabdi/Downloads/rstbl_d0_wt.vs.ko.txt")
d3_wt.vs.het <- read.delim("/Users/satabdi/Downloads/rstbl_d3_wt.vs.het.txt")
d3_wt.vs.ko <- read.delim("/Users/satabdi/Downloads/rstbl_d3_wt.vs.ko.txt")
d7_wt.vs.het <- read.delim("/Users/satabdi/Downloads/rstbl_d7_wt.vs.het.txt")
d7_wt.vs.ko <- read.delim("/Users/satabdi/Downloads/rstbl_d7_wt.vs.ko.txt")


# # 1. Venn
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")
install.packages("BioVenn")
library(BioVenn)

list_z <- NULL

# create gene lists for venn diagram. list the rownames of d0_wt.vs.het and d0_wt.vs.ko their padj is lower than 0.05 
d0_het_degenes <- rownames(subset(d0_wt.vs.het, padj < 0.05))
d0_ko_degenes <- rownames(subset(d0_wt.vs.ko, padj < 0.05))
d3_het_degenes <- rownames(subset(d3_wt.vs.het, padj < 0.05))
d3_ko_degenes <- rownames(subset(d3_wt.vs.ko, padj < 0.05))
d7_het_degenes <- rownames(subset(d7_wt.vs.het, padj < 0.05))
d7_ko_degenes <- rownames(subset(d7_wt.vs.ko, padj < 0.05))


d0.het.vs.ko <- draw.venn(d0_het_degenes, d0_ko_degenes, list_z, 
                          title = "DEG Overlap: Day 0", t_f = "Helvetica", 
                          subtitle = NULL, 
                          x_c = "#000096", xtitle = "RAI1-HET", xt_f = "Helvetica", 
                          y_c = "yellow", ytitle = "RAI1-KO", yt_f = "Helvetica", 
                          nrtype = "abs", nr_f = "Helvetica", 
                          output = "pdf", 
                          filename = sprintf("/Users/satabdi/Downloads/%s-d0_het.vs.ko-biovenn.pdf", Sys.Date()))

d3.het.vs.ko <- draw.venn(d3_het_degenes, d3_ko_degenes, list_z, 
                           title = "DEG Overlap: Day 3", t_f = "Helvetica", 
                           subtitle = NULL, 
                           x_c = "#000096", xtitle = "RAI1-HET", xt_f = "Helvetica", 
                           y_c = "yellow", ytitle = "RAI1-KO", yt_f = "Helvetica", 
                           nrtype = "abs", nr_f = "Helvetica", 
                           output = "pdf", 
                           filename = sprintf("/Users/satabdi/Downloads/%s-d3_het.vs.ko-biovenn.pdf", Sys.Date()))

d7.het.vs.ko <- draw.venn(d7_het_degenes, d7_ko_degenes, list_z, 
                           title = "DEG Overlap: Day 7", t_f = "Helvetica", 
                           subtitle = NULL, 
                           x_c = "#000096", xtitle = "RAI1-HET", xt_f = "Helvetica", 
                           y_c = "yellow", ytitle = "RAI1-KO", yt_f = "Helvetica", 
                           nrtype = "abs", nr_f = "Helvetica", 
                           output = "pdf", 
                           filename = sprintf("/Users/satabdi/Downloads/%s-d7_het.vs.ko-biovenn.pdf", Sys.Date()))



# chi square test: practice
data_frame <- read.csv("https://goo.gl/j6lRXD")
table(data_frame$treatment, data_frame$improvement)
chisq.test(data_frame$treatment, data_frame$improvement, correct = FALSE)

# chi square test: real data
# create a data frame with the DEG counts for each time point
get_padj <- function(df) {
  df[[grep("_padj$", colnames(df), value = TRUE)]]
}

genes <- rownames(d0_wt.vs.het)
d0_het <- get_padj(d0_wt.vs.het) < 0.05
d0_ko  <- get_padj(d0_wt.vs.ko)  < 0.05
d3_het <- get_padj(d3_wt.vs.het) < 0.05
d3_ko  <- get_padj(d3_wt.vs.ko)  < 0.05
d7_het <- get_padj(d7_wt.vs.het) < 0.05
d7_ko  <- get_padj(d7_wt.vs.ko)  < 0.05

df.chi <- data.frame(genes, d0_het, d0_ko, d3_het, d3_ko, d7_het, d7_ko)



head(df.chi)
# remove NA values
df.chi <- df.chi[complete.cases(df.chi), ]
d0chi <- chisq.test(df.chi$d0_het, df.chi$d0_ko, correct = FALSE)
d3chi <- chisq.test(df.chi$d3_het, df.chi$d3_ko, correct = FALSE)
d7chi <- chisq.test(df.chi$d7_het, df.chi$d7_ko, correct = FALSE)

library(ggplot2)
install.packages("ggstance")
library(ggstance)

# create the data frame
X <- c(d0chi$statistic, d3chi$statistic, d7chi$statistic)
p_value <- c(d0chi$p.value, d3chi$p.value, d7chi$p.value)
day <- c("0", "3", "7")
df.chi <- data.frame(day, X, p_value)

# Load required packages
library(ggplot2)
# install.packages("gridExtra")
library(gridExtra)

# create a table of chi square test results and save it as a pdf file
chi.table <- gridExtra::tableGrob(df.chi, rows = NULL, theme = ttheme_default(base_size = 12))
print(chi.table)

pdf(sprintf("/Users/satabdi/Downloads/%s-chisq.pdf", Sys.Date()), width = 8, height = 6)
grid.arrange(chi.table, ncol = 1)
dev.off()

# make scatter plots of log2 fold changes of wt.vs.het and wt.vs.ko

# prepare the dataframe with all restables combined
# add sample names to colnames
colnames(d0_wt.vs.het) <- paste0("d0_wt.vs.het_", colnames(d0_wt.vs.het))
colnames(d0_wt.vs.ko) <- paste0("d0_wt.vs.ko_", colnames(d0_wt.vs.ko))
colnames(d3_wt.vs.het) <- paste0("d3_wt.vs.het_", colnames(d3_wt.vs.het))
colnames(d3_wt.vs.ko) <- paste0("d3_wt.vs.ko_", colnames(d3_wt.vs.ko))
colnames(d7_wt.vs.het) <- paste0("d7_wt.vs.het_", colnames(d7_wt.vs.het))
colnames(d7_wt.vs.ko) <- paste0("d7_wt.vs.ko_", colnames(d7_wt.vs.ko))


# combine all restables into one data frame using row.names  
library(dplyr)
all_rstbl <- cbind(d0_wt.vs.het, d0_wt.vs.ko, d3_wt.vs.het, d3_wt.vs.ko, d7_wt.vs.het, d7_wt.vs.ko)
colnames(all_rstbl) 
head(all_rstbl)

# make a scatter plot of log2 fold changes of wt.vs.het and wt.vs.ko

# select DEGS for each time point
d0deg <- subset(all_rstbl, d0_wt.vs.het_padj < 0.05 | d0_wt.vs.ko_padj < 0.05)
d3deg <- subset(all_rstbl, d3_wt.vs.het_padj < 0.05 | d3_wt.vs.ko_padj < 0.05)
d7deg <- subset(all_rstbl, d7_wt.vs.het_padj < 0.05 | d7_wt.vs.ko_padj < 0.05)

library(ggplot2)
install.packages("hexbin")
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

d3_scatter <- ggplot(d3deg, aes(x = d3_wt.vs.het_log2FoldChange, y = d3_wt.vs.ko_log2FoldChange)) + 
  stat_binhex(bins = 60) + 
  geom_abline(intercept = 0, slope = 1, color = "gray", linetype = "dashed") + 
  geom_smooth(method = "lm", se = TRUE, color = "blue") + 
  scale_x_continuous(limits = c(-10, 10)) + 
  scale_y_continuous(limits = c(-10, 10)) +
  theme_bw() + 
  theme(text = element_text(family = "Helvetica")) +
  xlab("log2 fold change (wt vs het)") + 
  ylab("log2 fold change (wt vs ko)") + 
  ggtitle("d3")

d3_scatter

d7_scatter <- ggplot(d7deg, aes(x = d7_wt.vs.het_log2FoldChange, y = d7_wt.vs.ko_log2FoldChange)) + 
  stat_binhex(bins = 60) + 
  geom_abline(intercept = 0, slope = 1, color = "gray", linetype = "dashed") + 
  geom_smooth(method = "lm", se = TRUE, color = "blue") + 
  scale_x_continuous(limits = c(-10, 10)) + 
  scale_y_continuous(limits = c(-10, 10)) +
  theme_bw() + 
  theme(text = element_text(family = "Helvetica")) +
  xlab("log2 fold change (wt vs het)") + 
  ylab("log2 fold change (wt vs ko)") + 
  ggtitle("d7")

d7_scatter


# arrange the four plots in one pdf file
library(gridExtra)
pdf(sprintf("/Users/satabdi/Downloads/%s-log2fc_scatter.pdf", Sys.Date()), width = 10, height = 8)
grid.arrange(d0_scatter, d3_scatter, d7_scatter, ncol = 2)
dev.off()
