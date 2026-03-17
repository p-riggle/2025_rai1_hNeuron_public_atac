install.packages(c("dplyr", "ggplot2", "hexbin", "gridExtra"))

library(dplyr)
library(ggplot2)
library(hexbin)
library(gridExtra)

# read rstbls
astro_het <- read.csv("/Users/satabdi/Downloads//rstbl_Astro_wt.vs.het.csv")
astro_ko  <- read.csv("/Users/satabdi/Downloads/rstbl_Astro_wt.vs.ko.csv")

neuron_het <- read.csv("/Users/satabdi/Downloads/rstbl_Neuron_wt.vs.het.csv")
neuron_ko  <- read.csv("/Users/satabdi/Downloads/rstbl_Neuron_wt.vs.ko.csv")

# merge het and ko tables
astro_het  <- astro_het  %>% rename(gene = 1)
astro_ko   <- astro_ko   %>% rename(gene = 1)
neuron_het <- neuron_het %>% rename(gene = 1)
neuron_ko  <- neuron_ko  %>% rename(gene = 1)

astro_all <- astro_het %>%
  inner_join(astro_ko, by = "gene", suffix = c("_het", "_ko"))

neuron_all <- neuron_het %>%
  inner_join(neuron_ko, by = "gene", suffix = c("_het", "_ko"))

#Filter DEGs (padj < 0.05 in either comparison)
astro_deg <- astro_all %>%
  filter(padj_het < 0.05 | padj_ko < 0.05)

neuron_deg <- neuron_all %>%
  filter(padj_het < 0.05 | padj_ko < 0.05)


#Scatterplot function 
make_hex_scatter <- function(df, title) {
  ggplot(df, aes(x = log2FoldChange_het,
                 y = log2FoldChange_ko)) +
    stat_binhex(bins = 60) +
    geom_abline(intercept = 0, slope = 1,
                color = "gray", linetype = "dashed") +
    geom_smooth(method = "lm", se = TRUE, color = "blue") +
    scale_x_continuous(limits = c(-10, 10)) +
    scale_y_continuous(limits = c(-10, 10)) +
    theme_bw() +
    theme(text = element_text(family = "Helvetica")) +
    labs(
      x = "log2 fold change (wt vs het)",
      y = "log2 fold change (wt vs ko)",
      title = title
    )
}

#Generate plots
astro_scatter  <- make_hex_scatter(astro_deg,  "Astrocyte")
neuron_scatter <- make_hex_scatter(neuron_deg, "Neuron")

astro_scatter
neuron_scatter

#Save Plots
pdf(sprintf("/Users/satabdi/Downloads/%s-log2fc_scatter_Astro_Neuron.pdf", Sys.Date()),
    width = 12, height = 6)   

grid.arrange(astro_scatter, neuron_scatter, ncol = 2)

dev.off()



