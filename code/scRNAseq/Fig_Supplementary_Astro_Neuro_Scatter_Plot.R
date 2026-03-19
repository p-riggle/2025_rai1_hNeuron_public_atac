# ==================================================
# Load Libraries
# ==================================================

library(dplyr)
library(ggplot2)
library(hexbin)
library(gridExtra)

# ==================================================
# Read DESeq2 result tables
# ==================================================

astro_het <- read.csv("/Users/satabdi/Downloads/rstbl_Astro_wt.vs.het.csv")
astro_ko  <- read.csv("/Users/satabdi/Downloads/rstbl_Astro_wt.vs.ko.csv")

neuron_het <- read.csv("/Users/satabdi/Downloads/rstbl_Neuron_wt.vs.het.csv")
neuron_ko  <- read.csv("/Users/satabdi/Downloads/rstbl_Neuron_wt.vs.ko.csv")

# ==================================================
# Rename first column to gene
# ==================================================

astro_het  <- astro_het  %>% rename(gene = 1)
astro_ko   <- astro_ko   %>% rename(gene = 1)

neuron_het <- neuron_het %>% rename(gene = 1)
neuron_ko  <- neuron_ko  %>% rename(gene = 1)

# ==================================================
# Merge HET and KO tables
# ==================================================

astro_all <- astro_het %>%
  inner_join(astro_ko, by = "gene", suffix = c("_het", "_ko"))

neuron_all <- neuron_het %>%
  inner_join(neuron_ko, by = "gene", suffix = c("_het", "_ko"))

# ==================================================
# Filter DEGs (padj < 0.05 )
# ==================================================

astro_deg <- astro_all %>%
  filter(padj_het < 0.05 | padj_ko < 0.05)

neuron_deg <- neuron_all %>%
  filter(padj_het < 0.05 | padj_ko < 0.05)

# ==================================================
# Global axis (for uniform comparison)
# ==================================================

global_limit <- max(
  abs(astro_deg$log2FoldChange_het),
  abs(astro_deg$log2FoldChange_ko),
  abs(neuron_deg$log2FoldChange_het),
  abs(neuron_deg$log2FoldChange_ko),
  na.rm = TRUE
)

global_limit <- ceiling(global_limit)

# ==================================================
# Neuron-specific axis (zoomed view)
# ==================================================

neuron_limit <- max(
  abs(neuron_deg$log2FoldChange_het),
  abs(neuron_deg$log2FoldChange_ko),
  na.rm = TRUE
)

neuron_limit <- ceiling(neuron_limit)

# ==================================================
# Scatterplot function 
# ==================================================

make_hex_scatter <- function(df, title, limit, subtitle = NULL) {
  
  ggplot(df,
         aes(x = log2FoldChange_het,
             y = log2FoldChange_ko)) +
    
    stat_binhex(bins = 60) +
    
    geom_abline(
      intercept = 0,
      slope = 1,
      color = "gray40",
      linetype = "dashed"
    ) +
    
    
    geom_smooth(
      method = "lm",
      se = FALSE,
      color = "blue",
      linewidth = 0.8,
      alpha = 0.6
    ) +
    
    coord_fixed(
      xlim = c(-limit, limit),
      ylim = c(-limit, limit)
    ) +
    
    theme_bw() +
    
    theme(
      text = element_text(size = 14, family = "Helvetica")
    ) +
    
    labs(
      x = "log2 Fold Change (WT vs HET)",
      y = "log2 Fold Change (WT vs KO)",
      title = title,
      subtitle = subtitle
    )
}

# ==================================================
# Generate plots
# ==================================================

astro_scatter <- make_hex_scatter(
  astro_deg,
  "gRG",
  global_limit
)

neuron_scatter <- make_hex_scatter(
  neuron_deg,
  "Neuron",
  global_limit
)

neuron_zoom_scatter <- make_hex_scatter(
  neuron_deg,
  "Neuron (Zoomed)",
  neuron_limit,
  subtitle = "Axis scaled to neuron-specific range"
)

# ==================================================
# Display plots
# ==================================================

grid.arrange(
  astro_scatter,
  neuron_scatter,
  neuron_zoom_scatter,
  ncol = 3
)

# ==================================================
# Save plots
# ==================================================

pdf(
  sprintf("/Users/satabdi/Downloads/%s-log2fc_scatter_Astro_Neuron_zoom.pdf", Sys.Date()),
  width = 15,
  height = 5
)

grid.arrange(
  astro_scatter,
  neuron_scatter,
  neuron_zoom_scatter,
  ncol = 3
)

dev.off()
