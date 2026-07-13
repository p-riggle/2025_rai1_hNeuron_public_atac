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
# Filter DEGs (padj < 0.05)
# ==================================================

astro_deg <- astro_all %>%
  filter(padj_het < 0.05 | padj_ko < 0.05)

neuron_deg <- neuron_all %>%
  filter(padj_het < 0.05 | padj_ko < 0.05)

# ==================================================
# Scatterplot function
# X-axis = WT vs KO
# Y-axis = WT vs HET
# ==================================================

make_hex_scatter <- function(df, title) {
  
  # ----------------------------------------------
  # Linear regression
  # HET ~ KO
  # ----------------------------------------------
  
  fit <- lm(
    log2FoldChange_het ~ log2FoldChange_ko,
    data = df
  )
  
  # ----------------------------------------------
  # Extract regression statistics
  # ----------------------------------------------
  
  slope <- coef(fit)[2]
  intercept <- coef(fit)[1]
  r2 <- summary(fit)$r.squared
  pval <- summary(fit)$coefficients[2,4]
  
  # ----------------------------------------------
  # Print statistics to console
  # ----------------------------------------------
  
  cat("\n========================\n")
  cat(title, "\n")
  cat("========================\n")
  cat("Slope =", round(slope, 4), "\n")
  cat("Intercept =", round(intercept, 4), "\n")
  cat("R² =", round(r2, 4), "\n")
  cat("P-value =", signif(pval, 3), "\n")
  
  # ----------------------------------------------
  # Label for plot
  # ----------------------------------------------
  
  stats_label <- paste0(
    "Slope = ", round(slope, 3),
    "\nR² = ", round(r2, 3)
  )
  
  # ----------------------------------------------
  # Plot
  # ----------------------------------------------
  
  ggplot(
    df,
    aes(
      x = log2FoldChange_ko,
      y = log2FoldChange_het
    )
  ) +
    
    # Hex bins
    stat_binhex(binwidth = c(0.1, 0.1)) +
    
    # Diagonal reference line
    geom_abline(
      intercept = 0,
      slope = 1,
      color = "gray40",
      linetype = "dashed"
    ) +
    
    # Regression line
    geom_smooth(
      method = "lm",
      se = FALSE,
      color = "blue",
      linewidth = 0.8
    ) +
    
    # Add slope/R² text
    annotate(
      "text",
      x = -4.8,
      y = 4.5,
      label = stats_label,
      hjust = 0,
      size = 4
    ) +
    
    # Fixed axis limits
    coord_fixed(
      xlim = c(-5, 5),
      ylim = c(-5, 5),
      expand = FALSE
    ) +
    
    theme_bw() +
    
    labs(
      x = "log2 Fold Change (WT vs KO)",
      y = "log2 Fold Change (WT vs HET)",
      title = title
    )
}

# ==================================================
# Generate plots
# ==================================================

astro_scatter <- make_hex_scatter(
  astro_deg,
  "gRG"
)

neuron_scatter <- make_hex_scatter(
  neuron_deg,
  "Neuron"
)

# ==================================================
# Display plots
# ==================================================

combined_plot <- grid.arrange(
  astro_scatter,
  neuron_scatter,
  ncol = 2
)

# ==================================================
# Save plots
# ==================================================

ggsave(
  filename = sprintf(
    "/Users/satabdi/Downloads/%s_log2fc_gRG_neuro_scatterplot.pdf",
    Sys.Date()
  ),
  plot = combined_plot,
  width = 10,
  height = 5,
  dpi = 600
)