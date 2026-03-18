# ==========================================================
# 1️⃣ Libraries
# ==========================================================

library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

# ==========================================================
# 2️⃣ Gene list
# ==========================================================

genes_of_interest <- c(
  "SRGAP2C","SRGAP2P1","SRGAP2B","SRGAP2L","SYNGAP1","SRGAP2P2","PPP1R17",
  "POTEI","POTE2beta","POTEF","A26C1B","POTEACTIN","POTE2alpha","HARE5",
  "ZEB2","PDGFD","NOTCH2NLA","N2N","NOTCH2NL","NOTCH2NLB","NOTCH2NLC",
  "ETM6","NIID","OPDM3","ARHGAP11B","TBC1D3E","TMEM14B","CROCCP2","CROCCL1"
)

# ==========================================================
# 3️⃣ Load expression matrix
# ==========================================================

expr <- read.csv(
  "/Users/satabdi/Downloads/DEseq_rlogNorm.csv",
  row.names = 1
)

expr_df <- as.data.frame(expr)

# ==========================================================
# 4️⃣ Subset genes of interest
# ==========================================================

expr_df <- expr_df[rownames(expr_df) %in% genes_of_interest,]

expr_df <- expr_df[
  genes_of_interest[genes_of_interest %in% rownames(expr_df)],]

# ==========================================================
# 5️⃣ Safe replicate averaging
# ==========================================================

avg_reps <- function(col1,col2){
  
  cols <- intersect(c(col1,col2), colnames(expr_df))
  
  rowMeans(expr_df[,cols,drop=FALSE], na.rm = TRUE)
  
}

# ==========================================================
# 6️⃣ Day20 HET
# ==========================================================

for(i in 0:6){
  
  expr_df[[paste0("day20.het_",i)]] <-
    avg_reps(
      paste0("day20.het.rep1.het_",i),
      paste0("day20.het.rep2.het_",i)
    )
  
}

# ==========================================================
# 7️⃣ Day20 KO
# ==========================================================

for(i in 0:6){
  
  expr_df[[paste0("day20.ko_",i)]] <-
    avg_reps(
      paste0("day20.ko.rep1.ko_",i),
      paste0("day20.ko.rep2.ko_",i)
    )
  
}

# ==========================================================
# 8️⃣ Day20 WT
# ==========================================================

for(i in 0:6){
  
  expr_df[[paste0("day20.wt_",i)]] <-
    avg_reps(
      paste0("day20.wt.rep1.wt_",i),
      paste0("day20.wt.rep2.wt_",i)
    )
  
}

# ==========================================================
# 9️⃣ Day60 HET
# ==========================================================

for(i in c(0,1,2,3,5)){
  
  expr_df[[paste0("day60.het_",i)]] <-
    avg_reps(
      paste0("day60.het.rep1.het_",i),
      paste0("day60.het.rep2.het_",i)
    )
  
}

# ==========================================================
# 🔟 Day60 KO
# ==========================================================

for(i in c(0,1,2,3,5,6)){
  
  expr_df[[paste0("day60.ko_",i)]] <-
    avg_reps(
      paste0("day60.ko.rep1.ko_",i),
      paste0("day60.ko.rep2.ko_",i)
    )
  
}

# ==========================================================
# 1️⃣1️⃣ Day60 WT
# ==========================================================

for(i in c(0,1,2,3,5)){
  
  expr_df[[paste0("day60.wt_",i)]] <-
    avg_reps(
      paste0("day60.wt.rep1.wt_",i),
      paste0("day60.wt.rep2.wt_",i)
    )
  
}

# ==========================================================
# 1️⃣2️⃣ Select averaged columns
# ==========================================================

final_cols <- c(
  paste0("day20.het_",0:6),
  paste0("day20.ko_",0:6),
  paste0("day20.wt_",0:6),
  paste0("day60.het_",c(0,1,2,3,5)),
  paste0("day60.ko_",c(0,1,2,3,5,6)),
  paste0("day60.wt_",c(0,1,2,3,5))
)

heatmap_mat <- as.matrix(expr_df[,final_cols])

# ==========================================================
# 1️⃣3️⃣ Z-score normalization
# ==========================================================

z_mat <- t(scale(t(heatmap_mat)))
z_mat <- z_mat[complete.cases(z_mat),]

# ==========================================================
# 1️⃣4️⃣ Convert to long format for dot plot
# ==========================================================

dot_data <- as.data.frame(z_mat) %>%
  rownames_to_column("Gene") %>%
  pivot_longer(
    cols = -Gene,
    names_to = "Sample",
    values_to = "Zscore"
  )

dot_data$size <- abs(dot_data$Zscore)

# ==========================================================
# 1️⃣5️⃣ Dot plot
# ==========================================================
dot_plot <- ggplot(dot_data,
                   aes(x = Sample,
                       y = Gene)) +
  
  geom_point(aes(fill = Zscore,
                 size = size),
             shape = 21,          # allows border
             color = "black",     # border color
             stroke = 0.4) +      # border thickness
  
  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0,
    breaks = c(-4, -2, 0, 2, 4),
    name = "Z-score"
  ) +
  
  scale_size(
    range = c(1,6),
    name = "Z-score"
  ) +
  
  theme_minimal(base_size = 13) +
  
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 7),
    panel.grid.major = element_line(color = "grey85"),
    panel.grid.minor = element_blank()
  ) +
  
  labs(
    x = "Sample",
    y = "Gene",
    title = "Human-specific Gene Expression"
  )

dot_plot

# ==========================================================
# 1️⃣6️⃣ Save figure
# ==========================================================

ggsave(
  "/Users/satabdi/Downloads/human_specific_gene_dotplot_zscore.pdf",
  dot_plot,
  width = 14,
  height = 8
)
