#!/usr/bin/env Rscript
# Side-by-side volcano plots of HOMER knownResults.txt for up vs down peaks
# X-axis: only positive values (log2 fold enrichment >= 0)
# y = -log10(p-value)
#
# Usage:  Rscript atac_sidebyside_volcano_positive_xaxis.R [homer_dir] [out_dir]
# Needs:  ggplot2, ggrepel, dplyr, readr, patchwork

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(readr)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
HOMER_DIR <- if (length(args) >= 1) args[1] else
  "results/ATACseq/RNA_correlation/motifs"
OUT_DIR   <- if (length(args) >= 2) args[2] else file.path(HOMER_DIR, "volcano")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

RUNS <- c("up_vs_background", "down_vs_background")

# Human-readable description of each comparison
DESC <- list(
  up_vs_background   = c(fg = "up peaks",   bg = "background"),
  down_vs_background = c(fg = "down peaks", bg = "background")
)

Q_CUT      <- 0.05   # Benjamini q-value threshold for "significant"
FC_CUT     <- 0.25   # |log2 FC| threshold for coloring
N_LABEL    <- 10     # label top N motifs per plot
PSEUDO_PCT <- 0.1    # pseudocount added to percentages to avoid log(0)

read_known <- function(path) {
  d <- read_tsv(path, skip = 1, col_names = FALSE, show_col_types = FALSE,
                col_types = cols(.default = col_character()))
  stopifnot(ncol(d) >= 9)
  d <- d[, 1:9]
  names(d) <- c("motif", "consensus", "pvalue", "logp", "qvalue",
                "n_target", "pct_target", "n_bg", "pct_bg")
  d %>%
    mutate(
      logp       = as.numeric(logp),
      qvalue     = as.numeric(qvalue),
      pct_target = as.numeric(sub("%", "", pct_target)),
      pct_bg     = as.numeric(sub("%", "", pct_bg)),
      neglog10p  = -logp / log(10),
      log2fc     = log2((pct_target + PSEUDO_PCT) / (pct_bg + PSEUDO_PCT)),
      label      = sub("/.*$", "", motif)
    )
}

plot_sidebyside <- function() {
  # Read both datasets
  sets <- c(up_vs_background = "Up peaks", down_vs_background = "Down peaks")
  colors <- c("Up peaks" = "#C0392B", "Down peaks" = "#2471A3")
  
  d <- bind_rows(lapply(names(sets), function(run) {
    f <- file.path(HOMER_DIR, run, "knownResults.txt")
    if (!file.exists(f)) stop("Missing: ", f)
    read_known(f) %>% 
      mutate(
        set = sets[[run]],
        status = case_when(
          qvalue < Q_CUT & log2fc >= FC_CUT ~ "Enriched",
          qvalue < Q_CUT & log2fc <= -FC_CUT ~ "Depleted",
          TRUE ~ "Not significant"
        )
      )
  })) %>%
    mutate(
      set = factor(set, levels = unname(sets)),
      status = factor(status, levels = c("Enriched", "Depleted", "Not significant"))
    )
  
  # Create individual plots for each set
  plots <- lapply(names(sets), function(run) {
    set_name <- sets[[run]]
    d_subset <- d %>% filter(set == set_name)
    
    # Select labels: top N by p-value for enriched motifs
    lab <- d_subset %>%
      filter(status == "Enriched") %>%
      arrange(desc(neglog10p)) %>%
      distinct(label, .keep_all = TRUE) %>%
      slice_head(n = N_LABEL) %>%
      ungroup()
    
    ds <- DESC[[run]]
    
    p <- ggplot(d_subset, aes(log2fc, neglog10p, color = status)) +
      geom_vline(xintercept = c(-FC_CUT, FC_CUT), linetype = "dashed",
                 color = "grey60", linewidth = 0.3) +
      geom_vline(xintercept = 0, linetype = "solid",
                 color = "grey40", linewidth = 0.5) +
      geom_point(alpha = 0.6, size = 2) +
      geom_text_repel(data = lab, aes(label = label), 
                      size = 3, max.overlaps = Inf,
                      box.padding = 0.4, min.segment.length = 0, 
                      show.legend = FALSE,
                      segment.color = "grey50", segment.size = 0.3) +
      scale_color_manual(
        values = c("Enriched" = colors[set_name],
                   "Depleted" = "grey70",
                   "Not significant" = "grey90"),
        name = NULL, drop = FALSE
      ) +
      scale_x_continuous(limits = c(0, NA)) +  # Only positive x-axis
      labs(
        title = set_name,
        x = paste0("log2 fold enrichment (% ", ds["fg"], " / % ", ds["bg"], ")"),
        y = expression(-log[10](italic(P)))
      ) +
      theme_classic(base_size = 11) +
      theme(
        panel.grid.minor = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(face = "bold", hjust = 0.5)
      )
    
    return(p)
  })
  
  # Combine plots side-by-side
  combined <- plots[[1]] + plots[[2]] + 
    plot_layout(ncol = 2, guides = "collect") &
    theme(legend.position = "bottom")
  
  # Add overall title
  combined <- combined + 
    plot_annotation(
      title = "HOMER known motifs: Concordant peaks vs background (positive enrichment only)",
      theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5))
    )
  
  ggsave(file.path(OUT_DIR, "sidebyside_up_down_volcano_positive_xaxis.pdf"), 
         combined, width = 14, height = 6)
  ggsave(file.path(OUT_DIR, "sidebyside_up_down_volcano_positive_xaxis.png"), 
         combined, width = 14, height = 6, dpi = 300)
  
  message("Wrote side-by-side volcano plot (positive x-axis only)")
  invisible(combined)
}

plot_sidebyside()
message("Done. Output in ", OUT_DIR)
