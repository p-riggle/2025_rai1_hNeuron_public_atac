#!/usr/bin/env Rscript
# Volcano plots of HOMER knownResults.txt for the concordant-peak runs against
# background: up_vs_background and down_vs_background (individually and combined)
#
# x = log2 fold enrichment (% target with motif / % background with motif)
# y = -log10(p-value)   (computed from HOMER's natural-log "Log P-value" column,
#                        so very small p-values don't underflow to 0)
#
# Usage:  Rscript homer_known_volcano.R [homer_dir] [out_dir]
# Needs:  ggplot2, ggrepel, dplyr, readr   (install.packages(c("ggplot2","ggrepel","dplyr","readr")))
# need to download the HOMER knownResults.txt files from Greatlakes to motifs/down_vs_background and up_vs_background in Rstudio before running.

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
HOMER_DIR <- if (length(args) >= 1) args[1] else
  "results/ATACseq/RNA_correlation/motifs"
OUT_DIR   <- if (length(args) >= 2) args[2] else file.path(HOMER_DIR, "volcano")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

RUNS <- c("up_vs_background", "down_vs_background")

# Human-readable description of each comparison (used in titles / axis labels)
DESC <- list(
  up_vs_background   = c(fg = "up peaks",   bg = "background"),
  down_vs_background = c(fg = "down peaks", bg = "background")
)

Q_CUT      <- 0.05   # Benjamini q-value threshold for "significant"
FC_CUT     <- 0.25   # |log2 FC| threshold for coloring (set 0 to ignore)
N_LABEL    <- 10     # label top N motifs per side (by p-value)
PSEUDO_PCT <- 0.1    # pseudocount added to percentages to avoid log(0)

read_known <- function(path) {
  # Read by position: HOMER header text varies slightly between versions
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
      # "Name(Family)/Cell-Factor-ChIP-Seq/Source" -> "Name(Family)"
      label      = sub("/.*$", "", motif)
    )
}

plot_run <- function(run) {
  f <- file.path(HOMER_DIR, run, "knownResults.txt")
  if (!file.exists(f)) { warning("Missing: ", f); return(invisible(NULL)) }
  
  d <- read_known(f) %>%
    mutate(
      status = case_when(
        qvalue < Q_CUT & log2fc >=  FC_CUT ~ "Enriched in foreground",
        qvalue < Q_CUT & log2fc <= -FC_CUT ~ "Depleted in foreground",
        TRUE                               ~ "Not significant"
      ),
      status = factor(status, levels = c("Enriched in foreground",
                                         "Depleted in foreground",
                                         "Not significant"))
    )
  
  # Label the most significant motifs on each side; de-duplicate near-identical
  # motif names (HOMER has many ChIP-seq replicates of the same factor)
  lab <- d %>%
    filter(status != "Not significant") %>%
    arrange(desc(neglog10p)) %>%
    distinct(label, .keep_all = TRUE) %>%
    group_by(status) %>%
    slice_head(n = N_LABEL) %>%
    ungroup()
  
  ds <- DESC[[run]]
  p <- ggplot(d, aes(log2fc, neglog10p, color = status)) +
    geom_vline(xintercept = c(-FC_CUT, FC_CUT), linetype = "dashed",
               color = "grey60", linewidth = 0.3) +
    geom_point(alpha = 0.7, size = 1.8) +
    geom_text_repel(data = lab, aes(label = label), size = 3, max.overlaps = Inf,
                    box.padding = 0.4, min.segment.length = 0, show.legend = FALSE,
                    segment.color = "grey50", segment.size = 0.3) +
    scale_color_manual(values = c("Enriched in foreground" = "#C0392B",
                                  "Depleted in foreground" = "#2471A3",
                                  "Not significant"        = "grey70"),
                       name = NULL, drop = FALSE) +
    labs(
      title    = paste0("HOMER known motifs: ", ds["fg"], " vs ", ds["bg"]),
      subtitle = paste0("Colored: q < ", Q_CUT, " & |log2 FC| >= ", FC_CUT,
                        "   (", sum(d$status != "Not significant"), " of ", nrow(d), " motifs)"),
      x = paste0("log2 fold enrichment (% ", ds["fg"], " / % ", ds["bg"], ")"),
      y = expression(-log[10](italic(P)))
    ) +
    theme_classic(base_size = 12) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")
  
  ggsave(file.path(OUT_DIR, paste0(run, "_volcano.pdf")), p, width = 7, height = 6)
  ggsave(file.path(OUT_DIR, paste0(run, "_volcano.png")), p, width = 7, height = 6, dpi = 300)
  
  # Also save the parsed table for convenience
  write_csv(d %>% select(-label) %>% arrange(desc(neglog10p)),
            file.path(OUT_DIR, paste0(run, "_known_motifs_parsed.csv")))
  message("Wrote plots for ", run)
  invisible(p)
}

# ---- Combined plot: up-vs-background and down-vs-background on one volcano ----
# The two input files have the same name (knownResults.txt) but live in different
# directories (homer/up_vs_background/ and homer/down_vs_background/). The set label
# is taken from the directory name, so each point is tagged "Up peaks" or "Down peaks".
plot_combined <- function() {
  sets <- c(up_vs_background = "Up peaks", down_vs_background = "Down peaks")
  d <- bind_rows(lapply(names(sets), function(run) {
    f <- file.path(HOMER_DIR, run, "knownResults.txt")
    if (!file.exists(f)) stop("Missing: ", f)
    read_known(f) %>% mutate(set = sets[[run]])
  })) %>%
    mutate(
      set  = factor(set, levels = unname(sets)),
      sig  = qvalue < Q_CUT & abs(log2fc) >= FC_CUT,
      grp  = ifelse(sig, as.character(set), "Not significant"),
      grp  = factor(grp, levels = c(unname(sets), "Not significant"))
    )
  
  lab <- d %>%
    filter(sig) %>%
    arrange(desc(neglog10p)) %>%
    distinct(set, label, .keep_all = TRUE) %>%
    group_by(set) %>%
    slice_head(n = N_LABEL) %>%
    ungroup()
  
  p <- ggplot(d, aes(log2fc, neglog10p, color = grp)) +
    geom_vline(xintercept = c(-FC_CUT, FC_CUT), linetype = "dashed",
               color = "grey60", linewidth = 0.3) +
    geom_point(data = filter(d, !sig), alpha = 0.5, size = 1.5) +
    geom_point(data = filter(d, sig),  alpha = 0.85, size = 2) +
    geom_text_repel(data = lab, aes(label = label), size = 3, max.overlaps = Inf,
                    box.padding = 0.4, min.segment.length = 0, show.legend = FALSE,
                    segment.color = "grey50", segment.size = 0.3) +
    scale_color_manual(values = c("Up peaks" = "#C0392B", "Down peaks" = "#2471A3",
                                  "Not significant" = "grey75"), name = NULL, drop = FALSE) +
    labs(
      title    = "HOMER known motifs: concordant up and down peaks vs background",
      subtitle = paste0("Colored: q < ", Q_CUT, " & |log2 FC| >= ", FC_CUT,
                        "; each motif appears once per peak set"),
      x = "log2 fold enrichment (% foreground / % background)",
      y = expression(-log[10](italic(P)))
    ) +
    theme_classic(base_size = 12) +
    theme(panel.grid.minor = element_blank(), legend.position = "bottom")
  
  ggsave(file.path(OUT_DIR, "up_and_down_vs_background_volcano.pdf"), p, width = 8, height = 6.5)
  ggsave(file.path(OUT_DIR, "up_and_down_vs_background_volcano.png"), p, width = 8, height = 6.5, dpi = 300)
  message("Wrote combined up/down vs background plot")
  invisible(p)
}

invisible(lapply(RUNS, plot_run))   # individual up / down vs background plots
plot_combined()                     # plus the combined up + down vs background plot
message("Done. Output in ", OUT_DIR)