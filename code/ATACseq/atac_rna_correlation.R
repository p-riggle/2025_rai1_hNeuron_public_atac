# ATAC vs RNA-seq (KO vs WT): does chromatin accessibility track gene expression?
# Runs locally in RStudio. Inputs:
#   1) KO_vs_WT_all_peaks_annotated.csv  (from atac_annotate_peaks.R; has SYMBOL,
#      distanceToTSS, log2FC_shrunk, padj for every tested ATAC peak)
#   2) your RNA-seq DESeq2 results table (Ensembl IDs + log2FC)
#
# Outputs (in `outdir`), for each peak-to-gene distance scenario:
#   ATAC_vs_RNA_scatter_<scenario>.pdf   ATAC log2FC vs RNA log2FC, one point per gene
#   ATAC_vs_RNA_boxplot_<scenario>.pdf   RNA log2FC for genes near more / less / unchanged peaks
#   ATAC_RNA_merged_<scenario>.csv       the gene-level table behind the plots
#   ATAC_vs_RNA_summary.csv              correlations and concordance for all scenarios
#
# Install once if needed:
#   install.packages(c("ggplot2", "ggrepel"))
#   BiocManager::install(c("org.Hs.eg.db", "AnnotationDbi"))   # only to map Ensembl -> symbol

suppressPackageStartupMessages(library(ggplot2))
has_repel <- requireNamespace("ggrepel", quietly = TRUE)

# ======================= settings =======================
atac_file <- "results/ATACseq/annotation/KO_vs_WT_all_peaks_annotated.csv"
rna_file  <- "results/ineuron/deseq2/rstbl_d7_wt.vs.ko.txt"      # <- your RNA-seq table
outdir    <- "results/ATACseq/RNA_correlation"

# RNA table columns (check with names(read.csv(rna_file)) if unsure)
#   (skips the org.Hs.eg.db lookup)
rna_lfc_col      <- "log2FoldChange"
rna_padj_col     <- "padj"             # NULL if the table has no adjusted p-value
rna_basemean_col <- "baseMean"         # NULL if the table has no mean expression
flip_rna         <- FALSE            # TRUE if your RNA table is WT vs KO (ATAC here is KO vs WT)
rna_symbol_col <- NULL #if no gene names

# Thresholds
atac_padj_cut <- 0.05;  atac_lfc_cut <- 1     # "changed" ATAC peak (same as the strict set)
rna_padj_cut  <- 0.05;  rna_lfc_cut  <- 0.5   # "changed" RNA gene
min_basemean  <- 10                            # drop barely expressed genes (needs baseMean column)

# One gene can have several nearby peaks. "strongest": use a significant peak if the gene
# has one, otherwise the peak with the largest |log2FC|. "closest": peak nearest the TSS.
collapse <- "strongest"

# Peak-to-gene distance scenarios (bp from the nearest TSS)
scenarios <- c("promoter_1kb" = 1000, "within_3kb" = 3000, "all_nearest_gene" = Inf)
n_label <- 10          # genes labeled on the scatter plot
# ========================================================

dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

read_table <- function(f) {
  f <- path.expand(f)
  if (grepl("\\.(tsv|txt)$", f, ignore.case = TRUE)) read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
  else read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
}
need <- function(df, cols, what) {
  miss <- setdiff(cols, names(df))
  if (length(miss)) stop(what, " is missing column(s): ", paste(miss, collapse = ", "),
                         "\n  columns found: ", paste(names(df), collapse = ", "), call. = FALSE)
}

# ---------------- RNA-seq table ----------------
rna <- read_table(rna_file)
need(rna, c(rna_lfc_col,
            rna_padj_col, if (!is.null(rna_basemean_col)) rna_basemean_col), "RNA table")
has_padj <- !is.null(rna_padj_col)

rna$rna_lfc  <- as.numeric(rna[[rna_lfc_col]]) * (if (flip_rna) -1 else 1)
rna$rna_padj <- if (has_padj) as.numeric(rna[[rna_padj_col]]) else NA_real_
rna$rna_base <- if (!is.null(rna_basemean_col)) as.numeric(rna[[rna_basemean_col]]) else NA_real_
n_rna_start  <- nrow(rna)

rna$gene_id <- rownames(rna)
rna_id_col       <- "gene_id" 


if (!is.null(rna_symbol_col)) {
  rna$SYMBOL <- as.character(rna[[rna_symbol_col]])
} else {
  rna$ensembl <- sub("\\..*$", "", as.character(rna[[rna_id_col]]))
  if (!any(grepl("^ENSG", rna$ensembl)))
    warning("Column '", rna_id_col, "' does not look like human Ensembl gene IDs (ENSG...).")
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
    stop("Install org.Hs.eg.db (BiocManager::install('org.Hs.eg.db')) or set rna_symbol_col.")
  sym <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db,
                                                keys = unique(rna$ensembl), keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first"))
  rna$SYMBOL <- unname(sym[rna$ensembl])
}
n_mapped <- sum(!is.na(rna$SYMBOL) & rna$SYMBOL != "")
rna <- rna[!is.na(rna$SYMBOL) & rna$SYMBOL != "" & !is.na(rna$rna_lfc), ]
if (!all(is.na(rna$rna_base))) rna <- rna[!is.na(rna$rna_base) & rna$rna_base >= min_basemean, ]
# one row per symbol: keep the most highly expressed Ensembl ID
rna <- rna[if (all(is.na(rna$rna_base))) seq_len(nrow(rna)) else order(-rna$rna_base), ]
rna <- rna[!duplicated(rna$SYMBOL), ]
message("RNA: ", n_rna_start, " rows -> ", n_mapped, " with a gene symbol -> ",
        nrow(rna), " genes kept (non-NA log2FC, baseMean >= ", min_basemean, ")")

# ---------------- ATAC table ----------------
atac <- read_table(atac_file)
need(atac, c("peak", "SYMBOL", "distanceToTSS", "log2FC_shrunk", "padj"), "ATAC table")
atac <- atac[!is.na(atac$SYMBOL) & atac$SYMBOL != "" &
               !is.na(atac$distanceToTSS) & !is.na(atac$log2FC_shrunk), ]
atac$strict <- !is.na(atac$padj) & atac$padj < atac_padj_cut & abs(atac$log2FC_shrunk) > atac_lfc_cut
message("ATAC: ", nrow(atac), " peaks with a gene assignment, ", sum(atac$strict), " strict differential")

collapse_peaks <- function(d, how) {
  o <- if (how == "closest") order(abs(d$distanceToTSS)) else order(-d$strict, -abs(d$log2FC_shrunk))
  d <- d[o, ]
  d[!duplicated(d$SYMBOL), ]
}
fmt_p <- function(p) ifelse(is.na(p), "NA", format.pval(p, digits = 2, eps = 1e-16))

# ---------------- one scenario ----------------
run_one <- function(tag, max_dist) {
  d <- atac[abs(atac$distanceToTSS) <= max_dist, ]
  g <- collapse_peaks(d, collapse)
  m <- merge(g[, c("SYMBOL", "peak", "distanceToTSS", "log2FC_shrunk", "padj", "strict")],
             rna[, c("SYMBOL", "rna_lfc", "rna_padj", "rna_base")], by = "SYMBOL")
  names(m)[names(m) == "log2FC_shrunk"] <- "atac_lfc"
  names(m)[names(m) == "padj"]          <- "atac_padj"
  if (nrow(m) < 10) { message("[", tag, "] only ", nrow(m), " matched genes - skipped"); return(NULL) }
  
  m$rna_sig <- abs(m$rna_lfc) > rna_lfc_cut & (if (has_padj) (!is.na(m$rna_padj) & m$rna_padj < rna_padj_cut) else TRUE)
  both    <- m$strict & m$rna_sig
  concord <- sign(m$atac_lfc) == sign(m$rna_lfc)
  m$cat      <- ifelse(both & concord & m$atac_lfc > 0, "Concordant up", 
                       ifelse(both & concord & m$atac_lfc < 0, "Concordant down", 
                              ifelse(both, "Discordant", "Other")))
  m$atac_cat <- factor(ifelse(m$strict & m$atac_lfc > 0, "More open in KO",
                              ifelse(m$strict & m$atac_lfc < 0, "Less open in KO", "No sig. change")),
                       levels = c("Less open in KO", "No sig. change", "More open in KO"))
  
  sp <- suppressWarnings(cor.test(m$atac_lfc, m$rna_lfc, method = "spearman", exact = FALSE))
  pe <- cor.test(m$atac_lfc, m$rna_lfc, method = "pearson")
  n_both <- sum(both); n_con <- sum(both & concord)
  p_binom <- if (n_both > 0) binom.test(n_con, n_both, 0.5)$p.value else NA
  grp <- split(m$rna_lfc, m$atac_cat)
  wt <- function(a, b) if (length(grp[[a]]) > 2 && length(grp[[b]]) > 2)
    suppressWarnings(wilcox.test(grp[[a]], grp[[b]])$p.value) else NA
  p_ml <- wt("More open in KO", "Less open in KO")
  p_mn <- wt("More open in KO", "No sig. change")
  p_ln <- wt("Less open in KO", "No sig. change")
  
  title_txt <- sprintf("%s (peak within %s of TSS)", gsub("_", " ", tag),
                       if (is.finite(max_dist)) paste0(max_dist / 1000, " kb") else "any distance")
  
  # ---- scatter ----
  cols <- c("Concordant up" = "firebrick", "Concordant down" = "steelblue", Discordant = "#000000")
  cat_n <- sprintf("%s (n = %d)", names(cols), sapply(names(cols), function(k) sum(m$cat == k)))
  lab <- m[both, ]
  lab <- head(lab[order(-abs(lab$atac_lfc * lab$rna_lfc)), ], n_label)
  rng <- function(v, extra) range(c(quantile(v, c(0.001, 0.999)), extra))
  sc <- ggplot(m, aes(atac_lfc, rna_lfc)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(data = m[m$cat == "Other", ], color = "grey75", size = 0.6, alpha = 0.5) +
    geom_smooth(method = "lm", formula = y ~ x, color = "black", linewidth = 0.5, se = TRUE) +
    geom_point(data = m[m$cat != "Other", ], aes(color = cat), size = 1.5, alpha = 0.85) +
    scale_color_manual(values = cols, breaks = names(cols), labels = cat_n, name = NULL) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.05, vjust = 1.3, size = 3.4,
             label = sprintf("Spearman rho = %.2f (p %s)\nPearson r = %.2f\nn = %s genes",
                             unname(sp$estimate), fmt_p(sp$p.value), unname(pe$estimate),
                             formatC(nrow(m), big.mark = ",", format = "d"))) +
    coord_cartesian(xlim = rng(m$atac_lfc, lab$atac_lfc), ylim = rng(m$rna_lfc, lab$rna_lfc)) +
    labs(x = "ATAC log2FC (KO/WT, shrunk), peak nearest each gene",
         y = "RNA-seq log2FC (KO/WT)", title = title_txt,
         subtitle = sprintf("Both changed: %d genes, %.0f%% concordant (sign test p %s)",
                            n_both, ifelse(n_both > 0, 100 * sum(m$cat %in% c("Concordant up", "Concordant down")) / n_both, NA), fmt_p(p_binom))) +
    theme_classic()
  if (has_repel && nrow(lab) > 0)
    sc <- sc + ggrepel::geom_text_repel(data = lab, aes(label = SYMBOL), size = 3,
                                        max.overlaps = Inf, min.segment.length = 0)
  ggsave(file.path(outdir, paste0("ATAC_vs_RNA_scatter_", tag, ".pdf")), sc, width = 6.8, height = 6)
  
  # ---- boxplot: RNA log2FC by ATAC category ----
  n_lab <- table(m$atac_cat)
  xlabs <- setNames(sprintf("%s\n(n = %s)", levels(m$atac_cat),
                            formatC(as.integer(n_lab), big.mark = ",", format = "d")), levels(m$atac_cat))
  bp <- ggplot(m, aes(atac_cat, rna_lfc, fill = atac_cat)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    geom_violin(alpha = 0.35, color = NA, scale = "width") +
    geom_boxplot(width = 0.25, outlier.shape = NA, fill = "white") +
    scale_fill_manual(values = c("Less open in KO" = "steelblue", "No sig. change" = "grey70",
                                 "More open in KO" = "firebrick"), guide = "none") +
    scale_x_discrete(labels = xlabs) +
    coord_cartesian(ylim = quantile(m$rna_lfc, c(0.01, 0.99))) +
    labs(x = NULL, y = "RNA-seq log2FC (KO/WT)", title = title_txt,
         subtitle = sprintf("Wilcoxon p: more vs less %s | more vs unchanged %s | less vs unchanged %s",
                            fmt_p(p_ml), fmt_p(p_mn), fmt_p(p_ln))) +
    theme_classic()
  ggsave(file.path(outdir, paste0("ATAC_vs_RNA_boxplot_", tag, ".pdf")), bp, width = 5.5, height = 5)
  
  write.csv(m[order(-abs(m$atac_lfc)), ], file.path(outdir, paste0("ATAC_RNA_merged_", tag, ".csv")),
            row.names = FALSE)
  
  data.frame(scenario = tag, max_dist_bp = max_dist, n_genes = nrow(m),
             spearman_rho = unname(sp$estimate), spearman_p = sp$p.value,
             pearson_r = unname(pe$estimate), pearson_p = pe$p.value,
             n_both_changed = n_both, n_concordant = n_con,
             pct_concordant = ifelse(n_both > 0, 100 * n_con / n_both, NA), sign_test_p = p_binom,
             median_rna_lfc_more_open = median(grp[["More open in KO"]]),
             median_rna_lfc_less_open = median(grp[["Less open in KO"]]),
             median_rna_lfc_unchanged = median(grp[["No sig. change"]]),
             wilcox_p_more_vs_less = p_ml, row.names = NULL)
}

summ <- do.call(rbind, Map(run_one, names(scenarios), scenarios))
write.csv(summ, file.path(outdir, "ATAC_vs_RNA_summary.csv"), row.names = FALSE)
print(format(summ, digits = 3), row.names = FALSE)
message("Done. Outputs in: ", outdir)