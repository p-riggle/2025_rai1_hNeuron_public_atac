# GO enrichment for genes near peaks that are more / less accessible in KO.
# Runs locally in RStudio after atac_annotate_peaks.R.
#
# Install once if needed:
#   BiocManager::install(c("clusterProfiler", "enrichplot", "org.Hs.eg.db"))

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(enrichplot)
  library(org.Hs.eg.db)
  library(ggplot2)
})

indir  <- "results/ATACseq/annotation"          # <- output folder of atac_annotate_peaks.R
outdir <- file.path("results/ATACseq/GO")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# Only link a peak to a gene if it is within this distance of the gene's TSS.
# Nearest-gene assignment is unreliable for far-away distal peaks. Set to Inf
# to use every peak. The SAME filter is applied to the background below.
max_dist <- 50000

more <- read.csv(file.path(indir, "More_open_in_KO_annotated.csv"))
less <- read.csv(file.path(indir, "Less_open_in_KO_annotated.csv"))
all_ <- read.csv(file.path(indir, "KO_vs_WT_all_peaks_annotated.csv"))

genes_of <- function(d) {
  d <- d[!is.na(d$SYMBOL) & !is.na(d$distanceToTSS) & abs(d$distanceToTSS) <= max_dist, ]
  unique(d$SYMBOL)
}

gene_lists <- list("More open in KO" = genes_of(more),
                   "Less open in KO" = genes_of(less))
# Background: genes near ANY peak that was tested in DESeq2 (same distance filter).
universe <- genes_of(all_)

message("Genes: more open = ", length(gene_lists[[1]]),
        " | less open = ", length(gene_lists[[2]]),
        " | background = ", length(universe))

# ---- GO Biological Process, each set against the shared background ----
go_one <- function(g) {
  enrichGO(gene = g, universe = universe, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
           ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05,
           minGSSize = 10, maxGSSize = 500)
}
go <- lapply(gene_lists, go_one)
go <- lapply(go, function(x) if (is.null(x) || nrow(as.data.frame(x)) == 0) x else simplify(x))

# ---- tables ----
for (n in names(go)) {
  tab <- as.data.frame(go[[n]])
  message(n, ": ", nrow(tab), " enriched GO terms (after simplify)")
  write.csv(tab, file.path(outdir, paste0("GO_BP_", gsub(" ", "_", n), ".csv")),
            row.names = FALSE)
}

# ---- plots: top terms per set + a side-by-side comparison ----
for (n in names(go)) {
  if (!is.null(go[[n]]) && nrow(as.data.frame(go[[n]])) > 0) {
    pdf(file.path(outdir, paste0("GO_BP_dotplot_", gsub(" ", "_", n), ".pdf")),
        width = 7, height = 6)
    print(dotplot(go[[n]], showCategory = 15, title = paste("GO BP:", n)))
    dev.off()
  }
}

cc <- compareCluster(geneClusters = gene_lists, fun = "enrichGO",
                     universe = universe, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                     ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 500)
if (!is.null(cc) && nrow(as.data.frame(cc)) > 0) {
  cc <- simplify(cc)
  write.csv(as.data.frame(cc), file.path(outdir, "GO_BP_compare.csv"), row.names = FALSE)
  pdf(file.path(outdir, "GO_BP_compare_dotplot.pdf"), width = 8, height = 8)
  print(dotplot(cc, showCategory = 10, title = "GO BP: more vs less open in KO"))
  dev.off()
}

message("Done. Outputs in: ", outdir)
