# Annotate differential ATAC peaks (KO vs WT) with ChIPseeker.
# Runs fine locally in RStudio: download these from greatlakes
#   Data/Results/DESeq2/  ->  strict_more_open_in_KO.bed
#                             strict_less_open_in_KO.bed
#                             KO_vs_WT_all_peaks.csv
# into one folder and set `indir` below. Output goes to `outdir`.

# Install once if needed:
#   if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#   BiocManager::install(c("ChIPseeker", "TxDb.Hsapiens.UCSC.hg38.knownGene",
#                          "org.Hs.eg.db", "GenomicRanges"), version="3.23")
#  install.packages("ggplot2")

suppressPackageStartupMessages({
  library(ChIPseeker)
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
  library(org.Hs.eg.db)
  library(GenomicRanges)
  library(ggplot2)
})

indir  <- "results/ATACseq/DESeq2"            # <- folder with the three files above
outdir <- file.path("results/ATACseq/annotation")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

read_bed <- function(f) {
  d <- read.table(f, sep = "\t", header = FALSE, stringsAsFactors = FALSE)
  colnames(d)[1:4] <- c("chr", "start", "end", "peak")
  # BED is 0-based; GRanges is 1-based
  GRanges(d$chr, IRanges(d$start + 1, d$end), peak = d$peak)
}

# Reference = every peak that was tested in DESeq2 (background for comparison)
all_tab <- read.csv(file.path(indir, "KO_vs_WT_all_peaks.csv"))
gr_all  <- GRanges(all_tab$chr, IRanges(all_tab$start + 1, all_tab$end), peak = all_tab$peak)

peaks <- list(
  "More open in KO" = read_bed(file.path(indir, "strict_more_open_in_KO.bed")),
  "Less open in KO" = read_bed(file.path(indir, "strict_less_open_in_KO.bed")),
  "All tested peaks" = gr_all
)

# ---- annotate (promoter = +/- 3 kb of TSS) ----
anno <- lapply(peaks, annotatePeak, TxDb = txdb, tssRegion = c(-3000, 3000),
               annoDb = "org.Hs.eg.db", verbose = FALSE)

# ---- plots ----
pdf(file.path(outdir, "annotation_features.pdf"), width = 7, height = 3.5)
print(plotAnnoBar(anno))
dev.off()

pdf(file.path(outdir, "distance_to_TSS.pdf"), width = 7, height = 3.5)
print(plotDistToTSS(anno))
dev.off()

# ---- feature fractions as a table (for reporting) ----
feat <- do.call(rbind, lapply(names(anno), function(n) {
  s <- anno[[n]]@annoStat
  data.frame(set = n, feature = s$Feature, percent = round(s$Frequency, 1))
}))
write.csv(feat, file.path(outdir, "feature_percentages.csv"), row.names = FALSE)
print(feat)

# ---- annotated tables (peak -> nearest gene) ----
for (n in names(anno)[1:2]) {
  tab <- as.data.frame(anno[[n]])
  tab$start <- tab$start - 1   # back to BED coordinates
  write.csv(tab, file.path(outdir, paste0(gsub(" ", "_", n), "_annotated.csv")),
            row.names = FALSE)
}

# ---- merge gene names back into the full DESeq2 table ----
ann_all <- as.data.frame(anno[["All tested peaks"]])
keep <- ann_all[, c("peak", "SYMBOL", "annotation", "distanceToTSS")]
merged <- merge(all_tab, keep, by = "peak", all.x = TRUE)
merged <- merged[order(merged$padj), ]
write.csv(merged, file.path(outdir, "KO_vs_WT_all_peaks_annotated.csv"), row.names = FALSE)

message("Done. Outputs in: ", outdir)
