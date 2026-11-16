# ATAC-seq differential accessibility: KO vs WT (DESeq2)
# Run in the VS Code R terminal after atac_consensus_counts.sbatch finishes.

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})
# If DESeq2 is missing:
#   if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#   BiocManager::install(c("DESeq2", "apeglm"))
#   install.packages(c("ggplot2", "pheatmap"))

cons    <- file.path("data/ATACseq/consensus")
outdir  <- file.path("results/ATACseq/DESeq2")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

samples <- c("WT1", "WT2", "KO1", "KO2")

# ---- Load featureCounts output ----
fc <- read.delim(file.path(cons, "consensus_counts.txt"), comment.char = "#",
                 check.names = FALSE)
counts <- as.matrix(fc[, 7:10])
colnames(counts) <- samples
rownames(counts) <- fc$Geneid

peak_info <- data.frame(peak = fc$Geneid, chr = fc$Chr,
                        start = fc$Start - 1, end = fc$End)  # back to BED coords

coldata <- data.frame(
  row.names = samples,
  condition = factor(c("WT", "WT", "KO", "KO"), levels = c("WT", "KO"))  # WT = reference
)

# ---- FRiP from featureCounts summary ----
sm <- read.delim(file.path(cons, "consensus_counts.txt.summary"), check.names = FALSE)
tot      <- colSums(sm[, -1])
assigned <- unlist(sm[sm$Status == "Assigned", -1])
frip <- data.frame(sample = samples, fragments_in_peaks = as.numeric(assigned),
                   total_fragments = as.numeric(tot),
                   FRiP = as.numeric(assigned) / as.numeric(tot))
write.csv(frip, file.path(outdir, "FRiP.csv"), row.names = FALSE)
print(frip)

# ---- DESeq2 ----
dds <- DESeqDataSetFromMatrix(countData = counts, colData = coldata, design = ~ condition)
dds <- dds[rowSums(counts(dds) >= 10) >= 2, ]   # drop very low-count peaks
dds <- DESeq(dds)
message("Peaks tested: ", nrow(dds))

res <- results(dds, contrast = c("condition", "KO", "WT"), alpha = 0.05)
if (requireNamespace("apeglm", quietly = TRUE)) {
  res_shr <- lfcShrink(dds, coef = "condition_KO_vs_WT", type = "apeglm")
} else {
  res_shr <- lfcShrink(dds, coef = "condition_KO_vs_WT", type = "normal")
}

out <- cbind(peak_info[match(rownames(res), peak_info$peak), ],
             as.data.frame(res),
             log2FC_shrunk = res_shr$log2FoldChange)
out <- out[order(out$padj), ]
write.csv(out, file.path(outdir, "KO_vs_WT_all_peaks.csv"), row.names = FALSE)

sig <- subset(out, !is.na(padj) & padj < 0.05)
write.csv(sig, file.path(outdir, "KO_vs_WT_significant_padj0.05.csv"), row.names = FALSE)
# (BED files for downstream work come from the strict set written below)

message("Significant (padj<0.05): ", nrow(sig),
        "  | more open in KO: ", sum(sig$log2FoldChange > 0),
        "  | less open in KO: ", sum(sig$log2FoldChange < 0))

# ---- Strict differential set: padj < 0.05 AND |shrunk log2FC| > 1 ----
# Uses the shrunken fold change (less inflated for low-count peaks with n=2).
# These are the sets to use for annotation, motifs and heatmaps.
strict <- subset(out, !is.na(padj) & padj < 0.05 & abs(log2FC_shrunk) > 1)
message("Strict (padj<0.05, |shrunk log2FC|>1): ", nrow(strict),
        "  | more open in KO: ", sum(strict$log2FC_shrunk > 0),
        "  | less open in KO: ", sum(strict$log2FC_shrunk < 0))
write.csv(strict, file.path(outdir, "KO_vs_WT_strict_padj0.05_lfc1.csv"), row.names = FALSE)
write.table(subset(strict, log2FC_shrunk > 0)[, c("chr", "start", "end", "peak")],
            file.path(outdir, "strict_more_open_in_KO.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(subset(strict, log2FC_shrunk < 0)[, c("chr", "start", "end", "peak")],
            file.path(outdir, "strict_less_open_in_KO.bed"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# ---- Plots ----
vsd <- vst(dds, blind = TRUE)

# ---- QC summary: exact PCA variance, size factors, global shift, peaks tested ----
# plotPCA only uses the top 500 variable peaks, so compute variance on all peaks.
pca_all  <- prcomp(t(assay(vsd)))
pc_var   <- round(100 * pca_all$sdev^2 / sum(pca_all$sdev^2), 1)
qc_lines <- c(
  paste("Peaks tested:", nrow(dds)),
  paste("PCA variance explained (all peaks, PC1-PC4, %):", paste(pc_var, collapse = " / ")),
  paste("Size factors:", paste(names(sizeFactors(dds)), round(sizeFactors(dds), 3),
                               sep = "=", collapse = ", ")),
  paste("Median log2FC (should be ~0 if no global shift):",
        round(median(res$log2FoldChange, na.rm = TRUE), 4)),
  paste("Significant padj<0.05:", nrow(sig)),
  paste("Strict padj<0.05 & |shrunk log2FC|>1:", nrow(strict),
        "(more open in KO:", sum(strict$log2FC_shrunk > 0),
        "| less open in KO:", sum(strict$log2FC_shrunk < 0), ")")
)
writeLines(qc_lines, file.path(outdir, "qc_summary.txt"))
message(paste(qc_lines, collapse = "\n"))

pdf(file.path(outdir, "PCA.pdf"), width = 5, height = 4.5)
pca <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
pv  <- round(100 * attr(pca, "percentVar"))
print(ggplot(pca, aes(PC1, PC2, color = condition, label = name)) +
        geom_point(size = 4) + geom_text(vjust = -1, show.legend = FALSE) +
        xlab(paste0("PC1: ", pv[1], "% variance")) +
        ylab(paste0("PC2: ", pv[2], "% variance")) +
        theme_bw() + coord_cartesian(clip = "off"))
dev.off()

pdf(file.path(outdir, "sample_distance_heatmap.pdf"), width = 5, height = 4.5)
d <- as.matrix(dist(t(assay(vsd))))
heatmap(d, main = "Sample distances (VST)")
dev.off()

pdf(file.path(outdir, "MA_plot.pdf"), width = 5.5, height = 4.5)
plotMA(res_shr, ylim = c(-4, 4), main = "KO vs WT")
dev.off()

pdf(file.path(outdir, "volcano.pdf"), width = 5.5, height = 5)
vd <- out[!is.na(out$padj), ]
vd$status <- ifelse(vd$padj < 0.05 & vd$log2FoldChange > 0, "More open in KO",
                    ifelse(vd$padj < 0.05 & vd$log2FoldChange < 0, "Less open in KO", "NS"))
print(ggplot(vd, aes(log2FoldChange, -log10(padj), color = status)) +
        geom_point(size = 0.6, alpha = 0.6) +
        scale_color_manual(values = c("More open in KO" = "firebrick",
                                      "Less open in KO" = "steelblue", "NS" = "grey70")) +
        theme_bw() + labs(x = "log2 fold change (KO/WT)", y = "-log10 adjusted p"))
dev.off()

saveRDS(dds, file.path(outdir, "dds.rds"))
message("Done. Outputs in: ", outdir)

