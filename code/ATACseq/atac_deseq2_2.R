# ATAC-seq differential accessibility: KO vs WT (DESeq2)
# Run in the VS Code R terminal after atac_consensus_counts.sbatch finishes.

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(ggrastr)
})
# If DESeq2 is missing:
#   if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#   BiocManager::install(c("DESeq2", "apeglm"))
#   install.packages(c("ggplot2", "pheatmap", "ggrepel"))


cons    <- file.path("data/ATACseq/consensus")
outdir  <- file.path("results/ATACseq/DESeq2_2")
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
gc()  # free memory before the heaviest step
shrink_method <- "apeglm"
res_shr <- tryCatch({
  if (!requireNamespace("apeglm", quietly = TRUE)) stop("apeglm not installed")
  lfcShrink(dds, coef = "condition_KO_vs_WT", type = "apeglm")
}, error = function(e) {
  message("apeglm failed (", conditionMessage(e), "); falling back to type = 'normal'")
  shrink_method <<- "normal"
  lfcShrink(dds, coef = "condition_KO_vs_WT", type = "normal")
})
message("LFC shrinkage method used: ", shrink_method)

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
pheatmap(d, main = "Sample distances (VST)")
dev.off()

# ---- MA and volcano plots ----
# Cutoffs shown on the plots (same as the strict set above)
lfc_cut   <- 1
padj_cut  <- 0.05
n_label   <- 10          # genes to label per direction on the volcano plot
label_rank <- "lfc"     # pick labeled peaks by largest |shrunk log2FC| ("lfc") or smallest padj ("padj")
ma_require_sig <- TRUE  # TRUE: MA colors need padj < padj_cut too; FALSE: color by log2FC alone

# Gene names for labels come from the ChIPseeker annotation (atac_annotate_peaks.R).
# Copy KO_vs_WT_all_peaks_annotated.csv to the path below and re-run; until then,
# labels fall back to the peak ID.
anno_file <- file.path("results/ATACseq/annotation", "KO_vs_WT_all_peaks_annotated.csv")

plot_df <- out[!is.na(out$padj), ]
plot_df$label <- plot_df$peak
if (file.exists(anno_file)) {
  an  <- read.csv(anno_file)
  sym <- an$SYMBOL[match(plot_df$peak, an$peak)]
  plot_df$label <- ifelse(is.na(sym) | sym == "", plot_df$peak, sym)
} else {
  message("No annotation file at ", anno_file, " - volcano labels will use peak IDs.")
}

cols <- c("More open in KO" = "firebrick", "Less open in KO" = "steelblue", "NS" = "grey70")

# MA plot: red = log2FC > 1, blue = log2FC < -1, grey = otherwise
ma_df <- out[!is.na(out$baseMean) & !is.na(out$log2FC_shrunk), ]
ma_df$ma_status <- "NS"
ma_df$ma_status[ma_df$log2FC_shrunk >  lfc_cut] <- "More open in KO"
ma_df$ma_status[ma_df$log2FC_shrunk < -lfc_cut] <- "Less open in KO"
if (ma_require_sig) {
  ma_df$ma_status[is.na(ma_df$padj) | ma_df$padj >= padj_cut] <- "NS"
}
ma_df$ma_status <- factor(ma_df$ma_status, levels = c("NS", "Less open in KO", "More open in KO"))
ma_df <- ma_df[order(ma_df$ma_status), ]   # draw grey first so colored points stay visible

pdf(file.path(outdir, "MA_plot.pdf"), width = 6, height = 4.8)
print(ggplot(ma_df, aes(log10(baseMean), log2FC_shrunk, color = ma_status)) +
        geom_point_rast(size = 0.6, alpha = 0.6, dpi = 300) +
        scale_color_manual(values = cols, name = NULL) +
        coord_cartesian(ylim = c(-4, 4)) +
        theme_classic() + labs(x = "log10 mean normalized counts",
                          y = "log2 fold change (KO/WT, shrunk)", title = "KO vs WT"))
dev.off()

# Volcano plot: dashed lines at the significance cutoffs, top genes labeled
plot_df$status <- "NS"
plot_df$status[plot_df$padj < padj_cut & plot_df$log2FC_shrunk >  lfc_cut] <- "More open in KO"
plot_df$status[plot_df$padj < padj_cut & plot_df$log2FC_shrunk < -lfc_cut] <- "Less open in KO"
plot_df$neglog10p <- -log10(pmax(plot_df$padj, 1e-300))   # guard against padj = 0

top_genes <- function(d, group) {
  d <- d[d$status == group, ]
  d <- if (label_rank == "lfc") d[order(-abs(d$log2FC_shrunk)), ]
  else d[order(d$padj, -abs(d$log2FC_shrunk)), ]
  d <- d[!duplicated(d$label), ]    # one label per gene
  head(d, n_label)
}
lab_df <- rbind(top_genes(plot_df, "More open in KO"), top_genes(plot_df, "Less open in KO"))
write.csv(lab_df[, c("peak", "label", "log2FC_shrunk", "padj", "status")],
          file.path(outdir, "volcano_labeled_peaks.csv"), row.names = FALSE)

vp <- ggplot(plot_df, aes(log2FC_shrunk, neglog10p, color = status)) +
  geom_point_rast(size = 0.6, alpha = 0.6, dpi = 300) +
  geom_vline(xintercept = c(-lfc_cut, lfc_cut), linetype = "dashed", color = "grey30") +
  geom_hline(yintercept = -log10(padj_cut), linetype = "dashed", color = "grey30") +
  geom_point(data = lab_df, size = 1.8, shape = 21, color = "black",
             aes(fill = status), show.legend = FALSE) +
  scale_color_manual(values = cols, name = NULL) +
  scale_fill_manual(values = cols) +
  theme_classic() + labs(x = "log2 fold change (KO/WT, shrunk)", y = "-log10 adjusted p")
if (requireNamespace("ggrepel", quietly = TRUE)) {
  vp <- vp + ggrepel::geom_text_repel(data = lab_df, aes(label = label), color = "black",
                                      size = 3, max.overlaps = Inf, box.padding = 0.5,
                                      min.segment.length = 0)
} else {
  message("ggrepel not installed - labels may overlap. install.packages('ggrepel')")
  vp <- vp + geom_text(data = lab_df, aes(label = label), color = "black",
                       size = 3, vjust = -0.8, check_overlap = TRUE)
}
pdf(file.path(outdir, "volcano.pdf"), width = 6.5, height = 6)
print(vp)
dev.off()

saveRDS(dds, file.path(outdir, "dds.rds"))
message("Done. Outputs in: ", outdir)