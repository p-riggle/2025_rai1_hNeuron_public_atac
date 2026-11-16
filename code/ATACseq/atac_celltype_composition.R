# Cell-composition check: do whole cell-type marker sets shift in KO vs WT?
# Runs locally in RStudio. Uses the RNA-seq DESeq2 table and the annotated ATAC table.
#
# Logic: if KO and WT cultures differ in WHICH cells they contain (e.g. more mature
# neurons, fewer mesenchymal-like cells), the markers of each cell type should move
# together as a set, in both RNA and promoter accessibility. If only a few genes move,
# that points to gene-specific regulation instead. This is a marker-based check, not a
# deconvolution; confirm any composition difference with imaging or single-cell data.
#
# Outputs (in `outdir`):
#   marker_sets_summary.pdf   RNA and promoter-ATAC log2FC per marker set
#   marker_genes_rna.pdf      every marker gene's RNA log2FC, grouped by set
#   marker_sets_stats.csv     per-set statistics
#   marker_genes_table.csv    gene-level values behind the plots
#
# Needs only ggplot2 (+ org.Hs.eg.db if your RNA table has Ensembl IDs but no symbols).

suppressPackageStartupMessages(library(ggplot2))

# ======================= settings =======================
rna_file  <- "results/ineuron/deseq2/rstbl_d7_wt.vs.ko.txt"
atac_file <- "results/ATACseq/annotation/KO_vs_WT_all_peaks_annotated.csv"
outdir    <- "results/ATACseq/celltype_composition"

# RNA table columns (same as atac_rna_correlation.R)
rna_id_col       <- "gene_id"
rna_symbol_col   <- NULL
rna_lfc_col      <- "log2FoldChange"
rna_padj_col     <- "padj"
rna_basemean_col <- "baseMean"
flip_rna         <- FALSE        # TRUE if your RNA table is WT vs KO
min_basemean     <- 10           # genes below this count as "not expressed" for the statistics
atac_max_dist    <- 3000         # ATAC: use the peak closest to the TSS, within this distance (bp)

# Marker sets. Edit freely, these are a starting point from commonly used human markers.
# A gene is counted once (in the first set that lists it).
markers <- list(
  "Mature neuron"   = c("MAP2", "TUBB3", "RBFOX3", "SYP", "SNAP25", "SYT1", "STMN2", "NEFL", "NEFM",
                        "NEFH", "GAP43", "ELAVL4", "SLC17A6", "SLC17A7", "CAMK2A", "GRIN1", "GRIA2",
                        "DLG4", "NRXN1", "SCN2A"),
  "Immature neuron" = c("DCX", "NEUROD1", "NEUROD2", "NEUROD6", "MYT1L", "NHLH1", "NHLH2", "SOX11",
                        "NEUROG2", "TUBB2B"),
  "Neural progenitor" = c("NES", "SOX2", "PAX6", "SOX1", "HES1", "HES5", "FABP7", "ASCL1", "HMGA2",
                          "NOTCH1"),
  "Pluripotency"    = c("POU5F1", "NANOG", "LIN28A", "DPPA4", "DPPA2", "ZFP42", "DNMT3B",
                        "UTF1", "PRDM14"),
  "Mesenchymal / ECM" = c("COL1A1", "COL1A2", "COL3A1", "COL5A1", "FN1", "VIM", "THY1", "DCN", "LUM",
                          "FBN1", "ACTA2", "TAGLN", "S100A4", "PDGFRB", "CDH11", "SPARC"),
  "Glia"            = c("GFAP", "S100B", "AQP4", "ALDH1L1", "SLC1A2", "SOX9", "OLIG1", "OLIG2", "PLP1",
                        "MBP"),
  "Proliferation"   = c("MKI67", "TOP2A", "PCNA", "MCM2", "CCNB1", "CDK1", "AURKB", "CENPF", "BIRC5",
                        "E2F1")
)
rna_padj_cut <- 0.05
# ========================================================

rna_file <- path.expand(rna_file); atac_file <- path.expand(atac_file)
outdir <- path.expand(outdir); dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

read_table <- function(f) {
  if (grepl("\\.(tsv|txt)$", f, ignore.case = TRUE)) read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
  else read.csv(f, check.names = FALSE, stringsAsFactors = FALSE)
}
need <- function(df, cols, what) {
  miss <- setdiff(cols, names(df))
  if (length(miss)) stop(what, " is missing column(s): ", paste(miss, collapse = ", "),
                         "\n  columns found: ", paste(names(df), collapse = ", "), call. = FALSE)
}
fmt_p <- function(p) ifelse(is.na(p), "NA", format.pval(p, digits = 2, eps = 1e-16))

# ---------------- RNA ----------------
rna <- read_table(rna_file)

rna$gene_id <- rownames(rna)
rna_id_col       <- "gene_id" 

need(rna, c(rna_lfc_col, if (is.null(rna_symbol_col)) rna_id_col else rna_symbol_col,
            rna_padj_col, rna_basemean_col), "RNA table")
rna$rna_lfc  <- as.numeric(rna[[rna_lfc_col]]) * (if (flip_rna) -1 else 1)
rna$rna_padj <- if (!is.null(rna_padj_col)) as.numeric(rna[[rna_padj_col]]) else NA_real_
rna$rna_base <- if (!is.null(rna_basemean_col)) as.numeric(rna[[rna_basemean_col]]) else NA_real_
if (!is.null(rna_symbol_col)) {
  rna$SYMBOL <- as.character(rna[[rna_symbol_col]])
} else {
  rna$ensembl <- sub("\\..*$", "", as.character(rna[[rna_id_col]]))
  if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
    stop("Install org.Hs.eg.db (BiocManager::install('org.Hs.eg.db')) or set rna_symbol_col.")
  sym <- suppressMessages(AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = unique(rna$ensembl),
                                                keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first"))
  rna$SYMBOL <- unname(sym[rna$ensembl])
}
rna <- rna[!is.na(rna$SYMBOL) & rna$SYMBOL != "" & !is.na(rna$rna_lfc), ]
rna <- rna[order(-rna$rna_base), ]; rna <- rna[!duplicated(rna$SYMBOL), ]
rna$expressed <- !is.na(rna$rna_base) & rna$rna_base >= min_basemean
message("RNA: ", nrow(rna), " genes with a symbol and a log2FC (", sum(rna$expressed),
        " with baseMean >= ", min_basemean, ")")

# ---------------- ATAC: one promoter-proximal peak per gene ----------------
atac <- read_table(atac_file)
need(atac, c("SYMBOL", "distanceToTSS", "log2FC_shrunk", "padj"), "ATAC table")
atac <- atac[!is.na(atac$SYMBOL) & atac$SYMBOL != "" & !is.na(atac$distanceToTSS) &
               !is.na(atac$log2FC_shrunk) & abs(atac$distanceToTSS) <= atac_max_dist, ]
atac <- atac[order(abs(atac$distanceToTSS)), ]; atac <- atac[!duplicated(atac$SYMBOL), ]
atac$atac_sig <- !is.na(atac$padj) & atac$padj < 0.05
message("ATAC: ", nrow(atac), " genes with a peak within ", atac_max_dist, " bp of the TSS")

# ---------------- marker table ----------------
mk <- do.call(rbind, lapply(names(markers), function(s)
  data.frame(set = s, SYMBOL = markers[[s]], stringsAsFactors = FALSE)))
mk <- mk[!duplicated(mk$SYMBOL), ]
mk$set <- factor(mk$set, levels = names(markers))
i <- match(mk$SYMBOL, rna$SYMBOL)
mk$rna_lfc <- rna$rna_lfc[i]; mk$rna_padj <- rna$rna_padj[i]; mk$rna_base <- rna$rna_base[i]
mk$rna_expressed <- rna$expressed[i]
mk$rna_sig <- !is.na(mk$rna_padj) & mk$rna_padj < rna_padj_cut
j <- match(mk$SYMBOL, atac$SYMBOL)
mk$atac_lfc <- atac$log2FC_shrunk[j]; mk$atac_padj <- atac$padj[j]; mk$atac_sig <- atac$atac_sig[j]
write.csv(mk, file.path(outdir, "marker_genes_table.csv"), row.names = FALSE)
nf <- sum(is.na(mk$rna_lfc))
if (nf) message(nf, " marker genes not in the RNA table: ", paste(mk$SYMBOL[is.na(mk$rna_lfc)], collapse = ", "))

# ---------------- per-set statistics ----------------
# Compare each set with all OTHER expressed genes (RNA) / other promoter peaks (ATAC)
set_stat <- function(vals, bg) {
  vals <- vals[!is.na(vals)]
  if (length(vals) < 3) return(c(n = length(vals), median = NA, mean = NA, p = NA))
  c(n = length(vals), median = median(vals), mean = mean(vals),
    p = suppressWarnings(wilcox.test(vals, bg)$p.value))
}
rows <- lapply(levels(mk$set), function(s) {
  g <- mk$SYMBOL[mk$set == s]
  r <- mk[mk$set == s & mk$rna_expressed %in% TRUE, ]
  rb <- rna$rna_lfc[rna$expressed & !(rna$SYMBOL %in% g)]
  ab <- atac$log2FC_shrunk[!(atac$SYMBOL %in% g)]
  rs <- set_stat(r$rna_lfc, rb); as_ <- set_stat(mk$atac_lfc[mk$set == s], ab)
  data.frame(set = s, n_listed = length(g),
             rna_n_expressed = rs["n"], rna_median_lfc = rs["median"], rna_mean_lfc = rs["mean"],
             rna_n_up_sig = sum(r$rna_sig & r$rna_lfc > 0), rna_n_down_sig = sum(r$rna_sig & r$rna_lfc < 0),
             rna_p = rs["p"],
             atac_n_with_peak = as_["n"], atac_median_lfc = as_["median"], atac_mean_lfc = as_["mean"],
             atac_p = as_["p"], row.names = NULL)
})
st <- do.call(rbind, rows)
st$rna_padj  <- p.adjust(st$rna_p, "BH")
st$atac_padj <- p.adjust(st$atac_p, "BH")
write.csv(st, file.path(outdir, "marker_sets_stats.csv"), row.names = FALSE)
print(format(st[, c("set", "rna_n_expressed", "rna_median_lfc", "rna_n_up_sig", "rna_n_down_sig", "rna_padj",
                    "atac_n_with_peak", "atac_median_lfc", "atac_padj")], digits = 2), row.names = FALSE)

# ---------------- figure 1: set-level summary (RNA and ATAC side by side) ----------------
mod_rna  <- "RNA-seq"
mod_atac <- paste0("ATAC promoter peak (<= ", atac_max_dist / 1000, " kb)")
long <- rbind(
  data.frame(set = mk$set, SYMBOL = mk$SYMBOL, modality = mod_rna, lfc = mk$rna_lfc,
             sig = mk$rna_sig, shown = mk$rna_expressed %in% TRUE),
  data.frame(set = mk$set, SYMBOL = mk$SYMBOL, modality = mod_atac, lfc = mk$atac_lfc,
             sig = mk$atac_sig %in% TRUE, shown = TRUE))
long <- long[!is.na(long$lfc), ]
long$modality <- factor(long$modality, levels = c(mod_rna, mod_atac))
ann <- rbind(
  data.frame(set = st$set, modality = mod_rna, n = st$rna_n_expressed, padj = st$rna_padj),
  data.frame(set = st$set, modality = mod_atac, n = st$atac_n_with_peak, padj = st$atac_padj))
ann$set <- factor(ann$set, levels = names(markers)); ann$modality <- factor(ann$modality, levels = levels(long$modality))
ann$label <- sprintf("n=%d\nadj p %s", ann$n, fmt_p(ann$padj))

ylim <- range(c(long$lfc, 0), na.rm = TRUE)
p1 <- ggplot(long[long$shown, ], aes(set, lfc)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_boxplot(outlier.shape = NA, fill = NA, width = 0.55, color = "grey30") +
  geom_jitter(aes(color = sig), width = 0.15, height = 0, size = 1.6, alpha = 0.85) +
  geom_text(data = ann, aes(label = label), y = ylim[2], vjust = 1.1, size = 2.7, lineheight = 0.9) +
  scale_color_manual(values = c("TRUE" = "firebrick", "FALSE" = "grey60"), labels = c("TRUE" = "padj < 0.05", "FALSE" = "n.s."), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  facet_wrap(~ modality, ncol = 1) +
  labs(x = NULL, y = "log2FC (KO/WT)", title = "Cell-type marker sets in KO vs WT",
       subtitle = "Each point is one marker gene; adj p = Wilcoxon vs all other genes (BH across sets)") +
  theme_bw() + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(outdir, "marker_sets_summary.pdf"), p1, width = 7.5, height = 7.5)

# ---------------- figure 2: every marker gene, RNA ----------------
g <- mk[!is.na(mk$rna_lfc), ]
g <- g[order(g$set, g$rna_lfc), ]
g$SYMBOL <- factor(g$SYMBOL, levels = g$SYMBOL)
g$status <- ifelse(!g$rna_expressed, "low expression", ifelse(g$rna_sig, "padj < 0.05", "n.s."))
p2 <- ggplot(g, aes(rna_lfc, SYMBOL, fill = status)) +
  geom_col(width = 0.7) +
  geom_vline(xintercept = 0, color = "grey30") +
  scale_fill_manual(values = c("padj < 0.05" = "firebrick", "n.s." = "grey60", "low expression" = "grey88"), name = NULL) +
  facet_grid(set ~ ., scales = "free_y", space = "free_y") +
  labs(x = "RNA log2FC (KO/WT)", y = NULL, title = "Marker genes, RNA-seq",
       subtitle = paste0("Light bars: baseMean < ", min_basemean, " (barely expressed; excluded from the statistics)")) +
  theme_bw() + theme(strip.text.y = element_text(angle = 0), axis.text.y = element_text(size = 7))
ggsave(file.path(outdir, "marker_genes_rna.pdf"), p2, width = 7, height = max(5, 1.5 + 0.17 * nrow(g)))

message("Done. Outputs in: ", outdir)