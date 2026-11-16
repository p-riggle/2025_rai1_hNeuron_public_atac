# Build BED files of the ATAC peaks at CONCORDANT genes, for HOMER motif analysis.
# Runs locally in RStudio after atac_rna_correlation.R. Then upload the three BEDs to
# Great Lakes (Data/Results/Motifs_concordant/) and run atac_homer_concordant.sbatch.
#
# Concordant up   = peak more open in KO AND gene up in KO
# Concordant down = peak less open in KO AND gene down in KO
# One peak per gene: the one the correlation script used (strongest peak near the gene).

# ======================= settings =======================
merged_file <- "results/ATACseq/RNA_correlation/ATAC_RNA_merged_all_nearest_gene.csv"
peaks_file  <- "results/ATACseq/annotation/KO_vs_WT_all_peaks_annotated.csv"   # peak, chr, start, end
outdir      <- "results/ATACseq/RNA_correlation/motifs"
# ========================================================

merged_file <- path.expand(merged_file); peaks_file <- path.expand(peaks_file)
outdir <- path.expand(outdir)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

m  <- read.csv(merged_file, stringsAsFactors = FALSE)
pk <- read.csv(peaks_file,  stringsAsFactors = FALSE)
for (nm in list(list(m, c("peak", "cat", "atac_lfc"), "merged table"),
                list(pk, c("peak", "chr", "start", "end"), "peaks table"))) {
  miss <- setdiff(nm[[2]], names(nm[[1]]))
  if (length(miss)) stop(nm[[3]], " is missing: ", paste(miss, collapse = ", "), call. = FALSE)
}

conc <- m[m$cat == "Concordant", ]
up_ids   <- unique(conc$peak[conc$atac_lfc > 0])   # a peak shared by two genes is kept once
down_ids <- unique(conc$peak[conc$atac_lfc < 0])

to_bed <- function(ids) {
  d <- pk[match(ids, pk$peak), c("chr", "start", "end", "peak")]
  if (anyNA(d$chr)) stop("Some peak IDs were not found in the peaks table.", call. = FALSE)
  d$score <- 0; d$strand <- "+"
  d
}
write_bed <- function(d, f) write.table(d, file.path(outdir, f), sep = "\t", quote = FALSE,
                                        row.names = FALSE, col.names = FALSE)

up_bed   <- to_bed(up_ids)
down_bed <- to_bed(down_ids)
# Background: every other tested peak (all consensus peaks that went into DESeq2),
# minus the two foreground sets so targets never sit in their own background.
bg_bed <- to_bed(setdiff(pk$peak, c(up_ids, down_ids)))

write_bed(up_bed,   "concordant_up_peaks.bed")
write_bed(down_bed, "concordant_down_peaks.bed")
write_bed(bg_bed,   "background_peaks.bed")

message("Concordant genes: ", nrow(conc), " -> peaks: up = ", nrow(up_bed), ", down = ", nrow(down_bed),
        "; background = ", nrow(bg_bed))
message("Peak widths (bp): ", paste(range(up_bed$end - up_bed$start), collapse = "-"),
        " (HOMER will use the central 200 bp)")
if (min(nrow(up_bed), nrow(down_bed)) < 200)
  message("NOTE: fewer than ~200 peaks in a set; HOMER has little power with sets this small.")
message("Done. Upload the three BED files in ", outdir, " to Data/Results/Motifs_concordant/ on Great Lakes.")
