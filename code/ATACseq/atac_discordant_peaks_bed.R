# Build BED files of the ATAC peaks at DISCORDANT genes, for HOMER motif analysis.
# Runs after atac_rna_correlation.R. Then copy the three BEDs to Great Lakes
# (Data/Results/Motifs_discordant/) and run atac_homer_discordant.sbatch.
#
# Discordant open   = peak more open in KO AND gene down in KO
# Discordant closed = peak less open in KO AND gene up in KO
# One peak per gene: the one the correlation script used (strongest peak near the gene).

# ======================= settings =======================
merged_file <- "results/ATACseq/RNA_correlation/ATAC_RNA_merged_all_nearest_gene.csv"
peaks_file  <- "results/ATACseq/annotation/KO_vs_WT_all_peaks_annotated.csv"   # peak, chr, start, end
outdir      <- "results/ATACseq/RNA_correlation/motifs_discordant"
# ========================================================

merged_file <- path.expand(merged_file); peaks_file <- path.expand(peaks_file)
outdir <- path.expand(outdir)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

m  <- read.csv(merged_file, stringsAsFactors = FALSE)
pk <- read.csv(peaks_file,  stringsAsFactors = FALSE)
for (nm in list(list(m, c("peak", "cat", "atac_lfc", "rna_lfc"), "merged table"),
                list(pk, c("peak", "chr", "start", "end"), "peaks table"))) {
  miss <- setdiff(nm[[2]], names(nm[[1]]))
  if (length(miss)) stop(nm[[3]], " is missing: ", paste(miss, collapse = ", "), call. = FALSE)
}

disc <- m[m$cat == "Discordant", ]
open_ids   <- unique(disc$peak[disc$atac_lfc > 0 & disc$rna_lfc < 0])   # a peak shared by two genes is kept once
closed_ids <- unique(disc$peak[disc$atac_lfc < 0 & disc$rna_lfc > 0])
shared <- intersect(open_ids, closed_ids)
if (length(shared)) stop(length(shared), " peaks fall in both sets.", call. = FALSE)

to_bed <- function(ids) {
  d <- pk[match(ids, pk$peak), c("chr", "start", "end", "peak")]
  if (anyNA(d$chr)) stop("Some peak IDs were not found in the peaks table.", call. = FALSE)
  d$score <- 0; d$strand <- "+"
  d
}
write_bed <- function(d, f) write.table(d, file.path(outdir, f), sep = "\t", quote = FALSE,
                                        row.names = FALSE, col.names = FALSE)

open_bed   <- to_bed(open_ids)
closed_bed <- to_bed(closed_ids)
# Background: every other tested peak (all consensus peaks that went into DESeq2),
# minus the two foreground sets so targets never sit in their own background.
bg_bed <- to_bed(setdiff(pk$peak, c(open_ids, closed_ids)))

write_bed(open_bed,   "discordant_open_peaks.bed")
write_bed(closed_bed, "discordant_closed_peaks.bed")
write_bed(bg_bed,     "background_peaks.bed")

message("Discordant genes: ", nrow(disc), " -> peaks: open (gene down) = ", nrow(open_bed),
        ", closed (gene up) = ", nrow(closed_bed), "; background = ", nrow(bg_bed))
message("Peak widths (bp): ", paste(range(open_bed$end - open_bed$start), collapse = "-"),
        " (HOMER will use the central 200 bp)")
if (min(nrow(open_bed), nrow(closed_bed)) < 200)
  message("NOTE: fewer than ~200 peaks in a set; HOMER has little power with sets this small.")
message("Done. Copy the three BED files in ", outdir, " to Data/Results/Motifs_discordant/ on Great Lakes.")
