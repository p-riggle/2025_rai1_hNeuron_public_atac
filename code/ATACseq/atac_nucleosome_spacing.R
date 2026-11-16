# Nucleosome spacing from ATAC-seq fragment sizes: RAI1 KO vs WT
# Runs locally in RStudio on the small tables made by atac_fragment_sizes.sbatch.
#
# Download from  Data/Results/Nucleosome/  into one folder (set `indir`):
#   WT1.fragment_lengths.tsv  WT2.fragment_lengths.tsv  KO1.fragment_lengths.tsv  KO2.fragment_lengths.tsv
#
# Outputs (in `outdir`):
#   fragment_size_distribution_log10.pdf / _linear.pdf   all four samples, nucleosome classes shaded
#   fragment_size_KO_vs_WT_log2ratio.pdf                 where the size distribution shifts in KO
#   fragment_size_classes.pdf / .csv                     fraction of fragments in NFR / mono / di / tri
#   nucleosome_peak_positions.pdf / fragment_size_peaks.csv   mono / di / tri peak positions and spacing
#
# Only needs ggplot2:  install.packages("ggplot2")

suppressPackageStartupMessages(library(ggplot2))

# ======================= settings =======================
indir    <- "data/ATACseq/nucleosome"
outdir   <- file.path("results/ATACseq/nucleosome")
samples  <- c("WT1", "WT2", "KO1", "KO2")
max_len  <- 800        # x-axis limit for the plots (counts go to 1000 bp)
smooth_w <- 21         # running-mean width (bp) used for peak finding and the ratio plot;
#   21 averages over two turns of the ~10 bp DNA-helix ripple

# Fragment classes (ATACseqQC / Ou et al. 2018 boundaries)
classes <- data.frame(name = c("NFR", "Mono", "Di", "Tri"),
                      xmin = c(1, 180, 315, 558), xmax = c(99, 247, 473, 615))
# Windows for finding the nucleosome peaks (bp)
peak_windows <- list(Mono = c(150, 300), Di = c(320, 520), Tri = c(540, 720))

cols <- c(WT1 = "#2b6cb0", WT2 = "#63b3ed", KO1 = "#c53030", KO2 = "#fc8181")
gcol <- c(WT = "#2b6cb0", KO = "#c53030")
# ========================================================

dir.create(path.expand(outdir), recursive = TRUE, showWarnings = FALSE)
indir <- path.expand(indir); outdir <- path.expand(outdir)

# ---------------- load ----------------
frag <- do.call(rbind, lapply(samples, function(s) {
  f <- file.path(indir, paste0(s, ".fragment_lengths.tsv"))
  if (!file.exists(f)) stop("Missing file: ", f, call. = FALSE)
  d <- read.delim(f)
  stopifnot(all(c("length", "count") %in% names(d)))
  d <- d[order(d$length), ]
  d$sample   <- s
  d$genotype <- sub("[0-9]+$", "", s)
  d$frac     <- d$count / sum(d$count)
  d$smooth   <- as.numeric(stats::filter(d$frac, rep(1 / smooth_w, smooth_w), sides = 2))
  d
}))
frag$sample   <- factor(frag$sample, levels = samples)
frag$genotype <- factor(frag$genotype, levels = c("WT", "KO"))
totals <- tapply(frag$count, frag$sample, sum)
message("Fragments counted (<= 1000 bp): ", paste(names(totals), format(totals, big.mark = ","), collapse = " | "))

# Shaded class boxes. Finite y limits are needed (-Inf/Inf are dropped on a log axis);
# pass the range of the data being plotted so the boxes never change the axis range.
shade <- function(ymin, ymax) list(
  annotate("rect", xmin = classes$xmin, xmax = classes$xmax, ymin = ymin, ymax = ymax,
           fill = "grey85", alpha = 0.5),
  annotate("text", x = (classes$xmin + classes$xmax) / 2, y = ymax, label = classes$name,
           vjust = 1.4, size = 3, color = "grey30"))

# ---------------- 1. fragment size distribution ----------------
for (ys in c("log10", "linear")) {
  d <- frag[frag$length <= max_len & frag$frac > 0, ]
  p <- ggplot(d, aes(length, frac, color = sample)) +
    shade(if (ys == "log10") min(d$frac) else 0, max(d$frac)) +
    geom_line(linewidth = 0.5) +
    scale_color_manual(values = cols, name = NULL) +
    labs(x = "Fragment length (bp)", y = "Fraction of fragments",
         title = "ATAC-seq fragment size distribution", subtitle = "RAI1 KO vs WT") +
    theme_bw()
  if (ys == "log10") p <- p + scale_y_log10()
  ggsave(file.path(outdir, paste0("fragment_size_distribution_", ys, ".pdf")), p, width = 7, height = 4.5)
}

# ---------------- 2. KO / WT log2 ratio of the size distribution ----------------
wide <- reshape(frag[, c("sample", "length", "smooth")], idvar = "length", timevar = "sample", direction = "wide")
wt <- rowMeans(wide[, paste0("smooth.", grep("^WT", samples, value = TRUE))])
ko <- rowMeans(wide[, paste0("smooth.", grep("^KO", samples, value = TRUE))])
ratio <- data.frame(length = wide$length, log2ratio = log2(ko / wt))
ratio <- ratio[is.finite(ratio$log2ratio) & ratio$length <= max_len, ]
pr <- ggplot(ratio, aes(length, log2ratio)) +
  shade(min(ratio$log2ratio), max(ratio$log2ratio)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_line(linewidth = 0.6, color = "black") +
  labs(x = "Fragment length (bp)", y = "log2(KO / WT) fraction of fragments",
       title = "Shift in fragment size distribution",
       subtitle = paste0("Replicate means, ", smooth_w, " bp running mean; above 0 = KO enriched")) +
  theme_bw()
ggsave(file.path(outdir, "fragment_size_KO_vs_WT_log2ratio.pdf"), pr, width = 7, height = 4)

# ---------------- 3. fraction of fragments per class ----------------
cls <- do.call(rbind, lapply(samples, function(s) {
  d <- frag[frag$sample == s, ]
  data.frame(sample = s, genotype = d$genotype[1],
             class = classes$name,
             fraction = sapply(seq_len(nrow(classes)), function(i)
               sum(d$frac[d$length >= classes$xmin[i] & d$length <= classes$xmax[i]])))
}))
cls$class <- factor(cls$class, levels = classes$name)
cls$genotype <- factor(cls$genotype, levels = c("WT", "KO"))
write.csv(cls, file.path(outdir, "fragment_size_classes.csv"), row.names = FALSE)

pc <- ggplot(cls, aes(genotype, fraction, fill = genotype)) +
  stat_summary(fun = mean, geom = "col", width = 0.6, alpha = 0.6) +
  geom_point(aes(shape = sample), size = 2.5, position = position_jitter(width = 0.08, height = 0, seed = 1)) +
  scale_fill_manual(values = gcol, guide = "none") +
  facet_wrap(~ class, scales = "free_y", nrow = 1) +
  labs(x = NULL, y = "Fraction of fragments", shape = NULL,
       title = "Fragments per size class",
       subtitle = "NFR < 100 bp, mono 180-247, di 315-473, tri 558-615 bp (bars = mean, points = replicates)") +
  theme_bw()
ggsave(file.path(outdir, "fragment_size_classes.pdf"), pc, width = 8, height = 3.8)

# ---------------- 4. nucleosome peak positions and spacing ----------------
find_peak <- function(d, win) {
  w <- d[d$length >= win[1] & d$length <= win[2] & !is.na(d$smooth), ]
  i <- which.max(w$smooth)
  c(pos = w$length[i], edge = w$length[i] %in% range(w$length))
}
pk <- do.call(rbind, lapply(samples, function(s) {
  d <- frag[frag$sample == s, ]
  r <- lapply(peak_windows, function(w) find_peak(d, w))
  data.frame(sample = s, genotype = as.character(d$genotype[1]),
             mono_bp = r$Mono["pos"], di_bp = r$Di["pos"], tri_bp = r$Tri["pos"],
             spacing_di_minus_mono = r$Di["pos"] - r$Mono["pos"],
             spacing_tri_minus_di  = r$Tri["pos"] - r$Di["pos"],
             mono_at_window_edge = as.logical(r$Mono["edge"]),
             di_at_window_edge   = as.logical(r$Di["edge"]),
             tri_at_window_edge  = as.logical(r$Tri["edge"]),
             row.names = NULL)
}))
write.csv(pk, file.path(outdir, "fragment_size_peaks.csv"), row.names = FALSE)
print(pk, row.names = FALSE)
if (any(pk$di_at_window_edge | pk$tri_at_window_edge))
  message("NOTE: a di- or tri-nucleosome peak sits at the edge of its search window, i.e. there is no ",
          "clear peak there (weak, or fragments > ~500 bp missing from the alignment). Do not trust that spacing.")

long <- reshape(pk[, c("sample", "genotype", "mono_bp", "di_bp", "spacing_di_minus_mono")],
                direction = "long", idvar = c("sample", "genotype"),
                varying = c("mono_bp", "di_bp", "spacing_di_minus_mono"), v.names = "bp",
                timevar = "metric", times = c("Mono peak", "Di peak", "Di - mono spacing"))
long$metric   <- factor(long$metric, levels = c("Mono peak", "Di peak", "Di - mono spacing"))
long$genotype <- factor(long$genotype, levels = c("WT", "KO"))
pp <- ggplot(long, aes(genotype, bp, color = genotype)) +
  stat_summary(fun = mean, geom = "crossbar", width = 0.5, linewidth = 0.3) +
  geom_point(size = 2.8, position = position_jitter(width = 0.06, height = 0, seed = 1)) +
  scale_color_manual(values = gcol, guide = "none") +
  facet_wrap(~ metric, scales = "free_y", nrow = 1) +
  labs(x = NULL, y = "bp", title = "Nucleosome peak positions and spacing",
       subtitle = "Peak = maximum of the smoothed fragment-size distribution (points = replicates, bar = mean)") +
  theme_bw()
ggsave(file.path(outdir, "nucleosome_peak_positions.pdf"), pp, width = 7.5, height = 3.8)

message("Done. Outputs in: ", outdir)
