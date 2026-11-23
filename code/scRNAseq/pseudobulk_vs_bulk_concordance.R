# ============================================================================
# pseudobulk_vs_bulk_concordance.R
#
# Purpose
#   Build an all-cell ("in-silico bulk") pseudobulk for every scRNA-seq sample
#   (10x Flex; day 20 and day 60; WT / HET / KO; 2 clones each) and compare it
#   with the real bulk RNA-seq of the same cortical-neuron differentiation
#   (code/cn_rnaseq; d20 and d60; 4 replicates each) to assess
#     (1) gene-level agreement of expression levels (all genes, collagen
#         genes highlighted),
#     (2) gene-level agreement of genotype effects (log2 fold changes),
#     (3) batch effects (platform, clone / replicate) and replicate variability.
#
# Inputs (read-only; nothing is written under symlink_rai1_hesc/)
#   * ../symlink_rai1_hesc/data/scRNAseq/iwaselab/cts
#       dgCMatrix of raw RNA-assay counts, genes x (sample _ cell type),
#       produced by AggregateExpression(geo_so, assays = "RNA", slot = "counts",
#       group.by = c("orig.ident", "cell_type")) in code/scRNAseq/6_pseudobulk.R.
#       Summing its columns over cell types gives the all-cell pseudobulk.
#   * ../symlink_rai1_hesc/data/cn_rnaseq/blueberry_output/combine_counts/
#       gene_expected_count.txt and gene_TPM.txt (RSEM, gene symbols).
#   * results/cn_rnaseq/deseq2/rstbl_d{20,60}_wt.vs.{het,ko}.txt
#       Published bulk DESeq2 results (code/cn_rnaseq/deseq2/s1_deseq2.r).
#       Despite the file names, log2FoldChange = log2(mutant / WT).
#   * ../symlink_rai1_hesc/data/scRNAseq/10x_analysis_11534-BZ/.../
#       sample_filtered_feature_bc_matrix.h5 (Cell Ranger, one per sample; the
#       same matrices read by code/scRNAseq/1_first_steps.r). Used for per-cell
#       QC, to check that `cts` is reproduced exactly, and to build the
#       QC-filtered pseudobulk (PB_CELL_FILTER=qc).
#
# Usage (from the repository root 2025_rai1_hNeuron_public_atac/):
#   sbatch code/scRNAseq/pseudobulk_vs_bulk_concordance.sbatch
#   or, inside the single_cell_env_0.12.0.sif container:
#   Rscript code/scRNAseq/pseudobulk_vs_bulk_concordance.R
#   PB_CELL_FILTER=qc Rscript code/scRNAseq/pseudobulk_vs_bulk_concordance.R
#
# ----------------------------------------------------------------------------
# scRNA-seq CELL FILTER (PB_CELL_FILTER; default "none")
# ----------------------------------------------------------------------------
#   * "none" -> results/scRNAseq/pseudobulk_vs_bulk/
#       All cells in the published object, i.e. Cell Ranger-called cells with
#       >= 50 detected genes (CreateSeuratObject(min.features = 50)). This is
#       what the published pipeline used: 1_first_steps.r previews
#       subset(nFeature_RNA > 500 & percent.mt < 15) but never assigns it, and
#       results/scRNAseq/1_qc/1_cell_filtering_counts.csv has prefilter =
#       postfilter for every sample. Counts come from `cts`.
#   * "qc"   -> results/scRNAseq/pseudobulk_vs_bulk_cellqc/
#       The threshold 1_first_steps.r intended, applied per cell to the
#       Cell Ranger matrices: nFeature_RNA > 500 AND percent.mt < 15
#       (strict inequalities, as in the Seurat subset() call; percent.mt =
#       100 x counts on ^MT- genes / total counts, as PercentageFeatureSet()).
#       Rationale: < 500 detected genes flags empty/ambient droplets or
#       damaged cells; high mitochondrial fraction flags cells with lost
#       cytoplasmic RNA. Filtering acts before cell typing, so the all-cell
#       pseudobulk needs no cell-type labels.
#   * In both modes the all-cell reconstruction from the .h5 files is checked
#     against `cts` (must agree exactly); the QC mode stops if it does not.
#
# ----------------------------------------------------------------------------
# STATISTICAL PARAMETERS AND RATIONALE (values are set in CONFIG below)
# ----------------------------------------------------------------------------
#  Unit of replication
#   * The SAMPLE is the unit of replication (scRNA-seq: one hESC clone at one
#     timepoint; bulk: one replicate culture). Cells are never treated as
#     replicates; this avoids pseudoreplication.
#   * Bulk replicates (rep1-4) and scRNA-seq clones (rep1-2) are NOT assumed to
#     be paired: platforms are compared at the level of condition means.
#
#  Normalisation for expression-level comparison
#   * 10x Flex is probe-based, so its counts do not scale with transcript
#     length; full-length bulk RNA-seq counts do. Bulk is therefore expressed as
#     TPM (length-normalised) and pseudobulk as TMM-normalised CPM (edgeR).
#   * Both are re-scaled to 1e6 over the SHARED gene universe (genes in both the
#     Flex probe set and the bulk annotation), so both are "per million over the
#     same genes".
#   * log2(x + 1) is used for plotting and correlation (pseudocount = 1 per
#     million: stabilises variance for low values and keeps zeros finite).
#
#  Correlation statistics (expression and log2FC)
#   * Pearson r (linear association, as requested) with a 95% CI from
#     cor.test() (Fisher z). With ~15,000 genes every p-value is effectively 0
#     and genes are not independent, so effect sizes and CIs are reported,
#     not p-values, in figures. p-values are still written to the tables.
#   * Spearman rho (rank-based, robust to outliers / non-linearity) with a 95%
#     CI from Fisher z using SE = sqrt(1.06 / (n - 3)) (Fieller et al. 1957).
#   * OLS fit (y = pseudobulk, x = bulk): slope, intercept and R^2 describe the
#     linear relationship; the fitted line is drawn on each plot.
#   * Lin's concordance correlation coefficient (CCC) is reported for log2FC
#     only, where both axes are on the same scale and agreement with the
#     identity line is meaningful. It is not meaningful across platforms for
#     absolute expression (different units / capture), so it is omitted there.
#   * Sign concordance of log2FC: among genes significant in the bulk data
#     (padj < 0.05), proportion with the same direction in the pseudobulk;
#     exact binomial test against 0.5 (chance).
#
#  Gene filters
#   * Expression scatter: ALL shared genes are plotted. Statistics are reported
#     for all shared genes and for "expressed" genes (mean >= 1 per million in
#     EITHER platform). The "either" rule avoids selecting genes on agreement
#     between platforms, which would inflate correlation.
#   * log2FC scatter: genes with >= 10 total counts across the samples of the
#     contrast in the pseudobulk (DESeq2 pre-filter convention) and a non-NA
#     bulk log2FC. Collagen genes failing this filter are listed in the tables.
#
#  Differential expression (pseudobulk)
#   * DESeq2, design ~ conditions (day x genotype groups, both timepoints in one
#     model so dispersion is estimated from all 12 samples), betaPrior = TRUE
#     (normal prior on log2FC). This mirrors the published bulk model
#     (s1_deseq2.r: ~ conditions, betaPrior = TRUE), so the two sets of log2FC
#     are shrunk with the same estimator and are directly comparable.
#   * n = 2 clones per genotype: very low power. log2FC agreement, not the
#     number of significant genes, is the comparison of interest.
#
#  Collagen genes
#   * Collagen alpha-chain genes = HGNC symbols matching ^COL[0-9]+A[0-9]+$
#     (e.g. COL1A1, COL4A2). Excludes non-collagen genes that start with "COL"
#     (COLEC*, COLGALT*, COLQ).
#   * Collagen-specific deviation: an OLS line is fitted to NON-collagen genes;
#     the mean residual of collagen genes measures whether collagen genes sit
#     systematically above/below the genome-wide relationship (e.g. collagen
#     under-represented in pseudobulk, as expected if collagen-rich mesenchymal
#     cells are lost during dissociation). Significance: expression-matched
#     permutation test, 10,000 random gene sets of the same size, drawn from the
#     same decile bins of bulk expression as the collagen genes (expression
#     level strongly affects noise and probe/read capture). Two-sided empirical
#     p = (1 + #|null| >= |observed|) / (1 + n_perm). Benjamini-Hochberg FDR
#     within each family (expression conditions; log2FC contrasts).
#
#  Batch effects and variability
#   * Combined PCA (bulk d20/d60 + pseudobulk) on log2 per-million values of
#     expressed shared genes, top 2,000 most variable genes, centred, not scaled.
#     Shown before and after limma::removeBatchEffect(batch = platform,
#     design = ~ day * genotype). The correction is for visualisation only: it
#     removes the platform main effect while preserving biology, and asks
#     whether samples then group by day / genotype rather than platform.
#   * Cross-platform correlation: genes centred within platform (removes the
#     platform-specific gene-level offset, e.g. probe vs read capture), then
#     Pearson r between every pseudobulk and every bulk sample on the top 2,000
#     variable genes. A pseudobulk sample whose best match is a bulk sample of
#     the same day/genotype indicates concordant biology.
#   * Within-platform variance decomposition: DESeq2 vst (blind = TRUE) on each
#     platform's d20+d60 samples, PCA on the top 2,000 variable genes, and for
#     PCs 1-5 the R^2 of lm(PC ~ factor) for day, genotype and clone/replicate
#     (each factor separately). A replicate/clone factor explaining a large
#     share of a top PC indicates a replicate-level batch effect. The scRNA-seq
#     samples were multiplexed in one Flex pool (BC001-BC012, run 11534-BZ),
#     so there is no library-prep batch within scRNA-seq.
#   * Replicate variability: (a) Pearson r between replicates of the same
#     day x genotype group (expressed genes, log2 per million); (b) edgeR
#     biological coefficient of variation (BCV = sqrt(common dispersion)),
#     per platform x day, design ~ genotype, filterByExpr + TMM.
#
#  Bulk samples excluded
#   * Rai1-het_d34_rep1 and Rai1-het_d60_rep4, excluded as outliers in the
#     published bulk analysis (s1_deseq2.r). Kept identical for consistency.
#
#  Randomness: set.seed(42) before permutations.
# ============================================================================

# ============================================================================
# CONFIGURATION
# ============================================================================
pb_cell_filter <- match.arg(Sys.getenv("PB_CELL_FILTER", "none"), c("none", "qc"))

CONFIG <- list(
  # ---- Inputs (relative to the repository root) ----
  sc_cts_path      = "../symlink_rai1_hesc/data/scRNAseq/iwaselab/cts",
  bulk_counts_path = "../symlink_rai1_hesc/data/cn_rnaseq/blueberry_output/combine_counts/gene_expected_count.txt",
  bulk_tpm_path    = "../symlink_rai1_hesc/data/cn_rnaseq/blueberry_output/combine_counts/gene_TPM.txt",
  bulk_deseq_paths = c(
    d20_het = "results/cn_rnaseq/deseq2/rstbl_d20_wt.vs.het.txt",
    d20_ko  = "results/cn_rnaseq/deseq2/rstbl_d20_wt.vs.ko.txt",
    d60_het = "results/cn_rnaseq/deseq2/rstbl_d60_wt.vs.het.txt",
    d60_ko  = "results/cn_rnaseq/deseq2/rstbl_d60_wt.vs.ko.txt"
  ),
  # ---- scRNA-seq cell-level QC (see header "scRNA-seq CELL FILTER") ----
  sc_cell_filter = pb_cell_filter,
  # Per-cell QC from the Cell Ranger .h5 files (streamed; < 2 GB RAM, a few minutes).
  # Always run when sc_cell_filter = "qc". Skip with PB_RUN_CELL_QC=FALSE.
  run_cell_qc    = as.logical(Sys.getenv("PB_RUN_CELL_QC", "TRUE")),
  sc_h5_template = "../symlink_rai1_hesc/data/scRNAseq/10x_analysis_11534-BZ/Sample_11534-BZ-Pool01/per_sample_outs/%s/count/sample_filtered_feature_bc_matrix.h5",
  # orig.ident -> Cell Ranger sample, exactly as in 1_first_steps.r (sample_dirs)
  sc_h5_samples = c(day20_wt_rep1  = "11534-BZ-2",  day20_wt_rep2  = "11534-BZ-4",
                    day20_ko_rep1  = "11534-BZ-6",  day20_ko_rep2  = "11534-BZ-7",
                    day20_het_rep1 = "11534-BZ-10", day20_het_rep2 = "11534-BZ-12",
                    day60_wt_rep1  = "11534-BZ-14", day60_wt_rep2  = "11534-BZ-16",
                    day60_ko_rep1  = "11534-BZ-18", day60_ko_rep2  = "11534-BZ-20",
                    day60_het_rep1 = "11534-BZ-22", day60_het_rep2 = "11534-BZ-24"),
  seurat_min_features = 50,    # CreateSeuratObject(min.features = 50) in 1_first_steps.r
  qc_min_features     = 500,   # keep nFeature_RNA >  500 (1_first_steps.r)
  qc_max_percent_mt   = 15,    # keep percent.mt   <  15  (1_first_steps.r)
  mt_regex            = "^MT-",
  h5_chunk_cells      = 2000,  # cells read per block (memory control only)

  # ---- Outputs ----
  out_dir = if (pb_cell_filter == "qc") "results/scRNAseq/pseudobulk_vs_bulk_cellqc" else
    "results/scRNAseq/pseudobulk_vs_bulk",

  # ---- Design ----
  days            = c("d20", "d60"),                 # timepoints present on both platforms
  genotypes       = c("wt", "het", "ko"),
  contrasts       = list(het = c("het", "wt"), ko = c("ko", "wt")),  # c(test, reference)
  bulk_exclude    = c("Rai1-het_d34_rep1", "Rai1-het_d60_rep4"),

  # ---- Statistical parameters (see header for rationale) ----
  log_pseudocount        = 1,      # log2(per-million + 1)
  expressed_min_pm       = 1,      # "expressed": mean >= 1 per million in either platform
  min_total_counts_deseq = 10,     # pseudobulk pre-filter for log2FC comparison
  deseq_beta_prior       = TRUE,   # same estimator as the published bulk analysis
  padj_cutoff            = 0.05,
  n_hvg                  = 2000,   # variable genes for PCA / correlation heatmap
  n_pcs_decomp           = 5,
  n_perm                 = 10000,
  n_expr_bins            = 10,     # expression-matching bins for the permutation null
  conf_level             = 0.95,
  collagen_regex         = "^COL[0-9]+A[0-9]+$",
  random_seed            = 42,

  # ---- Plotting ----
  base_size        = 8,            # pt; journal figures are typically 6-8 pt
  col_collagen     = "#C0392B",
  col_other        = "grey65",
  genotype_colors  = c(WT = "#7b7b7b", HET = "#E3B505", KO = "#bd6565"),
  platform_shapes  = c(Bulk = 21, Pseudobulk = 24),
  day_colors       = c(d20 = "#4C9BE8", d60 = "#1B3A6B"),
  n_label_collagen = 8,
  raster_dpi       = 600,
  png_dpi          = 300
)

# ============================================================================
# PACKAGES
# ============================================================================
pkgs <- c("Matrix", "DESeq2", "edgeR", "limma", "dplyr", "tidyr", "tibble", "readr",
          "ggplot2", "ggrepel", "patchwork", "pheatmap", "RColorBrewer", "scales")
missing_pkgs <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) stop("Missing packages: ", paste(missing_pkgs, collapse = ", "))

suppressPackageStartupMessages({
  library(Matrix)
  library(DESeq2)
  library(edgeR)
  library(limma)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})
use_raster <- requireNamespace("ggrastr", quietly = TRUE)

# ============================================================================
# HELPERS
# ============================================================================
log_msg <- function(...) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), sprintf(...)))

fig_dir <- file.path(CONFIG$out_dir, "figures")
tab_dir <- file.path(CONFIG$out_dir, "tables")
for (d in c(fig_dir, tab_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

save_plot <- function(p, name, width, height) {
  ggsave(file.path(fig_dir, paste0(name, ".pdf")), p, width = width, height = height,
         units = "in", device = grDevices::pdf, useDingbats = FALSE)
  ggsave(file.path(fig_dir, paste0(name, ".png")), p, width = width, height = height,
         units = "in", dpi = CONFIG$png_dpi, bg = "white")
  log_msg("Saved figure %s", name)
}

save_pheatmap <- function(ph, name, width, height) {
  grDevices::pdf(file.path(fig_dir, paste0(name, ".pdf")), width = width, height = height,
                 useDingbats = FALSE)
  grid::grid.newpage(); grid::grid.draw(ph$gtable); grDevices::dev.off()
  grDevices::png(file.path(fig_dir, paste0(name, ".png")), width = width, height = height,
                 units = "in", res = CONFIG$png_dpi)
  grid::grid.newpage(); grid::grid.draw(ph$gtable); grDevices::dev.off()
  log_msg("Saved heatmap %s", name)
}

save_tbl <- function(df, name) readr::write_csv(df, file.path(tab_dir, paste0(name, ".csv")))

# Dense point layers are rasterised so PDFs stay small while text/lines stay vector.
rast <- function(layer) if (use_raster) ggrastr::rasterise(layer, dpi = CONFIG$raster_dpi) else layer

theme_pub <- function() {
  theme_classic(base_size = CONFIG$base_size, base_family = "sans") +
    theme(axis.text        = element_text(color = "black"),
          axis.line        = element_line(linewidth = 0.3),
          axis.ticks       = element_line(linewidth = 0.3),
          strip.background = element_blank(),
          strip.text       = element_text(face = "bold"),
          plot.title       = element_text(face = "bold", size = CONFIG$base_size + 1),
          legend.key.size  = unit(0.35, "cm"))
}

geno_label <- function(g) factor(toupper(g), levels = c("WT", "HET", "KO"))

lin_ccc <- function(x, y) {
  mx <- mean(x); my <- mean(y)
  2 * mean((x - mx) * (y - my)) / (mean((x - mx)^2) + mean((y - my)^2) + (mx - my)^2)
}

# Correlation / regression summary for y (pseudobulk) against x (bulk).
cor_stats <- function(x, y, with_ccc = FALSE) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; n <- length(x)
  if (n < 4) return(tibble(n = n))
  pt  <- cor.test(x, y, method = "pearson", conf.level = CONFIG$conf_level)
  rho <- suppressWarnings(cor(x, y, method = "spearman"))
  zq  <- qnorm(1 - (1 - CONFIG$conf_level) / 2)
  se  <- sqrt(1.06 / (n - 3))
  fit <- lm(y ~ x)
  out <- tibble(n = n,
                pearson_r = unname(pt$estimate),
                pearson_ci_low = pt$conf.int[1], pearson_ci_high = pt$conf.int[2],
                pearson_p = pt$p.value,
                spearman_rho = rho,
                spearman_ci_low = tanh(atanh(rho) - zq * se),
                spearman_ci_high = tanh(atanh(rho) + zq * se),
                ols_slope = unname(coef(fit)[2]), ols_intercept = unname(coef(fit)[1]),
                r_squared = summary(fit)$r.squared)
  if (with_ccc) out$lin_ccc <- lin_ccc(x, y)
  out
}

stat_label <- function(s, with_ccc = FALSE) {
  lab <- sprintf("r = %.2f [%.2f, %.2f]\nrho = %.2f\nR\u00b2 = %.2f, n = %s",
                 s$pearson_r, s$pearson_ci_low, s$pearson_ci_high, s$spearman_rho,
                 s$r_squared, format(s$n, big.mark = ","))
  if (with_ccc) lab <- paste0(lab, sprintf("\nCCC = %.2f", s$lin_ccc))
  lab
}

# Multi-line labels need hjust = 0 to stay left-aligned; leading spaces give the inset.
indent_lines <- function(x) gsub("(^|\n)", "\\1  ", x)

# Expression-matched permutation test for the mean residual of a gene set
# relative to an OLS line fitted on all other genes.
gene_set_residual_test <- function(x, y, in_set, strat, n_perm, n_bins) {
  ok <- is.finite(x) & is.finite(y) & is.finite(strat)
  x <- x[ok]; y <- y[ok]; in_set <- in_set[ok]; strat <- strat[ok]
  if (sum(in_set) < 3) return(tibble(n_set = sum(in_set)))
  fit   <- lm(y ~ x, subset = !in_set)
  resid <- y - (coef(fit)[1] + coef(fit)[2] * x)
  obs   <- mean(resid[in_set])
  brks  <- unique(quantile(strat, probs = seq(0, 1, length.out = n_bins + 1)))
  bins  <- cut(strat, brks, include.lowest = TRUE)
  need  <- table(bins[in_set])
  need  <- need[need > 0]
  pools <- lapply(names(need), function(b) resid[!in_set & bins == b])
  null  <- vapply(seq_len(n_perm), function(i) {
    mean(unlist(Map(function(pool, k) pool[sample.int(length(pool), k)], pools, as.integer(need))))
  }, numeric(1))
  tibble(n_set = sum(in_set), mean_residual = obs,
         null_mean = mean(null), null_sd = sd(null),
         z = (obs - mean(null)) / sd(null),
         p_perm = (1 + sum(abs(null - mean(null)) >= abs(obs - mean(null)))) / (1 + n_perm))
}

# Stream one Cell Ranger .h5 matrix (CSC, genes x cells) in blocks of cells.
# Returns per-cell QC metrics (computed as Seurat does) and per-gene count sums
# over (a) cells kept by CreateSeuratObject(min.features) and (b) cells that
# also pass the nFeature / percent.mt filter.
read_h5_pseudobulk <- function(path, sample) {
  f <- hdf5r::H5File$new(path, mode = "r")
  on.exit(f$close_all())
  genes  <- gsub("_", "-", make.unique(f[["matrix/features/name"]][]))   # as Read10X() + Seurat
  dims   <- as.numeric(f[["matrix/shape"]][])
  indptr <- as.numeric(f[["matrix/indptr"]][])
  n_cell <- dims[2]
  is_mt  <- grepl(CONFIG$mt_regex, genes)
  sum_kept <- sum_pass <- numeric(dims[1])
  cells <- vector("list", ceiling(n_cell / CONFIG$h5_chunk_cells))
  for (k in seq_along(cells)) {
    s <- (k - 1) * CONFIG$h5_chunk_cells + 1
    e <- min(s + CONFIG$h5_chunk_cells - 1, n_cell)
    lo <- indptr[s] + 1; hi <- indptr[e + 1]
    idx <- if (hi >= lo) f[["matrix/indices"]][lo:hi] else integer(0)
    val <- if (hi >= lo) f[["matrix/data"]][lo:hi]    else numeric(0)
    m <- Matrix::sparseMatrix(i = as.integer(idx) + 1L, p = as.integer(indptr[s:(e + 1)] - indptr[s]),
                              x = as.numeric(val), dims = c(dims[1], e - s + 1))
    n_feat <- diff(m@p)
    n_cnt  <- Matrix::colSums(m)
    pct_mt <- 100 * Matrix::colSums(m[is_mt, , drop = FALSE]) / pmax(n_cnt, 1)
    kept   <- n_feat >= CONFIG$seurat_min_features
    pass   <- kept & n_feat > CONFIG$qc_min_features & pct_mt < CONFIG$qc_max_percent_mt
    sum_kept <- sum_kept + Matrix::rowSums(m[, kept, drop = FALSE])
    sum_pass <- sum_pass + Matrix::rowSums(m[, pass, drop = FALSE])
    cells[[k]] <- tibble(nFeature_RNA = n_feat, nCount_RNA = n_cnt, percent_mt = pct_mt,
                         kept_seurat = kept, pass_qc = pass)
  }
  names(sum_kept) <- names(sum_pass) <- genes
  list(cells = bind_rows(cells) %>% mutate(sample = sample, .before = 1),
       sum_kept = sum_kept, sum_pass = sum_pass)
}

read_deseq_table <- function(path) {
  if (!file.exists(path)) stop("Bulk DESeq2 table not found: ", path)
  df <- utils::read.table(path, header = TRUE, sep = "\t", quote = "\"", check.names = FALSE,
                          stringsAsFactors = FALSE, comment.char = "")
  tibble(gene = rownames(df), bulk_baseMean = df$baseMean,
         bulk_log2fc = df$log2FoldChange, bulk_padj = df$padj)
}

# ============================================================================
# 1. LOAD DATA
# ============================================================================
set.seed(CONFIG$random_seed)
for (p in c(CONFIG$sc_cts_path, CONFIG$bulk_counts_path, CONFIG$bulk_tpm_path)) {
  if (!file.exists(p)) stop("Input not found: ", p, " (working dir: ", getwd(), ")")
}

## 1a. scRNA-seq: sum sample x cell-type pseudobulks into all-cell pseudobulks
log_msg("Reading scRNA-seq sample x cell-type counts: %s", CONFIG$sc_cts_path)
cts <- readRDS(CONFIG$sc_cts_path)
# AggregateExpression() turned "_" in orig.ident into "-", so columns are
# "<day>-<genotype>-<rep>_<cell type>"; cell-type names contain no "_".
col_sample <- sub("_.*$", "", colnames(cts))
sc_counts  <- as.matrix(cts %*% Matrix::t(Matrix::fac2sparse(factor(col_sample))))
sc_counts  <- round(sc_counts)
log_msg("Pseudobulk: %d genes x %d samples (from %d sample x cell-type columns)",
        nrow(sc_counts), ncol(sc_counts), ncol(cts))

## 1a-2. Cell-level QC from the Cell Ranger matrices; optionally re-build the
##       pseudobulk from QC-passing cells only (see header "scRNA-seq CELL FILTER")
if (CONFIG$run_cell_qc || CONFIG$sc_cell_filter == "qc") {
  if (!requireNamespace("hdf5r", quietly = TRUE)) stop("Cell QC needs the hdf5r package")
  h5 <- lapply(names(CONFIG$sc_h5_samples), function(s) {
    path <- sprintf(CONFIG$sc_h5_template, CONFIG$sc_h5_samples[[s]])
    if (!file.exists(path)) stop("Cell Ranger matrix not found: ", path)
    log_msg("Cell QC: streaming %s (%s)", s, CONFIG$sc_h5_samples[[s]])
    read_h5_pseudobulk(path, s)
  })
  # orig.ident uses "_", AggregateExpression() column names use "-"
  pb_names  <- gsub("_", "-", names(CONFIG$sc_h5_samples))
  h5_kept   <- do.call(cbind, lapply(h5, `[[`, "sum_kept")); colnames(h5_kept) <- pb_names
  h5_pass   <- do.call(cbind, lapply(h5, `[[`, "sum_pass")); colnames(h5_pass) <- pb_names
  cell_qc   <- bind_rows(lapply(h5, `[[`, "cells"))

  # Check: all Seurat-kept cells must reproduce cts exactly (same genes, same counts).
  # Genes absent from cts were dropped by CreateSeuratObject(min.cells = 1), i.e. zero.
  extra_genes <- setdiff(rownames(h5_kept), rownames(sc_counts))
  genes_ok  <- all(rownames(sc_counts) %in% rownames(h5_kept))
  max_diff  <- if (genes_ok) max(abs(h5_kept[rownames(sc_counts), colnames(sc_counts)] - sc_counts)) else Inf
  reproduces_cts <- genes_ok && max_diff < 0.5 && sum(h5_kept[extra_genes, ]) == 0
  log_msg("Cell QC: all-cell reconstruction reproduces cts exactly: %s (max |diff| = %g)",
          reproduces_cts, max_diff)

  qc_summary <- cell_qc %>%
    group_by(sample) %>%
    summarise(cellranger_cells = n(),
              cells_in_published_object = sum(kept_seurat),
              fail_nFeature = sum(kept_seurat & nFeature_RNA <= CONFIG$qc_min_features),
              fail_percent_mt = sum(kept_seurat & percent_mt >= CONFIG$qc_max_percent_mt),
              pass_qc = sum(pass_qc),
              pct_removed = 100 * (1 - pass_qc / cells_in_published_object),
              median_nFeature = median(nFeature_RNA[kept_seurat]),
              median_percent_mt = median(percent_mt[kept_seurat]),
              .groups = "drop") %>%
    left_join(tibble(sample = names(CONFIG$sc_h5_samples),
                     counts_published = colSums(h5_kept),
                     counts_qc_pass = colSums(h5_pass)), by = "sample") %>%
    mutate(pct_counts_removed = 100 * (1 - counts_qc_pass / counts_published))
  save_tbl(qc_summary, "sc_cell_qc_summary")
  readr::write_csv(cell_qc, file.path(tab_dir, "sc_cell_qc_per_cell.csv.gz"))
  log_msg("Cell QC: %d of %d cells pass nFeature > %d & percent.mt < %g (%.1f%% removed)",
          sum(qc_summary$pass_qc), sum(qc_summary$cells_in_published_object),
          CONFIG$qc_min_features, CONFIG$qc_max_percent_mt,
          100 * (1 - sum(qc_summary$pass_qc) / sum(qc_summary$cells_in_published_object)))

  sample_lab <- function(s) {
    p <- strsplit(s, "_")
    vapply(p, function(x) paste("PB", sub("^day", "d", x[1]), toupper(x[2]), x[3], sep = "_"), "")
  }
  qc_long <- cell_qc %>%
    filter(kept_seurat) %>%
    mutate(label = factor(sample_lab(sample), levels = sort(unique(sample_lab(sample))))) %>%
    select(label, pass_qc, `Detected genes (nFeature_RNA)` = nFeature_RNA,
           `Mitochondrial counts (%)` = percent_mt) %>%
    tidyr::pivot_longer(-c(label, pass_qc), names_to = "metric")
  thresholds <- tibble(metric = c("Detected genes (nFeature_RNA)", "Mitochondrial counts (%)"),
                       y = c(CONFIG$qc_min_features, CONFIG$qc_max_percent_mt))
  p_qc <- ggplot(qc_long, aes(label, value)) +
    geom_violin(fill = "grey85", colour = "grey40", linewidth = 0.25, scale = "width") +
    geom_hline(data = thresholds, aes(yintercept = y), colour = CONFIG$col_collagen,
               linetype = "dashed", linewidth = 0.4) +
    facet_wrap(~ metric, ncol = 1, scales = "free_y") +
    scale_y_continuous(trans = scales::pseudo_log_trans(sigma = 1),
                       breaks = c(0, 1, 5, 15, 50, 100, 500, 1000, 5000, 10000)) +
    labs(x = NULL, y = NULL,
         title = "scRNA-seq per-cell QC (cells in the published object)",
         caption = sprintf("Dashed: intended filter (keep nFeature_RNA > %d and percent.mt < %g). Pseudo-log y axis.",
                           CONFIG$qc_min_features, CONFIG$qc_max_percent_mt)) +
    theme_pub() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  save_plot(p_qc, "FigS0_cell_qc", 7.2, 5)

  if (CONFIG$sc_cell_filter == "qc") {
    if (!reproduces_cts) stop("QC filter needs an exact reconstruction of cts; check sample mapping")
    sc_counts <- round(h5_pass[rownames(sc_counts), colnames(sc_counts)])
    log_msg("Using QC-filtered pseudobulk (nFeature_RNA > %d & percent.mt < %g)",
            CONFIG$qc_min_features, CONFIG$qc_max_percent_mt)
  }
}
log_msg("scRNA-seq cell filter: %s", CONFIG$sc_cell_filter)

sc_meta <- tibble(sample = colnames(sc_counts)) %>%
  tidyr::separate(sample, into = c("day_raw", "genotype", "rep"), sep = "-", remove = FALSE) %>%
  mutate(day = sub("^day", "d", day_raw), platform = "Pseudobulk",
         label = paste("PB", day, toupper(genotype), rep, sep = "_")) %>%
  select(-day_raw)
stopifnot(all(sc_meta$day %in% CONFIG$days), all(sc_meta$genotype %in% CONFIG$genotypes))
colnames(sc_counts) <- sc_meta$label
celltype_composition <- tibble(column = colnames(cts), sample = col_sample,
                               cell_type = sub("^[^_]*_", "", colnames(cts)),
                               total_counts = Matrix::colSums(cts))
# Cell-type labels exist only in the published (unfiltered) object, so this
# table describes all cells in both filter modes.
save_tbl(celltype_composition, "sc_counts_per_sample_celltype")

## 1b. Bulk RNA-seq
read_rsem <- function(path, suffix) {
  x <- utils::read.delim(path, check.names = FALSE)
  rownames(x) <- x$gene_id; x$gene_id <- NULL
  colnames(x) <- sub(paste0("\\|", suffix, "$"), "", colnames(x))
  as.matrix(x)
}
log_msg("Reading bulk RSEM counts and TPM")
bulk_counts_all <- round(read_rsem(CONFIG$bulk_counts_path, "expected_count"))
bulk_tpm_all    <- read_rsem(CONFIG$bulk_tpm_path, "TPM")
stopifnot(identical(colnames(bulk_counts_all), colnames(bulk_tpm_all)))

bulk_meta <- tibble(sample = colnames(bulk_counts_all)) %>%
  mutate(genotype = sub("^Rai1-([a-z]+)_.*$", "\\1", sample),
         day      = sub("^.*_(d[0-9]+)_.*$", "\\1", sample),
         rep      = sub("^.*_(rep[0-9]+)$", "\\1", sample),
         platform = "Bulk",
         label    = paste("Bulk", day, toupper(genotype), rep, sep = "_")) %>%
  filter(day %in% CONFIG$days, !sample %in% CONFIG$bulk_exclude)
log_msg("Bulk samples used (d20/d60, outliers excluded): %d", nrow(bulk_meta))
bulk_counts <- bulk_counts_all[, bulk_meta$sample]; colnames(bulk_counts) <- bulk_meta$label
bulk_tpm    <- bulk_tpm_all[, bulk_meta$sample];    colnames(bulk_tpm)    <- bulk_meta$label

## 1c. Shared gene universe
shared_genes <- sort(intersect(rownames(sc_counts), rownames(bulk_tpm)))
collagen_all <- grep(CONFIG$collagen_regex, shared_genes, value = TRUE)
gene_universe <- tibble(
  source = c("Flex probe set", "Bulk annotation", "Shared", "Collagen (shared)"),
  n_genes = c(nrow(sc_counts), nrow(bulk_tpm), length(shared_genes), length(collagen_all)))
save_tbl(gene_universe, "gene_universe")
save_tbl(tibble(collagen_gene = grep(CONFIG$collagen_regex, union(rownames(sc_counts), rownames(bulk_tpm)),
                                     value = TRUE)) %>%
           mutate(in_flex = collagen_gene %in% rownames(sc_counts),
                  in_bulk = collagen_gene %in% rownames(bulk_tpm)) %>% arrange(collagen_gene),
         "collagen_gene_coverage")
log_msg("Shared genes: %d (collagen: %d)", length(shared_genes), length(collagen_all))

# ============================================================================
# 2. NORMALISATION ON THE SHARED GENE UNIVERSE
# ============================================================================
dge_sc <- edgeR::calcNormFactors(edgeR::DGEList(sc_counts[shared_genes, ]), method = "TMM")
sc_pm  <- edgeR::cpm(dge_sc, log = FALSE)
sc_pm  <- sweep(sc_pm, 2, colSums(sc_pm), "/") * 1e6
bulk_pm <- sweep(bulk_tpm[shared_genes, ], 2, colSums(bulk_tpm[shared_genes, ]), "/") * 1e6

log_pm <- function(m) log2(m + CONFIG$log_pseudocount)
sc_log   <- log_pm(sc_pm)
bulk_log <- log_pm(bulk_pm)
all_meta <- bind_rows(bulk_meta, sc_meta) %>%
  mutate(genotype_lab = geno_label(genotype),
         platform = factor(platform, levels = c("Bulk", "Pseudobulk")),
         day = factor(day, levels = CONFIG$days))
all_log  <- cbind(bulk_log, sc_log)[, all_meta$label]

is_expressed <- rowMeans(sc_pm) >= CONFIG$expressed_min_pm | rowMeans(bulk_pm) >= CONFIG$expressed_min_pm
log_msg("Expressed genes (mean >= %g per million in either platform): %d",
        CONFIG$expressed_min_pm, sum(is_expressed))

# ============================================================================
# 3. EXPRESSION-LEVEL CONCORDANCE (all genes, collagen highlighted)
# ============================================================================
group_mean <- function(logm, meta, day, geno) {
  cols <- meta$label[meta$day == day & meta$genotype == geno]
  rowMeans(logm[, cols, drop = FALSE])
}
conditions <- expand.grid(day = CONFIG$days, genotype = CONFIG$genotypes, stringsAsFactors = FALSE)
expr_long <- bind_rows(
  lapply(seq_len(nrow(conditions)), function(i) {
    d <- conditions$day[i]; g <- conditions$genotype[i]
    tibble(gene = shared_genes, day = d, genotype = g,
           bulk = group_mean(bulk_log, bulk_meta, d, g),
           pseudobulk = group_mean(sc_log, sc_meta, d, g))
  }),
  tibble(gene = shared_genes, day = "All", genotype = "all",
         bulk = rowMeans(bulk_log), pseudobulk = rowMeans(sc_log))
) %>%
  mutate(collagen = gene %in% collagen_all, expressed = gene %in% shared_genes[is_expressed],
         condition = ifelse(day == "All", "All samples", paste(day, toupper(genotype))))

expr_stats <- expr_long %>%
  group_by(condition, day, genotype) %>%
  group_modify(function(df, key) {
    bind_rows(
      cor_stats(df$bulk, df$pseudobulk) %>% mutate(gene_set = "all shared genes"),
      cor_stats(df$bulk[df$expressed], df$pseudobulk[df$expressed]) %>% mutate(gene_set = "expressed genes"),
      cor_stats(df$bulk[df$collagen], df$pseudobulk[df$collagen]) %>% mutate(gene_set = "collagen genes"))
  }) %>% ungroup() %>% relocate(gene_set, .after = genotype)
save_tbl(expr_stats, "expression_correlation_stats")
save_tbl(expr_long %>% select(condition, gene, collagen, expressed,
                              bulk_log2_tpm1 = bulk, pseudobulk_log2_cpm1 = pseudobulk),
         "expression_per_gene")

expr_scatter <- function(df, stats_df, facet = FALSE) {
  lab_df <- stats_df %>% filter(gene_set == "all shared genes")
  lab_df$label <- vapply(seq_len(nrow(lab_df)), function(i) stat_label(lab_df[i, ]), character(1))
  coll <- df %>% filter(collagen)
  to_label <- coll %>% group_by(condition) %>% slice_max(bulk, n = CONFIG$n_label_collagen) %>% ungroup()
  p <- ggplot(df, aes(bulk, pseudobulk)) +
    rast(geom_point(data = df %>% filter(!collagen), color = CONFIG$col_other,
                    size = 0.25, alpha = 0.35, stroke = 0)) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "black",
                linewidth = 0.4, linetype = "dashed") +
    geom_point(data = coll, color = CONFIG$col_collagen, size = 1.1, alpha = 0.9, stroke = 0) +
    geom_text_repel(data = to_label, aes(label = gene), size = 2, color = CONFIG$col_collagen,
                    segment.size = 0.2, min.segment.length = 0, max.overlaps = Inf,
                    box.padding = 0.25, seed = CONFIG$random_seed) +
    geom_text(data = lab_df, aes(x = -Inf, y = Inf, label = indent_lines(label)), hjust = 0, vjust = 1.1,
              size = 2, lineheight = 0.9, inherit.aes = FALSE) +
    labs(x = expression("Bulk RNA-seq, " * log[2] * "(TPM + 1)"),
         y = expression("scRNA-seq pseudobulk, " * log[2] * "(CPM + 1)")) +
    coord_cartesian(clip = "off") + theme_pub()
  if (facet) p <- p + facet_wrap(~condition, nrow = length(CONFIG$days), dir = "v")
  p
}

df_all <- expr_long %>% filter(condition == "All samples")
p_expr_all <- expr_scatter(df_all, expr_stats %>% filter(condition == "All samples")) +
  ggtitle("Pseudobulk vs bulk expression (all genes)") +
  labs(caption = sprintf("Red: collagen genes (n = %d). Dashed: OLS fit, all genes.", sum(df_all$collagen)))
save_plot(p_expr_all, "Fig1A_expression_scatter_all_samples", 3.6, 3.8)

df_cond <- expr_long %>% filter(condition != "All samples") %>%
  mutate(condition = factor(condition, levels = unique(paste(conditions$day, toupper(conditions$genotype)))))
p_expr_cond <- expr_scatter(df_cond, expr_stats %>% filter(condition != "All samples") %>%
                              mutate(condition = factor(condition, levels = levels(df_cond$condition))),
                            facet = TRUE)
save_plot(p_expr_cond, "FigS1_expression_scatter_by_condition", 7.2, 5.0)

# ============================================================================
# 4. LOG2 FOLD-CHANGE CONCORDANCE
# ============================================================================
log_msg("Pseudobulk DESeq2 (design ~ conditions, betaPrior = %s)", CONFIG$deseq_beta_prior)
sc_coldata <- as.data.frame(sc_meta %>% mutate(conditions = factor(paste(day, genotype, sep = "_"))))
rownames(sc_coldata) <- sc_coldata$label
dds_sc <- DESeqDataSetFromMatrix(sc_counts[, sc_coldata$label], colData = sc_coldata, design = ~ conditions)
dds_sc <- DESeq(dds_sc, betaPrior = CONFIG$deseq_beta_prior, quiet = TRUE)

lfc_long <- bind_rows(lapply(CONFIG$days, function(d) {
  bind_rows(lapply(names(CONFIG$contrasts), function(cn) {
    ct  <- CONFIG$contrasts[[cn]]
    key <- paste0(d, "_", cn)
    res <- results(dds_sc, contrast = c("conditions", paste0(d, "_", ct[1]), paste0(d, "_", ct[2])),
                   alpha = CONFIG$padj_cutoff)
    used <- sc_coldata$label[sc_coldata$conditions %in% paste0(d, "_", ct)]
    tot  <- rowSums(sc_counts[rownames(res), used])
    sc_df <- tibble(gene = rownames(res), sc_baseMean = res$baseMean, sc_log2fc = res$log2FoldChange,
                    sc_padj = res$padj, sc_total_counts = tot)
    read_deseq_table(CONFIG$bulk_deseq_paths[[key]]) %>%
      inner_join(sc_df, by = "gene") %>%
      mutate(day = d, contrast = paste(toupper(ct[1]), "vs", toupper(ct[2])),
             collagen = gene %in% collagen_all,
             passes_filter = sc_total_counts >= CONFIG$min_total_counts_deseq &
               is.finite(bulk_log2fc) & is.finite(sc_log2fc))
  }))
})) %>% mutate(panel = paste(day, contrast))
save_tbl(lfc_long, "log2fc_per_gene")
save_tbl(lfc_long %>% filter(collagen) %>% arrange(panel, gene), "log2fc_collagen_genes")

lfc_stats <- lfc_long %>% filter(passes_filter) %>%
  group_by(panel, day, contrast) %>%
  group_modify(function(df, key) {
    sig  <- df %>% filter(!is.na(bulk_padj), bulk_padj < CONFIG$padj_cutoff, sc_log2fc != 0)
    agree <- sum(sign(sig$bulk_log2fc) == sign(sig$sc_log2fc))
    bt <- if (nrow(sig) > 0) binom.test(agree, nrow(sig), p = 0.5) else NULL
    sign_tbl <- tibble(n_bulk_deg = nrow(sig), n_sign_agree = agree,
                       sign_concordance = if (nrow(sig) > 0) agree / nrow(sig) else NA_real_,
                       sign_binom_p = if (!is.null(bt)) bt$p.value else NA_real_)
    coll <- df %>% filter(collagen)
    bind_rows(
      bind_cols(cor_stats(df$bulk_log2fc, df$sc_log2fc, with_ccc = TRUE), sign_tbl) %>%
        mutate(gene_set = "all genes passing filter"),
      cor_stats(coll$bulk_log2fc, coll$sc_log2fc, with_ccc = TRUE) %>%
        mutate(gene_set = "collagen genes",
               mean_lfc_bulk = mean(coll$bulk_log2fc), mean_lfc_pseudobulk = mean(coll$sc_log2fc)))
  }) %>% ungroup() %>% relocate(gene_set, .after = contrast)
save_tbl(lfc_stats, "log2fc_correlation_stats")

lfc_df <- lfc_long %>% filter(passes_filter)
lfc_lab <- lfc_stats %>% filter(gene_set == "all genes passing filter")
lfc_lab$label <- paste0(vapply(seq_len(nrow(lfc_lab)), function(i) stat_label(lfc_lab[i, ], TRUE), character(1)),
                        sprintf("\nsign agree (bulk DEG) = %.0f%%", 100 * lfc_lab$sign_concordance))
lfc_coll <- lfc_df %>% filter(collagen)
lfc_to_label <- lfc_coll %>% group_by(panel) %>%
  slice_max(abs(bulk_log2fc), n = CONFIG$n_label_collagen) %>% ungroup()
lim <- quantile(abs(c(lfc_df$bulk_log2fc, lfc_df$sc_log2fc)), 0.999, na.rm = TRUE)
p_lfc <- ggplot(lfc_df, aes(bulk_log2fc, sc_log2fc)) +
  geom_hline(yintercept = 0, linewidth = 0.2, color = "grey40") +
  geom_vline(xintercept = 0, linewidth = 0.2, color = "grey40") +
  geom_abline(slope = 1, intercept = 0, linewidth = 0.3, color = "grey30", linetype = "dotted") +
  rast(geom_point(data = lfc_df %>% filter(!collagen), color = CONFIG$col_other,
                  size = 0.25, alpha = 0.35, stroke = 0)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, color = "black",
              linewidth = 0.4, linetype = "dashed") +
  geom_point(data = lfc_coll, color = CONFIG$col_collagen, size = 1.1, alpha = 0.9, stroke = 0) +
  geom_text_repel(data = lfc_to_label, aes(label = gene), size = 2, color = CONFIG$col_collagen,
                  segment.size = 0.2, min.segment.length = 0, max.overlaps = Inf, seed = CONFIG$random_seed) +
  geom_text(data = lfc_lab, aes(x = -Inf, y = Inf, label = indent_lines(label)), hjust = 0, vjust = 1.1,
            size = 1.9, lineheight = 0.9, inherit.aes = FALSE) +
  facet_wrap(~panel, nrow = 2) +
  coord_cartesian(xlim = c(-lim, lim), ylim = c(-lim, lim)) +
  labs(x = expression("Bulk RNA-seq " * log[2] * "FC"),
       y = expression("scRNA-seq pseudobulk " * log[2] * "FC"),
       title = "Genotype effects: pseudobulk vs bulk",
       caption = "Red: collagen genes. Dashed: OLS fit; dotted: identity. log2FC from DESeq2 (betaPrior) on both platforms.") +
  theme_pub()
save_plot(p_lfc, "Fig2_log2fc_scatter", 6.0, 6.0)

# ============================================================================
# 5. COLLAGEN-SPECIFIC DEVIATION (expression-matched permutation)
# ============================================================================
log_msg("Collagen residual permutation tests (%d permutations each)", CONFIG$n_perm)
set.seed(CONFIG$random_seed)
coll_expr <- expr_long %>% group_by(condition) %>%
  group_modify(~ gene_set_residual_test(.x$bulk, .x$pseudobulk, .x$collagen, .x$bulk,
                                        CONFIG$n_perm, CONFIG$n_expr_bins)) %>%
  ungroup() %>% mutate(analysis = "expression", fdr_bh = p.adjust(p_perm, "BH"))
coll_lfc <- lfc_df %>% group_by(panel) %>%
  group_modify(~ gene_set_residual_test(.x$bulk_log2fc, .x$sc_log2fc, .x$collagen,
                                        log10(.x$bulk_baseMean + 1), CONFIG$n_perm, CONFIG$n_expr_bins)) %>%
  ungroup() %>% rename(condition = panel) %>%
  mutate(analysis = "log2FC", fdr_bh = p.adjust(p_perm, "BH"))
coll_tests <- bind_rows(coll_expr, coll_lfc) %>% relocate(analysis)
save_tbl(coll_tests, "collagen_residual_permutation_tests")

p_coll_res <- ggplot(coll_tests, aes(mean_residual, condition)) +
  geom_vline(xintercept = 0, linewidth = 0.3, color = "grey40") +
  geom_linerange(aes(xmin = null_mean - 1.96 * null_sd, xmax = null_mean + 1.96 * null_sd),
                 linewidth = 0.8, color = "grey75") +
  geom_point(color = CONFIG$col_collagen, size = 1.8) +
  geom_text(aes(label = sprintf("FDR = %.2g", fdr_bh)), nudge_y = 0.35, size = 1.9) +
  facet_wrap(~analysis, scales = "free", ncol = 2) +
  scale_x_continuous(expand = expansion(mult = 0.2)) +
  coord_cartesian(clip = "off") +
  labs(x = "Mean collagen residual from genome-wide OLS (pseudobulk - fitted)", y = NULL,
       title = "Collagen genes vs genome-wide trend",
       caption = "Grey: 95% range of expression-matched random gene sets (10,000 permutations).") +
  theme_pub()
save_plot(p_coll_res, "Fig4B_collagen_residuals", 6.0, 3.2)

# ============================================================================
# 6. BATCH EFFECTS: COMBINED PCA, CROSS-PLATFORM CORRELATION
# ============================================================================
log_msg("Batch-effect diagnostics")
expr_log <- all_log[is_expressed, ]
pca_df <- function(m, n_hvg) {
  hv <- head(order(apply(m, 1, var), decreasing = TRUE), n_hvg)
  pc <- prcomp(t(m[hv, ]), center = TRUE, scale. = FALSE)
  ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  list(scores = as_tibble(pc$x[, 1:min(10, ncol(pc$x))], rownames = "label") %>% left_join(all_meta, by = "label"),
       var = ve)
}
pca_raw <- pca_df(expr_log, CONFIG$n_hvg)
design_bio <- model.matrix(~ day * genotype, data = all_meta)
expr_corr <- limma::removeBatchEffect(expr_log, batch = all_meta$platform, design = design_bio)
pca_cor <- pca_df(expr_corr, CONFIG$n_hvg)

pca_plot <- function(pp, title) {
  ggplot(pp$scores, aes(PC1, PC2, fill = genotype_lab, shape = platform, color = day)) +
    geom_point(size = 2.2, stroke = 0.8) +
    scale_fill_manual(values = CONFIG$genotype_colors, name = "Genotype") +
    scale_shape_manual(values = CONFIG$platform_shapes, name = "Platform") +
    scale_color_manual(values = CONFIG$day_colors, name = "Day") +
    guides(fill = guide_legend(override.aes = list(shape = 21, color = NA)),
           shape = guide_legend(override.aes = list(fill = "grey80"))) +
    labs(x = sprintf("PC1 (%.1f%%)", pp$var[1]), y = sprintf("PC2 (%.1f%%)", pp$var[2]), title = title) +
    theme_pub()
}
p_pca <- (pca_plot(pca_raw, "Uncorrected") | pca_plot(pca_cor, "Platform effect removed")) +
  plot_layout(guides = "collect") +
  plot_annotation(caption = sprintf("Top %d variable expressed genes; log2(per million + 1).\nCorrection (limma::removeBatchEffect, design ~ day*genotype) is for visualisation only.",
                                    CONFIG$n_hvg),
                  theme = theme(plot.caption = element_text(size = CONFIG$base_size - 1)))
save_plot(p_pca, "Fig3A_combined_PCA", 7.2, 3.2)
save_tbl(bind_rows(pca_raw$scores %>% mutate(version = "uncorrected"),
                   pca_cor$scores %>% mutate(version = "platform_removed")) %>%
           select(version, label, platform, day, genotype, rep, PC1:PC5), "combined_pca_scores")

# Cross-platform correlation on within-platform-centred profiles
centre_within <- function(m, meta) {
  for (pl in levels(meta$platform)) {
    cols <- meta$label[meta$platform == pl]
    m[, cols] <- m[, cols] - rowMeans(m[, cols])
  }
  m
}
cent <- centre_within(expr_log, all_meta)
hv_c <- head(order(apply(cent, 1, var), decreasing = TRUE), CONFIG$n_hvg)
bulk_ord <- all_meta %>% filter(platform == "Bulk") %>% arrange(day, genotype_lab, rep) %>% pull(label)
sc_ord   <- all_meta %>% filter(platform == "Pseudobulk") %>% arrange(day, genotype_lab, rep) %>% pull(label)
xcor <- cor(cent[hv_c, sc_ord], cent[hv_c, bulk_ord], method = "pearson")
save_tbl(as_tibble(xcor, rownames = "pseudobulk_sample"), "cross_platform_correlation_matrix")

best_match <- as_tibble(xcor, rownames = "label") %>%
  pivot_longer(-label, names_to = "bulk_sample", values_to = "r") %>%
  left_join(all_meta %>% select(label, day, genotype), by = "label") %>%
  left_join(all_meta %>% select(bulk_sample = label, bulk_day = day, bulk_genotype = genotype),
            by = "bulk_sample") %>%
  group_by(label, day, genotype, bulk_day, bulk_genotype) %>%
  summarise(mean_r = mean(r), .groups = "drop_last") %>%
  ungroup() %>% group_by(label) %>%
  mutate(rank = rank(-mean_r)) %>%
  summarise(best_bulk_group = paste(bulk_day[rank == 1], toupper(bulk_genotype[rank == 1])),
            best_mean_r = max(mean_r),
            matched_group_r = mean_r[bulk_day == day & bulk_genotype == genotype],
            matched_group_rank = rank[bulk_day == day & bulk_genotype == genotype], .groups = "drop")
save_tbl(best_match, "cross_platform_best_match")

ann_col <- as.data.frame(all_meta %>% filter(platform == "Bulk") %>% select(label, Day = day, Genotype = genotype_lab))
rownames(ann_col) <- ann_col$label; ann_col$label <- NULL
ann_row <- as.data.frame(all_meta %>% filter(platform == "Pseudobulk") %>% select(label, Day = day, Genotype = genotype_lab))
rownames(ann_row) <- ann_row$label; ann_row$label <- NULL
ann_colors <- list(Day = CONFIG$day_colors, Genotype = CONFIG$genotype_colors)
lim_x <- max(abs(xcor))
ph_x <- pheatmap::pheatmap(
  xcor, cluster_rows = FALSE, cluster_cols = FALSE,
  color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(101),
  breaks = seq(-lim_x, lim_x, length.out = 102),
  annotation_col = ann_col[bulk_ord, , drop = FALSE], annotation_row = ann_row[sc_ord, , drop = FALSE],
  annotation_colors = ann_colors, fontsize = CONFIG$base_size - 1, border_color = NA,
  gaps_row = sum(grepl("_d20_", sc_ord)), gaps_col = sum(grepl("_d20_", bulk_ord)),
  main = sprintf("Pearson r, within-platform-centred profiles (top %d variable genes)", CONFIG$n_hvg),
  silent = TRUE)
save_pheatmap(ph_x, "Fig3B_cross_platform_correlation_heatmap", 7.2, 3.8)

# ============================================================================
# 7. WITHIN-PLATFORM VARIANCE DECOMPOSITION (clone / replicate batch)
# ============================================================================
vst_safe <- function(counts) {
  dds <- DESeqDataSetFromMatrix(counts, colData = data.frame(row.names = colnames(counts), x = rep(1, ncol(counts))), design = ~ 1)
  dds <- dds[rowSums(counts(dds)) >= CONFIG$min_total_counts_deseq, ]
  tryCatch(assay(vst(dds, blind = TRUE)),
           error = function(e) assay(varianceStabilizingTransformation(dds, blind = TRUE)))
}
decomp <- function(counts, meta, platform_name) {
  v  <- vst_safe(counts[, meta$label])
  hv <- head(order(apply(v, 1, var), decreasing = TRUE), CONFIG$n_hvg)
  pc <- prcomp(t(v[hv, ]), center = TRUE, scale. = FALSE)
  ve <- 100 * pc$sdev^2 / sum(pc$sdev^2)
  sc <- as_tibble(pc$x, rownames = "label") %>% left_join(meta, by = "label")
  factors <- c(day = "day", genotype = "genotype", `clone / replicate` = "rep")
  bind_rows(lapply(seq_len(CONFIG$n_pcs_decomp), function(k) {
    bind_rows(lapply(names(factors), function(f) {
      fit <- lm(sc[[paste0("PC", k)]] ~ factor(sc[[factors[[f]]]]))
      tibble(platform = platform_name, PC = paste0("PC", k), pct_var = ve[k], factor = f,
             r_squared = summary(fit)$r.squared,
             p = anova(fit)[["Pr(>F)"]][1])
    }))
  }))
}
var_decomp <- bind_rows(
  decomp(bulk_counts, bulk_meta %>% mutate(label = label), "Bulk"),
  decomp(sc_counts, sc_meta, "Pseudobulk"))
save_tbl(var_decomp, "within_platform_pc_variance_decomposition")

p_decomp <- var_decomp %>%
  mutate(PC_lab = sprintf("%s\n(%.0f%%)", PC, pct_var),
         factor = factor(factor, levels = c("day", "genotype", "clone / replicate"))) %>%
  ggplot(aes(PC_lab, factor, fill = r_squared)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", r_squared)), size = 2) +
  scale_fill_gradient(low = "white", high = "#2C6FAC", limits = c(0, 1), name = expression(R^2)) +
  facet_wrap(~platform, scales = "free_x") +
  labs(x = NULL, y = NULL, title = "Variance of top PCs explained by each factor",
       caption = "R\u00b2 of lm(PC ~ factor), factors fitted separately; vst, top variable genes.") +
  theme_pub() + theme(axis.line = element_blank(), axis.ticks = element_blank())
save_plot(p_decomp, "Fig3C_pc_variance_decomposition", 6.0, 2.4)

# ============================================================================
# 8. REPLICATE VARIABILITY (replicate correlation, edgeR BCV)
# ============================================================================
rep_cor <- bind_rows(lapply(levels(all_meta$platform), function(pl) {
  m_pl <- all_meta %>% filter(platform == pl)
  pm   <- if (pl == "Bulk") bulk_pm else sc_pm
  keep <- rowMeans(pm) >= CONFIG$expressed_min_pm
  lm_  <- log_pm(pm[keep, ])
  m_pl %>% group_by(day, genotype_lab) %>%
    group_modify(function(g, key) {
      if (nrow(g) < 2) return(tibble())
      prs <- combn(g$label, 2)
      tibble(rep_a = prs[1, ], rep_b = prs[2, ],
             pearson_r = apply(prs, 2, function(p) cor(lm_[, p[1]], lm_[, p[2]])))
    }) %>% ungroup() %>% mutate(platform = pl, n_genes = sum(keep))
}))
save_tbl(rep_cor, "replicate_pairwise_correlation")

bcv <- bind_rows(lapply(levels(all_meta$platform), function(pl) {
  bind_rows(lapply(CONFIG$days, function(d) {
    m <- all_meta %>% filter(platform == pl, day == d)
    cnt <- (if (pl == "Bulk") bulk_counts else sc_counts)[shared_genes, m$label]
    y <- edgeR::DGEList(cnt, group = m$genotype)
    y <- y[edgeR::filterByExpr(y), , keep.lib.sizes = FALSE]
    y <- edgeR::calcNormFactors(y)
    y <- edgeR::estimateDisp(y, model.matrix(~ genotype, data = m))
    tibble(platform = pl, day = d, n_samples = ncol(y), n_genes = nrow(y),
           common_bcv = sqrt(y$common.dispersion),
           median_tagwise_bcv = median(sqrt(y$tagwise.dispersion)),
           median_tagwise_bcv_collagen = median(sqrt(y$tagwise.dispersion[rownames(y) %in% collagen_all])))
  }))
}))
save_tbl(bcv, "edgeR_bcv")

p_rep <- ggplot(rep_cor, aes(genotype_lab, pearson_r, color = platform, shape = platform)) +
  geom_point(position = position_dodge(width = 0.6), size = 1.6) +
  facet_wrap(~day) +
  scale_color_manual(values = c(Bulk = "black", Pseudobulk = CONFIG$col_collagen), name = NULL) +
  scale_shape_manual(values = c(Bulk = 16, Pseudobulk = 17), name = NULL) +
  labs(x = NULL, y = "Pearson r between replicates\n(log2 per million, expressed genes)",
       title = "Replicate concordance") + theme_pub()
p_bcv <- ggplot(bcv, aes(day, common_bcv, fill = platform)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  geom_text(aes(label = sprintf("%.2f", common_bcv)), position = position_dodge(width = 0.7),
            vjust = -0.3, size = 2) +
  scale_fill_manual(values = c(Bulk = "grey30", Pseudobulk = CONFIG$col_collagen), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Common BCV (edgeR)", title = "Biological coefficient of variation") + theme_pub()
save_plot(p_rep | p_bcv, "Fig3D_replicate_variability", 7.2, 2.8)

# ============================================================================
# 9. COLLAGEN GENES ACROSS ALL SAMPLES (variability view)
# ============================================================================
coll_show <- collagen_all[collagen_all %in% shared_genes[is_expressed]]
z_within <- function(m) t(scale(t(m)))
coll_z <- cbind(z_within(bulk_log[coll_show, bulk_ord, drop = FALSE]),
                z_within(sc_log[coll_show, sc_ord, drop = FALSE]))
# Genes constant within a platform give NaN z-scores and cannot be clustered.
coll_z <- coll_z[rowSums(!is.finite(coll_z)) == 0, , drop = FALSE]
ann_all <- as.data.frame(all_meta %>% select(label, Platform = platform, Day = day, Genotype = genotype_lab))
rownames(ann_all) <- ann_all$label; ann_all$label <- NULL
ph_c <- pheatmap::pheatmap(
  coll_z, cluster_cols = FALSE, cluster_rows = TRUE, clustering_distance_rows = "correlation",
  color = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(101),
  breaks = seq(-2.5, 2.5, length.out = 102),
  annotation_col = ann_all[c(bulk_ord, sc_ord), , drop = FALSE],
  annotation_colors = c(ann_colors, list(Platform = c(Bulk = "grey30", Pseudobulk = CONFIG$col_collagen))),
  gaps_col = length(bulk_ord), show_colnames = TRUE, fontsize = CONFIG$base_size - 1,
  border_color = NA, main = "Collagen genes: z-score within platform (log2 per million)",
  silent = TRUE)
save_pheatmap(ph_c, "Fig4A_collagen_heatmap_all_samples", 7.2, 4.6)

# ============================================================================
# 10. SESSION INFO
# ============================================================================
writeLines(c(capture.output(print(CONFIG)), "", capture.output(sessionInfo())),
           file.path(CONFIG$out_dir, "run_parameters_and_sessionInfo.txt"))
log_msg("Done. Outputs in %s", CONFIG$out_dir)
