# GO enrichment (Biological Process, Molecular Function, Cellular Component) for genes whose
# chromatin accessibility and expression change in the OPPOSITE direction in KO vs WT
# ("Discordant" genes from atac_rna_correlation.R). Runs locally in RStudio.
#
# Input: ATAC_RNA_merged_all_nearest_gene.csv, written by atac_rna_correlation.R.
#   For distance-restricted versions use the other merged tables from that script, e.g.
#   ATAC_RNA_merged_within_50kb.csv or ATAC_RNA_merged_promoter_3kb.csv.
#
# Install once if needed:
#   BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"))
#   install.packages("ggplot2")

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
})

# ======================= settings =======================
merged_file <- "results/ATACseq/RNA_correlation/ATAC_RNA_merged_all_nearest_gene.csv"
outdir      <- "results/ATACseq/RNA_correlation/GO_discordant"

ontologies  <- c(BP = "Biological Process", MF = "Molecular Function", CC = "Cellular Component")
padj_cut    <- 0.05
min_gs      <- 5      # smallest / largest GO term (genes in the background) to test
max_gs      <- 500
simplify_cutoff <- 0.8 # merge GO terms that are >80% semantically similar (keeps the best p-value)
n_show      <- 12      # terms shown per ontology in each figure
split_direction <- TRUE  # also run: more open + up in KO, and less open + down in KO
# ========================================================

merged_file <- path.expand(merged_file)
outdir      <- path.expand(outdir)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

m <- read.csv(merged_file, stringsAsFactors = FALSE)
miss <- setdiff(c("SYMBOL", "cat", "atac_lfc", "rna_lfc"), names(m))
if (length(miss)) stop("Merged table is missing: ", paste(miss, collapse = ", "), call. = FALSE)

# Background = every gene that could have been called concordant: it had an ATAC peak
# assigned AND was measured in the RNA-seq. Using this (not the whole genome) stops
# the result from simply reflecting "genes with accessible chromatin and expression".
universe <- unique(m$SYMBOL)
conc     <- m[m$cat == "Discordant", ]

sets <- list(all = list(label = "All discordant genes", genes = unique(conc$SYMBOL)))
if (split_direction) {
  sets$up   <- list(label = "Discordant: more open + up in KO", # note this is actually more open + down in KO
                    genes = unique(conc$SYMBOL[conc$atac_lfc > 0]))
  sets$down <- list(label = "Discordant: less open + down in KO", # note this is actually less open + up in KO
                    genes = unique(conc$SYMBOL[conc$atac_lfc < 0]))
}

message("Background: ", length(universe), " genes | discordant: ", length(sets$all$genes),
        if (split_direction) paste0(" (", length(sets$up$genes), " up, ", length(sets$down$genes), " down)"))
write.csv(data.frame(SYMBOL = conc$SYMBOL, atac_lfc = conc$atac_lfc, rna_lfc = conc$rna_lfc,
                     direction = ifelse(conc$atac_lfc > 0, "more open + up", "less open + down")),
          file.path(outdir, "discordant_genes.csv"), row.names = FALSE)

# ---------------- enrichment ----------------
go_one <- function(genes, ont) {
  if (length(genes) < 5) return(NULL)
  res <- tryCatch(
    enrichGO(gene = genes, universe = universe, OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
             ont = ont, pAdjustMethod = "BH", pvalueCutoff = padj_cut, qvalueCutoff = padj_cut,
             minGSSize = min_gs, maxGSSize = max_gs),
    error = function(e) { message("  enrichGO failed for ", ont, ": ", conditionMessage(e)); NULL })
  if (is.null(res) || nrow(as.data.frame(res)) == 0) return(NULL)
  res <- tryCatch(clusterProfiler::simplify(res, cutoff = simplify_cutoff, by = "p.adjust", select_fun = min),
                  error = function(e) res)
  d <- as.data.frame(res)
  ratio <- function(x) vapply(strsplit(x, "/"), function(v) as.numeric(v[1]) / as.numeric(v[2]), numeric(1))
  d$GeneRatioNum     <- ratio(d$GeneRatio)
  d$FoldEnrichment   <- d$GeneRatioNum / ratio(d$BgRatio)
  d$ontology         <- ont
  d
}

wrap <- function(x, w = 42) vapply(strwrap(x, width = w, simplify = FALSE), paste, character(1), collapse = "\n")

for (tag in names(sets)) {
  s <- sets[[tag]]
  message("\n== ", s$label, " (", length(s$genes), " genes) ==")
  res <- lapply(names(ontologies), function(o) go_one(s$genes, o))
  names(res) <- names(ontologies)
  for (o in names(res)) {
    message("  ", ontologies[[o]], ": ", if (is.null(res[[o]])) 0 else nrow(res[[o]]),
            " enriched terms (adj. p < ", padj_cut, ", after simplify)")
    if (!is.null(res[[o]]))
      write.csv(res[[o]], file.path(outdir, sprintf("GO_%s_%s.csv", tag, o)), row.names = FALSE)
  }
  all_res <- do.call(rbind, res[!vapply(res, is.null, logical(1))])
  if (is.null(all_res)) { message("  nothing enriched; no figure made"); next }
  write.csv(all_res, file.path(outdir, sprintf("GO_%s_all_ontologies.csv", tag)), row.names = FALSE)
  
  # ---- figure: top terms per ontology, one panel each ----
  top <- do.call(rbind, lapply(split(all_res, all_res$ontology), function(d) head(d[order(d$p.adjust), ], n_show)))
  top$ontology <- factor(unname(ontologies[top$ontology]), levels = unname(ontologies))
  top <- top[order(top$ontology, top$GeneRatioNum), ]
  top$term <- factor(wrap(top$Description), levels = unique(wrap(top$Description)))
  p <- ggplot(top, aes(GeneRatioNum, term, size = Count, color = p.adjust)) +
    geom_point() +
    facet_grid(ontology ~ ., scales = "free_y", space = "free_y") +
    scale_color_gradient(low = "firebrick", high = "steelblue", name = "adj. p") +
    scale_size_continuous(range = c(1.8, 6)) +
    labs(x = "Gene ratio", y = NULL, title = paste0("GO enrichment: ", s$label),
         subtitle = sprintf("%d genes vs %d-gene background; top %d terms per ontology",
                            length(s$genes), length(universe), n_show)) +
    theme_bw() + theme(strip.text.y = element_text(angle = 0))
  ggsave(file.path(outdir, sprintf("GO_%s_dotplot.pdf", tag)), p,
         width = 8, height = max(4, 1.6 + 0.3 * nrow(top)))
}

message("\nDone. Outputs in: ", outdir)
