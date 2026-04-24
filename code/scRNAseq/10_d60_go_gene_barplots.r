#########################################################################################
### Day 20 pseudobulk bar plot for postsynaptic genes across NPC and IPC cell types
#########################################################################################

library(tidyverse)
library(DESeq2)

postsynaptic_genes <- c("GRIN3A", "SLITRK2", "MAGEE1", "ITGA8", "CBLN1")
cell_types         <- c("NPC-G2/M", "NPC-G1", "Mes-WT", "Mes-KO", "IPC")
genotype_levels    <- c("wt", "het", "ko")
genotype_labels    <- c("wt" = "WT", "het" = "HET", "ko" = "KO")
genotype_colors    <- c("WT" = "#7b7b7b", "HET" = "#f1d174", "KO" = "#bd6565")

# Load Day 20 data
dds_d20         <- readRDS("data/scRNAseq/iwaselab/dds_d20")
norm_counts_d20 <- counts(dds_d20, normalized = TRUE)
coldata_d20     <- as.data.frame(colData(dds_d20)) %>%
	rownames_to_column("sample_id")

output_dir <- "results/scRNAseq/7_pseudobulk/select_gene_barplots"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Prepare expression data
postsynaptic_expr_long <- norm_counts_d20[postsynaptic_genes, , drop = FALSE] %>%
	as.data.frame() %>%
	rownames_to_column("gene") %>%
	pivot_longer(
		cols = -gene,
		names_to = "sample_id",
		values_to = "expression"
	) %>%
	left_join(coldata_d20, by = "sample_id") %>%
	filter(cluster_d20 %in% cell_types, genotype_d20 %in% genotype_levels) %>%
	mutate(
		gene      = factor(gene, levels = postsynaptic_genes),
		genotype  = factor(genotype_labels[genotype_d20], levels = c("WT", "HET", "KO")),
		cell_type = factor(cluster_d20, levels = cell_types)
	)

missing_genes <- setdiff(postsynaptic_genes, unique(as.character(postsynaptic_expr_long$gene)))
if (length(missing_genes) > 0) {
	stop(sprintf("Missing postsynaptic genes in day20 normalized counts: %s", paste(missing_genes, collapse = ", ")))
}

postsynaptic_bar_df <- postsynaptic_expr_long %>%
	group_by(gene, cell_type, genotype) %>%
	summarise(mean_expression = mean(expression), .groups = "drop")

plot_day20_postsynaptic_genes_norm_counts <- ggplot(
	postsynaptic_bar_df,
	aes(x = cell_type, y = mean_expression, fill = genotype)
) +
	geom_col(width = 0.75, position = position_dodge(width = 0.8)) +
	geom_point(
		data = postsynaptic_expr_long,
		aes(x = cell_type, y = expression, fill = genotype),
		inherit.aes = FALSE,
		position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
		size = 1.8,
		shape = 21,
		stroke = 0.35,
		color = "black"
	) +
	facet_wrap(~ gene, nrow = 1, scales = "free_y") +
	scale_fill_manual(values = genotype_colors, drop = FALSE) +
	labs(
		title = "Day 20 postsynaptic genes pseudobulk expression (NPC-G2/M, NPC-G1, Mes-WT, Mes-KO, IPC)",
		subtitle = "Bars show genotype mean normalized counts; points show replicate pseudobulk samples",
		x = NULL,
		y = "Normalized counts"
	) +
	theme_classic(base_size = 11) +
	theme(
		legend.position = "top",
		legend.title = element_blank(),
		strip.background = element_blank(),
		strip.text = element_text(face = "bold"),
		axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
		plot.title = element_text(face = "bold")
	)

ggsave(
	filename = file.path(output_dir, "day20_postsynaptic_genes_npc_ipc_normalized_counts_barplot.pdf"),
	plot = plot_day20_postsynaptic_genes_norm_counts,
	width = 16,
	height = 4.5
)

ggsave(
	filename = file.path(output_dir, "day20_postsynaptic_genes_npc_ipc_normalized_counts_barplot.png"),
	plot = plot_day20_postsynaptic_genes_norm_counts,
	width = 16,
	height = 4.5,
	dpi = 300
)

plot_day20_postsynaptic_genes_norm_counts

####################################################

### Day 20 collagen/ECM genes across NPC and IPC cell types

collagen_genes <- c("COL22A1", "C1QTNF7", "COL6A3", "EMILIN1", "CCBE1")

collagen_expr_long <- norm_counts_d20[collagen_genes, , drop = FALSE] %>%
	as.data.frame() %>%
	rownames_to_column("gene") %>%
	pivot_longer(
		cols = -gene,
		names_to = "sample_id",
		values_to = "expression"
	) %>%
	left_join(coldata_d20, by = "sample_id") %>%
	filter(cluster_d20 %in% cell_types, genotype_d20 %in% genotype_levels) %>%
	mutate(
		gene      = factor(gene, levels = collagen_genes),
		genotype  = factor(genotype_labels[genotype_d20], levels = c("WT", "HET", "KO")),
		cell_type = factor(cluster_d20, levels = cell_types)
	)

missing_collagen_genes <- setdiff(collagen_genes, unique(as.character(collagen_expr_long$gene)))
if (length(missing_collagen_genes) > 0) {
	stop(sprintf("Missing collagen genes in day20 normalized counts: %s", paste(missing_collagen_genes, collapse = ", ")))
}

collagen_bar_df <- collagen_expr_long %>%
	group_by(gene, cell_type, genotype) %>%
	summarise(mean_expression = mean(expression), .groups = "drop")

plot_day20_collagen_genes_norm_counts <- ggplot(
	collagen_bar_df,
	aes(x = cell_type, y = mean_expression, fill = genotype)
) +
	geom_col(width = 0.75, position = position_dodge(width = 0.8)) +
	geom_point(
		data = collagen_expr_long,
		aes(x = cell_type, y = expression, fill = genotype),
		inherit.aes = FALSE,
		position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
		size = 1.8,
		shape = 21,
		stroke = 0.35,
		color = "black"
	) +
	facet_wrap(~ gene, nrow = 1, scales = "free_y") +
	scale_fill_manual(values = genotype_colors, drop = FALSE) +
	labs(
		title = "Day 20 collagen/ECM genes pseudobulk expression (NPC-G2/M, NPC-G1, Mes-WT, Mes-KO, IPC)",
		subtitle = "Bars show genotype mean normalized counts; points show replicate pseudobulk samples",
		x = NULL,
		y = "Normalized counts"
	) +
	theme_classic(base_size = 11) +
	theme(
		legend.position = "top",
		legend.title = element_blank(),
		strip.background = element_blank(),
		strip.text = element_text(face = "bold"),
		axis.text.x = element_text(face = "bold", angle = 45, hjust = 1),
		plot.title = element_text(face = "bold")
	)

ggsave(
	filename = file.path(output_dir, "day20_collagen_genes_npc_ipc_normalized_counts_barplot.pdf"),
	plot = plot_day20_collagen_genes_norm_counts,
	width = 16,
	height = 4.5
)

ggsave(
	filename = file.path(output_dir, "day20_collagen_genes_npc_ipc_normalized_counts_barplot.png"),
	plot = plot_day20_collagen_genes_norm_counts,
	width = 16,
	height = 4.5,
	dpi = 300
)

plot_day20_collagen_genes_norm_counts

