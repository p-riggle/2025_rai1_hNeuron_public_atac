
#########################################################################################
### Day 60 neuron pseudobulk bar plot for selected genes
#########################################################################################

library(tidyverse)
library(DESeq2)

genes_of_interest <- c("NRGN", "ITGA8", "SHISA9", "GRM5", "GRIA4")
genotype_levels <- c("wt", "het", "ko")
genotype_labels <- c("wt" = "WT", "het" = "HET", "ko" = "KO")
genotype_colors <- c("WT" = "#7b7b7b", "HET" = "#f1d174", "KO" = "#bd6565")

dds_d60 <- readRDS("data/scRNAseq/iwaselab/dds_d60")
norm_counts_d60 <- counts(dds_d60, normalized = TRUE)

coldata_d60 <- as.data.frame(colData(dds_d60)) %>%
	rownames_to_column("sample_id")

dds_d20 <- readRDS("data/scRNAseq/iwaselab/dds_d20")
norm_counts_d20 <- counts(dds_d20, normalized = TRUE)

coldata_d20 <- as.data.frame(colData(dds_d20)) %>%
	rownames_to_column("sample_id")

norm_expr_long <- norm_counts_d60[genes_of_interest, , drop = FALSE] %>%
	as.data.frame() %>%
	rownames_to_column("gene") %>%
	pivot_longer(
		cols = -gene,
		names_to = "sample_id",
		values_to = "expression"
	) %>%
	left_join(coldata_d60, by = "sample_id") %>%
	filter(cluster_d60 == "Neuron", genotype_d60 %in% genotype_levels) %>%
	mutate(
		gene = factor(gene, levels = genes_of_interest),
		genotype = factor(genotype_labels[genotype_d60], levels = c("WT", "HET", "KO"))
	)

missing_norm_genes <- setdiff(genes_of_interest, unique(as.character(norm_expr_long$gene)))
if (length(missing_norm_genes) > 0) {
	stop(sprintf("Missing requested genes in day60 normalized counts: %s", paste(missing_norm_genes, collapse = ", ")))
}

norm_bar_df <- norm_expr_long %>%
	group_by(gene, genotype) %>%
	summarise(
		mean_expression = mean(expression),
		.groups = "drop"
	)

plot_day60_neuron_genes_norm_counts <- ggplot(norm_bar_df, aes(x = genotype, y = mean_expression, fill = genotype)) +
	geom_col(width = 0.72) +
	geom_point(
		data = norm_expr_long,
		aes(x = genotype, y = expression),
		inherit.aes = FALSE,
		position = position_jitter(width = 0.12, height = 0),
		size = 2.2,
		shape = 21,
		stroke = 0.35,
		fill = "white",
		color = "black"
	) +
	facet_wrap(~ gene, scales = "free_y", nrow = 1) +
	scale_fill_manual(values = genotype_colors, drop = FALSE) +
	labs(
		title = "Day 60 neuron pseudobulk expression",
		subtitle = "Bars show genotype mean normalized counts; points show replicate pseudobulk samples",
		x = NULL,
		y = "Normalized counts"
	) +
	theme_classic(base_size = 11) +
	theme(
		legend.position = "none",
		strip.background = element_blank(),
		strip.text = element_text(face = "bold"),
		axis.text.x = element_text(face = "bold"),
		plot.title = element_text(face = "bold")
	)

output_dir <- "results/scRNAseq/7_pseudobulk/select_gene_barplots"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

ggsave(
	filename = file.path(output_dir, "day60_neuron_top5_genes_normalized_counts_barplot.pdf"),
	plot = plot_day60_neuron_genes_norm_counts,
	width = 11,
	height = 3.8
)

ggsave(
	filename = file.path(output_dir, "day60_neuron_top5_genes_normalized_counts_barplot.png"),
	plot = plot_day60_neuron_genes_norm_counts,
	width = 11,
	height = 3.8,
	dpi = 300
)

plot_day60_neuron_genes_norm_counts

####################################################

### Day 60 gRG pseudobulk bar plot for selected genes

grg_genes_of_interest <- c("MT-ND5", "MT-ND2", "MT-ATP6", "MT-ND1", "MT-CO2")

grg_norm_expr_long <- norm_counts_d60[grg_genes_of_interest, , drop = FALSE] %>%
	as.data.frame() %>%
	rownames_to_column("gene") %>%
	pivot_longer(
		cols = -gene,
		names_to = "sample_id",
		values_to = "expression"
	) %>%
	left_join(coldata_d60, by = "sample_id") %>%
	filter(cluster_d60 == "gRG", genotype_d60 %in% genotype_levels) %>%
	mutate(
		gene = factor(gene, levels = grg_genes_of_interest),
		genotype = factor(genotype_labels[genotype_d60], levels = c("WT", "HET", "KO"))
	)

missing_grg_genes <- setdiff(grg_genes_of_interest, unique(as.character(grg_norm_expr_long$gene)))
if (length(missing_grg_genes) > 0) {
	stop(sprintf("Missing requested genes in day60 gRG normalized counts: %s", paste(missing_grg_genes, collapse = ", ")))
}

grg_norm_bar_df <- grg_norm_expr_long %>%
	group_by(gene, genotype) %>%
	summarise(
		mean_expression = mean(expression),
		.groups = "drop"
	)

plot_day60_grg_genes_norm_counts <- ggplot(grg_norm_bar_df, aes(x = genotype, y = mean_expression, fill = genotype)) +
	geom_col(width = 0.72) +
	geom_point(
		data = grg_norm_expr_long,
		aes(x = genotype, y = expression),
		inherit.aes = FALSE,
		position = position_jitter(width = 0.12, height = 0),
		size = 2.2,
		shape = 21,
		stroke = 0.35,
		fill = "white",
		color = "black"
	) +
	facet_wrap(~ gene, scales = "free_y", nrow = 1) +
	scale_fill_manual(values = genotype_colors, drop = FALSE) +
	labs(
		title = "Day 60 gRG pseudobulk expression",
		subtitle = "Bars show genotype mean normalized counts; points show replicate pseudobulk samples",
		x = NULL,
		y = "Normalized counts"
	) +
	theme_classic(base_size = 11) +
	theme(
		legend.position = "none",
		strip.background = element_blank(),
		strip.text = element_text(face = "bold"),
		axis.text.x = element_text(face = "bold"),
		plot.title = element_text(face = "bold")
	)

ggsave(
	filename = file.path(output_dir, "day60_gRG_top5_genes_normalized_counts_barplot.pdf"),
	plot = plot_day60_grg_genes_norm_counts,
	width = 11,
	height = 3.8
)

ggsave(
	filename = file.path(output_dir, "day60_gRG_top5_genes_normalized_counts_barplot.png"),
	plot = plot_day60_grg_genes_norm_counts,
	width = 11,
	height = 3.8,
	dpi = 300
)

plot_day60_grg_genes_norm_counts

####################################################

### Day 20 NPC-G1 pseudobulk bar plot for DNA replication genes

npcg1_genes_of_interest <- c("MCM2", "MCM4", "MCM5", "MCM6", "MCM7")

npcg1_norm_expr_long <- norm_counts_d20[npcg1_genes_of_interest, , drop = FALSE] %>%
	as.data.frame() %>%
	rownames_to_column("gene") %>%
	pivot_longer(
		cols = -gene,
		names_to = "sample_id",
		values_to = "expression"
	) %>%
	left_join(coldata_d20, by = "sample_id") %>%
	filter(cluster_d20 == "NPC-G1", genotype_d20 %in% genotype_levels) %>%
	mutate(
		gene = factor(gene, levels = npcg1_genes_of_interest),
		genotype = factor(genotype_labels[genotype_d60], levels = c("WT", "HET", "KO"))
	)

missing_npcg1_genes <- setdiff(npcg1_genes_of_interest, unique(as.character(npcg1_norm_expr_long$gene)))
if (length(missing_npcg1_genes) > 0) {
	stop(sprintf("Missing requested genes in day20 NPC-G1 normalized counts: %s", paste(missing_npcg1_genes, collapse = ", ")))
}

npcg1_norm_bar_df <- npcg1_norm_expr_long %>%
	group_by(gene, genotype) %>%
	summarise(
		mean_expression = mean(expression),
		.groups = "drop"
	)

plot_day20_npcg1_genes_norm_counts <- ggplot(npcg1_norm_bar_df, aes(x = genotype, y = mean_expression, fill = genotype)) +
	geom_col(width = 0.72) +
	geom_point(
		data = npcg1_norm_expr_long,
		aes(x = genotype, y = expression),
		inherit.aes = FALSE,
		position = position_jitter(width = 0.12, height = 0),
		size = 2.2,
		shape = 21,
		stroke = 0.35,
		fill = "white",
		color = "black"
	) +
	facet_wrap(~ gene, scales = "free_y", nrow = 1) +
	scale_fill_manual(values = genotype_colors, drop = FALSE) +
	labs(
		title = "Day 20 NPC-G1 pseudobulk expression",
		subtitle = "Bars show genotype mean normalized counts; points show replicate pseudobulk samples",
		x = NULL,
		y = "Normalized counts"
	) +
	theme_classic(base_size = 11) +
	theme(
		legend.position = "none",
		strip.background = element_blank(),
		strip.text = element_text(face = "bold"),
		axis.text.x = element_text(face = "bold"),
		plot.title = element_text(face = "bold")
	)

ggsave(
	filename = file.path(output_dir, "day20_NPC-G1_MCM_genes_normalized_counts_barplot.pdf"),
	plot = plot_day20_npcg1_genes_norm_counts,
	width = 11,
	height = 3.8
)

ggsave(
	filename = file.path(output_dir, "day20_NPC-G1_MCM_genes_normalized_counts_barplot.png"),
	plot = plot_day20_npcg1_genes_norm_counts,
	width = 11,
	height = 3.8,
	dpi = 300
)

plot_day20_npcg1_genes_norm_counts

