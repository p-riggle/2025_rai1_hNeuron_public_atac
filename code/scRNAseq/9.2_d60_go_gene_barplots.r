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

missing_requested_genes <- setdiff(genes_of_interest, rownames(norm_counts_d60))
if (length(missing_requested_genes) > 0) {
	stop(sprintf("Missing requested genes in day60 Neuron normalized counts: %s", paste(missing_requested_genes, collapse = ", ")))
}

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

plot_day60_neuron_genes_norm_counts

####################################################

### Day 60 gRG pseudobulk bar plot for selected genes

grg_genes_of_interest <- c("MT-ND5", "MT-ND2", "MT-ATP6", "MT-ND1", "MT-CO2")

missing_requested_genes <- setdiff(grg_genes_of_interest, rownames(norm_counts_d60))
if (length(missing_requested_genes) > 0) {
	stop(sprintf("Missing requested genes in day60 gRG normalized counts: %s", paste(missing_requested_genes, collapse = ", ")))
}

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

plot_day60_grg_genes_norm_counts

####################################################
####  Mitochondrial Unfolded Protein Response (UPRmt) genes: HSPD1, HSPE1, CLPP, ATF4, ATF5
############################################

uprmt_genes_of_interest <- c("HSPD1", "HSPE1", "CLPP", "ATF4", "ATF5")

missing_requested_genes <- setdiff(uprmt_genes_of_interest, rownames(norm_counts_d60))
if (length(missing_requested_genes) > 0) {
	stop(sprintf("Missing requested genes in day60 gRG normalized counts: %s", paste(missing_requested_genes, collapse = ", ")))
}

uprmt_norm_expr_long <- norm_counts_d60[uprmt_genes_of_interest, , drop = FALSE] %>%
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
		gene = factor(gene, levels = uprmt_genes_of_interest),
		genotype = factor(genotype_labels[genotype_d60], levels = c("WT", "HET", "KO"))
	)

missing_uprmt_genes <- setdiff(uprmt_genes_of_interest, unique(as.character(uprmt_norm_expr_long$gene)))
if (length(missing_uprmt_genes) > 0) {
	stop(sprintf("Missing requested genes in day60 gRG normalized counts: %s", paste(missing_uprmt_genes, collapse = ", ")))
}

uprmt_norm_bar_df <- uprmt_norm_expr_long %>%
	group_by(gene, genotype) %>%
	summarise(
		mean_expression = mean(expression),
		.groups = "drop"
	)

plot_day60_grg_uprmt_genes_norm_counts <- ggplot(uprmt_norm_bar_df, aes(x = genotype, y = mean_expression, fill = genotype)) +
	geom_col(width = 0.72) +
	geom_point(
		data = uprmt_norm_expr_long,
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
		title = "Day 60 gRG UPRmt gene expression",
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
	filename = file.path(output_dir, "day60_gRG_uprmt_genes_normalized_counts_barplot.pdf"),
	plot = plot_day60_grg_uprmt_genes_norm_counts,
	width = 11,
	height = 3.8
)

plot_day60_grg_uprmt_genes_norm_counts

#############################################
#### ROS Defense genes: SOD1, SOD2 and mitophagy genes: PINK1, PRKN, BNIP3, BNIP3L
###############################################

ros_defense_genes_of_interest <- c("SOD1", "SOD2")

missing_requested_genes <- setdiff(ros_defense_genes_of_interest, rownames(norm_counts_d60))
if (length(missing_requested_genes) > 0) {
	stop(sprintf("Missing requested genes in day60 gRG normalized counts: %s", paste(missing_requested_genes, collapse = ", ")))
}

ros_defense_norm_expr_long <- norm_counts_d60[ros_defense_genes_of_interest, , drop = FALSE] %>%
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
		gene = factor(gene, levels = ros_defense_genes_of_interest),
		genotype = factor(genotype_labels[genotype_d60], levels = c("WT", "HET", "KO"))
	)

missing_ros_defense_genes <- setdiff(ros_defense_genes_of_interest, unique(as.character(ros_defense_norm_expr_long$gene)))
if (length(missing_ros_defense_genes) > 0) {
	stop(sprintf("Missing requested genes in day60 gRG normalized counts: %s", paste(missing_ros_defense_genes, collapse = ", ")))
}

ros_defense_norm_bar_df <- ros_defense_norm_expr_long %>%
	group_by(gene, genotype) %>%
	summarise(
		mean_expression = mean(expression),
		.groups = "drop"
	)

plot_day60_grg_ros_defense_genes_norm_counts <- ggplot(ros_defense_norm_bar_df, aes(x = genotype, y = mean_expression, fill = genotype)) +
	geom_col(width = 0.72) +
	geom_point(
		data = ros_defense_norm_expr_long,
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
		title = "Day 60 gRG ROS defense gene expression",
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

plot_day60_grg_ros_defense_genes_norm_counts

missing_requested_genes <- setdiff(ros_defense_genes_of_interest, rownames(norm_counts_d60))
if (length(missing_requested_genes) > 0) {
	stop(sprintf("Missing requested genes in day60 Neuron normalized counts: %s", paste(missing_requested_genes, collapse = ", ")))
}

ros_defense_neuron_norm_expr_long <- norm_counts_d60[ros_defense_genes_of_interest, , drop = FALSE] %>%
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
		gene = factor(gene, levels = ros_defense_genes_of_interest),
		genotype = factor(genotype_labels[genotype_d60], levels = c("WT", "HET", "KO"))
	)

missing_ros_defense_neuron_genes <- setdiff(ros_defense_genes_of_interest, unique(as.character(ros_defense_neuron_norm_expr_long$gene)))
if (length(missing_ros_defense_neuron_genes) > 0) {
	stop(sprintf("Missing requested genes in day60 Neuron normalized counts: %s", paste(missing_ros_defense_neuron_genes, collapse = ", ")))
}

ros_defense_neuron_norm_bar_df <- ros_defense_neuron_norm_expr_long %>%
	group_by(gene, genotype) %>%
	summarise(
		mean_expression = mean(expression),
		.groups = "drop"
	)

plot_day60_neuron_ros_defense_genes_norm_counts <- ggplot(ros_defense_neuron_norm_bar_df, aes(x = genotype, y = mean_expression, fill = genotype)) +
	geom_col(width = 0.72) +
	geom_point(
		data = ros_defense_neuron_norm_expr_long,
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
		title = "Day 60 Neuron ROS defense gene expression",
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

pdf(
	file = file.path(output_dir, "day60_grg_neuron_ros_defense_genes_normalized_counts_barplot.pdf"),
	width = 6,
	height = 7.6
)
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = 2, ncol = 1)))
print(plot_day60_grg_ros_defense_genes_norm_counts, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
print(plot_day60_neuron_ros_defense_genes_norm_counts, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
dev.off()

plot_day60_neuron_ros_defense_genes_norm_counts


###############################################
## mitophagy genes: PINK1, PRKN, BNIP3L
################################################


mitophagy_genes_of_interest <- c("PINK1", "PRKN", "BNIP3L")

missing_requested_genes <- setdiff(mitophagy_genes_of_interest, rownames(norm_counts_d60))
if (length(missing_requested_genes) > 0) {
	stop(sprintf("Missing requested genes in day60 gRG normalized counts: %s", paste(missing_requested_genes, collapse = ", ")))
}

mitophagy_norm_expr_long <- norm_counts_d60[mitophagy_genes_of_interest, , drop = FALSE] %>%
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
		gene = factor(gene, levels = mitophagy_genes_of_interest),
		genotype = factor(genotype_labels[genotype_d60], levels = c("WT", "HET", "KO"))
	)

missing_mitophagy_genes <- setdiff(mitophagy_genes_of_interest, unique(as.character(mitophagy_norm_expr_long$gene)))
if (length(missing_mitophagy_genes) > 0) {
	stop(sprintf("Missing requested genes in day60 gRG normalized counts: %s", paste(missing_mitophagy_genes, collapse = ", ")))
}

mitophagy_norm_bar_df <- mitophagy_norm_expr_long %>%
	group_by(gene, genotype) %>%
	summarise(
		mean_expression = mean(expression),
		.groups = "drop"
	)

plot_day60_grg_mitophagy_genes_norm_counts <- ggplot(mitophagy_norm_bar_df, aes(x = genotype, y = mean_expression, fill = genotype)) +
	geom_col(width = 0.72) +
	geom_point(
		data = mitophagy_norm_expr_long,
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
		title = "Day 60 gRG mitophagy gene expression",
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

plot_day60_grg_mitophagy_genes_norm_counts

missing_requested_genes <- setdiff(mitophagy_genes_of_interest, rownames(norm_counts_d60))
if (length(missing_requested_genes) > 0) {
	stop(sprintf("Missing requested genes in day60 Neuron normalized counts: %s", paste(missing_requested_genes, collapse = ", ")))
}

mitophagy_neuron_norm_expr_long <- norm_counts_d60[mitophagy_genes_of_interest, , drop = FALSE] %>%
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
		gene = factor(gene, levels = mitophagy_genes_of_interest),
		genotype = factor(genotype_labels[genotype_d60], levels = c("WT", "HET", "KO"))
	)

missing_mitophagy_neuron_genes <- setdiff(mitophagy_genes_of_interest, unique(as.character(mitophagy_neuron_norm_expr_long$gene)))
if (length(missing_mitophagy_neuron_genes) > 0) {
	stop(sprintf("Missing requested genes in day60 Neuron normalized counts: %s", paste(missing_mitophagy_neuron_genes, collapse = ", ")))
}

mitophagy_neuron_norm_bar_df <- mitophagy_neuron_norm_expr_long %>%
	group_by(gene, genotype) %>%
	summarise(
		mean_expression = mean(expression),
		.groups = "drop"
	)

plot_day60_neuron_mitophagy_genes_norm_counts <- ggplot(mitophagy_neuron_norm_bar_df, aes(x = genotype, y = mean_expression, fill = genotype)) +
	geom_col(width = 0.72) +
	geom_point(
		data = mitophagy_neuron_norm_expr_long,
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
		title = "Day 60 Neuron mitophagy gene expression",
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

pdf(
	file = file.path(output_dir, "day60_grg_neuron_mitophagy_genes_normalized_counts_barplot.pdf"),
	width = 9,
	height = 7.6
)
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = 2, ncol = 1)))
print(plot_day60_grg_mitophagy_genes_norm_counts, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
print(plot_day60_neuron_mitophagy_genes_norm_counts, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
dev.off()

plot_day60_neuron_mitophagy_genes_norm_counts

################################################################
### nuclear mitochondrial genes: TFAM, ATP5F1A, ATP5F1B, NDUFA1, SDHB, COX4I1
################################################################
nuclear_mito_genes_of_interest <- c("TFAM", "ATP5F1A", "ATP5F1B", "NDUFA1", "SDHB", "COX4I1")

missing_requested_genes <- setdiff(nuclear_mito_genes_of_interest, rownames(norm_counts_d60))
if (length(missing_requested_genes) > 0) {
	stop(sprintf("Missing requested genes in day60 gRG normalized counts: %s", paste(missing_requested_genes, collapse = ", ")))
}

nuclear_mito_grg_norm_expr_long <- norm_counts_d60[nuclear_mito_genes_of_interest, , drop = FALSE] %>%
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
		gene = factor(gene, levels = nuclear_mito_genes_of_interest),
		genotype = factor(genotype_labels[genotype_d60], levels = c("WT", "HET", "KO"))
	)

missing_nuclear_mito_grg_genes <- setdiff(nuclear_mito_genes_of_interest, unique(as.character(nuclear_mito_grg_norm_expr_long$gene)))
if (length(missing_nuclear_mito_grg_genes) > 0) {
	stop(sprintf("Missing requested genes in day60 gRG normalized counts: %s", paste(missing_nuclear_mito_grg_genes, collapse = ", ")))
}

nuclear_mito_grg_norm_bar_df <- nuclear_mito_grg_norm_expr_long %>%
	group_by(gene, genotype) %>%
	summarise(
		mean_expression = mean(expression),
		.groups = "drop"
	)

plot_day60_grg_nuclear_mito_genes_norm_counts <- ggplot(nuclear_mito_grg_norm_bar_df, aes(x = genotype, y = mean_expression, fill = genotype)) +
	geom_col(width = 0.72) +
	geom_point(
		data = nuclear_mito_grg_norm_expr_long,
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
		title = "Day 60 gRG nuclear mitochondrial gene expression",
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

missing_requested_genes <- setdiff(nuclear_mito_genes_of_interest, rownames(norm_counts_d60))
if (length(missing_requested_genes) > 0) {
	stop(sprintf("Missing requested genes in day60 Neuron normalized counts: %s", paste(missing_requested_genes, collapse = ", ")))
}

nuclear_mito_neuron_norm_expr_long <- norm_counts_d60[nuclear_mito_genes_of_interest, , drop = FALSE] %>%
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
		gene = factor(gene, levels = nuclear_mito_genes_of_interest),
		genotype = factor(genotype_labels[genotype_d60], levels = c("WT", "HET", "KO"))
	)

missing_nuclear_mito_neuron_genes <- setdiff(nuclear_mito_genes_of_interest, unique(as.character(nuclear_mito_neuron_norm_expr_long$gene)))
if (length(missing_nuclear_mito_neuron_genes) > 0) {
	stop(sprintf("Missing requested genes in day60 Neuron normalized counts: %s", paste(missing_nuclear_mito_neuron_genes, collapse = ", ")))
}

nuclear_mito_neuron_norm_bar_df <- nuclear_mito_neuron_norm_expr_long %>%
	group_by(gene, genotype) %>%
	summarise(
		mean_expression = mean(expression),
		.groups = "drop"
	)

plot_day60_neuron_nuclear_mito_genes_norm_counts <- ggplot(nuclear_mito_neuron_norm_bar_df, aes(x = genotype, y = mean_expression, fill = genotype)) +
	geom_col(width = 0.72) +
	geom_point(
		data = nuclear_mito_neuron_norm_expr_long,
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
		title = "Day 60 Neuron nuclear mitochondrial gene expression",
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

pdf(
	file = file.path(output_dir, "day60_grg_neuron_nuclear_mito_genes_normalized_counts_barplot.pdf"),
	width = 12,
	height = 7.6
)
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = 2, ncol = 1)))
print(plot_day60_grg_nuclear_mito_genes_norm_counts, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
print(plot_day60_neuron_nuclear_mito_genes_norm_counts, vp = grid::viewport(layout.pos.row = 2, layout.pos.col = 1))
dev.off()

plot_day60_neuron_nuclear_mito_genes_norm_counts


