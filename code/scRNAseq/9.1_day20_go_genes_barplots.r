#########################################################################################
### Day 20 pseudobulk bar plot for integrin and cilia genes: NPC-G2/M and Mes-KO
#########################################################################################

library(tidyverse)
library(DESeq2)

cell_types      <- c("NPC-G2/M", "Mes-KO")
genotype_levels <- c("wt", "het", "ko")
genotype_labels <- c("wt" = "WT", "het" = "HET", "ko" = "KO")
genotype_colors <- c("WT" = "#7b7b7b", "HET" = "#f1d174", "KO" = "#bd6565")

# Load Day 20 data
dds_d20        <- readRDS("data/scRNAseq/iwaselab/dds_d20")
norm_counts_d20 <- counts(dds_d20, normalized = TRUE)
coldata_d20     <- as.data.frame(colData(dds_d20)) %>%
	rownames_to_column("sample_id")

output_dir <- "results/scRNAseq/7_pseudobulk/select_gene_barplots"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

#-------------------------------------------------------------------------------
# Integrin genes: NPC-G2/M and Mes-KO
#-------------------------------------------------------------------------------
integrin_genes <- c("ITGB2", "ITGA3", "ITGA9", "ITGA1", "ITGA5")

integrin_expr_long <- norm_counts_d20[integrin_genes, , drop = FALSE] %>%
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
		gene      = factor(gene, levels = integrin_genes),
		genotype  = factor(genotype_labels[genotype_d20], levels = c("WT", "HET", "KO")),
		cell_type = factor(cluster_d20, levels = cell_types)
	)

missing_integrin_genes <- setdiff(integrin_genes, unique(as.character(integrin_expr_long$gene)))
if (length(missing_integrin_genes) > 0) {
	stop(sprintf("Missing integrin genes in day20 normalized counts: %s", paste(missing_integrin_genes, collapse = ", ")))
}

integrin_bar_df <- integrin_expr_long %>%
	group_by(gene, cell_type, genotype) %>%
	summarise(mean_expression = mean(expression), .groups = "drop")

plot_day20_integrin_genes_norm_counts <- ggplot(integrin_bar_df, aes(x = cell_type, y = mean_expression, fill = genotype)) +
	geom_col(width = 0.75, position = position_dodge(width = 0.8)) +
	geom_point(
		data = integrin_expr_long,
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
		title = "Day 20 integrin genes pseudobulk expression (NPC-G2/M and Mes-KO)",
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
		axis.text.x = element_text(face = "bold"),
		plot.title = element_text(face = "bold")
	)

ggsave(
	filename = file.path(output_dir, "day20_integrin_genes_npcg2m_mesko_normalized_counts_barplot.pdf"),
	plot = plot_day20_integrin_genes_norm_counts,
	width = 13,
	height = 4
)

ggsave(
	filename = file.path(output_dir, "day20_integrin_genes_npcg2m_mesko_normalized_counts_barplot.png"),
	plot = plot_day20_integrin_genes_norm_counts,
	width = 13,
	height = 4,
	dpi = 300
)

print("Day 20 integrin genes bar plots created successfully!")


#-------------------------------------------------------------------------------
# Cilia genes: NPC-G2/M and Mes-KO
#-------------------------------------------------------------------------------
cilia_genes <- c("CFAP100", "CFAP47", "DNAI1", "DNAH6", "TEKT1")

cilia_expr_long <- norm_counts_d20[cilia_genes, , drop = FALSE] %>%
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
		gene      = factor(gene, levels = cilia_genes),
		genotype  = factor(genotype_labels[genotype_d20], levels = c("WT", "HET", "KO")),
		cell_type = factor(cluster_d20, levels = cell_types)
	)

missing_cilia_genes <- setdiff(cilia_genes, unique(as.character(cilia_expr_long$gene)))
if (length(missing_cilia_genes) > 0) {
	stop(sprintf("Missing cilia genes in day20 normalized counts: %s", paste(missing_cilia_genes, collapse = ", ")))
}

cilia_bar_df <- cilia_expr_long %>%
	group_by(gene, cell_type, genotype) %>%
	summarise(mean_expression = mean(expression), .groups = "drop")

plot_day20_cilia_genes_norm_counts <- ggplot(cilia_bar_df, aes(x = cell_type, y = mean_expression, fill = genotype)) +
	geom_col(width = 0.75, position = position_dodge(width = 0.8)) +
	geom_point(
		data = cilia_expr_long,
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
		title = "Day 20 cilia genes pseudobulk expression (NPC-G2/M and Mes-KO)",
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
		axis.text.x = element_text(face = "bold"),
		plot.title = element_text(face = "bold")
	)

ggsave(
	filename = file.path(output_dir, "day20_cilia_genes_npcg2m_mesko_normalized_counts_barplot.pdf"),
	plot = plot_day20_cilia_genes_norm_counts,
	width = 13,
	height = 4
)

ggsave(
	filename = file.path(output_dir, "day20_cilia_genes_npcg2m_mesko_normalized_counts_barplot.png"),
	plot = plot_day20_cilia_genes_norm_counts,
	width = 13,
	height = 4,
	dpi = 300
)

print("Day 20 cilia genes bar plots created successfully!")

#-------------------------------------------------------------------------------
# Postsynaptic genes: NPC-G2/M and Mes-KO
#-------------------------------------------------------------------------------
postsynaptic_genes <- c("ERBB4", "GRIK1", "LRRTM1", "GABRA3", "GSG1L")

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

missing_postsynaptic_genes <- setdiff(postsynaptic_genes, unique(as.character(postsynaptic_expr_long$gene)))
if (length(missing_postsynaptic_genes) > 0) {
	stop(sprintf("Missing postsynaptic genes in day20 normalized counts: %s", paste(missing_postsynaptic_genes, collapse = ", ")))
}

postsynaptic_bar_df <- postsynaptic_expr_long %>%
	group_by(gene, cell_type, genotype) %>%
	summarise(mean_expression = mean(expression), .groups = "drop")

plot_day20_postsynaptic_genes_norm_counts <- ggplot(postsynaptic_bar_df, aes(x = cell_type, y = mean_expression, fill = genotype)) +
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
		title = "Day 20 postsynaptic genes pseudobulk expression (NPC-G2/M and Mes-KO)",
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
		axis.text.x = element_text(face = "bold"),
		plot.title = element_text(face = "bold")
	)

ggsave(
	filename = file.path(output_dir, "day20_postsynaptic_genes_npcg2m_mesko_normalized_counts_barplot.pdf"),
	plot = plot_day20_postsynaptic_genes_norm_counts,
	width = 13,
	height = 4
)

ggsave(
	filename = file.path(output_dir, "day20_postsynaptic_genes_npcg2m_mesko_normalized_counts_barplot.png"),
	plot = plot_day20_postsynaptic_genes_norm_counts,
	width = 13,
	height = 4,
	dpi = 300
)

print("Day 20 postsynaptic genes bar plots created successfully!")

#-------------------------------------------------------------------------------
# Presynaptic membrane genes: NPC-G2/M and Mes-KO
#-------------------------------------------------------------------------------
presynaptic_genes <- c("ERBB4", "SLC1A6", "OPRK1", "KCNJ8", "ABCC8")

presynaptic_expr_long <- norm_counts_d20[presynaptic_genes, , drop = FALSE] %>%
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
		gene      = factor(gene, levels = presynaptic_genes),
		genotype  = factor(genotype_labels[genotype_d20], levels = c("WT", "HET", "KO")),
		cell_type = factor(cluster_d20, levels = cell_types)
	)

missing_presynaptic_genes <- setdiff(presynaptic_genes, unique(as.character(presynaptic_expr_long$gene)))
if (length(missing_presynaptic_genes) > 0) {
	stop(sprintf("Missing presynaptic genes in day20 normalized counts: %s", paste(missing_presynaptic_genes, collapse = ", ")))
}

presynaptic_bar_df <- presynaptic_expr_long %>%
	group_by(gene, cell_type, genotype) %>%
	summarise(mean_expression = mean(expression), .groups = "drop")

plot_day20_presynaptic_genes_norm_counts <- ggplot(presynaptic_bar_df, aes(x = cell_type, y = mean_expression, fill = genotype)) +
	geom_col(width = 0.75, position = position_dodge(width = 0.8)) +
	geom_point(
		data = presynaptic_expr_long,
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
		title = "Day 20 presynaptic membrane genes pseudobulk expression (NPC-G2/M and Mes-KO)",
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
		axis.text.x = element_text(face = "bold"),
		plot.title = element_text(face = "bold")
	)

ggsave(
	filename = file.path(output_dir, "day20_presynaptic_genes_npcg2m_mesko_normalized_counts_barplot.pdf"),
	plot = plot_day20_presynaptic_genes_norm_counts,
	width = 13,
	height = 4
)

ggsave(
	filename = file.path(output_dir, "day20_presynaptic_genes_npcg2m_mesko_normalized_counts_barplot.png"),
	plot = plot_day20_presynaptic_genes_norm_counts,
	width = 13,
	height = 4,
	dpi = 300
)

print("Day 20 presynaptic genes bar plots created successfully!")

#-------------------------------------------------------------------------------
# Vascular markers: NPC-G2/M and Mes-KO
#-------------------------------------------------------------------------------
vascular_genes <- c("PECAM1", "PLVAP", "MMRN2", "ABCG2", "CDH5")
vascular_expr_long <- norm_counts_d20[vascular_genes, , drop = FALSE] %>%
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
		gene      = factor(gene, levels = vascular_genes),
		genotype  = factor(genotype_labels[genotype_d20], levels = c("WT", "HET", "KO")),
		cell_type = factor(cluster_d20, levels = cell_types)
	)

missing_vascular_genes <- setdiff(vascular_genes, unique(as.character(vascular_expr_long$gene)))
if (length(missing_vascular_genes) > 0) {
	stop(sprintf("Missing vascular genes in day20 normalized counts: %s", paste(missing_vascular_genes, collapse = ", ")))
}		
vascular_bar_df <- vascular_expr_long %>%
	group_by(gene, cell_type, genotype) %>%
	summarise(mean_expression = mean(expression), .groups = "drop")
plot_day20_vascular_genes_norm_counts <- ggplot(vascular_bar_df, aes(x = cell_type, y = mean_expression, fill = genotype)) +
	geom_col(width = 0.75, position = position_dodge(width = 0.8)) +
	geom_point(
		data = vascular_expr_long,
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
		title = "Day 20 vascular genes pseudobulk expression (NPC-G2/M and Mes-KO)",
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
		axis.text.x = element_text(face = "bold"),
		plot.title = element_text(face = "bold")
	)
ggsave(
	filename = file.path(output_dir, "day20_vascular_genes_npcg2m_mesko_normalized_counts_barplot.pdf"),
	plot = plot_day20_vascular_genes_norm_counts,
	width = 13,
	height = 4
)	

print("Day 20 vascular genes bar plots created successfully!")

#-------------------------------------------------------------------------------
# immune signaling: SYK, TRBC1, TRBC2

immune_genes <- c("SYK", "TRBC1", "TRBC2")
immune_expr_long <- norm_counts_d20[immune_genes, , drop = FALSE] %>%
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
		gene      = factor(gene, levels = immune_genes),
		genotype  = factor(genotype_labels[genotype_d20], levels = c("WT", "HET", "KO")),
		cell_type = factor(cluster_d20, levels = cell_types)
	)
missing_immune_genes <- setdiff(immune_genes, unique(as.character(immune_expr_long$gene)))
if (length(missing_immune_genes) > 0) {
	stop(sprintf("Missing immune genes in day20 normalized counts: %s", paste(missing_immune_genes, collapse = ", ")))
}
immune_bar_df <- immune_expr_long %>%
	group_by(gene, cell_type, genotype) %>%
	summarise(mean_expression = mean(expression), .groups = "drop")
plot_day20_immune_genes_norm_counts <- ggplot(immune_bar_df, aes(x = cell_type, y = mean_expression, fill = genotype)) +
	geom_col(width = 0.75, position = position_dodge(width = 0.8)) +
	geom_point(
		data = immune_expr_long,
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
		title = "Day 20 immune genes pseudobulk expression (NPC-G2/M and Mes-KO)",
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
		axis.text.x = element_text(face = "bold"),
		plot.title = element_text(face = "bold")
	)
ggsave(
	filename = file.path(output_dir, "day20_immune_genes_npcg2m_mesko_normalized_counts_barplot.pdf"),
	plot = plot_day20_immune_genes_norm_counts,
	width = 8,
	height = 4
)	
print("Day 20 immune genes bar plots created successfully!")


#-------------------------------------------------------------------------------
# myogenic genes: ACTA1, TTN (Titin), ACTC1

myogenic_genes <- c("ACTA1", "TTN", "ACTC1")
myogenic_cell_types <- c("NPC-G2/M", "NPC-G1", "Mes-WT", "IPC")

myogenic_expr_long <- norm_counts_d20[myogenic_genes, , drop = FALSE] %>%
	as.data.frame() %>%
	rownames_to_column("gene") %>%
	pivot_longer(
		cols = -gene,
		names_to = "sample_id",
		values_to = "expression"
	) %>%
	left_join(coldata_d20, by = "sample_id") %>%
	filter(cluster_d20 %in% myogenic_cell_types, genotype_d20 %in% genotype_levels) %>%
	mutate(
		gene = factor(gene, levels = myogenic_genes),
		genotype = factor(genotype_labels[genotype_d20], levels = c("WT", "HET", "KO")),
		cell_type = factor(cluster_d20, levels = myogenic_cell_types)
	)

missing_myogenic_genes <- setdiff(myogenic_genes, unique(as.character(myogenic_expr_long$gene)))
if (length(missing_myogenic_genes) > 0) {
	stop(sprintf("Missing myogenic genes in day20 normalized counts: %s", paste(missing_myogenic_genes, collapse = ", ")))
}

myogenic_bar_df <- myogenic_expr_long %>%
	group_by(gene, cell_type, genotype) %>%
	summarise(mean_expression = mean(expression), .groups = "drop")

plot_day20_myogenic_genes_norm_counts <- ggplot(myogenic_bar_df, aes(x = cell_type, y = mean_expression, fill = genotype)) +
	geom_col(width = 0.75, position = position_dodge(width = 0.8)) +
	geom_point(
		data = myogenic_expr_long,
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
		title = "Day 20 myogenic genes pseudobulk expression (NPC-G2/M, NPC-G1, Mes-WT, IPC)",
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
	filename = file.path(output_dir, "day20_myogenic_genes_npcg2m_npcg1_meswt_ipc_normalized_counts_barplot.pdf"),
	plot = plot_day20_myogenic_genes_norm_counts,
	width = 11,
	height = 4
)

ggsave(
	filename = file.path(output_dir, "day20_myogenic_genes_npcg2m_npcg1_meswt_ipc_normalized_counts_barplot.png"),
	plot = plot_day20_myogenic_genes_norm_counts,
	width = 11,
	height = 4,
	dpi = 300
)

print("Day 20 myogenic genes bar plots created successfully!")


#-------------------------------------------------------------------------------
# het-immune genes:  CD8A, BST2

het_immune_genes <- c("CD8A", "BST2")
het_immune_cell_types <- c("NPC-G2/M", "NPC-G1", "Mes-WT", "IPC")

het_immune_expr_long <- norm_counts_d20[het_immune_genes, , drop = FALSE] %>%
	as.data.frame() %>%
	rownames_to_column("gene") %>%
	pivot_longer(
		cols = -gene,
		names_to = "sample_id",
		values_to = "expression"
	) %>%
	left_join(coldata_d20, by = "sample_id") %>%
	filter(cluster_d20 %in% het_immune_cell_types, genotype_d20 %in% genotype_levels) %>%
	mutate(
		gene = factor(gene, levels = het_immune_genes),
		genotype = factor(genotype_labels[genotype_d20], levels = c("WT", "HET", "KO")),
		cell_type = factor(cluster_d20, levels = het_immune_cell_types)
	)

missing_het_immune_genes <- setdiff(het_immune_genes, unique(as.character(het_immune_expr_long$gene)))
if (length(missing_het_immune_genes) > 0) {
	stop(sprintf("Missing het-immune genes in day20 normalized counts: %s", paste(missing_het_immune_genes, collapse = ", ")))
}

het_immune_bar_df <- het_immune_expr_long %>%
	group_by(gene, cell_type, genotype) %>%
	summarise(mean_expression = mean(expression), .groups = "drop")

plot_day20_het_immune_genes_norm_counts <- ggplot(het_immune_bar_df, aes(x = cell_type, y = mean_expression, fill = genotype)) +
	geom_col(width = 0.75, position = position_dodge(width = 0.8)) +
	geom_point(
		data = het_immune_expr_long,
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
		title = "Day 20 het-immune genes pseudobulk expression (NPC-G2/M, NPC-G1, Mes-WT, IPC)",
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
	filename = file.path(output_dir, "day20_het_immune_genes_npcg2m_npcg1_meswt_ipc_normalized_counts_barplot.pdf"),
	plot = plot_day20_het_immune_genes_norm_counts,
	width = 9,
	height = 4
)

ggsave(
	filename = file.path(output_dir, "day20_het_immune_genes_npcg2m_npcg1_meswt_ipc_normalized_counts_barplot.png"),
	plot = plot_day20_het_immune_genes_norm_counts,
	width = 9,
	height = 4,
	dpi = 300
)

print("Day 20 het-immune genes bar plots created successfully!")

#-------------------------------------------------------------------------------
# non-neuronal epithelial genes: CLDN1, CLDN4, CFTR, CLIC5

epithelial_genes <- c("CLDN1", "CLDN4", "CFTR", "CLIC5")
epithelial_cell_types <- c("NPC-G2/M", "NPC-G1", "Mes-WT", "IPC")

epithelial_expr_long <- norm_counts_d20[epithelial_genes, , drop = FALSE] %>%
	as.data.frame() %>%
	rownames_to_column("gene") %>%
	pivot_longer(
		cols = -gene,
		names_to = "sample_id",
		values_to = "expression"
	) %>%
	left_join(coldata_d20, by = "sample_id") %>%
	filter(cluster_d20 %in% epithelial_cell_types, genotype_d20 %in% genotype_levels) %>%
	mutate(
		gene = factor(gene, levels = epithelial_genes),
		genotype = factor(genotype_labels[genotype_d20], levels = c("WT", "HET", "KO")),
		cell_type = factor(cluster_d20, levels = epithelial_cell_types)
	)

missing_epithelial_genes <- setdiff(epithelial_genes, unique(as.character(epithelial_expr_long$gene)))
if (length(missing_epithelial_genes) > 0) {
	stop(sprintf("Missing non-neuronal epithelial genes in day20 normalized counts: %s", paste(missing_epithelial_genes, collapse = ", ")))
}

epithelial_bar_df <- epithelial_expr_long %>%
	group_by(gene, cell_type, genotype) %>%
	summarise(mean_expression = mean(expression), .groups = "drop")

plot_day20_epithelial_genes_norm_counts <- ggplot(epithelial_bar_df, aes(x = cell_type, y = mean_expression, fill = genotype)) +
	geom_col(width = 0.75, position = position_dodge(width = 0.8)) +
	geom_point(
		data = epithelial_expr_long,
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
		title = "Day 20 non-neuronal epithelial genes pseudobulk expression (NPC-G2/M, NPC-G1, Mes-WT, IPC)",
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
	filename = file.path(output_dir, "day20_non_neuronal_epithelial_genes_npcg2m_npcg1_meswt_ipc_normalized_counts_barplot.pdf"),
	plot = plot_day20_epithelial_genes_norm_counts,
	width = 12,
	height = 4
)

ggsave(
	filename = file.path(output_dir, "day20_non_neuronal_epithelial_genes_npcg2m_npcg1_meswt_ipc_normalized_counts_barplot.png"),
	plot = plot_day20_epithelial_genes_norm_counts,
	width = 12,
	height = 4,
	dpi = 300
)

print("Day 20 non-neuronal epithelial genes bar plots created successfully!")



