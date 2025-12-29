# RAI1 hESC Multi-omics Analysis Project

## Project Overview
This is a neurodevelopmental genomics study examining RAI1 gene function in human embryonic stem cells (hESCs) using multiple genomic approaches. The project analyzes the developmental trajectory from hESCs to neurons across three genotypes (WT, HET, KO) at multiple timepoints (day 0, 20, 34, 60).

## Core Analysis Workflows

### Environment Setup
- **Container**: Use Singularity container `single_cell_env_0.12.0.sif` for scRNA-seq analysis
- **Conda environments**: `snakemake` for pipeline execution, `blueberry` for bulk RNA-seq
- **SLURM clusters**: Jobs run on HPC with module loading (`module load Bioinformatics`, `module load star`)

### Data Architecture
- **Raw data**: `data/{cutrun,cn_rnaseq,scRNAseq,ineuron}/` - FASTQ files and 10X outputs
- **Processed data**: BAM files, bigwig tracks, count matrices
- **Results**: `results/{analysis_type}/` with structured subdirectories by analysis step
- **Code**: `code/{analysis_type}/` containing R scripts, Snakemake workflows, config files

### Key Analysis Pipelines

#### 1. scRNA-seq Workflow (Primary analysis)
**Sequential R scripts pattern**: `1_first_steps.r` → `2_clustering.r` → `3_celltype.r` → `10_monocle.r`
- Uses **Seurat + BPCells + Monocle3** stack for memory efficiency
- **Critical pattern**: Always convert sparse matrices with `as(counts, "sparseMatrix")` for Monocle3 compatibility
- **Trajectory analysis**: Manual root cell selection from WT NPC-G2/M cells in specific UMAP regions
- **Cell types**: NPC-G1/G2M, IPC, Neuron, Astrocyte, Mes-WT/KO with consistent color scheme

#### 2. Bulk RNA-seq (CN-RNA-seq)
**Snakemake + DESeq2 pattern**: `2_rai1_reads.smk` → `deseq2/s1_deseq2.r` → comparative analysis scripts
- **Config-driven**: YAML files specify sample BAM/BAI file pairs from RSEM-STAR output
- **DESeq2 conventions**: Result files named `rstbl_{timepoint}_{comparison}.txt`
- **Blueberry integration**: External pipeline at `~/Apricot/Snakefile_nobarcode` for preprocessing

#### 3. CUT&RUN Analysis
**Pattern**: `bam_bw.smk` → `browser_upload.{an,dm6}.py` for UCSC browser tracks
- **Dual-genome**: Human (hg38) and spike-in Drosophila (dm6) analysis
- **Sleep delays**: 10-second delays between jobs to prevent server overload

### Development Patterns

#### Configuration Management
- **Sample definitions**: YAML configs map sample names to file paths
- **Naming convention**: `{genotype}_{timepoint}_rep{n}` (e.g., `wt_day0_rep1`)
- **File paths**: Always use absolute paths from project root

#### R Analysis Conventions
- **Working directory**: Always `setwd("/nfs/turbo/umms-siwase/SHIGEKI2/rai1_hesc")`
- **Color schemes**: Consistent cell type colors defined in `my_cell_colors` vector
- **Plot saving**: Use loops for batch plot generation, consistent sizing (width=12, height=8)
- **Object persistence**: Monocle3 objects saved with `save_monocle_objects()`

#### HPC Job Submission
- **SBATCH headers**: Templates in `d1.sh` files with 16-24 cores typical
- **Module loading**: Bioinformatics modules required for most genomics tools
- **Conda activation**: Source conda profile before environment activation

## Critical Integration Points
- **Seurat→Monocle3**: Transfer UMAP coordinates and manually set clusters/partitions
- **BigWig generation**: Links to Google Cloud Storage for browser visualization
- **Cross-timepoint analysis**: DESeq2 results integrated across developmental stages
- **Genotype comparisons**: Systematic WT vs HET vs KO analysis at each timepoint

## Common Issues & Solutions
- **Memory management**: Use BPCells for large single-cell datasets, garbage collect regularly
- **Sparse matrix compatibility**: Convert matrices for Monocle3 with explicit casting
- **Pseudotime root selection**: Use regional filtering of WT NPC-G2/M cells, not all cells
- **File path consistency**: Always use project-relative paths in configs and scripts

When contributing to analysis scripts, follow the established numerical naming scheme and save intermediate objects for reproducibility across the complex multi-step workflows.