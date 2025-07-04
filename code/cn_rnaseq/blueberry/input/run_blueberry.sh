#!/bin/bash

snakemake \
    --snakefile ~/Apricot/Snakefile_nobarcode \
    --configfile /nfs/dataden/umms-siwase/SHIGEKI/rai1ko_hesc/code/blueberry/input/2022-01-17_config.yaml \
    --dryrun --quiet

nohup snakemake \
    --snakefile ~/Apricot/Snakefile_nobarcode \
    --configfile /nfs/dataden/umms-siwase/SHIGEKI/rai1ko_hesc/code/blueberry/input/2022-01-17_config.yaml \
    --resources mem_gb=105 \
    --cores 36 \
    -p &

# merge replicates of resulted bigwigs
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

conda install -c bioconda ucsc-bigwigmerge
# or
conda install -c "bioconda/label/cf201901" ucsc-bigwigmerge






