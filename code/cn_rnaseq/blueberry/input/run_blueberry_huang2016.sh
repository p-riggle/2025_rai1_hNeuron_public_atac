#!/bin/bash

conda activate blueberry

snakemake \
    --snakefile ~/Apricot/Snakefile_nobarcode \
    --configfile /nfs/dataden/umms-siwase/SHIGEKI/rai1ko_hesc/code/blueberry/input/2023-06-16_config_huang2016.yaml \
    --dryrun --quiet

nohup snakemake \
    --snakefile ~/Apricot/Snakefile_nobarcode \
    --configfile /nfs/dataden/umms-siwase/SHIGEKI/rai1ko_hesc/code/blueberry/input/2023-06-16_config_huang2016.yaml \
    --resources mem_gb=105 \
    --cores 36 \
    -p &

