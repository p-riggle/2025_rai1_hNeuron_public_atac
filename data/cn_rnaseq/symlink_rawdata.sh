#!/bin/bash

find /nfs/dataden/umms-siwase/SHIGEKI/rai1ko_hesc/data/fastq/usftp21.novogene.com/01.RawData -name '*.gz' -exec ln -s {} /nfs/dataden/umms-siwase/SHIGEKI/rai1ko_hesc/code/blueberry/input/fastq \;

