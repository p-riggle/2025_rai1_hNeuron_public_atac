#!/bin/bash

# Prefetch
cd /data/fastq/kaster_mesc
prefetch $(<SRR_Acc_List.txt)

nohup fasterq-dump $(<SRR_Acc_List.txt) &

gzip *.fastq