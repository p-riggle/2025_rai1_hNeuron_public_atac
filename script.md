2026.03.12

STAR --runThreadN 16 \
             --genomeDir /nfs/turbo/umms-siwase/GENOME/MM39 \
             --readFilesIn data/cn_rnaseq/fastq/usftp21.novogene.com/01.RawData/HYc1t2/HYc1t2_1.fq.gz \
             --readFilesCommand zcat \
             --outSAMtype BAM SortedByCoordinate \
             --outFileNamePrefix data/bam/contrl_shRNA_rep1
             --outSAMunmapped Within \
             --outSAMattributes Standard

