# specify config file

configfile: "code/cutrun/config.yaml"

# Define a global variable for the sleep duration (in seconds)
sleep_duration = 10  # Adjust as needed

rule all:
    input:
        expand("data/cutrun/bam/sort_{sample}.bam", sample=config["samples"]), 
        expand("data/cutrun/bam/sort_{sample}.bam.bai", sample=config["samples"]),
        expand("data/cutrun/bam/bt2_{sample}.log", sample = config["samples"]),
        expand("data/cutrun/bigwig/{sample}.bw", sample=config["samples"]),
        expand("data/cutrun/bam_dm6/sort_dm6.{sample}.bam", sample=config["samples"]),
        expand("data/cutrun/bam_dm6/sort_dm6.{sample}.bam.bai", sample=config["samples"]),
        expand("data/cutrun/bam_dm6/dm6.{sample}.log", sample = config["samples"]),
        expand("data/cutrun/bigwig_dm6/dm6.{sample}.bw", sample=config["samples"]), 
        "data/cutrun/bigwig/1_upload.txt",
        "data/cutrun/bigwig_dm6/1_upload.dm6.txt"
        # "bam_read_count.csv",
        # "bam_read_count_dm6.csv"

# read mapping and converting sam to bam
rule bowtie2:
    input:
        fq1=lambda wildcards: config["samples"][wildcards.sample]["fq1"],
        fq2=lambda wildcards: config["samples"][wildcards.sample]["fq2"]
    output:
        bam1 = "data/cutrun/bam/ori_{sample}.bam"
    params:
        index = "/nfs/turbo/umms-siwase/GENOME_INDICES/HG38/bt2_hg38-index" # change
    log:
        log = "data/cutrun/bam/bt2_{sample}.log"
    resources:
        mem_mb = 120000 # 120 GB
    shell:
        """
        module load Bioinformatics  gcc/10.3.0-k2osx5y
        module load bowtie2
        module load samtools
        bowtie2 -p 16 -q --local \
        -x {params.index} \
        -1 {input.fq1} \
        -2 {input.fq2} | samtools view -h -S -b -q 10 --threads 16 > {output.bam1} 
        samtools stats {output.bam1} > data/cutrun/bam/{wildcards.sample}_raw_stats.txt
        sleep {sleep_duration}  # Introduce a delay to avoid overloading the server
        """

# Define the rule for sorting BAM
rule bam_sort:
    input:
        bam2 = "data/cutrun/bam/ori_{sample}.bam"
    output:
        bam_sorted = "data/cutrun/bam/sort_{sample}.bam"
    resources:
        mem_mb = 240000 # 240 GB
    shell:
        """
        module load Bioinformatics  gcc/10.3.0-k2osx5y
        module load samtools
        samtools sort -o {output.bam_sorted} -T {wildcards.sample}_temp --threads 16 {input.bam2} -m 1G
        sleep {sleep_duration}  # Introduce a delay to avoid overloading the server
        """

# Define the rule for indexing BAM
rule bam_index:
    input:
        bam_sorted = "data/cutrun/bam/sort_{sample}.bam"
    output:
        bai1 = "data/cutrun/bam/sort_{sample}.bam.bai"
    shell:
        """
        module load Bioinformatics  gcc/10.3.0-k2osx5y
        module load samtools
        samtools index {input.bam_sorted}
        sleep {sleep_duration}  # Introduce a delay to avoid overloading the server
        """

rule bam_to_bigwig:
    input: 
        bam3 = "data/cutrun/bam/sort_{sample}.bam",
        bai2 = "data/cutrun/bam/sort_{sample}.bam.bai"
    output:
        bigwig = "data/cutrun/bigwig/{sample}.bw"
    shell: 
        """
        source ~/anaconda3/etc/profile.d/conda.sh  # change if needed
        conda activate deeptools
        bamCoverage -b {input.bam3} -of bigwig -o {output.bigwig} --normalizeUsing BPM -p 16
        sleep {sleep_duration}  # Introduce a delay to avoid overloading the server
        """
    
rule bowtie2_dm6:
    input:
        fq1=lambda wildcards: config["samples"][wildcards.sample]["fq1"],
        fq2=lambda wildcards: config["samples"][wildcards.sample]["fq2"]
    output:
        bam1 = "data/cutrun/bam_dm6/dm6.{sample}.bam"
    params:
        index = "/nfs/turbo/umms-siwase/GENOME_INDICES/DM6/dm" # change
    log:
        log = "data/cutrun/bam_dm6/dm6.{sample}.log"
    resources:
        mem_mb = 240000 # 240 GB
    shell:
        """
        module load Bioinformatics  gcc/10.3.0-k2osx5y
        module load bowtie2
        module load samtools
        bowtie2 -p 16 -q --local \
        -x {params.index} \
        -1 {input.fq1} \
        -2 {input.fq2} | samtools view -h -S -b -q 10 --threads 16 > {output.bam1} 
        samtools stats {output.bam1} > data/cutrun/bam_dm6/dm6.{wildcards.sample}_raw_stats.txt
        sleep {sleep_duration}  # Introduce a delay to avoid overloading the server
        """
rule bam_sort_dm6:
    input:
        bam2 = "data/cutrun/bam_dm6/dm6.{sample}.bam"
    output:
        bam_sorted = "data/cutrun/bam_dm6/sort_dm6.{sample}.bam"
    resources:
        mem_mb = 240000 # 240 GB
    shell:
        """
        module load Bioinformatics  gcc/10.3.0-k2osx5y
        module load samtools
        samtools sort -o {output.bam_sorted} -T dm6.{wildcards.sample}_temp --threads 16 {input.bam2} -m 1G
        sleep {sleep_duration}  # Introduce a delay to avoid overloading the server
        """
rule bam_index_dm6:
    input:
        bam_sorted = "data/cutrun/bam_dm6/sort_dm6.{sample}.bam"
    output:
        bai1 = "data/cutrun/bam_dm6/sort_dm6.{sample}.bam.bai"
    shell:
        """
        module load Bioinformatics  gcc/10.3.0-k2osx5y
        module load samtools
        samtools index {input.bam_sorted}
        sleep {sleep_duration}  # Introduce a delay to avoid overloading the server
        """
rule bam_to_bigwig_dm6:
    input: 
        bam3 = "data/cutrun/bam_dm6/sort_dm6.{sample}.bam",
        bai2 = "data/cutrun/bam_dm6/sort_dm6.{sample}.bam.bai"
    output:
        bigwig = "data/cutrun/bigwig_dm6/dm6.{sample}.bw"
    shell: 
        """
        source ~/anaconda3/etc/profile.d/conda.sh  # change if needed
        conda activate deeptools
        bamCoverage -b {input.bam3} -of bigwig -o {output.bigwig} --normalizeUsing BPM -p 16
        sleep {sleep_duration}  # Introduce a delay to avoid overloading the server
        """

rule uploading_txt: 
    output: 
        "data/cutrun/bigwig/1_upload.txt",
        "data/cutrun/bigwig_dm6/1_upload.dm6.txt"
    shell: 
        """
        python3 code/browser_upload.an.py
        python3 code/browser_upload.dm6.py
        """
# rule count_reads: 
#     output: 
#         "bam_read_count.csv",
#         "bam_read_count_dm6.csv"
#     shell: 
#         """
#         module load samtools
#         python3 code/readcounts.py
#         python3 code/readcounts_dm6.py
#         """
