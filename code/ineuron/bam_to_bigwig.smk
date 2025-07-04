# specify config file

configfile: "code/config.yaml"

# Define a global variable for the sleep duration (in seconds)
sleep_duration = 10  # Adjust as needed

rule all:
    input:
        expand("data/bigwig/{sample}.fwd.bw", sample=config["samples"]),
        expand("data/bigwig/{sample}.rev.bw", sample=config["samples"]),
        expand("data/bam_rai1/{sample}.bam", sample=config["samples"]),
        expand("data/bam_rai1/{sample}.sorted.bam", sample=config["samples"]),
        expand("data/bam_rai1/{sample}.sorted.bam.bai", sample=config["samples"]), 
        "data/bam_rai1/1_upload_urls.txt"
        

rule bam_to_bigwig:
    input: 
        fq1=lambda wildcards: config["samples"][wildcards.sample]["bam"],
        fq2=lambda wildcards: config["samples"][wildcards.sample]["bai"]
    output:
        bigwig_f = "data/bigwig/{sample}.fwd.bw",
        bigwig_r = "data/bigwig/{sample}.rev.bw"
    shell: 
        """
        source ~/anaconda3/etc/profile.d/conda.sh
        conda activate deeptools
        bamCoverage -b {input.fq1} -of bigwig --filterRNAstrand forward -o {output.bigwig_f} --normalizeUsing BPM -p 16
        bamCoverage -b {input.fq1} -of bigwig --filterRNAstrand reverse -o {output.bigwig_r} --normalizeUsing BPM -p 16
        sleep {sleep_duration}  # Introduce a delay to avoid overloading the server
        """

rule extract_reads_rai1: 
    input: 
        fq1=lambda wildcards: config["samples"][wildcards.sample]["bam"],
        fq2=lambda wildcards: config["samples"][wildcards.sample]["bai"]
    output:
        bam_rai1 = "data/bam_rai1/{sample}.bam"
    shell: 
        """
        module load Bioinformatics
        module load samtools
        samtools view -b {input.fq1} 17:17595154-17840025 > {output.bam_rai1}
        sleep {sleep_duration}  # Introduce a delay to avoid overloading the server
        """ 

rule sort_bam:
    input:
        bam_rai1 = "data/bam_rai1/{sample}.bam"
    output:
        sorted_bam = "data/bam_rai1/{sample}.sorted.bam"
    shell:
        """
        module load Bioinformatics
        module load samtools
        samtools sort -o {output.sorted_bam} {input.bam_rai1}
        sleep {sleep_duration}  # Introduce a delay to avoid overloading the server
        """

rule index_bam:
    input:
        sorted_bam = "data/bam_rai1/{sample}.sorted.bam"  
    output:
        index_bam = "data/bam_rai1/{sample}.sorted.bam.bai"
    shell:
        """
        module load Bioinformatics
        module load samtools
        samtools index {input.sorted_bam} {output.index_bam}
        sleep {sleep_duration}  # Introduce a delay to avoid overloading the server
        """
        
rule upload_bam_rai1:
    input:
        config = "code/config.yaml"
    output:
        txt = "data/bam_rai1/1_upload_urls.txt"
    shell: 
        """
        python3 data/bam_rai1/browser_upload2.py
        """

