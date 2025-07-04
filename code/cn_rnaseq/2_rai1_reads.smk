# specify config file

configfile: "code/2025-02-12_config.yaml"

# Define a global variable for the sleep duration (in seconds)
sleep_duration = 10  # Adjust as needed

rule all:
    input:
        expand("data/bam_rai1/{sample}.sorted.rai1.bam", sample=config["samples"]),
        expand("data/bam_rai1/{sample}.sorted.rai1.bam.bai", sample=config["samples"]), 
        "data/bam_rai1/1_upload_urls.txt"
        
rule extract_reads_rai1: 
    input: 
        fq1=lambda wildcards: config["samples"][wildcards.sample]["bam"]
    output:
        bam_rai1 = "data/bam_rai1/{sample}.rai1.bam"
    shell: 
        """
        module load Bioinformatics
        module load samtools
        samtools view -b {input.fq1} "chr17:17595154-17840025" > {output.bam_rai1}
        sleep {sleep_duration}  # Introduce a delay to avoid overloading the server
        """ 

rule sort_bam:
    input:
        bam_rai1 = "data/bam_rai1/{sample}.rai1.bam"
    output:
        sorted_bam = "data/bam_rai1/{sample}.sorted.rai1.bam"
    shell:
        """
        module load Bioinformatics
        module load samtools
        samtools sort -o {output.sorted_bam} {input.bam_rai1}
        sleep {sleep_duration}  # Introduce a delay to avoid overloading the server
        """

rule index_bam:
    input:
        sorted_bam = "data/bam_rai1/{sample}.sorted.rai1.bam"  
    output:
        index_bam = "data/bam_rai1/{sample}.sorted.rai1.bam.bai"
    shell:
        """
        module load Bioinformatics
        module load samtools
        samtools index {input.sorted_bam} {output.index_bam}
        sleep {sleep_duration}  # Introduce a delay to avoid overloading the server
        """
        
rule upload_bam_rai1:
    input:
        config = "code/2025-02-12_config.yaml"
    output:
        txt = "data/bam_rai1/1_upload_urls.txt"
    shell: 
        """
        python3 code/browser_upload2.py
        """

