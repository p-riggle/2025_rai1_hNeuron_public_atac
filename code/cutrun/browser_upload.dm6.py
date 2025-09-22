import yaml

# Load the config file
with open('code/config.yaml', 'r') as file:
    config = yaml.safe_load(file)

# Extract sample names from the config
samples = list(config['samples'].keys())
        
base_url_template = 'track type=bigWig name="{sample}" bigDataUrl=https://storage.googleapis.com/siwase/an_rai1_cutrun_pilot/4_rai1_cr_ko_spikein/dm6/dm6.{sample}.bw visibility=full color=235,158,52 autoScale=on alwaysZero=off gridDefault=off maxHeightPixels=60 graphType=bar viewLimits=0:5 windowingFunction=mean smoothingWindow=2\n'


# Create the list of URLs
urls = [base_url_template.format(sample=sample) for sample in samples]

# Write the URLs to a text file
with open('data/bigwig_dm6/1_upload.dm6.txt', 'w') as file:
    for url in urls:
        file.write(url + '\n')