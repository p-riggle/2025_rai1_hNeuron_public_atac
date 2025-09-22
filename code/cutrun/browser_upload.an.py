import yaml

# Load the config file
with open('code/config.yaml', 'r') as file:
    config = yaml.safe_load(file)

# Extract sample names from the config
samples = list(config['samples'].keys())
        
base_url_template = 'track type=bigWig name="{sample}" bigDataUrl=https://storage.googleapis.com/siwase/an_h3q5dop_cutrun_pilot/{sample}.bw visibility=full color=255,222,61 autoScale=on alwaysZero=off gridDefault=off maxHeightPixels=60 graphType=bar viewLimits=0:5 windowingFunction=mean smoothingWindow=2\n'


# Create the list of URLs
urls = [base_url_template.format(sample=sample) for sample in samples]

# Write the URLs to a text file
with open('data/bigwig/1_upload.txt', 'w') as file:
    for url in urls:
        file.write(url + '\n')