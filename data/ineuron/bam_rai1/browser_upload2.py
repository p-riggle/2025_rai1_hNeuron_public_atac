import yaml

# Load the config file
with open('code/config.yaml', 'r') as file:
    config = yaml.safe_load(file)

# Extract sample names from the config
samples = list(config['samples'].keys())

# Generating a textfile for uploading bam files to the UCSC browser
base_url = "https://storage.googleapis.com/siwase/rai1ko_hesc/ineuron/bam/"

# Create the list of URLs
urls = [f"{base_url}{sample}.sorted.bam" for sample in samples]

# Write the URLs to a text file
with open('data/bam_rai1/1_upload_urls.txt', 'w') as file:
    for url in urls:
        file.write(url + '\n')