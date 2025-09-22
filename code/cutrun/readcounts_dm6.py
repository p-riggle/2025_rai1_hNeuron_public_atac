import os
import csv
import subprocess

# Function to count uniquely mapped reads in a BAM file using samtools
def count_uniquely_mapped_reads_dm6(bam_file):
    # Run samtools command to count uniquely mapped reads
    result = subprocess.run(['samtools', 'view', '-F', '0x4', '-F', '0x100', '-c', bam_file], stdout=subprocess.PIPE, check=True)
    return int(result.stdout.strip())

# Directory containing BAM files
bam_directory = "data/bam_dm6/"

# List to store the results
results = []

# Iterate over each BAM file in the directory
for filename in os.listdir(bam_directory):
    if filename.startswith("sort_") and filename.endswith(".bam"):
        bam_file_path = os.path.join(bam_directory, filename)
        read_count = count_uniquely_mapped_reads_dm6(bam_file_path)
        results.append([filename, read_count])

# Write the results to a CSV file
with open("bam_read_counts_dm6.csv", "w", newline="") as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(["BAM File", "Uniquely Mapped Reads"])
    csvwriter.writerows(results)

print("The CSV file 'bam_read_counts_dm6.csv' has been created with the read counts.")
