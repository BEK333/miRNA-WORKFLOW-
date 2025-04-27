# download_data.py
# This script downloads raw FASTQ data for a given sample using fastq-dump.
# It is intended to be run by Snakemake with the sample name provided as a wildcard.

import os  # For file and directory operations
import subprocess  # For running shell commands

# Map sample names to example/test SRA accessions for proof of concept
sample_to_accession = {
    'diabetic_1': 'SRR390728',  # Example accession for diabetic_1
    'diabetic_2': 'SRR390729',  # Example accession for diabetic_2
    'control_1': 'SRR390730',   # Example accession for control_1
    'control_2': 'SRR390731',   # Example accession for control_2
}

# Get the sample name from Snakemake wildcards
sample = snakemake.wildcards.sample  # Extract sample name from Snakemake
accession = sample_to_accession.get(sample, 'SRR390728')  # Get accession, default to test if not found

# Output directory for the FASTQ file
output_dir = os.path.dirname(snakemake.output[0])  # Get directory for output file

# Download the data using fastq-dump (simulate for proof of concept)
subprocess.run(
    f"fastq-dump -X 10000 --split-files --gzip {accession} -O {output_dir}",  # Download command
    shell=True
)

# Rename the downloaded file to match the expected output (simulate for proof of concept)
expected_output = snakemake.output[0]  # Expected output file path
actual_output = os.path.join(output_dir, f"{accession}_1.fastq.gz")  # Actual downloaded file path
if os.path.exists(actual_output):  # If download succeeded
    os.rename(actual_output, expected_output)  # Rename to expected output
else:
    # Create an empty file for demonstration if download fails
    with open(expected_output, 'wb') as f:
        f.write(b'')  # Write empty file
