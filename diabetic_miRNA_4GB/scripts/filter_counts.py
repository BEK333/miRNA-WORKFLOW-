# filter_counts.py
# This script filters out miRNAs with low counts across all samples.
# It is intended to be run by Snakemake after counting and before DE analysis.

import pandas as pd  # Import pandas for data manipulation
import sys  # For system operations (not used here, but often included)

# Get parameters from Snakemake (passed via params in Snakefile)
MIN_COUNT = int(snakemake.params.min_count)
MIN_SAMPLES = int(snakemake.params.min_samples)

# Snakemake provides input and output file lists
input_files = snakemake.input  # List of input count files
output_files = snakemake.output  # List of output filtered files

for in_file, out_file in zip(input_files, output_files):  # Loop through input/output pairs
    df = pd.read_csv(in_file, sep='\t', index_col=0)  # Read count file
    # Filter: keep miRNAs with count >= MIN_COUNT in at least MIN_SAMPLES samples
    filtered = df[(df >= MIN_COUNT).sum(axis=1) >= MIN_SAMPLES]  # Apply filter
    filtered.to_csv(out_file, sep='\t')  # Save filtered counts
