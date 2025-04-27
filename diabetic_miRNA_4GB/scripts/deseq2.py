# deseq2.py
# This script performs differential expression analysis using pydeseq2.
# It is intended to be run by Snakemake, which provides the input count files and output path.

import pandas as pd  # Import pandas for data manipulation
from pydeseq2 import DeseqDataSet  # Import DESeq2 wrapper

# Get padj_threshold from Snakemake params
padj_threshold = float(snakemake.params.padj_threshold)

# Read all count files (one per sample) and concatenate them into a single DataFrame.
# Each file is expected to be a tab-separated table with gene IDs as rows.
counts = pd.concat(
    [pd.read_csv(f, sep='\t', index_col=0) for f in snakemake.input],  # Read each count file
    axis=1  # Concatenate along columns (samples)
)

# Create metadata for each sample, assigning 'diabetic' or 'control' based on sample name.
metadata = pd.DataFrame({
    'condition': ['diabetic' if 'diabetic' in col else 'control' for col in counts.columns]  # Assign condition
})

# Initialize DESeq2 dataset with counts and metadata.
dds = DeseqDataSet(counts=counts.T, metadata=metadata, design_factors='condition')  # Transpose counts for DESeq2

# Run the DESeq2 analysis.
dds.deseq2()  # Perform differential expression analysis

# Extract results (log2 fold change, p-values, etc.) and save to CSV.
results = dds.results()  # Get results DataFrame
# Filter by padj_threshold parameter
results = results[results['padj'] < padj_threshold]
results.to_csv(snakemake.output[0])  # Save results to output file
