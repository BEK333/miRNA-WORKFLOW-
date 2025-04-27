# pathways.py
# This script performs pathway enrichment analysis using gseapy.
# It is intended to be run by Snakemake, which provides the input DE results and output path.

import pandas as pd  # For data manipulation
import gseapy as gp  # For pathway enrichment
import os  # For file path operations

# Get padj_threshold from Snakemake params
padj_threshold = float(snakemake.params.padj_threshold)

# Load differential expression results from CSV.
de_results = pd.read_csv(snakemake.input[0])  # Read DE results file

# Select genes with adjusted p-value < padj_threshold as significant.
targets = de_results[de_results['padj'] < padj_threshold]['gene'].tolist()  # Get significant genes

# Run pathway enrichment analysis using Enrichr (KEGG and Reactome gene sets).
enr = gp.enrichr(
    gene_list=targets,  # List of significant genes
    gene_sets=['KEGG_2021', 'Reactome_2022'],  # Databases to use
    outdir=os.path.dirname(snakemake.output[0])  # Output directory for results
)

# Save the enrichment results as an HTML report.
enr.results.to_html(snakemake.output[0])  # Save results as HTML
