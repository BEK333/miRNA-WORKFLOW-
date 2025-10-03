# miRgeworkflow/README.md

# miRge3 Quantification Workflow (Experimental)

This workflow directory contains a miRge3-based quantification block for small RNA sequencing. It is experimental and intended for demonstration purposes.

## What it Does

- Downloads or expects a test FASTQ  
- Runs **miRge3.0** quantification  
- Produces per-sample `merged_counts.csv`  
- Can generate a summary counts table

## How to Run

1. Configure `config.yaml` (if present) — you may leave defaults for test fetch  
2. Ensure your conda environment includes miRge3  
3. From this folder:
   ```bash
   snakemake --use-conda --cores 2
   ```

## Notes & Caveats

- This is **not** a full differential expression pipeline  
- It is meant as **portfolio proof** of handling miRge3  
- Resource demands are modest for small test data — may not scale  
- You can extend later to DE, plotting, target prediction
