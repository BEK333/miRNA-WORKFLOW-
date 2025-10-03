# README

# miRNA Workflow Prototype



## Structure

```
/
├── README.md
├── Snakefile
├── config.yaml
├── environment.yaml
├── miRgeworkflow/
└── (other supporting files)
```

## How to Use Prototype

1. **Edit `config.yaml`** to point to test FASTQ or leave defaults for fetch step  
2. **Activate conda/env** that includes miRge3 and dependencies  
3. **Run pipeline**:
   ```bash
   snakemake --cores 2 --use-conda
   ```

4. Outputs:
   - `results/<sample>/merged_counts.csv`
   - `results/summary_counts.csv`

This is not a full production pipeline; it’s intended to **show your workflow skills** in your CV or GitHub.

---

## Tips & Notes

- It’s minimal — only fetch, quantification, summary.
- You can extend later to differential expression, plotting, target prediction, etc.
- Mark as prototype or experimental in the README to set expectations.

---

# config.yaml

mirge3:
  lib_dir: MiRge3_Lib/
  species: hsa
  params: "--umi --isomiR"

threads: 2

---

# Snakefile

import pandas as pd
from pathlib import Path

SAMPLES = ["demo"]

rule all:
    input:
        expand("results/{s}/merged_counts.csv", s=SAMPLES),
        "results/summary_counts.csv"

rule fetch_dummy:
    output:
        "data/demo.fastq.gz"
    shell:
        """
        mkdir -p data
        wget -O {output} https://sra-download.ncbi.nlm.nih.gov/traces/sra43/SRR/006868/SRR006868/SRR006868.fastq.gz
        """

rule miRge_quant:
    input:
        fastq="data/{s}.fastq.gz"
    output:
        counts="results/{s}/merged_counts.csv"
    threads: 2
    conda: "environment.yaml"
    shell:
        """
        mkdir -p results/{wildcards.s}
        mirge3-cli -s hsa -l MiRge3_Lib/ -1 {input.fastq} -2 {input.fastq} --output results/{wildcards.s} --umi --isomiR
        mv results/{wildcards.s}/counts_merged.csv {output.counts}
        """

rule summarize:
    input:
        expand("results/{s}/merged_counts.csv", s=SAMPLES)
    output:
        "results/summary_counts.csv"
    run:
        dfs = []
        for f in input:
            df = pd.read_csv(f, index_col=0)
            dfs.append(df.sum(axis=0))
        combo = pd.concat(dfs, axis=0)
        combo.to_csv(output[0], header=["total_counts"])

---

# environment.yaml

name: pipeline_env
channels:
  - defaults
  - conda-forge
  - bioconda
dependencies:
  - python=3.9
  - snakemake>=9.3.2
  - pandas
  - mirge3
  - wget
  - bowtie
  - samtools
  # Add other dependencies as needed
