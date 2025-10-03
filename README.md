# miRNA Differential Expression & Pathway Workflow

This repository contains two alternative miRNA-seq analysis workflows, implemented in Snakemake:

- **diabetic_miRNA_4GB/** — A lightweight, parameterized pipeline including trimming, alignment, quantification, filtering, differential expression, target prediction, and summary reporting. Designed for small test data and limited disk usage (< 4 GB).
- **miRgeworkflow/** — An extended version incorporating **miRge3.0** for quantification. This workflow is more resource intensive and is experimental.

---

## Directory Structure

