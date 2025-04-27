# miRNA Differential Expression and Pathway Analysis Workflow

This project provides a fully parameterized, Snakemake-based workflow for miRNA-seq data analysis, including adapter trimming, alignment, quantification, filtering, differential expression, target prediction, pathway enrichment, and summary reporting.

## Project Structure
- `Snakefile`: Main workflow definition (thoroughly commented and parameterized)
- `config/config.yaml`: All workflow parameters, sample names, and reference file locations
- `envs/mirna.yaml`: Conda environment for reproducibility
- `scripts/`: Python scripts for data download, filtering, DE analysis, target prediction, pathway analysis, and summary reporting
- `results/`: Output directory for results
- `data/raw/`: Raw FASTQ files (downloaded)

## Key Features
- **Parameterization:** All major workflow settings (adapter, threads, filtering, p-value thresholds, TargetScan options) are set in `config/config.yaml` for easy modification.
- **Filtering:** Low-count miRNAs are filtered out before DE analysis, with thresholds set in the config file.
- **Target Prediction:** Significant miRNAs are analyzed for mRNA targets using TargetScan.
- **Summary Reporting:** Generates both CSV and HTML summary reports of workflow execution and results.
- **Reproducibility:** Uses conda environments and Snakemake best practices.
- **Extensive Comments:** All scripts and the workflow are line-by-line commented for clarity and future maintainability.

## Setup Instructions
1. **Clone this repository** (if not already):
   ```sh
   git clone <your-repo-url>
   cd diabetic_miRNA_4GB
   ```
2. **Download miRBase index**:
   ```sh
   mkdir -p index
   wget ftp://mirbase.org/pub/mirbase/CURRENT/hairpin.fa.gz -O index/miRBase_v22.fa.gz
   bowtie-build index/miRBase_v22.fa.gz index/miRBase_v22
   ```
3. **Edit `config/config.yaml`** to set your samples, adapter, filtering, and analysis parameters as needed.
4. **Create the conda environment**:
   ```sh
   conda env create -n mirna-py39 -f envs/mirna.yaml
   conda activate mirna-py39
   ```
5. **Run the workflow**:
   ```sh
   snakemake --use-conda --cores 4
   ```

## Notes
- The workflow is designed for <4GB disk usage with partial test data.
- Update accession numbers in `scripts/download_data.py` as needed.
- All scripts and the Snakefile are thoroughly commented for learning and reproducibility.
- For real target prediction, ensure TargetScan and required UTR files are available and paths are set in `config/config.yaml`.
- For questions or improvements, see the comments in each script and the config file.
