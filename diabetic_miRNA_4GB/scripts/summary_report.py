# summary_report.py
# This script generates a detailed summary (CSV and HTML) of the workflow execution.
# It collects sample info, file paths, counts, DE results, and pathway results.

import os  # For file and directory operations
import pandas as pd  # For data manipulation
from datetime import datetime  # For timestamps

# Get sample names from config
samples = ["diabetic_1", "diabetic_2", "control_1", "control_2"]  # List of sample names

# Collect file paths
raw_fastq = [f"data/raw/{s}.fastq.gz" for s in samples]  # Paths to raw FASTQ files
counts = [f"results/counts/{s}.counts" for s in samples]  # Paths to count files
de_results = "results/DE_pathways/diabetic_vs_control.csv"  # Path to DE results
pathway_html = "results/DE_pathways/pathways.html"  # Path to pathway report

# Gather summary info
data = []  # List to hold summary rows
for i, s in enumerate(samples):  # Loop through samples
    row = {
        "Sample": s,  # Sample name
        "FASTQ": raw_fastq[i],  # FASTQ file path
        "Counts": counts[i],  # Count file path
        "FASTQ Exists": os.path.exists(raw_fastq[i]),  # Check if FASTQ exists
        "Counts Exists": os.path.exists(counts[i]),  # Check if counts exist
    }
    data.append(row)  # Add row to summary

summary_df = pd.DataFrame(data)  # Create DataFrame from summary rows

# Add DE and pathway info
summary_info = {
    "DE Results": de_results,  # DE results file path
    "DE Exists": os.path.exists(de_results),  # Check if DE results exist
    "Pathway Report": pathway_html,  # Pathway report file path
    "Pathway Exists": os.path.exists(pathway_html),  # Check if pathway report exists
    "Workflow Completed": datetime.now().isoformat(),  # Timestamp
}

# Save summary table as CSV
os.makedirs("results/summary", exist_ok=True)  # Ensure output directory exists
summary_csv = "results/summary/summary.csv"  # Path to summary CSV
summary_df.to_csv(summary_csv, index=False)  # Save summary DataFrame as CSV

# Save detailed HTML report
summary_html = "results/summary/summary.html"  # Path to summary HTML
with open(summary_html, "w") as f:  # Open file for writing
    f.write("<h1>miRNA Workflow Summary Report</h1>")  # Write report title
    f.write(f"<p><b>Workflow completed:</b> {summary_info['Workflow Completed']}</p>")  # Write timestamp
    f.write("<h2>Sample File Summary</h2>")  # Section header
    f.write(summary_df.to_html(index=False))  # Write summary table
    f.write("<h2>Key Results</h2>")  # Section header
    f.write(f"<p><b>Differential Expression Results:</b> {de_results} - {'Found' if summary_info['DE Exists'] else 'Missing'}</p>")  # DE results status
    f.write(f"<p><b>Pathway Analysis Report:</b> {pathway_html} - {'Found' if summary_info['Pathway Exists'] else 'Missing'}</p>")  # Pathway report status
    f.write("<h2>Parameters</h2>")  # Section header
    f.write("<ul><li>Samples: " + ", ".join(samples) + "</li>")  # List samples
    f.write("<li>DE method: pydeseq2</li><li>Pathway method: gseapy</li></ul>")  # List methods
