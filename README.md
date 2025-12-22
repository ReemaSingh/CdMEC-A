# CdMEC-A: Contextual mDNA Mobile Element Classifier - Analyzer

CdMEC-A is a bioinformatic pipeline designed to analyze the genomic proximity between Antibiotic Resistance Genes (ARGs) and Mobile Genetic Elements (MGEs) in Clostridioides difficile. It provides a high-throughput workflow from raw FASTA sequences to publication-quality genomic maps.

## License
This project is licensed under the **GNU General Public License v3.0**. See the [LICENSE](LICENSE) file for details.

## Prerequisites
- **Python 3.7+**
- **NCBI BLAST+** (blastn and blastx must be in your PATH)

## Quick Start

### 1. Installation
```bash
git clone [https://github.com/yourusername/CdMEC-A.git](https://github.com/yourusername/CdMEC-A.git)
cd CdMEC-A
pip install -r requirements.txt
```

### 2. Database Preparation
CdMEC-A requires specific BLAST databases to function. Follow these steps to prepare them:
#### 1. Antibiotic Resistance Gene (ARG) Database
- Download protein_fasta_protein_homolog_model.fasta from the CARD website.
- Build the database:
```bash
makeblastdb -in protein_fasta_protein_homolog_model.fasta -dbtype prot -out card_protein_homolog_db
```
#### 2. Mobile Genetic Element (MGE) Database
This pipeline uses a curated set of C. difficile MGEs.
- Combine your MGE fasta sequences into a single file:
```bash
cat *.fasta > combined_C_Diff_mge_nucl.fasta
```
- Build the database:
```bash
makeblastdb -in combined_C_Diff_mge_nucl.fasta -dbtype nucl -out combined_C_Diff_mge_nucl_db -title "C. Difficile MGE Nucleotide Database"
```

### 3. Run Analysis
Place your FASTA files in a folder and run:
```bash
python bin/cdmec_analyzer.py -i ./test_samples -o ./results
```
### 4. Generate Reports & Plots
```bash
# Merge results and create distribution plot
python plotting/master_Collector.py -i ./results -o Study_Summary
# Generate high-risk mobility stats
python bin/cdmec_stats_generator.py -i Study_Summary.csv -o Final_Report
```

## Visualizations
CdMEC-A includes a powerful visualization suite to transition from raw data to publication-ready figures

### 1. Global Proximity Distribution
Generated via master_Collector.py, this plot visualizes the spatial relationship between all identified ARGs and MGEs across your entire dataset. It highlights whether genes are predominantly "embedded" within elements or clustered at specific distances.
```bash
python plotting/master_Collector.py -i ./results -o Global_Summary.
```
### 2. Comparative Genomic Mapping
Generated via arg_mge_plotter_Map.py, this script creates high-resolution, fragmented maps of specific contigs. It uses a tiered visualization system to prevent overlapping labels and clearly distinguishes between ARGs (Red) and MGEs (Blue).
```bash
# Example: Map a specific 50kb window on a contig
python plotting/arg_mge_plotter_Map.py -i ./results --start 1000 --end 50000 -o ./maps
```
### 3. Sample-Level Histograms
For a quick look at individual sample results, use arg_mge_plotter.py to generate density plots showing upstream vs. downstream distributions.
```bash
python plotting/arg_mge_plotter.py -i ./results -o ./sample_plots
```

## Citation
If you use CdMEC-A in your research, please cite:

Reema Singh (2026). CdMEC-A: A pipeline for ARG-MGE contextual analysis.





