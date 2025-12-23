# CdMEC-A: Contextual mDNA Mobile Element Classifier - Analyzer

CdMEC-A is a bioinformatic pipeline designed to analyze the genomic proximity between Antibiotic Resistance Genes (ARGs) and Mobile Genetic Elements (MGEs) in *Clostridioides difficile*. It provides a high-throughput workflow from raw FASTA sequences to publication-quality genomic maps.

## Prerequisites
- **Python 3.7+**
- **NCBI BLAST+** (blastn and blastx must be in your PATH)

## Quick Start

### 1. Installation and Setup
#### **Clone the repository:**
```bash
git clone https://github.com/ReemaSingh/CdMEC-A.git
cd CdMEC-A
```
#### Create and activate a virtual environment (Recommended):
```bash
python -m venv venv
source venv/bin/activate
```
#### Install dependencies:
```bash
pip install -r requirements.txt
```

### 2. Database Preparation
CdMEC-A requires specific BLAST databases to function. Follow these steps to prepare them:
#### 1. Antibiotic Resistance Gene (ARG) Database
- Download the latest **Protein Homolog Model** fasta file from the [CARD Download Page](https://card.mcmaster.ca/latest/data).
  * Filename: `protein_fasta_protein_homolog_model.fasta`
- Place the file in your root directory and build the database:
```bash
makeblastdb -in protein_fasta_protein_homolog_model.fasta -dbtype prot -out card_protein_homolog_db
```
#### 2. Mobile Genetic Element (MGE) Database
This pipeline uses a curated set of *C. difficile* specific mobile elements.

**Option A: Use Provided Sequences (Recommended)**
I have included the individual FASTA files in the `data/mge_references/` folder. To build the database, run:

```bash
# 1. Combine the individual FASTA files into one master reference
cat data/mge_references/*.fasta > combined_C_Diff_mge_nucl.fasta

# 2. Build the BLAST database
makeblastdb -in combined_C_Diff_mge_nucl.fasta -dbtype nucl -out combined_C_Diff_mge_nucl_db -title "C. Difficile MGE Nucleotide Database"
```
Option B: Manual Data Sourcing If you wish to update the sequences or verify the sources, the curated reference set includes:

* **Transposons:** Tn916-like, Tn6194, Tn6218-like, Tn6000, Tn5397, Tn5398, Tn4453a-b, Tn6944.
* **Integrons/Islands:** CdISt1 (AJ579717.1), CFR_region, and various IS family elements (IS30, ISL3).
* **Accessions:** AM180356.2, AY350745.1, MG973074.1.

### 3. Run Analysis
Place your FASTA files in a folder and run:
```bash
python bin/cdmec_analyzer.py -i ./test_samples -o ./results
```
### Parallelization Options
CdMEC-A allows for dual-layer parallelization to optimize throughput:

* **`-t` (Sample Threads):** Number of genomic samples to process simultaneously. 
* **`-bt` (BLAST Threads):** Number of CPU cores assigned to *each* BLAST process (the `-num_threads` parameter in BLAST+).

**Example:**
To process 2 samples at a time, using 4 CPUs for each BLAST search (Total 8 CPUs):
```bash
python bin/cdmec_analyzer.py -i ./test_samples -o ./results -t 2 -bt 4
```
### 4. Generate Reports & Plots
```bash
# Merge results and create distribution plot
python plotting/master_Collector.py -i ./results -o Study_Summary
# Generate high-risk mobility stats
python bin/cdmec_stats_generator.py -i Study_Summary_Master.csv -o Final_Report
# Generate summary statistics for publications
python bin/cdmec_reporter.py -i ./results
```

## Visualizations
CdMEC-A includes a powerful visualization suite to transition from raw data to publication-ready figures

### 1. Global Proximity Distribution
Generated via master_Collector.py, this plot visualizes the spatial relationship between all identified ARGs and MGEs across your entire dataset. It highlights whether genes are predominantly "embedded" within elements or clustered at specific distances.
```bash
python plotting/master_Collector.py -i ./results -o Global_Summary
```
### 2. Comparative Genomic Mapping
Generated via arg_mge_plotter_Map.py, this script creates high-resolution, fragmented maps of specific contigs. It uses a tiered visualization system to prevent overlapping labels and clearly distinguishes between ARGs (Red) and MGEs (Blue).
```bash
# Example: Map a specific 50kb window on a contig
python plotting/arg_mge_plotter_Map.py -i ./results --start 1000 --end 50000 -o ./maps
```
**Note:** When mapping conserved regions (e.g., common transposon backbones), plots from different samples may appear identical if the genetic architecture is conserved across the selected coordinate range.
### 3. Sample-Level Histograms
For a quick look at individual sample results, use arg_mge_plotter.py to generate density plots showing upstream vs. downstream distributions.
```bash
python plotting/arg_mge_plotter.py -i ./results -o ./sample_plots
```

## Troubleshooting & Special Environments

### Running on HPC Clusters (e.g., Compute Canada/Slurm)
If you are running this on a cluster that restricts internet access on compute nodes, load the environment before installing:

```bash
module load python/3.10
pip install --no-index -r requirements.txt
```

## Citation
If you use CdMEC-A in your research, please cite:

Reema Singh (2026). CdMEC-A: A pipeline for ARG-MGE contextual analysis.

## License
This project is licensed under the **GNU General Public License v3.0**. See the [LICENSE](LICENSE) file for details.

## Contact
This is a personal side project. Feedback and contributions are welcome! Please open an issue or contact me at reemasingh.gencompbio@gmail.com






