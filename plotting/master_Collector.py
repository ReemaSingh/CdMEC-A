# CdMEC-A: Contextual mDNA Mobile Element Classifier - Analyzer
# Copyright (C) 2025 [Dr. Reema Singh]
# Licensed under the GNU General Public License v3.0

import matplotlib
matplotlib.use('Agg')  # Essential for Linux/HPC environments
import pandas as pd
import glob
import os
import argparse
import matplotlib.pyplot as plt

def collect_data(input_dir):
    print(f"--- [STEP 1/3] Scanning Directory: {input_dir} ---")
    tsv_files = glob.glob(os.path.join(input_dir, "*_summary.tsv"))
    
    if not tsv_files:
        print("!!! ERROR: No *_summary.tsv files found. !!!")
        return None

    print(f"Found {len(tsv_files)} files. Merging now...")
    all_df = [pd.read_csv(f, sep='\t') for f in tsv_files]
    return pd.concat(all_df, ignore_index=True)

def create_enhanced_plot(df, output_prefix):
    print(f"--- [STEP 2/3] Generating Publication-Quality Plot ---")
    
    plt.figure(figsize=(12, 7))
    
    # 1. Create the Histogram
    # We use more bins (100) for better resolution near 0
    plt.hist(df['Proximity_bp'], bins=100, color='#2c7bb6', edgecolor='white', alpha=0.8)
    
    # 2. Add the "Embedded" Red Line
    plt.axvline(0, color='#d7191c', linestyle='--', linewidth=2.5, label='Embedded (0 bp)')
    
    # 3. Add Annotation Text
    plt.text(500, plt.ylim()[1]*0.9, 'Embedded Genes', color='#d7191c', 
             fontweight='bold', fontsize=12)
    
    # 4. Professional Labeling
    plt.title('ARG-MGE Genomic Proximity Distribution', fontsize=16, fontweight='bold', pad=20)
    plt.xlabel('Proximity to Mobile Genetic Element (bp)', fontsize=13)
    plt.ylabel('Frequency (Number of Hits)', fontsize=13)
    
    # 5. Styling
    plt.grid(axis='y', linestyle=':', alpha=0.7)
    plt.xlim(-10000, 10000) # Keep the 10kb window
    plt.legend(loc='upper right')
    
    # Save the file
    plot_name = f"{output_prefix}_Distribution.png"
    plt.tight_layout()
    plt.savefig(plot_name, dpi=300) # High resolution for publication
    print(f"SAVED: {plot_name}")

def main():
    parser = argparse.ArgumentParser(description="CdMEC-A HPC Collector & Visualizer")
    parser.add_argument("-i", "--input", required=True, help="Input directory with TSVs")
    parser.add_argument("-o", "--output", default="MyReport", help="Output prefix")
    args = parser.parse_args()

    master_df = collect_data(args.input)
    
    if master_df is not None:
        # Save Master CSV
        master_df.to_csv(f"{args.output}_Master.csv", index=False)
        print(f"--- [STEP 3/3] Created {args.output}_Master.csv ---")
        
        # Create Plot
        create_enhanced_plot(master_df, args.output)
        print("\nProcess Complete!")

if __name__ == "__main__":
    main()
