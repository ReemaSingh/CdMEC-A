# CdMEC-A: Clostridioides difficile Mobile Element Context Analyzer
# Copyright (C) 2025-2026 [Dr. Reema Singh]
# Licensed under the GNU General Public License v3.0

import matplotlib
matplotlib.use('Agg')  # Essential for Linux/HPC environments
import pandas as pd
import glob
import os
import argparse
import numpy as np
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
    print(f"--- [STEP 2/3] Generating Dual-Panel Publication-Quality Plot ---")
    
    # Pre-process classifications and handle string variations or missing data safely
    df['Inferred_Status'] = df['Inferred_Status'].astype(str).str.strip()
    
    # Standardize names to match your cluster statistics output format precisely
    status_counts = {
        'Chromosomal (Unlinked)': df['Inferred_Status'].str.contains('Chromosomal', case=False, na=False).sum(),
        'Embedded within MGE': df['Inferred_Status'].str.contains('Embedded', case=False, na=False).sum(),
        'MGE-Associated (Upstream)': df['Inferred_Status'].str.contains('Upstream', case=False, na=False).sum(),
        'MGE-Associated (Downstream)': df['Inferred_Status'].str.contains('Downstream', case=False, na=False).sum()
    }
    
    # Create a 1-row, 2-column figure layout for a comprehensive genomic overview
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))
    
    # -------------------------------------------------------------------------
    # PANEL A: Categorical Bar Chart (Answers the demand for explicit counts)
    # -------------------------------------------------------------------------
    categories = list(status_counts.keys())
    counts = list(status_counts.values())
    colors = ['#7f7f7f', '#d7191c', '#2c7bb6', '#fdae61'] # Professional publication palette
    
    bars = ax1.bar(categories, counts, color=colors, edgecolor='black', alpha=0.85, width=0.6)
    
    # Add explicit data counts and percentages on top of every single bar
    total_args = sum(counts) if sum(counts) > 0 else 1
    for bar in bars:
        height = bar.get_height()
        percentage = (height / total_args) * 100
        ax1.annotate(f'{int(height)}\n({percentage:.1f}%)',
                     xy=(bar.get_x() + bar.get_width() / 2, height),
                     xytext=(0, 5),  # 5 points vertical offset
                     textcoords="offset points",
                     ha='center', va='bottom', fontsize=10, fontweight='bold')
                     
    ax1.set_title('A: Categorical Structural Classification Counts', fontsize=13, fontweight='bold', pad=15)
    ax1.set_ylabel('Total ARG Count (Frequency)', fontsize=12)
    
    # FIX: Explicitly set the tick positions first to silence the UserWarning permanently
    ax1.set_xticks(range(len(categories)))
    ax1.set_xticklabels(categories, rotation=25, ha='right', fontsize=10)
    
    ax1.grid(axis='y', linestyle=':', alpha=0.5)
    
    # Adjust upper limit to make room for annotations comfortably
    ax1.set_ylim(0, max(counts) * 1.15)

    # -------------------------------------------------------------------------
    # PANEL B: Flanking Distance Proximity Distribution (Fixes the Squashed Line)
    # -------------------------------------------------------------------------
    # Convert proximity column to numeric matrix elements, forcing unlinked tracking strings to NaN
    df['Proximity_bp_numeric'] = pd.to_numeric(df['Proximity_bp'], errors='coerce')
    
    # Biological Filtering Constraint: Exclude unlinked (NaN) AND embedded (0 bp) elements
    # This isolates flanking neighborhoods to show true spatial variations
    flanking_distances = df['Proximity_bp_numeric'].dropna()
    flanking_distances = flanking_distances[flanking_distances != 0]
    
    if len(flanking_distances) > 0:
        # Use absolute distance boundaries for flanking visualization
        # Set up signed values so upstream can be visualized on the left (-) and downstream on the right (+)
        signed_distances = []
        for _, row in df.iterrows():
            stat = row['Inferred_Status']
            val = row['Proximity_bp_numeric']
            if pd.notna(val) and val != 0:
                if 'Upstream' in stat:
                    signed_distances.append(-abs(val)) # Represent upstream on the left spectrum
                elif 'Downstream' in stat:
                    signed_distances.append(abs(val))  # Represent downstream on the right spectrum
        
        # Plot flanking density map using 40 fine-resolution bins
        n, bins, patches = ax2.hist(signed_distances, bins=40, color='#2c7bb6', edgecolor='white', alpha=0.8, label='Flanking Windows')
        ax2.set_xlim(-10000, 10000)
        ax2.axvline(0, color='black', linestyle='-', linewidth=1.2, alpha=0.7) # Operational window separator
         
        # This compresses the visual height of the data bars, forcing them below the text zone
        y_ceiling = ax2.get_ylim()[1]
        ax2.set_ylim(0, y_ceiling * 1.35)
        
        # Recalculate fixed ceiling coordinate after adjustment
        adjusted_y_ceiling = ax2.get_ylim()[1]
        
        # Using 0.95 ensures it hangs perfectly inside the clean top lane
        ax2.text(-9500, adjusted_y_ceiling * 0.95, '← Upstream Flanking', color='#2c7bb6', 
                  fontweight='bold', fontsize=11, va='top', ha='left')
        ax2.text(9500, adjusted_y_ceiling * 0.95, 'Downstream Flanking →', color='#fdae61', 
                  fontweight='bold', fontsize=11, va='top', ha='right')
        
    else:
        ax2.text(0.5, 0.5, 'No Non-Zero Flanking Loci\nDetected Within 10-kb Window', 
                 ha='center', va='center', transform=ax2.transAxes, fontsize=12, color='gray')

    ax2.set_title('B: Fine-Resolution Flanking Distance Spectrum', fontsize=12, fontweight='semibold', pad=10)
    ax2.set_xlabel('Spatial Distance to Reference Target (bp)', fontsize=12)
    ax2.set_ylabel('Frequency (Number of Hits)', fontsize=12)
    ax2.grid(axis='y', linestyle=':', alpha=0.5)   

    # Global Master Styling Layout
    plt.suptitle('ARG-MGE Genomic Proximity Distribution', fontsize=16, fontweight='bold', y=0.98)
    
    plot_name = f"{output_prefix}_Distribution.png"
    plt.tight_layout()
    plt.savefig(plot_name, dpi=300) # Crisp, high-resolution file for PeerJ print proofs
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
