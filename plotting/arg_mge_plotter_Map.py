# CdMEC-A: Contextual mDNA Mobile Element Classifier - Analyzer
# Copyright (C) 2025 [Dr. Reema Singh]
# Licensed under the GNU General Public License v3.0

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.lines import Line2D
import argparse
import os
import glob

def load_data(file_path):
    """Robustly loads data and identifies ARG vs MGE status."""
    records = []
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
            for line in lines[1:]:
                t = line.split()
                if len(t) >= 8:
                    # Logic to distinguish ARG from MGE based on 'Inferred_Status'
                    # or proximity columns typical in cdmec outputs
                    status = t[7].upper() 
                    is_mge = "MGE" in status or "MOBILE" in status
                    records.append({
                        'Name': t[2], 
                        'Start': int(t[3]), 
                        'End': int(t[4]),
                        'Type': 'MGE' if is_mge else 'ARG'
                    })
        df = pd.DataFrame(records)
        df['Label'] = df['Name'].apply(lambda x: str(x).split('|')[-1].split('_')[0])
        return df
    except Exception as e:
        print(f"Skipping {file_path}: {e}")
        return None

def plot_simultaneous(df, sample_id, start_bp, end_bp, output_path):
    """Generates a fragmented map for a specific window with a legend."""
    df_win = df[(df['End'] >= start_bp) & (df['Start'] <= end_bp)].copy()
    
    num_rows = 10 
    bp_per_row = (end_bp - start_bp) / num_rows
    
    fig, axes = plt.subplots(num_rows, 1, figsize=(24, 4 * num_rows))
    plt.subplots_adjust(hspace=1.4) 

    # Color mapping for the legend
    color_map = {'ARG': '#e74c3c', 'MGE': '#3498db'} # Red for ARG, Blue for MGE

    for i in range(num_rows):
        ax = axes[i]
        r_start = start_bp + (i * bp_per_row)
        r_end = r_start + bp_per_row
        
        row_df = df_win[(df_win['Start'] >= r_start - 100) & (df_win['End'] <= r_end + 100)]
        ax.hlines(50, r_start, r_end, color='#7f8c8d', linewidth=10, alpha=0.5)
        
        sorted_genes = row_df.sort_values('Start').to_dict('records')
        for j, row in enumerate(sorted_genes):
            # Up/Down Alternating Tiers (Total 8 height levels)
            is_up = (j % 2 == 0)
            tier = (j // 2) % 4
            y_offset = 20 + (tier * 18)
            y_text = 50 + y_offset if is_up else 50 - y_offset
            
            # Feature Box
            color = color_map.get(row['Type'], 'gray')
            ax.add_patch(patches.Rectangle((row['Start'], 48), max(100, row['End']-row['Start']), 4, 
                                          facecolor=color, edgecolor='black', zorder=5))
            
            # Anchor Line
            ax.plot([row['Start'], row['Start']], [50, y_text], 'k:', lw=1, alpha=0.3)
            
            # Label
            ax.text(row['Start'], y_text, row['Label'], fontsize=16, weight='bold',
                    ha='center', va='bottom' if is_up else 'top',
                    bbox=dict(facecolor='white', edgecolor='none', alpha=0.8, pad=0.5))

        ax.set_xlim(r_start, r_end)
        ax.set_ylim(0, 140)
        ax.set_yticks([])
        ax.xaxis.set_major_formatter(plt.FuncFormatter(lambda x, p: format(int(x), ',')))
        ax.set_title(f"Segment: {int(r_start):,} - {int(r_end):,} bp", loc='right', fontsize=12)

    # Add the legend similar to your attached figure
    legend_elements = [
        patches.Patch(facecolor='#e74c3c', edgecolor='black', label='ARG (Antibiotic Resistance Gene)'),
        patches.Patch(facecolor='#3498db', edgecolor='black', label='MGE (Mobile Genetic Element)')
    ]
    fig.legend(handles=legend_elements, loc='lower center', ncol=2, fontsize=20, bbox_to_anchor=(0.5, 0.02))

    plt.suptitle(f"COMPARATIVE GENOMIC MAP: {sample_id}\nWindow: {start_bp:,} - {end_bp:,} bp", 
                 fontsize=32, y=0.98, weight='bold')
    
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_dir", required=True, help="Directory with .txt/.tsv files")
    parser.add_argument("--start", type=int, required=True, help="Fixed Start (bp)")
    parser.add_argument("--end", type=int, required=True, help="Fixed End (bp)")
    parser.add_argument("-o", "--output_dir", default="Simultaneous_Results")
    args = parser.parse_args()

    if not os.path.exists(args.output_dir): os.makedirs(args.output_dir)
    
    files = glob.glob(os.path.join(args.input_dir, "*.txt")) + glob.glob(os.path.join(args.input_dir, "*.tsv"))
    
    for f in files:
        name = os.path.basename(f).split('.')[0]
        print(f"Generating simultaneous map for: {name}")
        data = load_data(f)
        if data is not None:
            out_path = os.path.join(args.output_dir, f"{name}_Comparison.png")
            plot_simultaneous(data, name, args.start, args.end, out_path)
