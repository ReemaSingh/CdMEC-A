# CdMEC-A: Contextual mDNA Mobile Element Classifier - Analyzer
# Copyright (C) 2025 [Dr. Reema Singh]
# Licensed under the GNU General Public License v3.0

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os

def visualize_arg_mge_robust(file_path, output_dir):
    """
    Loads ARG-MGE association data and saves plots to the specified output directory.
    """
    print(f"\n--- Processing: {os.path.basename(file_path)} ---")
    
    base_name = os.path.splitext(os.path.basename(file_path))[0]
    
    # Ensure the output directory exists
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        print(f"Created directory: {output_dir}")

    column_names = [
        'Sample_ID', 'Contig_ID', 'ARG_Name', 'ARG_Start', 'ARG_End',
        'MGE_Association', 'Proximity_bp', 'Inferred_Status'
    ]
    
    records = [] 
    try:
        with open(file_path, 'r') as f:
            lines = f.readlines()
            if not lines: return
            for line in lines[1:]: 
                tokens = line.split() 
                if len(tokens) >= 8:
                    record = tokens[:7]
                    inferred_status = ' '.join(tokens[7:])
                    record.append(inferred_status)
                    records.append(record)
        
        df = pd.DataFrame(records, columns=column_names)
    except Exception as e:
        print(f"Error processing {file_path}: {e}")
        return

    # --- Data Cleaning ---
    try:
        df['Proximity_bp'] = pd.to_numeric(df['Proximity_bp'], errors='coerce')
        df.dropna(subset=['Proximity_bp'], inplace=True)
        df['Gene_Name'] = df['ARG_Name'].apply(lambda x: x.split('|')[-1] if '|' in x else x)
    except Exception as e:
        print(f"Skipping {file_path} due to cleaning error: {e}")
        return

    # --- Visualization 1: Bar Chart ---
    if 'Inferred_Status' in df.columns and not df.empty:
        plt.figure(figsize=(10, 6))
        sns.countplot(
            y='Inferred_Status',
            data=df,
            order=df['Inferred_Status'].value_counts().index,
            palette='viridis'
        )
        plt.title(f'ARG-MGE Status: {base_name}', fontsize=14)
        plt.tight_layout()
        
        # Save to output folder
        out_path1 = os.path.join(output_dir, f"{base_name}_status_bar.png")
        plt.savefig(out_path1)
        plt.close()
        print(f"Saved: {out_path1}")

    # --- Visualization 2: Distribution ---
    if not df.empty:
        plt.figure(figsize=(10, 6))
        df['Proximity_Type'] = df['Proximity_bp'].apply(
            lambda x: 'Downstream (Positive)' if x > 0 else ('Embedded (Zero)' if x == 0 else 'Upstream (Negative)')
        )

        sns.histplot(
            df, x='Proximity_bp', hue='Proximity_Type',
            multiple='stack', bins=15, kde=False,
            palette={'Downstream (Positive)': 'skyblue', 'Upstream (Negative)': 'salmon', 'Embedded (Zero)': 'lightgreen'}
        )
        plt.axvline(x=0, color='grey', linestyle='--', linewidth=1)
        plt.title(f'Proximity Distribution: {base_name}', fontsize=14)
        plt.tight_layout()
        
        # Save to output folder
        out_path2 = os.path.join(output_dir, f"{base_name}_proximity_dist.png")
        plt.savefig(out_path2)
        plt.close()
        print(f"Saved: {out_path2}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Batch visualize ARG-MGE data and save to a folder.")
    parser.add_argument('-i', '--input_dir', type=str, required=True, help="Directory with input files.")
    parser.add_argument('-o', '--output_dir', type=str, required=True, help="Directory to save plots.")

    args = parser.parse_args()

    if os.path.isdir(args.input_dir):
        for filename in os.listdir(args.input_dir):
            full_path = os.path.join(args.input_dir, filename)
            if os.path.isfile(full_path):
                visualize_arg_mge_robust(full_path, args.output_dir)
    else:
        print(f"Error: {args.input_dir} is not a valid directory.")
