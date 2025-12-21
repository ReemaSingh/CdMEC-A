# CdMEC-A: Contextual mDNA Mobile Element Classifier - Analyzer
# Copyright (C) 2025 [Dr. Reema Singh]
# Licensed under the GNU General Public License v3.0

import pandas as pd
import os
import argparse

def generate_publication_report(input_csv, output_prefix):
    if not os.path.exists(input_csv):
        print(f"Error: {input_csv} not found!")
        return

    # Load data
    df = pd.read_csv(input_csv)
    
    # Strip any hidden spaces from headers just in case
    df.columns = df.columns.str.strip()

    # Define the columns based on your provided snippet
    arg_col = 'ARG_Name'
    mge_col = 'MGE_Association'  # Changed from MGE_Name to match your file
    dist_col = 'Proximity_bp'
    status_col = 'Inferred_Status'

    # 1. Filter for high-risk mobility (Distance within 1kb or Embedded)
    # Using .abs() ensures we capture both Upstream (-1000) and Downstream (+1000)
    mobile_hits = df[df[dist_col].abs() <= 1000].copy()

    # 2. Table 1: Top 10 Most Mobile ARGs
    top_args = mobile_hits[arg_col].value_counts().head(10).reset_index()
    top_args.columns = ['Resistance_Gene', 'Mobile_Occurrence_Count']

    # 3. Table 2: Top 10 MGE Carriers
    # Note: This counts unique MGE accessions found in the 'MGE_Association' column
    top_mges = mobile_hits[mge_col].value_counts().head(10).reset_index()
    top_mges.columns = ['MGE_Accession_Info', 'Total_Cargo_Genes']

    # 4. Summary Statistics for Abstract
    total_hits = len(df)
    embedded_count = len(df[df[status_col] == 'Embedded within MGE'])
    
    # --- PRINTING TO TERMINAL ---
    print("\n" + "="*60)
    print(f"   CdMEC-A SUMMARY REPORT: {output_prefix}")
    print("="*60)
    print(f"Total Associations Found: {total_hits}")
    print(f"Embedded (High Risk):     {embedded_count}")
    
    print("\n[TABLE 1: TOP MOBILE ARGs (<1kb distance)]")
    print(top_args.to_string(index=False))

    print("\n[TABLE 2: TOP MGE ASSOCIATIONS]")
    print(top_mges.to_string(index=False))

    # --- SAVE OUTPUTS ---
    arg_file = f"{output_prefix}_Top_ARGs.csv"
    mge_file = f"{output_prefix}_Top_MGEs.csv"
    
    top_args.to_csv(arg_file, index=False)
    top_mges.to_csv(mge_file, index=False)
    
    print("\n" + "="*60)
    print(f"CSV Tables saved: {arg_file}, {mge_file}")
    print("="*60)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", required=True, help="Input Master CSV")
    parser.add_argument("-o", "--output", default="Cdiff_Analysis", help="Output prefix")
    args = parser.parse_args()
    
    generate_publication_report(args.input, args.output)
