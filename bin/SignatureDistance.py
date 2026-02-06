# CdMEC-A: Contextual mDNA Mobile Element Classifier - Analyzer
# Copyright (C) 2025 [Dr. Reema Singh]
# Licensed under the GNU General Public License v3.0
import pandas as pd
import os
import glob
import argparse
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm
from scipy import stats 

def parse_args():
    parser = argparse.ArgumentParser(description="One Health Parallel Analysis (Signature Logic)")
    parser.add_argument("--porcine", "-p", default="animal_fna_cdmec")
    parser.add_argument("--environment", "-e", default="env_fna_cdmec")
    parser.add_argument("--human", "-u", default="human_fna_cdmec")
    parser.add_argument("--pattern", default="*.tsv")
    parser.add_argument("--workers", "-w", type=int, default=os.cpu_count())
    return parser.parse_args()

def process_single_file(file_path):
    try:
        # Load columns + file source to track prevalence correctly
        df = pd.read_csv(file_path, sep="\t", usecols=["ARG_Name", "Proximity_bp"])
        df.columns = [c.strip() for c in df.columns]
        # Track which file this came from
        df['File_Source'] = os.path.basename(file_path)
        return df
    except Exception:
        return pd.DataFrame()
def run_analysis(folders, pattern, max_workers):
    final_rows = []

    for host, path in folders.items():
        if not os.path.isdir(path):
            continue

        files = glob.glob(os.path.join(path, pattern))
        total_samples = len(files)
        if total_samples == 0:
            continue

        print(f"[*] Processing {host} ({total_samples} files)...")
        
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            results = list(tqdm(executor.map(process_single_file, files), total=total_samples))

        if not results: continue
        host_df = pd.concat(results, ignore_index=True)
        # --- REVISED CALCULATION LOGIC ---
        summary = host_df.groupby("ARG_Name").agg(
            # Distance Signature
            Signature_Distance_bp=("Proximity_bp", lambda x: stats.mode(x, keepdims=True).mode[0]),
            # Total Count (for Redundancy)
            Total_Hits=("ARG_Name", "size"),
            # Unique Files (for Prevalence)
            Samples_With_Gene=("File_Source", "nunique")
        ).reset_index()

        # 1. Prevalence %: Percentage of genomes containing the gene (Max 100%)
        summary["Prevalence_Pct"] = (summary["Samples_With_Gene"] / total_samples) * 100
        
        # 2. Avg Copies: The "Redundancy" / Dosage metric (e.g., 15.4 copies/genome)
        summary["Avg_Copies_Per_Genome"] = summary["Total_Hits"] / total_samples
        
        summary["Host"] = host
        summary.rename(columns={"ARG_Name": "Gene"}, inplace=True)
        
        final_rows.append(summary)

    if final_rows:
        output_df = pd.concat(final_rows, ignore_index=True)
        # This will now have the columns your Visualization script expects
        output_df.to_csv("one_health_spatial_signatures.csv", index=False)
        print("\n Success: 'one_health_spatial_signatures.csv' created with Prevalence and Redundancy metrics.")

if __name__ == "__main__":
    args = parse_args()
    target_folders = {"Porcine": args.porcine, "Environment": args.environment, "Human": args.human}
    run_analysis(target_folders, args.pattern, args.workers)

