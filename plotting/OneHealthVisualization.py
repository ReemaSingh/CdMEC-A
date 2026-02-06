# CdMEC-A: Contextual mDNA Mobile Element Classifier - Analyzer
# Copyright (C) 2025 [Dr. Reema Singh]
# Licensed under the GNU General Public License v3.0

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os

def main():
    parser = argparse.ArgumentParser(description="CdMEC-A One Health Visualization Suite")
    parser.add_argument("input", nargs="?", default="one_health_spatial_signatures.csv")
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"Error: {args.input} not found. Ensure you ran the updated SignatureDistance.py first.")
        return

    # Load data
    df = pd.read_csv(args.input)
    
    # Extract short gene name for cleaner labels (e.g., 'ranA')
    df['Gene_Short'] = df['Gene'].apply(lambda x: str(x).split('|')[-1] if '|' in str(x) else str(x))

    # --- SELECTION LOGIC ---
    # Identify genes common to all 3 hosts for the core comparison heatmaps
    common_genes = df.groupby('Gene_Short').filter(lambda x: x['Host'].nunique() == 3)['Gene_Short'].unique()
    
    if len(common_genes) == 0:
        print("Warning: No common genes found across all 3 hosts. Using top 25 by overall prevalence.")
        top_genes = df.groupby('Gene_Short')['Prevalence_Pct'].mean().sort_values(ascending=False).head(25).index
    else:
        # Sort common genes by their average redundancy (dosage) to highlight 'super-signatures'
        top_genes = (df[df['Gene_Short'].isin(common_genes)]
                     .groupby('Gene_Short')['Avg_Copies_Per_Genome']
                     .mean().sort_values(ascending=False).head(25).index)
        
  # ---------------- FIGURE: PREVALENCE HEATMAP (%) ----------------
    prev_pivot = df.pivot_table(index='Gene_Short', columns='Host', values='Prevalence_Pct').fillna(0)
    plt.figure(figsize=(10, 12))
    sns.heatmap(prev_pivot.loc[top_genes], annot=True, fmt=".1f", cmap="YlGnBu", vmax=100,
                cbar_kws={'label': 'Prevalence (%)'})
    plt.title('Top Conserved Signatures - Prevalence (%)')
    plt.tight_layout()
    plt.savefig('prevalence.png')

    # ---------------- FIGURE: REDUNDANCY HEATMAP (DOSAGE) ----------------
    redu_pivot = df.pivot_table(index='Gene_Short', columns='Host', values='Avg_Copies_Per_Genome').fillna(0)
    plt.figure(figsize=(10, 12))
    sns.heatmap(redu_pivot.loc[top_genes], annot=True, fmt=".1f", cmap="YlOrRd",
                cbar_kws={'label': 'Avg Copies per Genome'})
    plt.title('Top Conserved Signatures - Redundancy (Dosage)')
    plt.tight_layout()
    plt.savefig('redundancy.png')

    # ---------------- FIGURE: GENOMIC DISTANCE HEATMAP (bp) ----------------
    dist_pivot = df.pivot_table(index='Gene_Short', columns='Host', values='Signature_Distance_bp').fillna(0)
    plt.figure(figsize=(10, 12))
    # RdBu_r highlights 0bp (internalized) as center/white, distant as blue/red
    sns.heatmap(dist_pivot.loc[top_genes], annot=True, fmt=".0f", cmap="RdBu_r", center=0,
                cbar_kws={'label': 'Distance to MGE (bp)'})
    plt.title('Spatial Signatures - Genomic Distance (bp)')
    plt.tight_layout()
    plt.savefig('distance.png')
    
    # ---------------- FIGURE: KEY ONE HEALTH MARKER BARS ----------------
    key_markers = ['vanA', 'vanC', 'tet(S)', 'tet(T)', 'ranA', 'bcrA', 'poxtA']
    key_df = df[df['Gene_Short'].isin(key_markers)]

    if not key_df.empty:
        # Prevalence Comparison
        plt.figure(figsize=(12, 6))
        sns.barplot(data=key_df, x='Gene_Short', y='Prevalence_Pct', hue='Host', palette='viridis')
        plt.axhline(100, color='red', linestyle='--', alpha=0.6, label='100% Fixation')
        plt.title('Prevalence of Priority One Health Markers')
        plt.ylabel('Prevalence (%)')
        plt.ylim(0, 120)
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig('key_prevalence.png')

        # Redundancy Comparison
        plt.figure(figsize=(12, 6))
        sns.barplot(data=key_df, x='Gene_Short', y='Avg_Copies_Per_Genome', hue='Host', palette='magma')
        plt.axhline(1, color='black', linestyle='-', alpha=0.3, label='Single Copy Baseline')
        plt.title('Genomic Redundancy of Priority Markers')
        plt.ylabel('Avg Copies Per Genome')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig('key_redundancy.png')

    print("\n Visualization Suite Complete!")
    print("Generated: prevalence.png, redundancy.png, distance.png, key_prevalence.png, key_redundancy.png")

if __name__ == "__main__":
    main()
