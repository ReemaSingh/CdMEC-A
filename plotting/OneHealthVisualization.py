# CdMEC-A: Contextual mDNA Mobile Element Classifier - Analyzer
# Copyright (C) 2025 [Dr. Reema Singh]
# Licensed under the GNU General Public License v3.0
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("input", nargs="?", default="one_health_spatial_signatures.csv")
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"Error: {args.input} not found.")
        return

    df = pd.read_csv(args.input)
    # Extract short gene name (e.g., 'vanA')
    df['Gene_Short'] = df['Gene'].apply(lambda x: str(x).split('|')[-1] if '|' in str(x) else str(x))

    # Identify genes common to all 3 hosts for the 'Top 25' comparison
    common_genes = df.groupby('Gene_Short').filter(lambda x: x['Host'].nunique() == 3)['Gene_Short'].unique()
    
    if len(common_genes) == 0:
        print("Warning: No common genes found across all hosts. Using top available.")
        top_genes = df.groupby('Gene_Short')['Conservation_Pct'].mean().sort_values(ascending=False).head(25).index
    else:
        top_genes = (df[df['Gene_Short'].isin(common_genes)]
                     .groupby('Gene_Short')['Conservation_Pct']
                     .mean().sort_values(ascending=False).head(25).index)

    # --- PLOT 1: DISTANCE HEATMAP (The 'Signature' Logic) ---
    dist_pivot = df.pivot_table(index='Gene_Short', columns='Host', values='Signature_Distance_bp').fillna(0)
    plt.figure(figsize=(12, 14))
    sns.heatmap(dist_pivot.loc[top_genes], annot=True, fmt=".0f", cmap="RdBu_r", center=0,
                cbar_kws={'label': 'Distance (bp)'})
    plt.title('Top Conserved Spatial Signatures: Genomic Distance (bp)')
    plt.tight_layout()
    plt.savefig('spatial_distance_heatmap.png')

    # --- PLOT 2: CONSERVATION HEATMAP (The 'Redundancy' Logic) ---
    cons_pivot = df.pivot_table(index='Gene_Short', columns='Host', values='Conservation_Pct').fillna(0)
    plt.figure(figsize=(12, 14))
    # Setting vmax=200 highlights the 120-200% range you mentioned
    sns.heatmap(cons_pivot.loc[top_genes], annot=True, fmt=".1f", cmap="YlGnBu", vmax=200,
                cbar_kws={'label': 'Prevalence / Redundancy (%)'})
    plt.title('Top Conserved Spatial Signatures: Conservation Prevalence (%)')
    plt.tight_layout()
    plt.savefig('conservation_heatmap.png')

    # --- PLOT 3: KEY ONE HEALTH MARKERS BAR CHART ---
    key_markers = ['vanA', 'vanC', 'MexB', 'tetB(P)', 'tet(S)', 'tet(T)', 'tet(O)', 'tetO']
    key_df = df[df['Gene_Short'].isin(key_markers)]

    if not key_df.empty:
        plt.figure(figsize=(14, 7))
        sns.barplot(data=key_df, x='Gene_Short', y='Conservation_Pct', hue='Host', palette='viridis')
        plt.axhline(100, color='red', linestyle='--', alpha=0.5, label='Single Copy Threshold')
        plt.title('Conservation & Redundancy of Key "One Health" Genetic Signatures')
        plt.ylabel('Conservation (%)')
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.xticks(rotation=45)
        plt.tight_layout()
        plt.savefig('key_genes_conservation.png')

    print("[✔] Visualization Complete. Files saved: distance_heatmap.png, conservation_heatmap.png, key_genes_conservation.png")
    plt.show()

if __name__ == "__main__":
    main()
