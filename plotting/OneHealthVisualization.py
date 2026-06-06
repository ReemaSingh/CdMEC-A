#!/usr/bin/env python3
# CdMEC-A: Clostridioides difficile Mobile Element Context Analyzer
# Copyright (C) 2025-2026 [Dr. Reema Singh]
# Licensed under the GNU General Public License v3.0

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import argparse
import os

def main():
    parser = argparse.ArgumentParser(description="CdMEC-A One Health Visualization Suite")
    parser.add_argument("input", nargs="?", default="one_health_spatial_signatures.csv",
                        help="Path to the input CSV signature file (default: one_health_spatial_signatures.csv)")
    parser.add_argument("--markers", nargs="+", 
                        default=['vanA', 'vanC', 'tet(S)', 'tet(T)', 'ranA', 'bcrA', 'poxtA'],
                        help="List of high-priority One Health sentinel markers to explicitly plot in bar charts.")
    args = parser.parse_args()

    if not os.path.exists(args.input):
        print(f"Error: {args.input} not found. Ensure you ran the updated SignatureDistance.py first.")
        return

    # Load data
    df = pd.read_csv(args.input)
    
    # Extract short gene name for cleaner labels (e.g., 'ranA')
    df['Gene_Short'] = df['Gene'].apply(lambda x: str(x).split('|')[-1] if '|' in str(x) else str(x))

    # --- SELECTION LOGIC FOR HEATMAPS ---
    common_genes = df.groupby('Gene_Short').filter(lambda x: x['Host'].nunique() == 3)['Gene_Short'].unique()
    
    if len(common_genes) == 0:
        print("Warning: No common genes found across all 3 hosts. Using top 25 by overall prevalence.")
        top_genes = df.groupby('Gene_Short')['Prevalence_Pct'].mean().sort_values(ascending=False).head(25).index
    else:
        top_genes = (df[df['Gene_Short'].isin(common_genes)]
                     .groupby('Gene_Short')['Avg_Copies_Per_Genome']
                     .mean().sort_values(ascending=False).head(25).index)
        
    # ---------------- FIGURE: PREVALENCE HEATMAP (%) ----------------
    prev_pivot = df.pivot_table(index='Gene_Short', columns='Host', values='Prevalence_Pct').fillna(0)
    prev_pivot_filtered = prev_pivot.reindex(top_genes, fill_value=0)
    
    fig, ax = plt.subplots(figsize=(10, 12))
    sns.heatmap(prev_pivot_filtered, annot=True, fmt=".1f", cmap="YlGnBu", vmax=100,
                cbar_kws={'label': 'Prevalence (%)'}, ax=ax)
    ax.set_title('Top Conserved Signatures - Prevalence (%)')
    
    # Clean axis label and italicize
    ax.set_ylabel('Gene')
    for label in ax.get_yticklabels():
        label.set_style('italic')
    
    plt.tight_layout()
    plt.savefig('prevalence.png')
    plt.close()

    # ---------------- FIGURE: REDUNDANCY HEATMAP (DOSAGE) ----------------
    redu_pivot = df.pivot_table(index='Gene_Short', columns='Host', values='Avg_Copies_Per_Genome').fillna(0)
    redu_pivot_filtered = redu_pivot.reindex(top_genes, fill_value=0)
    
    fig, ax = plt.subplots(figsize=(10, 12))
    sns.heatmap(redu_pivot_filtered, annot=True, fmt=".1f", cmap="YlOrRd",
                cbar_kws={'label': 'Avg Copies per Genome'}, ax=ax)
    ax.set_title('Top Conserved Signatures - Redundancy (Dosage)')
    
    # Clean axis label and italicize
    ax.set_ylabel('Gene')
    for label in ax.get_yticklabels():
        label.set_style('italic')
    
    plt.tight_layout()
    plt.savefig('redundancy.png')
    plt.close()

    # ---------------- FIGURE: GENOMIC DISTANCE HEATMAP (bp) ----------------
    dist_pivot = df.pivot_table(index='Gene_Short', columns='Host', values='Signature_Distance_bp').fillna(0)
    dist_pivot_filtered = dist_pivot.reindex(top_genes, fill_value=0)
    
    fig, ax = plt.subplots(figsize=(10, 12))
    sns.heatmap(dist_pivot_filtered, annot=True, fmt=".0f", cmap="RdBu_r", center=0,
                cbar_kws={'label': 'Distance to MGE (bp)'}, ax=ax)
    ax.set_title('Spatial Signatures - Genomic Distance (bp)')
    
    # Clean axis label and italicize
    ax.set_ylabel('Gene')
    for label in ax.get_yticklabels():
        label.set_style('italic')
    
    plt.tight_layout()
    plt.savefig('distance.png')
    plt.close()
    
    # ---------------- FIGURE: KEY ONE HEALTH MARKER BARS ----------------
    key_markers_lower = [str(m).lower() for m in args.markers]
    key_df = df[df['Gene_Short'].str.lower().isin(key_markers_lower)]
    
    if key_df.empty:
        print(f"\n[Notice]: None of the specified priority markers {args.markers} were found.")
        print("Automatically defaulting bar plots to the top priority markers by prevalence in this dataset.")
        top_bar_genes = df.groupby('Gene_Short')['Prevalence_Pct'].mean().sort_values(ascending=False).head(7).index
        key_df = df[df['Gene_Short'].isin(top_bar_genes)]

    if not key_df.empty:
        # Prevalence Comparison
        fig, ax = plt.subplots(figsize=(12, 6))
        sns.barplot(data=key_df, x='Gene_Short', y='Prevalence_Pct', hue='Host', palette='viridis', ax=ax)
        plt.axhline(100, color='red', linestyle='--', alpha=0.6, label='100% Fixation')
        ax.set_title('Prevalence of Priority One Health Markers')
        ax.set_ylabel('Prevalence (%)')
        ax.set_ylim(0, 120)
        
        # Clean x-axis label from "Gene_Short" to "Gene"
        ax.set_xlabel('Gene')
        for label in ax.get_xticklabels():
            label.set_style('italic')
        
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig('key_prevalence.png')
        plt.close()

        # Redundancy Comparison
        fig, ax = plt.subplots(figsize=(12, 6))
        sns.barplot(data=key_df, x='Gene_Short', y='Avg_Copies_Per_Genome', hue='Host', palette='magma', ax=ax)
        plt.axhline(1, color='black', linestyle='-', alpha=0.3, label='Single Copy Baseline')
        ax.set_title('Genomic Redundancy of Priority Markers')
        ax.set_ylabel('Avg Copies Per Genome')
        
        # Clean x-axis label from "Gene_Short" to "Gene"
        ax.set_xlabel('Gene')
        for label in ax.get_xticklabels():
            label.set_style('italic')
        
        plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        plt.savefig('key_redundancy.png')
        plt.close()

    print("\nVisualization Suite Complete!")
    print("Generated: prevalence.png, redundancy.png, distance.png, key_prevalence.png, key_redundancy.png")

if __name__ == "__main__":
    main()
