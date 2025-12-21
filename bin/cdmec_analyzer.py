# CdMEC-A: Contextual mDNA Mobile Element Classifier - Analyzer
# Copyright (C) 2025 [Dr. Reema Singh]
# Licensed under the GNU General Public License v3.0

import argparse
import subprocess
import json
import csv
import sys
import os
import glob
from concurrent.futures import ProcessPoolExecutor, as_completed

# --- Configuration & Global Variables ---
BLAST_TOOL_NUCL = "blastn"
BLAST_TOOL_PROT = "blastx"

ARG_DB_PATH = "card_protein_homolog_db"
MGE_DB_PATH = "combined_C_Diff_mge_nucl_db"
CONTEXT_THRESHOLD = 10000

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="CdMEC-A: Hybrid Multithreaded Context Analyzer.",
        formatter_class=argparse.RawTextHelpFormatter 
    )
    parser.add_argument("-i", "--input_dir", required=True, help="Directory containing FASTA files.")
    parser.add_argument("-o", "--output_dir", default="./cdmec_analysis_reports", help="Output directory.")
    parser.add_argument("-t", "--workers", type=int, default=2, 
                        help="Number of samples to process in parallel (Python workers).")
    parser.add_argument("-bt", "--blast_threads", type=str, default="4", 
                        help="Number of threads per BLAST command (BLAST -num_threads).")
    return parser.parse_args()

def run_homology_search(query_fasta, db_prefix, hit_type, blast_threads):
    """Executes BLAST with the -num_threads parameter."""
    tool = BLAST_TOOL_PROT if hit_type == "ARG" else BLAST_TOOL_NUCL
        
    blast_cmd = [
        tool, 
        "-query", query_fasta,
        "-db", db_prefix,
        "-num_threads", blast_threads,  # Internal BLAST multithreading
        "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore",
        "-evalue", "1e-5" 
    ]
    try:
        result = subprocess.run(blast_cmd, capture_output=True, check=True, timeout=600) 
        stdout_str = result.stdout.decode('utf-8') 
        return [line for line in stdout_str.split('\n') if line.strip()]
    except Exception as e:
        sys.stderr.write(f"Error in {hit_type} search for {os.path.basename(query_fasta)}: {str(e)}\n")
        return []

# ... [parse_hits, calculate_distance, analyze_context, write_output functions remain the same] ...

def parse_hits(raw_output_lines, hit_type):
    parsed_hits = []
    for line in raw_output_lines:
        fields = line.strip().split() 
        if not fields or len(fields) < 10: continue 
        try:
            q_start, q_end = int(fields[4]), int(fields[5])
            parsed_hits.append({
                "Contig_ID": fields[0], "Hit_Name": fields[1],
                "Start": min(q_start, q_end), "End": max(q_start, q_end), "Hit_Type": hit_type
            })
        except ValueError: continue
    return parsed_hits 

def calculate_distance(arg_start, arg_end, mge_start, mge_end):
    if max(arg_start, mge_start) < min(arg_end, mge_end): return 0, "Overlapping"
    elif mge_end < arg_start: return -(arg_start - mge_end), "Upstream"
    elif arg_end < mge_start: return (mge_start - arg_end), "Downstream"
    return 0, "Internal Overlap"

def analyze_context(arg_hits, mge_hits):
    results = []
    mge_by_contig = {}
    for hit in mge_hits:
        mge_by_contig.setdefault(hit['Contig_ID'], []).append(hit)
    for arg in arg_hits:
        contig_id = arg['Contig_ID']
        if contig_id in mge_by_contig:
            nearest_mge, min_dist, final_signed, final_status = None, float('inf'), None, None
            for mge in mge_by_contig[contig_id]:
                s_dist, p_status = calculate_distance(arg['Start'], arg['End'], mge['Start'], mge['End'])
                if abs(s_dist) <= CONTEXT_THRESHOLD and abs(s_dist) < min_dist:
                    min_dist, nearest_mge, final_signed = abs(s_dist), mge, s_dist
                    if p_status in ["Overlapping", "Internal Overlap"]: final_status = "Embedded within MGE"
                    elif any(x in nearest_mge['Hit_Name'] for x in ["Transposase", "Integrase"]):
                        final_status = f"Transposon-Associated ({p_status})"
                    elif "rep" in nearest_mge['Hit_Name'] or "plasmid" in nearest_mge['Hit_Name'].lower():
                        final_status = f"Plasmid-Associated ({p_status})"
                    else: final_status = f"MGE-Associated ({p_status})"
            if nearest_mge:
                results.append({
                    "ARG_Name": arg['Hit_Name'], "Contig_ID": contig_id,
                    "ARG_Start": arg['Start'], "ARG_End": arg['End'],
                    "MGE_Association": f"{nearest_mge['Hit_Name']}:{nearest_mge['Start']}-{nearest_mge['End']}",
                    "Proximity_bp": final_signed, "Inferred_Status": final_status
                })
    return results

def write_output(output_dir, sample_id, results):
    os.makedirs(output_dir, exist_ok=True)
    json_path = os.path.join(output_dir, f"{sample_id}_cdmec.json")
    with open(json_path, 'w') as f: json.dump({"Sample_ID": sample_id, "ARG_Hits": results}, f, indent=4)
    tsv_path = os.path.join(output_dir, f"{sample_id}_cdmec_summary.tsv")
    with open(tsv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=["Sample_ID", "Contig_ID", "ARG_Name", "ARG_Start", "ARG_End", "MGE_Association", "Proximity_bp", "Inferred_Status"], delimiter='\t')
        writer.writeheader()
        for hit in results: writer.writerow({**hit, "Sample_ID": sample_id})

def process_sample(input_fasta_path, output_dir, blast_threads):
    sample_id = os.path.basename(input_fasta_path).split('.')[0]
    raw_arg = run_homology_search(input_fasta_path, ARG_DB_PATH, "ARG", blast_threads)
    raw_mge = run_homology_search(input_fasta_path, MGE_DB_PATH, "MGE", blast_threads)
    context_results = analyze_context(parse_hits(raw_arg, "ARG"), parse_hits(raw_mge, "MGE"))
    if context_results:
        write_output(output_dir, sample_id, context_results)
        return f"Done: {sample_id}"
    return f"Done: {sample_id} (No hits)"

if __name__ == "__main__":
    args = parse_arguments()
    fasta_files = sorted(list(set([f for ext in ["*.fa", "*.fasta", "*.fna"] for f in glob.glob(os.path.join(args.input_dir, ext))])))
    
    if not fasta_files:
        print("No FASTA files found."); sys.exit(1)

    print(f"Processing {len(fasta_files)} files.")
    print(f"Settings: {args.workers} samples at a time, each using {args.blast_threads} BLAST threads.")
    
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(process_sample, f, args.output_dir, args.blast_threads): f for f in fasta_files}
        for future in as_completed(futures):
            print(future.result())
