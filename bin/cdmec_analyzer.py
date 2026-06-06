# CdMEC-A: Clostridioides difficile Mobile Element Context Analyzer
# Copyright (C) 2025-2026 [Dr. Reema Singh]
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

# Absolute paths are evaluated here to ensure parallel background workers do not lose directory scope
ARG_DB_PATH = os.path.abspath("card_protein_homolog_db")
MGE_DB_PATH = os.path.abspath("combined_C_Diff_mge_nucl_db")
CONTEXT_THRESHOLD = 10000

# Strict thresholds to eliminate false positives and control alignment inflation
MIN_ARG_IDENTITY = 80.0   # Mandatory 80% amino acid identity for true ARGs
MIN_ARG_COVERAGE = 80.0   # Mandatory 80% coverage of the target reference protein
MAX_EVALUE = "1e-10"      # Highly stringent baseline E-value threshold

def parse_arguments():
    parser = argparse.ArgumentParser(
        description="CdMEC-A: Rigorous Context Analyzer with Validated Homology Filters.",
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
    """Executes BLAST with exact hit evaluation formats."""
    tool = BLAST_TOOL_PROT if hit_type == "ARG" else BLAST_TOOL_NUCL
    
    if not query_fasta or not os.path.exists(query_fasta):
        sys.stderr.write(f"Worker Error: Target query file does not exist: '{query_fasta}'\n")
        return []

    # Outfmt 6 including 'slen' (Subject Sequence Length) to calculate true target coverage
    blast_cmd = [
        tool, 
        "-query", os.path.abspath(query_fasta),
        "-db", os.path.abspath(db_prefix),
        "-num_threads", str(blast_threads).strip(),
        "-outfmt", "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore slen",
        "-evalue", MAX_EVALUE 
    ]
    try:
        # capture_output=True isolates standard error streams for clear debugging
        result = subprocess.run(blast_cmd, capture_output=True, check=True) 
        stdout_str = result.stdout.decode('utf-8') 
        return [line for line in stdout_str.split('\n') if line.strip()]
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"\n--- BLAST CRASH DETECTED FOR SAMPLE: {os.path.basename(query_fasta)} ---\n")
        sys.stderr.write(f"Command Attempted: {' '.join(blast_cmd)}\n")
        sys.stderr.write(f"BLAST Error Output: {e.stderr.decode('utf-8')}\n")
        sys.stderr.write("-----------------------------------------------------\n")
        return []
    except Exception as e:
        sys.stderr.write(f"General processing error in {hit_type} module: {str(e)}\n")
        return []

def parse_and_filter_hits(raw_output_lines, hit_type):
    """
    Parses alignment outputs, applies stringent stringency filters to ARGs using slen,
    and merges highly overlapping fragments to stop artificial copy-number duplication.
    """
    raw_hits = []
    for line in raw_output_lines:
        fields = line.strip().split() 
        if not fields or len(fields) < 11: 
            continue 
        try:
            pident = float(fields[2])
            length = float(fields[3])
            q_start, q_end = int(fields[4]), int(fields[5])
            slen = float(fields[10]) # True length of the target reference protein/sequence
            
            if hit_type == "ARG":
                # Fixed: Uses slen explicitly and safely to avoid NameError
                coverage = (length / slen) * 100 if slen > 0 else 0
                if pident < MIN_ARG_IDENTITY or coverage < MIN_ARG_COVERAGE:
                    continue
            else:
                # Fixed: Uses slen for the MGE nucleotide tracking safely
                coverage = (length / slen) * 100 if slen > 0 else 0
                if pident < 70.0 or coverage < 30.0:
                    continue

            raw_hits.append({
                "Contig_ID": fields[0], 
                "Hit_Name": fields[1],
                "Start": min(q_start, q_end), 
                "End": max(q_start, q_end), 
                "Hit_Type": hit_type
            })
        except ValueError: 
            continue

    if not raw_hits:
        return []
        
    # De-duplication Engine: Merge fragment overlays mapping to identical loci
    raw_hits.sort(key=lambda x: (x["Contig_ID"], x["Hit_Name"], x["Start"]))
    
    filtered_hits = []
    current_hit = raw_hits[0]
    
    for next_hit in raw_hits[1:]:
        if (next_hit["Contig_ID"] == current_hit["Contig_ID"] and 
            next_hit["Hit_Name"] == current_hit["Hit_Name"] and 
            next_hit["Start"] <= current_hit["End"]):
            current_hit["End"] = max(current_hit["End"], next_hit["End"])
        else:
            filtered_hits.append(current_hit)
            current_hit = next_hit
    filtered_hits.append(current_hit)
    
    return filtered_hits

def calculate_distance(arg_start, arg_end, mge_start, mge_end):
    if max(arg_start, mge_start) < min(arg_end, mge_end): 
        return 0, "Overlapping"
    elif mge_end < arg_start: 
        return -(arg_start - mge_end), "Upstream"
    elif arg_end < mge_start: 
        return (mge_start - arg_end), "Downstream"
    return 0, "Internal Overlap"

def analyze_context(arg_hits, mge_hits):
    results = []
    mge_by_contig = {}
    for hit in mge_hits:
        mge_by_contig.setdefault(hit['Contig_ID'], []).append(hit)
        
    for arg in arg_hits:
        contig_id = arg['Contig_ID']
        nearest_mge, min_dist, final_signed, final_status = None, float('inf'), None, "Chromosomal (Unlinked)"
        
        if contig_id in mge_by_contig:
            for mge in mge_by_contig[contig_id]:
                s_dist, p_status = calculate_distance(arg['Start'], arg['End'], mge['Start'], mge['End'])
                if abs(s_dist) <= CONTEXT_THRESHOLD and abs(s_dist) < min_dist:
                    min_dist, nearest_mge, final_signed = abs(s_dist), mge, s_dist
                    if p_status in ["Overlapping", "Internal Overlap"]: 
                        final_status = "Embedded within MGE"
                    elif any(x in nearest_mge['Hit_Name'] for x in ["Transposase", "Integrase"]):
                        final_status = f"Transposon-Associated ({p_status})"
                    elif "rep" in nearest_mge['Hit_Name'] or "plasmid" in nearest_mge['Hit_Name'].lower():
                        final_status = f"Plasmid-Associated ({p_status})"
                    else: 
                        final_status = f"MGE-Associated ({p_status})"
                        
        results.append({
            "ARG_Name": arg['Hit_Name'], "Contig_ID": contig_id,
            "ARG_Start": arg['Start'], "ARG_End": arg['End'],
            "MGE_Association": f"{nearest_mge['Hit_Name']}:{nearest_mge['Start']}-{nearest_mge['End']}" if nearest_mge else "None",
            "Proximity_bp": final_signed if nearest_mge else "N/A", 
            "Inferred_Status": final_status
        })
    return results

def write_output(output_dir, sample_id, results):
    os.makedirs(output_dir, exist_ok=True)
    json_path = os.path.join(output_dir, f"{sample_id}_cdmec.json")
    with open(json_path, 'w') as f: 
        json.dump({"Sample_ID": sample_id, "ARG_Hits": results}, f, indent=4)
        
    tsv_path = os.path.join(output_dir, f"{sample_id}_cdmec_summary.tsv")
    with open(tsv_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=["Sample_ID", "Contig_ID", "ARG_Name", "ARG_Start", "ARG_End", "MGE_Association", "Proximity_bp", "Inferred_Status"], delimiter='\t')
        writer.writeheader()
        for hit in results: 
            writer.writerow({**hit, "Sample_ID": sample_id})

def process_sample(input_fasta_path, output_dir, blast_threads):
    sample_id = os.path.basename(input_fasta_path).split('.')[0].strip()
    if not sample_id:
        return "Skipped: Malformed empty path string detected"
        
    raw_arg = run_homology_search(input_fasta_path, ARG_DB_PATH, "ARG", blast_threads)
    raw_mge = run_homology_search(input_fasta_path, MGE_DB_PATH, "MGE", blast_threads)
    
    parsed_args = parse_and_filter_hits(raw_arg, "ARG")
    parsed_mges = parse_and_filter_hits(raw_mge, "MGE")
    
    context_results = analyze_context(parsed_args, parsed_mges)
    
    if context_results:
        write_output(output_dir, sample_id, context_results)
        mge_linked = sum(1 for x in context_results if x["Inferred_Status"] != "Chromosomal (Unlinked)")
        return f"Done: {sample_id} ({len(context_results)} total ARGs mapped, {mge_linked} linked to MGEs)"
    return f"Done: {sample_id} (No hits)"

if __name__ == "__main__":
    args = parse_arguments()
    
    # Safely extract, clean, and deduplicate physical fasta coordinates before execution
    raw_files = [f for ext in ["*.fa", "*.fasta", "*.fna"] for f in glob.glob(os.path.join(args.input_dir, ext))]
    fasta_files = sorted(list(set([os.path.abspath(f) for f in raw_files if os.path.isfile(f) and os.path.getsize(f) > 0])))
    
    if not fasta_files:
        print("Error: No valid, non-empty FASTA files discovered in the input path.")
        sys.exit(1)

    print(f"Successfully validated {len(fasta_files)} sample files for analysis.")
    print(f"Parallel Profile: Using {args.workers} workers, allocating {args.blast_threads} thread(s) per BLAST instance.")
    
    with ProcessPoolExecutor(max_workers=args.workers) as executor:
        futures = {executor.submit(process_sample, f, args.output_dir, args.blast_threads): f for f in fasta_files}
        for future in as_completed(futures):
            print(future.result())
