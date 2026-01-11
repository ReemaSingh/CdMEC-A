## Core Analysis (bin/)

```bash cdmec_analyzer.py```: The engine of the pipeline. It performs multithreaded BLAST searches (blastx for ARGs, blastn for MGEs) and identifies pairs located within a 10kb window on the same contig.

```bash cdmec_stats_generator.py```: A post-processing tool that filters the raw results to identify "High-Risk" associations (defined as distance < 1kb or embedded). It outputs summary tables of the most mobile ARGs and common MGE carriers.

```bash SignatureDistance.py```: The script calculates the most frequent physical distance (Spatial Signature) and occurrence rate (Redundancy) between resistance genes and mobile genetic elements to track stable mobilization units across different host environments.
