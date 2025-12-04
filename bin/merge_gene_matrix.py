#!/usr/bin/env python3
import argparse
import glob
import os
import sys
from collections import defaultdict
from typing import Dict, Set

## Example usage
#python3 merge_gene_matrix.py \
#    -i "counts/*gene_summary.tsv" \
#    -o "Deduped_AMR_analytic_matrix.csv"

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Merge multiple gene summary TSV files into a single matrix. "
            "Rows are genes, columns are samples, values are counts."
        )
    )
    parser.add_argument(
        "-i", "--input-pattern",
        required=True,
        help="Glob pattern to match gene summary files (e.g., 'counts/*gene_summary.tsv')"
    )
    parser.add_argument(
        "-o", "--output",
        required=True,
        help="Output file path (e.g., 'Deduped_AMR_analytic_matrix.csv')"
    )
    parser.add_argument(
        "--fill-value",
        type=int,
        default=0,
        help="Value to use when a gene is not found in a sample (default: 0)"
    )
    return parser.parse_args()


def parse_gene_summary(filepath: str) -> tuple:
    """Parse a gene summary TSV file.
    
    Returns:
        tuple: (sample_id, dict of gene -> count)
    """
    gene_counts = {}
    sample_id = None
    
    with open(filepath, "r") as f:
        # Parse header to get sample ID
        header = f.readline().strip().split("\t")
        if len(header) != 2:
            sys.exit(f"[ERROR] Expected 2 columns in {filepath}, got {len(header)}")
        
        sample_id = header[1]  # Second column header is the sample ID
        
        # Parse data rows
        for line in f:
            line = line.strip()
            if not line:
                continue
            
            parts = line.split("\t")
            if len(parts) != 2:
                continue
            
            gene, count = parts
            gene_counts[gene] = int(count)
    
    return sample_id, gene_counts


def merge_gene_summaries(
    file_pattern: str,
    fill_value: int = 0
) -> tuple:
    """Merge multiple gene summary files into a matrix.
    
    Returns:
        tuple: (list of sample_ids, set of all genes, dict of sample_id -> {gene -> count})
    """
    files = sorted(glob.glob(file_pattern))
    
    if not files:
        sys.exit(f"[ERROR] No files found matching pattern: {file_pattern}")
    
    print(f"[INFO] Found {len(files)} files matching pattern")
    
    all_genes: Set[str] = set()
    sample_data: Dict[str, Dict[str, int]] = {}
    sample_ids = []
    
    for filepath in files:
        print(f"[INFO] Reading: {filepath}")
        sample_id, gene_counts = parse_gene_summary(filepath)
        
        if sample_id in sample_data:
            print(f"[WARNING] Duplicate sample ID '{sample_id}' found in {filepath}, skipping")
            continue
        
        sample_ids.append(sample_id)
        sample_data[sample_id] = gene_counts
        all_genes.update(gene_counts.keys())
    
    print(f"[INFO] Total samples: {len(sample_ids)}")
    print(f"[INFO] Total unique genes: {len(all_genes)}")
    
    return sample_ids, all_genes, sample_data


def write_matrix(
    sample_ids: list,
    all_genes: Set[str],
    sample_data: Dict[str, Dict[str, int]],
    output_path: str,
    fill_value: int = 0
) -> None:
    """Write the merged matrix to a CSV file."""
    
    sorted_genes = sorted(all_genes)
    
    with open(output_path, "w") as out:
        # Write header
        header = ["gene_accession"] + sample_ids
        out.write(",".join(header) + "\n")
        
        # Write data rows
        for gene in sorted_genes:
            row = [gene]
            for sample_id in sample_ids:
                count = sample_data[sample_id].get(gene, fill_value)
                row.append(str(count))
            out.write(",".join(row) + "\n")


def main():
    args = parse_args()
    
    sample_ids, all_genes, sample_data = merge_gene_summaries(
        args.input_pattern,
        args.fill_value
    )
    
    write_matrix(
        sample_ids,
        all_genes,
        sample_data,
        args.output,
        args.fill_value
    )
    
    print(f"[INFO] Wrote matrix to: {args.output}")
    print(f"[INFO] Matrix dimensions: {len(all_genes)} genes x {len(sample_ids)} samples")


if __name__ == "__main__":
    main()