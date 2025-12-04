#!/usr/bin/env python3
"""
plot_histogram.py

Plot a histogram of values from a specified column in a coverage stats TSV file.

Usage:
    python3 plot_histogram.py -i input_file.tsv --column num_all_alignments
    python3 plot_histogram.py -i input_file.tsv --column num_all_alignments --include-zeros
    python3 plot_histogram.py -i input_file.tsv --column num_all_alignments --bins 50
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Plot histogram of a specified column from coverage stats TSV file."
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input TSV file (e.g., F1F01c_combined_coverage_stats.tsv)"
    )
    parser.add_argument(
        "--column",
        required=True,
        help="Column name to plot (e.g., num_all_alignments, num_primary_alignments)"
    )
    parser.add_argument(
        "--include-zeros",
        action="store_true",
        help="Include zero values in the histogram (default: exclude zeros)"
    )
    parser.add_argument(
        "--bins",
        type=int,
        default=30,
        help="Number of bins for histogram (default: 30)"
    )
    parser.add_argument(
        "-o", "--output",
        help="Output file path (default: auto-generated from input filename and column)"
    )
    return parser.parse_args()


def main():
    args = parse_args()
    
    # Read the TSV file
    df = pd.read_csv(args.input, sep="\t")
    
    # Check if column exists
    if args.column not in df.columns:
        print(f"[ERROR] Column '{args.column}' not found in file.")
        print(f"Available columns: {', '.join(df.columns)}")
        return
    
    # Get values from specified column
    values = df[args.column]
    
    # Filter zeros if not including them
    if args.include_zeros:
        plot_values = values
        zero_label = "including zeros"
    else:
        plot_values = values[values > 0]
        zero_label = "excluding zeros"
    
    # Extract sample ID from filename
    input_path = Path(args.input)
    sample_id = input_path.stem.replace("_coverage_stats", "")
    
    # Generate output filename if not specified
    if args.output:
        output_path = args.output
    else:
        zeros_suffix = "_with_zeros" if args.include_zeros else "_no_zeros"
        output_path = f"{sample_id}_{args.column}{zeros_suffix}_histogram.png"
    
    # Plot histogram
    fig, ax = plt.subplots(figsize=(10, 6))
    
    ax.hist(plot_values, bins=args.bins, edgecolor="black", alpha=0.7)
    ax.set_xlabel(args.column, fontsize=12)
    ax.set_ylabel("Frequency", fontsize=12)
    ax.set_title(f"{sample_id}\n{args.column} ({zero_label})", fontsize=14)
    
    # Add summary statistics
    stats_text = (
        f"N = {len(plot_values):,}\n"
        f"Mean = {plot_values.mean():.2f}\n"
        f"Median = {plot_values.median():.2f}\n"
        f"Max = {plot_values.max():,}"
    )
    if not args.include_zeros:
        stats_text += f"\nZeros excluded = {(values == 0).sum():,}"
    
    ax.text(
        0.95, 0.95, stats_text,
        transform=ax.transAxes,
        fontsize=10,
        verticalalignment="top",
        horizontalalignment="right",
        bbox=dict(boxstyle="round", facecolor="wheat", alpha=0.5)
    )
    
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    
    print(f"[INFO] Saved histogram to: {output_path}")
    print(f"[INFO] Total reads: {len(values):,}")
    print(f"[INFO] Reads with zero: {(values == 0).sum():,}")
    print(f"[INFO] Reads plotted: {len(plot_values):,}")


if __name__ == "__main__":
    main()
