import pandas as pd
import numpy as np
import sys
import os

# Usage: python summarize_multiqc_stats.py <input_multiqc_general_stats.txt> <output_file>

def summarize_fastqc(input_file, output_file):
    # Get the project name from the output file name without the suffix
    project_name = os.path.splitext(os.path.basename(output_file))[0]

    # Read the input file into a pandas DataFrame
    df = pd.read_csv(input_file, sep='\t')
    
    # Check that the number of rows is even
    if len(df) % 2 != 0:
        raise ValueError("The input file must have an even number of rows for paired reads.")
    
    # Initialize variables
    paired_samples = len(df) // 2
    total_sequences = []
    low_read_samples = []

    # Process each pair of rows
    for i in range(0, len(df), 2):
        row1 = df.iloc[i]
        row2 = df.iloc[i + 1]

        # Check that the total sequences match for each pair
        if row1['FastQC_mqc-generalstats-fastqc-total_sequences'] != row2['FastQC_mqc-generalstats-fastqc-total_sequences']:
            raise ValueError(f"Total sequences do not match for sample pair: {row1['Sample']} and {row2['Sample']}")

        # Append total sequences value
        total_sequences.append(row1['FastQC_mqc-generalstats-fastqc-total_sequences'])
    
    # Calculate statistics
    total_sequences_array = np.array(total_sequences)
    mean_total_sequences = round(np.mean(total_sequences_array), 1)
    median_total_sequences = round(np.median(total_sequences_array), 1)
    min_total_sequences = round(np.min(total_sequences_array), 1)
    max_total_sequences = round(np.max(total_sequences_array), 1)
    quantiles = np.percentile(total_sequences_array, [0, 25, 50, 75, 100])
    twenty_fifth_quantile = round(quantiles[1], 1)

    # Identify and sort low read samples
    for i in range(0, len(df), 2):
        row1 = df.iloc[i]
        if row1['FastQC_mqc-generalstats-fastqc-total_sequences'] <= twenty_fifth_quantile:
            low_read_samples.append((row1['Sample'], row1['FastQC_mqc-generalstats-fastqc-total_sequences']))
    
    low_read_samples_sorted = sorted(low_read_samples, key=lambda x: x[1])
    low_read_samples_str = "|".join([f"{sample} ({reads})" for sample, reads in low_read_samples_sorted])
    quantiles_str = ", ".join([f"{q:.1f}" for q in quantiles])

    # Write the results to the output file
    with open(output_file, 'w') as f:
        f.write("ProjectName\tTotal Paired Samples\tMean Total Sequences\tMedian Total Sequences\tMin Total Sequences\tMax Total Sequences\tLow_read_samples\tQuantiles\n")
        f.write(f"{project_name}\t{paired_samples}\t{mean_total_sequences}\t{median_total_sequences}\t{min_total_sequences}\t{max_total_sequences}\t{low_read_samples_str}\t{quantiles_str}\n")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python summarize_multiqc_stats.py <input_multiqc_general_stats.txt> <output_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    summarize_fastqc(input_file, output_file)
