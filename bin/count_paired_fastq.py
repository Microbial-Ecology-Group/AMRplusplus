#!/usr/bin/env python3

import argparse
import gzip
import os

## Example
# python3 count_paired_reads.py -f forward.fastq.gz -r reverse.fastq.gz -s sample_name -c NonHost


def parse_cmdline_params(cmdline_params):
    parser = argparse.ArgumentParser(description="Count the number of sequences in paired FASTQ files.")
    parser.add_argument('-f', '--forward_fastq', required=True,
                        help='Path to the forward FASTQ file (can be gzipped).')
    parser.add_argument('-r', '--reverse_fastq', required=True,
                        help='Path to the reverse FASTQ file (can be gzipped).')
    parser.add_argument('-s', '--sample_name', required=True,
                        help='Sample name to include in the output file.')
    parser.add_argument('-c', '--column_name', required=True,
                        help='Column name for the output file.')
    return parser.parse_args(cmdline_params)

def count_reads(fastq_file):
    open_func = gzip.open if fastq_file.endswith('.gz') else open
    with open_func(fastq_file, 'rt') as f:
        return sum(1 for _ in f) // 4

def write_output(output_file, sample_name, column_name, read_count):
    with open(output_file, 'w') as o:
        o.write(f'Sample\t{column_name}\n')
        o.write(f'{sample_name}\t{read_count}\n')

if __name__ == "__main__":
    opts = parse_cmdline_params(sys.argv[1:])
    
    forward_count = count_reads(opts.forward_fastq)
    reverse_count = count_reads(opts.reverse_fastq)
    
    # Check that read numbers match
    if forward_count != reverse_count:
        print(f"Error: The number of reads in the forward file ({forward_count}) does not match the number of reads in the reverse file ({reverse_count}).", file=sys.stderr)
        sys.exit(1)

    paired_reads = forward_count  

    output_file = f"{opts.sample_name}_{opts.column_name}_pe_read_counts.txt"
    write_output(output_file, opts.sample_name, opts.column_name, paired_reads)
