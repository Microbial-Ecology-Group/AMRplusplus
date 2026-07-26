#!/usr/bin/env python3
"""
Faster drop-in replacement for kraken2_long_to_wide.py
Uses multiprocessing to parse report files in parallel.
"""

import sys
import argparse
import re
import numpy as np
from multiprocessing import Pool, cpu_count
from collections import defaultdict

taxa_levels = {'D':0,'K':1,'P':2,'C':3,'O':4,'F':5,'G':6,'S':7,'U':8}

def parse_cmdline_params(args):
    p = argparse.ArgumentParser()
    p.add_argument('-i', '--input_files', nargs='+', required=True)
    p.add_argument('-o', '--output_file', required=True)
    p.add_argument('--merged', action='store_true')
    p.add_argument('--threads', type=int, default=cpu_count())
    return p.parse_args(args)


def parse_one_report(file_path):
    """Parse a single kraken report — runs in a worker process."""
    sample_id = file_path.split('/')[-1].replace('.kraken.report', '')
    counts = {}
    unclassified = [0, 0, 0.0]
    taxon_list = ['NA'] * 8
    previous_level = 0

    with open(file_path) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            if len(parts) < 6:
                continue

            node_count  = int(parts[2])
            node_level  = parts[3].strip()
            node_name   = parts[5].strip()

            if node_level == 'U':
                unclassified = [node_count, unclassified[1] + node_count, float(parts[0])]
                continue
            if node_level == 'R':
                unclassified[1] += int(parts[1])
                continue
            if len(node_level) > 1:
                if node_level[0] in ('U', 'R'):
                    continue
                parent = node_level[0]
            else:
                parent = node_level

            if parent not in taxa_levels:
                continue
            this_level = taxa_levels[parent]

            if len(node_level) == 1:
                taxon_list[this_level] = node_name
            if this_level < previous_level:
                taxon_list[this_level + 1:] = ['NA'] * (7 - this_level)
            previous_level = this_level

            if node_count == 0:
                continue

            taxon_str = '|'.join(taxon_list[:this_level + 1])
            counts[taxon_str] = counts.get(taxon_str, 0) + node_count

    return sample_id, counts, unclassified


def dict_to_matrix(counts_dict):
    samples = list(counts_dict.keys())
    all_taxa = list({t for c in counts_dict.values() for t in c})
    taxa_idx = {t: i for i, t in enumerate(all_taxa)}
    M = np.zeros((len(all_taxa), len(samples)), dtype=float)
    for j, sid in enumerate(samples):
        for taxon, cnt in counts_dict[sid].items():
            M[taxa_idx[taxon], j] = cnt
    return M, all_taxa, samples


def write_output(outfile, M, taxa, samples, unclassifieds):
    with open(outfile, 'w') as out:
        out.write(','.join(['taxa'] + samples) + '\n')
        for i, taxon in enumerate(taxa):
            out.write(f'"{taxon.replace(",","")}",')
            out.write(','.join(str(x) for x in M[i]) + '\n')

    with open(f"unclassifieds_{outfile}", 'w') as u:
        u.write('SampleID,NumberUnclassified,Total,PercentUnclassified\n')
        for sid, nums in unclassifieds.items():
            u.write(f"{sid},{','.join(str(x) for x in nums)}\n")


if __name__ == '__main__':
    opts = parse_cmdline_params(sys.argv[1:])

    print(f"Parsing {len(opts.input_files)} files with {opts.threads} threads...",
          file=sys.stderr)

    with Pool(processes=opts.threads) as pool:
        results = pool.map(parse_one_report, opts.input_files)

    counts_dict = {sid: c for sid, c, _ in results}
    uncls_dict  = {sid: u for sid, _, u in results}

    M, taxa, samples = dict_to_matrix(counts_dict)
    write_output(opts.output_file, M, taxa, samples, uncls_dict)
    print(f"Done → {opts.output_file}", file=sys.stderr)