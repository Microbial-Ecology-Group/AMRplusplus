#!/usr/bin/env python3

__author__ = "Steven Lakin"
__copyright__ = ""
__credits__ = ["Steven Lakin"]
__version__ = ""
__maintainer__ = "lakinsm"
__email__ = "lakinsm@colostate.edu"
__status__ = "Cows go moo."

import argparse
import sys

amr_level_names = {0: 'Class', 1: 'Mechanism', 2: 'Group'}

def parse_cmdline_params(cmdline_params):
    info = ""
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input_files', nargs='+', required=True,
                        help='Use globstar to pass a list of files, (Ex: *.tsv)')
    parser.add_argument('-o', '--output_directory', required=True,
                        help='Output directory for writing the AMR_analytic_matrix.csv file')
    return parser.parse_args(cmdline_params)

def amr_load_data(file_name_list):
    samples = {}
    labels = set()
    for file in file_name_list:
        with open(file, 'r') as f:
            data = f.read().split('\n')[1:]
            for entry in data:
                if not entry:
                    continue
                entry = entry.split('\t')
                sample = entry[0].split('.')[0]
                count = float(entry[2])
                gene_name = entry[1]
                try:
                    samples[sample][gene_name] = count
                except KeyError:
                    try:
                        samples[sample].setdefault(gene_name, count)
                    except KeyError:
                        samples.setdefault(sample, {gene_name: count})
                labels.add(gene_name)
    return samples, labels

def output_amr_analytic_data(outdir, S, L):
    with open(outdir + '/AMR_analytic_matrix.csv', 'w') as amr:
        local_sample_names = []
        for sample, dat in S.items():
            local_sample_names.append(sample)
        amr.write(','.join(local_sample_names) + '\n')
        for label in L:
            local_counts = []
            amr.write(label + ',')
            for local_sample in local_sample_names:
                if label in S[local_sample]:
                    local_counts.append(str(S[local_sample][label]))
                else:
                    local_counts.append(str(0))
            amr.write(','.join(local_counts) + '\n')

if __name__ == '__main__':
    opts = parse_cmdline_params(sys.argv[1:])
    S, L = amr_load_data(opts.input_files)
    output_amr_analytic_data(opts.output_directory, S, L)
