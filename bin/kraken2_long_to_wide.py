#!/usr/bin/env python3

import sys
import argparse
import numpy as np

__author__ = 'Steven Lakin'
__maintainer__ = 'lakinsm'
__email__ = 'Steven.Lakin@colostate.edu'


taxa_levels = {
    'D': 0,
    'K': 1,
    'P': 2,
    'C': 3,
    'O': 4,
    'F': 5,
    'G': 6,
    'S': 7,
    'U': 8
}

taxa_level_names = {
    0: 'Domain',
    1: 'Kingdom',
    2: 'Phylum',
    3: 'Class',
    4: 'Order',
    5: 'Family',
    6: 'Genus',
    7: 'Species',
    8: 'Unclassified'
}


def parse_cmdline_params(cmdline_params):
    info = ""
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input_files', nargs='+', required=True,
                        help='Use globstar to pass a list of files, (Ex: *.tsv)')
    parser.add_argument('-o', '--output_file', required=True,
                        help='Output file name for writing the kraken_analytic_matrix.csv file')
    return parser.parse_args(cmdline_params)


def dict_to_matrix(D):
    ncol = len(D.keys())
    unique_nodes = []
    samples = []
    for sample, tdict in D.items():
        for taxon in tdict.keys():
            if taxon not in unique_nodes:
                unique_nodes.append(taxon)
    nrow = len(unique_nodes)
    return_values = np.zeros((nrow, ncol), dtype=float)
    for j, (sample, tdict) in enumerate(D.items()):
        samples.append(sample)
        for i, taxon in enumerate(unique_nodes):
            if taxon in tdict:
                return_values[i, j] = float(tdict[taxon])
    return return_values, unique_nodes, samples


def kraken2_load_analytic_data(file_name_list):
    return_values = {}
    unclassifieds = {}  # { sample: [unclassified, total, percent] }
    for file in file_name_list:
        sample_id = file.split('/')[-1].replace('.kraken.report', '')
        unclassifieds.setdefault(sample_id, [0, 0, 0])
        with open(file, 'r') as f:
            data = f.read().split('\n')
            taxon_list = ['NA'] * 8
            previous_taxon_level = 0
            for line in data:
                if not line:
                    continue
                entries = line.split('\t')
                node_count = int(entries[2])
                node_level = entries[3]
                node_name = entries[5].strip()
                if node_level == 'U':
                    unclassifieds[sample_id][0] = node_count
                    unclassifieds[sample_id][1] += node_count
                    unclassifieds[sample_id][2] = float(entries[0])
                    continue
                elif node_level == 'R':
                    unclassifieds[sample_id][1] += int(entries[1])
                    continue
                if len(node_level) > 1:
                    if node_level[0] in ('U', 'R'):
                        continue
                    parent_node_level = node_level[0]
                else:
                    parent_node_level = node_level
                this_taxon_level = taxa_levels[parent_node_level]
                if len(node_level) == 1:
                    taxon_list[this_taxon_level] = node_name
                if this_taxon_level < previous_taxon_level:
                    taxon_list[this_taxon_level + 1:] = ['NA'] * (7 - this_taxon_level)
                previous_taxon_level = this_taxon_level
                if node_count == 0:
                    continue
                this_taxonomy_string = '|'.join(taxon_list[:this_taxon_level + 1])
                try:
                    return_values[sample_id][this_taxonomy_string] += node_count
                except KeyError:
                    try:
                        return_values[sample_id].setdefault(this_taxonomy_string, node_count)
                    except KeyError:
                        return_values.setdefault(sample_id, {this_taxonomy_string: node_count})
    return dict_to_matrix(return_values), unclassifieds


def output_kraken2_analytic_data(outfile, M, m_names, n_names, unclassifieds):
    with open(outfile, 'w') as out, \
            open(f"unclassifieds_{outfile}", 'w') as u_out:
        out.write(','.join(n_names) + '\n')
        for i, row in enumerate(M):
            out.write('\"{}\",'.format(
                m_names[i].replace(',', '')
            ))
            out.write(','.join([str(x) for x in row]) + '\n')
        u_out.write('SampleID,NumberUnclassified,Total,PercentUnclassified\n')
        for sample, numbers in unclassifieds.items():
            u_out.write('{},{}\n'.format(
                sample,
                ','.join([str(x) for x in numbers])
            ))


if __name__ == '__main__':
    opts = parse_cmdline_params(sys.argv[1:])
    kraken2_load_analytic_data(opts.input_files)
    (K, m, n), u = kraken2_load_analytic_data(opts.input_files)
    output_kraken2_analytic_data(opts.output_file, K, m, n, u)
