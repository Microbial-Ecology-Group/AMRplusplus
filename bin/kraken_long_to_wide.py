#!/usr/bin/env python3

import argparse
import numpy as np
import sys

__author__ = "Steven Lakin"
__copyright__ = ""
__credits__ = ["Steven Lakin"]
__version__ = ""
__maintainer__ = "lakinsm"
__email__ = "lakinsm@colostate.edu"
__status__ = "Cows go moo."

taxa_level = {'D': 0, 'P': 1, 'C': 2, 'O': 3, 'F': 4, 'G': 5, 'S': 6}
taxa_level_names = {1: 'Domain', 2: 'Phylum', 3: 'Class', 4: 'Order',
                    5: 'Family', 6: 'Genus', 7: 'Species', 8: 'Unclassified'}

def parse_cmdline_params(cmdline_params):
    info = ""
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input_files', nargs='+', required=True,
                        help='Use globstar to pass a list of files, (Ex: *.tsv)')
    parser.add_argument('-o', '--output_directory', required=True,
                        help='Output directory for writing the kraken_analytic_matrix.csv file')
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
    ret = np.zeros((nrow, ncol), dtype=np.float)
    for j, (sample, tdict) in enumerate(D.items()):
        samples.append(sample)
        for i, taxon in enumerate(unique_nodes):
            if taxon in tdict:
                ret[i, j] = np.float(tdict[taxon])
    return ret, unique_nodes, samples

def kraken_load_analytic_data(file_name_list):
    ret = {}
    for file in file_name_list:
        sample_id = file.split('/')[-1].replace('.kraken.report', '')
        with open(file, 'r') as f:
            data = f.read().split('\n')
            assignment_list = [''] * 15
            taxon_list = ['NA'] * 7
            for entry in data:
                if not entry:
                    continue
                temp_name = entry.split('\t')[5]
                space_level = int((len(temp_name) - len(temp_name.lstrip(' '))) / 2) - 1
                if (space_level <= 0) and (''.join(entry.split()[5:]) not in ('Viruses', 'Bacteria', 'Archaea')):
                    continue
                if space_level < 0:
                    space_level = 0
                entry = entry.split()
                if entry[3] == 'U':
                    continue
                node_name = ' '.join(entry[5:])
                assignment_list[space_level] = node_name
                assignment_len = len(assignment_list) - assignment_list.count('')
                if (space_level + 1) < assignment_len:
                    assignment_list = assignment_list[:space_level + 1] + [''] * (14 - space_level)
                if entry[3] != '-':
                    node_level = taxa_level[entry[3]]
                    taxon_list[node_level] = node_name
                    taxon_name = '|'.join(taxon_list[:node_level+1])
                    taxon_len = len(taxon_list) - taxon_list.count('NA')
                    if (node_level + 1) < taxon_len:
                        taxon_list = taxon_list[:node_level + 1] + ['NA'] * (6 - node_level)
                temp_list = [x for x in taxon_list]
                if entry[3] == '-':
                    temp_list = [x for x in taxon_list]
                    while temp_list and temp_list[-1] == 'NA':
                        temp_list.pop()
                if (space_level + 1) < assignment_len:
                    iter_loc = space_level + 1
                    while True:
                        if iter_loc == 0:
                            break
                        try:
                            iter_loc = temp_list.index(assignment_list[iter_loc])
                            break
                        except ValueError:
                            iter_loc -= 1
                    temp_list = [x for x in taxon_list[:iter_loc + 1]]
                    taxon_list = taxon_list[:iter_loc + 1] + ['NA'] * (6 - iter_loc)
                    taxon_name = '|'.join(temp_list)
                if float(entry[2]) == 0.0:
                    continue
                try:
                    ret[sample_id][taxon_name] += float(entry[2])
                except KeyError:
                    try:
                        ret[sample_id].setdefault(taxon_name, float(entry[2]))
                    except KeyError:
                        ret.setdefault(sample_id, {taxon_name: float(entry[2])})
    return dict_to_matrix(ret)

def output_kraken_analytic_data(outdir, M, m_names, n_names):
    with open(outdir + '/kraken_analytic_matrix.csv', 'w') as out:
        out.write(','.join(n_names) + '\n')
        for i, row in enumerate(M):
            out.write('\"{}\",'.format(m_names[i].replace(',', '')))
            out.write(','.join([str(x) for x in row]) + '\n')

if __name__ == '__main__':
    opts = parse_cmdline_params(sys.argv[1:])
    K, m, n = kraken_load_analytic_data(opts.input_files)
    output_kraken_analytic_data(opts.output_directory, K, m, n)
