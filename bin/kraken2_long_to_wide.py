#!/usr/bin/env python3

import sys
import argparse
import numpy as np
import re

__authors__ = 'Steven Lakin & Enrique Doster'
__maintainer__ = 'Enrique Doster'
__email__ = 'enriquedoster@tamu.edu'


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
    parser.add_argument('--merged', action='store_true',
                        help='If set, combine paired samples ending with merged/unmerged into one *_combined column')
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

def _split_base_and_tag(sample_id):
    """
    Return (base, tag) where tag is 'merged' or 'unmerged' if present at the end,
    preceded optionally by '.', '_' or '-'. Otherwise tag is None.
    Examples:
      'S1.merged'   -> ('S1', 'merged')
      'S1-unmerged' -> ('S1', 'unmerged')
      'S1'          -> ('S1', None)
    """
    m = re.search(r'^(?P<base>.*?)(?:[._-]?)(?P<tag>merged|unmerged)$', sample_id)
    if m:
        return m.group('base'), m.group('tag')
    return sample_id, None


def combine_merged_unmerged(counts_dict, unclassifieds_dict):
    """
    counts_dict: { sample_id: {taxon: count, ...}, ... }
    unclassifieds_dict: { sample_id: [unclassified, total, percent], ... }

    Returns new dicts where paired (base.merged, base.unmerged) are replaced by base_combined.
    """
    from collections import defaultdict

    # Group sample IDs by their 'base'
    groups = defaultdict(dict)  # base -> {'merged': id, 'unmerged': id}
    for sid in counts_dict.keys():
        base, tag = _split_base_and_tag(sid)
        if tag in ('merged', 'unmerged'):
            groups[base][tag] = sid

    # Start with a copy; we’ll rebuild to avoid in-place confusion
    out_counts = {}
    out_uncls  = {}

    # First, mark all paired sids so we can skip copying them verbatim
    paired = set()
    for base, d in groups.items():
        if 'merged' in d and 'unmerged' in d:
            paired.add(d['merged'])
            paired.add(d['unmerged'])

    # Copy over everything that is not part of a complete pair
    for sid, tdict in counts_dict.items():
        if sid not in paired:
            out_counts[sid] = dict(tdict)
            if sid in unclassifieds_dict:
                out_uncls[sid] = list(unclassifieds_dict[sid])

    # Now add combined entries for complete pairs
    for base, d in groups.items():
        if 'merged' in d and 'unmerged' in d:
            s_merged   = d['merged']
            s_unmerged = d['unmerged']
            combined_sid = f"{base}_combined"

            # Sum taxon counts
            combo = {}
            for k, v in counts_dict[s_merged].items():
                combo[k] = combo.get(k, 0.0) + float(v)
            for k, v in counts_dict[s_unmerged].items():
                combo[k] = combo.get(k, 0.0) + float(v)
            out_counts[combined_sid] = combo

            # Sum unclassifieds; recompute percent
            u1 = unclassifieds_dict.get(s_merged,  [0, 0, 0.0])
            u2 = unclassifieds_dict.get(s_unmerged, [0, 0, 0.0])
            u = [0, 0, 0.0]
            u[0] = (u1[0] if isinstance(u1[0], int) else int(u1[0])) + (u2[0] if isinstance(u2[0], int) else int(u2[0]))
            u[1] = (u1[1] if isinstance(u1[1], int) else int(u1[1])) + (u2[1] if isinstance(u2[1], int) else int(u2[1]))
            u[2] = (100.0 * u[0] / u[1]) if u[1] else 0.0
            out_uncls[combined_sid] = u

    return out_counts, out_uncls


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
                if len(entries) >= 6:
                    node_level = (entries[3] or '').strip()
                    node_name  = entries[5].strip()
                    if not node_level and (node_name.lower() == 'root' or entries[4].strip() == '1'):
                        entries[3] = 'R'
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
    return return_values, unclassifieds


def output_kraken2_analytic_data(outfile, M, m_names, n_names, unclassifieds):
    with open(outfile, 'w') as out, \
            open(f"unclassifieds_{outfile}", 'w') as u_out:
        out.write(','.join(['taxa'] + n_names) + '\n')
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
    
    counts_by_sample, uncls = kraken2_load_analytic_data(opts.input_files)
    if opts.merged:
        counts_by_sample, uncls = combine_merged_unmerged(counts_by_sample, uncls)
    # Now build the matrix once
    K, m, n = dict_to_matrix(counts_by_sample)

    output_kraken2_analytic_data(opts.output_file, K, m, n, uncls)
