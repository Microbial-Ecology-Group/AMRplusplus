# -*- coding: utf-8 -*-
"""
Modified on Wed Jul 15 2026

@author: Jared Young with edits from Enrique Doster
"""

import csv
import os
import argparse

def sanitize_lines(file_handle):
    """Generator to strip NULL bytes (\x00) before feeding to csv.reader."""
    for line in file_handle:
        yield line.replace('\0', '')

def main():
    parser = argparse.ArgumentParser(
        description='Generates a concatenated count matrix and annotation file for the raw files given',
        usage='python %(prog)s [-h] [-s save_dir] file_dir fn_struc sample_file'
    )

    parser.add_argument('-s', type=str, default='test_results/Results/', help='Path in which to store the generated files')
    parser.add_argument('-file_dir', type=str, default='test_results/SNV_analysis_output/snpCaller/', help='Path where the raw files are stored')
    parser.add_argument('-fn_struc', type=str, default='called_SNPs', help='Naming structure uses in all raw files, e.g., "called_SNPs.best_split_"')
    parser.add_argument('-sample_file', type=str, default='test_results/SNV_analysis_output/all_samples', help='Path to the file with the names of all the samples')
    parser.add_argument('-aln_wf', type=str, default='Standard', help='Alignment type workflow - Standard (default) or Deduped')
    parser.add_argument('-matrix_out', type=str, default=None, help='Filename for the output analytic matrix (optional; default set by -aln_wf)')
    parser.add_argument('-ann_out', type=str, default='resistome_SNV_megares_v4_annotations.csv', help='Filename for the output annotations file')

    args = parser.parse_args()

    # Set default matrix output name if not explicitly passed
    if args.matrix_out is None:
        if args.aln_wf == "Deduped":
            args.matrix_out = 'dedup_resistome_SNV_analytic_matrix.csv'
        else:
            args.matrix_out = 'resistome_SNV_analytic_matrix.csv'

    data_frame = []

    for fn in os.listdir(args.file_dir):
        if args.fn_struc in fn:
            file_path = os.path.join(args.file_dir, fn)
            with open(file_path, 'r', encoding='utf-8', errors='replace') as file:
                reader = csv.reader(sanitize_lines(file), delimiter='\t')
                for row in reader:
                    if not row:
                        continue
                    if row[0].split('|')[-1] == "RequiresSNPConfirmation":
                        continue
                    extras_removed = row[0].split('|')[:5]
                    asc = extras_removed[0].split('_')[1]
                    grp_asc = '-'.join([extras_removed[4]] + [asc])
                    snp_ann = '_'.join([grp_asc] + [row[3].upper(), row[2]])
                    joined = '|'.join(extras_removed + [snp_ann])
                    data_frame.append([joined] + row[5].split(','))

    # Adjusting data frame structure
    for i in range(len(data_frame)):
        while len(data_frame[i]) > 2:
            data_frame.append([data_frame[i][0], data_frame[i][-1]])
            data_frame[i] = data_frame[i][:-1]

    for i in range(len(data_frame)):
        data_frame[i] = [data_frame[i][0]] + data_frame[i][1].split('|')
        data_frame[i] = [data_frame[i][0], data_frame[i][2]] + data_frame[i][4:]
        data_frame[i][0] = '_'.join([data_frame[i][0], data_frame[i][1].upper()])
        data_frame[i] = [data_frame[i][0]] + data_frame[i][2:]

    # Reading in sample names
    all_samples = []
    s_add = "S_"

    with open(args.sample_file, 'r', encoding='utf-8', errors='replace') as file:
        for line in sanitize_lines(file):
            line = line.rstrip()
            if not line:
                continue
            if "." in line and "/" in line:
                line = line.split("/")[-1].split(".")[0]
            elif "." in line:
                line = line.split(".")[0]
            elif "/" in line:
                line = line.split("/")[-1]

            if line and line[0].isdigit():
                line = s_add + line
            all_samples.append(line)

    # Adding title row to data frame
    data_frame.insert(0, ['SNV_accession'] + all_samples)

    # Writing final count matrix
    os.makedirs(args.s, exist_ok=True)
    matrix_path = os.path.join(args.s, args.matrix_out)

    with open(matrix_path, 'w', newline='', encoding='utf-8') as file:
        writer = csv.writer(file)
        writer.writerows(data_frame)

    # Generating annotation file data frame
    ann_data_frame = []
    for row in data_frame[1:]:
        ann_data_frame.append([row[0], *row[0].split('|')[1:]])

    title_row = ['SNV_accession', 'Type', 'Class', 'Mechanism', 'Group', 'SNV']
    ann_data_frame.insert(0, title_row)

    annotation_path = os.path.join(args.s, args.ann_out)
    with open(annotation_path, 'w', newline='', encoding='utf-8') as file:
        writer = csv.writer(file)
        writer.writerows(ann_data_frame)

if __name__ == '__main__':
    main()