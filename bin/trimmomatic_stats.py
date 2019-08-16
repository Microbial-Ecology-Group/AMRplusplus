#!/usr/bin/env python3

__author__ = "Chris Dean"
__copyright__ = ""
__credits__ = ["Chris Dean"]
__version__ = ""
__maintainer__ = "cdeanj"
__email__ = "cdean11@colostate.edu"
__status__ = "Cows go moo."

import argparse
import glob
import os
import re
import sys

def parse_cmdline_params(cmdline_params):
    info = "Parses a Trimmomatic log file to obtain the total number of input reads and dropped reads"
    parser = argparse.ArgumentParser(description=info)
    parser.add_argument('-i', '--input_files', nargs='+', required=True,
                        help='Use globstar to pass a list of files, (Ex: *.tsv)')
    parser.add_argument('-o', '--output_file', required=True,
                        help='Output file to write mapping results to')
    return parser.parse_args(cmdline_params)

def header(output_file):
    with open(output_file, 'a') as o:
        o.write('Sample\tNumberOfInputReads\tForwardOnlySurviving\tReverseOnlySurviving\tDropped\n')
    o.close()

def qc_stats(input_list, output_file):
    for f in input_list:
        total = 0
        forward_surviving = 0
        reverse_surviving = 0
        dropped = 0
        with open(f, 'r') as fp:
            sample_name = os.path.basename(str(fp.name)).split('.', 1)[0]
            for line in fp:
                total = re.search('Input Read Pairs: (\d+)', line)
                forward_surviving = re.search('Forward Only Surviving: (\d+)', line)
                reverse_surviving = re.search('Reverse Only Surviving: (\d+)', line)
                dropped = re.search('Dropped: (\d+)', line)
                if total:
                    break
        fp.close()
        with open(output_file, 'a') as ofp:
            ofp.write(sample_name + '\t' + total.group(1) + '\t' + forward_surviving.group(1) + '\t' + reverse_surviving.group(1) + '\t' + dropped.group(1) + '\n')
        ofp.close()

if __name__ == "__main__":
    opts = parse_cmdline_params(sys.argv[1:])
    header(opts.output_file)
    qc_stats(opts.input_files, opts.output_file)
