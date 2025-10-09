#!/usr/bin/env python3
"""
compare_fastqs.py

Compare two directories of FASTQ files (whole_set vs subset) and create a
shell script of `cp` commands to move unmatched whole_set files into an
output directory.

Usage
-----
python3 compare_fastqs.py \
    --whole_set_dir /path/to/whole_set_fastqs --whole_set_suffix .fastq.gz \
    --subset_dir /path/to/subset_fastqs --subset_suffix .subset.fastq.gz \
    --output_dir remaining_150
"""

import argparse
import os
from pathlib import Path


def parse_args():
    p = argparse.ArgumentParser(
        description="Generate cp commands to move whole_set FASTQs without a matching subset FASTQ."
    )
    p.add_argument("--whole_set_dir", required=True, help="Directory containing whole_set FASTQ files")
    p.add_argument("--whole_set_suffix", required=True, help="Suffix of whole_set FASTQ files (e.g. .fastq.gz)")
    p.add_argument("--subset_dir", required=True, help="Directory containing subset FASTQ files")
    p.add_argument("--subset_suffix", required=True, help="Suffix of subset FASTQ files (e.g. .subset.fastq.gz)")
    p.add_argument("--output_dir", required=True, help="Name of the target directory to copy missing samples into")
    return p.parse_args()


def main():
    args = parse_args()
    whole_set_dir = Path(args.whole_set_dir)
    subset_dir = Path(args.subset_dir)

    # collect basenames (without suffixes)
    whole_set_files = {
        f.name[:-len(args.whole_set_suffix)]
        for f in whole_set_dir.glob(f"*{args.whole_set_suffix}")
        if f.is_file() and f.name.endswith(args.whole_set_suffix)
    }

    subset_files = {
        f.name[:-len(args.subset_suffix)]
        for f in subset_dir.glob(f"*{args.subset_suffix}")
        if f.is_file() and f.name.endswith(args.subset_suffix)
    }

    # find files present in whole_set but not in subset
    missing = sorted(whole_set_files - subset_files)

    # write shell script
    out_script = f"move_samples_to_{args.output_dir}.sh"
    with open(out_script, "w") as out:
        out.write("#!/bin/bash\n\n")
        out.write(f"mkdir -p {args.output_dir}\n\n")
        for stem in missing:
            whole_set_path = whole_set_dir / f"{stem}{args.whole_set_suffix}"
            out.write(f"cp {whole_set_path} {args.output_dir}/\n")

    print(f"Wrote {len(missing)} cp commands to {out_script}")


if __name__ == "__main__":
    main()

