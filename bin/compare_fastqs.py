#!/usr/bin/env python3
from __future__ import annotations
"""
compare_fastqs.py

Compare two directories of FASTQ files (whole_set vs subset) and create a
shell script of `cp` commands to move unmatched whole_set files into an
output directory.

Usage
-----
python3 compare_fastqs.py \
    --whole_set_dir /path/to/whole_set_fastqs \
    --whole_set_pattern "*{merged,unmerged}.fastq.gz" \
    --subset_dir /path/to/subset_fastqs \
    --subset_pattern "*R{1,2}.fastq.gz" \
    --output_dir remaining_150

Patterns support brace expansion like bash:
  - "*{merged,unmerged}.fastq.gz" matches files ending in merged.fastq.gz or unmerged.fastq.gz
  - "*R{1,2}.fastq.gz" matches files ending in R1.fastq.gz or R2.fastq.gz
  - "*.fastq.gz" matches any file ending in .fastq.gz (sample ID = everything before .fastq.gz)
"""

import argparse
import re
from pathlib import Path
from itertools import product


def parse_args():
    p = argparse.ArgumentParser(
        description="Generate cp commands to move whole_set FASTQs without a matching subset FASTQ."
    )
    p.add_argument(
        "--whole_set_dir", required=True,
        help="Directory containing whole_set FASTQ files"
    )
    p.add_argument(
        "--whole_set_pattern", required=True,
        help="Pattern for whole_set FASTQ files (e.g. '*{merged,unmerged}.fastq.gz' or '*.fastq.gz')"
    )
    p.add_argument(
        "--subset_dir", required=True,
        help="Directory containing subset FASTQ files"
    )
    p.add_argument(
        "--subset_pattern", required=True,
        help="Pattern for subset FASTQ files (e.g. '*R{1,2}.fastq.gz' or '*.subset.fastq.gz')"
    )
    p.add_argument(
        "--output_dir", required=True,
        help="Name of the target directory to copy missing samples into"
    )
    return p.parse_args()


def expand_braces(pattern: str) -> list[str]:
    """
    Expand brace expressions in a pattern, similar to bash brace expansion.
    
    Examples:
        "*{merged,unmerged}.fastq.gz" -> ["*merged.fastq.gz", "*unmerged.fastq.gz"]
        "*R{1,2}.fastq.gz" -> ["*R1.fastq.gz", "*R2.fastq.gz"]
        "*.fastq.gz" -> ["*.fastq.gz"]
        "*{A,B}{1,2}.txt" -> ["*A1.txt", "*A2.txt", "*B1.txt", "*B2.txt"]
    """
    brace_pattern = re.compile(r'\{([^{}]+)\}')
    matches = list(brace_pattern.finditer(pattern))
    
    if not matches:
        return [pattern]
    
    # Extract alternatives for each brace group
    alternatives_list = []
    for match in matches:
        alternatives = match.group(1).split(',')
        alternatives_list.append(alternatives)
    
    # Generate all combinations
    expanded = []
    for combo in product(*alternatives_list):
        result = pattern
        # Replace braces from right to left to preserve indices
        for match, replacement in zip(reversed(matches), reversed(combo)):
            result = result[:match.start()] + replacement + result[match.end():]
        expanded.append(result)
    
    return expanded


def extract_sample_id(filename: str, pattern: str) -> str | None:
    """
    Extract the sample ID from a filename based on the pattern.
    The sample ID is the part matched by the first '*' in the pattern.
    """
    expanded = expand_braces(pattern)
    
    for exp in expanded:
        if '*' not in exp:
            # No wildcard - exact match, use filename as ID
            if filename == exp:
                return filename
            continue
        
        # Split pattern at first *
        prefix, rest = exp.split('*', 1)
        
        # Check if filename starts with prefix
        if not filename.startswith(prefix):
            continue
        
        # Get all possible suffixes from the rest of the pattern
        if '{' in rest:
            possible_suffixes = expand_braces(rest)
        else:
            possible_suffixes = [rest]
        
        for possible_suffix in possible_suffixes:
            # Handle additional wildcards in suffix by converting to regex
            suffix_regex = re.escape(possible_suffix).replace(r'\*', '.*')
            suffix_match = re.search(suffix_regex + '$', filename[len(prefix):])
            if suffix_match:
                # Sample ID is between the prefix and where the suffix starts
                sample_id = filename[len(prefix):len(prefix) + suffix_match.start()]
                return sample_id
    
    return None


def find_files_with_pattern(directory: Path, pattern: str) -> tuple[dict[str, list[Path]], int]:
    """
    Find all files matching the pattern and group them by sample ID.
    
    Returns:
        - dict mapping sample_id -> list of matching file paths
        - total number of files found
    """
    expanded_patterns = expand_braces(pattern)
    
    matching_files: dict[str, list[Path]] = {}
    total_files = 0
    
    for exp in expanded_patterns:
        for f in directory.glob(exp):
            if f.is_file():
                sample_id = extract_sample_id(f.name, pattern)
                if sample_id is not None:
                    if sample_id not in matching_files:
                        matching_files[sample_id] = []
                    matching_files[sample_id].append(f)
                    total_files += 1
    
    return matching_files, total_files


def main():
    args = parse_args()
    
    whole_set_dir = Path(args.whole_set_dir)
    subset_dir = Path(args.subset_dir)
    
    print("=" * 60)
    print("FASTQ Comparison Report")
    print("=" * 60)
    
    # Show expanded patterns for clarity
    whole_set_expanded = expand_braces(args.whole_set_pattern)
    subset_expanded = expand_braces(args.subset_pattern)
    
    print(f"\nWhole set directory: {whole_set_dir}")
    print(f"  Pattern: {args.whole_set_pattern}")
    if len(whole_set_expanded) > 1:
        print(f"  Expands to: {', '.join(whole_set_expanded)}")
    
    print(f"\nSubset directory: {subset_dir}")
    print(f"  Pattern: {args.subset_pattern}")
    if len(subset_expanded) > 1:
        print(f"  Expands to: {', '.join(subset_expanded)}")
    
    # Find files and extract sample IDs
    whole_set_files, whole_set_count = find_files_with_pattern(whole_set_dir, args.whole_set_pattern)
    subset_files, subset_count = find_files_with_pattern(subset_dir, args.subset_pattern)
    
    print("\n" + "-" * 60)
    print("Files Processed")
    print("-" * 60)
    print(f"Whole set: {whole_set_count} files")
    print(f"Subset:    {subset_count} files")
    
    print("\n" + "-" * 60)
    print("Sample IDs Found")
    print("-" * 60)
    print(f"Whole set: {len(whole_set_files)} unique sample IDs")
    print(f"Subset:    {len(subset_files)} unique sample IDs")
    
    # Find sample IDs present in whole_set but not in subset
    missing_ids = sorted(set(whole_set_files.keys()) - set(subset_files.keys()))
    common_ids = set(whole_set_files.keys()) & set(subset_files.keys())
    
    print("\n" + "-" * 60)
    print("Comparison Results")
    print("-" * 60)
    print(f"Sample IDs in both:              {len(common_ids)}")
    print(f"Sample IDs missing from subset:  {len(missing_ids)}")
    
    # Count total files to be copied
    files_to_copy = sum(len(whole_set_files[sid]) for sid in missing_ids)
    print(f"Files to be copied:              {files_to_copy}")
    
    # Write shell script
    out_script = f"move_samples_to_{args.output_dir}.sh"
    
    with open(out_script, "w") as out:
        out.write("#!/bin/bash\n\n")
        out.write(f"mkdir -p {args.output_dir}\n\n")
        
        for sample_id in missing_ids:
            for filepath in sorted(whole_set_files[sample_id]):
                out.write(f"cp {filepath} {args.output_dir}/\n")
    
    print("\n" + "=" * 60)
    print(f"Output: {out_script}")
    print(f"  Contains {files_to_copy} cp commands for {len(missing_ids)} samples")
    print("=" * 60)


if __name__ == "__main__":
    main()
