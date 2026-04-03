#!/usr/bin/env python3
"""
combine_seqkit_stats.py

Finds tab-delimited text files (e.g., seqkit stats output) matching a regex
pattern under a given directory, adds a 'data_group' column derived from each
file's stem, and writes a single combined TSV.

Usage:
    python combine_seqkit_stats.py \
        --path /path/to/results \
        --pattern ".*counts.*\\.txt$" \
        --output combined_stats.txt

    # Or search non-recursively:
    python combine_seqkit_stats.py \
        --path /path/to/results \
        --pattern "raw_counts\\.txt$" \
        --output combined_stats.txt \
        --no-recursive
"""

import argparse
import re
import sys
from pathlib import Path

import pandas as pd


def parse_args():
    parser = argparse.ArgumentParser(
        description="Combine seqkit stats (or similar) TSV files into one, "
                    "adding a 'data_group' column from each file's stem."
    )
    parser.add_argument(
        "--path", "-p",
        required=True,
        help="Root directory to search for matching files."
    )
    parser.add_argument(
        "--pattern", "-r",
        required=True,
        help="Regular expression matched against each file's full path. "
             "Example: '.*raw_counts\\.txt$' or '.*stats.*\\.txt$'"
    )
    parser.add_argument(
        "--output", "-o",
        required=True,
        help="Output file path for the combined TSV."
    )
    parser.add_argument(
        "--suffix",
        default=None,
        help="Explicit suffix to strip from filenames when building "
             "data_group (e.g. '.txt'). Defaults to stripping the last "
             "extension via Path.stem."
    )
    parser.add_argument(
        "--no-recursive",
        action="store_true",
        default=False,
        help="Search only the top-level directory (non-recursive). "
             "Default is recursive."
    )
    parser.add_argument(
        "--sep",
        default="\t",
        help="Field separator used in input files (default: tab)."
    )
    return parser.parse_args()


def find_files(root: Path, pattern: re.Pattern, recursive: bool) -> list[Path]:
    """Return all files under root whose full path matches pattern."""
    glob = root.rglob("*") if recursive else root.glob("*")
    matched = sorted(
        f for f in glob
        if f.is_file() and pattern.search(str(f))
    )
    return matched


def file_stem(filepath: Path, suffix: str | None) -> str:
    """Return the data_group label for a file."""
    if suffix:
        name = filepath.name
        return name[: -len(suffix)] if name.endswith(suffix) else name
    return filepath.stem          # strips only the last extension


def read_file(filepath: Path, sep: str) -> pd.DataFrame | None:
    """Read a TSV/CSV file into a DataFrame, skipping comment lines."""
    try:
        df = pd.read_csv(filepath, sep=sep, comment="#")
        if df.empty:
            print(f"  [warn] {filepath} is empty — skipping.", file=sys.stderr)
            return None
        return df
    except Exception as exc:
        print(f"  [warn] Could not read {filepath}: {exc} — skipping.",
              file=sys.stderr)
        return None


def main():
    args = parse_args()

    root = Path(args.path)
    if not root.is_dir():
        sys.exit(f"ERROR: --path '{args.path}' is not a valid directory.")

    try:
        pattern = re.compile(args.pattern)
    except re.error as exc:
        sys.exit(f"ERROR: Invalid regex '{args.pattern}': {exc}")

    recursive = not args.no_recursive
    files = find_files(root, pattern, recursive)

    if not files:
        sys.exit(
            f"ERROR: No files matched pattern '{args.pattern}' "
            f"under '{root}' ({'recursive' if recursive else 'non-recursive'})."
        )

    print(f"Found {len(files)} matching file(s):")
    frames = []
    for f in files:
        group = file_stem(f, args.suffix)
        print(f"  {f}  →  data_group='{group}'")
        df = read_file(f, args.sep)
        if df is None:
            continue
        df.insert(0, "data_group", group)   # prepend so it's the first column
        frames.append(df)

    if not frames:
        sys.exit("ERROR: No files could be read successfully.")

    combined = pd.concat(frames, ignore_index=True)

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    combined.to_csv(out, sep="\t", index=False)

    print(f"\nCombined {len(frames)} file(s) → {out}")
    print(f"Total rows : {len(combined):,}")
    print(f"Columns    : {list(combined.columns)}")


if __name__ == "__main__":
    main()