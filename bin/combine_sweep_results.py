#!/usr/bin/env python3
"""
combine_sweep_results.py
════════════════════════════════════════════════════════════════════════════════

Combines per-sample outputs from coverage_threshold_sweep.py (run one sample at
a time) into single matrices — one per output file type — with a sample_id
column injected from each file's name, so every sample and file type becomes one
queryable table, ready for plot_sweep_dropoff.R.

Scans a single directory for files matching the known suffixes:
    <prefix>_results.csv, <prefix>_gene_detail.csv,
    <prefix>_redundancy.csv, <prefix>_length_quantiles.csv
where <prefix> is the BAM filename minus ".bam". From that:
    sample_id  = the full prefix (BAM basename), preserved as-is.

OUTPUT
──────
combined_results.csv, combined_gene_detail.csv, combined_redundancy.csv,
combined_length_quantiles.csv (only the types that actually exist in the input
are produced) — each with a sample_id column prepended. Written progressively;
each output file is cleared once at the start of this run.

Usage:
    python combine_sweep_results.py sweep_results/ --out-dir combined/
"""

import argparse
import csv
import sys
from pathlib import Path
from typing import List, Optional, Tuple

FILE_TYPES = {
    "_results.csv":          "results",
    "_gene_detail.csv":      "gene_detail",
    "_redundancy.csv":       "redundancy",
    "_length_quantiles.csv": "length_quantiles",
}


def classify_filename(fname: str) -> Optional[Tuple[str, str]]:
    """Returns (prefix, file_type) if fname matches a known suffix, else None."""
    for suffix, ftype in FILE_TYPES.items():
        if fname.endswith(suffix):
            return fname[: -len(suffix)], ftype
    return None


def make_combiner(out_dir: Path):
    """
    Per-file-type streaming appender: header written on the first row for that
    file type, appended after. Clears any pre-existing combined_*.csv from a
    previous run the first time each file type is touched, so reruns don't
    silently double up data.
    """
    written = set()

    def append_rows(rows: List[dict], file_type: str):
        if not rows:
            return
        p = out_dir / f"combined_{file_type}.csv"
        first = file_type not in written
        if first and p.exists():
            p.unlink()
        with open(p, "a", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            if first:
                writer.writeheader()
            writer.writerows(rows)
        if first:
            written.add(file_type)
            print(f"  created combined_{file_type}.csv", flush=True)

    return append_rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input_dir", help="Directory holding the per-sample sweep CSVs")
    ap.add_argument("--out-dir", default="combined_sweep_results")
    args = ap.parse_args()

    root = Path(args.input_dir)
    if not root.is_dir():
        sys.exit(f"[ERROR] {root} is not a directory")

    sweep_files = sorted(
        f for f in root.iterdir() if f.is_file() and classify_filename(f.name)
    )
    if not sweep_files:
        sys.exit(f"[ERROR] No matching *_results.csv / *_gene_detail.csv / etc. files "
                 f"found in {root}")

    print(f"[INFO] Found {len(sweep_files)} sweep file(s) in {root}", flush=True)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    append_rows = make_combiner(out_dir)

    counts = {ft: 0 for ft in set(FILE_TYPES.values())}
    n_files = 0

    for f in sweep_files:
        prefix, file_type = classify_filename(f.name)
        sample_id = prefix  # full BAM basename; no token-splitting

        with open(f, newline="") as fh:
            reader = csv.DictReader(fh)
            if not reader.fieldnames:
                continue
            first_col = reader.fieldnames[0]
            rows = []
            for row in reader:
                # Defensive: skip an embedded header row, in case this file was
                # ever manually concatenated before reaching us.
                if row.get(first_col) == first_col:
                    continue
                rows.append({"sample_id": sample_id, **row})

        append_rows(rows, file_type)
        counts[file_type] += len(rows)
        n_files += 1

    print(f"\n[INFO] Combined {n_files} file(s):", flush=True)
    for ft, n in counts.items():
        if n:
            print(f"    {ft:<18} {n:,} rows -> combined_{ft}.csv", flush=True)
    print(f"\nOutput: {out_dir}/", flush=True)


if __name__ == "__main__":
    main()
