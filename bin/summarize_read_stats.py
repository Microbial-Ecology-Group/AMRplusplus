#!/usr/bin/env python3
"""
summarize_read_stats.py

Combine per-step seqkit stats into a single tidy spreadsheet with one row per
(sample, step, read-type). Handles single-end, paired-end, and merged+unmerged
layouts by inferring the read-type from each filename.

Typical use inside AMR++: each pipeline step (Raw, QC_trimmed, Deduped,
NonHost, ...) runs `seqkit stats -T -a` on its FASTQs and writes a
`<step>_seqkit.tsv`. This script is then run once on all of those to emit the
final combined workbook / CSV.

Output is a plain CSV (openpyxl / any spreadsheet library intentionally NOT used —
stdlib only, so it runs in any environment without extra dependencies). The CSV opens
directly in Excel / LibreOffice / Google Sheets.

Usage
-----
    # Combine several per-step seqkit tables:
    python summarize_read_stats.py \
        -i Raw_seqkit.tsv QC_trimmed_seqkit.tsv Deduped_seqkit.tsv NonHost_seqkit.tsv \
        -o read_stats_summary.csv
"""

import argparse
import os
import re
import sys
from collections import OrderedDict


# ── Read-type inference ──────────────────────────────────────────────────────
# Ordered: first pattern that matches a filename wins. Tweak here if your
# naming ever changes rather than touching the logic below.
READ_TYPE_PATTERNS = [
    (re.compile(r'\.extendedFrags\b|\bmerged\b(?!.*unmerged)', re.I), 'merged'),
    (re.compile(r'\.notCombined\b|\bunmerged\b',               re.I), 'unmerged'),
    (re.compile(r'_R1\b|\.1P\b|_1\b|\bR1\b',                    re.I), 'R1'),
    (re.compile(r'_R2\b|\.2P\b|_2\b|\bR2\b',                    re.I), 'R2'),
]


def infer_read_type(filename):
    """Return 'merged'|'unmerged'|'R1'|'R2'|'single' for a FASTQ filename."""
    base = os.path.basename(filename)
    for pattern, label in READ_TYPE_PATTERNS:
        if pattern.search(base):
            return label
    return 'single'


def infer_sample_id(filename):
    """
    Strip directory, common FASTQ extensions, read-type tags, and step tags
    to recover a stable sample id that is consistent across R1/R2/merged/etc.
    """
    base = os.path.basename(filename)
    # strip extensions
    base = re.sub(r'\.(fastq|fq)(\.gz)?$', '', base, flags=re.I)
    # strip read-type / pairing tags
    base = re.sub(r'(\.extendedFrags|\.notCombined)', '', base, flags=re.I)
    base = re.sub(r'([._-])(merged|unmerged)\b', '', base, flags=re.I)
    base = re.sub(r'([._-])(R?[12]P?)\b', '', base, flags=re.I)
    # strip common per-step suffixes so the same sample lines up across steps
    base = re.sub(r'\.(trimmed|dedup|non\.host|nonhost|host)\b', '', base, flags=re.I)
    return base


# ── seqkit table parsing ─────────────────────────────────────────────────────
def parse_seqkit_table(path, step_label):
    """
    Parse one `seqkit stats -T` (tab-separated) file.
    Returns a list of row dicts, one per input FASTQ line.
    """
    rows = []
    with open(path) as fh:
        header = fh.readline().rstrip('\n').split('\t')
        # seqkit -T columns: file format type num_seqs sum_len min_len avg_len max_len ...
        idx = {name: i for i, name in enumerate(header)}

        def col(parts, name, cast=str, default=''):
            if name in idx and idx[name] < len(parts):
                val = parts[idx[name]]
                try:
                    return cast(val)
                except (ValueError, TypeError):
                    return default
            return default

        for line in fh:
            line = line.rstrip('\n')
            if not line:
                continue
            parts = line.split('\t')
            fpath = col(parts, 'file')
            rows.append(OrderedDict([
                ('sample_id',  infer_sample_id(fpath)),
                ('step',       step_label),
                ('read_type',  infer_read_type(fpath)),
                ('num_reads',  col(parts, 'num_seqs', lambda x: int(float(x)), 0)),
                ('total_bp',   col(parts, 'sum_len',  lambda x: int(float(x)), 0)),
                ('min_len',    col(parts, 'min_len',  lambda x: int(float(x)), 0)),
                ('avg_len',    col(parts, 'avg_len',  float, 0.0)),
                ('max_len',    col(parts, 'max_len',  lambda x: int(float(x)), 0)),
                ('file',       os.path.basename(fpath)),
            ]))
    return rows


def step_label_from_filename(path):
    """Derive a step label from a '<step>_seqkit.tsv' filename."""
    base = os.path.basename(path)
    base = re.sub(r'_seqkit\.tsv$', '', base, flags=re.I)
    base = re.sub(r'\.(tsv|txt|csv)$', '', base, flags=re.I)
    return base or 'unknown'


# ── Aggregation: add per-sample merged+unmerged rollup ───────────────────────
def add_pair_rollups(rows):
    """
    For samples that have both 'merged' and 'unmerged' rows within a step,
    add a synthetic 'combined' row summing their read counts — so the
    spreadsheet has an at-a-glance total for merged-pipeline samples.
    """
    # index existing (sample, step, read_type)
    grouped = {}
    for r in rows:
        grouped.setdefault((r['sample_id'], r['step']), {})[r['read_type']] = r

    extra = []
    for (sample, step), by_type in grouped.items():
        if 'merged' in by_type and 'unmerged' in by_type:
            m, u = by_type['merged'], by_type['unmerged']
            extra.append(OrderedDict([
                ('sample_id', sample),
                ('step',      step),
                ('read_type', 'merged+unmerged (combined)'),
                ('num_reads', m['num_reads'] + u['num_reads']),
                ('total_bp',  m['total_bp']  + u['total_bp']),
                ('min_len',   min(m['min_len'], u['min_len'])),
                ('avg_len',   ''),   # avg not meaningful across two pools
                ('max_len',   max(m['max_len'], u['max_len'])),
                ('file',      ''),
            ]))
    return rows + extra


# ── Output writers ───────────────────────────────────────────────────────────
def write_csv(rows, out_path):
    import csv
    fieldnames = ['sample_id', 'step', 'read_type', 'num_reads',
                  'total_bp', 'min_len', 'avg_len', 'max_len', 'file']
    with open(out_path, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow(r)


# ── Sorting for readability ──────────────────────────────────────────────────
def sort_rows(rows, step_order):
    order_idx = {s: i for i, s in enumerate(step_order)}
    type_order = {'single': 0, 'R1': 1, 'R2': 2, 'merged': 3, 'unmerged': 4,
                  'merged+unmerged (combined)': 5}
    return sorted(
        rows,
        key=lambda r: (
            r['sample_id'],
            order_idx.get(r['step'], 999),
            type_order.get(r['read_type'], 9),
        )
    )


# ── Main ─────────────────────────────────────────────────────────────────────
def parse_args():
    p = argparse.ArgumentParser(
        description="Combine per-step seqkit stats into one read-stats spreadsheet."
    )
    p.add_argument('-i', '--input_files', nargs='+', required=True,
                   help="One or more `seqkit stats -T` TSV files, typically "
                        "named <step>_seqkit.tsv")
    p.add_argument('-o', '--output_file', required=True,
                   help="Output CSV path (opens directly in Excel/LibreOffice/Sheets).")
    p.add_argument('--step-order', nargs='+', default=None,
                   help="Optional explicit ordering of step labels for output "
                        "(e.g. --step-order Raw QC_trimmed Deduped NonHost).")
    return p.parse_args()


def main():
    args = parse_args()

    all_rows = []
    seen_steps = []
    for path in args.input_files:
        if not os.path.isfile(path):
            print(f"[WARN] Skipping missing file: {path}", file=sys.stderr)
            continue
        step = step_label_from_filename(path)
        if step not in seen_steps:
            seen_steps.append(step)
        all_rows.extend(parse_seqkit_table(path, step))

    if not all_rows:
        sys.exit("ERROR: No data parsed from any input file.")

    all_rows = add_pair_rollups(all_rows)

    step_order = args.step_order if args.step_order else seen_steps
    all_rows = sort_rows(all_rows, step_order)

    write_csv(all_rows, args.output_file)

    n_samples = len({r['sample_id'] for r in all_rows})
    print(f"[INFO] Wrote {len(all_rows)} rows for {n_samples} sample(s) "
          f"across {len(step_order)} step(s) → {args.output_file}",
          file=sys.stderr)


if __name__ == '__main__':
    main()
