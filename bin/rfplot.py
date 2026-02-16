#!/usr/bin/env python3
"""
Plot rarefaction curves and export raw counts.

Usage example
-------------
python rfplot.py \
       --dir results/Temp_rarefaction \
       --s --sd figs/ \
       --prefix myRun
"""
import argparse
import csv
import os
from pathlib import Path

import matplotlib.pyplot as plt

# ────────────────────────── CLI ──────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("--dir", type=str, required=True,
                    help="Directory containing the *.tsv rarefaction files")
parser.add_argument("--nd", action="store_true",
                    help="No display: do not open the Matplotlib windows")
parser.add_argument("--s", action="store_true",
                    help="Save the generated figures")
parser.add_argument("--sd", type=str, default="",
                    help="Directory in which to save figures (needs -s)")
parser.add_argument("--prefix", type=str, default="out",
                    help="Prefix for output figure names and counts table")
args = parser.parse_args()

indir  = Path(args.dir).expanduser().resolve()
savedir = Path(args.sd).expanduser().resolve() if args.sd else Path.cwd()
savedir.mkdir(parents=True, exist_ok=True)

# ───────────────────── storage for plotting ──────────────────────
levels = ("gene", "group", "mech", "class", "type")
x_dict = {lv: [] for lv in levels}   # list of lists → per-sample X
y_dict = {lv: [] for lv in levels}   # list of lists → per-sample Y
sample_names = []

# collect rows for the counts-table
counts_rows = []

# ─────────────────── ingest every TSV file ───────────────────────
for fn in sorted(indir.glob("*.tsv")):
    parts = fn.stem.split(".")            # e.g. SAMPLE.gene
    if len(parts) < 2 or parts[-1] not in levels:
        continue

    sample = parts[0]
    level  = parts[-1]

    if sample not in sample_names:
        sample_names.append(sample)

    x_dict[level].append([])
    y_dict[level].append([])

    with fn.open() as fh:
        tsv = csv.reader(fh, delimiter="\t")
        for row in tsv:
            subsample_pct = int(row[0])
            count         = int(row[1])

            # for plots
            x_dict[level][-1].append(subsample_pct)
            y_dict[level][-1].append(count)

            # for summary file
            counts_rows.append(
                (sample, level, count, subsample_pct)
            )

# ──────────────────── write counts table ─────────────────────────
counts_file = Path(f"rarefaction_{args.prefix}_counts.txt")
with counts_file.open("w", newline="") as fh:
    writer = csv.writer(fh, delimiter="\t")
    writer.writerow(["Sample", "Level", "Count", "SubsamplePercent"])
    writer.writerows(counts_rows)

print(f"✓ wrote count table → {counts_file}")

# ────────────────────────── plotting ─────────────────────────────
fig_axes = [plt.subplots() for _ in levels]   # list of (fig, ax) tuples

for lv_idx, level in enumerate(levels):
    fig, ax = fig_axes[lv_idx]
    for i, sample in enumerate(sample_names):
        if i < len(x_dict[level]):            # sample had that level
            ax.plot(x_dict[level][i],
                    y_dict[level][i],
                    label=sample)

    ax.set_title(f"{level.capitalize()} Subsampling Features")
    ax.set_xlabel("% of data subsampled")
    ax.set_ylabel("unique features identified")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.legend(bbox_to_anchor=(1.05, 1.0), loc="upper left")
    fig.set_facecolor("white")
    fig.tight_layout()

# ──────────────────── show or save figures ───────────────────────
if not args.nd:
    plt.show()

if args.s:
    for (fig, _), level in zip(fig_axes, levels):
        outfile = savedir / f"{args.prefix}_{level.capitalize()}.png"
        fig.savefig(outfile, dpi=300)
        print(f"✓ saved {outfile}")
