#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Filter an SNV count matrix to the genes that were also detected in the resistome
count matrix.

Matching is on the MEG_ accession, the first pipe-delimited field of each
identifier, so an SNV row is kept when its gene appears in the resistome matrix
regardless of the variant suffix appended to the SNV accession.

CHANGED FROM THE ORIGINAL
-------------------------
The original version hardcoded four input and output paths under
'test_results/Results/' and selected between them with -aln_wf. That works for a
fixed directory layout but not inside a Nextflow task, where inputs are staged
into an isolated work directory under whatever names the process declares.

Paths are now explicit arguments. -aln_wf is retained so existing calls keep
working: when the three path arguments are omitted it reproduces the original
behaviour exactly.

Usage
-----
    # explicit paths (used by the AMR++ SNV workflow)
    python3 Filter_SNV_matrix.py \\
        --snv-matrix resistome_SNV_analytic_matrix.csv \\
        --resistome-matrix AMR_analytic_matrix.csv \\
        --output filtered_resistome_SNV_analytic_matrix.csv

    # original behaviour, standard workflow
    python3 Filter_SNV_matrix.py -aln_wf Standard
"""

import argparse
import os
import sys

import pandas as pd


# Original hardcoded layout, kept so -aln_wf alone still works.
DEFAULTS = {
    "Standard": {
        "snv":  "test_results/Results/resistome_SNV_analytic_matrix.csv",
        "res":  "test_results/Results/AMR_analytic_matrix.csv",
        "out":  "test_results/Results/filtered_resistome_SNV_analytic_matrix.csv",
    },
    "Deduped": {
        "snv":  "test_results/Results/dedup_resistome_SNV_analytic_matrix.csv",
        "res":  "test_results/Results/dedup_AMR_analytic_matrix.csv",
        "out":  "test_results/Results/filtered_dedup_resistome_SNV_analytic_matrix.csv",
    },
}

SNV_ID_COLUMN = "SNV_accession"
RES_ID_COLUMN = "gene_accession"


def main():
    p = argparse.ArgumentParser(
        description="Filter an SNV matrix to ARGs present in the resistome matrix.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("-aln_wf", type=str, default="Standard",
                   choices=["Standard", "Deduped"],
                   help="Which default path set to use when explicit paths are "
                        "not given (default: Standard).")
    p.add_argument("--snv-matrix", default=None,
                   help="SNV count matrix CSV. Overrides the -aln_wf default.")
    p.add_argument("--resistome-matrix", default=None,
                   help="Resistome count matrix CSV. Overrides the -aln_wf default.")
    p.add_argument("--output", default=None,
                   help="Output CSV. Overrides the -aln_wf default.")
    args = p.parse_args()

    d = DEFAULTS[args.aln_wf]
    snv_path = args.snv_matrix      or d["snv"]
    res_path = args.resistome_matrix or d["res"]
    out_path = args.output           or d["out"]

    for label, path in (("SNV matrix", snv_path), ("resistome matrix", res_path)):
        if not os.path.exists(path):
            sys.exit(f"[ERROR] {label} not found: {path}")

    snv_df = pd.read_csv(snv_path)
    res_df = pd.read_csv(res_path)

    for label, df, col, path in (("SNV", snv_df, SNV_ID_COLUMN, snv_path),
                                 ("resistome", res_df, RES_ID_COLUMN, res_path)):
        if col not in df.columns:
            sys.exit(f"[ERROR] {label} matrix {path} has no '{col}' column. "
                     f"Found: {', '.join(df.columns[:8])}")

    # Match on the MEG_ accession, the first pipe-delimited field.
    allowed = set(res_df[RES_ID_COLUMN].astype(str).str.split("|").str[0])
    snv_ids = snv_df[SNV_ID_COLUMN].astype(str).str.split("|").str[0]
    filtered = snv_df[snv_ids.isin(allowed)]

    out_dir = os.path.dirname(out_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    filtered.to_csv(out_path, index=False)

    kept, total = len(filtered), len(snv_df)
    pct = kept / total * 100 if total else 0.0
    print(f"[INFO] {kept:,} of {total:,} SNV rows retained ({pct:.1f}%)")
    print(f"[INFO] {len(allowed):,} distinct genes in the resistome matrix")
    if kept == 0 and total > 0:
        print("[WARN] No SNV rows matched. Check that both matrices were "
              "generated from the same reference database.")
    print(f"[INFO] Wrote {out_path}")


if __name__ == "__main__":
    main()
