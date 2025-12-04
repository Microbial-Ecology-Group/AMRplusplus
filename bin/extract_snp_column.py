#!/usr/bin/env python3
import argparse, sys, os, csv, gzip

def smart_open(path, mode="rt"):
    """Open plain or gzipped files transparently."""
    if str(path).endswith(".gz"):
        return gzip.open(path, mode=mode, encoding=None if "b" in mode else "utf-8")
    return open(path, mode=mode, encoding=None if "b" in mode else "utf-8")

def main():
    parser = argparse.ArgumentParser(
        description="Extract first column and a given sample column from a SNP count matrix CSV and write as TSV."
    )
    parser.add_argument("-s", "--sample-id", required=True, help="Sample ID to find in header")
    parser.add_argument("-m", "--matrix", required=True, help="Path to SNP count matrix CSV (can be .gz)")
    parser.add_argument("-o", "--out-tsv", required=True, help="Output TSV path")
    args = parser.parse_args()

    # --- Handle '.non.host' suffix if present ---
    if ".non.host" in args.sample_id:
        original_id = args.sample_id
        args.sample_id = args.sample_id.split(".non.host")[0]
        print(f"[INFO] Detected '.non.host' in sample ID. "
              f"Using '{args.sample_id}' instead of '{original_id}'")

    # --- Read matrix and find the target column ---
    with smart_open(args.matrix, "rt") as f:
        reader = csv.reader(f)
        try:
            header = next(reader)
        except StopIteration:
            sys.exit(f"[ERROR] Empty file: {args.matrix}")

        if args.sample_id not in header:
            sys.exit(f"[ERROR] Sample '{args.sample_id}' not found in header of {args.matrix}")

        sample_col_idx = header.index(args.sample_id)
        first_col_name = header[0]

        # Create output directory if needed
        os.makedirs(os.path.dirname(args.out_tsv) or ".", exist_ok=True)

        with open(args.out_tsv, "w", newline="") as out_f:
            writer = csv.writer(out_f, delimiter="\t")
            writer.writerow([first_col_name, args.sample_id])  # Header
            for row in reader:
                if len(row) > sample_col_idx:
                    writer.writerow([row[0], row[sample_col_idx]])

if __name__ == "__main__":
    main()
