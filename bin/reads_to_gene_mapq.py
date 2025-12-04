#!/usr/bin/env python3
import argparse
import os
import sys
from collections import Counter
from typing import Dict, Any, Tuple

import pysam

## Example command
#    python3 reads_to_gene_mapq.py \
#        -i "$bam" \
#        -r "counts/${sample_id}_per_read.tsv" \
#        -g "counts/${sample_id}_gene_summary.tsv" \
#        --sample-id "$sample_id" \
#        --min-mapq 0 \
#        --count-mode alignment

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "For each read (read-end) in a BAM/SAM, output the gene(s) it maps to "
            "via primary alignments and the MAPQ, plus a summary of counts per gene."
        )
    )
    parser.add_argument(
        "-i", "--input",
        required=True,
        help="Input alignment file (.bam or .sam)"
    )
    parser.add_argument(
        "-r", "--read-output",
        required=True,
        help="Output TSV with per-read information"
    )
    parser.add_argument(
        "-g", "--gene-summary",
        required=True,
        help="Output TSV with per-gene summary counts"
    )
    parser.add_argument(
        "--min-mapq",
        type=int,
        default=0,
        help=(
            "Minimum mapping quality for a PRIMARY alignment to be considered "
            "classified (default: 0)."
        ),
    )
    parser.add_argument(
        "--sample-id",
        default=None,
        help=(
            "Sample ID to use as the column header in the gene summary. "
            "If not provided, inferred from the BAM/SAM filename."
        ),
    )
    parser.add_argument(
        "--count-mode",
        choices=["read_end", "alignment"],
        default="read_end",
        help=(
            "How to summarize counts per gene. "
            "'read_end' (default): each read-end contributes exactly 1 count "
            "to its best primary alignment gene (highest MAPQ). "
            "'alignment': each individual alignment record that passes filters "
            "contributes 1 count (primary + secondary + supplementary)."
        ),
    )
    return parser.parse_args()


def guess_mode(path: str) -> str:
    if path.lower().endswith(".bam"):
        return "rb"
    return "r"


def open_alignment(path: str) -> pysam.AlignmentFile:
    mode = guess_mode(path)
    try:
        return pysam.AlignmentFile(path, mode)
    except OSError as e:
        sys.exit(f"[ERROR] Could not open {path}: {e}")


def sample_id_from_path(path: str) -> str:
    base = os.path.basename(path)
    for ext in [".bam", ".sam"]:
        if base.endswith(ext):
            return base[: -len(ext)]
    return base


def read_end_id(read: pysam.AlignedSegment) -> str:
    qname = read.query_name
    if read.is_paired:
        if read.is_read1:
            return f"{qname}/1"
        if read.is_read2:
            return f"{qname}/2"
    return qname


def aggregate_per_read_and_alignment_counts(
    aln: pysam.AlignmentFile,
    min_mapq: int = 0,
) -> Tuple[Dict[str, Dict[str, Any]], Counter]:
    """Iterate over all alignments and build:
    - per_read: dict keyed by read_end_id with:
        * 'primary_genes': dict gene -> list of MAPQs for primary alignments passing filter
        * 'secondary_genes': dict gene -> list of MAPQs for secondary alignments
        * 'supplementary_genes': dict gene -> list of MAPQs for supplementary alignments
    - align_counts: Counter(gene -> number of alignment records) using all
      mapped alignments (primary + secondary + supplementary) with MAPQ >= min_mapq.
    """
    per_read: Dict[str, Dict[str, Any]] = {}
    align_counts: Counter = Counter()

    for read in aln.fetch(until_eof=True):
        rid = read_end_id(read)

        if rid not in per_read:
            per_read[rid] = {
                "primary_genes": {},       # gene -> list[MAPQ]
                "secondary_genes": {},     # gene -> list[MAPQ]
                "supplementary_genes": {}, # gene -> list[MAPQ]
            }

        # Skip unmapped alignments entirely
        if read.is_unmapped:
            continue

        ref_name = aln.get_reference_name(read.reference_id)

        # Alignment-level counts: count all mapped alignments above MAPQ threshold
        if read.mapping_quality >= min_mapq:
            align_counts[ref_name] += 1

        # Secondary alignments
        if read.is_secondary:
            gene_map = per_read[rid]["secondary_genes"]
            gene_map.setdefault(ref_name, []).append(read.mapping_quality)
            continue

        # Supplementary alignments
        if read.is_supplementary:
            gene_map = per_read[rid]["supplementary_genes"]
            gene_map.setdefault(ref_name, []).append(read.mapping_quality)
            continue

        # Primary alignment: only classify if MAPQ >= threshold
        if read.mapping_quality >= min_mapq:
            gene_map = per_read[rid]["primary_genes"]
            gene_map.setdefault(ref_name, []).append(read.mapping_quality)

    return per_read, align_counts


def format_gene_mapq_field(gene_mapq_dict: Dict[str, list]) -> str:
    """Format a gene -> [mapq, mapq, ...] dict as gene(mapq,mapq)/gene(mapq)"""
    if not gene_mapq_dict:
        return "-"
    
    parts = []
    for gene in sorted(gene_mapq_dict.keys()):
        mqs = gene_mapq_dict[gene]
        mq_str = ",".join(str(q) for q in mqs)
        parts.append(f"{gene}({mq_str})")
    return "/".join(parts)


def count_alignments(gene_mapq_dict: Dict[str, list]) -> int:
    """Count total number of alignments in a gene -> [mapq, ...] dict"""
    return sum(len(mqs) for mqs in gene_mapq_dict.values())


def write_read_output(
    per_read: Dict[str, Dict[str, Any]],
    out_path: str,
) -> None:
    """Write per-read TSV with columns:
        read_id
        num_all_alignments
        num_primary_secondary
        num_primary_alignments
        num_secondary_alignments
        num_supplementary_alignments
        primary_genes
        secondary_genes
        supplementary_genes
    
    Gene columns are formatted as: gene(mapq,mapq)/gene(mapq)
    If no alignments in a category, the field is '-'.
    """
    with open(out_path, "w") as out:
        out.write(
            "read_id\t"
            "num_all_alignments\t"
            "num_primary_secondary\t"
            "num_primary_alignments\t"
            "num_secondary_alignments\t"
            "num_supplementary_alignments\t"
            "primary_genes\t"
            "secondary_genes\t"
            "supplementary_genes\n"
        )

        for rid, info in per_read.items():
            primary_genes = info["primary_genes"]
            secondary_genes = info["secondary_genes"]
            supplementary_genes = info["supplementary_genes"]

            num_primary = count_alignments(primary_genes)
            num_secondary = count_alignments(secondary_genes)
            num_supplementary = count_alignments(supplementary_genes)
            num_primary_secondary = num_primary + num_secondary
            num_all = num_primary + num_secondary + num_supplementary

            primary_field = format_gene_mapq_field(primary_genes)
            secondary_field = format_gene_mapq_field(secondary_genes)
            supplementary_field = format_gene_mapq_field(supplementary_genes)

            out.write(
                f"{rid}\t"
                f"{num_all}\t"
                f"{num_primary_secondary}\t"
                f"{num_primary}\t"
                f"{num_secondary}\t"
                f"{num_supplementary}\t"
                f"{primary_field}\t"
                f"{secondary_field}\t"
                f"{supplementary_field}\n"
            )


def summarize_genes_read_end(
    per_read: Dict[str, Dict[str, Any]]
) -> Counter:
    """Build a Counter of how many read-ends are classified to each gene.
    
    For each read-end:
      - If it has one or more primary alignments, pick the one with the highest MAPQ
      - That gene gets +1 count
      - If there's a tie, pick alphabetically first gene (deterministic)
      - Each read-end contributes AT MOST 1 count total
    
    Result: Counter(gene_accession -> count_of_read_ends)
    """
    counts = Counter()

    for info in per_read.values():
        primary_genes = info["primary_genes"]  # dict: gene -> [mapq, ...]

        if not primary_genes:
            continue

        # Find the best primary alignment: highest MAPQ, then alphabetically first gene
        best_gene = None
        best_mapq = -1

        for gene, mapqs in primary_genes.items():
            max_mapq_for_gene = max(mapqs)
            if max_mapq_for_gene > best_mapq or (max_mapq_for_gene == best_mapq and (best_gene is None or gene < best_gene)):
                best_gene = gene
                best_mapq = max_mapq_for_gene

        if best_gene is not None:
            counts[best_gene] += 1

    return counts


def write_gene_summary(
    gene_counts: Counter,
    sample_id: str,
    out_path: str,
) -> None:
    with open(out_path, "w") as out:
        out.write(f"gene_accession\t{sample_id}\n")
        for gene, count in sorted(gene_counts.items()):
            out.write(f"{gene}\t{count}\n")


def main():
    args = parse_args()

    if args.sample_id is not None:
        sample_id = args.sample_id
    else:
        sample_id = sample_id_from_path(args.input)

    aln = open_alignment(args.input)

    per_read, align_counts = aggregate_per_read_and_alignment_counts(
        aln,
        min_mapq=args.min_mapq,
    )

    write_read_output(per_read, args.read_output)

    if args.count_mode == "read_end":
        gene_counts = summarize_genes_read_end(per_read)
    else:
        gene_counts = align_counts

    write_gene_summary(gene_counts, sample_id, args.gene_summary)

    print(f"[INFO] Wrote per-read file to: {args.read_output}")
    print(f"[INFO] Wrote per-gene summary to: {args.gene_summary}")
    print(f"[INFO] Sample ID used: {sample_id}")
    print(f"[INFO] Count mode: {args.count_mode}")


if __name__ == "__main__":
    main()
