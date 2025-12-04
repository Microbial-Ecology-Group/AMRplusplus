#!/usr/bin/env python3
import argparse
import os
import sys
from collections import Counter, defaultdict
from typing import Dict, Any, Tuple, Set, Optional

import pysam

## Example command
#    python3 alignment_analyzer.py \
#        -i "$bam" \
#        -r "counts/${sample_id}_per_read.tsv" \
#        -g "counts/${sample_id}_gene_summary.tsv" \
#        --sample-id "$sample_id" \
#        --min-mapq 0 \
#        --count-mode alignment \
#        --min-gene-fraction 0.8 \
#        --include-supplementary \
#        --cigar-aware-coverage \
#        --coverage-output "counts/${sample_id}_coverage_stats.tsv"

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
            "Minimum mapping quality for an alignment to be considered "
            "(default: 0)."
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
            "contributes 1 count (primary + secondary, optionally + supplementary)."
        ),
    )
    parser.add_argument(
        "--include-supplementary",
        action="store_true",
        default=False,
        help=(
            "Include supplementary alignments in counts. "
            "By default, only primary and secondary alignments are counted. "
            "This flag adds supplementary (chimeric/split) alignments to the counts."
        ),
    )
    parser.add_argument(
        "--min-gene-fraction",
        type=float,
        default=0.0,
        help=(
            "Minimum fraction of gene length that must be covered by at least one "
            "alignment for that gene to be counted (0.0-1.0, default: 0.0 = no filter). "
            "E.g., 0.8 means 80%% of the gene must be covered. "
            "Applies to both read_end and alignment count modes."
        ),
    )
    parser.add_argument(
        "--cigar-aware-coverage",
        action="store_true",
        default=False,
        help=(
            "Use CIGAR-aware coverage calculation. "
            "Only counts positions with M/=/X operations as covered. "
            "Excludes deletions (D) and skipped regions (N) from coverage. "
            "More accurate but slightly slower. Recommended for RNA-seq or "
            "alignments with large indels."
        ),
    )
    parser.add_argument(
        "--coverage-output",
        default=None,
        help=(
            "Optional output TSV with per-gene coverage statistics. "
            "If not provided, coverage stats are not written."
        ),
    )
    parser.add_argument(
        "--force",
        action="store_true",
        default=False,
        help=(
            "Force overwrite of existing output files. "
            "By default, if output files exist, the sample is skipped."
        ),
    )
    return parser.parse_args()


def check_outputs_exist(args: argparse.Namespace) -> bool:
    """Check if all output files already exist."""
    outputs_to_check = [args.read_output, args.gene_summary]
    if args.coverage_output:
        outputs_to_check.append(args.coverage_output)
    
    return all(os.path.exists(f) for f in outputs_to_check)


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


def get_reference_lengths(aln: pysam.AlignmentFile) -> Dict[str, int]:
    """Extract reference (gene) lengths from the BAM/SAM header."""
    ref_lengths = {}
    for ref_name, ref_length in zip(aln.references, aln.lengths):
        ref_lengths[ref_name] = ref_length
    return ref_lengths


def get_covered_positions_cigar_aware(read: pysam.AlignedSegment) -> Set[int]:
    """Get reference positions covered by the read using CIGAR operations.
    
    Only counts positions where the read actually aligns (M, =, X operations).
    Does not count deletions (D) or skipped regions (N).
    
    CIGAR operations:
        M (0): alignment match (can be sequence match or mismatch)
        I (1): insertion to reference (consumes query, not reference)
        D (2): deletion from reference (consumes reference, not query)
        N (3): skipped region from reference (e.g., intron)
        S (4): soft clipping (consumes query, not reference)
        H (5): hard clipping (consumes neither)
        P (6): padding (consumes neither)
        = (7): sequence match (consumes both)
        X (8): sequence mismatch (consumes both)
    """
    covered = set()
    
    if read.cigartuples is None:
        return covered
    
    ref_pos = read.reference_start
    
    for op, length in read.cigartuples:
        if op in (0, 7, 8):  # M, =, X - actual alignment
            for i in range(length):
                covered.add(ref_pos + i)
            ref_pos += length
        elif op in (2, 3):  # D, N - consumes reference but not covered
            ref_pos += length
        # I (1), S (4), H (5), P (6) - don't consume reference
    
    return covered


def get_covered_positions_simple(read: pysam.AlignedSegment) -> Set[int]:
    """Get reference positions covered using simple start-end range.
    
    Counts all positions between reference_start and reference_end as covered.
    This includes deletions and skipped regions.
    """
    covered = set()
    
    if read.reference_start is not None and read.reference_end is not None:
        for pos in range(read.reference_start, read.reference_end):
            covered.add(pos)
    
    return covered


def aggregate_per_read_and_alignment_counts(
    aln: pysam.AlignmentFile,
    min_mapq: int = 0,
    include_supplementary: bool = False,
    cigar_aware_coverage: bool = False,
) -> Tuple[Dict[str, Dict[str, Any]], Counter, Counter, Dict[str, Set[int]]]:
    """Iterate over all alignments and build:
    - per_read: dict keyed by read_end_id with:
        * 'primary_genes': dict gene -> list of MAPQs for primary alignments passing filter
        * 'secondary_genes': dict gene -> list of MAPQs for secondary alignments
        * 'supplementary_genes': dict gene -> list of MAPQs for supplementary alignments
    - align_counts_with_supp: Counter(gene -> count) for primary + secondary + supplementary
    - align_counts_no_supp: Counter(gene -> count) for primary + secondary only
    - gene_coverage: dict gene -> set of covered positions
    """
    per_read: Dict[str, Dict[str, Any]] = {}
    align_counts_with_supp: Counter = Counter()
    align_counts_no_supp: Counter = Counter()
    gene_coverage: Dict[str, Set[int]] = defaultdict(set)

    # Choose coverage function
    get_covered_positions = (
        get_covered_positions_cigar_aware if cigar_aware_coverage 
        else get_covered_positions_simple
    )

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

        # Track coverage using chosen method
        covered_positions = get_covered_positions(read)
        gene_coverage[ref_name].update(covered_positions)

        # Check MAPQ threshold
        passes_mapq = read.mapping_quality >= min_mapq

        # Secondary alignments
        if read.is_secondary:
            gene_map = per_read[rid]["secondary_genes"]
            gene_map.setdefault(ref_name, []).append(read.mapping_quality)
            if passes_mapq:
                align_counts_with_supp[ref_name] += 1
                align_counts_no_supp[ref_name] += 1
            continue

        # Supplementary alignments
        if read.is_supplementary:
            gene_map = per_read[rid]["supplementary_genes"]
            gene_map.setdefault(ref_name, []).append(read.mapping_quality)
            if passes_mapq:
                align_counts_with_supp[ref_name] += 1
                # NOT added to align_counts_no_supp
            continue

        # Primary alignment
        if passes_mapq:
            gene_map = per_read[rid]["primary_genes"]
            gene_map.setdefault(ref_name, []).append(read.mapping_quality)
            align_counts_with_supp[ref_name] += 1
            align_counts_no_supp[ref_name] += 1

    return per_read, align_counts_with_supp, align_counts_no_supp, gene_coverage


def calculate_gene_fractions(
    gene_coverage: Dict[str, Set[int]],
    ref_lengths: Dict[str, int]
) -> Dict[str, float]:
    """Calculate the fraction of each gene covered by alignments."""
    gene_fractions = {}
    for gene, covered_positions in gene_coverage.items():
        gene_length = ref_lengths.get(gene, 0)
        if gene_length > 0:
            gene_fractions[gene] = len(covered_positions) / gene_length
        else:
            gene_fractions[gene] = 0.0
    return gene_fractions


def get_passing_genes(
    gene_fractions: Dict[str, float],
    min_fraction: float
) -> Optional[Set[str]]:
    """Return set of genes that meet the minimum coverage fraction threshold."""
    if min_fraction <= 0.0:
        return None  # No filtering
    return {gene for gene, frac in gene_fractions.items() if frac >= min_fraction}


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
    per_read: Dict[str, Dict[str, Any]],
    passing_genes: Optional[Set[str]] = None
) -> Counter:
    """Build a Counter of how many read-ends are classified to each gene.
    
    For each read-end:
      - If it has one or more primary alignments, pick the one with the highest MAPQ
      - That gene gets +1 count
      - If there's a tie, pick alphabetically first gene (deterministic)
      - Each read-end contributes AT MOST 1 count total
      - If passing_genes is provided, only count genes in that set
    
    Result: Counter(gene_accession -> count_of_read_ends)
    """
    counts = Counter()

    for info in per_read.values():
        primary_genes = info["primary_genes"]  # dict: gene -> [mapq, ...]

        if not primary_genes:
            continue

        # Filter to only passing genes if threshold is set
        if passing_genes is not None:
            primary_genes = {g: m for g, m in primary_genes.items() if g in passing_genes}
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


def filter_align_counts(
    align_counts: Counter,
    passing_genes: Optional[Set[str]] = None
) -> Counter:
    """Filter alignment counts to only include passing genes."""
    if passing_genes is None:
        return align_counts
    return Counter({g: c for g, c in align_counts.items() if g in passing_genes})


def write_gene_summary(
    gene_counts: Counter,
    sample_id: str,
    out_path: str,
) -> None:
    with open(out_path, "w") as out:
        out.write(f"gene_accession\t{sample_id}\n")
        for gene, count in sorted(gene_counts.items()):
            out.write(f"{gene}\t{count}\n")


def write_coverage_output(
    gene_fractions: Dict[str, float],
    ref_lengths: Dict[str, int],
    gene_coverage: Dict[str, Set[int]],
    passing_genes: Optional[Set[str]],
    min_fraction: float,
    out_path: str,
) -> None:
    """Write per-gene coverage statistics."""
    with open(out_path, "w") as out:
        out.write(
            "gene_accession\t"
            "gene_length\t"
            "covered_bases\t"
            "coverage_fraction\t"
            "coverage_percent\t"
            "passes_threshold\n"
        )
        
        # Include all genes that have any coverage or are in reference
        all_genes = set(ref_lengths.keys()) | set(gene_coverage.keys())
        
        for gene in sorted(all_genes):
            gene_length = ref_lengths.get(gene, 0)
            covered_bases = len(gene_coverage.get(gene, set()))
            fraction = gene_fractions.get(gene, 0.0)
            percent = fraction * 100
            passes = "yes" if (passing_genes is None or gene in passing_genes) else "no"
            
            out.write(
                f"{gene}\t"
                f"{gene_length}\t"
                f"{covered_bases}\t"
                f"{fraction:.4f}\t"
                f"{percent:.2f}\t"
                f"{passes}\n"
            )


def main():
    args = parse_args()

    if args.sample_id is not None:
        sample_id = args.sample_id
    else:
        sample_id = sample_id_from_path(args.input)

    # Check if outputs already exist
    if not args.force and check_outputs_exist(args):
        print(f"[SKIP] Sample '{sample_id}' already processed. Output files exist. Use --force to overwrite.")
        sys.exit(0)

    # Validate min-gene-fraction
    if not 0.0 <= args.min_gene_fraction <= 1.0:
        sys.exit(f"[ERROR] --min-gene-fraction must be between 0.0 and 1.0, got {args.min_gene_fraction}")

    aln = open_alignment(args.input)

    # Get reference lengths from header
    ref_lengths = get_reference_lengths(aln)
    print(f"[INFO] Found {len(ref_lengths)} references in BAM header")

    per_read, align_counts_with_supp, align_counts_no_supp, gene_coverage = aggregate_per_read_and_alignment_counts(
        aln,
        min_mapq=args.min_mapq,
        include_supplementary=args.include_supplementary,
        cigar_aware_coverage=args.cigar_aware_coverage,
    )

    # Calculate gene coverage fractions
    gene_fractions = calculate_gene_fractions(gene_coverage, ref_lengths)

    # Determine which genes pass the threshold
    passing_genes = get_passing_genes(gene_fractions, args.min_gene_fraction)

    if passing_genes is not None:
        num_passing = len(passing_genes)
        num_total = len(gene_coverage)
        print(f"[INFO] Gene fraction filter: {num_passing}/{num_total} genes pass >= {args.min_gene_fraction:.1%} coverage")

    # Write per-read output (unfiltered - shows all alignments)
    write_read_output(per_read, args.read_output)

    # Write coverage output if requested
    if args.coverage_output:
        write_coverage_output(
            gene_fractions,
            ref_lengths,
            gene_coverage,
            passing_genes,
            args.min_gene_fraction,
            args.coverage_output
        )
        print(f"[INFO] Wrote coverage stats to: {args.coverage_output}")

    # Apply gene fraction filter to counts based on mode
    if args.count_mode == "read_end":
        gene_counts = summarize_genes_read_end(per_read, passing_genes)
    else:  # alignment mode
        if args.include_supplementary:
            gene_counts = filter_align_counts(align_counts_with_supp, passing_genes)
        else:
            gene_counts = filter_align_counts(align_counts_no_supp, passing_genes)

    write_gene_summary(gene_counts, sample_id, args.gene_summary)

    print(f"[INFO] Wrote per-read file to: {args.read_output}")
    print(f"[INFO] Wrote per-gene summary to: {args.gene_summary}")
    print(f"[INFO] Sample ID used: {sample_id}")
    print(f"[INFO] Count mode: {args.count_mode}")
    print(f"[INFO] Include supplementary: {args.include_supplementary}")
    print(f"[INFO] CIGAR-aware coverage: {args.cigar_aware_coverage}")
    print(f"[INFO] Min gene fraction: {args.min_gene_fraction:.1%}")


if __name__ == "__main__":
    main()
