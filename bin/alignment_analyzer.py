#!/usr/bin/env python3
"""
alignment_analyzer.py
════════════════════════════════════════════════════════════════════════════════

For each read (read-end) in a BAM/SAM, output the gene(s) it maps to via primary
alignments and the MAPQ, plus a summary of counts per gene.

CHANGES IN THIS VERSION
────────────────────────
1. --min-query-coverage flag (0.0-1.0, default 0.0). Read-level counterpart to
   --min-gene-fraction: an alignment must cover this fraction of the READ's
   length (aligned query bases / total read length) to be counted.

2. Gene-length breadth coverage (--min-gene-fraction) now respects the same
   quality filters as alignment counting by default. --coverage-ignore-filters
   restores the old unfiltered behaviour.

3. --coverage-output gains per-gene query-coverage depth stats (mean/median/n)
   alongside the existing breadth-based coverage_fraction, plus a "meg_id"
   join column (gene_accession split on the first '|').

4. NEW — retained-read tracking. Every run now reports, and optionally writes
   to a small stats file, how many primary alignments existed BEFORE any
   --min-mapq/--min-query-coverage filtering vs how many passed (and the %
   retained). Previously this number wasn't tracked or surfaced anywhere —
   the script silently filtered without reporting what fraction of the
   original classified reads survived.

## Example command
#    python3 alignment_analyzer.py \
#        -i "$bam" \
#        -r "counts/${sample_id}_per_read.tsv" \
#        -g "counts/${sample_id}_gene_summary.tsv" \
#        --sample-id "$sample_id" \
#        --min-mapq 0 \
#        --count-mode alignment \
#        --min-gene-fraction 0.8 \
#        --min-query-coverage 0.8 \
#        --include-supplementary \
#        --cigar-aware-coverage \
#        --coverage-output "counts/${sample_id}_coverage_stats.tsv"
"""

import argparse
import os
import statistics
import sys
from collections import Counter, defaultdict
from typing import Dict, Any, Tuple, Set, Optional, List

import pysam


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "For each read (read-end) in a BAM/SAM, output the gene(s) it maps to "
            "via primary alignments and the MAPQ, plus a summary of counts per gene."
        )
    )
    parser.add_argument("-i", "--input", required=True, help="Input alignment file (.bam or .sam)")
    parser.add_argument("-r", "--read-output", required=True, help="Output TSV with per-read information")
    parser.add_argument("-g", "--gene-summary", required=True, help="Output TSV with per-gene summary counts")
    parser.add_argument("--min-mapq", type=int, default=0,
                        help="Minimum mapping quality for an alignment to be considered (default: 0).")
    parser.add_argument(
        "--min-query-coverage", type=float, default=0.0,
        help=(
            "Minimum fraction of the READ that must be covered by the alignment "
            "(aligned query bases / total read length, 0.0-1.0, default: 0.0). "
            "Read-level counterpart to --min-gene-fraction: filters out "
            "alignments built from only a short partial overlap of the read."
        ),
    )
    parser.add_argument(
        "--coverage-ignore-filters", action="store_true", default=False,
        help=(
            "Compute gene-length breadth coverage from ALL alignments, ignoring "
            "--min-mapq/--min-query-coverage/--include-supplementary. Restores "
            "pre-this-version behaviour. By default breadth now only counts "
            "alignments that pass the same filters used for counting."
        ),
    )
    parser.add_argument("--sample-id", default=None,
                        help="Sample ID column header in the gene summary; inferred from filename if omitted.")
    parser.add_argument(
        "--count-mode", choices=["read_end", "alignment"], default="read_end",
        help=(
            "'read_end' (default): each read-end contributes exactly 1 count to "
            "its best primary alignment gene (highest MAPQ). "
            "'alignment': each individual passing alignment record contributes 1 count."
        ),
    )
    parser.add_argument("--include-supplementary", action="store_true", default=False,
                        help="Include supplementary alignments in counts (and, by default, in breadth coverage).")
    parser.add_argument(
        "--min-gene-fraction", type=float, default=0.0,
        help="Minimum fraction of gene length covered by passing alignments for that gene to count (0.0-1.0).",
    )
    parser.add_argument("--cigar-aware-coverage", action="store_true", default=False,
                        help="Use CIGAR-aware breadth coverage (excludes deletions/skips). Recommended.")
    parser.add_argument("--coverage-output", default=None,
                        help="Optional output TSV with per-gene breadth AND query-coverage depth statistics.")
    parser.add_argument(
        "--stats-output", default=None,
        help=(
            "Optional output TSV with overall retained-read statistics: total "
            "primary alignments seen before filtering, how many passed "
            "--min-mapq/--min-query-coverage, and the percent retained. "
            "If omitted but --coverage-output is set, defaults to "
            "'<coverage-output>.stats.tsv'."
        ),
    )
    parser.add_argument("--force", action="store_true", default=False, help="Force overwrite of existing output files.")
    return parser.parse_args()


def check_outputs_exist(args: argparse.Namespace) -> bool:
    outputs_to_check = [args.read_output, args.gene_summary]
    if args.coverage_output:
        outputs_to_check.append(args.coverage_output)
    return all(os.path.exists(f) for f in outputs_to_check)


def guess_mode(path: str) -> str:
    return "rb" if path.lower().endswith(".bam") else "r"


def open_alignment(path: str) -> pysam.AlignmentFile:
    try:
        return pysam.AlignmentFile(path, guess_mode(path))
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


def extract_meg_id(gene_accession: str) -> str:
    """First pipe-delimited field, e.g. 'MEG_4288' from 'MEG_4288|Drugs|...'."""
    return gene_accession.split('|')[0] if gene_accession else gene_accession


def get_reference_lengths(aln: pysam.AlignmentFile) -> Dict[str, int]:
    return dict(zip(aln.references, aln.lengths))


def get_covered_positions_cigar_aware(read: pysam.AlignedSegment) -> Set[int]:
    covered = set()
    if read.cigartuples is None:
        return covered
    ref_pos = read.reference_start
    for op, length in read.cigartuples:
        if op in (0, 7, 8):
            covered.update(range(ref_pos, ref_pos + length))
            ref_pos += length
        elif op in (2, 3):
            ref_pos += length
    return covered


def get_covered_positions_simple(read: pysam.AlignedSegment) -> Set[int]:
    if read.reference_start is None or read.reference_end is None:
        return set()
    return set(range(read.reference_start, read.reference_end))


def aggregate_per_read_and_alignment_counts(
    aln: pysam.AlignmentFile,
    min_mapq: int = 0,
    min_query_coverage: float = 0.0,
    include_supplementary: bool = False,
    cigar_aware_coverage: bool = False,
    coverage_ignore_filters: bool = False,
) -> Tuple[Dict[str, Dict[str, Any]], Counter, Counter, Dict[str, Set[int]], Dict[str, List[float]], int, int]:
    """
    Returns (per_read, align_counts_with_supp, align_counts_no_supp,
              gene_coverage, gene_query_coverages,
              total_primary_seen, total_primary_passing)

    total_primary_seen     = every mapped PRIMARY alignment encountered,
                              regardless of --min-mapq/--min-query-coverage.
    total_primary_passing  = how many of those passed both filters and were
                              actually counted as a classification.
    These two let main() report exactly what fraction of originally-mapped
    reads survive filtering — previously untracked.
    """
    per_read: Dict[str, Dict[str, Any]] = {}
    align_counts_with_supp: Counter = Counter()
    align_counts_no_supp: Counter = Counter()
    gene_coverage: Dict[str, Set[int]] = defaultdict(set)
    gene_query_coverages: Dict[str, List[float]] = defaultdict(list)
    total_primary_seen = 0
    total_primary_passing = 0

    get_covered_positions = (
        get_covered_positions_cigar_aware if cigar_aware_coverage
        else get_covered_positions_simple
    )

    min_qcov_pct = min_query_coverage * 100.0

    for read in aln.fetch(until_eof=True):
        rid = read_end_id(read)

        if rid not in per_read:
            per_read[rid] = {"primary_genes": {}, "secondary_genes": {}, "supplementary_genes": {}}

        if read.is_unmapped:
            continue

        ref_name = aln.get_reference_name(read.reference_id)

        qlen = read.query_length or 0
        aln_len = read.query_alignment_length or 0
        qcov_pct = round(aln_len / qlen * 100.0, 2) if qlen else 0.0

        passes_mapq = read.mapping_quality >= min_mapq
        passes_qcov = qcov_pct >= min_qcov_pct
        passes_filter = passes_mapq and passes_qcov

        include_for_breadth = (
            coverage_ignore_filters
            or (passes_filter and (include_supplementary or not read.is_supplementary))
        )
        if include_for_breadth:
            gene_coverage[ref_name].update(get_covered_positions(read))

        if read.is_secondary:
            gene_map = per_read[rid]["secondary_genes"]
            gene_map.setdefault(ref_name, []).append((read.mapping_quality, qcov_pct))
            if passes_filter:
                align_counts_with_supp[ref_name] += 1
                align_counts_no_supp[ref_name] += 1
            continue

        if read.is_supplementary:
            gene_map = per_read[rid]["supplementary_genes"]
            gene_map.setdefault(ref_name, []).append((read.mapping_quality, qcov_pct))
            if passes_filter:
                align_counts_with_supp[ref_name] += 1
            continue

        # ── primary alignment ──────────────────────────────────────────────
        total_primary_seen += 1
        if passes_filter:
            total_primary_passing += 1
            gene_map = per_read[rid]["primary_genes"]
            gene_map.setdefault(ref_name, []).append((read.mapping_quality, qcov_pct))
            align_counts_with_supp[ref_name] += 1
            align_counts_no_supp[ref_name] += 1
            gene_query_coverages[ref_name].append(qcov_pct)

    return (per_read, align_counts_with_supp, align_counts_no_supp,
            gene_coverage, gene_query_coverages,
            total_primary_seen, total_primary_passing)


def calculate_gene_fractions(gene_coverage: Dict[str, Set[int]], ref_lengths: Dict[str, int]) -> Dict[str, float]:
    fractions = {}
    for gene, covered in gene_coverage.items():
        length = ref_lengths.get(gene, 0)
        fractions[gene] = (len(covered) / length) if length > 0 else 0.0
    return fractions


def calculate_query_coverage_stats(gene_query_coverages: Dict[str, List[float]]) -> Dict[str, Tuple[float, float, int]]:
    stats = {}
    for gene, values in gene_query_coverages.items():
        if values:
            stats[gene] = (round(statistics.mean(values), 2), round(statistics.median(values), 2), len(values))
        else:
            stats[gene] = (0.0, 0.0, 0)
    return stats


def get_passing_genes(gene_fractions: Dict[str, float], min_fraction: float) -> Optional[Set[str]]:
    if min_fraction <= 0.0:
        return None
    return {gene for gene, frac in gene_fractions.items() if frac >= min_fraction}


def format_gene_alignment_field(gene_mapq_qcov: Dict[str, List[Tuple[int, float]]]) -> str:
    if not gene_mapq_qcov:
        return "-"
    parts = []
    for gene in sorted(gene_mapq_qcov.keys()):
        entry_str = ",".join(f"{mq}:{qc:.1f}" for mq, qc in gene_mapq_qcov[gene])
        parts.append(f"{gene}({entry_str})")
    return "/".join(parts)


def count_alignments(gene_mapq_qcov: Dict[str, List[Tuple[int, float]]]) -> int:
    return sum(len(entries) for entries in gene_mapq_qcov.values())


def write_read_output(per_read: Dict[str, Dict[str, Any]], out_path: str) -> None:
    with open(out_path, "w") as out:
        out.write(
            "read_id\tnum_all_alignments\tnum_primary_secondary\t"
            "num_primary_alignments\tnum_secondary_alignments\tnum_supplementary_alignments\t"
            "primary_genes\tsecondary_genes\tsupplementary_genes\n"
        )
        for rid, info in per_read.items():
            primary, secondary, supplementary = info["primary_genes"], info["secondary_genes"], info["supplementary_genes"]
            num_primary, num_secondary, num_supplementary = (
                count_alignments(primary), count_alignments(secondary), count_alignments(supplementary)
            )
            out.write(
                f"{rid}\t{num_primary + num_secondary + num_supplementary}\t"
                f"{num_primary + num_secondary}\t{num_primary}\t{num_secondary}\t{num_supplementary}\t"
                f"{format_gene_alignment_field(primary)}\t"
                f"{format_gene_alignment_field(secondary)}\t"
                f"{format_gene_alignment_field(supplementary)}\n"
            )


def summarize_genes_read_end(per_read: Dict[str, Dict[str, Any]], passing_genes: Optional[Set[str]] = None) -> Counter:
    counts = Counter()
    for info in per_read.values():
        primary_genes = info["primary_genes"]
        if not primary_genes:
            continue
        if passing_genes is not None:
            primary_genes = {g: v for g, v in primary_genes.items() if g in passing_genes}
            if not primary_genes:
                continue
        best_gene, best_mapq = None, -1
        for gene, entries in primary_genes.items():
            max_mapq_for_gene = max(mq for mq, _qc in entries)
            if max_mapq_for_gene > best_mapq or (max_mapq_for_gene == best_mapq and (best_gene is None or gene < best_gene)):
                best_gene, best_mapq = gene, max_mapq_for_gene
        if best_gene is not None:
            counts[best_gene] += 1
    return counts


def filter_align_counts(align_counts: Counter, passing_genes: Optional[Set[str]] = None) -> Counter:
    if passing_genes is None:
        return align_counts
    return Counter({g: c for g, c in align_counts.items() if g in passing_genes})


def write_gene_summary(gene_counts: Counter, sample_id: str, out_path: str) -> None:
    with open(out_path, "w") as out:
        out.write(f"gene_accession\tmeg_id\t{sample_id}\n")
        for gene, count in sorted(gene_counts.items()):
            out.write(f"{gene}\t{extract_meg_id(gene)}\t{count}\n")


def write_coverage_output(
    gene_fractions: Dict[str, float], ref_lengths: Dict[str, int], gene_coverage: Dict[str, Set[int]],
    query_coverage_stats: Dict[str, Tuple[float, float, int]], passing_genes: Optional[Set[str]],
    min_fraction: float, out_path: str,
) -> None:
    with open(out_path, "w") as out:
        out.write(
            "gene_accession\tmeg_id\tgene_length\tcovered_bases\t"
            "coverage_fraction\tcoverage_percent\tpasses_threshold\t"
            "n_primary_alignments\tmean_query_coverage_pct\tmedian_query_coverage_pct\n"
        )
        all_genes = set(ref_lengths.keys()) | set(gene_coverage.keys()) | set(query_coverage_stats.keys())
        for gene in sorted(all_genes):
            gene_length = ref_lengths.get(gene, 0)
            covered_bases = len(gene_coverage.get(gene, set()))
            fraction = gene_fractions.get(gene, 0.0)
            percent = fraction * 100
            passes = "yes" if (passing_genes is None or gene in passing_genes) else "no"
            mean_qc, median_qc, n_aln = query_coverage_stats.get(gene, (0.0, 0.0, 0))
            out.write(
                f"{gene}\t{extract_meg_id(gene)}\t{gene_length}\t{covered_bases}\t"
                f"{fraction:.4f}\t{percent:.2f}\t{passes}\t{n_aln}\t{mean_qc:.2f}\t{median_qc:.2f}\n"
            )


def write_stats_output(
    out_path: str, sample_id: str, total_primary_seen: int, total_primary_passing: int,
    min_mapq: int, min_query_coverage: float, n_genes_baseline: int, n_genes_passing: Optional[int],
    min_gene_fraction: float,
) -> None:
    """
    Overall retained-read diagnostic — answers "how many reads were included
    as classifications, and what % of the total classified reads survived
    filtering?" which was previously not tracked anywhere in this script.
    """
    pct_retained = round(total_primary_passing / total_primary_seen * 100, 2) if total_primary_seen else 0.0
    pct_genes_passing = (
        round(n_genes_passing / n_genes_baseline * 100, 2)
        if (n_genes_passing is not None and n_genes_baseline) else None
    )
    with open(out_path, "w") as out:
        out.write("metric\tvalue\n")
        out.write(f"sample_id\t{sample_id}\n")
        out.write(f"min_mapq\t{min_mapq}\n")
        out.write(f"min_query_coverage\t{min_query_coverage}\n")
        out.write(f"min_gene_fraction\t{min_gene_fraction}\n")
        out.write(f"total_primary_alignments_seen\t{total_primary_seen}\n")
        out.write(f"total_primary_alignments_passing\t{total_primary_passing}\n")
        out.write(f"pct_primary_alignments_retained\t{pct_retained}\n")
        out.write(f"n_genes_detected_baseline\t{n_genes_baseline}\n")
        out.write(f"n_genes_passing_fraction_threshold\t{n_genes_passing if n_genes_passing is not None else 'NA'}\n")
        out.write(f"pct_genes_passing_of_baseline\t{pct_genes_passing if pct_genes_passing is not None else 'NA'}\n")


def main():
    args = parse_args()
    sample_id = args.sample_id if args.sample_id is not None else sample_id_from_path(args.input)

    if not args.force and check_outputs_exist(args):
        print(f"[SKIP] Sample '{sample_id}' already processed. Use --force to overwrite.")
        sys.exit(0)

    if not 0.0 <= args.min_gene_fraction <= 1.0:
        sys.exit(f"[ERROR] --min-gene-fraction must be 0.0-1.0, got {args.min_gene_fraction}")
    if not 0.0 <= args.min_query_coverage <= 1.0:
        sys.exit(f"[ERROR] --min-query-coverage must be 0.0-1.0, got {args.min_query_coverage}")

    aln = open_alignment(args.input)
    ref_lengths = get_reference_lengths(aln)
    print(f"[INFO] Found {len(ref_lengths)} references in BAM header")

    (per_read, align_counts_with_supp, align_counts_no_supp, gene_coverage,
     gene_query_coverages, total_primary_seen, total_primary_passing) = (
        aggregate_per_read_and_alignment_counts(
            aln, min_mapq=args.min_mapq, min_query_coverage=args.min_query_coverage,
            include_supplementary=args.include_supplementary,
            cigar_aware_coverage=args.cigar_aware_coverage,
            coverage_ignore_filters=args.coverage_ignore_filters,
        )
    )

    gene_fractions = calculate_gene_fractions(gene_coverage, ref_lengths)
    query_coverage_stats = calculate_query_coverage_stats(gene_query_coverages)
    passing_genes = get_passing_genes(gene_fractions, args.min_gene_fraction)
    n_genes_baseline = len(gene_coverage)

    pct_retained = round(total_primary_passing / total_primary_seen * 100, 2) if total_primary_seen else 0.0
    print(f"[INFO] Primary alignments: {total_primary_passing:,} / {total_primary_seen:,} "
          f"passed --min-mapq {args.min_mapq} + --min-query-coverage {args.min_query_coverage:.0%} "
          f"({pct_retained:.2f}% retained)")

    if passing_genes is not None:
        print(f"[INFO] Gene fraction filter: {len(passing_genes)}/{n_genes_baseline} "
              f"genes pass >= {args.min_gene_fraction:.1%} coverage")

    write_read_output(per_read, args.read_output)

    if args.coverage_output:
        write_coverage_output(
            gene_fractions, ref_lengths, gene_coverage, query_coverage_stats,
            passing_genes, args.min_gene_fraction, args.coverage_output,
        )
        print(f"[INFO] Wrote coverage stats (breadth + depth) to: {args.coverage_output}")

    stats_path = args.stats_output or (f"{args.coverage_output}.stats.tsv" if args.coverage_output else None)
    if stats_path:
        write_stats_output(
            stats_path, sample_id, total_primary_seen, total_primary_passing,
            args.min_mapq, args.min_query_coverage, n_genes_baseline,
            len(passing_genes) if passing_genes is not None else None,
            args.min_gene_fraction,
        )
        print(f"[INFO] Wrote retained-read stats to: {stats_path}")

    if args.count_mode == "read_end":
        gene_counts = summarize_genes_read_end(per_read, passing_genes)
    else:
        gene_counts = filter_align_counts(
            align_counts_with_supp if args.include_supplementary else align_counts_no_supp, passing_genes,
        )

    write_gene_summary(gene_counts, sample_id, args.gene_summary)

    print(f"[INFO] Wrote per-read file to: {args.read_output}")
    print(f"[INFO] Wrote per-gene summary to: {args.gene_summary}")
    print(f"[INFO] Sample ID: {sample_id} | Count mode: {args.count_mode} | "
          f"Min gene fraction: {args.min_gene_fraction:.1%} | "
          f"Min query coverage: {args.min_query_coverage:.1%}")


if __name__ == "__main__":
    main()