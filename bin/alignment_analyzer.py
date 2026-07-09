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
import re
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
        "--count-mode", choices=["read_end", "alignment", "fragment"], default="fragment",
        help=(
            "'read_end' (default): each read-end contributes exactly 1 count to "
            "its best primary alignment gene (highest MAPQ). "
            "'alignment': each individual passing alignment record contributes 1 count. "
            "'fragment': mates are resolved together into one physical fragment "
            "(base_read_id groups R1/R2 and strips FLASH suffixes), so a paired-end "
            "fragment and a FLASH-merged fragment are weighted IDENTICALLY (1 hit each). "
            "This is the counting unit used by coverage_threshold_sweep.py and fixes the "
            "old Combined-BAM problem where unmerged (paired) reads were counted twice "
            "relative to merged reads."
        ),
    )
    parser.add_argument(
        "--group-aware", dest="group_aware", action="store_true", default=True,
        help=(
            "(DEFAULT ON) Only meaningful with --count-mode fragment. When two mates of a "
            "fragment disagree on gene_accession, check whether the two accessions share the "
            "same MEGARes Group before treating it as real discordance. Same-Group disagreement "
            "(the multi-mapping/redundant-accession artifact from BWA tie-breaking among "
            "near-identical DB entries) is collapsed to ONE hit for the higher-match_qcov mate "
            "rather than double-counted. Group is parsed directly from the MEGARes "
            "pipe-delimited reference name (no annotation file needed). "
            "Use --no-group-aware to disable."
        ),
    )
    parser.add_argument(
        "--no-group-aware", dest="group_aware", action="store_false",
        help="Disable group-aware fragment resolution (mate disagreements count to each gene).",
    )
    parser.add_argument("--include-supplementary", action="store_true", default=False,
                        help="Include supplementary alignments in counts (and, by default, in breadth coverage).")
    parser.add_argument(
        "--min-gene-fraction", type=float, default=0.0,
        help="Minimum fraction of gene length covered by passing alignments for that gene to count (0.0-1.0).",
    )
    parser.add_argument("--edge-aware-qcov", dest="edge_aware_qcov", action="store_true", default=True,
        help=(
            "(DEFAULT ON) Compute query-coverage against a coverage-anchored read length: "
            "reconstruct full read length from the CIGAR (soft AND hard clips) so "
            "the aligner's clip choice doesn't change qcov, and discount clip bases "
            "that overhang the reference edge (a read running off the end of a short "
            "gene is not penalized for the overhanging portion). Only affects the "
            "qcov (read-axis) metric; gene breadth/fraction is unaffected. "
            "Use --no-edge-aware-qcov to disable."
        ),
    )
    parser.add_argument("--no-edge-aware-qcov", dest="edge_aware_qcov", action="store_false",
        help="Disable edge-aware qcov (use raw read.query_length as the denominator).",
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


def edge_aware_query_length(read: pysam.AlignedSegment, ref_length: int) -> int:
    """
    Coverage-anchored read length for query-coverage (qcov).

    Two problems this fixes vs a raw read.query_length:
      1. SOFT vs HARD clip asymmetry. read.query_length includes soft-clipped
         bases but NOT hard-clipped ones, so the same physical alignment scores a
         different qcov depending purely on whether the aligner soft- or
         hard-clipped. We reconstruct the full physical length from the CIGAR
         (aligned + BOTH clip types, ops 4=S and 5=H) so the choice is irrelevant.
      2. REFERENCE-EDGE OVERHANG. When a read runs past the end of a short gene,
         the overhanging bases COULD NOT have aligned (the reference simply ended).
         Penalizing qcov for them is wrong. We subtract the portion of each
         terminal clip that falls outside [0, ref_length) -- i.e. the bases with no
         reference position to map to -- from the denominator.

    Returns the number of read bases that actually HAD A CHANCE to align. qcov is
    then aln_len / this, so a read fully explained by the gene that merely runs off
    the gene's edge scores ~100% regardless of soft/hard clipping.

    Never returns less than query_alignment_length (qcov <= 100%).
    """
    if read.cigartuples is None:
        return read.query_alignment_length or 0

    aln_len = read.query_alignment_length or 0
    ops = read.cigartuples

    lead_clip = ops[0][1] if ops and ops[0][0] in (4, 5) else 0
    trail_clip = ops[-1][1] if len(ops) > 1 and ops[-1][0] in (4, 5) else 0

    total_physical = aln_len + lead_clip + trail_clip

    # bases of each terminal clip that overhang the reference edge (couldn't align)
    ref_start = read.reference_start if read.reference_start is not None else 0
    ref_end = read.reference_end if read.reference_end is not None else 0
    overhang_5 = max(0, lead_clip - ref_start)
    overhang_3 = max(0, trail_clip - (ref_length - ref_end)) if ref_length else 0

    eff = total_physical - overhang_5 - overhang_3
    return max(eff, aln_len)


def aggregate_per_read_and_alignment_counts(
    aln: pysam.AlignmentFile,
    min_mapq: int = 0,
    min_query_coverage: float = 0.0,
    include_supplementary: bool = False,
    cigar_aware_coverage: bool = False,
    coverage_ignore_filters: bool = False,
    edge_aware_qcov: bool = False,
    ref_lengths: Optional[Dict[str, int]] = None,
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

        aln_len = read.query_alignment_length or 0
        if edge_aware_qcov:
            ref_len = (ref_lengths or {}).get(ref_name, 0)
            qlen = edge_aware_query_length(read, ref_len)
        else:
            qlen = read.query_length or 0
        qcov_pct = round(aln_len / qlen * 100.0, 2) if qlen else 0.0

        # match_qcov_pct = fraction of the WHOLE read explained by genuine matches
        # (aligned bases minus edit distance), over the same (edge-aware if enabled)
        # read-length denominator as qcov. Used as the group-aware tie-break metric.
        # Reads with no NM tag fall back to qcov_pct so they still tie-break sensibly.
        try:
            nm = read.get_tag('NM')
            match_qcov_pct = round((aln_len - nm) / qlen * 100.0, 2) if qlen else 0.0
        except KeyError:
            match_qcov_pct = qcov_pct

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
            gene_map.setdefault(ref_name, []).append((read.mapping_quality, qcov_pct, match_qcov_pct))
            if passes_filter:
                align_counts_with_supp[ref_name] += 1
                align_counts_no_supp[ref_name] += 1
            continue

        if read.is_supplementary:
            gene_map = per_read[rid]["supplementary_genes"]
            gene_map.setdefault(ref_name, []).append((read.mapping_quality, qcov_pct, match_qcov_pct))
            if passes_filter:
                align_counts_with_supp[ref_name] += 1
            continue

        # ── primary alignment ──────────────────────────────────────────────
        total_primary_seen += 1
        if passes_filter:
            total_primary_passing += 1
            gene_map = per_read[rid]["primary_genes"]
            gene_map.setdefault(ref_name, []).append((read.mapping_quality, qcov_pct, match_qcov_pct))
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


# ── Fragment-level counting (mirrors coverage_threshold_sweep.resolve_fragment_hits) ──

# Same regexes coverage_threshold_sweep.py uses, so fragment identity is identical
# between the two tools: strip FLASH suffixes and a trailing /1 /2 .1 .2.
_RE_FLASH = re.compile(r'\.(extendedFrags|notCombined)[^/]*$')
_RE_PAIR  = re.compile(r'[/\.][12]$')


def get_base_read_id(read_id: str) -> str:
    """Collapse a read-end id to its physical-fragment id (matches the sweep)."""
    base = _RE_FLASH.sub('', read_id)
    return _RE_PAIR.sub('', base)


def parse_megares_group(ref_name: str) -> str:
    """
    Extract the Group field from a MEGARes pipe-delimited accession.
    Header layout: MEG_ID|Type|Class|Mechanism|Group[|SNP]. Returns the Group
    (5th field, index 4) or '' if the name isn't in the expected format.
    Mirrors coverage_threshold_sweep.parse_gene_reference so both tools derive
    Group the same way — directly from the reference name, no annotation CSV needed.
    """
    parts = ref_name.split('|')
    return parts[4] if len(parts) >= 5 else ''


def summarize_genes_fragment(
    per_read: Dict[str, Dict[str, Any]],
    passing_genes: Optional[Set[str]] = None,
    group_aware: bool = False,
) -> Tuple[Counter, Dict[str, int]]:
    """
    Fragment-level counting. First reduce each read-end to a single best gene,
    then group those best-calls by base_read_id and resolve mates together.

    "Best gene per read-end" and the group-aware disagreement tie-break both use
    match_qcov_pct (fraction of the read explained by genuine matches), consistent
    with coverage_threshold_sweep.resolve_fragment_hits. Ties broken alphabetically.

      both mates same gene            -> 1 hit
      only one mate present           -> 1 hit
      mates disagree, group_aware off -> 1 hit to EACH gene (legacy behaviour)
      mates disagree, group_aware on:
        same Group, different accession -> 1 hit to the higher-match_qcov mate
                                           (redundant-DB-entry artifact collapsed)
        different Groups                -> 1 hit to EACH gene (real cross-mechanism
                                           fragment; both mechanisms counted)

    Returns (hits, diagnostics).
    """
    # Step 1: best gene per read-end, scored by match_qcov (3rd tuple element)
    frag_groups: Dict[str, List[Tuple[str, float]]] = defaultdict(list)
    for rid, info in per_read.items():
        primary_genes = info["primary_genes"]
        if not primary_genes:
            continue
        if passing_genes is not None:
            primary_genes = {g: v for g, v in primary_genes.items() if g in passing_genes}
            if not primary_genes:
                continue
        best_gene, best_score = None, -1.0
        for gene, entries in primary_genes.items():
            # entries are (mapq, qcov_pct, match_qcov_pct); score on match_qcov
            max_mqcov = max(mqc for _mq, _qc, mqc in entries)
            if max_mqcov > best_score or (max_mqcov == best_score and (best_gene is None or gene < best_gene)):
                best_gene, best_score = gene, max_mqcov
        if best_gene is not None:
            frag_groups[get_base_read_id(rid)].append((best_gene, best_score))

    # Step 2: resolve mates within each fragment
    hits: Counter = Counter()
    diagnostics: Dict[str, int] = defaultdict(int)

    for base_id, entries in frag_groups.items():
        if len(entries) == 1:
            hits[entries[0][0]] += 1
            diagnostics['one_mate_only'] += 1
        elif len(entries) == 2:
            (g1, s1), (g2, s2) = entries
            if g1 == g2:
                hits[g1] += 1
                diagnostics['same_gene'] += 1
            elif group_aware:
                grp1, grp2 = parse_megares_group(g1), parse_megares_group(g2)
                if grp1 and grp1 == grp2:
                    # Same Group, different accession = redundant-DB-entry artifact.
                    # Collapse to ONE hit for the higher-match_qcov mate (s1/s2);
                    # alphabetical on exact tie.
                    winner = g1 if (s1 > s2 or (s1 == s2 and g1 < g2)) else g2
                    hits[winner] += 1
                    diagnostics['same_group_diff_gene'] += 1
                else:
                    # Genuinely different Groups = a real cross-mechanism fragment.
                    # Count ONE hit to EACH gene (both mechanisms are evidenced).
                    hits[g1] += 1
                    hits[g2] += 1
                    diagnostics['cross_group_both_counted'] += 1
            else:
                hits[g1] += 1
                hits[g2] += 1
                diagnostics['legacy_discordant_both_counted'] += 1
        else:
            for g, _ in entries:
                hits[g] += 1
            diagnostics['unexpected_multi_mate'] += 1

    return hits, diagnostics


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
            edge_aware_qcov=args.edge_aware_qcov,
            ref_lengths=ref_lengths,
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
    elif args.count_mode == "fragment":
        gene_counts, frag_diag = summarize_genes_fragment(
            per_read, passing_genes, group_aware=args.group_aware,
        )
        diag_str = ", ".join(f"{k}={v}" for k, v in sorted(frag_diag.items())) or "none"
        print(f"[INFO] Fragment resolution diagnostics: {diag_str}")
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