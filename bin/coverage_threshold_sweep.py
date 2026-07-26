#!/usr/bin/env python3
"""
coverage_threshold_sweep.py
════════════════════════════════════════════════════════════════════════════════

Single-pass-then-sweep diagnostic for choosing --min-query-coverage and
--min-gene-fraction thresholds for alignment_analyzer.py.

WHY A SEPARATE SCRIPT
──────────────────────
alignment_analyzer.py applies one fixed pair of thresholds per run. Testing
N query-coverage values x M gene-fraction values by re-running it N*M times
would re-scan a multi-gigabyte BAM N*M times. This script reads every
PRIMARY mapped alignment's (gene, mapq, query_coverage_pct, reference span,
read_length, alignment_length) ONCE into a small intermediate CSV, then
evaluates the entire threshold grid as fast DuckDB queries + in-memory
interval merges — no further BAM reads.

WHAT'S IN THIS VERSION
─────────────────────────
1. Two "retained" columns, on purpose (see below) — n_alignments_retained
   (qcov filter alone) vs n_alignments_in_passing_genes (qcov + gene_fraction
   together).
2. --exclude-snp-confirmation: drops alignments to genes flagged
   "RequiresSNPConfirmation" in MEGARes. Coverage of these genes only tells
   you the (often universally-conserved) gene is present — NOT that the
   specific resistance-conferring point mutation was observed. Applied as a
   fixed pre-filter alongside --min-mapq, before any threshold sweep.
3. Class/Mechanism/Group pass-counts added directly to the main grid output,
   so you can see how taxonomy-level detection is affected without needing
   the separate batch/taxonomy pipeline.
4. --redundancy-output: quantifies multi-mapping redundancy — groups of
   near-identical reference accessions (same Group label, similar length,
   each propped up by only a couple of reads) that are very likely the same
   underlying signal split across redundant database entries by alignment
   tie-breaking, rather than independent evidence.
5. --gene-detail-output gains taxonomy columns and per-gene read_length /
   alignment_length quantiles (p10/median/p90).
6. --length-quantiles-output: dataset-wide read_length / alignment_length
   quantiles (p5/p10/p25/median/p75/p90/p95) per query-coverage threshold —
   a quick characterization of your read-length distribution independent of
   any one gene.

READING THE OUTPUT — TWO DIFFERENT "RETAINED" COLUMNS, ON PURPOSE
────────────────────────────────────────────────────────────────────
  n_alignments_retained / pct_alignments_retained_of_baseline
      Alignments passing --min-query-coverage ALONE. Identical across every
      min_gene_fraction row for a given min_query_coverage.
  n_alignments_in_passing_genes / pct_alignments_in_passing_genes_of_baseline
      Alignments belonging to a gene that ALSO clears min_gene_fraction.
  n_genes_at_qcov_threshold / n_genes_passing_both / pct_genes_passing_of_baseline_genes
      Gene COUNTS, not alignment counts.
  n_classes_passing / n_mechanisms_passing / n_groups_passing
      Distinct taxonomy categories with >=1 gene passing both filters.

It's normal for alignment-retention to stay high while gene/category counts
drop sharply — alignment counts per gene are typically very skewed.

BREADTH APPROXIMATION
───────────────────────
Gene-length breadth is the UNION OF [reference_start, reference_end)
INTERVALS of qualifying alignments — not full CIGAR-aware position tracking
like alignment_analyzer.py's --cigar-aware-coverage. Close approximation for
typical short-read alignments; use this to narrow down thresholds, then
confirm the exact number for your final pair with --cigar-aware-coverage.

7. Fragment-aware hit counting (n_fragment_hits): a paired fragment where
   both mates agree on the same gene counts as ONE hit, not two. If only one
   mate is classified, that's still one hit. If the mates disagree, each
   gene gets its own hit (two total) — disagreement is itself informative,
   so it's not collapsed away. This applies uniformly to PE fragments and to
   the not-combined ("unmerged") subset of a merge workflow, since both
   share the same base_read_id structure; a true FLASH-merged read has no
   mate at all and trivially counts as one hit through the same code path.
   This sits alongside the existing alignment-record-based counts (read_count,
   n_alignments_retained, etc.) rather than replacing them, since those
   remain meaningful on their own — n_fragment_hits is what you want when
   comparing PE and merged workflows on equal footing, since raw alignment
   counts otherwise double-count concordant PE pairs relative to a single
   merged read covering the same physical fragment.

Usage:
    python coverage_threshold_sweep.py -i sample.bam -o sweep_results.csv \\
        --min-mapq 0 \\
        --exclude-snp-confirmation \\
        --query-coverage-sweep 0,0.5,0.6,0.7,0.8,0.9,0.95 \\
        --gene-fraction-sweep 0,0.1,0.25,0.5,0.8 \\
        --gene-detail-output gene_detail.csv \\
        --redundancy-output redundancy.csv \\
        --length-quantiles-output length_quantiles.csv
"""

import argparse
import csv
import os
import re
import statistics
import sys
import tempfile
from collections import defaultdict
from typing import Dict, List, Tuple, Optional

import duckdb
import pysam

# Same stripping convention used throughout this pipeline: a FLASH-merged
# read's name ends in '.extendedFrags...' with no mate suffix; a not-combined
# PE mate ends in '/1', '/2', '.1', or '.2'. Stripping both yields one
# base_read_id per physical fragment regardless of which case applies.
_RE_FLASH = re.compile(r'\.(extendedFrags|notCombined)[^/]*$')
_RE_PAIR  = re.compile(r'[/\.][12]$')


def get_base_read_id(read_id: str) -> str:
    base = _RE_FLASH.sub('', read_id)
    return _RE_PAIR.sub('', base)


def edge_aware_query_length(read, ref_length: int) -> int:
    """
    Coverage-anchored read length for query-coverage (qcov). Identical logic to
    alignment_analyzer.edge_aware_query_length so the two tools agree.

    Fixes two things vs raw read.query_length:
      1. Soft- vs hard-clip asymmetry (query_length includes soft clips but not
         hard clips). Full physical length is rebuilt from the CIGAR (aligned +
         both clip types, ops 4=S, 5=H) so the aligner's clip choice is irrelevant.
      2. Reference-edge overhang: clip bases that fall outside [0, ref_length)
         couldn't have aligned (the gene simply ended) and are discounted from the
         denominator, so a read running off a short gene's edge isn't penalized.

    Returns read bases that HAD A CHANCE to align; never below query_alignment_length.
    """
    if read.cigartuples is None:
        return read.query_alignment_length or 0
    aln_len = read.query_alignment_length or 0
    ops = read.cigartuples
    lead_clip = ops[0][1] if ops and ops[0][0] in (4, 5) else 0
    trail_clip = ops[-1][1] if len(ops) > 1 and ops[-1][0] in (4, 5) else 0
    total_physical = aln_len + lead_clip + trail_clip
    ref_start = read.reference_start if read.reference_start is not None else 0
    ref_end = read.reference_end if read.reference_end is not None else 0
    overhang_5 = max(0, lead_clip - ref_start)
    overhang_3 = max(0, trail_clip - (ref_length - ref_end)) if ref_length else 0
    eff = total_physical - overhang_5 - overhang_3
    return max(eff, aln_len)


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-i", "--input", required=True, nargs="+",
                    help=(
                        "Input BAM/SAM. Accepts one or more paths (e.g. "
                        "-i merged.bam unmerged.bam) to treat several files as ONE logical "
                        "sample — every alignment from every file given is combined into "
                        "the same per-gene breadth/fragment-hit computation. Use this for a "
                        "merge-workflow sample split across a merged BAM and a not-combined "
                        "(unmerged) BAM: gene_fraction needs the union of intervals from "
                        "BOTH to be correct, since breadth from two non-overlapping subsets "
                        "doesn't average — it has to be computed from the combined raw "
                        "alignments, not reconstructed from two separate breadth numbers "
                        "after the fact. Safer than pointing this at an externally-merged "
                        "Combined/ BAM, since this pipeline has previously found Combined/ "
                        "BAMs to triple-count reads depending on what went into the merge — "
                        "passing the original source files here sidesteps that risk entirely."
                    ))
    ap.add_argument("-o", "--output", required=True, help="Output sweep-grid CSV")
    ap.add_argument("--min-mapq", type=int, default=0,
                    help="Fixed MAPQ floor applied to every grid cell (default 0)")
    ap.add_argument("--edge-aware-qcov", dest="edge_aware_qcov", action="store_true", default=True,
                    help="(DEFAULT ON) Compute qcov (and match_qcov) against a coverage-anchored "
                         "read length: rebuild full length from the CIGAR so soft/hard "
                         "clips are equivalent, and discount clip bases that overhang "
                         "the reference edge (a read running off a short gene's end is "
                         "not penalized). Only affects read-axis qcov metrics; gene "
                         "breadth is unaffected. Matches alignment_analyzer.py's flag. "
                         "Use --no-edge-aware-qcov to disable.")
    ap.add_argument("--no-edge-aware-qcov", dest="edge_aware_qcov", action="store_false",
                    help="Disable edge-aware qcov (use raw read.query_length as denominator).")
    ap.add_argument("--min-aln-length", type=int, default=0,
                    help=(
                        "Fixed ABSOLUTE aligned-length floor (in bp), applied alongside "
                        "every swept --query-coverage-sweep value rather than replacing it. "
                        "Catches reads that don't have ENOUGH absolute matching sequence "
                        "even if their percentage technically passes (most relevant for "
                        "short reads). Default 0 (disabled)."
                    ))
    ap.add_argument("--max-unaligned-length", type=int, default=None,
                    help=(
                        "Fixed ABSOLUTE cap (in bp) on read_length - aln_length, applied "
                        "alongside every swept --query-coverage-sweep value. Catches the "
                        "OTHER direction: a long read (e.g. a 300bp FLASH-merged read) can "
                        "clear a 50%% query-coverage bar while still leaving 150bp "
                        "unexplained, far more absolute unaligned sequence than the same "
                        "50%% would tolerate on a 150bp read. Combined with the qcov_pct "
                        "sweep, this naturally makes longer reads satisfy a STRICTER "
                        "effective percentage without needing a different --query-coverage- "
                        "sweep value per read-length stratum: at a fixed cap, a 300bp read's "
                        "effective bar tightens automatically relative to a 150bp read's, "
                        "exactly compensating for read length rather than ignoring it. "
                        "Default None (disabled, identical to previous behavior)."
                    ))
    ap.add_argument("--exclude-snp-confirmation", action="store_true", default=False,
                    help=(
                        "Exclude alignments to genes flagged 'RequiresSNPConfirmation' "
                        "in MEGARes (e.g. RPOB, GYRA, 16S/23S rRNA mutation markers). "
                        "These are often universally-conserved housekeeping genes where "
                        "coverage alone does NOT confirm the specific resistance "
                        "mutation was observed. Applied as a fixed pre-filter, like "
                        "--min-mapq, before the threshold sweep."
                    ))
    ap.add_argument("--group-aware", "--group-aware-discordant", dest="group_aware",
                    action="store_true", default=True,
                    help=(
                        "(DEFAULT ON) Changes how a fragment whose two mates classify to "
                        "DIFFERENT gene_accessions is counted (see resolve_fragment_hits() "
                        "docstring for full reasoning). With this ON, mates that disagree on "
                        "gene_accession but share the same MEGARes Group are treated as the "
                        "SAME underlying signal split across redundant database entries (the "
                        "multi-mapping artifact tracked via --redundancy-output) and "
                        "counted as ONE hit, credited to whichever mate has the higher "
                        "match_qcov_pct. Genuine cross-Group disagreement (mates hit "
                        "truly different Groups) is instead counted as ONE hit to EACH "
                        "gene, since both mechanisms are real evidence; tracked "
                        "separately in the console diagnostics so you can see "
                        "how often each case actually occurs. Group is parsed from the "
                        "MEGARes pipe-delimited reference name. Use --no-group-aware "
                        "to restore legacy behavior (both genes get 1 hit each, 2 total). "
                        "(--group-aware-discordant is a backward-compatible alias.)"
                    ))
    ap.add_argument("--no-group-aware", "--no-group-aware-discordant", dest="group_aware",
                    action="store_false",
                    help="Disable group-aware resolution; mates that disagree on "
                         "gene_accession each get their own hit (legacy behavior).")
    ap.add_argument("--query-coverage-sweep", default="0,0.5,0.6,0.7,0.8,0.9,0.95",
                    help="Comma-separated min-query-coverage values to test (0-1). "
                         "Set to '0' to disable qcov filtering and use identity-sweep alone.")
    ap.add_argument("--gene-fraction-sweep", default="0,0.1,0.25,0.5,0.8",
                    help="Comma-separated min-gene-fraction values to test (0-1)")
    ap.add_argument("--identity-sweep", default="0",
                    help=(
                        "Comma-separated minimum percent-identity values to sweep (0-100). "
                        "Identity = (query_alignment_length - NM) / query_alignment_length * 100, "
                        "measuring how similar the ALIGNED PORTION of each read is to the "
                        "reference, independent of how much of the read aligned (qcov handles "
                        "that). Default '0' (no identity filtering, same as previous behavior). "
                        "Typical sweep: '0,80,90,95,97'. "
                        "Reads missing an NM tag pass at identity=0 but are excluded at any "
                        "higher threshold. Can be used alongside --query-coverage-sweep to "
                        "produce a 3-axis grid, or set --query-coverage-sweep 0 to sweep "
                        "identity × gene-fraction only."
                    ))
    ap.add_argument("--match-qcov-sweep", default="0",
                    help=(
                        "Comma-separated minimum 'match-only query-coverage' values to sweep "
                        "(0-100). match_qcov_pct = (query_alignment_length - NM) / "
                        "query_length * 100 — what fraction of the WHOLE read is explained "
                        "by genuine matches, with both clipped bases AND mismatches/indels "
                        "excluded from the numerator. Unlike plain --query-coverage-sweep "
                        "(which only checks clipping, blind to mismatches) this also "
                        "penalizes a fully-aligned-but-noisy read. Unlike --identity-sweep "
                        "alone (which only checks the aligned portion's quality, blind to "
                        "extent) this also penalizes a short-but-clean fragment — and "
                        "because NM >= 0, match_qcov_pct can NEVER exceed plain qcov_pct, so "
                        "a short fragment is capped at its own length-driven ceiling no "
                        "matter how clean it is; high identity cannot 'rescue' a short hit "
                        "the way it might first seem to. Default '0' (disabled, same as "
                        "previous behavior). Typical sweep for filtering on this alone: "
                        "'0,50,60,70,80,90'. Designed to be used INSTEAD of separately "
                        "sweeping --query-coverage-sweep and --identity-sweep, when you want "
                        "extent and quality combined into a single filtering decision rather "
                        "than tracked as two independent diagnostic axes — set "
                        "--query-coverage-sweep 0 --identity-sweep 0 and sweep this instead."
                    ))
    ap.add_argument("--tmp-csv", default=None,
                    help="Where to write the intermediate per-alignment CSV "
                         "(default: temp file, deleted after run)")
    ap.add_argument("--keep-tmp-csv", action="store_true", default=False,
                    help="Don't delete the intermediate per-alignment CSV after the run")
    ap.add_argument("--gene-detail-output", default=None,
                    help="Optional CSV: per-gene breadth, read count, taxonomy, and "
                         "read/alignment-length quantiles at each swept qcov threshold.")
    ap.add_argument("--redundancy-output", default=None,
                    help="Optional CSV: per-Group multi-mapping redundancy metrics at "
                         "each swept qcov threshold (see module docstring point 4).")
    ap.add_argument("--length-quantiles-output", default=None,
                    help="Optional CSV: dataset-wide read_length/alignment_length "
                         "quantiles at each swept qcov threshold.")
    return ap.parse_args()


def parse_gene_reference(ref_name: str) -> Tuple[str, str, str, str, str, str]:
    """MEGARes pipe-delimited name -> (MEG_ID, Type, Class, Mechanism, Group, SNP)."""
    parts = ref_name.split('|')
    if len(parts) >= 5:
        return (parts[0], parts[1], parts[2], parts[3], parts[4],
                parts[5] if len(parts) > 5 else '')
    return (ref_name, '', '', '', '', '')


def percentile(sorted_vals: List[float], q: float) -> float:
    """Linear-interpolation percentile (q in 0-100), robust down to n=1."""
    if not sorted_vals:
        return 0.0
    if len(sorted_vals) == 1:
        return float(sorted_vals[0])
    k = (len(sorted_vals) - 1) * (q / 100.0)
    f = int(k)
    c = min(f + 1, len(sorted_vals) - 1)
    if f == c:
        return float(sorted_vals[f])
    return sorted_vals[f] * (c - k) + sorted_vals[c] * (k - f)


def stream_primary_alignments_to_csv(
    bam_paths: List[str], out_csv: str, min_mapq: int,
    exclude_snp_confirmation: bool, edge_aware_qcov: bool = False
) -> Tuple[int, Dict[str, int], int]:
    """
    Single pass over EACH BAM in bam_paths (e.g. a sample's merged BAM and
    its not-combined/unmerged BAM together), writing every PRIMARY MAPPED
    alignment passing --min-mapq AND (if requested) not flagged
    RequiresSNPConfirmation into ONE shared per-alignment CSV. Both filters
    are fixed for the whole sweep, applied here once.

    Combining multiple files here (rather than relying on a pre-merged BAM)
    means gene_fraction/breadth downstream is computed from the TRUE union
    of intervals across every file — base_read_id collisions across files
    aren't a concern for files belonging to the same sample, since a given
    original fragment is FLASH-assigned to exactly one of merged/unmerged,
    never both, so read IDs don't collide between the two.

    Returns (total_primary_seen_before_any_filter, ref_lengths, n_excluded_by_snp).
    """
    total_seen = 0
    n_excluded_by_snp = 0
    out_dir = os.path.dirname(out_csv)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)

    ref_lengths: Dict[str, int] = {}

    with open(out_csv, "w", newline="") as fh:
        writer = csv.writer(fh)
        # 'blocks' replaces the old ref_start/ref_end pair. get_blocks()
        # returns the CIGAR-aware list of gapless aligned segments, correctly
        # splitting around deletions (D/N) — a read with CIGAR 40M10D50M
        # produces "start:50;60:end" rather than one span "start:end" that
        # silently counts the 10bp deletion gap as covered. This matches
        # ResistomeAnalyzer's per-base CIGAR walk for deletions, and fixes
        # the one blind spot the old interval-span approach had.
        writer.writerow(["gene_accession", "mapq", "qcov_pct",
                         "read_length", "aln_length", "base_read_id",
                         "pct_identity", "match_qcov_pct", "matched_length", "blocks"])

        for bam_path in bam_paths:
            with pysam.AlignmentFile(bam_path, "rb", check_sq=False) as bam:
                ref_lengths.update(dict(zip(bam.references, bam.lengths)))

                for read in bam.fetch(until_eof=True):
                    if read.flag & 0x900:
                        continue
                    if read.is_unmapped:
                        continue
                    total_seen += 1

                    if read.mapping_quality < min_mapq:
                        continue

                    gene = bam.get_reference_name(read.reference_id)

                    if exclude_snp_confirmation:
                        snp_field = gene.split('|')[-1] if '|' in gene else ''
                        if snp_field == 'RequiresSNPConfirmation':
                            n_excluded_by_snp += 1
                            continue

                    qlen = read.query_length or 0
                    aln_len = read.query_alignment_length or 0
                    if edge_aware_qcov:
                        _ref_len = ref_lengths.get(gene, 0)
                        qlen = edge_aware_query_length(read, _ref_len)
                    qcov_pct = round(aln_len / qlen * 100.0, 2) if qlen else 0.0
                    base_id = get_base_read_id(read.query_name)

                    # NM tag (edit distance) = mismatches + inserted bases +
                    # deleted reference bases. Absent in some poorly-formed SAMs;
                    # we write NULL rather than silently assigning 0% or 100%.
                    # pct_identity = (aln_len - NM) / aln_len * 100
                    #   -> 100% = every aligned base matches perfectly
                    #   -> lower values = sum of all edit events as a fraction
                    #      of the aligned (non-clipped) query length.
                    # Soft clips are correctly excluded from both aln_len and NM
                    # so clipping alone doesn't penalize identity — only real
                    # mismatches/indels inside the aligned region do.
                    try:
                        nm = read.get_tag('NM')
                        pct_identity = round((aln_len - nm) / aln_len * 100.0, 2) if aln_len else 0.0
                        # match_qcov_pct = what fraction of the WHOLE read is
                        # explained by genuine matches (clipped bases AND
                        # mismatches/indels both excluded from the numerator).
                        # Mathematically: match_qcov_pct = qcov_pct * pct_identity / 100,
                        # and since NM >= 0, match_qcov_pct can never exceed
                        # qcov_pct — a short fragment is capped at its own
                        # qcov_pct ceiling regardless of how clean it is, so
                        # this metric can't be "rescued" by high identity the
                        # way a naive intuition might worry about. It directly
                        # answers "how much of this read is confidently,
                        # correctly explained by this reference," combining
                        # extent and quality into one number.
                        match_qcov_pct = round((aln_len - nm) / qlen * 100.0, 2) if qlen else 0.0
                        # matched_length = absolute bp of TRUE matches, i.e.
                        # aln_length with mismatches/indels (NM) subtracted
                        # out. The exact integer underlying both pct_identity
                        # (matched_length / aln_length) and match_qcov_pct
                        # (matched_length / read_length) — stored directly so
                        # downstream length distributions don't have to
                        # re-derive it from already-rounded percentages.
                        matched_length = aln_len - nm
                    except KeyError:
                        pct_identity = None
                        match_qcov_pct = None
                        matched_length = None

                    raw_blocks = read.get_blocks()
                    blocks_str = (
                        ";".join(f"{s}:{e}" for s, e in raw_blocks)
                        if raw_blocks
                        else f"{read.reference_start}:{read.reference_end}"
                    )

                    writer.writerow([gene, read.mapping_quality, qcov_pct,
                                     qlen, aln_len, base_id,
                                     pct_identity, match_qcov_pct, matched_length, blocks_str])

    return total_seen, ref_lengths, n_excluded_by_snp


def resolve_fragment_hits(
    rows: List[Tuple[str, str, float]],
    taxonomy_cache: Optional[Dict[str, Tuple[str, str, str, str, str, str]]] = None,
    group_aware: bool = False,
) -> Tuple[Dict[str, int], Dict[str, int]]:
    """
    rows: list of (gene_accession, base_read_id, match_qcov_pct) for
    alignments passing the current threshold combination. match_qcov_pct is
    only used as a tie-breaker when group_aware=True; pass 0.0 if unused.

    DEFAULT (group_aware=False — exact previous behavior, unchanged):
      both mates agree on gene_accession -> 1 hit
      only one mate present              -> 1 hit
      mates disagree on gene_accession   -> 1 hit to EACH gene (2 total)

    group_aware=True (opt-in): before treating a gene-level disagreement as
    real discordance, check whether the two accessions share the same
    MEGARes Group. If they do, this is almost certainly the multi-mapping/
    redundant-accession artifact this pipeline has tracked all along via
    the redundancy diagnostics (BWA's tie-breaking among near-identical
    database entries) — both mates agree on what's biologically there, they
    just landed on different specific accessions. Counted as ONE hit,
    credited to whichever mate has the higher match_qcov_pct (the
    established combined extent+quality confidence metric).

    Genuine cross-Group disagreement (the two accessions belong to
    DIFFERENT Groups entirely — e.g. one mate hits a beta-lactamase, the
    other an aminoglycoside-modifying enzyme) is a real, more concerning
    signal: a chimeric fragment, PCR artifact, contamination, or a bad
    FLASH merge. Also resolved to ONE hit (the higher-match_qcov_pct mate),
    rather than doubling the count for what's likely noise either way —
    but tracked separately in diagnostics so you can see how often this
    actually happens, since it's a meaningfully different situation from
    same-Group disagreement.

    A true single/merged read has no mate sharing its base_read_id, so it
    falls into the "only one present" case trivially — same code path, no
    special-casing needed for merged vs PE.

    Returns (hits, diagnostics):
      hits        — gene_accession -> fragment-hit count
      diagnostics — counts of each resolution case actually encountered:
                    'same_gene', 'one_mate_only',
                    'same_group_diff_gene' (group_aware only),
                    'cross_group_both_counted' (group_aware only),
                    'legacy_discordant_both_counted' (group_aware=False only)
    """
    frag_groups: Dict[str, List[Tuple[str, float]]] = defaultdict(list)
    for gene, base_id, mq in rows:
        frag_groups[base_id].append((gene, mq))

    hits: Dict[str, int] = defaultdict(int)
    diagnostics: Dict[str, int] = defaultdict(int)

    for base_id, entries in frag_groups.items():
        if len(entries) == 1:
            hits[entries[0][0]] += 1
            diagnostics['one_mate_only'] += 1
        elif len(entries) == 2:
            (g1, mq1), (g2, mq2) = entries
            if g1 == g2:
                hits[g1] += 1
                diagnostics['same_gene'] += 1
            elif group_aware:
                group1 = (taxonomy_cache or {}).get(g1, (g1, '', '', '', '', ''))[4]
                group2 = (taxonomy_cache or {}).get(g2, (g2, '', '', '', '', ''))[4]
                if group1 and group1 == group2:
                    # Same Group, different accession = redundant-DB-entry artifact.
                    # Collapse to ONE hit for the higher-match_qcov mate (mq1/mq2).
                    winner = g1 if mq1 >= mq2 else g2
                    hits[winner] += 1
                    diagnostics['same_group_diff_gene'] += 1
                else:
                    # Genuinely different Groups = a real cross-mechanism fragment.
                    # Count ONE hit to EACH gene (both mechanisms evidenced).
                    hits[g1] += 1
                    hits[g2] += 1
                    diagnostics['cross_group_both_counted'] += 1
            else:
                hits[g1] += 1          # mates disagree -> one hit each (legacy)
                hits[g2] += 1
                diagnostics['legacy_discordant_both_counted'] += 1
        else:
            # Shouldn't happen — each base_read_id should have at most one
            # primary alignment per mate (R1, R2, or single). If it does,
            # fall back to counting each entry separately rather than guessing.
            for g, _ in entries:
                hits[g] += 1
            diagnostics['unexpected_3plus'] += 1
    return hits, diagnostics


def merged_interval_length(intervals: List[Tuple[int, int]]) -> int:
    if not intervals:
        return 0
    intervals = sorted(intervals)
    total = 0
    cur_start, cur_end = intervals[0]
    for s, e in intervals[1:]:
        if s <= cur_end:
            cur_end = max(cur_end, e)
        else:
            total += cur_end - cur_start
            cur_start, cur_end = s, e
    total += cur_end - cur_start
    return total


def compute_redundancy_table(
    threshold_labels: Dict[str, float],
    gene_data: Dict[str, dict],
    gene_fragment_hits: Dict[str, int],
    ref_lengths: Dict[str, int],
    taxonomy_cache: Dict[str, Tuple[str, str, str, str, str, str]],
) -> List[dict]:
    """
    Groups gene accessions by their taxonomy Group field and quantifies
    multi-mapping redundancy: many accessions, each backed by only a few
    reads, all roughly the same length, with no single accession dominating
    the read count -> strong sign these are the SAME underlying signal split
    across near-duplicate reference entries by alignment tie-breaking,
    rather than independent distinct hits.

    threshold_labels: dict of column_name -> value, prepended to every output
    row (e.g. {"min_query_coverage": 0.0, "min_identity": 90.0,
    "min_match_qcov": 0.0}) — generalized from a single qcov argument so the
    redundancy table stays labeled correctly across however many threshold
    axes are actually being swept.

    Reports both the raw alignment-record count (total_reads, as before) and
    the fragment-hit count (total_fragment_hits). A concordant PE pair
    inflates total_reads to 2 for what is really one piece of evidence —
    total_fragment_hits is the fairer number when comparing a PE workflow's
    redundancy against a merged workflow's.
    """
    by_group: Dict[str, List[Tuple[str, int, int, int]]] = defaultdict(list)
    for gene, d in gene_data.items():
        n = len(d['read_lengths'])
        if n == 0:
            continue
        _, _, _, _, group, _ = taxonomy_cache.get(gene, (gene, '', '', '', '', ''))
        if not group:
            continue
        by_group[group].append((gene, n, ref_lengths.get(gene, 0), gene_fragment_hits.get(gene, 0)))

    rows = []
    for group, entries in by_group.items():
        n_accessions = len(entries)
        total_reads = sum(e[1] for e in entries)
        max_reads = max(e[1] for e in entries)
        total_hits = sum(e[3] for e in entries)
        max_hits = max(e[3] for e in entries) if entries else 0
        lengths = [e[2] for e in entries if e[2] > 0]
        mean_length = statistics.mean(lengths) if lengths else 0.0
        length_cv = (
            statistics.pstdev(lengths) / mean_length
            if (len(lengths) > 1 and mean_length > 0) else 0.0
        )
        pct_top = round(max_reads / total_reads * 100, 2) if total_reads else 0.0
        pct_top_hits = round(max_hits / total_hits * 100, 2) if total_hits else 0.0

        # Heuristic flag — adjust thresholds to taste; raw columns are kept
        # so you can re-sort/re-filter with your own criteria instead.
        likely_redundant = (n_accessions >= 5 and pct_top < 50 and length_cv < 0.10)

        rows.append({
            **threshold_labels,
            "group": group,
            "n_accessions": n_accessions,
            "total_reads": total_reads,
            "max_accession_reads": max_reads,
            "pct_reads_in_top_accession": pct_top,
            "total_fragment_hits": total_hits,
            "max_accession_fragment_hits": max_hits,
            "pct_fragment_hits_in_top_accession": pct_top_hits,
            "mean_gene_length": round(mean_length, 1),
            "gene_length_cv": round(length_cv, 4),
            "likely_redundant": likely_redundant,
        })
    return rows


def main():
    args = parse_args()

    qcov_values       = sorted(float(x) for x in args.query_coverage_sweep.split(","))
    gf_values         = sorted(float(x) for x in args.gene_fraction_sweep.split(","))
    identity_values   = sorted(float(x) for x in args.identity_sweep.split(","))
    match_qcov_values = sorted(float(x) for x in args.match_qcov_sweep.split(","))

    for values, name, lo, hi in [
        (qcov_values,       "query-coverage", 0.0, 1.0),
        (gf_values,         "gene-fraction",  0.0, 1.0),
        (identity_values,   "identity",       0.0, 100.0),
        (match_qcov_values, "match-qcov",     0.0, 100.0),
    ]:
        for x in values:
            if not lo <= x <= hi:
                sys.exit(f"[ERROR] {name} sweep values must be {lo}-{hi}, got {x}")

    tmp_csv = args.tmp_csv or tempfile.mktemp(suffix=".csv")
    if len(args.input) > 1:
        print(f"[INFO] Scanning {len(args.input)} input files as ONE sample -> {tmp_csv}", flush=True)
        for p in args.input:
            print(f"         {p}", flush=True)
    else:
        print(f"[INFO] Scanning {args.input[0]} once -> {tmp_csv}", flush=True)
    if args.exclude_snp_confirmation:
        print("[INFO] --exclude-snp-confirmation is ON: genes flagged "
              "RequiresSNPConfirmation will be dropped entirely from this analysis.",
              flush=True)
    if args.min_aln_length > 0:
        print(f"[INFO] --min-aln-length {args.min_aln_length}bp is ON: applied alongside "
              f"every swept threshold combination.", flush=True)
    if args.max_unaligned_length is not None:
        print(f"[INFO] --max-unaligned-length {args.max_unaligned_length}bp is ON.", flush=True)
    if identity_values != [0.0]:
        print(f"[INFO] Identity sweep: {identity_values}", flush=True)
    else:
        print("[INFO] Identity sweep: off (all reads pass; set --identity-sweep to enable)", flush=True)
    if match_qcov_values != [0.0]:
        print(f"[INFO] Match-qcov sweep: {match_qcov_values}", flush=True)
    else:
        print("[INFO] Match-qcov sweep: off (all reads pass; set --match-qcov-sweep to enable)", flush=True)

    total_seen, ref_lengths, n_excluded_by_snp = stream_primary_alignments_to_csv(
        args.input, tmp_csv, args.min_mapq, args.exclude_snp_confirmation,
        edge_aware_qcov=args.edge_aware_qcov
    )
    print(f"[INFO] {total_seen:,} primary mapped alignments seen "
          f"(before any filtering)", flush=True)
    if args.exclude_snp_confirmation:
        print(f"[INFO] {n_excluded_by_snp:,} alignments excluded as "
              f"RequiresSNPConfirmation (after --min-mapq, before query-coverage sweep)",
              flush=True)

    con = duckdb.connect(":memory:")
    con.execute(f"CREATE TABLE aln AS SELECT * FROM read_csv('{tmp_csv}', header=true)")

    baseline_n = con.execute("SELECT COUNT(*) FROM aln").fetchone()[0]
    baseline_genes = con.execute("SELECT COUNT(DISTINCT gene_accession) FROM aln").fetchone()[0]

    distinct_genes = [r[0] for r in con.execute("SELECT DISTINCT gene_accession FROM aln").fetchall()]
    taxonomy_cache = {g: parse_gene_reference(g) for g in distinct_genes}

    def distinct_taxonomy_count(idx: int) -> int:
        return len({taxonomy_cache[g][idx] for g in distinct_genes if taxonomy_cache[g][idx]})

    baseline_classes = distinct_taxonomy_count(2)
    baseline_mechanisms = distinct_taxonomy_count(3)
    baseline_groups = distinct_taxonomy_count(4)

    print(f"[INFO] Baseline (--min-mapq {args.min_mapq}"
          f"{', SNP-confirmation excluded' if args.exclude_snp_confirmation else ''}): "
          f"{baseline_n:,} alignments, {baseline_genes:,} genes, "
          f"{baseline_classes:,} classes, {baseline_mechanisms:,} mechanisms, "
          f"{baseline_groups:,} groups\n", flush=True)

    rows_out = []
    gene_detail_rows = []
    redundancy_rows = []
    length_quantile_rows = []
    baseline_fragment_hits: Optional[int] = None   # captured on the first (smallest) qcov below

    QUANTILES_OVERALL = [5, 10, 25, 50, 75, 90, 95]

    for qcov in qcov_values:
      for identity in identity_values:
       for match_qcov in match_qcov_values:
        # Build the SQL WHERE clause. NULL handling for identity AND
        # match_qcov: at threshold=0, reads missing the underlying NM tag
        # pass (treated as unfiltered); at any positive threshold, NULL
        # values are excluded — we can't verify them, so they shouldn't
        # silently pass.
        base_where = (
            "qcov_pct >= ? AND aln_length >= ? "
            "AND (? = 0 OR (pct_identity IS NOT NULL AND pct_identity >= ?)) "
            "AND (? = 0 OR (match_qcov_pct IS NOT NULL AND match_qcov_pct >= ?))"
        )
        params_base = [qcov * 100.0, args.min_aln_length, identity, identity,
                       match_qcov, match_qcov]

        if args.max_unaligned_length is not None:
            sub = con.execute(
                "SELECT gene_accession, read_length, aln_length, base_read_id, "
                f"COALESCE(pct_identity, 0.0) AS pct_identity, "
                f"COALESCE(matched_length, 0) AS matched_length, blocks FROM aln "
                f"WHERE {base_where} AND (read_length - aln_length) <= ?",
                params_base + [args.max_unaligned_length],
            ).fetchall()
        else:
            sub = con.execute(
                "SELECT gene_accession, read_length, aln_length, base_read_id, "
                f"COALESCE(pct_identity, 0.0) AS pct_identity, "
                f"COALESCE(matched_length, 0) AS matched_length, blocks FROM aln WHERE {base_where}",
                params_base,
            ).fetchall()

        n_retained = len(sub)
        pct_retained = round(n_retained / baseline_n * 100, 2) if baseline_n else 0.0

        gene_data: Dict[str, dict] = defaultdict(lambda: {
            'intervals': [], 'read_lengths': [], 'aln_lengths': [],
            'identities': [], 'matched_lengths': []
        })
        for gene, rl, al, base_id, pid, mlen, blocks_str in sub:
            d = gene_data[gene]
            for block in blocks_str.split(';'):
                s_str, e_str = block.split(':')
                d['intervals'].append((int(s_str), int(e_str)))
            d['read_lengths'].append(rl)
            d['aln_lengths'].append(al)
            d['identities'].append(pid)
            d['matched_lengths'].append(mlen)

        gene_fragment_hits, discordant_diagnostics = resolve_fragment_hits(
            [(gene, base_id, (mlen / rl * 100.0 if rl else 0.0))
             for gene, rl, _, base_id, _, mlen, _ in sub],
            taxonomy_cache=taxonomy_cache,
            group_aware=args.group_aware,
        )
        n_fragment_hits_retained = sum(gene_fragment_hits.values())
        if baseline_fragment_hits is None:
            baseline_fragment_hits = n_fragment_hits_retained
        pct_fragment_hits_retained = (
            round(n_fragment_hits_retained / baseline_fragment_hits * 100, 2)
            if baseline_fragment_hits else 0.0
        )

        gene_fraction_at_filter: Dict[str, float] = {}
        gene_aln_count: Dict[str, int] = {}
        for gene, d in gene_data.items():
            length = ref_lengths.get(gene, 0)
            covered = merged_interval_length(d['intervals'])
            gene_fraction_at_filter[gene] = (covered / length) if length else 0.0
            gene_aln_count[gene] = len(d['read_lengths'])

        n_genes_at_filter = len(gene_fraction_at_filter)

        threshold_labels = {
            "min_query_coverage": qcov,
            "min_identity": identity,
            "min_match_qcov": match_qcov,
        }

        # ── redundancy and length-quantile outputs — written once per
        # (qcov, identity, match_qcov), not once per gf (a downstream filter) ──
        if args.redundancy_output:
            redundancy_rows.extend(
                compute_redundancy_table(threshold_labels, gene_data, gene_fragment_hits,
                                         ref_lengths, taxonomy_cache)
            )

        if args.length_quantiles_output:
            all_rl = sorted(rl for d in gene_data.values() for rl in d['read_lengths'])
            all_al = sorted(al for d in gene_data.values() for al in d['aln_lengths'])
            all_ml = sorted(ml for d in gene_data.values() for ml in d['matched_lengths'])
            all_id = sorted(pid for d in gene_data.values() for pid in d['identities'])
            # TRUE per-read match_qcov_pct quantiles: matched_lengths and
            # read_lengths are parallel lists (same index = same read) within
            # each gene's accumulator, so pairing them with zip() and dividing
            # PER READ before taking quantiles is mathematically correct —
            # unlike dividing matched_length_p{q} by read_length_p{q} (the
            # quantile of a ratio is NOT the ratio of quantiles in general).
            all_mq = sorted(
                (ml / rl * 100.0 if rl else 0.0)
                for d in gene_data.values()
                for rl, ml in zip(d['read_lengths'], d['matched_lengths'])
            )
            qrow = {**threshold_labels, "n_alignments": len(all_rl)}
            for q in QUANTILES_OVERALL:
                qrow[f"read_length_p{q}"] = round(percentile(all_rl, q), 1)
            for q in QUANTILES_OVERALL:
                qrow[f"aln_length_p{q}"] = round(percentile(all_al, q), 1)
            # matched_length: same as aln_length above but with mismatches/
            # indels (NM) subtracted out — the "with vs without accounting
            # for matches" comparison. Plot aln_length_p* and matched_length_p*
            # side by side to see how much of the reported aligned length is
            # actually composed of true matches vs noise.
            for q in QUANTILES_OVERALL:
                qrow[f"matched_length_p{q}"] = round(percentile(all_ml, q), 1)
            for q in QUANTILES_OVERALL:
                qrow[f"pct_identity_p{q}"] = round(percentile(all_id, q), 2)
            for q in QUANTILES_OVERALL:
                qrow[f"match_qcov_pct_p{q}"] = round(percentile(all_mq, q), 2)
            length_quantile_rows.append(qrow)

        if args.gene_detail_output:
            for gene, frac in gene_fraction_at_filter.items():
                d = gene_data[gene]
                rl_sorted = sorted(d['read_lengths'])
                al_sorted = sorted(d['aln_lengths'])
                ml_sorted = sorted(d['matched_lengths'])
                meg_id, gtype, gclass, mech, group, snp = taxonomy_cache.get(
                    gene, (gene, '', '', '', '', '')
                )
                gene_detail_rows.append({
                    **threshold_labels,
                    "gene_accession": gene,
                    "meg_id": meg_id, "type": gtype, "class": gclass,
                    "mechanism": mech, "group": group, "snp": snp,
                    "gene_length": ref_lengths.get(gene, 0),
                    "covered_bases": merged_interval_length(d['intervals']),
                    "gene_fraction": round(frac, 4),
                    "read_count": gene_aln_count[gene],
                    "n_fragment_hits": gene_fragment_hits.get(gene, 0),
                    "pct_of_baseline_alignments": round(
                        gene_aln_count[gene] / baseline_n * 100, 4) if baseline_n else 0.0,
                    "identity_p10": round(percentile(sorted(d['identities']), 10), 2),
                    "identity_median": round(percentile(sorted(d['identities']), 50), 2),
                    "identity_p90": round(percentile(sorted(d['identities']), 90), 2),
                    "read_length_p10": round(percentile(rl_sorted, 10), 1),
                    "read_length_median": round(percentile(rl_sorted, 50), 1),
                    "read_length_p90": round(percentile(rl_sorted, 90), 1),
                    "aln_length_p10": round(percentile(al_sorted, 10), 1),
                    "aln_length_median": round(percentile(al_sorted, 50), 1),
                    "aln_length_p90": round(percentile(al_sorted, 90), 1),
                    # matched_length: aln_length with mismatches/indels (NM)
                    # subtracted — "alignment length accounting for only
                    # matches." Compare against aln_length_p* above to see
                    # how much of the reported alignment is genuine match
                    # vs noise, per gene.
                    "matched_length_p10": round(percentile(ml_sorted, 10), 1),
                    "matched_length_median": round(percentile(ml_sorted, 50), 1),
                    "matched_length_p90": round(percentile(ml_sorted, 90), 1),
                })

        # ── main grid: (qcov × identity × match_qcov × gf) ───────────────────
        gf_group_summary = []  # (gf, n_groups_passing) — surfaced in the console line below,
                                # since the per-threshold print previously only showed counts
                                # from BEFORE the gf filter (n_genes_at_filter), making it look
                                # like gene-fraction wasn't being swept at all even though it was.
        for gf in gf_values:
            passing_genes = [g for g, f in gene_fraction_at_filter.items() if f >= gf]
            n_genes_passing = len(passing_genes)
            pct_genes_passing = round(n_genes_passing / baseline_genes * 100, 2) if baseline_genes else 0.0

            n_aln_in_passing_genes = sum(gene_aln_count[g] for g in passing_genes)
            pct_aln_in_passing_genes = round(n_aln_in_passing_genes / baseline_n * 100, 2) if baseline_n else 0.0

            n_hits_in_passing_genes = sum(gene_fragment_hits.get(g, 0) for g in passing_genes)
            pct_hits_in_passing_genes = (
                round(n_hits_in_passing_genes / baseline_fragment_hits * 100, 2)
                if baseline_fragment_hits else 0.0
            )

            passing_classes    = {taxonomy_cache[g][2] for g in passing_genes if taxonomy_cache[g][2]}
            passing_mechanisms = {taxonomy_cache[g][3] for g in passing_genes if taxonomy_cache[g][3]}
            passing_groups     = {taxonomy_cache[g][4] for g in passing_genes if taxonomy_cache[g][4]}
            gf_group_summary.append((gf, len(passing_groups)))

            rows_out.append({
                **threshold_labels,
                "min_gene_fraction":                          gf,
                "n_alignments_retained":                      n_retained,
                "pct_alignments_retained_of_baseline":        pct_retained,
                "n_alignments_in_passing_genes":              n_aln_in_passing_genes,
                "pct_alignments_in_passing_genes_of_baseline": pct_aln_in_passing_genes,
                "n_fragment_hits_retained":                   n_fragment_hits_retained,
                "pct_fragment_hits_retained_of_baseline":     pct_fragment_hits_retained,
                "n_fragment_hits_in_passing_genes":           n_hits_in_passing_genes,
                "pct_fragment_hits_in_passing_genes_of_baseline": pct_hits_in_passing_genes,
                "n_genes_at_filter":                          n_genes_at_filter,
                "n_genes_passing_both":                       n_genes_passing,
                "pct_genes_passing_of_baseline_genes":        pct_genes_passing,
                "n_classes_passing":                          len(passing_classes),
                "pct_classes_passing":  round(len(passing_classes)    / baseline_classes    * 100, 2) if baseline_classes    else 0.0,
                "n_mechanisms_passing":                       len(passing_mechanisms),
                "pct_mechanisms_passing": round(len(passing_mechanisms) / baseline_mechanisms * 100, 2) if baseline_mechanisms else 0.0,
                "n_groups_passing":                           len(passing_groups),
                "pct_groups_passing":   round(len(passing_groups)     / baseline_groups     * 100, 2) if baseline_groups     else 0.0,
            })

        gf_console_summary = " ".join(f"gf>={gf:.2g}:{n:,}grp" for gf, n in gf_group_summary)

        discordant_suffix = ""
        if args.group_aware:
            n_same_group = discordant_diagnostics.get('same_group_diff_gene', 0)
            n_cross_group = discordant_diagnostics.get('cross_group_both_counted', 0)
            n_disagreeing = n_same_group + n_cross_group
            discordant_suffix = (
                f" | {n_disagreeing:,} disagreeing pairs "
                f"({n_same_group:,} same-Group collapsed, {n_cross_group:,} cross-Group counted-both)"
            )

        print(f"  qcov>={qcov:.0%} identity>={identity:.0f}% match_qcov>={match_qcov:.0f}%: "
              f"{n_retained:,} alignments ({pct_retained:.1f}% of baseline), "
              f"{n_fragment_hits_retained:,} fragment hits ({pct_fragment_hits_retained:.1f}%), "
              f"{n_genes_at_filter:,} genes detectable (pre-gf) | {gf_console_summary}{discordant_suffix}",
              flush=True)

    # ── write outputs ────────────────────────────────────────────────────────────
    out_dir = os.path.dirname(args.output)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows_out[0].keys()))
        writer.writeheader()
        writer.writerows(rows_out)
    print(f"\n[INFO] Wrote {len(rows_out)} grid rows to {args.output}", flush=True)
    print("[INFO] Compare n_alignments_retained (qcov filter only) against "
          "n_alignments_in_passing_genes (qcov + gene_fraction together) to see "
          "how much alignment volume sits in genes with poor breadth.", flush=True)

    if args.gene_detail_output and gene_detail_rows:
        detail_dir = os.path.dirname(args.gene_detail_output)
        if detail_dir:
            os.makedirs(detail_dir, exist_ok=True)
        with open(args.gene_detail_output, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(gene_detail_rows[0].keys()))
            writer.writeheader()
            writer.writerows(gene_detail_rows)
        print(f"[INFO] Wrote {len(gene_detail_rows):,} per-gene rows to {args.gene_detail_output}",
              flush=True)
        print("[INFO] Long tail: sort by read_count asc + gene_fraction asc. "
              "Suspicious conserved-motif genes: sort by read_count DESC + gene_fraction asc.",
              flush=True)

    if args.redundancy_output and redundancy_rows:
        red_dir = os.path.dirname(args.redundancy_output)
        if red_dir:
            os.makedirs(red_dir, exist_ok=True)
        with open(args.redundancy_output, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(redundancy_rows[0].keys()))
            writer.writeheader()
            writer.writerows(redundancy_rows)
        print(f"[INFO] Wrote {len(redundancy_rows):,} group-redundancy rows to {args.redundancy_output}",
              flush=True)

        flagged = [r for r in redundancy_rows if r["min_query_coverage"] == qcov_values[0] and r["likely_redundant"]]
        flagged.sort(key=lambda r: -r["n_accessions"])
        if flagged:
            print(f"[INFO] Top likely-redundant groups at qcov={qcov_values[0]:.0%} "
                  f"(many near-identical accessions, no dominant single hit):", flush=True)
            for r in flagged[:10]:
                print(f"    {r['group']:<20} n_accessions={r['n_accessions']:<4} "
                      f"total_reads={r['total_reads']:<6} "
                      f"top_accession_share={r['pct_reads_in_top_accession']:.1f}% "
                      f"length_cv={r['gene_length_cv']:.3f}", flush=True)

    if args.length_quantiles_output and length_quantile_rows:
        lq_dir = os.path.dirname(args.length_quantiles_output)
        if lq_dir:
            os.makedirs(lq_dir, exist_ok=True)
        with open(args.length_quantiles_output, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(length_quantile_rows[0].keys()))
            writer.writeheader()
            writer.writerows(length_quantile_rows)
        print(f"[INFO] Wrote {len(length_quantile_rows)} length-quantile rows to "
              f"{args.length_quantiles_output}", flush=True)

    if not args.tmp_csv and not args.keep_tmp_csv:
        os.remove(tmp_csv)
    elif args.keep_tmp_csv:
        print(f"[INFO] Intermediate per-alignment CSV kept at: {tmp_csv}", flush=True)


if __name__ == "__main__":
    main()