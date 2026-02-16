#!/bin/bash
#SBATCH -J resistome_counts_80_alignment -o log_resistome_analyzer_80_alignment.out -t 12:00:00 --mem=30G --nodes=1 --ntasks=1 --cpus-per-task=1

# Input directory containing BAM files
input_dir="UM_uncat_AMR++_results/Alignment/BAM_files/Deduped"

# Output folder name
out_dir="deduped_resistome_counts_80gf_alignment"

# Create output directory if it doesn't exist
mkdir -p "$out_dir"

# Loop through all BAM files matching the pattern
for bam in "${input_dir}"/*_alignment_dedup.bam; do
    # Extract filename without path
    filename=$(basename "$bam")
    
    # Extract sample ID (everything before _alignment_dedup.bam)
    sample_id="${filename%_alignment_dedup.bam}"
    
    echo "Processing: $bam -> Sample ID: $sample_id"
    
    python3 alignment_analyzer.py \
        -i "$bam" \
        -r "${out_dir}/${sample_id}_per_read.tsv" \
        -g "${out_dir}/${sample_id}_gene_summary.tsv" \
        --sample-id "$sample_id" \
        --min-mapq 0 \
        --count-mode alignment \
        --min-gene-fraction 0.8 \
        --include-supplementary \
        --cigar-aware-coverage \
        --coverage-output "${out_dir}/${sample_id}_coverage_stats.tsv"
done

echo "Done! Results in ${out_dir}"