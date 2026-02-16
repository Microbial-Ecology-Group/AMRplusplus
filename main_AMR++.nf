nextflow.enable.dsl=2
// Example command:
// module load python/3.9.3_anaconda2021.11_mamba
// nextflow run main_AMR++.nf -profile conda --pipeline demo
// nextflow run main_AMR++.nf -profile conda --pipeline demo --kraken_db /mnt/c/Users/enriq/Dropbox/minikraken_8GB_20200312/

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
*/

// Determine pipeline type and input source
def pipelineType = "Paired-end"
def inputParam = "--reads"
def inputPath = params.reads

if (params.pipeline?.startsWith("se_") || params.pipeline == "se_AMR") {
    pipelineType = "Single-end"
    inputParam = "--reads"
    inputPath = params.reads
} else if (params.pipeline?.startsWith("merged_") || params.pipeline == "merged_AMR" || params.pipeline == "merged_AMR_wKraken") {
    pipelineType = "Merged reads"
    inputParam = "--merged_reads"
    inputPath = params.merged_reads ?: params.reads
} else if (params.pipeline?.startsWith("bam_")) {
    pipelineType = "BAM files"
    inputParam = "--bam_files"
    inputPath = params.bam_files ?: "Not specified"
}

log.info """\
 A M R + +    N F   P I P E L I N E
 ===================================
 Pipeline             : ${params.pipeline ?: 'demo'}
 Pipeline type        : ${pipelineType}
 Input (${inputParam.padRight(12)}) : ${inputPath}
 Output folder        : ${params.output}
 ===================================
"""


def helpMessage = """\
    ===============================================================================
                        AMR++ Nextflow Pipeline Help
    ===============================================================================

    USAGE:
        nextflow run main_AMR++.nf --pipeline <pipeline_name> [options]

    -------------------------------------------------------------------------------
                              MAIN PIPELINES
    -------------------------------------------------------------------------------
    
    demo                    Run a demonstration of AMR++ using test data
    standard_AMR            Standard AMR++ pipeline (QC → trim → host removal → resistome)
    fast_AMR                Fast AMR++ pipeline without host removal
    standard_AMR_wKraken    Standard pipeline with Kraken microbiome analysis

    -------------------------------------------------------------------------------
                           MERGED READ PIPELINES
    -------------------------------------------------------------------------------
    For analyzing reads that have been merged with FLASH or similar tools.

    merged_AMR              Standard pipeline for merged paired-end reads
    merged_AMR_wKraken      Merged read pipeline with Kraken analysis

    -------------------------------------------------------------------------------
                         SINGLE-END (SE) PIPELINES  
    -------------------------------------------------------------------------------
    For analyzing single-end sequencing data.

    se_AMR                  Standard AMR++ for single-end reads with Kraken
    se_AMR_wKraken          Single-end read pipeline with Kraken analysis

    -------------------------------------------------------------------------------
                         PAIRED-END SUBWORKFLOWS
    -------------------------------------------------------------------------------
    Individual analysis steps that can be run independently.

    eval_qc                 Run FastQC quality assessment
    trim_qc                 Run trimming and quality control (Trimmomatic)
    rm_host                 Remove host reads (BWA alignment to host genome)
    resistome               Perform resistome analysis (alignment + counting)
    align                   Perform alignment to MEGARes database only
    kraken                  Perform Kraken taxonomic classification
    qiime2                  Perform QIIME 2 16S rRNA analysis

    -------------------------------------------------------------------------------
                        SINGLE-END SUBWORKFLOWS
    -------------------------------------------------------------------------------
    Individual analysis steps for single-end data.

    se_eval_qc              FastQC for single-end reads
    se_trim_qc              Trimming for single-end reads
    se_rm_host              Host removal for single-end reads
    se_resistome            Resistome analysis for single-end reads
    se_kraken               Kraken analysis for single-end reads

    -------------------------------------------------------------------------------
                        MERGED READ SUBWORKFLOWS
    -------------------------------------------------------------------------------
    Individual analysis steps for merged read data.

    merge_reads             Merge paired-end reads with FLASH
    merged_rm_host          Host removal for merged reads
    merged_resistome        Resistome analysis for merged reads
    merged_kraken           Kraken analysis for merged reads

    -------------------------------------------------------------------------------
                           BAM FILE WORKFLOWS
    -------------------------------------------------------------------------------
    For analyzing pre-aligned BAM files.

    bam_resistome           Resistome analysis from BAM files
    bam_resistome_counts    Resistome counting from BAM files

    -------------------------------------------------------------------------------
                              OPTIONS
    -------------------------------------------------------------------------------

    Input/Output:
        --reads             Path to input FASTQ files (default: params.config)
        --output            Output directory (default: params.config)
        --bam_files         Path to BAM files (for bam_* pipelines)
        --merged_reads      Path to merged FASTQ files (for merged_* pipelines)

    Reference Databases:
        --host              Path to host genome for host removal
        --amr               Path to AMR database (MEGARes)
        --annotation        Path to AMR annotation file
        --kraken_db         Path to Kraken database
        --dada2_db          Path to DADA2 database (for qiime2)

    Analysis Options:
        --snp N             Skip SNP confirmation analysis (default: Y)
        --deduped Y         Include deduplicated count analysis
                            (increases runtime and storage requirements)

    -------------------------------------------------------------------------------
                              PROFILES
    -------------------------------------------------------------------------------
    Use -profile to specify your computing environment:

        local               Software already in \$PATH (default)
        local_slurm         Local installation with SLURM job control
        conda               Uses mamba to create conda environment
        conda_slurm         Conda with SLURM job control
        singularity         Uses Singularity container
        singularity_slurm   Singularity with SLURM job control
        apptainer           Uses Apptainer (.sif) container
        docker              Uses Docker container

    -------------------------------------------------------------------------------
                              EXAMPLES
    -------------------------------------------------------------------------------

    # Run demonstration
    nextflow run main_AMR++.nf -profile conda --pipeline demo

    # Standard analysis with host removal
    nextflow run main_AMR++.nf -profile conda --pipeline standard_AMR \\
        --reads 'data/*_R{1,2}.fastq.gz' --host /path/to/host.fasta

    # Fast analysis (no host removal)
    nextflow run main_AMR++.nf -profile conda --pipeline fast_AMR \\
        --reads 'data/*_R{1,2}.fastq.gz'

    # Single-end analysis
    nextflow run main_AMR++.nf -profile conda --pipeline se_AMR \\
        --reads 'data/*.fastq.gz'

    # Analyze existing BAM files
    nextflow run main_AMR++.nf -profile conda --pipeline bam_resistome \\
        --bam_files 'alignments/*.bam'

    # Run with Kraken microbiome analysis
    nextflow run main_AMR++.nf -profile conda --pipeline standard_AMR_wKraken \\
        --reads 'data/*_R{1,2}.fastq.gz' --kraken_db /path/to/kraken_db

    ===============================================================================
    For more information, visit: https://github.com/Microbial-Ecology-Group/AMRplusplus
    ===============================================================================
    """

// Helper function for consistent pipeline logging
def logPipelineStart(String pipelineName, String description) {
    log.info """
    ===============================================================================
                         AMR++ Pipeline: ${pipelineName}
    ===============================================================================
    ${description}
    -------------------------------------------------------------------------------
    """
}

Channel
    .fromFilePairs( params.reads , size: (params.reads =~ /\{/) ? 2 : 1)
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .map { id, files -> 
        def modified_baseName = id.split('\\.')[0]
        tuple(modified_baseName, files)
    }
    .set {fastq_files}

// Load main pipeline workflows
include { STANDARD_AMRplusplus } from './subworkflows/AMR++_standard.nf' 
include { FAST_AMRplusplus } from './subworkflows/AMR++_fast.nf'
include { STANDARD_AMRplusplus_wKraken } from './subworkflows/AMR++_standard_wKraken.nf'

// Load merged read workflows
include { STANDARD_merged_AMRplusplus } from './subworkflows/AMR++_merged_standard.nf'
include { STANDARD_merged_AMRplusplus_wKraken } from './subworkflows/AMR++_merged_standard_wKraken.nf'
include { FASTQ_MERGE_WF } from "$baseDir/subworkflows/fastq_merging.nf"
include { MERGED_FASTQ_RM_HOST_WF } from "$baseDir/subworkflows/fastq_host_removal.nf" 
include { MERGED_FASTQ_RESISTOME_WF } from "$baseDir/subworkflows/fastq_resistome.nf"
include { MERGED_FASTQ_KRAKEN_WF } from "$baseDir/subworkflows/fastq_microbiome.nf"

// Load SE read workflows
include { SE_AMRplusplus } from './subworkflows/AMR++_SE_standard.nf'
include { SE_AMRplusplus_wKraken } from './subworkflows/AMR++_SE_standard_wKraken.nf'
include { FASTQ_QC_SE_WF } from './subworkflows/fastq_information.nf'
include { FASTQ_TRIM_SE_WF } from './subworkflows/fastq_QC_trimming.nf'
include { FASTQ_RM_HOST_SE_WF   } from './subworkflows/fastq_host_removal.nf'
include { FASTQ_RESISTOME_SE_WF } from './subworkflows/fastq_resistome.nf'
include { FASTQ_KRAKEN_SE_WF    } from './subworkflows/fastq_microbiome.nf'


// Load subworkflows
include { FASTQ_QC_WF } from './subworkflows/fastq_information.nf'
include { FASTQ_TRIM_WF } from './subworkflows/fastq_QC_trimming.nf'
include { FASTQ_ALIGN_WF } from './subworkflows/fastq_align.nf'
include { FASTQ_RM_HOST_WF } from './subworkflows/fastq_host_removal.nf' 
include { FASTQ_RESISTOME_WF } from './subworkflows/fastq_resistome.nf'
include { FASTQ_KRAKEN_WF } from './subworkflows/fastq_microbiome.nf'
include { FASTQ_QIIME2_WF } from './subworkflows/fastq_16S_qiime2.nf'

// Load BAM subworkflows
include { BAM_RESISTOME_WF } from './subworkflows/bam_resistome.nf'
include { BAM_RESISTOME_COUNTS_WF } from './subworkflows/bam_resistome_counts.nf'



workflow {
    // =========================================================================
    //                          HELP / DEFAULT
    // =========================================================================
    if (params.pipeline == null || params.pipeline == "help") {
        println helpMessage
        FAST_AMRplusplus(fastq_files, params.amr, params.annotation)
    }

    // =========================================================================
    //                         MAIN PIPELINES
    // =========================================================================
    else if(params.pipeline == "demo") {
        logPipelineStart("Demo", 
            "Running AMR++ demonstration with test data.\n    This pipeline runs a fast analysis without host removal.")
        FAST_AMRplusplus(fastq_files, params.amr, params.annotation)
    } 
    else if(params.pipeline == "standard_AMR") {
        logPipelineStart("Standard AMR++",
            "Full pipeline: QC → Trimming → Host Removal → Resistome Analysis\n    Includes SNP confirmation for relevant genes.")
        STANDARD_AMRplusplus(fastq_files, params.host, params.amr, params.annotation)
    } 
    else if(params.pipeline == "fast_AMR") {
        logPipelineStart("Fast AMR++",
            "Streamlined pipeline without host removal step.\n    Use when host contamination is not a concern.")
        FAST_AMRplusplus(fastq_files, params.amr, params.annotation)
    } 
    else if(params.pipeline == "standard_AMR_wKraken") {
        logPipelineStart("Standard AMR++ with Kraken",
            "Full pipeline with Kraken taxonomic classification.\n    Combines resistome and microbiome analysis.")
        STANDARD_AMRplusplus_wKraken(fastq_files, params.host, params.amr, params.annotation, params.kraken_db)
    } 

    // =========================================================================
    //                    SINGLE-END PIPELINES
    // =========================================================================
    else if(params.pipeline == "se_AMR") {
        logPipelineStart("Single-End AMR++",
            "Standard AMR++ pipeline for single-end reads.\n    Includes host removal, resistome, and Kraken analysis.")
        Channel
            .fromPath(params.reads)
            .map { f -> tuple(f.name.replaceFirst(/\.f(ast)?q(\.gz)?$/, ''), f) }
            .set { read_se_ch }
        SE_AMRplusplus_wKraken( read_se_ch , params.host, params.amr, params.annotation )
    }
    else if(params.pipeline == "se_AMR_wKraken") {
        logPipelineStart("Single-End AMR++",
            "Standard AMR++ pipeline for single-end reads.\n    Includes host removal, resistome, and Kraken analysis.")
        Channel
            .fromPath(params.reads)
            .map { f -> tuple(f.name.replaceFirst(/\.f(ast)?q(\.gz)?$/, ''), f) }
            .set { read_se_ch }
        SE_AMRplusplus_wKraken( read_se_ch , params.host, params.amr, params.annotation )
    }

    // =========================================================================
    //                      MERGED READ PIPELINES
    // =========================================================================
    else if(params.pipeline == "merged_AMR") {
        logPipelineStart("Merged Reads AMR++",
            "Standard pipeline for FLASH-merged paired-end reads.\n    Handles both merged and unmerged read fractions.")
        STANDARD_merged_AMRplusplus(fastq_files, params.host, params.amr, params.annotation)
    } 
    else if(params.pipeline == "merged_AMR_wKraken") {
        logPipelineStart("Merged Reads AMR++ with Kraken",
            "Merged read pipeline with Kraken taxonomic analysis.\n    Combines resistome and microbiome analysis for merged reads.")
        STANDARD_merged_AMRplusplus_wKraken(fastq_files, params.host, params.amr, params.annotation)
    } 

    // =========================================================================
    //                     PAIRED-END SUBWORKFLOWS
    // =========================================================================
    else if(params.pipeline == "eval_qc") {
        logPipelineStart("Quality Control Evaluation",
            "Running FastQC analysis on input reads.\n    Generates quality reports for read assessment.")
        FASTQ_QC_WF( fastq_files )
    } 
    else if(params.pipeline == "trim_qc") {
        logPipelineStart("Trimming & QC",
            "Running Trimmomatic for adapter removal and quality trimming.\n    Generates trimmed reads and QC reports.")
        FASTQ_TRIM_WF( fastq_files )
    }
    else if(params.pipeline == "rm_host") {
        logPipelineStart("Host Removal",
            "Removing host-derived reads using BWA alignment.\n    Host genome: ${params.host}")
        FASTQ_RM_HOST_WF(params.host, fastq_files )
    } 
    else if(params.pipeline == "resistome") {
        logPipelineStart("Resistome Analysis",
            "Performing resistome analysis with MEGARes database.\n    Includes alignment, counting, and SNP confirmation.")
        FASTQ_RESISTOME_WF( fastq_files, params.amr, params.annotation )
    }  
    else if(params.pipeline == "align") {
        logPipelineStart("Alignment Only",
            "Aligning reads to MEGARes database with BWA.\n    Generates BAM files for downstream analysis.")
        FASTQ_ALIGN_WF( fastq_files, params.amr)
    }  
    else if(params.pipeline == "kraken") {
        logPipelineStart("Kraken Taxonomic Classification",
            "Running Kraken2 for taxonomic classification.\n    Database: ${params.kraken_db}")
        FASTQ_KRAKEN_WF( fastq_files )
    }
    else if(params.pipeline == "qiime2") {
        logPipelineStart("QIIME 2 Analysis",
            "Running QIIME 2 for 16S rRNA analysis with DADA2.\n    Database: ${params.dada2_db}")
        Channel
            .fromFilePairs( params.reads, flat: true )
            .ifEmpty { exit 1, "Read pair files could not be found: ${params.reads}" }
            .map { name, forward, reverse -> [ forward.drop(forward.findLastIndexOf{"/"})[0], forward, reverse ] }
            .map { name, forward, reverse -> [ name.toString().take(name.toString().indexOf("_")), forward, reverse ] }
            .map { name, forward, reverse -> [ name +","+ forward + ",forward\n" + name +","+ reverse +",reverse" ] }
            .flatten()
            .collectFile(name: 'manifest.txt', newLine: true, storeDir: "${params.output}/demux", seed: "sample-id,absolute-filepath,direction")
            .set { ch_manifest }
        
        FASTQ_QIIME2_WF( ch_manifest , params.dada2_db)
    }

    // =========================================================================
    //                    SINGLE-END SUBWORKFLOWS
    // =========================================================================
    else if(params.pipeline == "se_eval_qc") {
        logPipelineStart("Single-End QC Evaluation",
            "Running FastQC analysis on single-end reads.")
        Channel
            .fromPath(params.reads)
            .map { f -> tuple(f.baseName, f) }
            .set { read_se_ch }
        FASTQ_QC_SE_WF( read_se_ch )
    }
    else if(params.pipeline == "se_trim_qc") {
        logPipelineStart("Single-End Trimming",
            "Running Trimmomatic on single-end reads.")
        Channel
            .fromPath(params.reads)
            .map { f -> tuple(f.baseName, f) }
            .set { read_se_ch }
        FASTQ_TRIM_SE_WF( read_se_ch )
    }
    else if(params.pipeline == "se_rm_host") {
        logPipelineStart("Single-End Host Removal",
            "Removing host reads from single-end data.\n    Host genome: ${params.host}")
        Channel
            .fromPath(params.reads)
            .map { f -> tuple(f.baseName, f) }
            .set { read_se_ch }
        FASTQ_RM_HOST_SE_WF( params.host, read_se_ch )
    }
    else if(params.pipeline == "se_resistome") {
        logPipelineStart("Single-End Resistome",
            "Performing resistome analysis on single-end reads.")
        Channel
            .fromPath(params.reads)
            .map { f -> tuple(f.baseName, f) }
            .set { read_se_ch }
        FASTQ_RESISTOME_SE_WF( read_se_ch , params.amr, params.annotation )
    }
    else if(params.pipeline == "se_kraken") {
        logPipelineStart("Single-End Kraken",
            "Running Kraken2 on single-end reads.\n    Database: ${params.kraken_db}")
        Channel
            .fromPath(params.reads)
            .map { f -> tuple(f.baseName, f) }
            .set { read_se_ch }
        FASTQ_KRAKEN_SE_WF( read_se_ch )
    }

    // =========================================================================
    //                    MERGED READ SUBWORKFLOWS
    // =========================================================================
    else if(params.pipeline == "merge_reads") {
        logPipelineStart("Read Merging",
            "Merging paired-end reads with FLASH.\n    Generates merged and unmerged read files.")
        FASTQ_MERGE_WF( fastq_files )
    }  
    else if(params.pipeline == "merged_rm_host") {
        logPipelineStart("Merged Reads Host Removal",
            "Removing host reads from merged paired-end data.\n    Host genome: ${params.host}")
        Channel
            .fromPath( params.merged_reads, glob:true )
            .ifEmpty { error "No FASTQs match: ${params.merged_reads}" }
            .map { f ->
                def m = (f.name =~ /(.+?).(extendedFrags|notCombined)\.fastq\.gz$/)
                if( !m ) error "Unrecognised FLASH name: ${f.name}"
                def sid   = m[0][1]
                def rtype = m[0][2]
                tuple( sid, tuple(rtype, f) )
            }
            .groupTuple()
            .map { sid, list ->
                def merged_fq   = list.find { it[0] == 'extendedFrags' }?.getAt(1)
                def unmerged_fq = list.find { it[0] == 'notCombined'   }?.getAt(1)
                if( !merged_fq || !unmerged_fq )
                    error "Sample ${sid} is missing merged or unmerged FASTQ"
                tuple( sid, merged_fq, unmerged_fq )
            }
            .set { to_host_rm_ch }
        MERGED_FASTQ_RM_HOST_WF(params.host, to_host_rm_ch)
    }  
    else if(params.pipeline == "merged_resistome") {
        logPipelineStart("Merged Reads Resistome",
            "Performing resistome analysis on merged reads.\n    Analyzes both merged and unmerged fractions.")
        Channel
          .fromFilePairs( params.merged_reads, glob: true )
          .ifEmpty { error "No FASTQ files match: ${params.merged_reads}" }
          .map { sample_id, files ->
            def merged   = files.find { it.name.contains('merged')   }
            def unmerged = files.find { it.name.contains('unmerged') }
            assert merged && unmerged : "Sample $sample_id missing one of merged/unmerged"
            tuple( sample_id, merged, unmerged )
          }
          .set { to_resistome_ch }
        MERGED_FASTQ_RESISTOME_WF(to_resistome_ch, params.amr, params.annotation)
    }  
    else if(params.pipeline == "merged_kraken") {
        logPipelineStart("Merged Reads Kraken",
            "Running Kraken2 on merged reads.\n    Database: ${params.kraken_db}")
        Channel
          .fromFilePairs( params.merged_reads, glob: true )
          .ifEmpty { error "No FASTQ files match: ${params.merged_reads}" }
          .map { sample_id, files ->
            def merged   = files.find { it.name.contains('merged')   }
            def unmerged = files.find { it.name.contains('unmerged') }
            assert merged && unmerged : "Sample $sample_id missing one of merged/unmerged"
            tuple( sample_id, merged, unmerged )
          }
          .set { to_microbiome_ch }
        MERGED_FASTQ_KRAKEN_WF(to_microbiome_ch)
    }  

    // =========================================================================
    //                       BAM FILE WORKFLOWS
    // =========================================================================
    else if(params.pipeline == "bam_resistome"){
        logPipelineStart("BAM Resistome Analysis",
            "Performing resistome analysis on pre-aligned BAM files.\n    Input: ${params.bam_files}\n    Use --output to specify output directory.")
        Channel
            .fromPath(params.bam_files)
            .ifEmpty { exit 1, "BAM files could not be found: ${params.bam_files}" }
            .map { file ->
                def modified_baseName = file.baseName.split('\\.')[0]
                tuple(modified_baseName, file)
            }
            .set {bam_files_ch}
        BAM_RESISTOME_WF( bam_files_ch , params.amr, params.annotation )
    }
    else if(params.pipeline == "bam_resistome_counts"){
        logPipelineStart("BAM Resistome Counts",
            "Generating resistome counts from pre-aligned BAM files.\n    Input: ${params.bam_files}\n    Use --output to specify output directory.")
        Channel
            .fromPath(params.bam_files)
            .ifEmpty { exit 1, "BAM files could not be found: ${params.bam_files}" }
            .map { file ->
                def modified_baseName = file.baseName.split('\\.')[0]
                tuple(modified_baseName, file)
            }
            .set {bam_files_ch}
        BAM_RESISTOME_COUNTS_WF( bam_files_ch , params.amr, params.annotation )
    }

    // =========================================================================
    //                         ERROR HANDLING
    // =========================================================================
    else {
        log.error """
    ===============================================================================
                                   ERROR
    ===============================================================================
    Unknown pipeline: '${params.pipeline}'
    
    Please specify a valid pipeline using --pipeline <name>
    
    Run with --pipeline help to see all available options.
    ===============================================================================
        """
        println helpMessage
        System.exit(1)  
    }
}


workflow.onComplete {
    log.info """
    ===============================================================================
                            Pipeline Complete
    ===============================================================================
    Pipeline        : ${params.pipeline}
    Started at      : ${workflow.start}
    Finished at     : ${workflow.complete}
    Duration        : ${workflow.duration}
    Success         : ${workflow.success}
    Work directory  : ${workflow.workDir}
    Output directory: ${params.output}
    ===============================================================================
    """
}
