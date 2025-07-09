nextflow.enable.dsl=2
// Example command:
// module load python/3.9.3_anaconda2021.11_mamba
// nextflow run main_AMR++.nf -profile conda --pipeline demo
// nextflow run main_AMR++.nf -profile conda --pipeline demo --kraken_db /mnt/c/Users/enriq/Dropbox/minikraken_8GB_20200312/

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
*/

log.info """\
 A M R + +    N F   P I P E L I N E
 ===================================
 Reads being analyzed : ${params.reads}
 With the pipeline    : ${params.pipeline}
 Output folder        : ${params.output}

===================================
 Running the ${params.pipeline} pipeline
===================================
"""


def helpMessage = """\
    AMR++ Nextflow Pipeline Help
    =============================

    Available pipelines:
        - demo: Run a demonstration of AMR++
        - standard_AMR: Run the standard AMR++ pipeline
        - fast_AMR: Run the fast AMR++ pipeline without host removal.
        - standard_AMR_wKraken: Run the standard AMR++ pipeline with Kraken
    Available pipeline subworkflows:
        - eval_qc: Run FastQC analysis
        - trim_qc: Run trimming and quality control
        - rm_host: Remove host reads
        - resistome: Perform resistome analysis
        - align: Perform alignment to MEGARes database
        - kraken: Perform Kraken analysis
        - qiime2: Perform QIIME 2 analysis
        - bam_resistome: Perform resistome analysis on BAM files

    To run a specific pipeline/subworkflow, use the "--pipeline" option followed by the pipeline name:
        nextflow run main_AMR++.nf --pipeline <pipeline_name> [other_options]

    To analyze your samples or otherwise change how AMR++ runs, modify the "params.config" file 
    or add more parameters to the command line.

    Finally, consider your computing environment and modify the "-profile" option. By default,
    AMR++ assumes all software dependencies are in your \$PATH, as in the "local" profile. Here are 
    the other options:
        - local: Assumes all sofware is already in your \$PATH
        - local_slurm: Local installation and adds control over slurm job submission.
        - conda: Uses "mamba" to install a conda environment. 
        - conda_slurm: Uses "mamba" and adds control over slurm job submission.
        - singularity: Uses a "singularity" image container.
        - apptainer: Uses a "Singularity" (.sif) image container with the apptainer container system. 
        - singularity_slurm: Singularity image and adds control over slurm job submission.
        - docker: Uses a docker image container.

    Analysis for genes requiring SNP confirmation occurs by default, add `--snp N` to your command or modify the "params.txt" if you wish to skip this analysis. 
    
    To include deduplicated count analysis, add `--deduped Y` to your command. 
    Please be aware that adding deduplicated counts will significantly increase run time and temp file storage requirements.

    """

Channel
    .fromFilePairs( params.reads , size: (params.reads =~ /\{/) ? 2 : 1)
    .ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .map { id, files -> 
        def modified_baseName = id.split('\\.')[0]
        tuple(modified_baseName, files)
    }
    .set {fastq_files}

// Load null pipeline
params.pipeline = null

// Load main pipeline workflows
include { STANDARD_AMRplusplus } from './subworkflows/AMR++_standard.nf' 
include { FAST_AMRplusplus } from './subworkflows/AMR++_fast.nf'
include { STANDARD_AMRplusplus_wKraken } from './subworkflows/AMR++_standard_wKraken.nf'

// Load merged read workflows
include { STANDARD_merged_AMRplusplus } from './subworkflows/AMR++_merged_standard.nf'
include { FASTQ_MERGE_WF } from "$baseDir/subworkflows/fastq_merging.nf"
include { MERGED_FASTQ_RM_HOST_WF } from "$baseDir/subworkflows/fastq_host_removal.nf" 
include { MERGED_FASTQ_RESISTOME_WF } from "$baseDir/subworkflows/fastq_resistome.nf"

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
    if (params.pipeline == null || params.pipeline == "help") {
        println helpMessage

        //run with demo params, use params.config
        FAST_AMRplusplus(fastq_files, params.amr, params.annotation)
    }
    else if(params.pipeline == "demo") {
        log.info """\
===================================
Running a demonstration of AMR++
===================================
        """
        //run with demo params, use params.config
        FAST_AMRplusplus(fastq_files, params.amr, params.annotation)
    } 
    else if(params.pipeline == "standard_AMR") {


        STANDARD_AMRplusplus(fastq_files,params.host, params.amr, params.annotation)
    } 
    else if(params.pipeline == "fast_AMR") {

        FAST_AMRplusplus(fastq_files, params.amr, params.annotation)
    } 
    else if(params.pipeline == "standard_AMR_wKraken") {

        STANDARD_AMRplusplus_wKraken(fastq_files,params.host, params.amr, params.annotation, params.kraken_db)
    } 
    else if(params.pipeline == "eval_qc") {
        log.info """\
===================================
Running the ${params.pipeline} subworkflow
===================================
        """
        FASTQ_QC_WF( fastq_files )
    } 
    else if(params.pipeline == "trim_qc") {
        log.info """\
===================================
Running the ${params.pipeline} subworkflow
===================================
        """
        FASTQ_TRIM_WF( fastq_files )
    }
    else if(params.pipeline == "rm_host") {
        log.info """\
===================================
Running the ${params.pipeline} subworkflow
===================================
        """
        FASTQ_RM_HOST_WF(params.host, fastq_files )
    } 
    else if(params.pipeline == "resistome") {
        log.info """\
===================================
Running the ${params.pipeline} subworkflow
===================================
        """
        FASTQ_RESISTOME_WF( fastq_files, params.amr, params.annotation )
    }  
    else if(params.pipeline == "align") {
        log.info """\
===================================
Running the ${params.pipeline} subworkflow
===================================
        """
        FASTQ_ALIGN_WF( fastq_files, params.amr)
    }  
    else if(params.pipeline == "kraken") {
        log.info """\
===================================
Running the ${params.pipeline} subworkflow
===================================
        """
        FASTQ_KRAKEN_WF( fastq_files , params.kraken_db)
    }
    else if(params.pipeline == "merged_AMR") {
        STANDARD_merged_AMRplusplus(fastq_files,params.host, params.amr, params.annotation)
    } 
    else if(params.pipeline == "merge_reads") {
        FASTQ_MERGE_WF( fastq_files )
    }  
    else if(params.pipeline == "merged_rm_host") {
        Channel
            .fromPath( params.merged_reads , glob:true )
            .ifEmpty { error "No FASTQs match: ${params.merged_reads}" }
            .map { Path f ->
                // capture sample ID and read-type
                def m = (f.name =~ /(.+?)_(merged|unmerged)\.dedup\.fastq\.gz$/)
                if( !m ) error "Bad name for deduped reads: ${f.name}"
                tuple( m[0][1], tuple(m[0][2], f) )      // (sid , (type , file))
            }
            .groupTuple()                                // (sid , [ (type,file) , … ])
            .map { sid, list ->
                def merged_fq   = list.find { it[0] == 'merged'   }?.getAt(1)
                def unmerged_fq = list.find { it[0] == 'unmerged' }?.getAt(1)
                if( !merged_fq || !unmerged_fq )
                    error "Sample '${sid}' missing merged or unmerged FASTQ"
                tuple( sid, merged_fq, unmerged_fq )      // final 3-element tuple
            }
            .set { to_host_rm_ch }
        MERGED_FASTQ_RM_HOST_WF(params.host, to_host_rm_ch)
    }  
    else if(params.pipeline == "merged_resistome") {
        Channel
          .fromFilePairs( params.merged_reads, glob: true )
          .ifEmpty { error "No FASTQ files match: ${params.merged_reads}" }
          .map { sample_id, files ->
            //
            // files will be e.g.
            //   [ Path(…/S1_test_merged.dedup.fastq.gz),
            //     Path(…/S1_test_unmerged.dedup.fastq.gz) ]
            //
            def merged   = files.find { it.name.contains('merged')   }
            def unmerged = files.find { it.name.contains('unmerged') }
            assert merged && unmerged : "Sample $sample_id missing one of merged/unmerged"
            tuple( sample_id, merged, unmerged )
          }
          .set { to_resistome_ch }

        MERGED_FASTQ_RESISTOME_WF(to_resistome_ch, amr,annotation)
    
    }  


    else if(params.pipeline == "bam_resistome"){
        log.info """\
                    =======================================
                    Running resistome analysis on bam files
                    with bwa alignments to the MEGARes db.

                    Use the --bam_files argument and change
                    the --output flag to keep track of your
                    standard vs deduped data.
                    =======================================
                """
        Channel
            .fromPath(params.bam_files)
            .ifEmpty { exit 1, "bam files could not be found: ${params.bam_files}" }
            .map { file ->
            def modified_baseName = file.baseName.split('\\.')[0]
              tuple(modified_baseName, file)
                }
            .set {bam_files_ch}
        BAM_RESISTOME_WF( bam_files_ch , params.amr, params.annotation )
    }
    else if(params.pipeline == "bam_resistome_counts"){
        log.info """\
                    =======================================
                    Running resistome analysis on bam files
                    with bwa alignments to the MEGARes db.

                    Use the --bam_files argument and change
                    the --output flag to keep track of your
                    standard vs deduped data.
                    =======================================
                """
        Channel
            .fromPath(params.bam_files)
            .ifEmpty { exit 1, "bam files could not be found: ${params.bam_files}" }
            .map { file ->
            def modified_baseName = file.baseName.split('\\.')[0]
              tuple(modified_baseName, file)
                }
            .set {bam_files_ch}
        BAM_RESISTOME_COUNTS_WF( bam_files_ch , params.amr, params.annotation )
    }
    else {
            println "ERROR ################################################################"
            println "Please choose a pipeline!!!" 
            println ""
            println "To test the pipeline, use the \"demo\" pipeline or omit the pipeline flag:"
            println ""
            println "ERROR ################################################################"
            println helpMessage
            println "Exiting ..."
            System.exit(0)  
    }
}


workflow.onComplete {
    println "Pipeline completed!"
    println "Started at  $workflow.start"
    println "Finished at $workflow.complete"
    println "Time elapsed: $workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
