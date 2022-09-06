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
 reads        : ${params.reads}
 output       : ${params.output}
 """

Channel
.fromFilePairs( params.reads , size: ("${params.reads}" =~ /\{/) ? 2 : 1)
.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
.set { fastq_files }

// Default is pipeline is null to warn users below
params.pipeline = null

// Load main pipeline workflows

include { STANDARD_AMRplusplus } from './subworkflows/AMR++_standard.nf' 
include { FAST_AMRplusplus } from './subworkflows/AMR++_fast.nf'
include { STANDARD_AMRplusplus_wKraken } from './subworkflows/AMR++_standard_wKraken.nf'

// Load subworkflows
include { FASTQ_QC_WF } from './subworkflows/fastq_information.nf'
include { FASTQ_TRIM_WF } from './subworkflows/fastq_QC_trimming.nf'
include { FASTQ_RM_HOST_WF } from './subworkflows/fastq_host_removal.nf' 
include { FASTQ_RESISTOME_WF } from './subworkflows/fastq_resistome.nf'
include { FASTQ_KRAKEN_WF } from './subworkflows/fastq_microbiome.nf'


workflow {
    
    if (params.pipeline == "demo") {

        //run with demo params, use params.config
        FAST_AMRplusplus(fastq_files, params.amr, params.annotation)
        
    } else if(params.pipeline == "standard_AMR") {

        STANDARD_AMRplusplus(fastq_files,params.reference, params.amr, params.annotation)
        
    } else if(params.pipeline == "fast_AMR") {

        FAST_AMRplusplus(fastq_files, params.amr, params.annotation)
    } 
    else if(params.pipeline == "standard_AMR_wKraken") {

        STANDARD_AMRplusplus_wKraken(fastq_files,params.reference, params.amr, params.annotation, params.kraken_db)
    } 
    else if(params.pipeline == "multiqc") {

        FASTQ_QC_WF( fastq_files )
    } 
    else {
            println "ERROR ################################################################"
            println "Please choose a pipeline!!!" 
            println ""
            println "To test the pipeline, use the \"demo\" pipeline :"
            println ""
            println "ERROR ################################################################"
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