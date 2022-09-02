nextflow.enable.dsl=2
// Example command:
// nextflow run main_amr++.nf -profile standard --multiqc "/home/enrique/Dropbox/Projects/AMR++_v3_update/nfcamp-tutorial/multiqc/" -resume

/*
 * Default pipeline parameters. They can be overriden on the command line eg.
 * given `params.foo` specify on the run command line `--foo some_value`.
 */

//params.reads = "$baseDir/data/ggal/ggal_gut_{1,2}.fq"
// nextflow run main_amr++.nf -profile local

log.info """\
 A M R + +    N F   P I P E L I N E
 ===================================
 reads        : ${params.reads}
 output       : ${params.output}
 """

// Load main pipeline workflows

include { STANDARD_AMRplusplus } from './subworkflows/AMR++_standard.nf' 
include { FAST_AMRplusplus } from './subworkflows/AMR++_fast.nf'

// Load subworkflows
include { FASTQ_QC_WF } from './subworkflows/fastq_information.nf'
include { FASTQ_TRIM_WF } from './subworkflows/fastq_QC_trimming.nf'
include { FASTQ_RM_HOST_WF } from './subworkflows/fastq_host_removal.nf' 
include { FASTQ_RESISTOME_WF } from './subworkflows/fastq_resistome.nf'

Channel
.fromFilePairs( params.reads , size: ("${params.reads}" =~ /\{/) ? 2 : 1)
.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
.set { fastq_files }


workflow {
    /// This works for individual processes
    //FASTQC_QC_WF(fastq_files)
    //FASTQ_TRIM_WF(fastq_files)
    //FASTQ_RM_HOST_WF(params.reference , TRIMMING_WF.out.trimmed_reads)
    //FASTQ_RESISTOME_WF(fastq_files, params.amr, params.annotation)


    /// STD_AMRplusplus workflow is the standard analysis which includes
    // fastqc, QC trimming, host removal, followed by alignment to MEGARes and SNP confirmation
    
    STANDARD_AMR++(fastq_files,params.reference, params.amr, params.annotation)
    
    //FAST_AMR++(fastq_files, params.amr, params.annotation)
}

//workflow {
//    if (params.pipeline != "" && params.fastq == "") {
        //run with demo params, use params.config     
//    } else if(params.pipeline == "standard_AMR" && params.fastq != "") {

//    STANDARD_AMR++(fastq_files,params.reference, params.amr, params.annotation)

        
//    } else {
  //          println "ERROR ################################################################"
    //        println "Please choose one between fast5 and fastq as input!!!" 
      //      println "ERROR ################################################################"
        //    println "Exiting ..."
          //  System.exit(0)  
    //}
//}

workflow.onComplete {
    println "Pipeline completed!"
    println "Started at  $workflow.start"
    println "Finished at $workflow.complete"
    println "Time elapsed: $workflow.duration"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}