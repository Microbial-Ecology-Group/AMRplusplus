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


// Load null pipeline
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
include { FASTQ_QIIME2_WF } from './subworkflows/fastq_16S_qiime2.nf'

// Load BAM subworkflows
include { BAM_RESISTOME_WF } from './subworkflows/bam_resistome.nf'



workflow {
    if (params.pipeline == null || params.pipeline == "demo") {
        log.info """\
        Running a demonstration of AMR++
        ===================================
        To include SNP analysis, add `--snp Y` to your command.
        ===================================        
        To include deduplicated count analysis, add `--deduped Y` to your command. 
        Please be aware that adding deduplicated counts will significantly increase run time and temp file storage requirements.
        ===================================
        """
        //run with demo params, use params.config
        FAST_AMRplusplus(fastq_files, params.amr, params.annotation)
        
    } else if(params.pipeline == "standard_AMR") {

        STANDARD_AMRplusplus(fastq_files,params.host, params.amr, params.annotation)
        
    } else if(params.pipeline == "fast_AMR") {

        FAST_AMRplusplus(fastq_files, params.amr, params.annotation)
    } 
    else if(params.pipeline == "standard_AMR_wKraken") {

        STANDARD_AMRplusplus_wKraken(fastq_files,params.host, params.amr, params.annotation, params.kraken_db)
    } 
    else if(params.pipeline == "eval_qc") {

        FASTQ_QC_WF( fastq_files )
    } 
    else if(params.pipeline == "trim_qc") {

        FASTQ_TRIM_WF( fastq_files )
    }
    else if(params.pipeline == "rm_host") {

        FASTQ_RM_HOST_WF(params.host, fastq_files )
    } 
    else if(params.pipeline == "resistome") {

        FASTQ_RESISTOME_WF( fastq_files, params.amr, params.annotation )
    }  
    else if(params.pipeline == "kraken") {

        FASTQ_KRAKEN_WF( fastq_files , params.kraken_db)
    }
    else if(params.pipeline == "qiime2") {
        Channel
            .fromFilePairs( params.reads, flat: true )
            .ifEmpty { exit 1, "Read pair files could not be found: ${params.reads}" }
            .map { name, forward, reverse -> [ forward.drop(forward.findLastIndexOf{"/"})[0], forward, reverse ] } //extract file name
            .map { name, forward, reverse -> [ name.toString().take(name.toString().indexOf("_")), forward, reverse ] } //extract sample name
            .map { name, forward, reverse -> [ name +","+ forward + ",forward\n" + name +","+ reverse +",reverse" ] } //prepare basic synthax
            .flatten()
            .collectFile(name: 'manifest.txt', newLine: true, storeDir: "${params.output}/demux", seed: "sample-id,absolute-filepath,direction")
            .set { ch_manifest }
        
        FASTQ_QIIME2_WF( ch_manifest , params.dada2_db)
    }
    else if(params.pipeline == "bam_resistome"){
        Channel
        .fromPath(params.bam_files)
        .map { file -> [ id:file.baseName,file:file] }
        .set {bam_files_ch}
        
        BAM_RESISTOME_WF( bam_files_ch , params.amr, params.annotation )

    }

    else {
            println "ERROR ################################################################"
            println "Please choose a pipeline!!!" 
            println ""
            println "To test the pipeline, use the \"demo\" pipeline or omit the pipeline flag:"
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
