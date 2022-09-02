// Load modules
include { index as amr_index ; index as host_index } from '../modules/Alignment/bwa' addParams(EXTRAPARS: "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36")
include { bwa_align } from '../modules/Alignment/bwa' addParams(EXTRAPARS: "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36")
include { fastqc ; multiqc } from '../modules/Fastqc/fastqc' addParams(EXTRAPARS: "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36")
include { runqc } from '../modules/Trimming/trimmomatic' addParams(EXTRAPARS: "LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:36")

// contaminant removal
include { bwa_rm_contaminant_fq } from '../modules/Alignment/bwa' addParams(EXTRAPARS: "test")

// resistome
include { runresistome } from '../modules/Resistome/resistome' addParams(EXTRAPARS: "test")
include { runsnp } from '../modules/Resistome/resistome' addParams(EXTRAPARS: "test")
include { resistomeresults } from '../modules/Resistome/resistome' addParams(EXTRAPARS: "test")
include { runrarefaction } from '../modules/Resistome/resistome' addParams(EXTRAPARS: "test")


workflow STANDARD_AMRplusplus {
    take: 
        read_pairs_ch
        hostfasta
        amr
        annotation

    main:
        // fastqc
        fastqc( read_pairs_ch )
        multiqc(fastqc.out.collect(), params.multiqc )
        // runqc trimming
        runqc(read_pairs_ch)
        // make indices for host and amr
        amr_index(amr)
        host_index(hostfasta)
        // remove host DNA
        bwa_rm_contaminant_fq(hostfasta,host_index.out, runqc.out.paired_fastq )
        // AMR alignment
        bwa_align(amr, amr_index.out, bwa_rm_contaminant_fq.out.nonhost_reads )
        runresistome(bwa_align.out.bwa_sam,amr, annotation )
        //runsnp(bwa_align.out.bwa_sam )
        resistomeresults(runresistome.out.resistome_counts.collect())
        runrarefaction(bwa_align.out.bwa_sam, annotation, amr)

    emit:
        fastqc = fastqc.out   
        multiqc = multiqc.out
        trim_reads = runqc.out.paired_fastq
        non_host_reads = bwa_rm_contaminant_fq.out.nonhost_reads
}
