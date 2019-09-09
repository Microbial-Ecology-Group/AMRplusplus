#!/usr/bin/env nextflow

/*
vim: syntax=groovy
-*- mode: groovy;-*-
*/

if (params.help ) {
    return help()
}
if( params.host_index ) {
    host_index = Channel.fromPath(params.host_index).toSortedList()
    //if( host_index.isEmpty() ) return index_error(host_index)
}
if( params.host ) {
    host = file(params.host)
    if( !host.exists() ) return host_error(host)
}
if( params.amr ) {
    amr = file(params.amr)
    if( !amr.exists() ) return amr_error(amr)
}
if( params.adapters ) {
    adapters = file(params.adapters)
    if( !adapters.exists() ) return adapter_error(adapters)
}
if( params.annotation ) {
    annotation = file(params.annotation)
    if( !annotation.exists() ) return annotation_error(annotation)
}

kraken_db = params.kraken_db
threads = params.threads

threshold = params.threshold

min = params.min
max = params.max
skip = params.skip
samples = params.samples

leading = params.leading
trailing = params.trailing
slidingwindow = params.slidingwindow
minlen = params.minlen

Channel
    .fromFilePairs( params.reads, flat: true )
    .ifEmpty { exit 1, "Read pair files could not be found: ${params.reads}" }
    .set { reads }

process RunQC {
    tag { sample_id }

    publishDir "${params.output}/RunQC", mode: 'copy', pattern: '*.fastq',
        saveAs: { filename ->
            if(filename.indexOf("P.fastq") > 0) "Paired/$filename"
            else if(filename.indexOf("U.fastq") > 0) "Unpaired/$filename"
            else {}
        }

    input:
        set sample_id, file(forward), file(reverse) from reads

    output:
        set sample_id, file("${sample_id}.1P.fastq.gz"), file("${sample_id}.2P.fastq.gz") into (paired_fastq)
        set sample_id, file("${sample_id}.1U.fastq.gz"), file("${sample_id}.2U.fastq.gz") into (unpaired_fastq)
        file("${sample_id}.trimmomatic.stats.log") into (trimmomatic_stats)

    """
     ${JAVA} -jar ${TRIMMOMATIC}/trimmomatic.jar \
      PE \
      -threads ${threads} \
      $forward $reverse ${sample_id}.1P.fastq ${sample_id}.1U.fastq ${sample_id}.2P.fastq ${sample_id}.2U.fastq \
      ILLUMINACLIP:${adapters}:2:30:10:3:TRUE \
      LEADING:${leading} \
      TRAILING:${trailing} \
      SLIDINGWINDOW:${slidingwindow} \
      MINLEN:${minlen} \
      2> ${sample_id}.trimmomatic.stats.log

      gzip *fastq
    """
}

trimmomatic_stats.toSortedList().set { trim_stats }

process QCStats {
    tag { sample_id }

    publishDir "${params.output}/RunQC", mode: 'copy',
        saveAs: { filename ->
            if(filename.indexOf(".stats") > 0) "Stats/$filename"
            else {}
        }

    input:
        file(stats) from trim_stats

    output:
	file("trimmomatic.stats")

    """
    ${PYTHON3} $baseDir/bin/trimmomatic_stats.py -i ${stats} -o trimmomatic.stats
    """
}

if( !params.host_index ) {
    process BuildHostIndex {
        publishDir "${params.output}/BuildHostIndex", mode: "copy"

        tag { host.baseName }

        input:
            file(host)

        output:
            file '*' into (host_index)

        """
        ${BWA} index ${host}
        """
    }
}

process AlignReadsToHost {
    tag { sample_id }

    publishDir "${params.output}/AlignReadsToHost", mode: "copy"

    input:
        set sample_id, file(forward), file(reverse) from paired_fastq
        file index from host_index
        file host

    output:
        set sample_id, file("${sample_id}.host.sam") into (host_sam)

    """
    ${BWA} mem ${host} ${forward} ${reverse} -t ${threads} > ${sample_id}.host.sam
    """
}

process RemoveHostDNA {
    tag { sample_id }

    publishDir "${params.output}/RemoveHostDNA", mode: "copy", pattern: '*.bam',
	saveAs: { filename ->
            if(filename.indexOf(".bam") > 0) "NonHostBAM/$filename"
        }

    input:
        set sample_id, file(sam) from host_sam

    output:
        set sample_id, file("${sample_id}.host.sorted.removed.bam") into (non_host_bam)
        file("${sample_id}.samtools.idxstats") into (idxstats_logs)

    """
    ${SAMTOOLS} view -bS ${sam} | ${SAMTOOLS} sort -@ ${threads} -o ${sample_id}.host.sorted.bam
    ${SAMTOOLS} index ${sample_id}.host.sorted.bam && ${SAMTOOLS} idxstats ${sample_id}.host.sorted.bam > ${sample_id}.samtools.idxstats
    ${SAMTOOLS} view -h -f 4 -b ${sample_id}.host.sorted.bam -o ${sample_id}.host.sorted.removed.bam
    """
}

idxstats_logs.toSortedList().set { host_removal_stats }

process HostRemovalStats {
    tag { sample_id }

    publishDir "${params.output}/RemoveHostDNA", mode: "copy",
        saveAs: { filename ->
            if(filename.indexOf(".stats") > 0) "HostRemovalStats/$filename"
        }

    input:
        file(stats) from host_removal_stats

    output:
        file("host.removal.stats")

    """
    ${PYTHON3} $baseDir/bin/samtools_idxstats.py -i ${stats} -o host.removal.stats
    """
}

process BAMToFASTQ {
    tag { sample_id }

    publishDir "${params.output}/BAMToFASTQ", mode: "copy"

    input:
        set sample_id, file(bam) from non_host_bam

    output:
        set sample_id, file("${sample_id}.non.host.R1.fastq"), file("${sample_id}.non.host.R2.fastq") into (non_host_fastq, non_host_fastq_kraken)

    """
    ${BEDTOOLS}  \
       bamtofastq \
      -i ${bam} \
      -fq ${sample_id}.non.host.R1.fastq \
      -fq2 ${sample_id}.non.host.R2.fastq
    """
}

if( !params.amr_index ) {
    process BuildAMRIndex {
        tag { amr.baseName }

        input:
            file(amr)

        output:
            file '*' into (amr_index)

        """
        ${BWA} index ${amr}
        """
    }
}

process AlignToAMR {
     tag { sample_id }

     publishDir "${params.output}/AlignToAMR", mode: "copy"

     input:
         set sample_id, file(forward), file(reverse) from non_host_fastq
         file index from amr_index
         file amr

     output:
         set sample_id, file("${sample_id}.amr.alignment.sam") into (resistome_sam, rarefaction_sam, snp_sam)

     """
     ${BWA} mem ${amr} ${forward} ${reverse} -t ${threads} > ${sample_id}.amr.alignment.sam
     """
}

process RunResistome {
    tag { sample_id }

    publishDir "${params.output}/RunResistome", mode: "copy"

    input:
        set sample_id, file(sam) from resistome_sam
        file annotation
        file amr

    output:
        file("${sample_id}.gene.tsv") into (resistome)

    """
    ${RESISTOME} -ref_fp ${amr} \
      -annot_fp ${annotation} \
      -sam_fp ${sam} \
      -gene_fp ${sample_id}.gene.tsv \
      -group_fp ${sample_id}.group.tsv \
      -class_fp ${sample_id}.class.tsv \
      -mech_fp ${sample_id}.mechanism.tsv \
      -t ${threshold}
    """
}

process RunRarefaction {
    tag { sample_id }

    publishDir "${params.output}/RunRarefaction", mode: "copy"

    input:
        set sample_id, file(sam) from rarefaction_sam
        file annotation
        file amr

    output:
        set sample_id, file("*.tsv") into (rarefaction)

    """
    ${RAREFACTION} \
      -ref_fp ${amr} \
      -sam_fp ${sam} \
      -annot_fp ${annotation} \
      -gene_fp ${sample_id}.gene.tsv \
      -group_fp ${sample_id}.group.tsv \
      -class_fp ${sample_id}.class.tsv \
      -mech_fp ${sample_id}.mech.tsv \
      -min ${min} \
      -max ${max} \
      -skip ${skip} \
      -samples ${samples} \
      -t ${threshold}
    """
}

process RunSNPFinder {
    tag { sample_id }

    publishDir "${params.output}/RunSNPFinder", mode: "copy"

    input:
        set sample_id, file(sam) from snp_sam
        file amr

    output:
        set sample_id, file("*.tsv") into (snp)

    """
    ${SNPFINDER} \
      -amr_fp ${amr} \
      -sampe ${sam} \
      -out_fp ${sample_id}.tsv
    """
}

process RunKraken {
    tag { sample_id }

    publishDir "${params.output}/RunKraken", mode: "copy"

    input:
       set sample_id, file(forward), file(reverse) from non_host_fastq_kraken

    output:
       file("${sample_id}.kraken.filtered.report") into
       kraken_report

    """
    ${KRAKEN2} --preload --db ${kraken_db} --paired ${forward} ${reverse} --threads ${threads} --report ${sample_id}.kraken.report > ${sample_id}.kraken.raw
    ${KRAKEN2} --preload --db ${kraken_db} --confidence 1 --paired ${forward} ${reverse} --threads ${threads} --report ${sample_id}.kraken.filtered.report > ${sample_id}.kraken.raw
    """
}

resistome.toSortedList().set { amr_l_to_w }

process AMRLongToWide {
    tag { }

    publishDir "${params.output}/AMRLongToWide", mode: "copy"

    input:
        file(resistomes) from amr_l_to_w

    output:
        file("AMR_analytic_matrix.csv") into amr_master_matrix

    """
    mkdir ret
    ${PYTHON3} $baseDir/bin/amr_long_to_wide.py -i ${resistomes} -o ret
    mv ret/AMR_analytic_matrix.csv .
    """
}

kraken_report.toSortedList().set { kraken_l_to_w }

process KrakenLongToWide {
    tag { }

    publishDir "${params.output}/KrakenLongToWide", mode: "copy"

    input:
        file(kraken_reports) from kraken_l_to_w

    output:
        file("kraken_analytic_matrix.csv") into kraken_master_matrix

    """
    mkdir ret
    ${PYTHON3} $baseDir/bin/kraken2_long_to_wide.py -i ${kraken_reports} -o ret
    mv ret/kraken_analytic_matrix.csv .
    """
}

def nextflow_version_error() {
    println ""
    println "This workflow requires Nextflow version 0.25 or greater -- You are running version $nextflow.version"
    println "Run ./nextflow self-update to update Nextflow to the latest available version."
    println ""
    return 1
}

def adapter_error(def input) {
    println ""
    println "[params.adapters] fail to open: '" + input + "' : No such file or directory"
    println ""
    return 1
}

def amr_error(def input) {
    println ""
    println "[params.amr] fail to open: '" + input + "' : No such file or directory"
    println ""
    return 1
}

def annotation_error(def input) {
    println ""
    println "[params.annotation] fail to open: '" + input + "' : No such file or directory"
    println ""
    return 1
}

def fastq_error(def input) {
    println ""
    println "[params.reads] fail to open: '" + input + "' : No such file or directory"
    println ""
    return 1
}

def host_error(def input) {
    println ""
    println "[params.host] fail to open: '" + input + "' : No such file or directory"
    println ""
    return 1
}

def index_error(def input) {
    println ""
    println "[params.host_index] fail to open: '" + input + "' : No such file or directory"
    println ""
    return 1
}

def help() {
    println ""
    println "Program: AmrPlusPlus"
    println "Documentation: https://github.com/colostatemeg/amrplusplus/blob/master/README.md"
    println "Contact: Christopher Dean <cdean11@colostate.edu>"
    println ""
    println "Usage:    nextflow run main.nf [options]"
    println ""
    println "Input/output options:"
    println ""
    println "    --reads         STR      path to FASTQ formatted input sequences"
    println "    --adapters      STR      path to FASTA formatted adapter sequences"
    println "    --host          STR      path to FASTA formatted host genome"
    println "    --host_index    STR      path to BWA generated index files"
    println "    --amr           STR      path to AMR resistance database"
    println "    --annotation    STR      path to AMR annotation file"
    println "    --output        STR      directory to write process outputs to"
    println "    --KRAKENDB      STR      path to kraken database"
    println ""
    println "Trimming options:"
    println ""
    println "    --leading       INT      cut bases off the start of a read, if below a threshold quality"
    println "    --minlen        INT      drop the read if it is below a specified length"
    println "    --slidingwindow INT      perform sw trimming, cutting once the average quality within the window falls below a threshold"
    println "    --trailing      INT      cut bases off the end of a read, if below a threshold quality"
    println ""
    println "Algorithm options:"
    println ""
    println "    --threads       INT      number of threads to use for each process"
    println "    --threshold     INT      gene fraction threshold"
    println "    --min           INT      starting sample level"
    println "    --max           INT      ending sample level"
    println "    --samples       INT      number of sampling iterations to perform"
    println "    --skip          INT      number of levels to skip"
    println ""
    println "Help options:"
    println ""
    println "    --help                   display this message"
    println ""
    return 1
}
