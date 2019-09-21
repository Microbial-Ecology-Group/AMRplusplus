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
if(params.kraken_db) {
    kraken_db = file(params.kraken_db)
}


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
     ${JAVA} -jar ${TRIMMOMATIC} \
      PE \
      -threads ${threads} \
      $forward $reverse ${sample_id}.1P.fastq.gz ${sample_id}.1U.fastq.gz ${sample_id}.2P.fastq.gz ${sample_id}.2U.fastq.gz \
      ILLUMINACLIP:${adapters}:2:30:10:3:TRUE \
      LEADING:${leading} \
      TRAILING:${trailing} \
      SLIDINGWINDOW:${slidingwindow} \
      MINLEN:${minlen} \
      2> ${sample_id}.trimmomatic.stats.log
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

process NonHostReads {
    tag { sample_id }

    publishDir "${params.output}/NonHostReads", mode: "copy"

    input:
        set sample_id, file(bam) from non_host_bam

    output:
        set sample_id, file("${sample_id}.non.host.R1.fastq"), file("${sample_id}.non.host.R2.fastq") into (non_host_fastq_megares, non_host_fastq_dedup,non_host_fastq_kraken)

    """
    ${BEDTOOLS}  \
       bamtofastq \
      -i ${bam} \
      -fq ${sample_id}.non.host.R1.fastq \
      -fq2 ${sample_id}.non.host.R2.fastq
    """
}

/*
-
--
---
---- nonhost reads for megares and kraken2
---
--
-
*/


/*
---- Run Kraken2
*/



process RunKraken {
    tag { sample_id }

    publishDir "${params.output}/RunKraken", mode: 'copy',
        saveAs: { filename ->
            if(filename.indexOf(".kraken.raw") > 0) "Standard/$filename"
            else if(filename.indexOf(".kraken.report") > 0) "Standard_report/$filename"
            else if(filename.indexOf(".kraken.filtered.report") > 0) "Filtered_report/$filename"
            else if(filename.indexOf(".kraken.filtered.raw") > 0) "Filtered/$filename"
            else {}
        }

    input:
       set sample_id, file(forward), file(reverse) from non_host_fastq_kraken

   output:
      file("${sample_id}.kraken.report") into (kraken_report,kraken_extract_taxa)
      set sample_id, file("${sample_id}.kraken.raw") into kraken_raw
      file("${sample_id}.kraken.filtered.report") into kraken_filter_report
      file("${sample_id}.kraken.filtered.raw") into kraken_filter_raw

     """
     ${KRAKEN2} --db ${kraken_db} --paired ${forward} ${reverse} --threads ${threads} --report ${sample_id}.kraken.report > ${sample_id}.kraken.raw
     ${KRAKEN2} --db ${kraken_db} --confidence 1 --paired ${forward} ${reverse} --threads ${threads} --report ${sample_id}.kraken.filtered.report > ${sample_id}.kraken.filtered.raw
     """
}

kraken_report.toSortedList().set { kraken_l_to_w }
kraken_filter_report.toSortedList().set { kraken_filter_l_to_w }

process KrakenResults {
    tag { }

    publishDir "${params.output}/KrakenResults", mode: "copy"

    input:
        file(kraken_reports) from kraken_l_to_w

    output:
        file("kraken_analytic_matrix.csv") into kraken_master_matrix

    """
    ${PYTHON3} $baseDir/bin/kraken2_long_to_wide.py -i ${kraken_reports} -o kraken_analytic_matrix.csv
    """
}

process FilteredKrakenResults {
    tag { sample_id }

    publishDir "${params.output}/FilteredKrakenResults", mode: "copy"

    input:
        file(kraken_reports) from kraken_filter_l_to_w

    output:
        file("filtered_kraken_analytic_matrix.csv") into filter_kraken_master_matrix

    """
    ${PYTHON3} $baseDir/bin/kraken2_long_to_wide.py -i ${kraken_reports} -o filtered_kraken_analytic_matrix.csv
    """
}



/*
---- Run alignment to MEGAres
*/

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
         set sample_id, file(forward), file(reverse) from non_host_fastq_megares
         file index from amr_index
         file amr

     output:
         set sample_id, file("${sample_id}.amr.alignment.sam") into (megares_resistome_sam, megares_rarefaction_sam, megares_snp_sam , megares_snpfinder_sam, megares_RGI_sam)
         set sample_id, file("${sample_id}.amr.alignment.dedup.sam") into (megares_dedup_resistome_sam,megares_dedup_RGI_sam)
         set sample_id, file("${sample_id}.amr.alignment.dedup.bam") into (megares_dedup_resistome_bam)


     """
     ${BWA} mem ${amr} ${forward} ${reverse} -t ${threads} -R '@RG\\tID:${sample_id}\\tSM:${sample_id}' > ${sample_id}.amr.alignment.sam
     ${SAMTOOLS} view -S -b ${sample_id}.amr.alignment.sam > ${sample_id}.amr.alignment.bam
     ${SAMTOOLS} sort -n ${sample_id}.amr.alignment.bam -o ${sample_id}.amr.alignment.sorted.bam
     ${SAMTOOLS} fixmate ${sample_id}.amr.alignment.sorted.bam ${sample_id}.amr.alignment.sorted.fix.bam
     ${SAMTOOLS} sort ${sample_id}.amr.alignment.sorted.fix.bam -o ${sample_id}.amr.alignment.sorted.fix.sorted.bam
     ${SAMTOOLS} rmdup -S ${sample_id}.amr.alignment.sorted.fix.sorted.bam ${sample_id}.amr.alignment.dedup.bam
     ${SAMTOOLS} view -h -o ${sample_id}.amr.alignment.dedup.sam ${sample_id}.amr.alignment.dedup.bam
     #rm ${sample_id}.amr.alignment.bam
     #rm ${sample_id}.amr.alignment.sorted*.bam
     """
}

process RunResistome {
    tag { sample_id }

    publishDir "${params.output}/RunResistome", mode: "copy"

    input:
        set sample_id, file(sam) from megares_resistome_sam
        file annotation
        file amr

    output:
        file("${sample_id}.gene.tsv") into (megares_resistome_counts, SNP_confirm_long)

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

megares_resistome_counts.toSortedList().set { megares_amr_l_to_w }

process ResistomeResults {
    tag { }

    publishDir "${params.output}/ResistomeResults", mode: "copy"

    input:
        file(resistomes) from megares_amr_l_to_w

    output:
        file("AMR_analytic_matrix.csv") into amr_master_matrix

    """
    ${PYTHON3} $baseDir/bin/amr_long_to_wide.py -i ${resistomes} -o AMR_analytic_matrix.csv
    """
}


/* samtools deduplication of megares alignment */
process SamDedupRunResistome {
    tag { sample_id }

    publishDir "${params.output}/SamDedupRunResistome", mode: "copy"

    input:
        set sample_id, file(sam) from megares_dedup_resistome_sam
        file annotation
        file amr

    output:
        file("${sample_id}.gene.tsv") into (megares_dedup_resistome_counts)

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

megares_dedup_resistome_counts.toSortedList().set { megares_dedup_amr_l_to_w }

process SamDedupResistomeResults {
    tag { }

    publishDir "${params.output}/SamDedup_ResistomeResults", mode: "copy"

    input:
        file(resistomes) from megares_dedup_amr_l_to_w

    output:
        file("SamDedup_AMR_analytic_matrix.csv") into megares_dedup_amr_master_matrix

    """
    ${PYTHON3} $baseDir/bin/amr_long_to_wide.py -i ${resistomes} -o SamDedup_AMR_analytic_matrix.csv
    """
}

process RunRarefaction {
    tag { sample_id }

    publishDir "${params.output}/RunRarefaction", mode: "copy"

    input:
        set sample_id, file(sam) from megares_rarefaction_sam
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



/*
---- Confirmation of alignments to genes that require SNP confirmation with RGI
*/

process ExtractSNP {
     tag { sample_id }

     publishDir "${params.output}/ExtractMegaresSNPs", mode: "copy",
         saveAs: { filename ->
             if(filename.indexOf(".snp.fasta") > 0) "SNP_fasta/$filename"
             else if(filename.indexOf("gene.tsv") > 0) "Gene_hits/$filename"
             else {}
         }

     input:
         set sample_id, file(sam) from megares_RGI_sam
         file annotation
         file amr

     output:
         set sample_id, file("*.snp.fasta") into megares_snp_fasta
         set sample_id, file("${sample_id}*.gene.tsv") into (resistome_hits)

     """
     awk -F "\\t" '{if (\$1!="@SQ" && \$1!="@RG" && \$1!="@PG" && \$1!="@HD" && \$3="RequiresSNPConfirmation" ) {print ">"\$1"\\n"\$10}}' ${sam} | tr -d '"'  > ${sample_id}.snp.fasta
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

process RunRGI {
     tag { sample_id }
     errorStrategy 'ignore'
	
     publishDir "${params.output}/RunRGI", mode: "copy"

     input:
         set sample_id, file(fasta) from megares_snp_fasta

     output:
         set sample_id, file("${sample_id}_rgi_output.txt") into rgi_results

     """
     alias diamond='echo "${DIAMOND}"'
     cp ${fasta} ${fasta}.temp.contig.fsa
     ${RGI} main --low_quality --input_sequence ${fasta} --output_file ${sample_id}_rgi_output -a diamond -n ${threads} --clean
     """
}


process SNPconfirmation {
     tag { sample_id }
     errorStrategy 'ignore'

     publishDir "${params.output}/SNPConfirmation", mode: "copy",
         saveAs: { filename ->
             if(filename.indexOf("_rgi_perfect_hits.csv") > 0) "Perfect_RGI/$filename"
             else if(filename.indexOf("_rgi_strict_hits.csv") > 0) "Strict_RGI/$filename"
             else if(filename.indexOf("_rgi_loose_hits.csv") > 0) "Loose_RGI/$filename"
             else {}
         }

     input:
         set sample_id, file(rgi) from rgi_results

     output:
         set sample_id, file("${sample_id}_rgi_strict_hits.csv") into strict_snp_long_hits
     """
     ${PYTHON3} $baseDir/bin/RGI_aro_hits.py ${rgi} ${sample_id}
     """
}

process Confirmed_AMR_hits {
     tag { sample_id }

     publishDir "${params.output}/SNP_confirmed_counts", mode: "copy"

     input:
         set sample_id, file(megares_counts) from resistome_hits
         set sample_id, file(strict_rgi_counts) from strict_snp_long_hits

     output:
         file("${sample_id}*strict_SNP_confirmed_counts") into strict_confirmed_counts

     """
     ${PYTHON3} $baseDir/bin/RGI_long_combine.py ${strict_rgi_counts} ${megares_counts} ${sample_id}.strict_SNP_confirmed_counts ${sample_id}
     """
}


strict_confirmed_counts.toSortedList().set { strict_confirmed_amr_l_to_w }

process Confirmed_ResistomeResults {
     tag {}

     publishDir "${params.output}/Confirmed_ResistomeResults", mode: "copy"

     input:
         file(strict_confirmed_resistomes) from strict_confirmed_amr_l_to_w

     output:
         file("strict_SNP_confirmed_AMR_analytic_matrix.csv") into strict_confirmed_matrix

     """
     ${PYTHON3} $baseDir/bin/amr_long_to_wide.py -i ${strict_confirmed_resistomes} -o strict_SNP_confirmed_AMR_analytic_matrix.csv
     """
}

/*
---- Confirmation of deduped alignments to genes that require SNP confirmation with RGI.
*/


process ExtractDedupSNP {
     tag { sample_id }

     publishDir "${params.output}/ExtractDedupMegaresSNPs", mode: "copy",
         saveAs: { filename ->
             if(filename.indexOf(".snp.fasta") > 0) "SNP_fasta/$filename"
             else if(filename.indexOf("gene.tsv") > 0) "Gene_hits/$filename"
             else {}
         }

     input:
         set sample_id, file(sam) from megares_dedup_RGI_sam
         file annotation
         file amr

     output:
         set sample_id, file("*.snp.fasta") into dedup_megares_snp_fasta
         set sample_id, file("${sample_id}*.gene.tsv") into (dedup_resistome_hits)

     """
     awk -F "\\t" '{if (\$1!="@SQ" && \$1!="@RG" && \$1!="@PG" && \$1!="@HD" && \$3="RequiresSNPConfirmation" ) {print ">"\$1"\\n"\$10}}' ${sam} | tr -d '"'  > ${sample_id}.snp.fasta
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

process RunDedupRGI {
     tag { sample_id }
     errorStrategy 'ignore'
     publishDir "${params.output}/RunDedupRGI", mode: "copy"

     input:
         set sample_id, file(fasta) from dedup_megares_snp_fasta

     output:
         set sample_id, file("${sample_id}_rgi_output.txt") into dedup_rgi_results

     """
     alias diamond='echo "${DIAMOND}"'
     cp ${fasta} ${fasta}.temp.contig.fsa
     ${RGI} main --low_quality --input_sequence ${fasta} --output_file ${sample_id}_rgi_output -a diamond -n ${threads} --clean
     """
}


process DedupSNPconfirmation {
     tag { sample_id }
     errorStrategy 'ignore'

     publishDir "${params.output}/DedupSNPConfirmation", mode: "copy",
         saveAs: { filename ->
             if(filename.indexOf("_rgi_perfect_hits.csv") > 0) "Perfect_RGI/$filename"
             else if(filename.indexOf("_rgi_strict_hits.csv") > 0) "Strict_RGI/$filename"
             else if(filename.indexOf("_rgi_loose_hits.csv") > 0) "Loose_RGI/$filename"
             else {}
         }

     input:
         set sample_id, file(rgi) from dedup_rgi_results

     output:
         set sample_id, file("${sample_id}_rgi_strict_hits.csv") into dedup_strict_snp_long_hits
     """
     ${PYTHON3} $baseDir/bin/RGI_aro_hits.py ${rgi} ${sample_id}
     """
}

process ConfirmDedupAMRHits {
     tag { sample_id }

     publishDir "${params.output}/SNP_confirmed_counts", mode: "copy"

     input:
         set sample_id, file(megares_counts) from dedup_resistome_hits
         set sample_id, file(strict_rgi_counts) from dedup_strict_snp_long_hits

     output:
         file("${sample_id}*strict_SNP_confirmed_counts") into dedup_strict_confirmed_counts

     """
     ${PYTHON3} $baseDir/bin/RGI_long_combine.py ${strict_rgi_counts} ${megares_counts} ${sample_id}.strict_SNP_confirmed_counts ${sample_id}
     """
}


dedup_strict_confirmed_counts.toSortedList().set { dedup_strict_confirmed_amr_l_to_w }

process DedupSNPConfirmed_ResistomeResults {
     tag {}

     publishDir "${params.output}/Confirmed_ResistomeResults", mode: "copy"

     input:
         file(strict_confirmed_resistomes) from dedup_strict_confirmed_amr_l_to_w

     output:
         file("strict_SNP_confirmed_AMR_analytic_matrix.csv") into dedup_strict_confirmed_matrix

     """
     ${PYTHON3} $baseDir/bin/amr_long_to_wide.py -i ${strict_confirmed_resistomes} -o strict_SNP_confirmed_dedup_AMR_analytic_matrix.csv
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
