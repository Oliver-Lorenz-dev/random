#!/usr/bin/env nextflow

params.fastq  = "$baseDir/nftest/*.fastq.gz"
reads = Channel
                .fromPath(params.fastq)
                .map { file -> tuple(file.baseName, file) }

reads.into { reads_qc; reads_fastp }

params.transcripts = "$baseDir/nftest/gencode.v38.transcripts.fa.gz"
transcriptome = Channel.fromPath(params.transcripts)
params.outdir = "rna-seq-results"

process qc {

    input:
    set readID, file(readFile) from reads_qc

    output:
    file("fastqc_${readID}_logs") into fastqc_ch

    script:
    """
    mkdir fastqc_${readID}_logs
    fastqc -o fastqc_${readID}_logs -f fastq -q ${readFile}
    """
}

process preprocessing {

    input:
    set readID, file(readFile) from reads_fastp

    output:
    file("fastpout/${readID}") into processed_reads

    """
    mkdir fastpout
    fastp -i ${readFile} -o fastpout/${readID}
    """
}

process create_index {

    input:
    file x from transcriptome

    output:
    file("index") into index

    """
    salmon index -t $x -i index --gencode
    """
}

process quantification {

    input:
    file idx from index
    file read from processed_reads

    output:
    file("read_counts/${read}/quant.sf")

    """
    mkdir read_counts
    salmon quant -i $idx -l A -r ${read} -o read_counts/${read}
    """
}

process multiqc {
    publishDir params.outdir, mode:'copy'

    input:
    file '*' from fastqc_ch.collect()

    output:
    file 'multiqc_report.html'

    script:
    """
    multiqc .
    """
}
