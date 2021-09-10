#!/usr/bin/env nextflow

params.fastq  = '*.fastq.gz'
reads = Channel
                .fromPath(params.fastq)
                .map { file -> tuple(file.baseName, file) }

reads.into { reads_qc; reads_fastp }

params.transcripts = 'gencode.v38.transcripts.fa'
transcriptome = Channel.fromPath(params.transcripts)


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
    set readID, file(readFile) from processed_reads

    output:
    file("read_counts/$read/quant.sf") into read_counts

    """
    mkdir read_counts
    salmon quant -i $idx -l A -r $read -o read_counts/$read
    """
}
