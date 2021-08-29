#!/usr/bin/env nextflow

/*
 * This is the nextflow script to unpack downloaded TCGA files and then
 * repackage them into fastq.gz files for downstream RNA Seq processing
 */

// Declare Inputs

// Declare References 

// Tools

// Create channel stream
Channel.fromPath("*/*.tar.gz")
  .set{ ch_folder }
  
process unTar {

    label 'untar'

    input:
        file(tarball) from ch_folder
    output:
        file("*.fastq") into ch_fastq

    script:
    """
    tar -xvzf ${tarball}
    """
}


process gzip {

    label 'untar'

    input:
        tuple file(fq1), file(fq2) from ch_fastq
    output:
        file("*.fastq.gz") into ch_gz

    publishDir path: './fastqs', mode: 'copy'

    script:
    """
    gzip ${fq1}
    """
}