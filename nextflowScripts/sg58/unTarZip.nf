#!/usr/bin/env nextflow

/*
 * This is the nextflow script to unpack downloaded TCGA files and then
 * repackage them into fastq.gz files for downstream RNA Seq processing
 */

// Declare Inputs

// Declare References 

// Tools
pigzModule = "pigz/2.3.4"

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

    label 'vardict_comp'

    input:
        tuple file(fq1), file(fq2) from ch_fastq
    output:
        file("*.fastq.gz") into ch_gz

    publishDir path: './fastqs', mode: 'copy'

    module pigzModule

    script:
    """
    pigz -f -p ${task.cpus} ${fq1}
    pigz -f -p ${task.cpus} ${fq2}
    """
}