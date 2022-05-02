#!/usr/bin/env nextflow

/*
 * This is to run sequenza processing for samples and CNV analysis
 * 
 */

// Declare Inputs
refFolder = file("/projects/tn26/referenceFiles/human/ucschg19")

// Declare References 
refBase        = "${refFolder}/genome"
ref            = "${refBase}.fa"

// Tools
singularityModule = 'singularity/3.7.1'
bowtie2Module = 'bowtie2/2.2.9'

// Create channel stream
Channel.fromFilePairs("fastq_trimmed/*-pair{1,2}.fastq")
  .set{ ch_hlahd }
  
process hlahd {
     
    label 'hlahd'

    input:
        set sampName, file(vardict) from ch_hlahd
    output:
        set sampName, file("${sampName}*") into ch_hlahdDone

    publishDir path: './hlahd', mode: 'copy'

    module bowtie2Module

    script:
    """
    /home/nwong/rc78_scratch/angela.pizzola/hlahd.1.4.0/bin/hlahd.sh  \
      
    """
}
