#!/usr/bin/env nextflow

/*
 * Another attempt at a nextflow script 2020-01-15
 */

// Declare Inputs
refFolder = file("refFiles")
inputDirectory = file('rawFiles/catFastq')

// Declare References 
refBase = "$refFolder/genome"
ref     = "${refBase}.fa"

// Create channel stream
Channel.fromFilePairs("rawFiles/catFastq/*_R{1,2}_001.fastq.gz")
  .set{ ch_fastaIn }

process align_bwa {

   input:
     set sampName, file(fastqs) from ch_fastaIn

   output:
     set sampName, file("${sampName}.sam") into ch_mappedBams

   publishDir path: './out_bam', mode: 'copy'

   script:
   """
   bwa mem -t ${task.cpus} -R "@RG\tID:${sampName}\tPU:${sampName}\tSM:${sampName}\tPL:ILLUMINA\tLB:rhAmpSeq" \
       $ref ${fastqs[0]} ${fastqs[1]} > ${sampName}.sam
   """
}
