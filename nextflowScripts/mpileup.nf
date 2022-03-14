#!/usr/bin/env nextflow

//* This script is to create an mpileup according to a BED file of locations of interest
//* The file will be placed into the refFolder

// Declare Inputs
refFolder = file("/home/nwong/ls25_scratch/nick.wong/refFiles")
refFile = "${refFolder}/loi.bed"

// Load modules
samtoolsModule = 'samtools/1.9-gcc5'

// Create channel stream
Channel.fromPath("fc*/bams/*.bam")
  .set{ ch_bamIn }

process mpileup {
   
   label 'fastqc'


   input:
     file(bams) from ch_bamIn

   output:
     file("${bams.getSimpleName()}.mpileup") into ch_outMpileup


   publishDir path: './mpileup', mode: 'copy'

   script:
   """
   
   samtools mpileup  \
      -l /scratch/ls25/nick.wong/refFiles/loi.bed \
      -f /scratch/ls25/nick.wong/refFiles/genome.fa \
      -o ${bams.getSimpleName()}.mpileup \
      -d 0 \
      ${bams}
   """
}
