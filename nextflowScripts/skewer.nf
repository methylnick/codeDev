#!/usr/bin/env nextflow

/*
 * Another attempt at a nextflow script 2020-01-18
 * This is a test and dev script so it is going straight to the short queue
 */

// Declare Inputs
refFolder = file("/home/nwong/ls25_scratch/nick.wong/refFiles")
inputDirectory = file('fastq')
AF_THR          = 0.0025

// Declare References 
refBase        = "${refFolder}/genome"
ref            = "${refBase}.fa"
target         = "${refFolder}/chip_v2_amplicon_bed_target_20200324.bed"
picardInsert   = "${refFolder}/chip_v2_amplicon_bed_target_20200324.bed.interval"
picardAmplicon = "${refFolder}/chip_v2_amplicon_bed_amplicon_20200324.bed.interval"
vardictAmp     = "${refFolder}/chip_v2_amplicon_bed_vardict_20200324.bed"

// Tools
picardModule   = 'picard/2.9.2'
bwaModule      = 'bwa/0.7.17-gcc5'
samtoolsModule = 'samtools/1.9-gcc5'
bedtoolsModule = 'bedtools/2.27.1-gcc5'
bcftoolsModule = 'bcftools/1.8'
RModule        = 'R/3.6.0-mkl'
fastqcModule         = 'fastqc/0.11.7'
skewerModule	= 'skewer/20170212'

// Create channel stream
Channel.fromFilePairs("fastq/*_R{1,2}_001.fastq.gz")
  .set{ ch_fastaIn }
  
Channel.fromFilePairs("fastq/*_R{1,2}_001.fastq.gz")
  .set{ ch_fastQC }

process skewer {
   
   label 'fastqc'

   input:
     set sampName, file(fastqs) from ch_fastQC

   output:
     set sampName, file("*-trimmed-pair1.fastq.gz"), file("*-trimmed-pair2.fastq.gz"), file("*-trimmed.log") into ch_outSkewer


   publishDir path: './trimmedFastq', mode: 'copy'

    module      skewerModule

   script:
   """
   skewer -t ${task.cpus} -z  ${fastqs[0]} ${fastqs[1]}
   """
}

process fastqc {
   
   label 'fastqc'

   input:
     set sampName, file(fastqs) from ch_outSkewer

   output:
     set sampName, file("*.zip"), file("*.html") into ch_outFastQC


   publishDir path: './trimmedFastq/fastqc', mode: 'copy'

    module      fastqcModule

   script:
   """
   fastqc -t ${task.cpus} ${fastqs[0]} ${fastqs[1]}
   """
}
