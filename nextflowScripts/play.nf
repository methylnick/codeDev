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

// Tools
picardJar      = '~/picard.jar'
bwaModule      = 'bwa/0.7.17-gcc5'
samtoolsModule = 'samtools/1.9'

// Global Resource Configuration Options
globalExecutor    = 'slurm'
globalStageInMode = 'symlink'
globalCores       = 1
bwaCores	  = 4
vepCores          = 4
globalMemoryS     = '6 GB'
globalMemoryM     = '8 GB'
globalMemoryL     = '64 GB'
globalTimeS       = '8m'
globalTimeM       = '3h'
globalTimeL       = '24h'
globalQueueS      = 'short'
globalQueueL      = 'comp'

// Create channel stream
Channel.fromFilePairs("rawFiles/catFastq/*_R{1,2}_001.fastq.gz")
  .set{ ch_fastaIn }

process align_bwa {

   input:
     set sampName, file(fastqs) from ch_fastaIn

   output:
     set sampName, file("${sampName}.sam") into ch_mappedBams

   publishDir path: './out_bam', mode: 'copy'

    executor    globalExecutor
    stageInMode globalStageInMode
    module      bwaModule
    module      samtoolsModule
    cpus        bwaCores
    memory      globalMemoryM
    time        globalTimeM
    queue       globalQueueL



   script:
   """
   bwa mem -t ${task.cpus} -R "@RG\tID:${sampName}\tPU:${sampName}\tSM:${sampName}\tPL:ILLUMINA\tLB:rhAmpSeq" \
       $ref ${fastqs[0]} ${fastqs[1]} > ${sampName}.sam
   """
}
