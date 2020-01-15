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
picardModule   = 'picard/2.9.2'
bwaModule      = 'bwa/0.7.17-gcc5'
samtoolsModule = 'samtools/1.9-gcc5'
RModule        = 'R/3.6.0-mkl(default)'

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
     set sampName, file("${sampName}.sorted.bam"), file("${sampName}.sorted.bam.bai") into ch_mappedBams

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
   bwa mem -t ${task.cpus} -R "@RG\\tID:${sampName}\\tPU:${sampName}\\tSM:${sampName}\\tPL:ILLUMINA\\tLB:rhAmpSeq" \
       $ref ${fastqs[0]} ${fastqs[1]} | samtools view -u -h -q 1 - \
       | samtools sort -@ $bwaCores -o "${sampName}.sorted.bam"
   samtools index "${sampName}.sorted.bam" "${sampName}.sorted.bam.bai"
   """
}

process bam_qc {

   input:
     set sampName, file(bam), file(bai) from ch_mappedBams

   output:
     set sampName, file("${sampName}_alignment_metrics.txt"), file("${sampName}_multiple_metrics.txt") into ch_outQC

   publishDir path: './qc_out', mode: 'copy'

    executor    globalExecutor
    stageInMode globalStageInMode
    module      picardModule
    module      samtoolsModule
    cpus        globalCores
    memory      globalMemoryM
    time        globalTimeM
    queue       globalQueueL

   script:
   """
   picard CollectAlignmentSummaryMetrics   \
     I="${bam}"  \
     O="${sampName}_alignment_metrics.txt"  \
     R="${ref}" 

   picard CollectMultipleMetrics \
     I="${bam}"  \
     O="${sampName}_multiple_metrics.txt"  \
     R="${ref}" 
   """
}
