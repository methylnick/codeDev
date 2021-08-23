#!/usr/bin/env nextflow

/*
 * 2021-07-26 this is a nf script for extracting fastq files from a BAM from 
 * TCGA, then recreating the bam after alignment with appropriate reference
 * genome
 */

// Declare Inputs

// Declare reference files 

// Tools
samtoolsModule = 'samtools/1.9-gcc5'
bedtoolsModule = 'bedtools/2.27.1-gcc5'

// Create channel stream
Channel.fromPath("*.bam")
  .set{ ch_bam }


process bamtofastq {
   
   label 'fastqc'

   input:
    file(bam) from ch_bam

   output:
    file("*.zip"), file("*.html") into ch_outFastQC

   publishDir path: './fastqs', mode: 'copy'

   module samtoolsModule
   module bedtoolsModule

   script:
   """
   samtools view -u -f 1 -F 12 ${bam} | \
     samtools sort -@ ${task.cpus} -T ${bam.getSimpleName()} -n ${bam} | \
     samtools fastq -@ ${task.cpus} -1 ${bam.getSimpleName()}_mapped.R1.fastq.gz \
     -2 ${bam.getSimpleName()}_mapped.R2.fastq.gz \
     -0 /dev/null -s /dev/null -n
   
   samtools view -u -f 4 -F 264 ${bam} | \
     samtools sort -@ ${task.cpus} -T ${bam.getSimpleName()} -n ${bam} | \
     samtools fastq -@ ${task.cpus} -1 ${bam.getSimpleName()}_u_m.R1.fastq.gz \
     -2 ${bam.getSimpleName()}_u_m.R2.fastq.gz \
     -0 /dev/null -s /dev/null -n

   samtools view -u -f 8 -F 260 ${bam} | \
     samtools sort -@ ${task.cpus} -T ${bam.getSimpleName()} -n ${bam} | \
     samtools fastq -@ ${task.cpus} -1 ${bam.getSimpleName()}_m_u.R1.fastq.gz \
     -2 ${bam.getSimpleName()}_m_u.R2.fastq.gz \
     -0 /dev/null -s /dev/null -n

   samtools view -u -f 12 -F 256 ${bam} | \
     samtools sort -@ ${task.cpus} -T ${bam.getSimpleName()} -n ${bam} | \
     samtools fastq -@ ${task.cpus} -1 ${bam.getSimpleName()}_u_u.R1.fastq.gz \
     -2 ${bam.getSimpleName()}_u_u.R2.fastq.gz \
     -0 /dev/null -s /dev/null -n
   """
}