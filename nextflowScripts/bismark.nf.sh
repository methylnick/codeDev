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
refBase        = "${refFolder}/referenceFiles"
ref            = "${refBase}.fa"

/* Target regions commented out for future inclusion
 *
 * target         = "${refFolder}/chip_v2_amplicon_bed_target_20200324.bed"
 * picardInsert   = "${refFolder}/chip_v2_amplicon_bed_target_20200324.bed.interval"
 * picardAmplicon = "${refFolder}/chip_v2_amplicon_bed_amplicon_20200324.bed.interval"
 * vardictAmp     = "${refFolder}/chip_v2_amplicon_bed_vardict_20200324.bed"
 */


// Tools
picardModule   = 'picard/2.9.2'
bwaModule      = 'bwa/0.7.17-gcc5'
bowtieModule   = 'bowtie2/2.3.5'
samtoolsModule = 'samtools/1.9-gcc5'
bedtoolsModule = 'bedtools/2.27.1-gcc5'
bcftoolsModule = 'bcftools/1.8'
RModule        = 'R/3.6.0-mkl'
fastqcModule         = 'fastqc/0.11.7'
skewerModule   = 'skewer/20170212'

// Create channel stream
Channel.fromFilePairs("fastq/*_R{1,2}_001.fastq.gz")
  .set{ ch_fastaIn }
  
Channel.fromFilePairs("fastq/*_R{1,2}_001.fastq.gz")
  .set{ ch_fastQC }

process fastqc {
   
   label 'fastqc'

   input:
     set sampName, file(fastqs) from ch_fastQC

   output:
     set sampName, file("*.zip"), file("*.html") into ch_outFastQC


   publishDir path: './qc_out/fastqc', mode: 'copy'

    module      fastqcModule

   script:
   """
   fastqc -t ${task.cpus} ${fastqs[0]} ${fastqs[1]}
   """
}

process skewer {
   label 'bwa'   

   input:
	 set sampName, file(rawfqs) from ch_fastaIn
   
   output:
     set sampName, file("*-trimmed-pair1.fastq.gz"), file("*-trimmed-pair2.fastq.gz") into (ch_fastaToBwa, ch_fastaToFastqc)
     set sampName, file("*-trimmed-pair1.fastq.gz"), file("*-trimmed-pair2.fastq.gz"), file("*-trimmed.log") into ch_skewEnd
     
   publishDir path: './fastq/trimmed', mode: 'copy'
     module skewerModule
   
   script:
   """
   skewer -t ${task.cpus} -q 20 ${rawfqs[0]} ${rawfqs[1]} -z -o "${sampName}"
   """
}

process fastqc_skew {

   label 'fastqc'

   input:
     set sampName, file(fastqs) from ch_fastaToFastqc

   output:
     set sampName, file("*.zip"), file("*.html") into ch_outFastQC2


   publishDir path: './qc_out/fastqc_trimmed', mode: 'copy'

    module      fastqcModule

   script:
   """
   fastqc -t ${task.cpus} ${fastqs[0]} ${fastqs[1]}
   """
}

process align_bismark {

   label 'bwa_genomics_12'

   input:
     set sampName, file(fastqs) from ch_fastaToBwa

   output:
     set sampName, file("${sampName}_bismark_bt2_pe.bam"), file("${sampName}._bismark_bt2_PE_report.html"), file("${sampName}._bismark_bt2_PE_report.txt") into (ch_mappedBams, ch_mappedBams2, ch_mappedBams3, ch_mappedBams4, ch_mappedBams5, ch_mappedBams6)
     
   publishDir path: './out_bismark', mode: 'copy'

    module      bowtieModule
    module      samtoolsModule

   script:
   """
   /home/nwong/bin/Bismark_v0.19.0/bismark ${refBase}  \
   -1 ${fastqs[0]} \
   -2 ${fastqs[1]} \
   -q \
   -X 1000 \
   --multicore 6 \
   --un \
   -o out_bismark \
   --non_directional \
   --score_min L,0,-0.6 
   """
}

