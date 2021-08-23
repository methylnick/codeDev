#!/usr/bin/env nextflow

/*
 * 2021-08-24 This script is to extract read depth coverage of bam files given
 * a bed file of locations of interest. This was made for samples with no variants
 * to see if there is sufficient coverage and likely to be reference base. 
 */

// Declare Inputs
refFolder = file("/home/nwong/ls25_scratch/nick.wong/refFiles")

// Declare References 
refBase        = "${refFolder}/genome"
ref            = "${refBase}.fa"
target         = "${refFolder}/2021-08-24-BedFeatures.bed"

// Tools
picardModule   = 'picard/2.9.2'
bwaModule      = 'bwa/0.7.17-gcc5'
samtoolsModule = 'samtools/1.9-gcc5'
bedtoolsModule = 'bedtools/2.27.1-gcc5'
bcftoolsModule = 'bcftools/1.8'
RModule        = 'R/3.6.0-mkl'
fastqcModule         = 'fastqc/0.11.7'
skewerModule   = 'skewer/20170212'
singularityModule = 'singularity/3.7.1'

// Create channel stream
Channel.fromFilePath("chip/out_bam/*.bam")
  .set{ ch_bamIn }


process coverage_qc {
	
	label 'genomics_qc'
	
	input:
	  file(bam) from ch_bamIn
	  
	output:
	  file("${bam.getSimpleName()}.bedtools.coverage") into ch_coverageQC
	  
	publishDir path: './qc_out/bedtools', mode: 'copy'
	
    module      bedtoolsModule
    
    script:
    """
    bedtools coverage -a ${target} \
      -b ${bam} \
      -hist > ${bam.getSimpleName()}.bedtools.coverage
    """
	
}

