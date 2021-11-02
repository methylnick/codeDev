#!/usr/bin/env nextflow

/*
 * This is the nextflow script to run GATK picard HS metrics of
 * sorted mduped bam files. 
 */

// Declare Inputs
refFolder = file("/projects/rc78/refFiles")

// Declare References 
refBase        = "${refFolder}/hg38/bwa_0.7.17-gcc5/hg38_ref_noAlt"
ref            = "${refBase}.fa"
target         = "${refFolder}/hg38/agilentClinicalResearchExomeV2/S30409818_Covered.interval"
baits          = "${refFolder}/hg38/agilentClinicalResearchExomeV2/S30409818_Padded.interval"

// Tools
singularityModule = 'singularity/3.7.1'

// Create channel stream
Channel.fromPath("/home/nwong/sg58/2021-07-21-WES/bamFiles/*.bam")
  .set{ ch_bam }
  
process tumNorm {
    
    label 'genomics_qc'
    
    input:
      file(bam) from ch_bam
	  
    output:
      file("${bam.getSimpleName()}.HS_metrics.txt") into ch_gatkOut2
	  
    publishDir path: './qc_out', mode: 'copy'

    script:
    """
    module purge
    module load samtools
    export PATH=/home/nwong/bin/gatk-4.2.0.0:$PATH
    gatk CollectHsMetrics  \
     -I ${bam} \
     -BI ${baits} \
     -TI ${target} \
     -O ${bam.getSimpleName()}.HS_metrics.txt
    """
}
