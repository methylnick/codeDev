#!/usr/bin/env nextflow

/*
 * This is the nextflow script to run GATK Mutect in Tumour Normal Mode
 * on filtered VCF files. 
 */

// Declare Inputs
refFolder = file("/projects/rc78/refFiles")

// Declare References 
refBase        = "${refFolder}/hg38/bwa_0.7.17-gcc5/hg38_ref_noAlt"
ref            = "${refBase}.fa"

// Tools
singularityModule = 'singularity/3.7.1'

// Create channel stream
Channel.fromPath("/fs03/sg58/angela.pizzola/2021-07-21-WES/bamFiles/*.bam")
  .set{ ch_bam }
  
process tumNorm {
    
    label 'vep'
    
    input:
      file(bam) from ch_bam
	  
    output:
      tuple file("${bam.getSimpleName()}.vcf.gz"), file("${bam.getSimpleName()}.vcf.gz.stats"), file("${bam.getSimpleName()}.vcf.gz.tbi") into ch_gatkOut2
	  
    publishDir path: './tumourNormal/', mode: 'copy'

    script:
    """
    module purge
    module load samtools
    export PATH=/home/nwong/bin/gatk-4.2.0.0:$PATH
    gatk Mutect2  \
     -R ${ref} \
     -I ${bam} \
     -I  /fs03/sg58/angela.pizzola/2021-07-21-WES/IR1_1_7_HFJ25DSX2_AGAGTCAA_L002_sorted.mdups.bam \
     --germline-resource ${refFolder}/hg38/af-only-gnomad.hg38.vcf.gz \
     -O ${bam.getSimpleName()}.vcf.gz \
    """
}
