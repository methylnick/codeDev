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
Channel.fromPath("bamFiles/*.bam")
  .set{ ch_vcf }
  
process tumNorm {
    
    label 'vardict_genomics'
    
    input:
      set sampName, file(vcf) from ch_vcf
	  
    output:
      set sampName, file("${sampName}.vcf.gz"), file("${sampName}.vcf.gz.stats"), file("${sampName}.vcf.gz.tbi") into ch_gatkOut2
	  
    publishDir path: './tumourNormal/', mode: 'copy'

    script:
    """
    module purge
    module load samtools
    export PATH=/home/nwong/bin/gatk-4.2.0.0:$PATH
    gatk Mutect2  \
     -R ${ref} \
     -I ${bam} \
     -I  IR1_1_7_HFJ25DSX2_AGAGTCAA_L002_sorted.mdups.bam \
     --germline-resource ${refFolder}/af-only-gnomad.hg38.vcf.gz \
     -O ${sampName}.vcf.gz \
    """
}