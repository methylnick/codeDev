#!/usr/bin/env nextflow

/*
 * 2021-07-14 this is a nf script for WGS alignment and processing to get to basic sequencing
 * metrics for gender typing
 */

// Declare Inputs
refFolder = file("/scratch/tn26/nick.wong/refFiles/hg38")

// Declare References 
ref            = "${refFolder}/bwa_0.7.17-gcc5/hg38_ref_noAlt.fa"

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
Channel.fromPath("out_bam/*.sorted.bam")
  .set{ ch_bamIn }

process picardMarkDup {
	
   label 'vardict_genomics'

   input:
     set file(bam) from ch_bamIn

   output:
     set sampName, file("${bam.getSimpleName()}_sorted.mdups.bam"), file("${bam.getSimpleName()}_sorted.mdups.metrics") into (ch_outMdups, ch_outMdups2)

   publishDir path: './bamFiles', mode: 'copy'

    module      picardModule

   script:
   """
   picard MarkDuplicates \
       TMP_DIR=${bam.getSimpleName()}_tmp \
       VALIDATION_STRINGENCY=LENIENT \
       INPUT=${bam} \
       OUTPUT=${bam.getSimpleName()}_sorted.mdups.bam \
       METRICS_FILE=${bam.getSimpleName()}_sorted.mdups.metrics
   """
}

process mutect2 {
    
    label 'vardict_genomics'
    
    input:
      set file(bam), file(bai) from ch_outMdups
	  
    output:
      set file("${bam.getSimpleName()}.vcf.gz") into ch_gatkOut
	  
    publishDir path: './tumourOnly/', mode: 'copy'

    script:
    """
    module purge
    module load samtools
    export PATH=/home/nwong/bin/gatk-4.2.0.0:$PATH
    gatk Mutect2  \
     -R ${ref} \
     -I ${bam} \
     -O ${sbam.getSimpleName()}.vcf.gz \
    """
}

process tumNorm {
    
    label 'vardict_genomics'
    
    input:
      set file(bam), file(bai) from ch_outMdups2
	  
    output:
      set file("${bam.getSimpleName()}.vcf.gz") into ch_gatkOut2
	  
    publishDir path: './tumourNormal/', mode: 'copy'

    script:
    """
    module purge
    module load samtools
    export PATH=/home/nwong/bin/gatk-4.2.0.0:$PATH
    gatk Mutect2  \
     -R ${ref} \
     -I ${bam} \
     -I  bamFiles/HFJ25DSX2_AGAGTCAA_sorted.mdups.bam \
     --germline-resource ${refFolder}/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
     -O ${bam.getSimpleName()}.g.vcf.gz \
    """
}
