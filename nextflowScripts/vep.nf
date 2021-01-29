#!/usr/bin/env nextflow

/*
 * 2021-01-12 Working up and testing the nf script for production
 * This is a test and dev script so it is going straight to the short queue
 * Include AMELX loss of Y processing and vardict calling 
 * BAM filter employed to extract alignments most likely to be associated with
 * assay proper. Vardict calling the insertion deployed to call Loss of Y
 */

// Declare Inputs
refFolder = file("/home/nwong/ls25_scratch/nick.wong/refFiles")
inputDirectory = file('fastq')
AF_THR          = 0.0025

// Declare References 
refBase        = "${refFolder}/genome"
ref            = "${refBase}.fa"
loyRefBase     = "${refFolder}/loy/hg19_chrX.fa.gz"
target         = "${refFolder}/chip_v2_amplicon_bed_target_20200324.bed"
picardInsert   = "${refFolder}/chip_v2_amplicon_bed_target_20200324.bed.interval"
picardAmplicon = "${refFolder}/chip_v2_amplicon_bed_amplicon_20200324.bed.interval"
vardictAmp     = "${refFolder}/chip_v2_amplicon_bed_vardict_20200324.bed"
loyVardictAmp  = "${refFolder}/chip_v2_loy_amplicon_bed_vardict_20210113.bed"

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
Channel.fromPath("*.vcf")
  .set{ ch_vardict }
  
process vep {

    containerOptions "-B ${refFolder}/vep:/opt/vep/.vep"

    input:
        file(vardict) from ch_vardict
    output:
        file("${vardict.getSimpleName()}.vep*") into ch_VEPDone

    publishDir path: './chip/vep', mode: 'copy'

    script:
    """
    vep -i ${vardict} \
        -o ${vardict.getSimpleName()}.vep.vcf \
        --everything \
        --fork ${task.cpus} \
        --dir_cache /opt/vep/.vep \
        -offline
    """
}
