#!/usr/bin/env nextflow

/*
 * This is the nextflow script to run a singularity container immcantation
 * on sequencing using presto NEB/AbSeq pipeline. Source of the material is
 * here: https://immcantation.readthedocs.io/en/stable/docker/pipelines.html#abseqpipeline
 */

// Declare Inputs
params.index = '/scratch/rc78/angela.pizzola/2021-07-21-WES/bamFiles/bamsToSplit/manifest.csv'

// Tools
singularityModule = 'singularity/3.7.1'

// Create channel stream
Channel
    .fromPath(params.index)
    .splitCsv(header:true)
    .map{ row -> tuple(row.chr, file(row.norm), file(row.samp1), file(row.samp2), file(row.samp3), file(row.samp4), file(row.samp5)) }
    .set { samples_ch }
  
process presto {

    label 'presto'

    containerOptions "-B /scratch/mx82:/data"

    input:
        set chr, file(sampN), file(samp1), file(samp2), file(samp3), file(samp4), file(samp5) from samples_ch
    output:
        file("${chr}.multisnv.vcf") into ch_Done

    module singularityModule

    publishDir path: './multisnv', mode: 'copy'

    script:
    """
    /multisnv/multiSNV  \
       -N 6  \
       --fasta /opt/hg38_ref_noAlt.fa  \
       --bam  /opt/${sampN} \
              /opt/${samp1} \
              /opt/${samp2} \
              /opt/${samp3} \
              /opt/${samp4} \
              /opt/${samp5} \
  -f /opt/${chr}.multisnv.vcf
    """
}
