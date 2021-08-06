#!/usr/bin/env nextflow

/*
 * This is the nextflow script to run a singularity container to run VEP
 * Testing only, using VEP 102 and need to point to 102 annotations too
 */

// Declare Inputs
refFolder = file("/projects/rc78/refFiles")

// Declare References 
refBase        = "${refFolder}/hg38/bwa_0.7.17-gcc5/hg38_ref_noAlt"
ref            = "${refBase}.fa"

// Tools
singularityModule = 'singularity/3.7.1'

// Create channel stream
Channel.fromPath("*.vcf.gz")
  .set{ ch_vardict }
  
process vep {

    containerOptions "-B ${refFolder}/vep_homo_sapiens:/opt/vep/.vep"

    input:
        file(vardict) from ch_vardict
    output:
        file("${vardict.getSimpleName()}.vep*") into ch_VEPDone

    publishDir path: './vep', mode: 'copy'

    script:
    """
    vep -i ${vardict} \
        -o ${vardict.getSimpleName()}.vep.vcf \
        --everything \
        --fork ${task.cpus} \
        --dir /opt/vep/.vep \
        --offline
    """
}
