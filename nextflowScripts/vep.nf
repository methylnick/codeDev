#!/usr/bin/env nextflow

/*
 * This is the nextflow script to run a singularity container to run VEP
 * Testing only, using VEP 102 and need to point to 102 annotations too
 */

// Declare Inputs
refFolder = file("/home/nwong/ls25_scratch/nick.wong/refFiles")

// Declare References 
refBase        = "${refFolder}/genome"
ref            = "${refBase}.fa"

// Tools
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
        --assembly 102_GRCh37 \
        --everything \
        --fork ${task.cpus} \
        --dir_cache /opt/vep/.vep \
        -offline
    """
}
