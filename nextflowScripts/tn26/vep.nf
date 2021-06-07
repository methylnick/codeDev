#!/usr/bin/env nextflow

/*
 * This is the nextflow script to run a singularity container to run VEP
 * for project tn26, using VEP 102 and need to point to 102 annotations
 */

// Declare Inputs
refFolder = file("/projects/tn26/referenceFiles/human/ucschg19")

// Declare References 
refBase        = "${refFolder}/genome"
ref            = "${refBase}.fa"

// Tools
singularityModule = 'singularity/3.7.1'

// Create channel stream
Channel.fromPath("*.vcf.gz")
  .set{ ch_vardict }
  
process vep {

    containerOptions "-B ${refFolder}/vep:/opt/vep/.vep"

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
