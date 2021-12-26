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

    label 'vep'

    containerOptions "-B ${refFolder}/vep:/opt/vep/.vep"

    input:
        file(vardict) from ch_vardict
    output:
        file("${vardict.getSimpleName()}.vep*") into ch_VEPDone

    module singularityModule

    publishDir path: './chip/vep', mode: 'copy'

    script:
    """
    vep -i ${vardict} \
        -o ${vardict.getSimpleName()}.vep.vcf \
        --everything \
        --fork ${task.cpus} \
        --dir /opt/vep/.vep \
        --offline \
        --vcf
    """
}

process maf {

    label 'maf'

    containerOptions "-B ${refFolder}:/refFiles"

    input:
        set sampName, file(vepInFile) from ch_VEPDone
    output:
        set sampName, file("${sampName}.maf"), file("${sampName}*.vcf")  into ch_MAFDone
    
    module singularityModule

    publishDir path: './chip/maf', mode: 'copy'

    script:
    """
    perl /opt/vep/src/ensembl-vep/vcf2maf-main/vcf2maf.pl \
    --input-vcf ${vepInFile} \
    --output-maf ${sampName}.maf \
    --cache-version=102 \
    --ref-fasta /refFiles/genome.fa \
    --vep-path /opt/vep/src/ensembl-vep \
    --vep-data /refFiles/vep
    """
}