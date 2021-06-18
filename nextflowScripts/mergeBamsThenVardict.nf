#!/usr/bin/env nextflow

/*
 * 2021-06-15 This is a working script to merge bam files from two flowcells
 * Then rerun of the QC and Vardict processing on the merged BAM
 * requires some initial BASH scripting to rename the original BAM files into
 * Something akin to a fastq file going into the script.
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
Channel.fromFilePairs("flowcell1/*_{1,2}.bam")
  .set{ ch_bam }

process mergeBam {

   label 'bwa_genomics'

   input:
     set sampName, file(fq1), file(fq2) from ch_bam

   output:
     set sampName, file("${sampName}.merged.bam"), file("${sampName}.merged.bam.bai") into (ch_mappedBams, ch_mappedBams2, ch_mappedBams3, ch_mappedBams4, ch_mappedBams5, ch_mappedBams6)
     
   publishDir path: './output/out_bam', mode: 'copy'

    module      samtoolsModule

   script:
   """
   samtools merge -@ ${task.cpus} ${sampName}.merged.bam ${fq1} ${fq2} 

   samtools index "${sampName}.merged.bam" "${sampName}.merged.bam.bai"
   """
}

process bam_stats {
   label 'fastqc'
	
   input:
     set sampName, file(bam), file(bai) from ch_mappedBams6

   output:
      set sampName, file("${sampName}.samtools.stats"), file("${sampName}.idxstats") into ch_outSAMStats
   
   publishDir path: './qc_out/samtools', mode: 'copy'
   
    module		samtoolsModule
    
   script:
   """
   samtools stats -@ ${task.cpus} ${bam} > ${sampName}.samtools.stats
   samtools idxstats ${bam} > ${sampName}.idxstats
   """
}

process bam_qc1 {
	
   label 'genomics_qc'

   input:
     set sampName, file(bam), file(bai) from ch_mappedBams

   output:
     set sampName, file("${sampName}_alignment_metrics.txt") into ch_outQC

   publishDir path: './qc_out/picard', mode: 'copy'

    module      picardModule
    module      samtoolsModule
    module      RModule

   script:
   """
   picard CollectAlignmentSummaryMetrics   \
     I="${bam}"  \
     O="${sampName}_alignment_metrics.txt"  \
     R="${ref}" 
   """
}

process bam_qc2 {
	
   label 'genomics_qc'

   input:
     set sampName, file(bam), file(bai) from ch_mappedBams2

   output:
     set sampName, file("${sampName}.mm*") into ch_outQC2

   publishDir path: './qc_out/picard', mode: 'copy'

    module      picardModule
    module      samtoolsModule
    module      RModule

   script:
   """
   picard CollectMultipleMetrics \
     I="${bam}"  \
     O="${sampName}.mm"  \
     R="${ref}" 
   """
}

process bam_qc3 {
	
   label 'genomics_qc'

   input:
     set sampName, file(bam), file(bai) from ch_mappedBams3

   output:
     set sampName, file("${sampName}_pcr_metrics.txt") into ch_outQC3

   publishDir path: './qc_out/picard', mode: 'copy'

    module      picardModule
    module      samtoolsModule
    module      RModule

   script:
   """
   picard CollectTargetedPcrMetrics   \
     I="${bam}"  \
     O="${sampName}_pcr_metrics.txt"  \
     R="${ref}" \
     AMPLICON_INTERVALS="${picardAmplicon}" \
     TARGET_INTERVALS="${picardAmplicon}"
   """
}

process coverage_qc {
	
	label 'genomics_qc'
	
	input:
	  set sampName, file(bam), file(bai) from ch_mappedBams4
	  
	output:
	  set sampName, file("${sampName}.bedtools.coverage") into ch_coverageQC
	  
	publishDir path: './qc_out/bedtools', mode: 'copy'
	
    module      bedtoolsModule
    
    script:
    """
    bedtools coverage -a ${target} \
      -b ${bam} \
      -hist > ${sampName}.bedtools.coverage
    """
	
}

process vardict {
    
    label 'vardict_comp'
    
    input:
      set sampName, file(bam), file(bai) from ch_mappedBams5
	  
    output:
      set sampName, file("${sampName}.tsv") into ch_vcfMake
	  
    publishDir path: './chip/raw_variants/', mode: 'copy'

    script:
    """
    module purge
    module load samtools
    export PATH=/home/nwong/bin/VarDict-1.7.0/bin:$PATH
    VarDict -G ${ref} -f ${AF_THR} -N ${sampName} -b ${bam} -th ${task.cpus} \
      ${vardictAmp}  \
      >  "${sampName}.tsv"
    """
}

process makeVCF {
    
    label 'vardict_loy'

    input:
        set sampName, file(tsv) from ch_vcfMake.filter{ sampName, file -> !file.isEmpty() }
    output:
        set sampName, file("${sampName}.vardict.vcf") into ch_vardict
    
    publishDir path: './chip/vardict', mode: 'copy'

    module RModule

    script:
    """
    module purge
    module load R/3.6.0-mkl
    cat ${tsv} | /home/nwong/bin/VarDict-1.7.0/bin/teststrandbias.R | \
        /home/nwong/bin/VarDict-1.7.0/bin/var2vcf_valid.pl -N "${sampName}" \
        -f ${AF_THR} -E > "${sampName}.vardict.vcf"
    """
}

process vep {

    label 'vep'

    containerOptions "-B ${refFolder}/vep:/opt/vep/.vep"

    input:
        set sampName, file(vardict) from ch_vardict
    output:
        set sampName, file("${sampName}.vep*") into ch_VEPDone
    
    module singularityModule

    publishDir path: './chip/vep', mode: 'copy'

    script:
    """
    vep -i ${vardict} \
        -o ${sampName}.vep.vcf \
        --everything \
        --fork ${task.cpus} \
        --dir /opt/vep/.vep \
        --offline \
        --vcf
    """
}









