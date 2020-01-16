#!/usr/bin/env nextflow

/*
 * Another attempt at a nextflow script 2020-01-17
 * This is a test and dev script so it is going straight to the short queue
 */

// Declare Inputs
refFolder = file("refFiles")
inputDirectory = file('fastq')
AF_THR           = 0.0025

// Declare References 
refBase = "${refFolder}/genome"
ref     = "${refBase}.fa"
target  = "${refFolder}/rhAMPTargetROI_hg19.sorted.bed"
vardictAmp = "${refFolder}/vardict_amplicon_rhAMPSeq.sorted.bed"

// Tools
picardModule   = 'picard/2.9.2'
bwaModule      = 'bwa/0.7.17-gcc5'
samtoolsModule = 'samtools/1.9-gcc5'
RModule        = 'R/3.6.0-mkl(default)'
bedtoolsModule = 'bedtools/2.27.1-gcc5'
bcftoolsModule = 'bcftools/1.8'

// Global Resource Configuration Options
globalExecutor    = 'slurm'
globalStageInMode = 'symlink'
globalCores       = 1
bwaCores	      = 4
vepCores          = 4
globalMemoryS     = '6 GB'
globalMemoryM     = '8 GB'
globalMemoryL     = '64 GB'
globalTimeS       = '30m'
globalTimeM       = '3h'
globalTimeL       = '24h'
globalQueueS      = 'short'
globalQueueL      = 'comp'

// Create channel stream
Channel.fromFilePairs("rawFiles/catFastq/*_R{1,2}_001.fastq.gz")
  .set{ ch_fastaIn }

process align_bwa {

   input:
     set sampName, file(fastqs) from ch_fastaIn

   output:
     set sampName, file("${sampName}.sorted.bam"), file("${sampName}.sorted.bam.bai") into ch_mappedBams

   publishDir path: './out_bam', mode: 'copy'

    executor    globalExecutor
    stageInMode globalStageInMode
    module      bwaModule
    module      samtoolsModule
    cpus        bwaCores
    memory      globalMemoryM
    time        globalTimeS
    queue       globalQueueS

   script:
   """
   bwa mem -t ${task.cpus} -R "@RG\\tID:${sampName}\\tPU:${sampName}\\tSM:${sampName}\\tPL:ILLUMINA\\tLB:rhAmpSeq" \
       $ref ${fastqs[0]} ${fastqs[1]} | samtools view -u -h -q 1 - \
       | samtools sort -@ $bwaCores -o "${sampName}.sorted.bam"
   samtools index "${sampName}.sorted.bam" "${sampName}.sorted.bam.bai"
   """
}

process bam_qc {

   input:
     set sampName, file(bam), file(bai) from ch_mappedBams

   output:
     set sampName, file("${sampName}_alignment_metrics.txt"), file("${sampName}_multiple_metrics.txt") into ch_outQC

   publishDir path: './qc_out/picard', mode: 'copy'

    executor    globalExecutor
    stageInMode globalStageInMode
    module      picardModule
    module      samtoolsModule
    module      RModule
    cpus        globalCores
    memory      globalMemoryM
    time        globalTimeS
    queue       globalQueueS

   script:
   """
   picard CollectAlignmentSummaryMetrics   \
     I="${bam}"  \
     O="${sampName}_alignment_metrics.txt"  \
     R="${ref}" 

   picard CollectMultipleMetrics \
     I="${bam}"  \
     O="${sampName}_multiple_metrics.txt"  \
     R="${ref}" 
     
   picard picard CollectTargetedPcrMetrics   \
     I="${bam}"  \
     O="${sampName}_pcr_metrics.txt"  \
     R="${ref}" \
     AMPLICON_INTERVALS="${refFolder}/rhAMPTargetROI_hg19.interval_list" \
     TARGET_INTERVALS="${refFolder}/CHIPMutations_hg19_050219.extend.interval_list"
   """
}

process coverage_qc {
	input:
	  set sampName, file(bam), file(bai) from ch_mappedBams
	  
	output:
	  set sampName, file("${sampName}.bedtools.coverage") into ch_coverageQC
	  
	publishDir path: '.qc_out/bedtools', mode: 'copy'
	
	executor    globalExecutor
    stageInMode globalStageInMode
    module      bedtoolsModule
    cpus        globalCores
    memory      globalMemoryM
    time        globalTimeS
    queue       globalQueueS
    
    script:
    """
    bedtools coverage -a ${target}. \
      -b ${bam} \
      -hist > ${sampName}.bedtools.coverage
    """
	
}

process vardict {
	input:
	  set sampName, file(bam), file(bai) from ch_mappedBams
	  
	output:
	  set sampName, file("sampName}.vcf") into ch_vardict
	  
	publishDir path: '.vardict/', mode: 'copy'
	
	executor    globalExecutor
    stageInMode globalStageInMode
    module      bwaModule
    cpus        globalCores
    memory      globalMemoryM
    time        globalTimeS
    queue       globalQueueS
    
    script:
    """
    export PATH=/home/nwong/bin/VarDict-master:$PATH
    vardict -G ${ref} -f ${AF_THR} -N ${sampName} -b ${bam} \
      -z 0 -c 1 -S 2 -E 3 -g 4 -y \
      ${vardictAmp} | teststrandbias.R  \
      | var2vcf_valid.pl \
      -N ${sampName} -E -f ${AF_THR} > "${sampName}.vardict.vcf"

    """
}

process sortVCFS {

    input:
        set sampName, file(vcf) from ch_vardict
    output:
        set sampName, file("${sampName}.sorted.vcf.gz") into ch_sortedVCF

    publishDir path: './variants_raw_out', mode: 'copy'                                    
    
    module      bcftoolsModule
    executor    globalExecutor                                                    
    stageInMode globalStageInMode                                                 
    module      bwaModule
    memory      globalMemoryM 
    time        globalTimeS
    queue       globalQueueS

    script:
    """
    bcftools sort -o "${sampName}.sorted.vcf.gz" -O z ${vcf}
    """
}

process indexVCFS {
    input:
        set baseName, file(vcf) from ch_sortedVCF
    output:
        set baseName, file(vcf), file("${sampName}.sorted.vcf.gz.tbi") into ch_indexedVCF

    publishDir path: './variants_raw_out', mode: 'copy'                                    
    
    module      bcftoolsModule
    executor    globalExecutor                                                    
    stageInMode globalStageInMode                                                 
    module      bwaModule
    memory      globalMemoryM 
    time        globalTimeS
    queue       globalQueueS

    script:
    """
    bcftools index -f --tbi ${vcf} -o ${sampName}.sorted.vcf.gz.tbi
    """
}

