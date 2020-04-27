#!/usr/bin/env nextflow

/*
 * Another attempt at a nextflow script 2020-01-18
 * This is a test and dev script so it is going straight to the short queue
 */

// Declare Inputs
refFolder = file("/home/nwong/ls25_scratch/nick.wong/refFiles")
inputDirectory = file('fastq')
AF_THR          = 0.0025

// Declare References 
refBase        = "${refFolder}/genome"
ref            = "${refBase}.fa"
target         = "${refFolder}/chip_v2_amplicon_bed_target_20200324.bed"
picardInsert   = "${refFolder}/chip_v2_amplicon_bed_target_20200324.bed.interval"
picardAmplicon = "${refFolder}/chip_v2_amplicon_bed_amplicon_20200324.bed.interval"
vardictAmp     = "${refFolder}/chip_v2_amplicon_bed_vardict_20200324.bed"

// Tools
picardModule   = 'picard/2.9.2'
bwaModule      = 'bwa/0.7.17-gcc5'
samtoolsModule = 'samtools/1.9-gcc5'
bedtoolsModule = 'bedtools/2.27.1-gcc5'
bcftoolsModule = 'bcftools/1.8'
RModule        = 'R/3.6.0-mkl'
fastqcModule         = 'fastqc/0.11.7'

// Create channel stream
Channel.fromFilePairs("fastq/*_R{1,2}_001.fastq.gz")
  .set{ ch_fastaIn }
  
Channel.fromFilePairs("fastq/*_R{1,2}_001.fastq.gz")
  .set{ ch_fastQC }

process fastqc {
   
   label 'fastqc'

   input:
     set sampName, file(fastqs) from ch_fastQC

   output:
     set sampName, file("*.zip"), file("*.html") into ch_outFastQC


   publishDir path: './qc_out/fastqc', mode: 'copy'

    module      fastqcModule

   script:
   """
   fastqc -t ${task.cpus} ${fastqs[0]} ${fastqs[1]}
   """
}

process align_bwa {

   label 'bwa'

   input:
     set sampName, file(fastqs) from ch_fastaIn

   output:
     set sampName, file("${sampName}.sorted.bam"), file("${sampName}.sorted.bam.bai") into ch_mappedBams
     set sampName, file("${sampName}.sorted.bam"), file("${sampName}.sorted.bam.bai") into ch_mappedBams2
     set sampName, file("${sampName}.sorted.bam"), file("${sampName}.sorted.bam.bai") into ch_mappedBams3
     set sampName, file("${sampName}.sorted.bam"), file("${sampName}.sorted.bam.bai") into ch_mappedBams4
     set sampName, file("${sampName}.sorted.bam"), file("${sampName}.sorted.bam.bai") into ch_mappedBams5
     set sampName, file("${sampName}.sorted.bam"), file("${sampName}.sorted.bam.bai") into ch_mappedBams6
     
   publishDir path: './out_bam', mode: 'copy'

    module      bwaModule
    module      samtoolsModule

   script:
   """
   bwa mem -t ${task.cpus} -R "@RG\\tID:${sampName}\\tPU:${sampName}\\tSM:${sampName}\\tPL:ILLUMINA\\tLB:rhAmpSeq" \
       $ref ${fastqs[0]} ${fastqs[1]} | samtools view -u -h -q 1 - \
       | samtools sort -@ $task.cpus -o "${sampName}.sorted.bam"
   samtools index "${sampName}.sorted.bam" "${sampName}.sorted.bam.bai"
   """
}

process bam_stats {
   label 'fastqc'
	
   input:
     set sampName, file(bam), file(bai) from ch_mappedBams6

   output:
      set sampName, file("${sampName}.samtools.stats"), file("${sampName}.idxstats") into ch_outSAMStats
   
   publishDir path: '.qc_out/samtools', mode: 'copy'
   
    module		samtoolsModule
    
   script:
   """
   samtools stats -@ ${task.cpus} ${bam} > ${sampName}.samtools.stats
   samtools idxstats ${bam} > ${sampName}.idxstats
   """
}

process bam_qc1 {
	
   label 'medium_1h'

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
	
   label 'medium_1h'

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
	
   label 'medium_1h'

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
	
	label 'medium_1h'
	
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
    
    label 'vardict'
    
    input:
      set sampName, file(bam), file(bai) from ch_mappedBams5
	  
    output:
      set sampName, file("${sampName}.tsv") into ch_vcfMake
	  
    publishDir path: './raw_variants/', mode: 'copy'

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
    
    label 'small_1'

    input:
        set sampName, file(tsv) from ch_vcfMake
    output:
        set sampName, file("${sampName}.vardict.vcf") into ch_vardict
    
    publishDir path: './vardict', mode: 'copy'

    script:
    """
    module purge
    module load R
    export PATH=/home/nwong/bin/VarDict-1.7.0/bin:$PATH
    cat ${tsv} | teststrandbias.R | \
        var2vcf_valid.pl -N "${sampName}" \
        -f ${AF_THR} -E > "${sampName}.vardict.vcf"
    """
}

process sortVCFS {

    label 'medium_1h'

    input:
        set sampName, file(vcf) from ch_vardict
    output:
        set sampName, file("${sampName}.sorted.vcf.gz") into ch_sortedVCF

    publishDir path: './variants_raw_out', mode: 'copy'                                    
    
    module      bcftoolsModule                                               
    module      bwaModule

    script:
    """
    bcftools sort -o "${sampName}.sorted.vcf.gz" -O z ${vcf}
    """
}

process indexVCFS {
    
    label 'small_1'
    
    input:
        set sampName, file(vcf) from ch_sortedVCF
    output:
        set sampName, file(vcf), file("${sampName}.sorted.vcf.gz.tbi") into ch_indexedVCF

    publishDir path: './variants_raw_out', mode: 'copy'                                    
    
    module      bcftoolsModule                                                
    module      bwaModule

    script:
    """
    bcftools index -f --tbi ${vcf} -o "${sampName}.sorted.vcf.gz.tbi"
    """
}


