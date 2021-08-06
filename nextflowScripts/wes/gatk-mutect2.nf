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
Channel.fromFilePairs("fastqs/*_{1,2}.fastq.gz")
  .set{ ch_fastaIn }

  Channel.fromFilePairs("fastqs/*_{1,2}.fastq.gz")
  .set{ ch_fastQC}

process fastqc {
   
   label 'fastqc'
   errorStrategy 'ignore'

   input:
     set sampName, file(fastqs) from ch_fastQC

   output:
     set sampName, file("*.zip"), file("*.html") into ch_outFastQC


   publishDir path: './qc_out/fastqc', mode: 'copy'

   script:
   """
   export PATH=/home/nwong/bin/FastQC:$PATH
   fastqc -t ${task.cpus} ${fastqs[0]} ${fastqs[1]}
   """
}

process skewer {
   label 'fastqc'   

   input:
	 set sampName, file(rawfqs) from ch_fastaIn
   
   output:
     set sampName, file("${sampName}-trimmed-pair1.fastq.gz"), file("${sampName}-trimmed-pair2.fastq.gz") into (ch_fastaToBwa, ch_fastaToFastqc)
     
   publishDir path: './fastq_trimmed', mode: 'copy'
     module skewerModule
   
   script:
   """
   skewer -t ${task.cpus} -q 20 ${rawfqs[0]} ${rawfqs[1]} -z -o "${sampName}"
   """
}

process fastqc_skew {

   label 'fastqc'
   errorStrategy 'ignore'

   input:
     set sampName, file(fq1), file(fq2) from ch_fastaToFastqc

   output:
     set sampName, file("*.zip"), file("*.html") into ch_outFastQC2

   publishDir path: './qc_out/fastqc_trimmed', mode: 'copy'

   script:
   """
   export PATH=/home/nwong/bin/FastQC:$PATH
   fastqc -t ${task.cpus} ${fq1} ${fq2}
   """
}

process align_bwa {

   label 'bwa_genomics'

   input:
     set sampName, file(fq1), file(fq2) from ch_fastaToBwa

   output:
     set sampName, file("${sampName}.sorted.bam"), file("${sampName}.sorted.bam.bai") into (ch_mappedBams, ch_mappedBams2, ch_mappedBams3, ch_mappedBams4)
     
   publishDir path: './chip/out_bam', mode: 'copy'

    module      bwaModule
    module      samtoolsModule

   script:
   """
   bwa mem -t ${task.cpus} -R "@RG\\tID:${sampName}\\tPU:${sampName}\\tSM:${sampName}\\tPL:ILLUMINA\\tLB:rhAmpSeq" \
       $ref ${fq1} ${fq2} | samtools view -u -h -q 1 - \
       | samtools sort -@ ${task.cpus} -o "${sampName}.sorted.bam"
   samtools index "${sampName}.sorted.bam" "${sampName}.sorted.bam.bai"
   """
}


process bam_stats {
   label 'fastqc'
	
   input:
     set sampName, file(bam), file(bai) from ch_mappedBams3

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


process picardMarkDup {
	
   label 'vardict_genomics'

   input:
     set sampName, file(bam), file(bai) from ch_mappedBams4

   output:
     set sampName, file("${sampName}_sorted.mdups.bam"), file("${sampName}_sorted.mdups.bam.bai"), file("${sampName}_sorted.mdups.metrics") into (ch_outMdups, ch_outMdups2)

   publishDir path: './bamFiles', mode: 'copy'

    module      picardModule

   script:
   """
   picard MarkDuplicates \
       TMP_DIR=${sampName}_tmp \
       VALIDATION_STRINGENCY=LENIENT \
       INPUT=${bam} \
       OUTPUT=${sampName}_sorted.mdups.bam \
       METRICS_FILE=${sampName}_sorted.mdups.metrics
   """
}

process mutect2 {
    
    label 'vardict_genomics'
    
    input:
      set sampName, file(bam), file(bai) , file(metrics) from ch_outMdups
	  
    output:
      set sampName, file("${sampName}.vcf.gz"), file("${sampName}.vcf.gz.tbi"), file("${sampName}.vcf.gz.stats") into ch_gatkOut
	  
    publishDir path: './tumourOnly/', mode: 'copy'

    script:
    """
    module purge
    module load samtools
    export PATH=/home/nwong/bin/gatk-4.2.0.0:$PATH
    gatk Mutect2  \
     -R ${ref} \
     -I ${bam} \
     -O ${sampName}.vcf.gz 
    """
}

process filterMutect {
    
    label 'vardict_genomics'
    
    input:
      set sampName, file(vcf), file(tbi) , file(stats) from ch_gatkOut
	  
    output:
      set sampName, file("${sampName}_filtered.vcf.gz"), file("${sampName}_filtered.vcf.gz.tbi"), file("${sampName}_filtered.vcf.gz.filteringStats.tsv") into ch_filterOutOut
	  
    publishDir path: './tumourOnly/', mode: 'copy'

    script:
    """
    module purge
    module load samtools
    export PATH=/home/nwong/bin/gatk-4.2.0.0:$PATH
    gatk FilterMutectCalls  \
     -R ${ref} \
     -V ${vcf} \
     -O ${sampName}_filtered.vcf.gz 
    """
}

process tumNorm {
    
    label 'vardict_genomics'
    
    input:
      set sampName, file(bam), file(bai), file(metrics) from ch_outMdups2
	  
    output:
      set sampName, file("${sampName}.vcf.gz") into ch_gatkOut2
	  
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
     -O ${sampName}.g.vcf.gz \
    """
}
