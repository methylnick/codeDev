#!/usr/bin/env nextflow

/*
 * 2021-07-14 this is a nf script for WGS alignment and processing to get to basic sequencing
 * metrics for gender typing
 */

// Declare Inputs
refFolder = file("/scratch/yu81/nick.wong/refFiles/ncbi-genomes-2022-04-27")
vepFolder = file("/scratch/yu81/nick.wong/refFiles/VEP")

// Declare References 
ref            = "${refFolder}/GCA_000002655.1_ASM265v1_genomic.fna"


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
     set sampName, file("${sampName}.sorted.bam"), file("${sampName}.sorted.bam.bai") into (ch_mappedBams, ch_mappedBams2, ch_mappedBams3, ch_mappedBams4, ch_mappedBams5)
     
   publishDir path: './out_bam', mode: 'copy'

    module      bwaModule
    module      samtoolsModule

   script:
   """
   bwa mem -t ${task.cpus} -R "@RG\\tID:${sampName}\\tPU:${sampName}\\tSM:${sampName}\\tPL:ILLUMINA\\tLB:rhAmpSeq" \
       ${ref} ${fq1} ${fq2} | samtools view -u -h -q 1 - \
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

process gatk {
    
    label 'vardict_genomics'
    
    input:
      set sampName, file(bam), file(bai) from ch_mappedBams4
	  
    output:
      set sampName, file("${sampName}.g.vcf.gz"), file("${sampName}.g.vcf.gz.tbi") into ch_gatkOut
	  
    publishDir path: './gvcfs/', mode: 'copy'

    script:
    """
    module purge
    module load samtools
    export PATH=/home/nwong/bin/gatk-4.2.0.0:$PATH
    gatk HaplotypeCaller  \
     -R ${ref} \
     -I ${bam} \
     -O ${sampName}.g.vcf.gz \
     -ERC GVCF \
     -ploidy 1 \
     -stand-call-conf 20.0 \
     --dont-use-soft-clipped-bases true
    """
}

process gatk_gt {
    
    label 'vardict_genomics'
    
    input:
      set sampName, file(vcf), file(tbi) from ch_gatkOut
	  
    output:
      set sampName, file("${sampName}.vcf.gz"), file("${sampName}.vcf.gz.tbi") into ch_gatkGTOut
	  
    publishDir path: './vcfs/', mode: 'copy'

    script:
    """
    module purge
    module load samtools
    export PATH=/home/nwong/bin/gatk-4.2.0.0:$PATH
    gatk GenotypeGVCFs  \
     -R ${ref} \
     -V ${vcf} \
     -O ${sampName}.vcf.gz 
    """
}

process vep {
     
    label 'vep'

    containerOptions "-B ${vepFolder}/vep:/opt/vep/.vep"

    input:
        file(vardict) from ch_gatkGTOut
    output:
        file("${sampName}.vep*") into ch_VEPDone

    publishDir path: './vep', mode: 'copy'

    module singularityModule

    script:
    """
    vep -i ${vardict} \
        -o ${sampName}.vep.vcf \
        --everything \
        --fork ${task.cpus} \
        --dir /opt/vep/.vep \
        --offline \
        --vcf \
        --species aspergillus_fumigatus \
        --cache_version 53
    """
}