#!/usr/bin/env nextflow

/*
 * RNASeq NextFlow Pipeline
 * Human GRCh38 Ensembl Reference
 */

// Declare Inputs
refFolder = file("/home/nwong/df22_scratch/references/iGenomes")
inputDirectory = file('fastq')

// Declare References 
starRef        = "${refFolder}/Homo_sapiens/Ensembl/GRCh38/Sequence/STARIndex"
ref            = "${refFolder}/Homo_sapiens/Ensembl/GRCh38/Sequence/WholeGenomeFasta/genome.fa"
genes.         = "${refFolder}/Homo_sapiens/Ensembl/GRCh38/Annotation/Genes/genes.gtf"

// Tools
picardModule   = 'picard/2.9.2'
bwaModule      = 'bwa/0.7.17-gcc5'
samtoolsModule = 'samtools/1.9-gcc5'
bedtoolsModule = 'bedtools/2.27.1-gcc5'
bcftoolsModule = 'bcftools/1.8'
RModule        = 'R/3.6.0-mkl'
fastqcModule   = 'fastqc/0.11.7'
starModule     = 'star/2.5.2b'
subreadModule  = 'subread/1.5.1'


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

process align_STAR {

   label 'bwa'

   input:
     set sampName, file(fastqs) from ch_fastaIn

   output:
     set sampName, file("${sampName}_Aligned.out.bam"), file("${sampName}_Log.final.out"), file(${sampName}_Log.out"), file("${sampName}_Log.progress.out"), file("${sampName}_SJ.out.tab")  into ch_mappedBams
     

    module      starModule
    module      samtoolsModule

   script:
   """
   STAR --runThreadN ${task.cpus} \
       --runMode alignReads \
       --genomeDir ${starRef} \
       --outSAMattrRGline ID:${sampName} PU:${sampName} SM:${sampName} PL:ILLUMINA LB:RNASeq \
       --readFilesIn ${fastqs[0]} ${fastqs[1]} \
       --outFileNamePrefix ${sampName}_ \
       --outSAMtype BAM Unsorted --outSAMunmapped Within --readFilesCommand zcat
   """
}

process sort_BAM {
   label 'fastqc'
	
   input:
     set sampName, file(bam), file(Log[1-4]) from ch_mappedBams

   output:
      set sampName, file("${sampName}_sorted.bam"), file("${sampName}_sorted.bam.bai") into ch_outBAMSorted
   
   publishDir path: './raw_bam', mode: 'copy'
   
    module		samtoolsModule
    
   script:
   """
   samtools sort -@ ${task.cpus} \
       -T temp_${sampName} \
       -o ${sampName}_sorted.bam \
       -O bam \
       ${bam} 
   """
}

process picardMarkDup {
	
   label 'medium_6h'

   input:
     set sampName, file(bam) from ch_outBAMSorted

   output:
     set sampName, file("${sampName}_sorted.mdups.bam"), file("${sampName}_sorted.mdups.bam.bai"), file("${sampName}_sorted.mdups.metrics") into ch_outMdups
     set sampName, file("${sampName}_sorted.mdups.bam"), file("${sampName}_sorted.mdups.bam.bai"), file("${sampName}_sorted.mdups.metrics") into ch_outMdups2
     set sampName, file("${sampName}_sorted.mdups.bam"), file("${sampName}_sorted.mdups.bam.bai"), file("${sampName}_sorted.mdups.metrics") into ch_outMdups3
     set sampName, file("${sampName}_sorted.mdups.bam"), file("${sampName}_sorted.mdups.bam.bai"), file("${sampName}_sorted.mdups.metrics") into ch_outMdups4
     set sampName, file("${sampName}_sorted.mdups.bam"), file("${sampName}_sorted.mdups.bam.bai"), file("${sampName}_sorted.mdups.metrics") into ch_outMdups5
     set sampName, file("${sampName}_sorted.mdups.bam"), file("${sampName}_sorted.mdups.bam.bai"), file("${sampName}_sorted.mdups.metrics") into ch_outMdups6
     set sampName, file("${sampName}_sorted.mdups.bam"), file("${sampName}_sorted.mdups.bam.bai"), file("${sampName}_sorted.mdups.metrics") into ch_outMdups7
     set sampName, file("${sampName}_sorted.mdups.bam"), file("${sampName}_sorted.mdups.bam.bai"), file("${sampName}_sorted.mdups.metrics") into ch_outMdups8
     set sampName, file("${sampName}_sorted.mdups.bam"), file("${sampName}_sorted.mdups.bam.bai"), file("${sampName}_sorted.mdups.metrics") into ch_outMdups9

   publishDir path: './bamFiles', mode: 'copy'

    module      picardModule

   script:
   """
   picard MarkDuplicates \
       TMP_DIR=${sampName}_tmp \
       VALIDATION_STRINGENCY=LENIENT \
       CREATE_INDEX=true \
       INPUT=${bam} \
       OUTPUT=${sampName}_sorted.mdups.bam \
       METRICS_FILE=${sampName}_sorted.mdups.metrics
   """
}

process bam_stats1 {
   label 'fastqc'
	
   input:
     set sampName, file(bam), file(bai) from ch_outMdups

   output:
      set sampName, file("${sampName}_sorted.mdups.flagstat") into ch_outSAMStats
   
   publishDir path: '.qc_out/samtools', mode: 'copy'
   
    module		samtoolsModule
    
   script:
   """
   samtools flagstat -@ ${task.cpus} ${bam} > ${sampName}_sorted.mdups.flagstat
   """
}

process bam_stats2 {
   label 'fastqc'
	
   input:
     set sampName, file(bam), file(bai) from ch_outMdups2

   output:
      set sampName, file("${sampName}_sorted.mdups.idxstats") into ch_outSAMStats2
   
   publishDir path: '.qc_out/samtools', mode: 'copy'
   
    module		samtoolsModule
    
   script:
   """
   samtools idxstats ${bam} > ${sampName}_sorted.mdups.idxstats
   """
}

process bam_stats3 {
   label 'fastqc'
	
   input:
     set sampName, file(bam), file(bai) from ch_outMdups3

   output:
      set sampName, file("${sampName}_sorted.mdups.stats") into ch_outSAMStats3
   
   publishDir path: '.qc_out/samtools', mode: 'copy'
   
    module		samtoolsModule
    
   script:
   """
   samtools stats -@ ${task.cpus} ${bam} > ${sampName}_sorted.mdups.stats
   """
}

process bam_qc1 {
	
   label 'medium_1h'

   input:
     set sampName, file(bam), file(bai) from ch_outMdups4

   output:
     set sampName, file("${sampName}_alignment_metrics.txt") into ch_outQC

   publishDir path: './qc_out/picard', mode: 'copy'

    module      picardModule
    module      samtoolsModule
    module      RModule

   script:
   """
   picard CollectAlignmentSummaryMetrics \
     VALIDATION_STRINGENCY=LENIENT \
     I="${bam}"  \
     O="${sampName}_alignment_metrics.txt"  \
     R="${ref}" 
   """
}

process bam_qc2 {
	
   label 'medium_1h'

   input:
     set sampName, file(bam), file(bai) from ch_outMdups5

   output:
     set sampName, file("${sampName}_sorted_mdups_insert_size.metrics"), file("${sampName}.insertSize.pdf") into ch_outQC2

   publishDir path: './qc_out/picard', mode: 'copy'

    module      picardModule
    module      samtoolsModule
    module      RModule

   script:
   """
   picard CollectInsertSizeMetrics \
     VALIDATION_STRINGENCY=LENIENT \
     I="${bam}"  \
     O="${sampName}_sorted_mdups_insert_size.metrics"  \
     HISTOGRAM_FILE="${sampName}.insertSize.pdf"
     R="${ref}" 
   """
}

process bam_qc3 {
	
   label 'medium_1h'

   input:
     set sampName, file(bam), file(bai) from ch_outMdups6

   output:
     set sampName, file("${sampName}_sorted_mdups_gc.metrics"), file("${sampName}_sorted_mdups_gc.pdf"), file("${sampName}_sorted_mdups_gc_summary.metrics") into ch_outQC3

   publishDir path: './qc_out/picard', mode: 'copy'

    module      picardModule
    module      samtoolsModule
    module      RModule

   script:
   """
   picard CollectGcBiasMetrics \
     VALIDATION_STRINGENCY=LENIENT   \
     I="${bam}"  \
     O="${sampName}_sorted_mdups_gc.metrics"  \
     CHART_OUTPUT="${sampName}_sorted_mdups_gc.pdf" \
     SUMMARY_OUTPUT="${sampName}_sorted_mdups_gc_summary.metrics" \
     R="${ref}"
   """
}

process bedtools {
	
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


