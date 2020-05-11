#!/usr/bin/env nextflow

/*
 * RNASeq NextFlow Pipeline
 * Human GRCh38 Ensembl Reference
 */

// Declare Inputs
refFolder = file("/home/nwong/df22_scratch/references/iGenomes")
inputDirectory = file('fastq')

// Declare References 
starRef        = "/home/nwong/df22_scratch/nick.wong/2020-05-10-PanCancer/starIndex"
ref            = "${refFolder}/Homo_sapiens/Ensembl/GRCh38/Sequence/WholeGenomeFasta/genome.fa"
genes         = "${refFolder}/Homo_sapiens/Ensembl/GRCh38/Annotation/Genes/genes.gtf"

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
     set sampName, file("${sampName}_Aligned.out.bam"), file("${sampName}_Log.final.out"), file("${sampName}_Log.out"), file("${sampName}_Log.progress.out"), file("${sampName}_SJ.out.tab")  into ch_mappedBams
     

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
   label 'bwa'
	
   input:
     set sampName, file(bam), file(Log1), file(Log2), file(Log3), file(Log4) from ch_mappedBams

   output:
      set sampName, file("${sampName}_sorted.bam") into ch_outBAMSorted
   
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
	
   label 'bwa'

   input:
     set sampName, file(bam) from ch_outBAMSorted

   output:
     set sampName, file("${sampName}_sorted.mdups.bam"), file("${sampName}_sorted.mdups.metrics") into ch_outMdups
     set sampName, file("${sampName}_sorted.mdups.bam"), file("${sampName}_sorted.mdups.metrics") into ch_outMdups2
     set sampName, file("${sampName}_sorted.mdups.bam"), file("${sampName}_sorted.mdups.metrics") into ch_outMdups3
     set sampName, file("${sampName}_sorted.mdups.bam"), file("${sampName}_sorted.mdups.metrics") into ch_outMdups4
     set sampName, file("${sampName}_sorted.mdups.bam"), file("${sampName}_sorted.mdups.metrics") into ch_outMdups5
     set sampName, file("${sampName}_sorted.mdups.bam"), file("${sampName}_sorted.mdups.metrics") into ch_outMdups6
     set sampName, file("${sampName}_sorted.mdups.bam"), file("${sampName}_sorted.mdups.metrics") into ch_outMdups7
     set sampName, file("${sampName}_sorted.mdups.bam"), file("${sampName}_sorted.mdups.metrics") into ch_outMdups8
     set sampName, file("${sampName}_sorted.mdups.bam"), file("${sampName}_sorted.mdups.metrics") into ch_outMdups9
     set sampName, file("${sampName}_sorted.mdups.bam"), file("${sampName}_sorted.mdups.metrics") into ch_outMdups10
     set sampName, file("${sampName}_sorted.mdups.bam"), file("${sampName}_sorted.mdups.metrics") into ch_outMdups11

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
	  set sampName, file(bam), file(bai) from ch_outMdups7
	  
	output:
	  set sampName, file("${sampName}.bedGraph") into ch_bigWig
	  
	publishDir path: './coverageFiles', mode: 'copy'
	
    module      bedtoolsModule
    
    script:
    """
    bedtools genomecov -bga \
      -split \
      -ibam ${bam} \
      -g ${ref} > ${sampName}.bedGraph.tmp
    LC_COLLATE=C  sort -k1,1 -k2,2n ${sampName}.bedGraph.tmp > ${sampName}.bedGraph
    rm ${sampName}.bedGraph.tmp
    """
	
}

process qualimap {
	
	label 'medium_1h'
	
	input:
	  set sampName, file(bam), file(bai) from ch_outMdups8
	  
	output:
	  set sampName, file("${sampName}/*") into ch_outQualimap
	  
	publishDir path: './qc_out/qualimap/${sampName}', mode: 'copy'
    
    script:
    """
    export PATH="$PATH:/home/nwong/bin/qualimap_v2.2.1"
    qualimap rnaseq -bam ${bam} \
       --paired \
       --sorted \
       -gtf ${genes} \
       -outdir ${sampName} \
       -oc ${sampName}_counts.txt \
       --sequencing-protocol strand-specific-reverse \
       --java-mem-size=16G
    """	
}

process fc_ustrand {
	
	label 'vardict'
	
	input:
	  set sampName, file(bam), file(bai) from ch_outMdups9.collect()
	  
	output:
	  set sampName, file("*.txt") into ch_outFeatureCounts1
	  
	publishDir path: './counts', mode: 'copy'
	
	module      subreadModule
    
    script:
    """
    featureCounts -T ${task.cpus} \
       -s 0 \
       -a ${genes} \
       -p \
       -o NonStrandedCounts.txt \
       ${bam}
    """
	
}

process fc_strand {
	
	label 'vardict'
	
	input:
	  set sampName, file(bam), file(bai) from ch_outMdups10.collect()
	  
	output:
	  set sampName, file("*.txt") into ch_outFeatureCounts2
	  
	publishDir path: './counts', mode: 'copy'
	
	module      subreadModule
    
    script:
    """
    featureCounts -T ${task.cpus} \
       -s 1 \
       -a ${genes} \
       -p \
       -o StrandedCounts.txt \
       ${bam}
    """
	
}

process fc_revStrand {
	
	label 'vardict'
	
	input:
	  set sampName, file(bam), file(bai) from ch_outMdups11.collect()
	  
	output:
	  set sampName, file("*.txt") into ch_outFeatureCounts3
	  
	publishDir path: './counts', mode: 'copy'
	
	module      subreadModule
    
    script:
    """
    featureCounts -T ${task.cpus} \
       -s 2 \
       -a ${genes} \
       -p \
       -o ReverseStrandedCounts.txt \
       ${bam}
    """
}
