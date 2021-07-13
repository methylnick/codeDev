#!/usr/bin/env nextflow

/*
 * 2021-03-31 Working up and testing the nf script for production
 * This is a test and dev script so it is going straight to the short queue
 * Include AMELX loss of Y processing and vardict calling 
 * BAM filter employed to extract alignments most likely to be associated with
 * assay proper. Vardict calling the insertion deployed to call Loss of Y.
 * also now runs VEP on the samples properly with singularity container and
 * stops when samples/files with no variants called (an empty .tsv file) results
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
sidVardictAmp  = "${refFolder}/sid_v2_amplicon_bed_vardict_20210628.bed.intervals.bed"

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
Channel.fromFilePairs("fastqs/*_R{1,2}_001.fastq.gz")
  .set{ ch_fastaIn }
  
Channel.fromFilePairs("fastqs/*_R{1,2}_001.fastq.gz")
  .set{ ch_fastQC }

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
     set sampName, file("${sampName}-trimmed-pair1.fastq.gz"), file("${sampName}-trimmed-pair2.fastq.gz") into (ch_fastaToBwa, ch_fastaToFastqc, ch_loyBwa)
     
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
     set sampName, file("${sampName}.sorted.bam"), file("${sampName}.sorted.bam.bai") into (ch_mappedBams, ch_mappedBams2, ch_mappedBams3, ch_mappedBams4, ch_mappedBams5, ch_mappedBams6, ch_mappedBams7)
     
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

process loy_bwa {

   label 'bwa_genomics'

   input:
     set sampName, file(fq1), file(fq2) from ch_loyBwa

   output:
     set sampName, file("${sampName}.bam") into (ch_loyBams, ch_loyBams2)

   publishDir path: './loy/out_bam', mode: 'copy'

    module      bwaModule
    module      samtoolsModule

   script:
   """
   bwa mem -t ${task.cpus} -R "@RG\\tID:${sampName}\\tPU:${sampName}\\tSM:${sampName}\\tPL:ILLUMINA\\tLB:rhAmpSeq" \
       ${loyRefBase} ${fq1} ${fq2} | samtools view -u -h -q 1 - \
       | samtools sort -@ ${task.cpus} -o "${sampName}.bam"
   """
}


process loy_bam {

   label 'genomics_qc'

   input:
     set sampName, file(bam) from ch_loyBams

   output:
     set sampName, file("${sampName}_filtered.bam"), file("${sampName}_filtered.bam.bai") into (ch_loyVardict)

   publishDir path: './loy/out_bam', mode: 'copy'

    module      samtoolsModule

   shell:
   '''
    samtools view -h !{bam} | awk 'substr($0,1,1)=="@" || ($9>=103 && $9<=105) || ($9<=-103 && $9>=-105)' | \
       samtools view -b > !{sampName}_104.bam
    samtools view -h !{bam} | awk 'substr($0,1,1)=="@" || ($9>=132 && $9<=134) || ($9<=-132 && $9>=-134)' | \
       samtools view -b > !{sampName}_133.bam
    samtools view -h !{bam} | awk 'substr($0,1,1)=="@" || ($9>=127 && $9<=129) || ($9<=-127 && $9>=-129)' | \
       samtools view -b > !{sampName}_128.bam
    samtools view -h !{bam} | awk 'substr($0,1,1)=="@" || ($9>=144 && $9<=145) || ($9<=-144 && $9>=-145)' | \
       samtools view -b > !{sampName}_144.bam
    samtools merge -f -h !{sampName}_104.bam !{sampName}_filtered.bam !{sampName}_104.bam \
       !{sampName}_133.bam !{sampName}_128.bam !{sampName}_144.bam
    samtools index !{sampName}_filtered.bam
   '''
}

process loy_bam2 {

   label 'genomics_qc'

   input:
     set sampName, file(bam) from ch_loyBams2

   output:
     set sampName, file("${sampName}_filtered.bam"), file("${sampName}_filtered.bam.bai") into (ch_loyVardict2)

   publishDir path: './loy/out_bam_133', mode: 'copy'

    module      samtoolsModule

   shell:
   '''
    samtools view -h !{bam} | awk 'substr($0,1,1)=="@" || ($9>=132 && $9<=134) || ($9<=-132 && $9>=-134)' | \
       samtools view -b > !{sampName}_133.bam
    samtools view -h !{bam} | awk 'substr($0,1,1)=="@" || ($9>=127 && $9<=129) || ($9<=-127 && $9>=-129)' | \
       samtools view -b > !{sampName}_128.bam
    samtools merge -f -h !{sampName}_133.bam !{sampName}_filtered.bam !{sampName}_133.bam !{sampName}_128.bam
    samtools index !{sampName}_filtered.bam
   '''
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
    
    label 'vardict_genomics'
    
    input:
      set sampName, file(bam), file(bai) from ch_mappedBams5
	  
    output:
      set sampName, file("${sampName}.tsv") into ch_vcfMake
	  
    publishDir path: './chip/raw_variants/', mode: 'copy'

    script:
    """
    module purge
    module load samtools
    export PATH=/home/nwong/bin/VarDict-1.8.0/bin:$PATH
    VarDict -G ${ref} -f ${AF_THR} -N ${sampName} -b ${bam} -th ${task.cpus} \
      ${vardictAmp}  \
      >  "${sampName}.tsv"
    """
}

process gatk {
    
    label 'vardict_genomics'
    
    input:
      set sampName, file(bam), file(bai) from ch_mappedBams7
	  
    output:
      set sampName, file("${sampName}.g.vcf.gz"), file("${sampName}.g.vcf.gz.tbi") into ch_gatkOut
	  
    publishDir path: './sampID/gvcfs/', mode: 'copy'

    script:
    """
    module purge
    module load samtools
    export PATH=/home/nwong/bin/gatk-4.2.0.0:$PATH
    gatk HaplotypeCaller  \
     -R ${ref} \
     -L ${sidVardictAmp} \
     -I ${bam} \
     -O ${sampName}.g.vcf.gz \
     -ERC GVCF
    """
}

process gatk_gt {
    
    label 'vardict_genomics'
    
    input:
      set sampName, file(vcf), file(tbi) from ch_gatkOut
	  
    output:
      set sampName, file("${sampName}.vcf.gz"), file("${sampName}.vcf.gz.tbi") into ch_gatkGTOut
	  
    publishDir path: './sampID/vcfs/', mode: 'copy'

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
    cat ${tsv} | /home/nwong/bin/VarDict-1.8.0/bin/teststrandbias.R | \
        /home/nwong/bin/VarDict-1.8.0/bin/var2vcf_valid.pl -N "${sampName}" \
        -f ${AF_THR} -E > "${sampName}.vardict.vcf"
    """
}

process vep {

    label 'vep'

    containerOptions "-B ${refFolder}/vep:/opt/vep/.vep -B ${refFolder}:/refFolder"

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
        --fasta /refFolder/genome.fa \
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

process loy_vardict {

    label 'vardict_loy'

    input:
      set sampName, file(bam), file(bai) from ch_loyVardict

    output:
      set sampName, file("${sampName}.tsv") into ch_loyVcf

    publishDir path: './loy/raw_variants/', mode: 'copy'

    script:
    """
    module purge
    module load samtools
    export PATH=/home/nwong/bin/VarDict-1.8.0/bin:$PATH
    VarDict -G ${ref} -f ${AF_THR} -N ${sampName} -b ${bam} -th ${task.cpus} \
      -c 1 -S 2 -E 3 -g 4 ${loyVardictAmp}  \
      >  "${sampName}.tsv"
    """
}

process loy_vardict2 {

    label 'vardict_loy'

    input:
      set sampName, file(bam), file(bai) from ch_loyVardict2

    output:
      set sampName, file("${sampName}.tsv") into ch_loyVcf2

    publishDir path: './loy/raw_variants_133', mode: 'copy'

    script:
    """
    module purge
    module load samtools
    export PATH=/home/nwong/bin/VarDict-1.8.0/bin:$PATH
    VarDict -G ${ref} -f ${AF_THR} -N ${sampName} -b ${bam} -th ${task.cpus} \
      -c 1 -S 2 -E 3 -g 4 ${loyVardictAmp}  \
      >  "${sampName}.tsv"
    """
}

process loymakeVCF {

    label 'vardict_loy'

    input:
        set sampName, file(tsv) from ch_loyVcf.filter{ sampName, file -> !file.isEmpty() }
    output:
        set sampName, file("${sampName}.vardict.vcf") into ch_loyDone

    publishDir path: './loy/vardict', mode: 'copy'

    module RModule

    script:
    """
    module purge
    module load R/3.6.0-mkl
    cat ${tsv} | /home/nwong/bin/VarDict-1.8.0/bin/teststrandbias.R | \
        /home/nwong/bin/VarDict-1.7.0/bin/var2vcf_valid.pl -N "${sampName}" \
        -f ${AF_THR} -E > "${sampName}.vardict.vcf"
    """
}

process loymakeVCF2 {

    label 'vardict_loy'

    input:
        set sampName, file(tsv) from ch_loyVcf2.filter{ sampName, file -> !file.isEmpty() }
    output:
        set sampName, file("${sampName}.vardict.vcf") into ch_loyDone2

    publishDir path: './loy/vardict_133', mode: 'copy'

    module RModule

    script:
    """
    module purge
    module load R/3.6.0-mkl
    cat ${tsv} | /home/nwong/bin/VarDict-1.8.0/bin/teststrandbias.R | \
        /home/nwong/bin/VarDict-1.7.0/bin/var2vcf_valid.pl -N "${sampName}" \
        -f ${AF_THR} -E > "${sampName}.vardict.vcf"
    """
}

