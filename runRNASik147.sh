#!/bin/bash
# Load RNAsik Module
ml load RNAsik-pipe/1.4.7

#Run the pipeline, paired end RNA Seq data
RNAsik -align star -fastaRef refFiles/genome.fa -fqDir fastq/ -count -gtfFile refFiles/genes.gtf -prePro -samplesSheet samplesheet.txt -fastqc -multiqc -exonicRate -extn "_L001_R1.fastq.gz" -threads 24 

mail -s "RNAsik Analysis Complete" nick.wong@monash.edu <<< 'Ben Shields analysis is complete'
