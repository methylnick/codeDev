#!/bin/bash

# This script runs RNASik pipeline v1.4.7
# Assumes there is a samplesheet.txt file
# Results are the works

# Load RNAsik Module
ml load RNAsik-pipe/1.4.7

#Run the pipeline, paired end RNA Seq data 24 CPUS/threads
RNAsik -align star \
  -fastaRef refFiles/genome.fa \
  -fqDir fastq/ \
  -count \
  -gtfFile refFiles/genes.gtf \
  -prePro \
  -samplesSheet samplesheet.txt \
  -fastqc \
  -multiqc \
  -exonicRate \
  -extn "_L001_R1.fastq.gz" \
  -threads 24 

mail -s "RNAsik Analysis Complete" nick.wong@monash.edu <<< 'RNASik v1.4.7 is complete'
