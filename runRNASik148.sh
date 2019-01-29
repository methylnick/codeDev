#!/bin/bash

# This script runs RNASik v1.4.8

# Load RNAsik Module
ml load RNAsik-pipe/1.4.8

#Run the pipeline, paired end RNA Seq data
RNAsik -align star \
  -fastaRef referenceFiles/genome.fa \
  -fqDir fastq/ \
   -counts \
   -gtfFile referenceFiles/genes.gtf \
   -prePro \
   -fastqc \
   -multiqc \
   -exonicRate \
   -extn "_001.fastq.gz" \
   -threads 24

mail -s "RNAsik Analysis Complete" nick.wong@monash.edu <<< 'RNASik v1.4.8 is complete'





