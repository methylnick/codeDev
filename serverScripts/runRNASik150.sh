#!/bin/bash

# This script runs RNASik v1.5.0 with the works

ml load RNAsik-pipe/1.5.0

RNAsik -fqDir fastq \
       -align star  \
       -fastaRef referenceFiles/genome.fa  \
       -gtfFile referenceFiles/genes.gtf  \
       -counts  \
       -all

mail -s "RNAsik Analysis Complete" nick.wong@monash.edu <<< 'RNASik v1.5.0 is complete'