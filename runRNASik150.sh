#!/bin/bash

ml load RNAsik-pipe/1.5.0

RNAsik -fqDir fastq \
       -align star  \
       -fastaRef referenceFiles/genome.fa  \
       -gtfFile referenceFiles/genes.gtf  \
       -counts  \
       -all
