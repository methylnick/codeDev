#!/usr/bin/env nextflow

/*
 * This is my first script in processing amplicon seq data
 * Created on 10 July 2019
 *
 */



params.in = "fastq/*.fastq.gz"
fqfiles = Channel.fromPath(params.in)

/*
 * run fastq on raw files
 */
 process Fastqc {
	 input:
	 file input from fqfiles
	 
	 output: 
	 file "fastqc/${input.baseName}" into fastqcOut
	 
	 script:
	 """
	 fastqc $input
	 """
 }
 
