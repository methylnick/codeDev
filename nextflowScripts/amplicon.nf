#!/usr/bin/env nextflow

/*
 * This is my first script in processing amplicon seq data
 * Created on 10 July 2019
 *
 */

/*
 * Declare parameters
 */
 
params.in = "fastq/*.fastq.gz"
fqfiles = Channel
                 .fromPath(params.in)
                 .map { file -> tuple(file.simpleName, file)}
                 
fqfiles.into { fqfiles_readsplit; fqfiles_fastqc; fqfiles_align}
results_path = \$PWD
fqPairs = Channel 
				 .fromFilePairs('fastq/*_{R1,R2}.fastq.gz', 'fastq/*_{1,2}.fq.gz')
				 .println()


/*
 * Operate on the fastq file to extract machine and flow cell
 * information for populating the RG line in a SAM output
 */
process extractFqHeader {
	input:
	file fastq from fqfiles_readsplit
	
	output:
	stdout readout
	val fastq into header
	
	script:
	"""
	zcat $fastq | head -1 | awk '{print \$1}'
	"""
}

readout.subscribe { print "This is read header $it" }

process extractFields {
	input:
	val head from header
	
	output:
	stdout headout
	
	script:
	"""
	IFS=: read ins run flo lan til xps yps <<< $head
	echo \$read
	echo \$ins
	echo \$run
	echo \$flo
	echo \$lan
	echo \$til
	echo \$xps
	echo \$yps
	"""
}

headout.subscribe { print "These are the values extracted $it"}



/*
 * run fastq on raw files
 */
 process Fastqc {
	 publishDir "$results_path/nf_fastqc/", mode: 'copy', pattern: "_fastqc.html"
	 publishDir "$results_path/nf_fastqc/", mode: 'copy', pattern: "_fastqc.zip"
	 
	 input:
	 file input from fqfiles_fastqc
	 
	 output:
	 file "*_fastqc.html" into fastqcReport
	 file "*_fastqc.zip" into fastqcArchive

	 script:
	 """
	 fastqc $input
	 """
 }
 
