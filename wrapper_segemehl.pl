#!/usr/bin/perl

#==================#
#   HEADER START   #
#==================#
### Name: wrapper_segemehl.pl
### Created: Apr 12, 2013
### Modified: Apr 12, 2013
### Author: Christina Herrmann
### Company: Zavolan Group, Biozentrum, University of Basel
### Version: v0.1
### Adapted from: wrapper_sra_to_fastq.pl
### Requirements: n/a
#==================#
### Description: Writes and submits jobs for segemehl
### Output: Job files
### Usage: perl ./wrapper_segemehl.pl; for information on required and optional arguments
#==================#
#    HEADER END    #
#==================#


#==========================#
#   PRE-REQUISITES START   #
#==========================#
#---> PRAGMAS / PACKAGES <---#
use strict;
use warnings;
#==========================#
#    PRE-REQUISITES END    #
#==========================#


#================#
#   MAIN START   #
#================#
#---> MAIN VARIABLES <---#
my $job_name_base_pl = '/import/bc2/home/zavolan/cherrmann/ALTSPLICE/GSE12946/jobs/segemehl/segemehl';

#---> BODY <---#
# Open directory
opendir DIR, '/import/bc2/home/zavolan/cherrmann/ALTSPLICE/GSE12946/input/fastq';
# Read FASTQ files from directory
my @files = map {s/\.[^.]+$//;$_} grep {/\.fastq$/} readdir DIR;
# Close directory
closedir DIR;
# Traverse through @files array
foreach my $file (@files) {
	# Print job
	&print_job($file);
}

#---> STATUS MESSAGE <---#
print "\nDone.\n";

#---> PROGRAM EXIT <---#
exit;
#================#
#    MAIN END    #
#================#


#=======================#
#   SUBROUTINES START   #
#=======================#
sub print_job {
#---> PASS ARGUMENTS ---#
my $file_base = shift;

#---> BODY <---#
# Build job name
my $job_name_pl =  $job_name_base_pl . '_' . $file_base . '.job';
# Open output file
open OUT, ">$job_name_pl";

## Print job
print OUT
'#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q fs_long@@x3755
#$ -j y
#$ -l mem_total=64000M
#$ -cwd
#$ -l sjpn=1
#$ -o /import/bc2/home/zavolan/cherrmann/ALTSPLICE/GSE12946/jobs/segemehl
# previous line specifies location for job log files

## GENERAL INFO
# Author: Alexander Kanitz
# Created: 15-APR-2013
# Modified: 15-APR-2013

## STATUS MESSAGE
echo "JOB STARTED."

## VARIABLES
# File basename
base=' . $file_base . '
# Paths to script files
s_dir=/import/bc2/home/zavolan/grubera/snoRNA-paper/split-mapping/segemehl
# Genome file
gen=/import/bc2/home/zavolan/GROUP/tmp/BWA_benchmarking/indices/hg19.fa
# Index file
idx=/import/bc2/home/zavolan/GROUP/tmp/BWA_benchmarking/indices/index_hg19_seg/index_hg19.idx
# Path to input file
i_dir=/import/bc2/home/zavolan/cherrmann/ALTSPLICE/GSE12946/input/fastq
# Path to output files
o_dir=/import/bc2/home/zavolan/cherrmann/ALTSPLICE/GSE12946/input/sam

## COMMAND
# Map sequences
time $s_dir/segemehl.x -i $idx -d $gen -q $i_dir/$base.fastq --differences 1 --accuracy 85 --threads 1 --silent -o $o_dir/segemehl_${base}.sam -u $o_dir/segemehl_${base}_unmatched

## STATUS MESSAGE
echo "JOB COMPLETED."
	
';	

# Close output file
close OUT;

# Submit job
system("qsub", $job_name_pl);
}
#=======================#
#    SUBROUTINES END    #
#=======================#
