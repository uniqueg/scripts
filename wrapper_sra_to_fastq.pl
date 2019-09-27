#!/usr/bin/perl

#==================#
#   HEADER START   #
#==================#
### Name: wrapper_sra_to_fastq.pl
### Created: Apr 12, 2013
### Modified: Apr 12, 2013
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
### Version: v0.1
### Adapted from: wrapper_summarize_overlaps.pl
### Requirements: n/a
#==================#
### Description: Writes and submits jobs for sra fastq-dump
### Output: Job files
### Usage: perl ./wrapper_sra_to_fastq.pl; for information on required and optional arguments
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
my $job_name_base_pl = '/import/bc2/home/zavolan/cherrmann/ALTSPLICE/GSE12946/jobs/sra_to_fastq';
my $queue_pl = 'fs_short';
my $job_dir_pl = '/import/bc2/home/zavolan/cherrmann/ALTSPLICE/GSE12946/jobs';
my $date_cr_pl = '12-APR-2013';
my $date_mod_pl = '12-APR-2013';
my $s_dir_pl = '/import/bc2/home/zavolan/bilebi00/bin';
my $i_dir_pl = '/import/bc2/home/zavolan/cherrmann/ALTSPLICE/GSE12946/input/sra';
my $o_dir_pl = '/import/bc2/home/zavolan/cherrmann/ALTSPLICE/GSE12946/input/fastq';

#---> BODY <---#
# Open directory
opendir DIR, $i_dir_pl;
# Read SRA files from directory
my @files = map {s/\.[^.]+$//;$_} grep {/\.sra$/} readdir DIR;
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
#$ -q ' . $queue_pl . '
#$ -j y
#$ -cwd
#$ -o ' . $job_dir_pl . '
	
## GENERAL INFO
# Author: Alexander Kanitz
# Created: ' . $date_cr_pl . '
# Modified: ' . $date_mod_pl . '

## STATUS MESSAGE
echo "JOB STARTED."

## VARIABLES
# Path to script
s_dir=' . $s_dir_pl . '
# Path to input files
i_dir=' . $i_dir_pl . '
# Path to output files
o_dir=' . $o_dir_pl . '

## COMMANDS
time $s_dir/fastq-dump --split-spot --skip-technical --minReadLen 25 --outdir $o_dir/ $i_dir/' . $file_base . '.sra

## STATUS MESSAGE
echo "JOB COMPLETED."';	

# Close output file
close OUT;

# Submit job
system("qsub", $job_name_pl);
}
#=======================#
#    SUBROUTINES END    #
#=======================#

