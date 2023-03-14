#!/usr/bin/perl

#==================#
#   HEADER START   #
#==================#
### Name: wrapper_RD_submission.pl
### Created: Mar 26, 2014
### Modified: Mar 26, 2014
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
### Version: v0.1
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
my $job_name_base_pl = '/import/bc2/home/zavolan/GROUP/ISO_BENCH/data/benchmarking_raw_output/RD/JOBS/step_2/RD_quant_step_2';
my $queue_pl = 'fs_long';
my $job_dir_pl = "'/import/bc2/home/zavolan/GROUP/ISO_BENCH/data/benchmarking_raw_output/RD/JOBS/step_2'";
my $date_cr_pl = '26-MAR-2014';
my $infile_dir = '/import/bc2/home/zavolan/GROUP/ISO_BENCH/data/benchmarking_raw_output/RD/step_1_bins';

#---> BODY <---#
# Open directory
opendir DIR, $infile_dir;
# Read files from directory
my @files = grep { /^\./ && -f "$infile_dir/$_" } readdir DIR;
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
my $outfile_base = $file_base;
$outfile_base =~ s/step_1/step_2/;

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
#$ -o ' . $job_dir_pl . '

## GENERAL INFO
# Author: Alexander Kanitz
# Created: ' . $date_cr_pl . '

## VARIABLES
# Executable
RD="/import/bc2/home/zavolan/krini/bin/Rscript /import/bc2/home/zavolan/krini/soft/RD/r-RefSeq-isoform.Rscript"
# Input filename
input="/import/bc2/home/zavolan/GROUP/ISO_BENCH/data/benchmarking_raw_output/RD/step_1_bins/' . $file_base . '"
# Output filename
output="/import/bc2/home/zavolan/GROUP/ISO_BENCH/data/benchmarking_raw_output/RD/step_2_bins/' . $outfile_base . '"
# Minimum reads per gene
threshold=100

## STATUS MESSAGE
echo "JOB STARTED."

## COMMAND
# Map sequences
time $RD --input $input --output $output --threshold $threshold

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