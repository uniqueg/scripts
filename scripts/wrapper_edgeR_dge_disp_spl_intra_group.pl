#!/usr/bin/perl

#==================#
#   HEADER START   #
#==================#
### Name: wrapper_edgeR_dge_disp_splice.pl
### Created: Apr 25, 2013
### Modified: Apr 28, 2013
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
### Version: v0.1
### Adapted from: n/a
### Requirements: n/a
#==================#
### Description: Writes and submits jobs for edgeR_dge_disp_splice.R
### Output: Job files for SGE queue
### Usage: perl ./wrapper_edgeR_dge_disp_splice.pl
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
my $in_file = '/import/bc2/home/zavolan/cherrmann/ALTSPLICE/GSE12946/input/other/inter_group_comparisons';
my $job_name_base = '/import/bc2/home/zavolan/cherrmann/ALTSPLICE/GSE12946/jobs/edgeR/';

#---> BODY <---#
# Generate 1-line input files from multiline input file
my $names_array_ref = &one_file_per_line($in_file);
# Initialize serial number for unique output filename generation
my $sn = 1;
## Traverse through each element $file of array @$names_array_ref
foreach my $file ( @$names_array_ref ) {
    # Generate job filename
    my $job_file = $job_name_base . 'edgeR_dge_disp_splice.job.' . $sn;
    # Print job with current file
	&print_job($file, $job_file);
	# Submit job
	system "qsub", $job_file;
	# Increase serial number count
	$sn++;
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
	my ($file, $job_file) = @_;
	#---> BODY <---#
	# Open output file handle
	open OUT, ">$job_file";
	# Print job
	print OUT
'#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q fs_long@@x3755
#$ -j y
#$ -cwd
#$ -o /import/bc2/home/zavolan/cherrmann/ALTSPLICE/GSE12946/jobs/edgeR

## GENERAL INFO
# Author: Alexander Kanitz
# Created: 25-APR-2013
# Modified: 25-APR-2013

## STATUS MESSAGE
echo "JOB STARTED."

## EXPORT
# Export R library path
export R_LIBS_USER=/import/bc2/home/zavolan/kanitz/R/R-2.15.2/library

## VARIABLES
# Path to Rscript
Rscript="/import/bc2/home/zavolan/kanitz/R/R-2.15.2/bin/Rscript"
# Path to script file
s_dir="/import/bc2/home/zavolan/cherrmann/ALTSPLICE/GSE12946/scripts"
# P value cutoff
p_value=0.05
# Output folder
o_dir="/import/bc2/home/zavolan/cherrmann/ALTSPLICE/GSE12946/output/edgeR"
# Prefix for output files
o_prefix="first_run_"

## COMMAND
time $Rscript $s_dir/edgeR_dge_disp_splice.R ' . $file . ' $p_value $splice $o_dir $o_prefix

## STATUS MESSAGE
echo "JOB COMPLETED."';
	# Close output file handle
	close OUT;	
}
#-----------------------#
sub one_file_per_line {
### Function: Reads file line by line and writes one file for each line; new line characters are removed; a serial number '.#' is appended to the original filename (with # from 1..number of rows in input file)
### Accepts: 1. Path to text file [STRING]
### Returns: 1. Reference to array containing the names of the generated files
### Dependencies: n/a
### Type: Generic
	## Pass arguments
	my $in_file = shift;
	# Declare array of output filenames
	my @out_names;
	# Initialize serial number for unique output filename generation
	my $sn = 1;
	# Open input file handle
	open IN, $in_file;
	## Traverse through input file line by line
	while (<IN>) {
		# Assign line content to dedicated variable
		my $line = $_;
		# Remove trailing newline character
		chomp $line;
		# Generate output filename
		my $out_file = $in_file . '.' . $sn;
		# Add output filename to array of output filenames
		push @out_names, $out_file;
		# Open output file handle
		open OUT, ">$out_file";
		# Write line to output file handle
		print OUT $line;
		# Increase serial number count
		$sn++;
		# Close output file handle
		close OUT;	
	}
	# Close output file handle
	close OUT;
	# Return array reference
	return \@out_names;
}
#=======================#
#    SUBROUTINES END    #
#=======================#