#!/usr/bin/perl

#==================#
#   HEADER START   #
#==================#
### Name: wrapper_test.pl
### Created: Mar 28, 2013
### Modified: Mar 28, 2013
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
### Version: v0.1
### Adapted from: n/a
### Requirements: n/a
#==================#
### Description: Writes and submits jobs for summarize_overlaps_genes.R
### Output: Job files for SGE queue
### Usage: perl ./wrapper_test.pl
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
my $job_name_base_pl = '/import/bc2/home/zavolan/kanitz/AS/GSE41716/sra/jobs/summarize_overlaps_genes';
my $queue_pl = 'fs_long@@x3755';
my $job_dir_pl = "'/import/bc2/home/zavolan/kanitz/AS/GSE41716/sra/jobs'";
my $date_cr_pl = '28-MAR-2013';
my $date_mod_pl = '28-MAR-2013';
my $r_lib_dir_pl = "'/import/bc2/home/zavolan/kanitz/R/R-2.15.2/library'";
my $Rscript_pl = "'/import/bc2/home/zavolan/kanitz/R/R-2.15.2/bin/Rscript'";
my $s_dir_pl = "'/import/bc2/home/zavolan/kanitz/AS/GSE41716/sra/scripts'";
my $bam_dir_pl = "'/import/bc2/home/zavolan/kanitz/AS/GSE41716/sra/input/bam";
my $bed_dir_pl = "'/import/bc2/home/zavolan/kanitz/resources'";
my $gen_pl = "'hg19'";
my $mode_pl = "'IntersectionStrict'";
my $o_dir_pl = "'/import/bc2/home/zavolan/kanitz/AS/GSE41716/sra/output'";

#---> BODY <---#
# Traverse through each BAM dir
for my $dir (1..25) {
	# Build job name
	my $job_name_pl =  $job_name_base_pl . $dir . '.job';
	# Open output file
	open OUT, ">$job_name_pl";
	# Print job
	&print_job($dir);
	# Close output file
	close OUT;
	# Submit job
	system("qsub", $job_name_pl);
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
my $dir = shift;
#---> BODY <---#
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

## EXPORT
# Export R library path
export R_LIBS_USER=' . $r_lib_dir_pl . '

## VARIABLES
# Path to Rscript
Rscript=' . $Rscript_pl . '
# Path to script
s_dir=' . $s_dir_pl . '
# Path to BAM files
bam_dir=' . $bam_dir_pl . '/' . $dir . "'" . '
# Path to BED file
bed_dir=' . $bed_dir_pl . '
# Genome information
gen=' . $gen_pl . '
# Overlap mode
mode=' . $mode_pl . '
# Path to output file
o_dir=' . $o_dir_pl . '

## COMMANDS
time $Rscript $s_dir/summarize_overlaps_genes.R $bam_dir $bed_dir/hg19_exon_coordinates.bed $gen $mode $o_dir/summarized_overlaps_hg19_genes_abyzov_' . $dir . '.R

## STATUS MESSAGE
echo "JOB COMPLETED."';	
}
#=======================#
#    SUBROUTINES END    #
#=======================#

