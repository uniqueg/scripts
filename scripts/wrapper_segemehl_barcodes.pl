#!/usr/bin/perl

#==================#
#   HEADER START   #
#==================#
### Name: wrapper_segemehl_barcodes.pl
### Created: Nov 19, 2013
### Modified: Nov 19, 2013
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
### Version: v1.0
### Adapted from: n/a
### Requirements: GetOpt::Long
#==================#
#    HEADER END    #
#==================#


#==========================#
#   PRE-REQUISITES START   #
#==========================#
#---> PRAGMAS / PACKAGES <---#
use strict;
use warnings;
use Getopt::Long;


#---> USAGE <---#
my $usage_info = &usage;

#---> OPTIONS / ARGUMENTS <---#
my $usage = '';
my $quiet = '';
my $barcodes = '';
my $out_dir = '';
my $options_result = GetOptions (
	'usage|help' => \$usage,
	'quiet' => \$quiet,
	#-----------------------#
	'barcodes=s' => \$barcodes,
	'out_dir=s' => \$out_dir
);
die $usage_info if $usage || !$options_result;
die $usage_info if !$barcodes || !$out_dir; 

#==========================#
#    PRE-REQUISITES END    #
#==========================#


#================#
#   MAIN START   #
#================#
#---> STATUS MESSAGE <---#
print "Starting '$0'...\n\n" unless $quiet;

#---> MAIN VARIABLES <---#
my $barcodes_array_ref;

#---> BODY <---#
	
	#---> Read barcodes <---#
	$barcodes_array_ref = line_to_array($barcodes);

	#---> Print and submit segemehl jobs <---#	
	&print_submit_segemehl_barcode_jobs($barcodes_array_ref);

#---> STATUS MESSAGE <---#
print "Done.\n" unless $quiet;

#---> PROGRAM EXIT <---#
exit 0;
#================#
#    MAIN END    #
#================#


#=======================#
#   SUBROUTINES START   #
#=======================#
sub usage {
### Function: Returns usage information for current script
### Accepts: n/a
### Returns: String with usage information
### Type: Specialized
'Usage: perl ./wrapper_segemehl_barcodes.pl [OPTIONS] --barcodes [FILE] --out_dir [PATH]

Description: Writes and submits SGE jobs for segemehl given a number of barcodes.

==================================================
Required arguments:
--barcodes	Text file with one barcode per line
--out_dir	Directory where jobs are written
==================================================
Optional arguments:
--usage|help	Show this information and die
--quiet	Shut up!
';
}
#-----------------------#
sub line_to_array {
### Function: Reads text file line by line and pushes each line (minus the newline character) into an array
### Accepts: 1. Path to text file [STRING]
### Returns: 1. Reference to array
### Dependencies: n/a
### Type: Generic

	#---> PASS ARGUMENTS ---#
	my $file = shift;

	#---> STATUS MESSAGE <---#
	print STDERR "Processing file '$file'..." . "\n" unless $quiet;

	#---> SUBROUTINE VARIABLES <---#
	my @lines;

	#---> BODY <---#

		#---> Open file <---#
		open FILE, "<", $file;

		#---> Push line to array <---#
		while (<FILE>) {
			chomp;
			push @lines, $_;
		}

		#---> Close file <---#
		close FILE;

	#---> STATUS MESSAGE <---#
	print STDERR "File '$file' processed." . "\n" unless $quiet;
	
	#---> RETURN VALUE <---#
	return \@lines;
	
}
#-----------------------#
sub print_submit_segemehl_barcode_jobs {

	#---> PASS ARGUMENTS ---#
	my $barcodes_array_ref = shift;

	#---> STATUS MESSAGE <---#
	print STDERR "Printing and submitting jobs..." . "\n" unless $quiet;

	#---> BODY <---#

		#---> Iterate over barcodes <---#
		foreach my $barcode (@$barcodes_array_ref) {
			
			#---> Build job name <---#
			my $job_name =  'segemehl_' . $barcode . '.job';
	
			#---> Open job file <---#
			open OUT, ">$job_name";

			#---> Print job <---#		
			print OUT '
#!/bin/bash
#$ -S /bin/bash
#$ -P project_zavolan
#$ -q fs_long@@x3755
#$ -l sjpn=1
#$ -j y
#$ -o "/import/bc2/home/zavolan/GROUP/KEI/eleni/jobs"

## VARIABLES
# Filename and basename
filename="' . $barcode . '.fa"
basename=${filename%.*}
# Segemehl executable
seg="/import/bc2/home/zavolan/kanitz/soft/bin/segemehl.x"
# Index file
idx="/import/bc2/home/zavolan/kanitz/resources/annotations/genomes/mouse/mm10/compiled/mm10_segemehl_index.idx"
# Reference file (FASTA)
ref="/import/bc2/home/zavolan/kanitz/resources/annotations/genomes/mouse/mm10/compiled/mm10.fa"
# Path to read file (FASTA|FASTQ)
reads="/import/bc2/home/zavolan/GROUP/KEI/eleni/files/time_course/processed/${filename}"
# Output file prefix
prefix="/import/bc2/home/zavolan/GROUP/KEI/eleni/files/time_course/processed/segemehl_${basename}"

## COMMAND
# Map sequences
time $seg -i $idx -d $ref -q $reads --differences 1 --accuracy 90 --threads 1 --silent -o ${prefix}.sam -u ${prefix}_unmatched
';

			#---> Close job file <---#
			close OUT;
	
			#---> Submit job file <---#
			system("qsub", $job_name);

		}

	#---> STATUS MESSAGE <---#
	print STDERR "Jobs printed and submitted." . "\n" unless $quiet;

}
#=======================#
#    SUBROUTINES END    #
#=======================#