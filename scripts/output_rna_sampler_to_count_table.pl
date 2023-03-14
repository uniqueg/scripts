#!/usr/bin/perl

#==================#
#   HEADER START   #
#==================#
### Name: output_rna_sampler_to_count_table.pl
### Created: Jul 26, 2013
### Modified: Jul 26, 2013
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
### Version: v0.1
### Adapted from: n/a
### Requirements: n/a
#==================#
### Description: Converts the output from the 'rna' tool written by Cedrid Gobet and Felix Naef (to estimate isoform expression) into a count table format
### Output: The output is of the form: Feature ID \t Count (= transcript fraction times total count of reads for the corresponding gene)
### Usage: perl ./output_rna_sampler_to_count_table.pl for information on required and optional arguments
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
my $in_file = '';
my $out_file = '';
my $options_result = GetOptions (
	'usage|help' => \$usage,
	'quiet|no_verbose' => \$quiet,
	#-----------------------#
	'in_file=s' => \$in_file,
	'out_file=s' => \$out_file
);
die $usage_info if $usage || !$options_result;
die $usage_info if !$in_file || !$out_file; 

#---> GLOBAL VARIABLES <---# 

#==========================#
#    PRE-REQUISITES END    #
#==========================#


#================#
#   MAIN START   #
#================#
#---> STATUS MESSAGE <---#
print "Starting 'output_rna_sampler_to_count_table.pl'...\n\n" unless $quiet;

#---> BODY <---#
convert_rna_output_to_count_table($in_file, $out_file);

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
'Usage: perl ./output_rna_sampler_to_count_table.pl [OPTIONS] --in_file [FILE] --out_file [FILE] --arg3_name [arg3_type]

Description: Converts the output from the "rna" tool written by Cedrid Gobet and Felix Naef (to estimate isoform expression) into a count table of the form: Feature ID \t Count (= transcript fraction times total count of reads for the corresponding gene)

==================================================
Required arguments:
--in_file	Input filename (output of program "rna") [FILE]
--out_file	Output filename [FILE]
==================================================
Optional arguments:
--usage|help	Show this information and die
--quiet	Shut up!
';
}
#-----------------------#
sub convert_rna_output_to_count_table {
### Function: Converts the output from the "rna" tool written by Cedrid Gobet and Felix Naef (to estimate isoform expression) into a count table of the form: Feature ID \t Count (= transcript fraction times total count of reads for the corresponding gene)
### Accepts: 1. Input file (output file from program 'rna'); 2. output filename
### Returns: n/a (converted file is directly written to output file)
### Dependencies: n/a
### Type: Specialized
	#---> PASS ARGUMENTS ---#
	my ($in_file, $out_file) = @_;
	
	#---> STATUS MESSAGE <---#
	print STDOUT "Converting input file..." . "\n" unless $quiet;

	#---> BODY <---#
	# Open file handle IN
	open IN, $in_file;
	
	# Open file handle OUT
	open OUT, ">$out_file";
	
	## Traverse through file handle IN line by line
	while (<IN>) {
		# Remove trailing entry separator
		chomp;
		# Next if empty string
		next if $_ eq "";
		# Split line by tab
		my ($gene_id, $trx_id, $gene_count, $trx_fraction) = split /\t/;
		# Calculate count
		my $trx_count = $trx_fraction * $gene_count;
		# Print line to output file
		print OUT $gene_id . "|" . $trx_id . "\t" . $trx_count . "\n";
	}
		
	# Close file handle OUT
	close OUT;
		
	# Close file handle IN
	close IN;
	
	#---> STATUS MESSAGE <---#
	print STDOUT "Conversion successful." . "\n\n" unless $quiet;
}
#=======================#
#    SUBROUTINES END    #
#=======================#