#!/usr/bin/perl

### Author: Alexander Kanitz
### Created: 06-FEB-2013
### Modified: 06-FEB-2013
### Adapted from: N/A
### Description: Converts mRNA/exon/intron coordinate TAB file derived from '/import/bc2/home/zavolan/rodak/auxiliaryData/processedData/*/coordinates/' and subfolders (with * being 1 for human, 2 for mouse and 3 for worms) to BED file 
### Arguments: 1. input filename; 2. output filename 
### Output: 1. BED file derived from input TAB file
### Usage: perl ./clipz_coordinates_to_bed.pl coordinates.tab coodinates.bed

### A. Pre-requisites
## Pragmas
use strict;
use warnings;
## Command-line arguments / initialization
my $in_file = shift;
my $out_file = shift;
###

### B. Main program
# 
rearrangeTab($in_file, $out_file);
# Print status message
print "Done.\n";
# Exit
exit;
###

### C. Subroutines
sub rearrangeTab {
### Accepts: 1. tab-delimited input file; 2. output file name
### Returns: 2. tab-delimited output file (rearranged)
	## Pass arguments
	my $in_file = shift;
	my $out_file = shift;
	# Open input file handle
	open IN, $in_file;
	# Open output file handle
	open OUT, ">$out_file";
	## Traverse through input file line by line
	while (<IN>) {
		# Remove trailing newline character
		chomp;
		# Split line by tab;
		my @in_line = split /\t/;
		# Rearrange
		push my @out_line, $in_line[1], $in_line[3], $in_line[4], $in_line[0], "1", $in_line[2];
		# Join output line array elements by tab
		print OUT join("\t", @out_line) . "\n";
	}
	# Close output file handle
	close OUT;
	# Close input file handle
	close IN;
	## Status messages
	print "Processed file '$in_file'.\n";
	print "Output written to file '$out_file'.\n";
}
###
