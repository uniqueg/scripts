#!/usr/bin/perl

### Author: Alexander Kanitz
### Created: 15-JAN-2013/ Modified: 15-JAN-2013
### Description: Remove from a tab file of patterns and distances (output from distance_bed_bed.pl / concatenate_files.pl) all patterns that do have more or less than 2 CG di-nucleotides 
### Usage: perl /path/to/remove_all_but_2CG.pl </path/to/input_file> <path/to/output_file>

### A. Pre-requisites
## Pragmas
use strict;
use warnings;
## Command-line arguments / initialization
my $in_file = shift;
my $out_file = shift;
###

### B. Main program
# Read and write files
&readWrite($in_file, $out_file);
# Exit
exit;
###

### C. Subroutines

sub readWrite {
### Reads text line by line, transforms each line and re-writes to output file if transformation successful
### Reuse: Replace function called for line transformation according to needs
	## Pass arguments
	my $in = shift;
	my $out = shift;
	# Open output file
	open OUT, ">$out";
	# Open input file
	open IN, $in;
	# Read input file line by line
	while (<IN>) {
		# Line transformation
		my $transformed_line = &removeUnwantedPatterns($_);
		# Print if transformation successful
		print OUT $transformed_line if $transformed_line;
	}
	# Close input bed file
	close IN;
	# Close output file
	close OUT;
}

sub removeUnwantedPatterns {
### Concatenates files in array and writes to output file
	# Pass arguments
	my $line = shift;
	# Calculate number of CG matches
	my $count =()= $line =~ /CG/gi;
	# Traverse through @files array
	if ($count == 2) {
		return $line;
	}
	# Return FALSE/0
	return 0;
}