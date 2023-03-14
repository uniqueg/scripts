#!/usr/bin/perl

#######
### AUTHOR INFO & CREATION DATE:
### ----------------------------
### Author:		Alexander Kanitz
### Created: 	15-JAN-2013
### Modified:	15-JAN-2013
#######

#######
### FUNCTION:
### ---------
### The script converts a single file in oligomap mapping format into a tab-separated value file that can be analyzed with 'compare_tsv.pl'.
#######

#######
### ARGUMENTS:
### ----------
### 1. Input file name (including full path if not in current working directory)
### General Notes: Compare the examples in USAGE for further information.
#######

#######
### OUTPUT:
### -------
### For each entry in the oligomap mapping file, one line per alignment with the following entries is written to <FILE>.tsv (in the same directory as the input file): 1. sequence ID, 2. chromosome, 3. strand, 4. start position, 5. end position, 6. sequence.
#######

#######
### USAGE:
### ------
### perl /path/to/oligomap2tsv.pl /path/to/gem_mapping_file
### Examples:
### 1. perl ~/oligomap2tsv.pl ~/test.map (converts test.map to test.map.tsv; mapping file and script in home directory)
### 2. perl oligomap2tsv.pl test.map (converts test.map to test.map.tsv; mapping file and script in current working directory) 
#######

#######
### OVERVIEW:
### ---------
### A. Pre-requisites
### B. File processing
### C. Clean-up
#######

### A. Pre-requisites
## Pragmas
use warnings;
use strict;
# Command-line argument
my $file = shift;
###

### B. Main program
# Open input file handle
open IN, $file;
# Open output file handle
open OUT, ">$file.tsv";
## Crawl through input file entry by entry
while (<IN>) {
	# Remove trailing newline character
	chomp;
	# Add line to dedicated variable
	my $line = $_;
	# Test whether line is first line of alignment entry
	if ($_ =~ m/^seq/) {
		my @line = split(/\t/, $line);
		my ($id) = split(/ /, $line[0], 2);
		my $chr = $line[1];
		my ($start, $stop) = split(/\.\./, $line[2]);
		my $str = <IN>;
		chomp $str;
		$str =~ s/errors(.)*orientation:.//;
		my $seq = <IN>;
		chomp $seq;
		print OUT "$id\t$chr\t$str\t$start\t$stop\t$seq\n";
	}
}
# Close input file handle
close OUT;
# Close output file handle
close IN;
###

### C. Clean-up
# Exit
exit;
###