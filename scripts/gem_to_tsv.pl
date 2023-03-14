#!/usr/bin/perl

#######
### AUTHOR INFO & CREATION DATE:
### ----------------------------
### Author:		Alexander Kanitz
### Created: 	17-DEC-2012
### Modified:	15-JAN-2013
#######

#######
### FUNCTION:
### ---------
### The script converts a single file in GEM mapping format (.map) into a tab-separated value file that can be analyzed with 'compare_tsv.pl'.
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
### For each line in the GEM mapping file, one line per alignment with the following entries is written to <FILE>.tsv (in the same directory as the input file): 1. sequence ID, 2. chromosome, 3. strand, 4. start position, 5. end position, 6. sequence.
#######

#######
### USAGE:
### ------
### perl /path/to/gem2tsv.pl /path/to/gem_mapping_file
### Examples:
### 1. perl ~/gem2tsv.pl ~/test.map (converts test.map to test.map.tsv; mapping file and script in home directory)
### 2. perl gem2tsv.pl test.map (converts test.map to test.map.tsv; mapping file and script in current working directory) 
#######

#######
### OVERVIEW:
### ---------
### A. Pre-requisites
### B. File processing
### C. Clean-up
#######

### A. Pre-requisites
use warnings;
use strict;
my $usage = 'Usage: perl /path/to/gem2tsv.pl /path/to/gem_mapping_file\n';
my $file = shift or die "$usage";
###

### B. File processing
#!# Open file, split lines, print (ID => 'rest')
open (GEM, $file);
open (OUT, ">$file.tsv");
while (<GEM>) {
	chomp;
	my @cols = split("\t");
	# if applicable, split multiple coordinates
	my @map_no = split(",", $cols[3]);
	# traverse through each set of coordinates
	foreach (@map_no) {
		# split coordinates into individual components
		my @xy = split(":");
		# calculate length
		my $len = length($cols[1]);
		# test whether coordinates are defined
		if ( $xy[0] ne "-" ) {
			# calculate end positon
			my $end = $xy[2] + $len - 1;
			print OUT "$cols[0]\t$xy[0]\t$xy[1]\t$xy[2]\t$end\t$len\n";
		}
		# if no coordinates are available, print "NA"
		else {
			print OUT "$cols[0]\tNA\tNA\tNA\tNA\t$len\n";
       	}
	}
}
close (OUT);
###

### C. Clean-up
exit;
###