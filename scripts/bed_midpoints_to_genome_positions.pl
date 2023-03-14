#!/usr/bin/perl

### Author: Alexander Kanitz
### Created: 30-JAN-2013
### Modified: 30-JAN-2013
### Description: Maps the midpoints of regions in a BED-file to genome positions; strand information disregarded
### Arguments: Accepts a BED file and an output file name
### Output: Returns a table indicating the chromosome name and midpoint of each region in the BED file
### Usage: perl bed_midpoints_to_genome_positions.pl input.bed output.tab

### A. Pre-requisites
## Pragmas
use strict;
use warnings;
## Command-line arguments / initialization
my $bed = "s_8_sequence_BARC4_trim.hg19-hgV02-aln2.3.m100.bedweight.shifted";#shift;
my $out = "s_8_sequence_BARC4_trim.out";#shift;
###

### B. Main program
# Load BED file to hash of arrays
my $bed_HoA_ref = &bedMidRegionToHashOfSortedArrays($bed);
# Calculate and print continuous genome positions regions in BED file
&bedGenomePositions($bed_HoA_ref, $out);
# Exit
exit;
###

### C. Subroutines
sub bedMidRegionToHashOfSortedArrays {
### Accepts file in bed format
### Returns a hash of of arrays containing region midpoints grouped by chromosome
### Arrays are sorted by size in ascending order
	## Pass argument
	my $bed_file = shift;
	# Declare/initialize variables
	my %HoA;
	# Open bed file handle
	open BED, $bed_file;
	## Crawl through bed file line by line
	while (<BED>) {
		# Remove trailing new line characters
		chomp;
		# Split line into array @line by tabs
		my @line = split /\t/, $_;
		# Extract chromosome (= key for outer hash) from @line array
		my $chr = $line[0];
		# Calculate midpoint of region from start/end positions
		my $mid = ($line[1] + $line[2]) / 2;
		## IF hash entry with key $chr does not exist, initialize it
		unless (exists $HoA{$chr}) {
			## Add value for each key of inner hash
   			$HoA{$chr} = [ $mid ];
		}
		## ELSE push to existing values (i.e. arrays) of hash
		else {
			## For each key of inner hash, push values to existing ones (array!)
   			push @{$HoA{$chr}}, $mid;
		}
	}
	# Close bed file handle
	close BED;
	## Sort each array of %HoA
	foreach my $array_ref (values %HoA) {
		@{$array_ref} = sort { $a <=> $b } @{$array_ref};
	}
	# Print status message
	print "File $bed_file processed.\n";
	# Return bed file hash of hashes of arrays
	return \%HoA;
}
sub bedGenomePositions {
### Writes continuous genome positions for each region to the indicated output file
### Requires a reference to a hash of arrays containing mid-points of regions sorted by chromosome, a reference to a hash containing continuous genomic starting positions for each chromosome and an output file name 
	## Pass arguments
	my %bed = %{shift()};
	my $out = shift;
	# Print status message
	print "Calculating positions...\n";
	# Open output file handle 
	open OUT, ">$out";
	## Crawl through chromosomes (= outer hash of HoA) one by one
	foreach my $chr (keys %bed) {
		## Crawl through values in HoA arrays (i.e. regions in target bed file)
		foreach my $region (@{$bed{$chr}}) {	
			# Print chromosome and region midpoint to output file
			print OUT "$chr\t$region\n";
		}
	}
	# Close output file handle
	close OUT;
	# Screen output
	print "Positions written to file: '$out'\n";
}
###