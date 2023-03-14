#!/usr/bin/perl

#######
### AUTHOR INFO & CREATION DATE:
### ----------------------------
### Author:		Alexander Kanitz
### Created: 	07-JAN-2013
### Modified:	15-JAN-2013
#######

#######
### FUNCTION:
### ---------
### The program merges overlapping or nearby regions of a BED file. A genome file in FASTA format is required to extract the lengths of chromosomes. The ranges over which regions are merged are indicated by the 'flanking parameter', e.g. by setting it to 100 will result in merging of regions that are 200 nucleotides or less apart. Note that only the first three columns of the BED file (chromosome, start and stop positions) are regarded, all other information is discarded and NOT included in the output. 
#######

#######
### ARGUMENTS:
### ----------
### 1. Input file name (include path if not in current working directory)
### 2. Genome file name (include path if not in current working directory)
### 3. Flanking region (integer)
### 4. Output file name (path optional)
### General Notes: Compare the examples in USAGE for further information.
#######

#######
### OUTPUT:
### -------
### A BED file containing chromosome, start and stop positions of the merged regions. The value of the flanking parameter is subtracted/added from/to the start/stop positions (min start position = 1; max stop position = last nt of chromosome).
#######

#######
### USAGE:
### ------
### perl /path/to/merge_regions_from_bed.pl <path/to/bed/file> <path/to/genome/fasta/file> <flanking parameter/integer> <path/to/output/file>
### Examples:
### 1. perl ./merge_regions_from_bed.pl ./test.bed ./hg19.fa 0	out.bed (only overlapping regions are merged)
### 2. perl ./merge_regions_from_bed.pl ./test.bed ./hg19.fa 500 out_flank500.bed (regions are flanked with 500 nts upstream and downstream before merging)
#######

#######
### OVERVIEW:
### ---------
### A. Pre-requisites
### B. Get chromosome sizes
### C. Read input bed file
### D. Merge regions and print to file
### E. Clean-up
#######

### A. Pre-requisites
use warnings;
use strict;
my $in_file = shift;
my $gen_file = shift;
my $flank = shift;
my $out_file = shift;
my %chr_size;
my %chr_reg;
###

### B. Get chromosome sizes			## BETTER USE hg19_sizes.fa & sub fastaToHash (see distance_bed_bed)
# Open genome file handle
open GEN, "$gen_file";
# Initialize chromosome identifier variable $id
my $id;
# Crawl through genome fasta file line by line
while (<GEN>) {
	# Remove trailing new line character from input line
	chomp;
	## If line starts with ">" (identifier line), extract identifier and initialize new element of hash %chr_size with 0 length; else, add length of line to value of the latest key 
	if (substr($_, 0, 1) eq ">") {
		$id = substr $_, 1;
		$chr_size{$id} = 0;
	}
	else {
		$chr_size{$id} = $chr_size{$id} + length($_);
	}
}
# Close genome file handle
close GEN;
# Print status message
print "Chromosome sizes lookup table created.\n";
###

### C. Read input bed file
# Open input file handle
open IN, "$in_file";
# Crawl through input BED file line by line
while (<IN>) {
	# Split lines by tabulators into chromosome ($chr), start ($start) and stop ($end) positions; remaining information is discarded
	my ($chr, $start, $end) = split /\t/, $_, 4;
	## Start and stop positions of each region in the BED file are stored in hash %chr_reg (keys = $chr; values = start/stop positions of each region; start and stop values within a region are separated by ':', between regions by ',')
	## If hash entry of the current chromosome $chr exists, append; else, initialize new element
	if (exists $chr_reg{$chr}) {
		$chr_reg{$chr} = $chr_reg{$chr} . ",$start:$end";
	}
	else {
		$chr_reg{$chr} = "$start:$end";
	}
}
# Close input file handle
close IN;
# Print status message
print "BED file read.\n";
###

### D. Merge regions and print to file
# Open output file handle
open OUT, ">$out_file";
# Crawl through each chromosome (hash %chr_reg)
foreach my $key (keys %chr_reg) {
	# Split different regions into array @reg_orig
	my @reg_orig = split(/,/, $chr_reg{$key});
	# Sort array by start position
	my @reg_sort = map { $_->[0] } sort { $a->[1] <=> $b->[1] } map { [$_, /(\d+):/] } @reg_orig;
	# Split coordinates of first entry (lowest start position) of the current chromosome into array @coord
	my @coord = split(/:/, shift @reg_sort);
	# Crawl through remaining entries of the current chromosome
	foreach (@reg_sort) {
		# Split coordinates of next region (next highest start position) into array @tmp
		my @tmp = split(/:/, $_);
		## If the start position of the current (@tmp) region (minus the flanking parameter if applicable) is bigger than the stop position of the previous (@coord) region (plus the flanking parameter if applicable), print previous region and then replace with current region; else, replace stop position of previous region with that of current region (-> extension of region!)    
		if ($tmp[0] - $flank > $coord[1] + $flank) {
			# Set start position to 1 if out of bounds because of flanking parameter
			$coord[0] = ($coord[0] - $flank) < 1 ? 1 : $coord[0] - $flank;
			# Set start position to last nucleotide of chromosome if out of bounds because of flanking parameter
			$coord[1] = ($coord[1] + $flank) > $chr_size{$key} ? $chr_size{$key} : $coord[1] + $flank;
			print OUT "$key\t$coord[0]\t$coord[1]\n";
			@coord = @tmp;
		}
		else {
			$coord[1] = $tmp[1];
		}
	}
	## Print last entry of each chromosome
	# Set start position to 1 if out of bounds because of flanking parameter 
	$coord[0] = ($coord[0] - $flank) < 1 ? 1 : $coord[0] - $flank;
	# Set start position to last nucleotide of chromosome if out of bounds because of flanking parameter
	$coord[1] = ($coord[1] + $flank) > $chr_size{$key} ? $chr_size{$key} : $coord[1] + $flank;
	print OUT "$key\t$coord[0]\t$coord[1]\n";
}
# Close output file handle
close OUT;
# Print status message
print "Overlapping regions merged.\n";
###

### E. Clean-up
# Print status message
print "Results written to file: $out_file\n";
# Exit program
exit;
###