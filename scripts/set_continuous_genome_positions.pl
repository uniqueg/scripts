#!/usr/bin/perl

### Author: Alexander Kanitz
### Created: 21-JAN-2013
### Modified: 21-JAN-2013
### Description: (Arbitrarily) defines continuous genomic positions for genomes; disregards strand information
### Arguments: Accepts a chromosome size table in FASTA format (chromosome name in identifier line) and an output file name
### Output: Returns a table in FASTA format indicating (arbitrary) continuous genomic starting positions of each chromosome (chromosome name in identifier line); chromsomes are sorted in the following way: chr1-22, chrX, chrY, chrM, unmapped fragments (chr known), unknown fragments (chr not known)
### Usage: perl set_continuous_genome_positions.pl hg19_sizes.fa hg19_abs_positions.fa

### A. Pre-requisites
## Pragmas
use strict;
use warnings;
## Command-line arguments / initialization
my $in = shift;
my $out = shift;
###

### B. Main program
# Obtain hash of chromosome names and sizes from input file
my $hash_ref = &fastaToHash($in);
# Sort chromosome hash
$hash_ref = &setAbsGenPosition($hash_ref);
# Print hash to output file in fasta format
&printHashToFasta($hash_ref, $out);
# Exit
exit;
###

### C. Subroutines
sub fastaToHash {
### Pass file in fasta format and load into hash
### Keys: ID lines (minus leading '>')
### Values: Sequences (multiple lines are concatenated, separating newline characters - if present - are removed) 
	# Pass arguments
	my $file = shift;
	# Screen output
	print "Processing file \"$file\".\n";
	# Declare initialize variables
	my %hash;
	my $seq = "";
	my $id;
	# Open fasta file handle
	open FA, "$file";
	## Crawl through fasta file line by line
	while (<FA>) {
		# Remove trailing new line characters
		chomp;
		# Pass line to dedicated variable
		my $line = $_;
		# Execute code block if identifier line (contains leading '>')
		if ($line =~ s/^>//) {
			# Unless empty string, pass sequence ($seq) to hash (key = $id)
			$hash{$id} = $seq if $seq ne "";
			# Set identifier variable ($id)
			$id = $line;
			# Empty sequence string
			$seq = "";
		}
		## Else concatenate line to existing sequence 
		else {
			$seq = $seq.$line;
		}
	}
	# Add last entry to hash
	$hash{$id} = $seq;
	# Close fasta file handle
	close FA;
	# Screen output
	print "File processed.\n";
	# Return hash reference
	return \%hash;
}
sub setAbsGenPosition{
	# Pass arguments
	my $hash_ref = shift;
	# Declare variables
	my $pos = 1;
	## Crawl through sorted hash keys/chromosomes	
	foreach my $chr (map { $_->[0] } sort sortChrHash2 sort sortChrHash1 map { [$_, /chr(X|Y|M|Un|\d+)(.*)/] } keys %{$hash_ref}) {
		# Assign size of current chromosome to dedicated variable
		my $chr_size = ${$hash_ref}{$chr};
		# Assign current total continuous position $pos to starting position of current chromosome
		${$hash_ref}{$chr} = $pos;
		# Add size of current chromosome to total continuous position $pos
		$pos += $chr_size;
	}
	# Return hash reference
	return $hash_ref;
}
sub printHashToFasta {
### Accepts a sequence hash and prints it to the indicated output file in the sorting order indicated by sortChrHash1 and 2
	# Pass arguments
	my $hash_ref = shift;
	my $out_file = shift;
	# Open output file handle
	open OUT, ">$out_file";
	## Loop through sequence hash sorted by keys (sequence IDs)
	foreach my $chr (map { $_->[0] } sort sortChrHash2 sort sortChrHash1 map { [$_, /chr(X|Y|M|Un|\d+)(.*)/] } keys %{$hash_ref}) {
		# Print in fasta format
		print OUT ">$chr\n${$hash_ref}{$chr}\n";
	}
	# Close output file handle
	close OUT;
	# Screen output
	print "Sequence sizes written to file \"$out_file\".\n";
}
sub sortChrHash1{
	## Sort according to capture group [1] ("1"-"22", "M", "X", "Y" or "Un")
	if 		($a->[1] eq 'Un') 	{ return  1 }						## 2a. 'Un' last
	elsif 	($b->[1] eq 'Un') 	{ return -1 }
	elsif 	($a->[1] eq 'M') 	{ return  1 }						## 2b. 'M' second last
	elsif 	($b->[1] eq 'M') 	{ return -1 }
	elsif 	($a->[1] eq 'Y') 	{ return  1 }						## 2c. 'Y' third last
	elsif 	($b->[1] eq 'Y') 	{ return -1 }
	elsif 	($a->[1] eq 'X') 	{ return  1 }						## 2d. 'X' fourth last
	elsif 	($b->[1] eq 'X') 	{ return -1 }
	else 						{ return $a->[1] <=> $b->[1] } 		#  2e. Remainder numerical
}
sub sortChrHash2{
	## Sort ASCII-betical according to capture group [2] ("" or "_xxxxxx")
	$a->[2] cmp $b->[2]
}
###