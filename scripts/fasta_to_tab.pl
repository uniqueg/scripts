#!/usr/bin/perl

### Author: Alexander Kanitz
### Created: 21-JAN-2013
### Modified: 21-JAN-2013
### Description: Maps the midpoints of regions in a BED-file to (arbitrarily) defined continuous genomc positions (i.e. no need for chromosome name); strand information disregarded
### Arguments: Accepts a FASTA file and an output file name
### Output: Returns a tab-separated table indicating the identifier and sequence/value
### Usage: perl fasta_to_tab.pl input.fa output.tab

### A. Pre-requisites
## Pragmas
use strict;
use warnings;
## Command-line arguments / initialization
my $fa = "hg19_sizes.fa";#shift;
my $out = "hg19_sizes.tab";#shift;
###

### B. Main program
# Load FASTA file
my $hash_ref = &fastaToHash($fa);
# Print hash to output file in tab format
&printHashToTab($hash_ref, $out);
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
	# Declare/initialize variables
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
	print "File $file processed.\n";
	# Return hash reference
	return \%hash;
}
sub printHashToTab {
### Accepts an id/sequence hash reference and prints it to the indicated tab-delimited output file (id line, seq)
	# Pass arguments
	my $hash_ref = shift;
	my $out_file = shift;
	# Open output file handle
	open OUT, ">$out_file";
	## Loop through sequence hash sorted by keys (sequence IDs)
	foreach my $id (keys %{$hash_ref}) {
		# Print in fasta format
		print OUT "$id\t${$hash_ref}{$id}\n";
	}
	# Close output file handle
	close OUT;
	# Screen output
	print "Sequence sizes written to file \"$out_file\".\n";
}
###