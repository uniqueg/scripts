#!/usr/bin/perl

### Author: Alexander Kanitz
### Created: 14-JAN-2013/ Modified: 14-JAN-2013
### Description: The program writes a fasta file indicating the sizes of sequences in a fasta file   
### Usage: perl /path/to/fasta_seq_to_fasta_size.pl ./in.fa ./out.fa

### A. Pre-requisites
## Pragmas
use warnings;
use strict;
## Command-line arguments / initialization
my $in_file = shift;
my $out_file = shift;
###

### B. Main program 
# Load input FASTA file to hash
my $in_hash_ref = &fastaToHash($in_file);
# Calculate sequence sizes
$in_hash_ref = &seqSizesHash($in_hash_ref);
# Print sizes to output FASTA file
&printFastaHash($in_hash_ref, $out_file);
# Exit program
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

sub seqSizesHash {
	# Pass arguments
	my %seq_hash = %{shift()};
	# Screen output
	print "Calculating sizes... ";
	## Loop through sequence hash sorted by keys (sequence IDs)
	foreach my $seq_id (sort keys %seq_hash) {
		$seq_hash{$seq_id} = length $seq_hash{$seq_id};
	}
	# Screen output
	print "Done.\n";
	# Return hash reference
	return \%seq_hash;			
}

sub printFastaHash {
	# Pass arguments
	my %seq_hash = %{shift()};
	my $out_file = shift;
	# Open output file handle
	open OUT, ">$out_file";
	## Loop through sequence hash sorted by keys (sequence IDs)
	foreach my $seq_id (sort keys %seq_hash) {
		print OUT ">$seq_id\n$seq_hash{$seq_id}\n";
	}
	# Close output file handle
	close OUT;
	# Screen output
	print "Sequence sizes written to file \"$out_file\".\n";
}
