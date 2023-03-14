#!/usr/bin/perl

#######
### AUTHOR INFO & CREATION DATE:
### ----------------------------
### Author:		Alexander Kanitz
### Created: 	10-JAN-2013
### Modified:	14-JAN-2013
#######

#######
### FUNCTION:
### ---------
### Search for one or more patterns across one or more sequences (e.g. genome). Both patterns and sequences have to be in fasta file format. 
#######

#######
### ARGUMENTS:
### ----------
### 1. Input pattern file in fasta format
### 2. Input sequence file in fasta format
### 3. Reverse complement required? 0 for NO, any other integer for YES 
### 4. Output file
### General Notes: Compare the examples in USAGE for further information.
#######

#######
### OUTPUT:
### -------
### A BED file listing each occurrence of each pattern in the pattern file (with or without reverse complement) across all sequenced in the sequence file (identifier line of sequence file assumed to contain only chromsome information) 
#######

#######
### USAGE:
### ------
### perl /path/to/pat_search.pl <path/to/pat_fa_file> </path/to/seq_fa_file> <reverse complement/integer> <path/to/output_file>
### Examples:
### 1. perl pat_search.pl patterns.fa sequences.fa 0 out.bed (reverse complement not considered) 
### 2. perl pat_search.pl patterns.fa sequences.fa 1 out.bed (adds reverse complement for each pattern in patterns.fa)
#######

#######
### OVERVIEW:
### ---------
### A. Pre-requisites
### B. Main program
### C. Subroutines
#######


### A. Pre-requisites

## Pragmas
use strict;
use warnings;
## Command-line/initialization variables
my $in_pat = shift;
my $in_seq = shift;
my $rc = shift;
my $out_file = shift;


### B. Main program

# Read pattern fasta file to array
my $pat_array_ref = &fastaToArray($in_pat);
# Add reverse complements for each pattern if desired ($rc <> 0)
$pat_array_ref = &revCompl($pat_array_ref) if ($rc);
# Read sequence fasta file entry by entry, search for patterns and write to output file (bed file if sequence fasta ID line contains only chromosome!)
&patInGenome($pat_array_ref, $in_seq, $out_file);
# Exit program
exit;


### C. Subroutines

sub fastaToArray {
### Pass file in fasta format and load into array
### Array elements: Sequences (multiple lines are concatenated, separating newline characters - if present - are removed) 
	# Pass arguments
	my $file = shift;
	# Screen output
	print "Processing fasta file \"$file\".\n";
	# Declare initialize variables
	my @array;
	my $seq = "";
	# Open fasta file handle
	open FA, "$file";
	## Crawl through fasta file line by line
	while (<FA>) {
		# Remove trailing new line characters
		chomp;
		# Pass line to dedicated variable
		my $line = $_;
		# Execute code block IF line is identifier (contains leading '>')
		if ($line =~ s/^>//) {
			# Unless empty string, push sequence to array
			push @array, $seq if $seq ne "";
			# Empty sequence string
			$seq = "";
		}
		## ELSE concatenate line to existing sequence 
		else {
			$seq = $seq.$line;
		}
	}
	# Push last entry to array
	push @array, $seq;
	# Close fasta file handle
	close FA;
	# Screen output
	print "File \"$file\" processed.\n";
	# Return array reference
	return \@array;
}

sub revCompl {
### Pass array reference containing patterns
### For each entry, add another entry with the reverse complement (suffix "_comp" added to key/id)
### Only A, C, G and T (both upper- and lowercase) are considered
	# Pass arguments
	my @pat = @{shift()};
	my $rev;
	my $rev_compl;
	## Crawl through all pattern entries
	foreach (@pat) {
		# Compute reverse of current pattern
		$rev = reverse($_);
		# Compute complement
		$rev_compl =~ tr/AaCcGgTt/TtGgCcAa/;
		# Push to array @pat
		push @pat, $rev_compl;
	}
	# Return array reference
	return \@pat;
}

sub patInGenome {
### Search for one or more patterns across one or more sequences
### Patterns are passed as an array reference
### An output file in BED format is written (chrom, chrom start, chrom end, found pattern)
### For proper output, sequence file should contain relevant/meaningful IDs (i.e. chromosome name)
	## Pass arguments
	my @pat = @{shift()};
	my $seq_file = shift();
	my $out_file = shift();
	my $seq = "";
	my $chr;
	# Open output file handle
	open OUT, ">$out_file";
	# Open sequence file handle;
	open SEQ, "$seq_file";
	## Traverse through all sequences in sequence file
	while (<SEQ>) {
		# Remove trailing new line characters
		chomp;
		# Pass line to dedicated variable
		my $line = $_;
		# Execute code block if identifier line (contains leading '>')
		if ($line =~ s/^>//) {
			## Unless empty string, search in sequence
			if ($seq ne "") {
				# Traverse through all patterns in array
				foreach (@pat) {
					# Search for all occurrences of current pattern
					while($seq =~ /$_/ig) {
						# For each found occurence, write chromosome, start/end position and "ID:found pattern" to output file
						print OUT $chr . "\t" . (pos($seq) - length($&) + 1) . "\t" . pos($seq) . "\t" . "$&" . "\n";
					}
				}
				# Screen output
				print "Chromosome $chr completed.\n";					
			}		
			# Set identifier variable ($id)
			$chr = $line;
			# Empty sequence string
			$seq = "";
		}
		## Else concatenate line to existing sequence 
		else {
			$seq = $seq.$line;
		}
	}
	## Search last sequence
	# Traverse through all patterns in array
	foreach (@pat) {
		# Search for all occurrences of current pattern
		while($seq =~ /$_/ig) {
			# For each found occurence, write chromosome, start/end position and "ID:found pattern" to output file
			print OUT $chr . "\t" . (pos($seq) - length($&) + 1) . "\t" . pos($seq) . "\t" . "$&" . "\n";
		}
	}
	# Screen output
	print "Chromosome $chr completed.\n";			
	# Close sequence file handle;
	close SEQ;
	# Close output file handle
	close OUT;
	# Screen output
	print "Search completed.\n";
}