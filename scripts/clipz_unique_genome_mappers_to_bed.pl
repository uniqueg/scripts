#!/usr/bin/perl

### Author: Alexander Kanitz
### Created: 24-JAN-2013
### Modified: 24-JAN-2013
### Adapted from: ag-filter-mappings-from-clipz.pl by Andreas R. Gruber
### Description: The scripts evaluates genome mappings of a CLIPZ experiment, filters reads and writes relevant information to a BED file. Reads are disregarded if they are not unique genome mappers or map against bacterial, fungal, viral, adaptor- or vector-derived sequences
### Arguments: 1. CLIPZ experiment file "mapped_sequences"; 2. CLIPZ experiment file "genome_mappings"; 3. output file
### Output: Tab-delimited BED file containing the following columns: chr, chr_start, chr_end, seq_id, count, strand
### Usage: perl clipz_unique_genome_mappers_to_bed.pl mapped_sequences genome_mappings unique_genome_mappers.bed

### A. Pre-requisites
## Pragmas
use strict;
use warnings;
## Command-line arguments / initialization
my $mapped_sequences = shift;
my $genome_mappings = shift;
my $out_file = shift;
###

### B. Main program
# 
parseMappedSequences($mapped_sequences, $genome_mappings, $out_file);
# Print status message
print "Done.\n";
# Exit
exit;
###

### C. Subroutines
sub parseMappedSequences { 
###
	# Pass arguments
	my $mapped_sequences = shift;
	my $genome_mappings = shift;
	my $out_file = shift;
	# Declare hash variable to store relevant info
	my $count_ref;
	# Print status message
	print "Processing file '$mapped_sequences'...\n";
	# Open filehandle for mapped reads
	open MAP, "$mapped_sequences";
	## Crawl through mapped reads file line by line
	while (<MAP>) {
		# Remove trailing newline character
		chomp;
		# Skip header line
		next if $_ =~ m/^id/;
		# Split line by tab into @columns array
		my @cols = split /\t/;
		## Skip if read is...
		next if $cols[-1]					||		# ...putative adaptor
				$cols[-2] eq 'bacterial'	||		# ...bacterial
				$cols[-2] eq 'fungus'		||		# ...fungal
				$cols[-2] eq 'vector'		||		# ...part of a vector
				$cols[-2] eq 'viral'		||		# ...viral
				$cols[-5] ne '1';					# ...not a unique mapper
		# Save count
		$$count_ref{$cols[0]} = $cols[2];
		# Write to file every 10^6 lines
		$count_ref = &crossRefWithGenomeSequencesFileAndWriteToBed($count_ref, $genome_mappings, $out_file) if scalar keys %$count_ref > 1000000;
	}
	# Close filehandle for mapped reads
	close MAP;
	# Write remaining lines to file
	&crossRefWithGenomeSequencesFileAndWriteToBed($count_ref, $genome_mappings, $out_file);
}
sub crossRefWithGenomeSequencesFileAndWriteToBed {
###
	# Pass arguments
	my $count_ref = shift; 
	my $genome_mappings = shift;
	my $out_file = shift;
	# Print status message
	print "\tCross-referencing batch of reads against file '$genome_mappings' and writing to file '$out_file'...\n";
	# Open filehandle for reads mapped against the genome
	open GEN, "$genome_mappings";
	# Open output filehandle and enable appending
	open OUT, ">>$out_file";
	## Crawl through reads mapped against the genome file line by line
	while (<GEN>) {
		# Remove trailing newline character
		chomp;
		# Skip header line
		next if $_ =~ m/^id/;
		# Split line by tab into @columns array
		my @cols = split /\t/;
		# Skip line if no count/annotation info available for current read
		next unless exists $$count_ref{$cols[1]};
		## For each count
		foreach (1..$$count_ref{$cols[1]}) {
			# Print relevant info to BED file $out_file
			print OUT "$cols[2]\t$cols[6]\t$cols[7]\tseq$cols[1]\t1\t$cols[5]\n";
		}
	}
	# Close output filehandle
	close OUT;
	# Close filehandle for reads mapped against the genome
	close(GEN);
	# Empty hash
	%$count_ref = ();
	# Return hash reference
	return $count_ref;
}
###