#!/usr/bin/perl

### Author: Alexander Kanitz
### Created: 05-FEB-2013
### Modified: 07-FEB-2013
### Adapted from: N/A 
### Description: Determines which genomic crosslink positions (as determined by genome site extraction based on mutations on CLIPZ) overlap with regions of interest indicated in a BED file
### Arguments: 1. Crosslink site file (from CLIPZ analysis); 2. BED file with regions of interest; 3. output file
### Output: Output file (tab-delimited) indicating feature name, feature chromosome, feature start, feature end, feature strand, start best crosslink, end best crosslink, rank best crosslink, score best crosslink, total crosslink score, number of crosslinks
### Usage: perl ./clipz_overlap_cl_sites_bed.pl cl_sites.tab regions_of_interest.bed output_file


### A. Pre-requisites
## Pragmas
use warnings;
use strict;
## Command-line arguments / initialization
my $cl_sites_file = "/home/kanitz/Dropbox/Work/Eclipse/general/input/1288_cl_sites_C2T";#shift;
my $bed_file = "/home/kanitz/Dropbox/Work/Eclipse/general/input/mm9_exon_coordinates.bed";#shift;
my $out_file = "/home/kanitz/Dropbox/Work/Eclipse/general/output/1288_cl_sites_in_exons";#shift;
###


### B. Main program 
# Load CLIPZ crosslinking sites file to hash of arrays
my $cl_sites_HoAoA_ref = &tabToHashOfArraysOfArrays($cl_sites_file, 9, 1, 1);
# Crawl through target BED file and determine overlaps; update window position counts in hash of arrays of arrays
my $results_array_ref = &overlapBEDwithCLIPZclSites($cl_sites_HoAoA_ref, $bed_file);
# Print output file
&printArrayRef($results_array_ref, $out_file, "\n");
## Status message
print "Done.\n";
## Exit program
exit;
###


### C. Subroutines
sub tabToHashOfArraysOfArrays {
### Accepts: 1. tab-delimited file; 2. total number of columns; 3. number of grouping column; 4. 0 (no header) or 1 (header present)
### Returns: Hash of arrays of arrays, grouped by indicated column (hash); each row constitutes an entry to an (outer) array which holds the other column values (inner array) in the order they appear
	## Pass arguments
	my $tab_file = shift;
	my $cols = shift;
	my $first = shift;
	my $header = shift;
	# Declare/initialize variables
	my $HoAoA_ref;
	# Open bed file handle
	open TAB, $tab_file;
	# Discard header if present
	<TAB> if $header;
	## Crawl through tab file line by line
	while (<TAB>) {
		# Remove trailing newline character
		chomp;
		# Split line into array @line by tabs
		my @line = split /\s+/;
		# Skip line, if number of columns does not match specified number
		next if @line != $cols;
		# Extract grouping column
		my $key = splice @line, ($first-1), 1;
		# Initiate outer array if not yet present 
		$$HoAoA_ref{$key} = [] unless exists $$HoAoA_ref{$key};
		# Push values in line to hash of arrays of arrays		
		push $$HoAoA_ref{$key}, [@line];
	}
	# Close bed file handle
	close TAB;
	# Print status message
	print "File '$tab_file' processed.\n";
	# Return bed file hash of hashes of arrays
	return $HoAoA_ref;
}
sub overlapBEDwithCLIPZclSites {
### Function: Crawl through target BED file line by line and determine overlap with query CLIPZ crosslink sites
### Accepts: 1. crosslink site file hash (chromosomes) of arrays (lines of file) of arrays (start, end, strand, score, rank, class, mutations, copies) from CLIPZ genome site extraction based on crosslinking enrichment, 2. BED file
### Returns: Array reference of results: feature name, feature chromosome, feature start, feature end, feature strand, start best crosslink, end best crosslink, rank best crosslink, score best crosslink, total crosslink score, number of crosslinks; separated by TAB
	## Arguments / declarations
	my $cl_sites_HoAoA_ref = shift;
	my $bed_file = shift;
	# Declare results array variable
	my @results;
	# Build header
	my $header = 	"Feature name"			.	"\t"	.
					"Chromosome feature"	.	"\t"	.
					"Start feature"			.	"\t"	.
					"End feature"			.	"\t"	.
					"Strand feature"		.	"\t"	.
					"Start best crosslink"	.	"\t"	.
					"End best crosslink"	.	"\t"	.
					"Best crosslink rank"	.	"\t"	.
					"Best crosslink score"	.	"\t"	.
					"Total crosslink score"	.	"\t"	.					
					"Number of crosslinks";					
	# Push header to results array
	push @results, $header;					
	# Open positions file handle
	open BED, "$bed_file";
	## Crawl through fasta file line by line
	while (<BED>) {
		# Remove trailing new line characters
		chomp;
		# Split line into relevant coordinates by tabs
		my ($chr, $start, $end, $name, $score, $strand) = split /\s+/;
		# Declare variables that keep track of lowest crosslink rank, lowest crosslink score, total crosslink score and number of crosslinks 
		my ($min_rank, $min_score, $total_score, $no_overlaps, $start_best_cl, $end_best_cl) = 0; 
		## Crawl through each crosslink site
		foreach my $array_ref (@{$$cl_sites_HoAoA_ref{$chr}}) {
			# Get crosslink site coordinates, score, rank and class from current array
			my $start_cl	=	$$array_ref[0];
			my $end_cl		=	$$array_ref[1];
			my $strand_cl	=	$$array_ref[2];
			my $score_cl	=	$$array_ref[3];
			my $rank_cl		=	$$array_ref[4];
			my $class_cl	=	$$array_ref[5];
			## IF the current crosslinking site...
			if (
				$start		<	$start_cl	&&		## ...overlaps with the current region of interest AND
				$start		<	$end_cl		&&
				$end		>	$start_cl	&&
				$end		>	$end_cl		&&
				$strand		eq	$strand_cl	&&		# ...is on the same strand as the current region of interest AND
				$class_cl	!~	m/Not/i				# ...is not classified as "Not crosslinked" ...
				)
			{	# ...assign/adjust lowest crosslink rank, lowest crosslink score,  and 
				# ...increase number of crosslinks
				$no_overlaps++;
				# ...adjust total crosslink score
				$total_score += $score_cl;
				# ...and adjust variables for best crosslink if current crosslink site has lowest rank so far
				if ($rank_cl < $min_rank || $min_rank == 0) {
					$min_rank = $rank_cl;
					$min_score = $score_cl;
					$start_best_cl = $start_cl;
					$end_best_cl = $end_cl;
				} 
			}
		}
		## IF at least one overlapping crosslink site was found for the current region of interest... 
		if ($no_overlaps) {
			# ...build line for output file
			my $line = 	$name			.	"\t"	.
						$chr			.	"\t"	.
						$start			.	"\t"	.
						$end			.	"\t"	.
						$strand			.	"\t"	.
						$start_best_cl	.	"\t"	.
						$end_best_cl	.	"\t"	.
						$min_rank		.	"\t"	.
						$min_score		.	"\t"	.
						$total_score	.	"\t"	.
			 			$no_overlaps;
			# ...push line to results array
			push @results, $line;		
		}
	}
	# Close fasta file handle
	close BED;
	# Screen output
	print "File '$bed_file' processed.\n";
	# Return results array reference
	return \@results;
}
sub printArrayRef {
### Accepts: 1. array reference; 2. output file name; 3. separator
### Returns: The referenced array, separated by the specified separator, is written to the indicated output file
	## Pass arguments
	my $array_ref = shift;
	my $out_file = shift;
	my $sep = shift;
	# Open output file handle
	open OUT, ">$out_file";
	## Print each element of referenced array to output file
	foreach my $element (@$results_array_ref) {
		print OUT $element . $sep;
	}
	# Screen output
	print "Output written to file '$out_file'.\n";
	# Close output file handle
	close OUT;
}