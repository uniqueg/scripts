#!/usr/bin/perl

### Author: Alexander Kanitz
### Created: 15-MAY-2013
### Modified: 15-MAY-2013
### Adapted from: overlap_counts_window_around_pos_vs_bed.pl 
### Description: Counts the overlaps of start sites of regions of a BED file for any position of a window of an indicated size around one or more input coordinates
### Arguments: 1. File with input coordinates (column1: chromsome, column2: coordinate); 2. reference/target BED file; 3. window size; 4. output file
### Output: Output file with one row for each query position, elements separated by SPACE; rownames: <chromosome:start_pos-end_pos>; rest: window position counts
### Usage: perl ./overlap_start_site_counts_window_around_pos_vs_bed.pl query.pos [TAB] target.bed [BED] window_size [INTEGER] output_file [FILE]

### A. Pre-requisites
## Pragmas
use strict;
use warnings;
## Command-line arguments / initialization
my $query = shift;
my $target = shift;
my $win_size = shift;
my $out_file = shift;
###

### B. Main program
# Set window around each position and fill with zeros; return hash (chromosomes) of arrays (positions) of arrays (start, end position, window position counts) reference
my $pos_HoAoA_ref = &windowAroundPositionsToHashOfSortedArraysOfArrays($query, $win_size);
# Crawl through target BED file and calculate overlaps between BED region start sites and windows; update window position counts in hash of arrays of arrays
$pos_HoAoA_ref = &overlapBEDStartWithWindowAroundPositionsToHashOfSortedArraysOfArrays($pos_HoAoA_ref, $target);
# Print results to file
&printForWindowAroundPositionsToHashOfSortedArraysOfArrays($pos_HoAoA_ref, $out_file);
# Print status message
print "Done.\n";
# Exit
exit;
###

### C. Subroutines
sub windowAroundPositionsToHashOfSortedArraysOfArrays {
### Function: Loads positions into hash (key = chromosomes) of arrays (positions) of arrays (start, end position, window position counts = 0)
### Accepts: Positions file in the following BED-derived format: chromosome TAB position 
### Returns: Hash (key = chr) of arrays (positions) of arrays (start, end positions, 0s for each window position); the outer arrays are sorted by ascending start, then end positions in ascending numerical order
	## Arguments / declaration
	my $pos_file = shift;
	my $win_size = shift;
	my %HoAoA;
	# Open positions file handle
	open POS, "$pos_file";
	## Crawl through fasta file line by line
	while (<POS>) {
		# Remove trailing new line characters
		chomp;
		# Split line into relevant coordinates by tabs
		my ($chr, $pos) = split /\t/;
		# Compute window start and end positions (round down to avoid problems with .5 'mid-point' coordinates)
		my $start = int($pos - ($win_size / 2));
		my $end = int($pos + ($win_size / 2));
		# Use current chromosome value as hash key to array of positions (outer array); push pair of start/end positions (inner array) to outer array
		push @{$HoAoA{$chr}}, [$start, $end];
		# Push one zero for each window position to inner array
	}
	# Close positions file handle
	close POS;
	## Sort and push zeros
	foreach my $AoA_ref (values %HoAoA) {
		# Sort each outer array of %HoAoA by first, then second elements in inner arrays (number sort)
		my @sortedAoA = sort { $a->[0] <=> $b->[0] || $a->[1] <=> $b->[1] } @{$AoA_ref};
		## Push one zero for each window position to inner array
		foreach my $A_ref (@sortedAoA) {
			for ($$A_ref[0]..$$A_ref[1]) {
				push @$A_ref, 0;
			}
		}
		# Assign sorted AoA reference to original AoA reference 
		$AoA_ref = \@sortedAoA;
	}
	# Screen output
	print "File '$pos_file' processed.\n";
	# Return positions file hash of arrays of arrays reference
	return \%HoAoA;
}
sub overlapBEDStartWithWindowAroundPositionsToHashOfSortedArraysOfArrays {
### Function: Crawl through target BED file line by line and determine BED region start site overlap with query positions
### Accepts: 1. positions file hash (chromosomes) of arrays (positions) of arrays (start, end position, window position counts) from sub windowAroundPositionsToHashOfSortedArraysOfArrays, 2. BED file
### Returns: Hash of arrays of arrays reference with updated window position counts
	## Arguments / declarations
	my $query_HoAoA_ref = shift;
	my $target = shift;
	# Open positions file handle
	open BED, "$target";
	## Crawl through BED file line by line
	while (<BED>) {
		# Remove trailing new line characters
		chomp;
		# Split line into relevant coordinates by tabs
		my ($chr, $start, $end, $name, $score, $strand) = split /\t/; 
		# Swap start and end coordinates if region is on minus strand
		if ($strand eq "-") {
			($start, $end) = ($end, $start);
		}
		# Declare variable that keeps track of last/previous distance
		my $last_distance;
		## Crawl through each query region array
		foreach my $array_ref (@{$$query_HoAoA_ref{$chr}}) {
			# Get query start/end positions from current array
			my $start_win = $$array_ref[0];
			my $end_win = $$array_ref[1];
			# Declare variable for overlap
			my $overlap;
			# Calculate absolute distance between target position and query region midpoint
			my $distance = &DistAbsPosFromMidRegion($start, $start_win, $end_win);
			## Determine overlap position if start of target lies between current query region
			if ($start >= $start_win && $start <= $end_win) {
				# Calculate overlap position relative to start of window
				$overlap = $start - $start_win - 1;
			}
			## Break inner loop (foreach) IF distance between target and query regions is increasing
			elsif (defined $last_distance) {
				last if $distance > $last_distance;
			}
			# Update last with current distance
			$last_distance = $distance;
			## Increase count at overlap position IF an overlap was found  (add 2 to index because first two array elements contain window start/end positions)
			$$array_ref[$overlap+2]++ if defined $overlap;
		} 
	}
	# Close fasta file handle
	close BED;
	# Screen output
	print "File '$target' processed.\n";
	# Return query hash of arrays of arrays reference, updated counts
	return $query_HoAoA_ref;
}
sub printForWindowAroundPositionsToHashOfSortedArraysOfArrays {
### Function: Prints hash (chromosomes) of arrays (positions) of arrays (start, end position, window position counts) from sub windowAroundPositionsToHashOfSortedArraysOfArrays
### Accepts: Positions file hash (chromosomes) of arrays (positions) of arrays (start, end position, window position counts) from sub windowAroundPositionsToHashOfSortedArraysOfArrays
### Returns: Output file with one row for each query position, elements separated by SPACE; rownames: <chromosome:start_pos-end_pos>; rest: window position counts
	## Arguments / declarations
	my $query_HoAoA_ref = shift;
	my $out_file = shift;
	# Open output file handle
	open OUT, ">$out_file";	
	## Crawl through outer array references
	foreach my $key (keys %$pos_HoAoA_ref) {
		## Crawl through inner array references
		foreach my $A_ref (@{$$pos_HoAoA_ref{$key}}) {
			# Print "rownames": chromosome:start_pos-end_pos
			print OUT "<$key:$$A_ref[0]-$$A_ref[1]>";
			## Discard first two array elements
			shift @$A_ref;
			shift @$A_ref;
			## Print window position counts
			foreach (@$A_ref) {
				print OUT " $_";
			}		
		# Print line break
		print OUT "\n";
		}
	}
	# Close fasta file handle
	close OUT;	
	# Screen output
	print "Results printed to file '$out_file'.\n";	
}
sub DistAbsPosFromMidRegion {
### Accepts: Position (pos) and start and end points of a region (start, end)
### Returns: Absolute distance between position and region midpoint 
	abs( $_[0] - ( ($_[2] + $_[3]) / 2 ) );
}
###