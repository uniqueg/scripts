#!/usr/bin/perl

#######
### AUTHOR INFO & CREATION DATE:
### ----------------------------
### Author:		Alexander Kanitz
### Created: 	11-JAN-2013
### Modified:	21-JAN-2013
#######

#######
### FUNCTION:
### ---------
### Calculates the relative distance for each region of a query BED file to the closest region in a target BED file
#######

#######
### ARGUMENTS:
### ----------
### 1. Query BED file name
### 2. Target BED file name
### 3. Genome sizes file name
### 4. Output file name
### General Notes: Compare the examples in USAGE for further information.
#######

#######
### OUTPUT:
### -------
### A tab-separated file containing the IDs (4th column!) of the query BED file and the closest distance to any of the regions in the target BED file; negative values indicate that the query is upstream of the target
#######

#######
### USAGE:
### ------
### perl /path/to/distance_bed_bed.pl <path/to/query_bed_file> <path/to/target_bed_file> <path/to/genome_sizes_fasta_file> <path/to/output_file>
### Example: perl ./distance_bed_bed.pl ./query.bed ./target.bed ./hg19_sizes.fa ./out
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
use warnings;
use strict;
## Command-line arguments / initialization
my $query_file = shift;
my $target_file = shift;
my $gen_file = shift;
my $out_file = shift;
###


### B. Main program 
# Load query BED file to hash of arrays
my $query_HoHoA_ref = &tabToHashOfHashesOfArrays($query_file, 4, 1, "start", "stop", "rest");
# Load target BED file to hash of arrays
my $target_HoA_ref = &bedMidRegionToHashOfSortedArrays($target_file);
# Create chromosome sizes lookup hash
my $chr_size_hash_ref = &fastaToHash($gen_file);
# Calculate distances
&bedDistance($query_HoHoA_ref, $target_HoA_ref, $chr_size_hash_ref);
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

sub tabToHashOfHashesOfArrays {
### Accepts file in bed format (or other tab-separated file)
### Returns a hash of hashes of arrays containing all column values (addressable by column name) grouped by a common column
### The total number of columns, the number of the grouping column (= key for the outer hash), and the names of the remaining columns (= keys for the inner hashes) are passed when calling (e.g. 3, 1, "start", "stop", "strand")
### Values of the non-grouping columns are used to populate the values of the inner hashes (i.e. arrays)
	## Pass arguments
	my $bed_file = shift;
	my $no = shift;
	my $first = shift;
	## Pass hash key arguments
	my @keys;
	foreach (@_) {
		push @keys, $_;
	}
	# Declare/initialize variables
	my %HoHoA;
	# Open bed file handle
	open BED, $bed_file;
	## Crawl through bed file line by line
	while (<BED>) {
		# Remove trailing new line characters
		chomp;
		# Split line into array @line by tabs
		my @line = split /\t/, $_, $no;
		# Extract lookup/key value for outer hash from @line array
		my $value = splice @line, ($first - 1), 1;
		## IF outer hash entry with key $value does not exist, create both outer and inner hash
		unless (exists $HoHoA{$value}) {
			## Add value for each key of inner hash
			foreach (@keys) {
    			$HoHoA{$value}{$_} = [ shift @line ];
			}
		}
		## ELSE push to values (i.e. arrays) of inner hash
		else {
			## For each key of inner hash, push values to existing ones (array!)
			foreach (@keys) {
    			push @{ $HoHoA{$value}{$_} } , shift @line;
			}
		}
	}
	# Close bed file handle
	close BED;
	# Print status message
	print "File $bed_file processed.\n";
	# Return bed file hash of hashes of arrays
	return \%HoHoA;
}

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

sub bedDistance {
###
###
###
	## Pass arguments
	my %query = %{shift()};
	my %target = %{shift()};
	my %chr_size = %{shift()};
	# Print status message
	print "Calculating distances...\n";
	# Open output file handle 
	open OUT, ">$out_file";
	## Crawl through chromosomes (= outer hash of HoHoA) one by one
	foreach my $chr (keys %query) {
		## Repeat the following as long as there are values in HoHoA array "start"
		while (@{$query{$chr}{"start"}}) {
			## Declare/initialize distance variables ($min_distance gets chromosome size as initial value)
			my $dist;
			my $min_dist = $chr_size{$chr};
			## Populate query-related variables by shifting from HoHoA arrays
			my $start_query = shift @{$query{$chr}{"start"}};
			my $stop_query = shift @{$query{$chr}{"stop"}};
			my $mid_query = ($start_query + $stop_query) / 2;
			my $id_query = shift @{$query{$chr}{"rest"}};
			## Declare target midpoint variable
			my $mid_target;
			## Test whether current chromosome (= outer hash of HoHoA) exists in target HoHoA; else skip chromosome and print nothing!
			if (exists $target{$chr}) {
				## Crawl through values in HoA arrays (i.e. regions in target bed file)
				for (my $i = 0; $i < @{$target{$chr}}; $i++) {	
					# Initialize target midpoint by shifting from HoA array
					$mid_target = ${$target{$chr}}[$i];
					# Calculate distance between mid-point of current query and target ranges
					$dist = $mid_query - $mid_target;
					# Proceed to next query range if absolute distance increases compared to minimum distance (i.e. last distance)
					last if abs $dist >= abs $min_dist;
					# Set minimum distance to current distance
					$min_dist = $dist;
				}
				# Print ID of query and minimum distance to output file
				print OUT "$id_query\t$min_dist\n";
			}
		}
	}
	# Close output file handle
	close OUT;
	# Screen output
	print "Minimum distances written to file: $out_file\n";
}