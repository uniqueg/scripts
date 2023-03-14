#!/usr/bin/perl

#==================#
#   HEADER START   #
#==================#
### Name: find_mates.pl
### Created: Sep 6, 2013
### Modified: Sep 6, 2013
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
### Version: v1.0
### Adapted from: n/a
### Requirements: n/a
#==================#
### Description: From two BED files containing first mate and second mate reads/regions (following the naming format: read_id\1 and read_id\2), extract all read pairs that are on the same chromosome and opposite strands and whose merged region is in a specified size window.
### Output: BED file with valid mates, merged into single regions ("\1" and "\2" dropped from IDs)
### Usage: perl ./find_mates.pl for information on required and optional arguments
#==================#
#    HEADER END    #
#==================#


#==========================#
#   PRE-REQUISITES START   #
#==========================#
#---> PRAGMAS / PACKAGES <---#
use strict;
use warnings;
use Getopt::Long;

#---> USAGE <---#
my $usage_info = &usage;

#---> OPTIONS / ARGUMENTS <---#
my $usage = '';
my $quiet = '';
my $mate_one = '';
my $mate_two = '';
my $out_file = '';
my $min_dist = 50;
my $max_dist = 5000;
my $head = 0;
my $options_result = GetOptions (
	'usage|help' => \$usage,
	'quiet' => \$quiet,
	'min_distance=i' => \$min_dist,
	'max_distance=i' => \$max_dist,
	#-----------------------#
	'mate_one=s' => \$mate_one,
	'mate_two=s' => \$mate_two,
	'out_file=s' => \$out_file
);
die $usage_info if $usage || !$options_result;
die $usage_info if !$mate_one || !$mate_two || !$out_file; 

#==========================#
#    PRE-REQUISITES END    #
#==========================#


#================#
#   MAIN START   #
#================#
#---> STATUS MESSAGE <---#
print "Starting '$0'...\n" unless $quiet;

#---> MAIN VARIABLES <---#
my $id_two_hoa_ref;

#---> BODY <---#
$id_two_hoa_ref = &tab_to_hoa_of_sorted_arrays($mate_two, 4);
&strip_mate_id($id_two_hoa_ref);
&validate_merge_and_print_mates_to_BED($mate_one, $id_two_hoa_ref, $out_file, $min_dist, $max_dist);

#---> STATUS MESSAGE <---#
print "Done.\n" unless $quiet;

#---> PROGRAM EXIT <---#
exit 0;
#================#
#    MAIN END    #
#================#


#=======================#
#   SUBROUTINES START   #
#=======================#
sub usage {
### Function: Returns usage information for current script
### Accepts: n/a
### Returns: String with usage information
### Type: Specialized
'Usage: perl ./find_mates.pl [OPTIONS] --read_file [FILE|BED|SAM] --id_file [FILE|TXT] --out_file [FILE|BED|SAM]

Description: From two BED files containing first mate and second mate reads/regions (following the naming format: read_id\1 and read_id\2), extract all read pairs that are on the same chromosome and opposite strands and whose merged region is in a specified size window.

==================================================
Required arguments:
--mate_one	Read file in BED (default) or SAM (with --sam switch) format
--mate_two	Flat text file with read IDs (one per line)
--out_file	Output file in same format as input file
==================================================
Optional arguments:
--min_distance [INT]	Maximum required distance between (start of) first and (end of) second mate (default: 50)
--max_distance [INT]	Maximum allowed distance between (start of) first and (end of) second mate (default: 5000)
--usage|help	Show this information and die
--quiet	Shut up!
';
}
#-----------------------#
sub tab_to_hoa_of_sorted_arrays {
### Function: Reads a TAB file into a hash of arrays of sorted arrays, grouped by the specified column
### Accepts: 1. TAB file; grouping column (1-based) [INTEGER]; 3. (OPTIONAL) one or more sets of sorting information (multiples of 3 in the order of priority!): A. Reference column for sorting [INTEGER]; B. Sorting type [STRING | ALLOWED: "num", "alpha"]; C. Sorting order [STRING | ALLOWED: "asc", "desc"]; DEFAULT: NO SORTING!
### Returns: 2. Reference to hash of arrays of sorted arrays (hash keys = chromosomes; elements outer arrays = TAB file rows; elements inner arrays = values of single TAB file row)
### Dependencies: n/a
### Type: Generic
	#---> PASS ARGUMENTS ---#
	my ($in_file, $group_col, @sort) = @_;	
	
	#---> STATUS MESSAGE <---#
	print "Reading second mates from file '$in_file' into hash..." . "\n" unless $quiet;

	#---> SUBROUTINE VARIABLES <---#
	my %HoAoA;	
	
	#---> BODY <---#
	# Open file handle TAB
	open TAB, $in_file;
	#-----------------------#
	## Crawl through TAB file line by line
	while (<TAB>) {
		# Remove trailing new line characters
		chomp;
		# Split line into array @line by tabs
		my @line = split /\t/, $_;
		# Extract grouping column value from @line array
		my $group = $line[$group_col - 1];
		# Push line values array @line to respective hash value / outer array
		push @{$HoAoA{$group}}, [ @line ];
	}		
	#-----------------------#
	# Close TAB file handle
	close TAB;
	#-----------------------#
	## Sort inner arrays of %HoAoA
	# Reverse sorting information array @sort, so that last entries have least priority 
	@sort = reverse @sort;
	# As long as there are entries in sorting information array @sort
	while (@sort) {
		## Shift information for one sorting round from array (sorting with least priority is executed first!)
		my $order = shift(@sort);
		my $type = shift(@sort);
		my $col = shift(@sort) - 1;
		## Check for different combinations of sorting type and order
		## Numerically in ascending order
		if ( $type eq "num" && $order eq "asc" ) {
			# Crawl through each inner array individually
			foreach my $array_ref (values %HoAoA) {
				# Sort...
				@{$array_ref} = sort { $a->[$col] <=> $b->[$col] } @{$array_ref};
			}
		}
		# Alphanumerically in ascending order
		elsif ( $type eq "alpha" && $order eq "asc" ) {
			# Crawl through each inner array individually
			foreach my $array_ref (values %HoAoA) {
				# Sort...
				@{$array_ref} = sort { $a->[$col] cmp $b->[$col] } @{$array_ref};
			}				
		}
		# Numerically in descending order
		elsif ( $type eq "num" && $order eq "desc" ) {
			# Crawl through each inner array individually
			foreach my $array_ref (values %HoAoA) {
				# Sort...
				@{$array_ref} = sort { $b->[$col] <=> $a->[$col] } @{$array_ref};
			}
		}
		# Alphanumerically in descending order
		elsif ( $type eq "alpha" && $order eq "desc" ) {
			# Crawl through each inner array individually
			foreach my $array_ref (values %HoAoA) {
				# Sort...
				@{$array_ref} = sort { $b->[$col] cmp $a->[$col] } @{$array_ref};
			}				
		}
		# Illegal sorting type and/or order 
		else {
			print "$order\t$col\t$type\n";
			print STDERR "Invalid value for sorting type (only 'num' and 'alpha' allowed!) or order (only 'asc' and 'desc' allowed!)" . "\n";
		}
	}
		
	#---> STATUS MESSAGE <---#
	print "Second mates read." . "\n" unless $quiet;
	
	#---> RETURN VALUE <---#
	return \%HoAoA;
}
#-----------------------#
sub strip_mate_id {
### Function: Strips off mate number in read IDs (in place substitution!)
### Accepts: Reference to hash of read IDs without terminal mate IDs (of the form "/1" and "/2")
### Returns: n/a
### Dependencies: n/a
### Type: Generic
	#---> PASS ARGUMENTS ---#
	my $read_id_hash_ref = shift;
	
	#---> STATUS MESSAGE <---#
	print STDOUT "Stripping mate identifiers." . "\n" unless $quiet;

	#---> BODY <---#
		
		## Traverse through each key $id of hash %$read_id_hash_ref
		foreach my $id ( keys %$read_id_hash_ref ) {
		    # Save hash value
		   	my $value = $read_id_hash_ref->{$id} if defined $read_id_hash_ref->{$id};
		    # Delete hash element
		    delete $read_id_hash_ref->{$id};
		    ## In hash key, strip off "/1" and "/2" and re-add the hash key to the hash
			if ($id =~ /\/1$/) {
				($id = $id) =~ s/\/1$//;
				$read_id_hash_ref->{$id} = defined ($value) ? $value : undef;
			}
			else {
				($id = $id) =~ s/\/2$//;
				$read_id_hash_ref->{$id} = defined ($value) ? $value : undef;
			}
		}
	
	#---> STATUS MESSAGE <---#
	print STDOUT "Mate identifiers stripped." . "\n" unless $quiet;
}
#-----------------------#
sub validate_merge_and_print_mates_to_BED {
### Function: function
### Accepts: 1. Input file in SAM or BED format; 2. Reference to hash containing read IDs as keys (empty strings as values); 3. Output file; 4. Flag whether header should be printed (only valid for SAM files)
### Returns: n/a
### Dependencies: n/a
### Type: Generic
	#---> PASS ARGUMENTS ---#
	my ($mate_one, $mate_two_hoa_ref, $out_file, $min_dist, $max_dist) = @_;
	
	#---> STATUS MESSAGE <---#
	print STDOUT "Crawling through first mates in file '$mate_one'..." . "\n" unless $quiet;

	#---> BODY <---#
	
		## Open file handles
		open FIRST_MATE, "<", $mate_one or die "[ERROR] Could not open file '$mate_one'!\n";
		open OUT, ">", $out_file or die "[ERROR] Could not open file '$out_file'!\n";
		
		#---> Process firt mate regions and compare to second mates in hash of arrays reference <---#
		while (<FIRST_MATE>) {
			# Split line into field values for first mate
			my ($refseq_1, $start_1, $stop_1, $name_1, $score_1, $strand_1) = split "\t";
			# Remove mate identifier from name
			($name_1 = $name_1) =~ s/\/\d$//; 	
			# Skip if second mate cannot be found
			if (exists $mate_two_hoa_ref->{$name_1}) {
				# Split referenced array into field values for second array
		
		###################### CHECK THIS FROM SCRIPT ON BC2 ##########################
		
				my ($refseq_2, $start_2, $stop_2, $name_2, $score_2, $strand_2) = @{$mate_two_hoa_ref->{$name_1}};
				# Skip entry if sequences are on different reference sequences
				next unless $refseq_1 eq $refseq_2;
				# Skip entry if sequences are NOT on different strands
				next unless $strand_1 ne $strand_2;
				# If first mate is on Watson strand...
				if ($strand_1 eq "+") {
					# Skip entry if distance is negative or too big
					next if ( ( $stop_2 - $start_1 ) > $max_dist ) || ( ( $stop_2 - $start_1 ) < $min_dist );
					# Print
					print OUT $refseq_1 . "\t" . $start_1 . "\t" . $stop_2 . "\t" . $name_1 . "\t" . $score_1 . ";" . $score_2 . "\t" . $strand_1 . "\n";
				}
				# If first mate is on Crick strand...
				if ($strand_1 eq "-") {
					# Skip entry if distance is negative or too big
					next if ( ( $stop_1 - $start_2 ) > $max_dist ) || ( ( $stop_1 - $start_2 ) < $min_dist );
					# Print
					print OUT $refseq_1 . "\t" . $start_2 . "\t" . $stop_1 . "\t" . $name_1 . "\t" . $score_1 . ";" . $score_2 . "\t" . $strand_1 . "\n";
				}
			}
		}
	
		## Close file handles
		close OUT;
		close FIRST_MATE;
	
	#---> STATUS MESSAGE <---#
	print STDOUT "Mates merged. Output written to '$out_file'." . "\n" unless $quiet;
}
#=======================#
#    SUBROUTINES END    #
#=======================#