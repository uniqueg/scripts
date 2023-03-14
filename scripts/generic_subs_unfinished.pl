sub tabMergeDuplicatesByColumnToSortedAoARef {
### Accepts: 1. name of tab-delimited output file to process; 2. number of header lines (default: 0); 3. column index of values to merge (default: 1)
### Returns: Reference to array (one line/row per merged value) of arrays (individual elements of that line/row); duplicate entries are merged according to the 'mergeRules' subroutine; the returned array of arrays is sorted according to the merged values in alphabetic; the header lines (if present) are  
### Type: Generic
### Dependencies: mergeRules subroutine (define rules for merging duplicate entries; leave blank to always keep the first entry and remove all duplicates)
	## Pass parameters
	my $in_file = shift;
	my $header = @_ ? shift : 0;
	my $by = @_ ? shift : 1;
	# Decrease $by by 1 to match array indices
	$by--;
	# Declare temporary hash of arrays and results array
	my %HoA;
	my @AoA;
	# Open input file handle
	open IN, $in_file;
	## While there are header lines...
	while ($header) {
		# Decrease header count by 1
		$header--;
		# Read line to dedicated variable
		my $line = <IN>;
		# Remove trailing newline character
		chomp($line);
		# Split header line by tab
		my @header = split /\t/, $line;
		push @AoA, [@header];
	}
	## Traverse through each line of input file
	while (<IN>) {
		# Remove trailing newline character
		chomp;
		# Split line by tabs
		my @line = split /\t/;
		## IF hash %HoA entry with same name (= key) exists...
		if (exists $HoA{$line[$by]}) {
			# Call subrouting defining merging rules
			$HoA{$line[$by]} = &mergeRules($HoA{$line[$by]}, \@line);
		} 
		## ELSE create new hash %HoA entry
		else {
			# Initialize empty array as %HoA value for current name (= key)
			$HoA{$line[$by]} = [];
			# Push line into 
			push $HoA{$line[$by]}, @line;
		}
	}
	# Close input file handle
	close IN;
	## Sort hash values alphanumerically and...
	foreach my $key (sort { $a cmp $b } keys %HoA) {
		# push to @results array
		push @AoA, $HoA{$key};
	} 
	# Screen output
	print "Processed file '$in_file'.\n";
	# Return @AoA reference
	return \@AoA;
}
#-----------------------#
sub mergeRules {
### Accepts: 1. reference to array containing current/reference value(s); 2. reference to array containing new value(s) 
### Returns: Reference to array with updated values determined by the merging rules
### Type: Specific to merging duplicates from 'clipz_overlap_cl_sites.pl' output
### Dependencies: Called by the generic subroutine 'tabMergeDuplicatesByColumnToSortedAoARef' to apply specific rules governing the merging of duplicate values
	## Pass arguments
	my $ref_array_ref = shift;
	my $new_array_ref = shift;
	## Do comparisons between existing and new entry
	## GENERAL MERGING RULES:
	# Increase total number of crosslinks
	$$ref_array_ref[10] += $$new_array_ref[10];
	# Calculate current total score
	$$ref_array_ref[9] += $$new_array_ref[9];
	## CONDITION-DEPENDENT MERGING RULES:
	## IF rank of new line is lower (i.e. better) than current/previous rank...
	if ($$new_array_ref[7] < $$ref_array_ref[7]) {
		# ...set new start position of best crosslink
		$$ref_array_ref[5] = $$new_array_ref[5];
		# ...set new end position of best crosslink
		$$ref_array_ref[6] = $$new_array_ref[6];
		# ...set new rank of best crosslink
		$$ref_array_ref[7] = $$new_array_ref[7];
		# ...set new score of best crosslink
		$$ref_array_ref[8] = $$new_array_ref[8];
	}
	## IF start coordinate of current exon is lower than current/previous one...
	if ($$new_array_ref[2] < $$ref_array_ref[2]) {
		# ...set to newer, lower value
		$$ref_array_ref[2] = $$new_array_ref[2];
	}
	## IF end coordinate of current exon is higher than current/previous one...
	if ($$new_array_ref[3] > $$ref_array_ref[3]) {
		# ...set to newer, lower value
		$$ref_array_ref[3] = $$new_array_ref[3];
	}
	# Return array reference
	return $ref_array_ref;	
}
#-----------------------#
sub rearrangeTab {
### Accepts: 1. tab-delimited input file; 2. reference to array holding new order; 3. output file name; 4; default value for unavailable values (default: "N/A")
### Returns: 2. tab-delimited output file (rearranged)
### Type: Generic
	## Pass arguments
	my $in_file = shift;
	my $new_order_array_ref = shift;
	my $out_file = shift;
	my $default_NA = @_ ? shift : "N/A";
	# Open input file handle
	open IN, $in_file;
	# Open output file handle
	open OUT, ">$out_file";
	## Traverse through input file line by line
	while (<IN>) {
		# Remove trailing newline character
		chomp;
		# Split line by tab;
		my @in_line = split /\t/;
		# Declare variable for output line array
		my @out_line;
		## Traverse through each element in the array referenced by $new_order_array_ref
		for my $index (0..$#$new_order_array_ref) {
			## IF a non-numerical character or an entry surrounded by '/' for the value of the current index of the new order array, use as string (i.e. write found value, minus surrounding slashes if available, to @out_line element)
			if ($$new_order_array_ref[$index] =~ m:^0-9]: | $$new_order_array_ref[$index] =~ m:^/.*/$:) {
				# Assign value of the current index of the new order array (i.e. a string) to a helper variable
				my $string = $$new_order_array_ref[$index];
				# Remove surrounding slashes '/' if present
				$string =~ s:^/(.*)/$:$1:;
				# Assign string to value of the current index of the @out_line array
				$out_line[$index] = $string;
			}
			## ELSE... 
			else {
				## ...check IF value in column specified by new order array reference is present (mind Perl array index numbering!)
				if (defined $in_line[$$new_order_array_ref[$index] - 1]) {
					## If so, assign to @out_line array element with the indicated column index (mind Perl array index numbering!)
					$out_line[$index] = $in_line[$$new_order_array_ref[$index] - 1];
				}
				## ELSE print warning message and assign default value for missing values
				else {
					# Get column number (mind Perl array index numbering!)
					my $column = $index + 1;
					# Print warning message to STDOUT
					print "Warning: No value in row $., column $column, writing '$default_NA' instead.\n";
					# Assign default value for missing values to @out_line array element with the indicated column index
					$out_line[$index] = $default_NA;
				}
			}
		}
		# Join output line array elements by tab
		print OUT join("\t", @out_line) . "\n";
	}
	# Close output file handle
	close OUT;
	# Close input file handle
	close IN;
	## Status messages
	print "Processed file '$in_file'.\n";
	print "Output written to file '$out_file'.\n";
}
#-----------------------#
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
#-----------------------#
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
#-----------------------#
sub bedContinuousPositions {
### Writes continuous genome positions for each region to the indicated output file
### Requires a reference to a hash of arrays containing mid-points of regions sorted by chromosome, a reference to a hash containing continuous genomic starting positions for each chromosome and an output file name 
	## Pass arguments
	my %bed = %{shift()};
	my %cont = %{shift()};
	my $out = shift;
	# Print status message
	print "Calculating positions...\n";
	# Open output file handle 
	open OUT, ">$out";
	## Crawl through chromosomes (= outer hash of HoA) one by one
	foreach my $chr (keys %bed) {
		## Crawl through values in HoA arrays (i.e. regions in target bed file)
		foreach my $region (@{$bed{$chr}}) {	
			# Calculate continuous genomic position for current region
			$region += $cont{$chr} - 1;
			# Print to output file
			print OUT "$region\n";
		}
	}
	# Close output file handle
	close OUT;
	# Screen output
	print "Positions written to file: '$out'\n";
}
#-----------------------#

#-----------------------#

#-----------------------#

#-----------------------#

#-----------------------#

#-----------------------#

#-----------------------#



#===================#
#   SANDBOX START   #
#===================#

sub fasta_to_hash {
### Function: Read FASTA file into hash; multiple sequence lines per entry, if present, are concatenated
### Accepts: 1. FASTA filename
### Returns: Reference to sequence hash
### Dependencies: n/a
### Type: Generic
	#---> PASS ARGUMENTS ---#
	my $in_file = shift;	
	
	#---> SUBROUTINE VARIABLES <---#
	my (%seq_hash, $id, $last_id_line);
	my $seq = "";
	
	#---> BODY <---#
	# Open file handle FA
	open FA, "$in_file";	

	## Traverse through file handle FA line by line
	while ($line = <FA>) {
		# Remove trailing newline character "\n"
		chomp $line;
		## Check whether current line is identifier line
		if ( $line =~ s/^>// ) {
			## If so, check whether hash key with current identifier $id already exists
			if ( defined $hash{$id} ) {
				## If so, check whether $id/$seq pair is a duplicate (i.e. hash value of current identifier $id matches current $seq)
				if ( $hash{$id} eq $seq ) {
					# If so, print warning message to STDERR...
					print STDERR "Sequence entry starting in line $last_id_line is a duplicate. Only one instance kept."
				}
				## Else, print warning message to STDERR...
				else {
					print STDERR "Identifier in line $last_id_line of '$in_file' not unique. Adding suffix.";
					# Try to add suffix by initializing a serial number $suf_nr... 
					my $suf_nr = 0;
					# ...and increasing it as long as a hash key made of the current serial number $suf_nr appended to $id via an underscore linker exists
					do $suf_nr++ while defined $hash{$id . "_" . $suf_nr};
					# Append latest serial number $suf_nr to $id via an underscore linker
					$id += "_" . $suf_nr;
				}
			}
			# Set $id to current line, $last_id_line to current line number in $in_file and $seq to empty string ""	
			($id, $last_id_line, $seq) = ($., $line, $seq);
			# Assign sequence $seq to hash key $id unless $seq is an empty string (very first instance!) 
			$hash{$id} = $seq unless $seq eq "";
		}
		## Else...
		else {
		    # Append current line to existing sequence string $seq
		    $seq = $seq.$line;
		}
	}
	# Add last entry to hash
	$hash{$id} = $seq;
	
	# Close file handle FA
	close FA;
		
	#---> STATUS MESSAGE <---#
	print "FASTA File $in_file read.\n";	
	
	#---> RETURN VALUE <---#
	return \%seq_hash;
}