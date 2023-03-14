#-----------------------#
sub bed_to_hoa_of_sorted_arrays {
### Function: Reads a BED file into a hash of arrays of sorted arrays, grouped by chromosome; numerical sorting of inner arrays by start then stop position in ascending order
### Accepts: 1. BED file
### Returns: 2. Reference to hash of arrays of sorted arrays (hash keys = chromosomes; elements outer arrays = BED file rows; elements inner arrays = values of single BED file row)
### Dependencies: n/a
### Type: Generic
	#---> PASS ARGUMENTS ---#
	my ($in_file, @sort) = @_;	
	
	#---> SUBROUTINE VARIABLES <---#
	my %HoAoA;	
	
	#---> BODY <---#
	# Open file handle BED
	open BED, $in_file;
	#-----------------------#
	## Crawl through BED file line by line
	while (<BED>) {
		# Remove trailing new line characters
		chomp;
		# Split line into array @line by tabs
		my @line = split /\t/, $_;
		# Extract chromosome (= key for outer hash) from @line array
		my $chr = $line[0];
		# Push line values array @line to respective hash value / outer array
		push @{$HoAoA{$chr}}, [ @line ];
	}		
	#-----------------------#
	# Close BED file handle
	close BED;
	#-----------------------#
	## Sort inner arrays of %HoAoA according to stop position
	# Crawl through each inner array individually
	foreach my $array_ref (values %HoAoA) {
		# Sort...
		@{$array_ref} = sort { $a->[2] <=> $b->[2] } @{$array_ref};
	}
	#-----------------------#
	## Sort inner arrays of %HoAoA according to start position
	# Crawl through each inner array individually
	foreach my $array_ref (values %HoAoA) {
		# Sort...
		@{$array_ref} = sort { $a->[1] <=> $b->[1] } @{$array_ref};
	}
		
	#---> STATUS MESSAGE <---#
	print "Read BED file '$in_file'" . "\n" unless $quiet;
	
	#---> RETURN VALUE <---#
	return \%HoAoA;
}
#-----------------------#
sub cat_adapter_to_seq_array {
### Function: Concatenate an adapter sequence to the 5'- or 3'-end of a set of sequences stored in an array
### Accepts: 1. Reference to sequence array; 2. Adapter sequence [STRING]; 3. Orientation [INT; ALLOWED = 3, 5 | DEFAULT = 5]
### Returns: Reference to array of sequences
### Dependencies: n/a
### Type: Generic
	#--- PASS ARGUMENTS ---#
	my ($seq_array_ref, $adapter, $orientation) = (@_, "");
	## Check whether orientation argument is illegal
	if ($orientation ne "3" && $orientation ne "5" ) {
		## If so, print status message and use default value
		print "Illegal value for orientation of adapter concatentation. Default value (5'-end addition) is used.\n" unless $quiet;
		$orientation = 5;		
	} 

	#--- SUBROUTINE VARIABLES ---#
	my @seqs;
	
	#--- BODY ---#
	## Traverse through each element $seq of array @$seq_array_ref
	foreach my $seq ( @$seq_array_ref ) {
	    # Push concatenated sequence to array (order depending on value in $orientation) 
		push @seqs, ($orientation == 5) ? $adapter.$seq : $seq.$adapter;
	}	
	
	#---> STATUS MESSAGE <---#
	print "Added sequence '$adapter' to the $orientation'-end of each sequence.\n" unless $quiet;	
	
	#--->  RETURN VALUE  <---#
	return \@seqs;	
}
#-----------------------#
sub check_key_value_pair {
### Function: Tests whether a specified key already exists in a hash; if so, and if the values match (i.e. entry is a duplicate), prints warning message to STDERR, but only one instance is kept; if so, and if the values do not match, prints warning message to STDERR and appends (via an optional specified linker) to key
### Accepts: 1. Hash key to test; 2. Hash value to test; 3. Hash reference; 4. Linker [STRING | DEFAULT = "_"]
### Returns: 1. Error code (0: original key unique ('normal' case); 1: key/value pair is a duplicate; 2: original key not unique, linker and serial number appended); 2. Unique key
### Dependencies: n/a
### Type: Generic
	#---> PASS ARGUMENTS ---#
	my ($key, $value, $hash_ref, $linker) = (@_, "_");
	
	#---> BODY <---#
	## Check whether key with current identifier $key exists in hash %$hash_ref
	if ( defined $$hash_ref{$key} ) {
		## If so, check whether $key/$value pair is a duplicate (i.e. hash value of current identifier $key matches current $value)
		if ( $$hash_ref{$key} eq $value ) {
			# If so, return error code 1 and original key $key
			return (1, $key);
		}
		## Else...
		else {
			# Try to add suffix by initializing a serial number $suf_nr... 
			my $suf_nr = 1;
			# ...and increasing it as long as a hash key made of the current serial number $suf_nr appended to $key via linker $linker exists
			$suf_nr++ while defined $$hash_ref{$key . "_" . $suf_nr};
			# Append latest serial number $suf_nr to $key via linker $linker
			$key .= $linker . $suf_nr;
			# Return error code 2 and new, unique key $key
			return (2, $key);
		}
	}
	# Else...
	else {
		# Return error code 0 and original key $key
		return (0, $key);	
	}	
}
#-----------------------#
sub dec_to_32d_bin {
### Function: Convert integer to 32-digit binary number
### Accepts: 1. Integer
### Returns: 32 digit binary representation of specified integer
### Dependencies: n/a
### Type: Generic
	#---> PASS ARGUMENTS ---#
	my $int = shift;
	
	#---> BODY <---#
	my $bin = unpack("B32", pack("N", $int));	
	
	#---> RETURN VALUE <---#
	return $bin;
}
#-----------------------#
sub fasta_to_hash {
### Function: Read FASTA file into hash; multiple sequence lines per entry, if present, are concatenated; no alphabet checking is performed
### Accepts: 1. FASTA filename
### Returns: Reference to sequence hash
### Dependencies: Subroutine 'check_key_value_pair'
### Type: Generic
	#---> PASS ARGUMENTS ---#
	my $in_file = shift;	
	
	#---> SUBROUTINE VARIABLES <---#
	my (%seq_hash, $id, $last_id_line);
	my $seq = "";
	
	#---> BODY <---#
	# Open file handle FA
	open FA, "$in_file";	
	#-----------------------#
	## Traverse through file handle FA line by line
	while (my $line = <FA>) {
		# Remove trailing newline character "\n"
		chomp $line;
		## Check whether current line is identifier line
		if ( $line =~ s/^>// ) {
			## Check whether identifier $id is defined and sequence $seq is not empty (the very first case!); if so... 
			if (defined $id && $seq ne "") {
				# Test whether key or key/value sequence pair is duplicate and generate unique $id (saved in $$result_array_ref[1]) if necessary
				my ($err_code, $new_id) = &check_key_value_pair($id, $seq, \%seq_hash);
				# Write warning message to STDERR if error code in $err_code is 1
				print STDERR "Sequence entry starting in line $last_id_line is a duplicate. Only one instance kept." . "\n" if $err_code == 1;
				# Write warning message to STDERR if error code in $err_code is 2				
				print STDERR "Identifier in line $last_id_line of '$in_file' not unique. Adding suffix." . "\n" if $err_code == 2;
				# Assign sequence $seq to hash key $new_id
				$seq_hash{$new_id} = $seq;				
			}
			# Set $id to current line, $last_id_line to current line number in $in_file and $seq to empty string ""	
			($id, $last_id_line, $seq) = ($line, $., "");
		}
		## Else...
		else {
		    # Append current line to existing sequence string $seq
		    $seq = $seq.$line;
		}
	}
	#-----------------------#
	## Add last entry to hash
	# Test whether key or key/value sequence pair is duplicate and generate unique $id (saved in $$result_array_ref[1]) if necessary
	my ($err_code, $new_id) = &check_key_value_pair($id, $seq, \%seq_hash);
	# Write warning message to STDERR if error code in $err_code is 1
	print STDERR "Sequence entry starting in line $last_id_line is a duplicate. Only one instance kept." . "\n" if $err_code == 1;
	# Write warning message to STDERR if error code in $err_code is 2				
	print STDERR "Identifier in line $last_id_line of '$in_file' not unique. Adding suffix." . "\n" if $err_code == 2;
	# Assign sequence $seq to hash key $new_id
	$seq_hash{$new_id} = $seq;		
	#-----------------------#
	# Close file handle FA
	close FA;
		
	#---> STATUS MESSAGE <---#
	print "FASTA File $in_file read.\n" unless $quiet;	
	
	#---> RETURN VALUE <---#
	return \%seq_hash;
}
#-----------------------#
sub files_regex {
### Function: Finds all files in a specified folder that contain a specified pattern in their filenames
### Accepts: 1. Path to folder; 2. Pattern [STRING]
### Returns: Reference to array of filenames
### Dependencies: n/a
### Type: Generic
	#---> PASS ARGUMENTS ---#
	my ($path, $pattern) = @_;
	
	#---> SUBROUTINE VARIABLES <---#
	my @filenames;
	
	#---> BODY <---#
	# Open directory handle DIR or die
	opendir (DIR, $path) or die $!;
	#-----------------------#
	## Traverse through files in DIR
	while (my $file = readdir(DIR)) {
 		# Push file to array @filenames if filename contains $pattern 
		push @filenames, $file if $file =~ m/$pattern/;
	}
	#-----------------------#
	# Close directory handle DIR
	closedir (DIR);
	
	#---> STATUS MESSAGE <---#
	print "Directory '$path' processed. " . scalar @filenames . " file(s) matched the indicated pattern '$pattern'." . "\n" unless $quiet;
	
	#---> RETURN VALUE <---#
	return \@filenames;
}
#-----------------------#
sub find_unique_n_mers_in_seq {
### Function: Extracts all possible, unique subsequences of a specified size from a a given sequence
### Accepts: 1. Query sequence [STRING]; 2. N-mer size [INT] 
### Returns: Reference to array of subsequences
### Dependencies: n/a
### Type: Generic
	#--- PASS ARGUMENTS ---#
	my ($seq, $n) = @_;

	#--- SUBROUTINE VARIABLES ---#
	my @seqs;

	#--- BODY ---#
	# Die if query sequence is shorter than N-mer size
	die "N-mer size in '--nt' bigger than length of sequence in '--query'!\nNo output can be generated.\nDied" if length($seq) < $n;
	## Generate subsequences by looping through query sequence nucleotide by nucleotide
	for ( my $start_pos = 0 ; $start_pos <= length($seq) - $n + 1; $start_pos++ ) {
	    # Extract subsequence and push to array
	    push @seqs, substr $seq, $start_pos, $n;
	}
	my $seqs_sorted_array_ref = &unique_array_elements(@seqs);

	#---> STATUS MESSAGE <---#
	print "Extracted all possible, unique subsequences from sequence '$seq'.\n" unless $quiet;

	#--->  RETURN VALUE  <---#
	return $seqs_sorted_array_ref;
}
#-----------------------#
sub name_array_elements_with_serial_number_and_add_to_hash {
### Function: Uniquely name elements of an array with serial numbers and an optional pre- or suffix; load to hash
### Accepts: 1. Reference to array; 2. Pre- or -suffix [STRING; DEFAULT = EMPTY STRING ""]; 3. Orientation [STRING; ALLOWED = "prefix", "suffix" | DEFAULT = "prefix"]
### Returns: Hash reference
### Dependencies: n/a
### Type: Generic
	#--- PASS ARGUMENTS ---#
	my ($array_ref, $fix, $orientation) = (@_, "", "");
	## Check whether orientation argument is illegal
	if ($orientation ne "prefix" && $orientation ne "suffix" ) {
		## If so, print status message and use default value
		print "Illegal value for orientation of concatentation. Default value ('prefix') is used.\n" unless $quiet;
		$orientation = "prefix";		
	} 

	#--- SUBROUTINE VARIABLES ---#
	my %hash;
	
	#--- BODY ---#
	# Determine necessary number of digits for unique numbering
	my $digits = length($#{$array_ref});
	## 
	for ( my $i = 1; $i <= $#{$array_ref}; ++$i ) {
	    # Generate unique name $id
	    my $id = ($orientation eq "prefix") ? $fix.sprintf("%0${digits}d", $i) : sprintf("%0${digits}d", $i).$fix;
		# Load hash (key: $id; value: array element)
		$hash{$id} = $$array_ref[$i-1];  
	}
	
	#---> STATUS MESSAGE <---#
	print "Added unique names.\n" unless $quiet;	
	
	#--->  RETURN VALUE  <---#
	return \%hash;	
}
#-----------------------#
sub aoa_to_file {
### Description: Writes the indicated, referenced array of arrays to the specified output file; the default settings produce a tab-separated file, but the intra- and interarray separators may be specified manually
### Accepts: 1. Reference to array of arrays; 2. Output file name; 3. String used to separate elements of inner arrays [DEFAULT = TAB "\t"]; 4. String used to separate outer array elements [DEFAULT = NEWLINE "\n"]; NOTE: If any one of the separators is set manually, always specify both!
### Returns: n/a 
### Dependencies: n/a
### Type: Generic
	#--- PASS ARGUMENTS ---#
	my ($AoA_ref, $out_file, $sep_intra, $sep_inter) = (@_, "\t", "\n");

	#--- BODY ---#
	# Open file handle OUT
	open OUT, ">$out_file";
	#-----------------------#
	## Loop through references to inner arrays (elements of outer array)
	foreach my $array_ref (@$AoA_ref) {
		# Print current inner array to output file; array elements are separated by $sep_intra (e.g. tab), different arrays by $sep_inter (e.g. newline)
		print OUT join($sep_intra,@$array_ref) . $sep_inter;
	}
	#-----------------------#
	# Close file handle OUT
	close OUT;

	#---> STATUS MESSAGE <---#
	print "Output written to file '$out_file'.\n" unless $quiet;	
}
#-----------------------#
sub hash_to_fasta_file {
### Function: Write sequence hash in FASTA format to indicated output fileUniquely name elements of an array with serial numbers and an optional pre- or suffix; load to hash
### Accepts: 1. Reference to sequence hash; 2. Output filename [STRING]
### Returns: n/a
### Dependencies: n/a
### Type: Generic
	#--- PASS ARGUMENTS ---#
	my ($hash_ref, $out_file) = @_;

	#--- BODY ---#
	# Open file handle OUT
	open OUT, ">$out_file";
	## Loop through sequence hash sorted by keys (sequence IDs)
	foreach my $seq (sort {$a cmp $b} keys %{$hash_ref}) {
		# Print in fasta format
		print OUT ">" . $seq . "\n" . ${$hash_ref}{$seq} . "\n";
	}		
	# Close file handle OUT
	close OUT;
	
	#---> STATUS MESSAGE <---#
	print "Sequences written to file '$out_file'.\n" unless $quiet;	
}
#-----------------------#
sub tab_to_aoa {
### Function: Reads a TAB file into an array of arrays
### Accepts: TAB file
### Returns: Reference to array of arrays (elements outer arrays = TAB file rows; elements inner arrays = values of single TAB file row)
### Dependencies: n/a
### Type: Generic
	#---> PASS ARGUMENTS ---#
	my $in_file = shift;	
	
	#---> SUBROUTINE VARIABLES <---#
	my @AoA;	
	
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
		# Push line values array @line to respective hash value / outer array
		push @AoA, [ @line ];
	}		
	#-----------------------#
	# Close TAB file handle
	close TAB;
	#-----------------------#
		
	#---> STATUS MESSAGE <---#
	print "Read TAB file '$in_file'" . "\n" unless $quiet;
	
	#---> RETURN VALUE <---#
	return \@AoA;
}
#-----------------------#
sub tab_to_hoa_of_sorted_arrays {
### Function: Reads a TAB file into a hash of arrays of sorted arrays, grouped by the specified column
### Accepts: 1. TAB file; grouping column [INTEGER]; 3. (OPTIONAL) one or more sets of sorting information (multiples of 3 in the order of priority!): A. Reference column for sorting [INTEGER]; B. Sorting type [STRING | ALLOWED: "num", "alpha"]; C. Sorting order [STRING | ALLOWED: "asc", "desc"]; DEFAULT: NO SORTING!
### Returns: 2. Reference to hash of arrays of sorted arrays (hash keys = chromosomes; elements outer arrays = TAB file rows; elements inner arrays = values of single TAB file row)
### Dependencies: n/a
### Type: Generic
	#---> PASS ARGUMENTS ---#
	my ($in_file, $group_col, @sort) = @_;	
	
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
	print "Read TAB file '$in_file'" . "\n" unless $quiet;
	
	#---> RETURN VALUE <---#
	return \%HoAoA;
}
#-----------------------#
sub unique_array_elements {
### Function: Removes duplicates of array; preserves order 
### Accepts: Array
### Returns: Array reference
### Dependencies: n/a
### Type: Generic
	#--- PASS ARGUMENTS ---#
	my @seqs = @_;

	#--- SUBROUTINE VARIABLES ---#
	my %uniq;
	
	#--- BODY ---#	
	my @sorted = grep !$uniq{$_}++, @_;

	#--->  RETURN VALUE  <---#
	return \@sorted;
}
#-----------------------#
sub one_file_per_line {
### Function: Reads file line by line and writes one file for each line; new line characters are removed; a serial number '.#' is appended to the original filename (with # from 1..number of rows in input file)
### Accepts: 1. Path to text file [STRING]
### Returns: 1. Reference to array containing the names of the generated files
### Dependencies: n/a
### Type: Generic
	## Pass arguments
	my $in_file = shift;
	# Declare array of output filenames
	my @out_names;
	# Initialize serial number for unique output filename generation
	my $sn = 1;
	# Open input file handle
	open IN, $in_file;
	## Traverse through input file line by line
	while (<IN>) {
		# Assign line content to dedicated variable
		my $line = $_;
		# Remove trailing newline character
		chomp $line;
		# Generate output filename
		$out_file = $in_file . ".$sn";
		# Add output filename to array of output filenames
		push @out_names, $out_file;
		# Open output file handle
		open OUT, ">$out_file";
		# Write line to output file handle
		print OUT $line;
		# Increase serial number count
		$sn++;
		# Close output file handle
		close OUT;	
	}
	# Close output file handle
	close OUT;
	# Return array reference
	return \@out_names;
}
#-----------------------#
sub total_length_bed {
### Function: Calculates the sum of the lengths of all regions in a BED file.
### Accepts: 1. Name of input BED file [FILE]; 2. Boolean (0 or 1) indicating whether the input BED file is in half/end-open format (default: 1)
### Returns: Sum of region lengths [INTEGER]
### Dependencies: n/a
### Type: Generic
	#---> PASS ARGUMENTS ---#
	my $in_file = shift;
	my $open = shift;

	#---> STATUS MESSAGE <---#
	print "Processing BED file '$in_file'...\n" unless $quiet;			

	#---> SUBROUTINE VARIABLES <---#
	my $sum = 0;
	
	#---> BODY <---#
	# Open file handle BED
	open BED, "$in_file";
	
	## Traverse through file handle BED line by line
	while (<BED>) {
		# Split line by tabs
		my ($chr, $start, $end) = split /\t/, $_;
		# Add current region length to sum of lengths
		$sum += $end - $start + 1 - $open;	
	}
		
	# Close file handle BED
	close BED;
	
	#---> STATUS MESSAGE <---#
	print "BED file '$in_file' processed!\nTotal length of all regions: $sum\n\n" unless $quiet;		
	
	#---> RETURN VALUE <---#
	return $sum;
}
#-----------------------#
sub are_hash_values_integers {
### Function: Verifies if all values of the hash passed in the arguments are integers. This is not very robust and will miss certain cases (scientific notation etc.).
### Accepts: 1. Hash reference
### Returns: '1' if all values of the specified hash are integers, else '0'
### Dependencies: n/a
### Type: Generic
	#---> PASS ARGUMENTS ---#
	my $hash_ref = shift;
	
	#---> STATUS MESSAGE <---#
	print STDERR "Testing whether hash contains integers only...\n" unless $quiet;	
	
	#---> BODY <---#
	## Traverse through each value $value of hash %$hash_ref
	foreach my $value ( values %$hash_ref ) {
    	## Use regular expression to test whether value is integer 
		unless (defined $value && $value =~ /^[+-]?\d+$/) {
			# Print status message
			print STDERR "Hash contains one or more non-integers!\n\n" unless $quiet;	
			# If not, return '0'
			return 0;
		}
	}
	
	#---> STATUS MESSAGE <---#
	print STDERR "Hash contains only integers!\n\n" unless $quiet;	
	
	#---> RETURN VALUE <---#
	return 1;
}
#-----------------------#
sub line_to_array {
### Function: Reads text file line by line and pushes each line (minus the newline character) into an array
### Accepts: 1. Path to text file [STRING]
### Returns: 1. Reference to array
### Dependencies: n/a
### Type: Generic

	#---> PASS ARGUMENTS ---#
	my $file = shift;

	#---> STATUS MESSAGE <---#
	print STDERR "Processing file '$file'..." . "\n" unless $quiet;

	#---> SUBROUTINE VARIABLES <---#
	my @lines;

	#---> BODY <---#

		#---> Open file <---#
		open FILE, "<", $file;

		#---> Push line to array <---#
		while (<FILE>) {
			chomp;
			push @lines, $_;
		}

		#---> Close file <---#
		close FILE;

	#---> STATUS MESSAGE <---#
	print STDERR "File '$file' processed." . "\n" unless $quiet;
	
	#---> RETURN VALUE <---#
	return \@lines;
	
}
#-----------------------#
sub bed_to_hoa {
### Function: Reads a BED file into a hash of arrays, with a specified field being the key (default: 4th/name column); a warning is issues if non-unique hash keys are found
### Accepts: 1. BED file; 2. field/column to use as hash key
### Returns: Reference to hash of arrays (hash keys = specified field; values = arrays containing the whole BED lines)
### Dependencies: n/a
### Type: Generic
	#---> PASS ARGUMENTS ---#
	my $in_file = shift;
	my $col = length scalar @_ ? shift : 4; 	

	#---> STATUS MESSAGE <---#
	print "Processing file '$in_file'...\n" unless $quiet;
	
	#---> SUBROUTINE VARIABLES <---#
	my %HoA;	
	
	#---> BODY <---#

		#---> Open file handle <---#
		open BED, $in_file;
	
		#---> Traverse through each line of file <---#
		while (<BED>) {

			#---> Remove trailing record separator <---#
			chomp;

			#---> Split line <---#
			my @line = split /\t/, $_;

			#---> Verify correct file format <---#
			die "[ERROR] File '$in_file' does not look like a valid input file in line $.\nExecution aborted.\n" unless defined $line[$col - 1];

			#---> Extract hash key <---#
			my $key = $line[$col - 1];

			#---> Warn if hash key not unique <---#
			warn "[WARNING] Entry '$key' not unique in line $.!\nDuplicate entry discared." and next if exists $HoA{$key};

			#---> Add line to hash <---#
			$HoA{$key} = \@line;
		}		
	
		#---> Close file handle <---#
		close BED;

	#---> STATUS MESSAGE <---#
	print "File '$in_file' processed.\n" unless $quiet;
	
	#---> RETURN VALUE <---#
	return \%HoA;
}