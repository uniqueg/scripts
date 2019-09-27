#!/usr/bin/perl

#==================#
#   HEADER START   #
#==================#
### Name: generate_GENCODE_exon_coordinates_table.pl
### Created: Oct 4, 2013
### Modified: Oct 4, 2013
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
### Version: v2.0
### Adapted from: generate_GENCODE_gene_transcript_lookup_table.pl
### Requirements: Getopt::Long
#==================#
### Description: Builds a BED file (0-based!) of exon junctions from "exon" entries of a GTF file (e.g. GENCODE); exon junctions can be printed in genome or transcriptome coordinates
### Output: Returns a tab-separated file of the following format: transcript TAB id start TAB stop
### Usage: perl ./perl generate_GENCODE_exon_junctions_bed.pl for information on required and optional arguments
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
## Set default values for options and arguments
my $usage = '';
my $quiet = '';
my $gtf = '';
my $out_file = '';
my $trx_coord = 0;
my $score = 0;
## Parse options
my $options_result = GetOptions (
	'usage|help' => \$usage,
	'quiet|no_verbose' => \$quiet,
	'trx' => \$trx_coord,
	'score=i' => \$score,
	#-----------------------#
	'gtf=s' => \$gtf,
	'out_file=s' => \$out_file
);
# Die with usage information if --usage option is set or if options parsing returned FALSE
die $usage_info if $usage || !$options_result;
# Die with usage information if required arguments are not set
die $usage_info if !$gtf || !$out_file; 
#==========================#
#    PRE-REQUISITES END    #
#==========================#


#================#
#   MAIN START   #
#================#
#---> STATUS MESSAGE <---#
print "Starting 'generate_GENCODE_gene_transcript_lookup_table.pl'...\n\n" unless $quiet;

#---> MAIN VARIABLES <---#
my ($exons_AoH_ref, $trx_coord_HoA_ref);

#---> BODY <---#
$exons_AoH_ref = &gencode_GTF_to_AoH($gtf, "exon");
$trx_coord_HoA_ref = &gencode_GTF_exons_AoH_to_trx_coord_HoA_ref($exons_AoH_ref);
&print_HoA_values_keys_alphanum_sorted($out_file, $trx_coord_HoA_ref);

#---> STATUS MESSAGE <---#
print "Done.\n" unless $quiet;

#---> PROGRAM EXIT <---#
exit;
#================#
#    MAIN END    #
#================#


#=======================#
#   SUBROUTINES START   #
#=======================#
sub usage {
### Function: Returns usage information for current program
### Accepts: n/a
### Returns: String with usage information
### Dependencies: n/a
### Type: Specialized
'Usage: perl ./generate_GENCODE_gene_transcript_lookup_table.pl [OPTIONS] --gtf [FILE] --out_file [FILE]
==================================================
Required arguments:
--gtf			GTF input file [FILE]
--out_file		Output filename [FILE]
==================================================
Optional arguments:
--trx	Print out coordinates relative to the transcript
--score	Score variable to be set in BED output (only relevant if --trx not set; default = 0)
--usage|help		Show this information and die
--quiet			Shut up!
';
}
#-----------------------#
sub gencode_GTF_to_AoH {
### Function: Reads a GENCODE GTF file into an array of hashes with the field values of one line populating one hash each; the used keys are: "chr", "source", "type", "start", "end", "score", "strand", "phase" as well as all the attribute keys found within a particular line ("e.g. "transcript_id", "gene_name"); attribute values are stored as array references in each hash to account for the possible occurrence of multiple values associated with a specific attribute key
### Accepts: 1. [REQUIRED] GTF file; 2+. feature types to be processed (e.g. "transcript", "gene", "exon") [DEFAULT: "", i.e. all feature types]; if indicated, all other feature types are discarded
### Returns: Reference to array of hashes; attribute 
### Dependencies: n/a
### Type: Generic

	#---> PASS ARGUMENTS<---#
	my $gtf_file = shift;
	my @selected_types = defined($_[0]) ? @_ : "";

	#---> STATUS MESSAGE <---#	
	print STDOUT "Processing GTF file '$gtf_file'...\n" unless $quiet;

	#---> SUBROUTINE VARIABLES <---#
	my @selected;	

	#---> BODY <---#

	# Open GTF filehandle
	open GTF, $gtf_file or die "Can't open $gtf.\n";
	## Iterate over GTF file line by line
	while (<GTF>) {

		# Ignore header
		next if(/^##/);
		# Remove trailing entry separator		
		chomp;
		
		# Split line by tabs
		my ($chr, $source, $type, $start, $end, $score, $strand, $phase, $attributes) = split("\t");

		# Skip line if not of desired type
		unless ($selected_types[0] eq "") {
			next unless grep /$type/, @selected_types;
		}

		# Store all fields but attributes in hash
		my %fields = (
			chr				=> $chr,
			source		 => $source,
			type			 => $type,
			start			=> $start,
			end				=> $end,
			score			=> $score,
			strand		 => $strand,
			phase			=> $phase,
		);
		
		# Split attributes by semicolon
		my @add_attributes = split ";", $attributes;

		## Store attribute keys and values in hash
		for (my $i = 0; $i < scalar @add_attributes; $i++) {
			 $add_attributes[$i] =~ /^(.+)\s(.+)$/;
			 my $c_type	= $1;
			 my $c_value = $2;
			 $c_type =~ s/^\s//;
			 if($c_type	&& $c_value){
				 if(!exists($fields{$c_type})){
					 $fields{$c_type} = [];
				 }
				 push(@{ $fields{$c_type} }, $c_value);
			 }
		}
	
		# Push hash to array of selected entries
		push @selected, \%fields;

	}		

	# Close GTF filehandle
	close GTF;
	
	#---> STATUS MESSAGE <---#
	print STDOUT "GTF file '$gtf_file' processed.\n\n" unless $quiet;

	#--->	RETURN VALUE	<---#
	return \@selected;

}
#-----------------------#
sub gencode_GTF_exons_AoH_to_trx_coord_HoA_ref {
### Function: Creates a hash of arrays of the form: transcript ID > exon string (in minimal BED format) for each exon of a GENCODE GTF file
### Accepts: 1. Reference to GENCODE GTF exons array of hashes (one exon per line)
### Returns: Reference to hash of arrays: transcript ID > exon coordinates (transcript ID, start, stop)
### Dependencies: n/a
### Type: Specialized
	#---> PASS ARGUMENTS<---#
	my @exons_AoH = @{shift()};

	#---> STATUS MESSAGE <---#	
	print STDOUT "Building transcript ID -> exon coordinates hash...\n" unless $quiet;

	#---> SUBROUTINE VARIABLES <---#
	my %HoA_exon_coord;

	#---> BODY <---#

	#---> BUILD TRANSCRIPT ID -> ARRAY OF EXON LENGTHS <---#
	## Iterate over every element $gtf_line (hash reference) of array @exons_trx_AoH
	foreach my $gtf_line (@exons_AoH) {
	    ## Skip line if "type" not "exon"
	    next unless $gtf_line->{"type"} eq "exon";
	    ## Issue warning and skip line if there is more than one ENSG or ENST associated with this line
	    if ( scalar @{$gtf_line->{"gene_id"}} != 1 || scalar @{$gtf_line->{"transcript_id"}} != 1 ) {
	    	print STDERR "[WARNING] Exon entry in GTF file associated with more than one gene or transcript ID:\n";
	    	## Print gene IDs
	    	print STDERR "Gene IDs:\n";
	    	foreach (@{$gtf_line->{"gene_id"}}) {
	    		print STDERR "$_\n";
	    	}
	    	## Print transcript IDs
	    	print STDERR "Transcript IDs:\n";
	    	foreach (@{$gtf_line->{"transcript_id"}}) {
	    		print STDERR "$_\n";
	    	}
	    	# Skip line
	    	print STDERR "Line skipped!\n";
	    	next;
	    };
	    # Extract transcript ID
	    my $transcript_id = ${$gtf_line->{"transcript_id"}}[0];
	    # Remove version number and quotes
	    $transcript_id =~ s/\A\"(.*)\.\d+\"\z/$1/;
	    # Initialize coordinates variable;
	    my $coords;
	    ## Build hash value string in minimal BED format (transcript coordinates)
	    if ($trx_coord) {
		    $coords = $transcript_id . "\t" . ($gtf_line->{"start"} -1) . "\t" . $gtf_line->{"end"};
	    }
	    ## Build hash value string in BED format (genomic coordinates)
	    else {
	    	$coords = $gtf_line->{"chr"} . "\t" . ($gtf_line->{"start"} - 1) . "\t" . $gtf_line->{"end"} . "\t" . $transcript_id . "\t" .  $score . "\t" . $gtf_line->{"strand"};
	    }
		# Push length of current feature to array reference associated with %HoA key 'transcript_id'	
		push @{$HoA_exon_coord{$transcript_id}}, $coords;
	}

	#---> SORT EXONS <---#
	## Interate over every set of exons
	foreach my $exon_bed_string_array_ref (values %HoA_exon_coord) {
	    ## Sort differently for transcript or genomic coordinates...
	    if ($trx_coord) {
			# Sort exons by transcript name, then transcript start, then end positions
		    $exon_bed_string_array_ref = &tab_string_array_ref_to_array_of_sorted_arrays($exon_bed_string_array_ref, 1, "alpha", "asc", 2, "num", "asc", 3, "num", "asc");
	 	    # Initiate sum of exon lengths
		    my $len_sum = 0;
	  	    ## Traverse through each element $exon_len of array @$trx_len_array_ref
		    foreach my $exon_bed_string ( @$exon_bed_string_array_ref ) {
				# Split by tab
				my @coords = split "\t", $exon_bed_string;
				# Add length of exon to $sum
				$len_sum += $coords[2] - $coords[1];
				## Replace genomic with transcript start and stop coordinates
				$coords[1] = $len_sum -1;
				$coords[2] = $len_sum;
				$exon_bed_string = join "\t", @coords;
			}
			# Remove final exon (no junction!)
			pop @$exon_bed_string_array_ref;
	    }
	    else {
			# Sort exons by transcript name, then genomic start, then end positions
		    $exon_bed_string_array_ref = &tab_string_array_ref_to_array_of_sorted_arrays($exon_bed_string_array_ref, 4, "alpha", "asc", 2, "num", "asc", 3, "num", "asc");
	    }
    }

	#---> STATUS MESSAGE <---#
	print STDOUT "Transcript ID -> exon coordinates hash built.\n\n" unless $quiet;
	
	#--->  RETURN VALUE  <---#
	return \%HoA_exon_coord;	
}
#-----------------------#
sub tab_string_array_ref_to_array_of_sorted_arrays {
### Function: Reads a reference to an array of tab-separated strings into an arrays of sorted arrays
### Accepts: 1. reference to array of tab-seaparated strings; 2-4. (or more; OPTIONAL!) one or more sets of sorting information (multiples of 3) in the order of priority: A. Reference column for sorting [INTEGER]; B. Sorting type [STRING | ALLOWED: "num", "alpha"]; C. Sorting order [STRING | ALLOWED: "asc", "desc"]; DEFAULT: NO SORTING!
### Returns: 2. Reference to sorted array of tab-separated strings
### Dependencies: n/a
### Type: Generic
	#---> PASS ARGUMENTS ---#
	my ($array_ref, @sort) = @_;	
	
	#---> SUBROUTINE VARIABLES <---#
	my @AoA;	
	
	#---> BODY <---#
	## Crawl through array line by line
	foreach my $line (@$array_ref) {
		# Split line into array @line by tabs
		my @line = split /\t/, $line;
		# Push line values array @line to respective hash value / outer array
		push @AoA, [ @line ];
	}		
	#-----------------------#
	## Sort inner arrays of @AoA
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
			# Sort...
			@AoA = sort { $a->[$col] <=> $b->[$col] } @AoA;
		}
		# Alphanumerically in ascending order
		elsif ( $type eq "alpha" && $order eq "asc" ) {
			# Sort...
			@AoA = sort { $a->[$col] cmp $b->[$col] } @AoA;
		}
		# Numerically in descending order
		elsif ( $type eq "num" && $order eq "desc" ) {
			# Sort...
			@AoA = sort { $b->[$col] <=> $a->[$col] } @AoA;
		}
		# Alphanumerically in descending order
		elsif ( $type eq "alpha" && $order eq "desc" ) {
			# Sort...
			@AoA = sort { $b->[$col] cmp $a->[$col] } @AoA;
		}
		# Illegal sorting type and/or order 
		else {
			print "$order\t$col\t$type\n";
			print STDERR "Invalid value for sorting type (only 'num' and 'alpha' allowed!) or order (only 'asc' and 'desc' allowed!)" . "\n";
		}
	}
	#-----------------------#
	## Join inner array elements by TABs
	foreach my $inner_array_ref (@AoA) {
    	$inner_array_ref = join "\t", @$inner_array_ref;
	}

	#---> RETURN VALUE <---#
	return \@AoA;
}
#-----------------------#
sub print_HoA_values_keys_alphanum_sorted {
### Function: Writes values (i.e. arrays) of a hash of arrays to output file, with one array element per line; order: alphanumerical sorting of hash keys 
### Accepts: 1. Output file; 2. HoA reference
### Returns: n/a
### Dependencies: n/a
### Type: Generic
	#---> PASS ARGUMENTS<---#
	my $out_file = shift;
	my $hash_ref = shift;

	#---> STATUS MESSAGE <---#	
	print STDOUT "Writing output to file '$out_file'...\n" unless $quiet;

	#---> BODY <---#
	# Open GTF filehandle
	open OUT, ">" . $out_file;
	## Sort keys and iterate over each key of the hash of array references
	foreach my $key (sort { $a cmp $b } keys %$hash_ref) {
			## Iterate over each array element of the current inner array
			foreach my $element (@{$$hash_ref{$key}}) {
				# Print element, followed by a newline character
				print OUT $element . "\n";
			}
	}
	# Close OUT filehandle
	close OUT;
	
	#---> STATUS MESSAGE <---#
	print STDOUT "File '$out_file' written.\n\n" unless $quiet;
}
#=======================#
#    SUBROUTINES END    #
#=======================#
