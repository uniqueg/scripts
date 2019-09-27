#!/usr/bin/perl

#==================#
#   HEADER START   #
#==================#
### Name: generate_GENCODE_gene_transcript_lookup_table.pl
### Created: Jun 14, 2013
### Modified: Jul 31, 2013
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
### Version: v2.0
### Adapted from: n/a
### Requirements: Getopt::Long
#==================#
### Description: Builds gene -> transcript lookup table from a GENCODE GTF file (only considers features of types "exon" and "transcript")
### Output: Returns a tab-separated file of the following format: gene_id TAB gene_id TAB number_of_transcripts TAB trx_1_id|trx_1_length,trx_2_id|trx_2_length, ...,trx_n_id|trx_n_length
### Output example line: ENSG00000223972	ENSG00000223972	ENST00000456328|1657,ENST00000515242|1653,ENST00000518655|1483
### Usage: perl ./perl generate_GENCODE_gene_transcript_lookup_table.pl for information on required and optional arguments
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
## Parse options
my $options_result = GetOptions (
	'usage|help' => \$usage,
	'quiet|no_verbose' => \$quiet,
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
my ($exons_trx_AoH_ref, $trx_length_hash_ref, $gene_trx_hash_ref);

#---> BODY <---#
$exons_trx_AoH_ref = &gencode_gtf_to_AoH($gtf, "exon", "transcript");
$trx_length_hash_ref = &gencode_gtf_AoH_to_trx_length_hash_ref($exons_trx_AoH_ref);
$gene_trx_hash_ref = &gencode_gtf_AoH_to_gene_trx_hash_ref($exons_trx_AoH_ref, $trx_length_hash_ref);
&print_hash_keys_alphanum_sorted($out_file, $gene_trx_hash_ref);

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
--usage|help		Show this information and die
--quiet			Shut up!
';
}
#-----------------------#
sub gencode_gtf_to_AoH {
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
sub gencode_gtf_AoH_to_trx_length_hash_ref {
### Function: Creates a hash with "gene_id TAB gene_name" as keys and arrays of corresponding "transcript_id"s as values
### Accepts: 1. [REQUIRED] GTF file; 2. [REQUIRED] Feature type to be processed (e.g. "transcript", "gene", "exon") [DEFAULT: "exon"]; all other types are discarded
### Returns: Reference to hash of arrays
### Dependencies: n/a
### Type: Generic
	#---> PASS ARGUMENTS<---#
	my @exons_trx_AoH = @{shift()};

	#---> STATUS MESSAGE <---#	
	print STDOUT "Building transcript -> length hash...\n" unless $quiet;

	#---> SUBROUTINE VARIABLES <---#
	my %HoA_trx_len;

	#---> BODY <---#

	#---> BUILD TRANSCRIPT ID -> ARRAY OF EXON LENGTHS HASH <---#
	## Iterate over every element $gtf_line (hash reference) of array @exons_trx_AoH
	foreach my $gtf_line (@exons_trx_AoH) {
	    ## Skip line if "type" not "transcript"
	    next unless $$gtf_line{"type"} eq "exon";
	    ## Issue warning and skip line if there is more than one ENSG or ENST associated with this line
	    if ( scalar @{$$gtf_line{"gene_id"}} != 1 || scalar @{$$gtf_line{"transcript_id"}} != 1 ) {
	    	print STDERR "[WARNING] Exon entry in GTF file associated with more than one gene or transcript ID:\n";
	    	## Print gene IDs
	    	print STDERR "Gene IDs:\n";
	    	foreach (@{$$gtf_line{"gene_id"}}) {
	    		print STDERR "$_\n";
	    	}
	    	## Print transcript IDs
	    	print STDERR "Transcript IDs:\n";
	    	foreach (@{$$gtf_line{"transcript_id"}}) {
	    		print STDERR "$_\n";
	    	}
	    	# Skip line
	    	next;
	    };
	    # Extract gene ID
	    my $gene_id = ${$$gtf_line{"gene_id"}}[0];
	    # Extract transcript ID
	    my $transcript_id = ${$$gtf_line{"transcript_id"}}[0];
	    ## Remove version numbers and quotes
	    $gene_id =~ s/\A\"(.*)\.\d+\"\z/$1/;
	    $transcript_id =~ s/\A\"(.*)\.\d+\"\z/$1/;
	    # Build unique transcript ID
	    my $uniq_trx_id = join "@", $gene_id, $transcript_id;
	    # Calculate exon length
	    my $exon_len = $$gtf_line{"end"} - $$gtf_line{"start"} + 1;
		## Issue warning and skip line if $exon_len is not an integer
		unless ($exon_len =~ /\A[-+]?\d+\z/) {
			print STDERR "[WARNING] Value '$exon_len' in transcript '$' is not an integer! Value skipped.\n";
			next;
		}
		# Push length of current feature to array reference associated with %HoA key 'transcript_id'	
		push @{$HoA_trx_len{$uniq_trx_id}}, $exon_len;
	}

	#---> REPLACE ARRAY OF EXON LENGTHS WITH SUM OF EXON LENGTHS <---#	
	## Traverse through each value $trx_len_array_ref of hash %HoA_trx_len
	foreach my $trx_len_array_ref ( values %HoA_trx_len ) {
	    # Initiate sum of exon lengths
	    my $len_sum = 0;
	    ## Traverse through each element $$exon_len of array @$trx_len_array_ref
	    foreach my $exon_len ( @$trx_len_array_ref ) {
			# Add it to $sum
			$len_sum += $exon_len;
		}
		# Replace reference to array of exon lengths with sum of exon lengths
		$trx_len_array_ref = $len_sum;
    }
		
	#---> STATUS MESSAGE <---#
	print STDOUT "Transcript -> length hash built.\n\n" unless $quiet;
	
	#--->  RETURN VALUE  <---#
	return \%HoA_trx_len;	
}
#-----------------------#
sub gencode_gtf_AoH_to_gene_trx_hash_ref {
### Function: Generates 'exon_id' -> transcripts hash from GENCODE GTF array of hashes
### Accepts: 1. Reference to GENCODE GTF array of hashes
### Returns: Reference to exon -> transcripts hash
### Dependencies: n/a
### Type: Specialized
	#---> PASS ARGUMENTS ---#
	my @exons_trx_AoH = @{shift()};
	my %trx_length_hash = %{shift()};
		
	#---> STATUS MESSAGE <---#
	print STDOUT "Building gene -> transcripts hash..." . "\n" unless $quiet;

	#---> SUBROUTINE VARIABLES <---#
	my %gene_trx_HoA;
	
	#---> BODY <---#

	#---> BUILD GENE ID -> ARRAY OF TRANSCRIPT IDS HASH <---#
	## Iterate over every element $gtf_line (hash reference) of array @exons_trx_AoH
	foreach my $gtf_line (@exons_trx_AoH) {
	    ## Skip line if "type" not "exon"
	    next unless $$gtf_line{"type"} eq "transcript";
	    ## Issue warning and skip line if there is more than one ENSG or ENST associated with this line
	    if ( scalar @{$$gtf_line{"gene_id"}} != 1 || scalar @{$$gtf_line{"transcript_id"}} != 1 ) {
	    	print STDERR "[WARNING] Exon entry in GTF file associated with more than one gene or transcript ID:\n";
	    	## Print gene IDs
	    	print STDERR "Gene IDs:\n";
	    	foreach (@{$$gtf_line{"gene_id"}}) {
	    		print STDERR "$_\n";
	    	}
	    	## Print transcript IDs
	    	print STDERR "Transcript IDs:\n";
	    	foreach (@{$$gtf_line{"transcript_id"}}) {
	    		print STDERR "$_\n";
	    	}
	    	# Skip line
	    	next;
	    };
	    # Extract gene ID
	    my $gene_id = ${$$gtf_line{"gene_id"}}[0];
	    # Extract transcript ID
	    my $transcript_id = ${$$gtf_line{"transcript_id"}}[0];
	    ## Remove version numbers and quotes
	    $gene_id =~ s/\A\"(.*)\.\d+\"\z/$1/;
	    $transcript_id =~ s/\A\"(.*)\.\d+\"\z/$1/;
		# Push $transcript_id of current line to array referenced by %gene_trx_HoA hash key $gene_id
		push @{$gene_trx_HoA{$gene_id}}, $transcript_id;
	}

	#---> REPLACE ARRAY OF EXON LENGTHS WITH TAB-SEPARATED STRING OF GENE ID, TRANSCRIPT COUNT & TRANSCRIPT ID|LENGTH... <---#	
	
	## Traverse through each key $gene_id of hash %gene_trx_HoA
	foreach my $gene_id ( keys %gene_trx_HoA ) {
	    # Count transcripts
		my $trx_count = scalar @{$gene_trx_HoA{$gene_id}};
		# Declare array to hold transcript|length strings
		my @trx_len;
		## Traverse through each element $transcript_id of array @{$gene_trx_HoA{$gene_id}}
		foreach my $transcript_id ( @{$gene_trx_HoA{$gene_id}} ) {
			# Declare switch that decides when/if to delete a missing/misbehaving transcript
			my $skip = 0;
		    # Generate lookup value for transcript -> length hash
	    	my $trx_len_hash_key = join "@", $gene_id, $transcript_id;
			# Declare string to hold transcript length
			my $transcript_len;
			# Obtain length of current transcript from transcript -> length hash
	    	if (defined $trx_length_hash{$trx_len_hash_key}) {
	    		$transcript_len = $trx_length_hash{$trx_len_hash_key};
	    	}
	    	else {
	    		print STDERR "[WARNING] No length found for the following transcript:\n$trx_len_hash_key\nTranscript skipped.\n";
	    		$skip = 1;
	    	}
	    	# Append current transcript ID and length to transcript|length string	    
			push @trx_len, ($transcript_id . "|" . $transcript_len) unless $skip;
		}
		# Join transcript|length array to transcript|length string
		my $trx_len_str = join ",", @trx_len;
		# Replace reference to array of transcripts with tab-separated string of the following format: gene ID, transcript count & transcript1 ID|length, transcript2 ID|length, ..., transcriptN ID|length
		$gene_trx_HoA{$gene_id} = $gene_id . "\t" . $trx_count . "\t" . $trx_len_str;
	}

	#---> STATUS MESSAGE <---#
	print STDOUT "Gene -> transcripts hash built." . "\n\n" unless $quiet;
	
	#---> RETURN VALUE <---#
	return \%gene_trx_HoA;
}
#-----------------------#
sub print_hash_keys_alphanum_sorted {
### Function: Writes key -> value hash pairs to output file (one pair per line)
### Accepts: 1. Output file; 2. Hash reference; 3. Separator [DEFAULT: TAB]
### Returns: n/a
### Dependencies: n/a
### Type: Generic
	#---> PASS ARGUMENTS<---#
	my $out_file = shift;
	my $hash_ref = shift;
	my $sep = defined($_[0]) ? $_ : "\t";

	#---> STATUS MESSAGE <---#	
	print STDOUT "Writing output to file '$out_file'...\n" unless $quiet;

	#---> BODY <---#
	# Open GTF filehandle
	open OUT, ">" . $out_file;
	## Iterate over each key value pair of 'exon -> transcript' hash of array references
	foreach my $key (sort { $a cmp $b } keys %$hash_ref) {
		# Print hash key TAB hash value
		print OUT $key . $sep . $$hash_ref{$key} . "\n";
	}
	# Close OUT filehandle
	close OUT;
	
	#---> STATUS MESSAGE <---#
	print STDOUT "File '$out_file' written.\n\n" unless $quiet;
}
#=======================#
#    SUBROUTINES END    #
#=======================#
