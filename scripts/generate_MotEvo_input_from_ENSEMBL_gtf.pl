#!/usr/bin/perl

#==================#
#   HEADER START   #
#==================#
### Name: generate_MotEvo_input_regions.pl
### Created: Aug 1, 2013
### Modified: Aug 1, 2013
### Author: Alexander Kanitz, Stamatina Stamelaki
### Company: Zavolan Group, Biozentrum, University of Basel
### Version: v1.0
### Adapted from: n/a
### Requirements: Getopt::Long
#==================#
### Description: Generates a MotEvo compatible regions file from an ENSEMBL GTF annotation file; feature types (e.g. "exon", "transcript") can be selected as well as flanking regions etc. selected
### Output: Returns a tab-separated file of the following format: chromosome TAB start TAB end TAB strand TAB extra_info
### Output example line: chr1	11869	12227	+	ENSG00000223972_chr1_+_11869_12227
### Usage: perl ./perl generate_MotEvo_input_from_ENSEMBL_gtf.pl --usage for information on required and optional arguments
#==================#
#    HEADER END	   #
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
my $type = '';
my $up = 0;
my $down = 0;
my $left = 0;
my $right = 0;
my $sep;
my $gtf;
my $prefix = '';
## Parse options
my $options_result = GetOptions (
	'usage|help' => \$usage,
	'quiet|no_verbose' => \$quiet,
	'type=s' => \$type,
	'up=s' => \$up,
	'down=s' => \$down,
	'left=s' => \$left,
	'right=s' => \$right,
	'sep|separately' => \$sep,
	#-----------------------#
	'gtf=s' => \$gtf,
	'prefix=s' => \$prefix
);
# Die with usage information if --usage option is set or if options parsing returned FALSE
die $usage_info if $usage || !$options_result;
# Die with usage information if required arguments are not set
die "[ERROR] Required arguments missing! Execution halted.\n\n$usage_info" if !$gtf || !$prefix; 
# Die with warning if --up or --down are used together with --left or --right and --sep is not set
die "[ERROR] Cannot use any of --up and --down together with any of --left and --right if --sep is not set! Execution halted.\n\n$usage_info" if ! $sep && ($up || $down) && ($left || $right);
# Split up $type parameter
my @types = split /\|/, $type;
#==========================#
#	PRE-REQUISITES END	#
#==========================#


#================#
#   MAIN START   #
#================#
#---> STATUS MESSAGE <---#
print "Starting 'generate_MotEvo_input_from_ENSEMBL_gtf.pl'...\n\n" unless $quiet;

#---> MAIN VARIABLES <---#
my ($region_AoH_ref, $region_type_HoH_ref);

#---> BODY <---#
$region_AoH_ref = ENSEMBL_gtf_to_AoH($gtf, @types);

print $region_AoH_ref . "\n";

foreach my $hash_ref (@$region_AoH_ref) {
	foreach my $key (keys %$hash_ref) {
		print "KEY: $key\tVALUE: " . $$hash_ref->{$key} . "\n";
	}
}

$region_type_HoH_ref = ENSEMBL_gtf_AoH_to_region_type_HoH_ref($region_AoH_ref);
&print_output($region_type_HoH_ref);

#---> STATUS MESSAGE <---#
print "Done.\n" unless $quiet;

#---> PROGRAM EXIT <---#
exit;
#================#
#	MAIN END	#
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
'Usage: perl ./generate_MotEvo_input_from_ENSEMBL_gtf.pl [OPTIONS] --gtf [FILE] --out_file [FILE]

Generates a MotEvo compatible regions file from a ENSEMBL GTF annotation file.

==================================================
Required arguments:
--gtf [GTF]		GTF input filename
--prefix [PATH/STRING]	Output file prefix
==================================================
Optional arguments:
--type [QUOTED STRING]	Select which feature types to select, e.g. "exon", "gene" etc. Multiple values are accepted, separated by the pipe character "|". (DEFAULT: empty string, here meaning :all)
--up [INT]		Extends each region for INT bases at the 5\' end ("upstream"). Cannot be used together with --left or --rigt unless --sep is set. Negative values shrink the region at the specified side (regions with length < 1 are skipped). Features not annotated on the "+" or "-" strand (e.g, "*" strand) are skipped with --up and --down.
--down [INT]		Extends each region for INT bases at the 3\' end ("downstream"). See --up for more details.
--left [INT]		Extends each region for INT bases to the left ("upstream" on the "+", "downstream" on the "-" strand). See --up for more detailed.
--right [INT]		Extends each region for INT bases to the right ("downstream" on the "+", "upstream" on the "-" strand). See --left for more details.
--sep|separately	Together with --up, --down, --left and/or --right writes out separate output files for the original region and the selected flanking regions
--usage|help		Show this information and die
--quiet			Shut up!
';
}
#-----------------------#
sub ENSEMBL_gtf_to_AoH {
### Function: Reads an ENSEMBL GTF file into an array of hashes with the field values of one line populating one hash each; the used keys are: "chr", "source", "type", "start", "end", "score", "strand", "phase" as well as all the attribute keys found within a particular line ("e.g. "transcript_id", "gene_name"); attribute values are stored as array references in each hash to account for the possible occurrence of multiple values associated with a specific attribute key
### Accepts: 1. [REQUIRED] GTF file; 2+. feature types to be processed (e.g. "transcript", "gene", "exon") [DEFAULT: "", i.e. all feature types]; if indicated, all other feature types are discarded
### Returns: Reference to array of hashes
### Dependencies: n/a
### Type: Generic

	#---> PASS ARGUMENTS <---#
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
		next if(/^\#\!/);
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
			chr			=> $chr,
			source		=> $source,
			type		=> $type,
			start		=> $start,
			end			=> $end,
			score		=> $score,
			strand		=> $strand,
			phase		=> $phase,
		);
		
		# Split attributes by semicolon
		my @attributes = split /\s*;\s*/, $attributes;

		## Store attribute keys and values in hash
		for my $attribute (@attributes) {
			warn "[WARNING] Non-standard attribute format in line " . $. . ". Entry skipped.\n" and next unless $attribute =~ /^\s*(\S+)\s+(\"[^\"]*\")\s*$/;
			my $attr_type = $1;
			my $attr_value = $2;
			$attr_value =~ s/\"//g;
			push(@{ $fields{$attr_type} }, $attr_value);
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
sub ENSEMBL_gtf_AoH_to_region_type_HoH_ref {
### Function: Determines IDs and coordinates for each selected region type ("main", "up", "down", "left", "right")
### Accepts: 1. Reference to ENSEMBL GTF array of hashes
### Returns: Region type (keys: "main", "up", "down", "left", "right") hash of hashes (ID -> coordinates)
### Dependencies: n/a
### Type: Specialized
	#---> PASS ARGUMENTS ---#
	my @features_AoH = @{shift()};
		
	#---> STATUS MESSAGE <---#
	print STDOUT "Building ID -> coordinates hashes..." . "\n" unless $quiet;

	#---> SUBROUTINE VARIABLES <---#
	my %region_type_HoH;
	
	#---> BODY <---#
	## Traverse through each element $gtf_line of array @features_AoH
	foreach my $gtf_line ( @features_AoH ) {
		## Make coordinates 0-based and open-ended (BED-like);
		$$gtf_line{"start"} -= 1;
		## Issue warning and skip line if there is more than one ENSG or ENST associated with this line
		if ( scalar @{$$gtf_line{"gene_id"}} != 1 || scalar @{$$gtf_line{"transcript_id"}} != 1 ) {
			print STDERR "[WARNING] Entry in GTF file associated with more than one gene or transcript ID:\n";
			## Print gene IDs
			print STDERR "Gene IDs:\n";
			foreach (@{$$gtf_line{"gene_id"}}) {
				print STDERR "$_\n";
			}
		};
		# Extract gene ID
		my $gene_id = ${$$gtf_line{"gene_id"}}[0];
		# Remove version numbers and quotes
		$gene_id =~ s/\A\"(.*)\.\d+\"\z/$1/;
		## Build inner hash pairs (ID -> coordinates) and add to %region_type_HoH
		if ($sep) {
			my $hash_pair_up_array_ref = &extract_upstream($gene_id, $$gtf_line{"chr"}, $$gtf_line{"start"}, $$gtf_line{"end"}, $$gtf_line{"strand"}) if $up;
			${$region_type_HoH{"up"}}{$$hash_pair_up_array_ref[0]} = $$hash_pair_up_array_ref[1] if $hash_pair_up_array_ref;
			my $hash_pair_down_array_ref = &extract_downstream($gene_id, $$gtf_line{"chr"}, $$gtf_line{"start"}, $$gtf_line{"end"}, $$gtf_line{"strand"}) if $down;
			${$region_type_HoH{"down"}}{$$hash_pair_down_array_ref[0]} = $$hash_pair_down_array_ref[1] if $hash_pair_down_array_ref;;
			my $hash_pair_left_array_ref = &extract_left($gene_id, $$gtf_line{"chr"}, $$gtf_line{"start"}, $$gtf_line{"end"}, $$gtf_line{"strand"}) if $left;
			${$region_type_HoH{"left"}}{$$hash_pair_left_array_ref[0]} = $$hash_pair_left_array_ref[1] if $hash_pair_left_array_ref;;
			my $hash_pair_right_array_ref = &extract_right($gene_id, $$gtf_line{"chr"}, $$gtf_line{"start"}, $$gtf_line{"end"}, $$gtf_line{"strand"}) if $right;
			${$region_type_HoH{"right"}}{$$hash_pair_right_array_ref[0]} = $$hash_pair_right_array_ref[1] if $hash_pair_right_array_ref;;
			my $hash_pair_region_only_array_ref = &extract_region($gene_id, $$gtf_line{"chr"}, $$gtf_line{"start"}, $$gtf_line{"end"}, $$gtf_line{"strand"});
			${$region_type_HoH{"main"}}{$$hash_pair_region_only_array_ref[0]} = $$hash_pair_region_only_array_ref[1] if $hash_pair_region_only_array_ref;;
		}
		elsif ($up || $down) {
			my $hash_pair_union_up_down_array_ref = &union_up_down($gene_id, $$gtf_line{"chr"}, $$gtf_line{"start"}, $$gtf_line{"end"}, $$gtf_line{"strand"});
			${$region_type_HoH{"main"}}{$$hash_pair_union_up_down_array_ref[0]} = $$hash_pair_union_up_down_array_ref[1] if $hash_pair_union_up_down_array_ref;
		}
		elsif ($left || $right) {
			my $hash_pair_union_left_right_array_ref = &union_left_right($gene_id, $$gtf_line{"chr"}, $$gtf_line{"start"}, $$gtf_line{"end"}, $$gtf_line{"strand"});
			${$region_type_HoH{"main"}}{$$hash_pair_union_left_right_array_ref[0]} = $$hash_pair_union_left_right_array_ref[1] if $hash_pair_union_left_right_array_ref;
		}
		else {
			my $hash_pair_region_only_array_ref = &extract_region($gene_id, $$gtf_line{"chr"}, $$gtf_line{"start"}, $$gtf_line{"end"}, $$gtf_line{"strand"});
			${$region_type_HoH{"main"}}{$$hash_pair_region_only_array_ref[0]} = $$hash_pair_region_only_array_ref[1];		
		}
	}
	
	#---> STATUS MESSAGE <---#
	print STDOUT "ID -> coordinates hashes built." . "\n\n" unless $quiet;
	
	#---> RETURN VALUE <---#
	return \%region_type_HoH;
}
#-----------------------#
sub extract_region {
### Function: Extracts region from GTF file
### Accepts: 1. Gene ID; 2. Chromosome; 3. Start position; 4. End position; 5. Strand
### Returns: Reference to array holding a hash key [0] value [1] pair of the form: gene_id_chr_str_start_end ->gene_id \t chr \t start \t end \t str
### Dependencies: n/a
### Type: Specialized
	#---> PASS ARGUMENTS <---#
	my ($gene_id, $chr, $start, $end, $strand) = @_;

	#---> DECLARE VARIABLES <---
	my @return;

	#---> BODY <---#
	$return[0] = join "_", $gene_id, $chr, $strand, $start, $end;
	$return[1] = join "\t", $gene_id, $chr, $start, $end, $strand;
	
	#---> RETURN VALUE <---#
	if ( ! $end > $start ) {
		return 0;
	}
	else {
		return \@return;
	}
}
#-----------------------#
sub extract_upstream {
### Function: Extracts upstream region defined by --up
### Accepts: 1. Gene ID; 2. Chromosome; 3. Start position; 4. End position; 5. Strand
### Returns: Reference to array holding a hash key [0] value [1] pair of the form: gene_id_chr_str_start_end ->gene_id \t chr \t start \t end \t str
### Dependencies: n/a
### Type: Specialized
	#---> PASS ARGUMENTS <---#
	my ($gene_id, $chr, $start, $end, $strand) = @_;

	#---> DECLARE VARIABLES <---
	my @return;

	#---> BODY <---#
	if ($strand eq "+") {
		$end = $start;
		$start = $start - $up;
	}
	elsif ($strand eq "-") {
		$start = $end;
		$end = $end + $up;
	}
	else {
		print STDERR "[WARNING]	Entry is not annotated on the \"+\" or \"-\" strand! Parameter --up makes no sense without strand information.\nEntry skipped.";
	}

	$return[0] = join "_", $gene_id, $chr, $strand, $start, $end, $up, "up";
	$return[1] = join "\t", $gene_id, $chr, $start, $end, $strand;
	
	#---> RETURN VALUE <---#
	if ( ! $end > $start ) {
		return 0;
	}
	else {
		return \@return;
	}
}
#-----------------------#
sub extract_downstream {
### Function: Extracts downstream region defined by --down
### Accepts: 1. Gene ID; 2. Chromosome; 3. Start position; 4. End position; 5. Strand
### Returns: Reference to array holding a hash key [0] value [1] pair of the form: gene_id_chr_str_start_end ->gene_id \t chr \t start \t end \t str
### Dependencies: n/a
### Type: Specialized
	#---> PASS ARGUMENTS <---#
	my ($gene_id, $chr, $start, $end, $strand) = @_;

	#---> DECLARE VARIABLES <---
	my @return;

	#---> BODY <---#
	if ($strand eq "+") {
		$start = $end;
		$end = $end + $down;
	}
	elsif ($strand eq "-") {
		$end = $start;
		$start = $start - $down;
	}
	else {
		print STDERR "[WARNING]	Entry is not annotated on the \"+\" or \"-\" strand! Parameter --up makes no sense without strand information.\nEntry skipped.";
	}

	$return[0] = join "_", $gene_id, $chr, $strand, $start, $end, $down, "down";
	$return[1] = join "\t", $gene_id, $chr, $start, $end, $strand;
	
	#---> RETURN VALUE <---#
	if ( ! $end > $start ) {
		return 0;
	}
	else {
		return \@return;
	}
}
#-----------------------#
sub extract_left {
### Function: Extracts upstream region defined by --left
### Accepts: 1. Gene ID; 2. Chromosome; 3. Start position; 4. End position; 5. Strand
### Returns: Reference to array holding a hash key [0] value [1] pair of the form: gene_id_chr_str_start_end ->gene_id \t chr \t start \t end \t str
### Dependencies: n/a
### Type: Specialized
	#---> PASS ARGUMENTS <---#
	my ($gene_id, $chr, $start, $end, $strand) = @_;

	#---> DECLARE VARIABLES <---
	my @return;

	#---> BODY <---#
	$end = $start;
	$start = $start - $left;

	$return[0] = join "_", $gene_id, $chr, $strand, $start, $end, $left, "left";
	$return[1] = join "\t", $gene_id, $chr, $start, $end, $strand;
	
	#---> RETURN VALUE <---#
	if ( ! $end > $start ) {
		return 0;
	}
	else {
		return \@return;
	}
}
#-----------------------#
sub extract_right {
### Function: Extracts upstream region defined by --right
### Accepts: 1. Gene ID; 2. Chromosome; 3. Start position; 4. End position; 5. Strand
### Returns: Reference to array holding a hash key [0] value [1] pair of the form: gene_id_chr_str_start_end ->gene_id \t chr \t start \t end \t str
### Dependencies: n/a
### Type: Specialized
	#---> PASS ARGUMENTS <---#
	my ($gene_id, $chr, $start, $end, $strand) = @_;

	#---> DECLARE VARIABLES <---
	my @return;

	#---> BODY <---#
	$start = $end;
	$end = $end + $right;

	$return[0] = join "_", $gene_id, $chr, $strand, $start, $end, $right, "right";
	$return[1] = join "\t", $gene_id, $chr, $start, $end, $strand;
	
	#---> RETURN VALUE <---#
	if ( ! $end > $start ) {
		return 0;
	}
	else {
		return \@return;
	}
}
#-----------------------#
sub union_up_down {
### Function: Calculates union of region defined by the original region, --up and/or --down
### Accepts: 1. Gene ID; 2. Chromosome; 3. Start position; 4. End position; 5. Strand
### Returns: Reference to array holding a hash key [0] value [1] pair of the form: gene_id_chr_str_start_end ->gene_id \t chr \t start \t end \t str
### Dependencies: n/a
### Type: Specialized
	#---> PASS ARGUMENTS <---#
	my ($gene_id, $chr, $start, $end, $strand) = @_;

	#---> DECLARE VARIABLES <---
	my @return;

	#---> BODY <---#
	if ($strand eq "+") {
		$start = $start - $up;
		$end = $end + $down;
	}
	elsif ($strand eq "-") {
		$start = $start - $down;
		$end = $end + $up;
	}
	else {
		print STDERR "[WARNING]	Entry is not annotated on the \"+\" or \"-\" strand! Parameter --up makes no sense without strand information.\nEntry skipped.";
	}

	$return[0] = join "_", $gene_id, $chr, $strand, $start, $end, $up, "up", $down, "down";
	$return[1] = join "\t", $gene_id, $chr, $start, $end, $strand;
	
	#---> RETURN VALUE <---#
	if ( ! $end > $start ) {
		return 0;
	}
	else {
		return \@return;
	}
}
#-----------------------#
sub union_left_right {
### Function: Calculates union of region defined by the original region, --left and/or --right
### Accepts: 1. Gene ID; 2. Chromosome; 3. Start position; 4. End position; 5. Strand
### Returns: Reference to array holding a hash key [0] value [1] pair of the form: gene_id_chr_str_start_end ->gene_id \t chr \t start \t end \t str
### Dependencies: n/a
### Type: Specialized
	#---> PASS ARGUMENTS <---#
	my ($gene_id, $chr, $start, $end, $strand) = @_;

	#---> DECLARE VARIABLES <---
	my @return;

	#---> BODY <---#
	$start = $start - $left;
	$end = $end + $right;

	$return[0] = join "_", $gene_id, $chr, $strand, $start, $end, $left, "left", $right, "right";
	$return[1] = join "\t", $gene_id, $chr, $start, $end, $strand;
	
	#---> RETURN VALUE <---#
	if ( ! $end > $start ) {
		return 0;
	}
	else {
		return \@return;
	}
}
#-----------------------#
sub print_output {
### Function: Generates the output files
### Accepts: 1. Reference to region type hash of hashes generated by subroutine 'ENSEMBL_gtf_AoH_to_region_type_HoH_ref'
### Returns: n/a
### Dependencies: n/a
### Type: Specific
	#---> PASS ARGUMENTS <---#
	my $region_type_HoH_ref = shift;

	#---> BODY <---#
	## Call generic print subroutine for each region type
	&print_hash_inverted_keys_alphanum_sorted("${prefix}.motevo", $$region_type_HoH_ref{"main"});
	&print_hash_inverted_keys_alphanum_sorted("${prefix}_${up}_up.motevo", $$region_type_HoH_ref{"up"}) if $up && $sep;
	&print_hash_inverted_keys_alphanum_sorted("${prefix}_${down}_down.motevo", $$region_type_HoH_ref{"down"}) if $down && $sep;
	&print_hash_inverted_keys_alphanum_sorted("${prefix}_${left}_left.motevo", $$region_type_HoH_ref{"left"}) if $left && $sep;
	&print_hash_inverted_keys_alphanum_sorted("${prefix}_${right}_right.motevo", $$region_type_HoH_ref{"right"}) if $right && $sep;
}
#-----------------------#
sub print_hash_inverted_keys_alphanum_sorted {
### Function: Writes value -> key hash pairs to output file (one pair per line)
### Accepts: 1. Output file; 2. Hash reference; 3. Separator [DEFAULT: TAB]
### Returns: n/a
### Dependencies: n/a
### Type: Generic
	#---> PASS ARGUMENTS <---#
	my $out_file = shift;
	my $hash_ref = shift;
	my $sep = defined($_[0]) ? $_ : "\t";

	#---> STATUS MESSAGE <---#	
	print STDOUT "Writing output to file '$out_file'...\n" unless $quiet;

	#---> BODY <---#
	# Open GTF filehandle
	open OUT, ">" . $out_file;
	## Iterate over each key value pair of 'ID -> transcript' hash of array references
	foreach my $key (sort { $a cmp $b } keys %$hash_ref) {
		# Print hash key TAB hash value
		print OUT $$hash_ref{$key} . $sep . $key . "\n";
	}
	# Close OUT filehandle
	close OUT;
	
	#---> STATUS MESSAGE <---#
	print STDOUT "File '$out_file' written.\n\n" unless $quiet;
}
#=======================#
#	SUBROUTINES END	#
#=======================#
