#!/usr/bin/perl

#==================#
#   HEADER START   #
#==================#
### Name: generate_GENCODE_gene_transcript_lookup_table.pl
### Created: Jul 31, 2013
### Modified: Jul 31, 2013
### Author: Alexander Kanitz, Stamatina Stamelaki
### Company: Zavolan Group, Biozentrum, University of Basel
### Version: v1.0
### Adapted from: n/a
### Requirements: Getopt::Long
#==================#
### Description: Builds Exon > Transcript lookup table; it considers only those features/lines of the GTF file that are annotated as exons
### Output: Returns a tab-separated file of the following format: Exon_name (Ensembl gene ID + genomic coordinates, separated by underscores: ENSG_X_chr_strand_start_end) TAB ENSG_X@ENST_1|ENSG_X@ENST_2|ENSG_X@ENST_N (where ENSG_X is the Ensembl gene ID and ENST_N are the Ensembl transcript IDs)
### Output example line: ENSG00000223972_chr1_+_12345_12456	\t ENSG00000223972@ENST00000456328|ENSG00000223972@ENST00000515242|ENSG00000223972@ENST00000518655
### Usage: perl ./perl generate_GENCODE_exon_transcript_lookup_table.pl for information on required and optional arguments
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
my ($exons_AoH_ref, $exon_trx_hash_ref);

#---> BODY <---#
$exons_AoH_ref = gencode_gtf_to_AoH($gtf, "exon");
$exon_trx_hash_ref = gencode_gtf_AoH_to_exon_trx_hash_ref($exons_AoH_ref);
&print_hash_keys_alphanum_sorted($out_file, $exon_trx_hash_ref);

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
'Usage: perl ./generate_GENCODE_exon_transcript_lookup_table.pl [OPTIONS] --gtf [FILE] --out_file [FILE]
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
sub gencode_gtf_AoH_to_exon_trx_hash_ref {
### Function: Generates 'exon_id' -> transcripts hash from GENCODE GTF array of hashes
### Accepts: 1. Reference to GENCODE GTF array of hashes
### Returns: Reference to exon -> transcripts hash
### Dependencies: n/a
### Type: Specialized
	#---> PASS ARGUMENTS ---#
	my @exons_AoH = @{shift()};
		
	#---> STATUS MESSAGE <---#
	print STDOUT "Building exon -> transcripts hash..." . "\n" unless $quiet;

	#---> SUBROUTINE VARIABLES <---#
	my %exon_trx_hash;
	
	#---> BODY <---#
	## Traverse through each element $gtf_line of array @exons_AoH
	foreach my $gtf_line ( @exons_AoH ) {
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
	    ## Build exon ID
	    my $exon_id = join "_", $gene_id, $$gtf_line{"chr"}, $$gtf_line{"strand"}, $$gtf_line{"start"}, $$gtf_line{"end"};
	    ## Build unique transcript ID
	    my $uniq_trx_id = join "@", $gene_id, $transcript_id; 
		## Add key value pair to hash; if key does not yet exist
        if (! defined $exon_trx_hash{$exon_id}) {
                # Add pair
                $exon_trx_hash{$exon_id} = $uniq_trx_id;
        } else {
        		# Else append value to existing value
                $exon_trx_hash{$exon_id} .= "|" . $uniq_trx_id;
        }		
	}
	
	#---> STATUS MESSAGE <---#
	print STDOUT "Exon -> transcripts hash built." . "\n\n" unless $quiet;
	
	#---> RETURN VALUE <---#
	return \%exon_trx_hash;
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
