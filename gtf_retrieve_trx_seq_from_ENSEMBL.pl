#!/usr/bin/perl

#==================#
#   HEADER START   #
#==================#
### Name: gtf_retrieve_trx_seq_from_ENSEMBL.pl
### Created: Jun 17, 2013
### Modified: Jun 17, 2013
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
### Version: v0.1
### Adapted from: n/a
### Requirements: ENSEMBL Perl API
#==================#
### Description: Retrieves sequences of transcripts by 'transcript_id' meta information from a GTF file
### Output: FASTA file with 'transcript ID' in identifier line
### Usage: perl ./gtf_retrieve_trx_seq_from_ENSEMBL.pl for information on required and optional arguments
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
use lib '/usr/lib/perl/modules/bioperl-1.2.3';
use lib '/usr/lib/perl/modules/ensembl/modules';
use Bio::EnsEMBL::Registry;

#---> USAGE <---#
my $usage_info = &usage;

#---> OPTIONS / ARGUMENTS <---#
my $usage = '';
my $quiet = '';
my $gtf = '';
my $out_file = '';
my $organism = 'Human';
my $options_result = GetOptions (
	'usage|help' => \$usage,
	'quiet|no_verbose' => \$quiet,
	'organism=s' => $organism,
	#-----------------------#
	'gtf=s' => \$gtf,
	'out_file=s' => \$out_file
);
die $usage_info if $usage || !$options_result;
die $usage_info if !$gtf || !$out_file; 
#==========================#
#    PRE-REQUISITES END    #
#==========================#


#================#
#   MAIN START   #
#================#
#---> STATUS MESSAGE <---#
print "Starting 'gtf_retrieve_trx_seq_from_ENSEMBL.pl'...\n\n" unless $quiet;

#---> MAIN VARIABLES <---#
my $ids;

#---> BODY <---#
$ids = &gtfTrxIDToArrayRef($gtf, "transcript");
&writeEnsemblTrxSeq($out_file, $ids, $organism);

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
'Usage: perl ./gtf_retrieve_trx_seq_from_ENSEMBL.pl [OPTIONS] --gtf [FILE] --out_file [FILE]
==================================================
Required arguments:
--gtf			GTF input file [FILE]
--out_file		Output filename [FILE]
==================================================
Optional arguments:
--organism	Organism [DEFAULT: "Human"]
--usage|help		Show this information and die
--quiet			Shut up!
';
}
#-----------------------#
sub gtfTrxIDToArrayRef {
### Function: Retrieves a list of unique 'transcript_id' identifiers from a GTF file 
### Accepts: 1. Input GTF file; 2. Type of entries to select from GTF file [DEFAULT: "transcript"]
### Returns: Reference to array of identifiers
### Dependencies: n/a
### Type: Generic
	#---> PASS ARGUMENTS ---#
	my $gtf_file = shift;
	my $type = defined($_[0]) ? shift : "transcript";

	#---> STATUS MESSAGE <---#
	print STDOUT "Processing file '$gtf_file' to extract transcript identifiers..." . "\n" unless $quiet;	
	
	#---> SUBROUTINE VARIABLES <---#
	my %transcripts;
	my @transcripts;
	
	#---> BODY <---#
	# Open GTF filehandle
	open GTF, $gtf_file;	
	## Iterate over GTF file line by line
	while (<GTF>) {
		# Assign line to dedicated variable
		my $line = $_;		
		# Skip line IF line starts with '#' (header line)
		next if $line =~ /\A#/;
		# Split line by tabs into array @gtf_line
		my @gtf_line = split /\t/, $line;
		# Skip line IF 'type field' does not match value in $type
		next unless $gtf_line[2] eq $type; 
		# Split 'metadata field' by semicolon into array @gtf_meta
		my @gtf_meta = split /;/, $gtf_line[8];
		# Initialize hash %gtf_meta;
		my %gtf_meta;
		## Iterate over elements of array @gtf_meta
		foreach my $meta (@gtf_meta) {
			## IF element matches pattern: zero or more whitespace + [CAPTURE 1: word characters] + one or more whitespace + optional '"' + [CAPTURE2: word characters] + optional '.' + optional word characters + optional '"'
			if ($meta =~ /\s*(\w*)\s+"?(\w*)\.?\w*"?/) {
				# Add capture groups to hash %gtf_meta (capture group 1 = key; capture group 2 = value)
				$gtf_meta{$1} = $2;	
			}
			## ELSE...
			else {
				# ...issue warning
				print STDERR "[WARNING] Meta information entry '$meta' skipped because it does not match the pre-defined mask. Please verify the syntax/format of the following line of your input data:\n$line\n"; 
			}
		}
		## IF 'transcript_id' exist as keys in hash %gtf_meta
		if (exists($gtf_meta{"gene_id"}) && exists($gtf_meta{"gene_name"}) && exists($gtf_meta{"transcript_id"})) {
			# Add to hash %transcripts as key (value empty; this is to avoid duplicates)
			$transcripts{$gtf_meta{"transcript_id"}} = "";
		}
		## ELSE...
		else {
			# ...issue warning
			print STDERR "[WARNING] Field 'transcript_id' not available. The following line was skipped:\n$line\n"; 
		}
	}
	# Extract unique transcript IDs (hash keys) and sort alphanumerically
	@transcripts = sort { $a cmp $b } keys %transcripts;
	# Close GTF filehandle
	close GTF;	
	
	#---> STATUS MESSAGE <---#
	print STDOUT "Transcript identifiers extracted." . "\n\n" unless $quiet;	
	
	#---> RETURN VALUE <---#
	return \@transcripts;
}
#-----------------------#
sub writeEnsemblTrxSeq {
### Function: Retrieves transcript sequence by 'transcript_id' from ENSEMBL and writes them to FASTA file
### Accepts: 1. Filename for FASTA output file [FILE]; 2. Reference to array of transcript identifiers [ARRAY REF]; 3. Organism [DEFAULT: "Human"]
### Returns: n/a
### Dependencies: n/a
### Type: Generic
	#---> PASS ARGUMENTS ---#
	my $out_file = shift;
	my $ids = shift;
	my $org = defined($_[0]) ? shift : "Human";	
	
	#---> STATUS MESSAGE <---#
	print STDOUT "Retrieving " . @$ids . " transcript sequences from ENSEMBL..." . "\n" unless $quiet;

	#---> SUBROUTINE VARIABLES <---#
	my $registry = 'Bio::EnsEMBL::Registry';
	my $index = 0;	
	
	#---> BODY <---#
	# Connect to Ensembl Registry
	$registry->load_registry_from_db(
	    -host => 'ensembldb.ensembl.org', # alternatively 'useastdb.ensembl.org'
	    -user => 'anonymous'
	);	

	# Set up adaptor to human transcript database
	my $trx_adaptor = $registry->get_adaptor( $org, 'Core', 'Transcript' );

	# Open OUT filehandle
	open OUT, ">" . $out_file;

	## Iterate through transcript identifiers in $ids array reference
	foreach my $id (@$ids) {
		# Save a copy of the identifier
		my $id_fetch = $id;
		# Replace ENSTR with ENST0 (this is done because ENSTR transcripts are not found!)
		$id_fetch =~ s/ENSTR/ENST0/;
		# Fetch transcript by identifier
		my $trx = $trx_adaptor->fetch_by_stable_id($id_fetch);
		# Fetch transcript sequence
		my $trx_seq = $trx->spliced_seq();
		# Increase index;
		$index++;
		# Print identifier and sequence in FASTA format
		print OUT ">" . $id . "\n" . $trx_seq . "\n";
		## IF index reaches the next thousand 1000
		if ($index % 1000 == 0) {
			# Print status message
			print STDOUT "$index sequences retrieved..." . "\n" unless $quiet;
		}
	}

	# Close OUT filehandle
	close OUT;
	
	#---> STATUS MESSAGE <---#
	print STDOUT "Sequences retrieved and written to FASTA file '$out_file'." . "\n\n" unless $quiet;
	
	#---> RETURN VALUE <---#
	return ;
}
#=======================#
#    SUBROUTINES END    #
#=======================#