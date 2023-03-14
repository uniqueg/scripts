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
### Description: From a BED or SAM file, find the mates of a list of read IDs. Naming convention: read_id\1 and read_id\2, for first and second mates respectively.
### Output: All lines of the input BED or SAM file that contain the mates of reads with the specified IDs.
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
my $read_file = '';
my $id_file = '';
my $out_file = '';
my $sam = 0;
my $head = 0;
my $options_result = GetOptions (
	'usage|help' => \$usage,
	'quiet' => \$quiet,
	'sam' => \$sam,
	'print_header' => \$head,
	#-----------------------#
	'read_file=s' => \$read_file,
	'id_file=s' => \$id_file,
	'out_file=s' => \$out_file
);
die $usage_info if $usage || !$options_result;
die $usage_info if !$read_file || !$id_file || !$out_file; 
$head = 0 if $sam == 0;

#==========================#
#    PRE-REQUISITES END    #
#==========================#


#================#
#   MAIN START   #
#================#
#---> STATUS MESSAGE <---#
print "Starting '$0'...\n" unless $quiet;

#---> MAIN VARIABLES <---#
my $id_hash_ref;

#---> BODY <---#
$id_hash_ref = &lines_to_hash($id_file);
&reverse_mate_id($id_hash_ref);
&print_BED_or_SAM_lines_with_specific_IDs($read_file, $id_hash_ref, $out_file, $sam, $head);

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

Description: From a BED or SAM file, find the mates of a list of read IDs. Naming convention: read_id\1 and read_id\2, for first and second mates respectively.

==================================================
Required arguments:
--read_file	Read file in BED (default) or SAM (with --sam switch) format
--id_file	Flat text file with read IDs (one per line)
--out_file	Output file in same format as input file
==================================================
Optional arguments:
--sam	Read file has SAM format (default: BED)
--print_header	Print the header for the SAM file (only if --sam is set)
--usage|help	Show this information and die
--quiet	Shut up!
';
}
#-----------------------#
sub lines_to_hash {
### Function: Read flat file into hash (each line being a key; no values)
### Accepts: 1. Filename
### Returns: Reference to hash
### Dependencies: n/a
### Type: Generic
	#---> PASS ARGUMENTS ---#
	my $in_file = shift;	
	
	#---> STATUS MESSAGE <---#
	print STDOUT "Reading file '$in_file'..." . "\n" unless $quiet;

	#---> SUBROUTINE VARIABLES <---#
	my %hash;
	
	#---> BODY <---#
	
		# Open file handle FA
		open FILE, "$in_file";	
		## Traverse through file handle FA line by line
		while (<FILE>) {
			# Remove trailing newline character "\n"
			chomp;
			## Check whether current line is identifier line
			$hash{$_} = undef;				
		}
		# Close file handle FA
		close FILE;
		
	#---> STATUS MESSAGE <---#
	print "File $in_file read." . "\n" unless $quiet;	
	
	#---> RETURN VALUE <---#
	return \%hash;
}
#-----------------------#
sub reverse_mate_id {
### Function: Reverses mate number in read IDs (in place substitution!)
### Accepts: Reference to hash of read IDs containing terminal mate IDs of the form "/1" and "/2"
### Returns: n/a
### Dependencies: n/a
### Type: Generic
	#---> PASS ARGUMENTS ---#
	my $read_id_hash_ref = shift;
	
	#---> STATUS MESSAGE <---#
	print STDOUT "Reversing mate identifiers." . "\n" unless $quiet;

	#---> BODY <---#
		
		## Traverse through each key $id of hash %$read_id_hash_ref
		foreach my $id ( keys %$read_id_hash_ref ) {
		    # Save hash value
		   	my $value = $read_id_hash_ref->{$id} if defined $read_id_hash_ref->{$id};
		    # Delete hash element
		    delete $read_id_hash_ref->{$id};
		    ## In hash key, replace "/1" with "/2" (and vice verse) and re-add the hash key to the hash
			if ($id =~ /\/1$/) {
				($id = $id) =~ s/\/1$/\/2/;
				$read_id_hash_ref->{$id} = defined ($value) ? $value : undef;
			}
			else {
				($id = $id) =~ s/\/2$/\/1/;
				$read_id_hash_ref->{$id} = defined ($value) ? $value : undef;
			}
		}
	
	#---> STATUS MESSAGE <---#
	print STDOUT "Mate identifiers reversed." . "\n" unless $quiet;
}
#-----------------------#
sub print_BED_or_SAM_lines_with_specific_IDs {
### Function: function
### Accepts: 1. Input file in SAM or BED format; 2. Reference to hash containing read IDs as keys (empty strings as values); 3. Output file; 4. Flag whether input is in SAM format; 5. Flag whether header should be printed (only valid/relevant for SAM files)
### Returns: n/a
### Dependencies: n/a
### Type: Generic
	#---> PASS ARGUMENTS ---#
	my ($read_file, $id_hash_ref, $out_file, $sam, $head) = @_;
	
	#---> STATUS MESSAGE <---#
	print STDOUT "Crawling through reads file '$read_file'..." . "\n" unless $quiet;

	#---> SUBROUTINE VARIABLES <---#
	
	#---> BODY <---#
	
		## Open file handles
		open READS, "<", $read_file or die "[ERROR] Could not open file '$read_file'!\n";
		open OUT, ">", $out_file or die "[ERROR] Could not open file '$out_file'!\n";
		
		#---> If $head is set ($sam only!), print header to OUT; once non-header line is met, exit <---# 
		{
			if ($head) {
	
				## Set block variables
				my $regex_header = '^\@[A-Za-z][A-Za-z](\t[A-Za-z][A-Za-z0-9]:[ -~]+)+$';
				my $regex_comment = '/^\@CO\t.*/';
				my $last_line;
				
				## Traverse through file handle READS line by line
				while (<READS>) {
					## If line is a header line, print to filehandle HEAD (if $head is set)
					if ( /$regex_header/ || /$regex_comment/ ) {
						print OUT if $head;
					}
					## Else exit loop
					else {
						$last_line = $_;
						last;
					}
				}
	
				## Set the correct filehandle position...
				{
					use bytes;
					seek READS, -length($last_line), 1;			
				}
	
			}
		}
	
		#---> Print line of READS to OUT if its identifier is present in hash <---#
		while (<READS>) {
			my @line = split "\t";
			my $id = $sam ? $line[0] : $line[3];
			print OUT if exists $$id_hash_ref{$id};
		}
	
		## Close file handles
		close OUT;
		close READS;
	
	#---> STATUS MESSAGE <---#
	print STDOUT "Read file processed. Output written to '$out_file'." . "\n" unless $quiet;
}
#=======================#
#    SUBROUTINES END    #
#=======================#