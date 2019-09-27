#!/usr/bin/perl

#    FASTX-toolkit - FASTA/FASTQ preprocessing tools.
#    Copyright (C) 2009  A. Gordon (gordon@cshl.edu)
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU Affero General Public License as
#   published by the Free Software Foundation, either version 3 of the
#   License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU Affero General Public License for more details.
#
#    You should have received a copy of the GNU Affero General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

use strict;
use warnings;
use File::Basename;
use IO::Handle;
use Data::Dumper;
use Getopt::Long;
use Carp;

## This program splits a FASTQ/FASTA file into several smaller files,
## Based on barcode matching.
##
## run with "--help" for usage information
##
## Assaf Gordon <gordon@cshl.edu> , 11sep2008

#------------------------------------------------#
#   Adapted by Alexander Kanitz on 10-JUL-2013   #
#------------------------------------------------#

# Forward declarations
sub load_barcodes ($);
sub parse_command_line ;
sub match_sequences ;
sub mismatch_count($$) ;
sub print_results;
sub open_and_detect_input_format;
sub read_record;
sub write_record($$);
sub usage();

# Global flags and arguments, 
# Set by command line argumens
my $in_file;
my $out_dir;
my $out_dir_unmatched = "";
my $barcodes;
my $barcodes_at_bol = 1;
my $trim = 0;
my $allow_partial_overlap = 0;
my $allowed_mismatches;
my $suffix = "barcode";
my $quiet = 0;
my $debug = 0;
my $fastq_format = 1;

# Global variables 
# Populated by 'create_output_files' and other
my %filenames;
my %files;
my $basename;
my $extension;
my %counts = ( 'unmatched' => 0 );
my $barcodes_length;
my @barcodes;
my $input_file_io;
my $fh;

# The Four lines per record in FASTQ format.
# (when using FASTA format, only the first two are used)
my $seq_name;
my $seq_bases;
my $seq_name2;
my $seq_qualities;



#
# Start of Program
#
parse_command_line;

load_barcodes($barcodes);

open_and_detect_input_format;

match_sequences;

print_results unless $quiet;
#
# End of program
#



#
# Parse command line
#
sub parse_command_line {
	my $help;

	usage() if (scalar @ARGV==0);

	my $result = GetOptions ( 
				  "infile=s" => \$in_file,
				  "bc=s" => \$barcodes,
				  "outdir=s" => \$out_dir,
				  "outdir_um=s" => \$out_dir_unmatched,
				  "suffix=s" => \$suffix,
				  "eol"  => sub { $barcodes_at_bol = 0 },
				  "trim" => \$trim,
				  "quiet" => \$quiet, 
				  "partial=i" => \$allow_partial_overlap,
				  "debug" => \$debug,
				  "mismatches=i" => \$allowed_mismatches,
				  "help" => \$help
				  ) ;

	usage() if ($help);

	die "[ERROR] Input FASTA/FASTQ file not specified (use '--infile [FILENAME]')!\n" unless defined $in_file;
	die "[ERROR] No barcodes specified (use '--bc [BARCODE_1|BARCODE_2|...|BARCODE_N]')!\n" unless defined $barcodes;	
	die "[ERROR] Output folder not specified (use '--outdir [PATH]')!\n" unless defined $out_dir;
	die "[ERROR] Bad suffix value '$suffix' (must be alphanumeric or underscore)!\n" unless $suffix =~ m/^[\w_]+$/;
	die "[ERROR] Option '--partial' allows only positive integers and 0!\n" if $allow_partial_overlap < 0;
	die "[ERROR] Option '--mismatches' allows only positive integers and 0!\n" if defined $allowed_mismatches && $allowed_mismatches < 0;
	die "[ERROR] Value for '--partial' ($allow_partial_overlap) cannot be greater than value for --mismatches ($allowed_mismatches)!\n" if defined $allowed_mismatches && $allow_partial_overlap > $allowed_mismatches;

    $out_dir_unmatched = $out_dir if $out_dir_unmatched eq "";

	$out_dir =~ s | ^(.)*\/$ | $1 |;
	$out_dir_unmatched =~ s | ^(.)*\/$ | $1 |; 
	
	exit unless $result;
}

#
# Read the barcode
#
sub load_barcodes ($) {

	# Pass parameter
	my $bc_string = shift;
	
	# Split barcode string by pipe ("|") character
	my @bc_seqs = split /\|/, $bc_string;

	# Set barcode counter
	my $count = 1;

	# Traverse over each barcode sequence
	foreach my $barcode (@bc_seqs) {

		# Characters in barcode to uppercase
		$barcode = uc($barcode);

		# Set identifier via specified suffix and current barcode count
		my $ident = $barcode;

		# Set barcode length from first instance (different lengths are not allowed, see below)
		$barcodes_length = length($barcode) unless defined $barcodes_length;
		
		# Set default allowed mismatches (one for every three nucleotides) unless manually specified
		$allowed_mismatches = int($barcodes_length / 3) unless defined $allowed_mismatches;

		# Sanity checks on the barcodes
		die "[ERROR] Bad barcode value '$barcode'. Only DNA alphabet [ACGT] allowed.\n" unless $barcode =~ m/^[AGCT]+$/;
		die "[ERROR] Found barcodes of different lengths (not supported yet).\n" unless $barcodes_length == length($barcode);
		die "[ERROR] Barcode '$barcode' is shorter or equal to allowed mismatches ($allowed_mismatches). Specify fewer mismatches.\n" if length($barcode) <= $allowed_mismatches;

	 	# Add barcode and ID to dedicated array
	 	push @barcodes, [$ident, $barcode];
	 	
		# Add modified versions of barcode if partial overlaps are allowed
		if ($allow_partial_overlap > 0) {
			foreach my $i (1 .. $allow_partial_overlap) {
				substr $barcode, ($barcodes_at_bol) ? 0 : -$i, $i, '';
	 			push @barcodes, [$ident, $barcode];
			}
		}
		
		# Increase barcode count
		$count++;

	}
	
	# DEBUGGING
	if ($debug) {
		print STDERR "barcode\tsequence\n";
		foreach my $barcoderef (@barcodes) {
			my ($ident, $seq) = @{$barcoderef};
			print STDERR $ident,"\t", $seq ,"\n";
		}
	}
}

# Create one output file for each barcode.
# (Also create a file for the dummy 'unmatched' barcode)
sub create_output_files {
	my %barcodes = map { $_->[0] => 1 } @barcodes; #generate a uniq list of barcode identifiers;
	$barcodes{"unmatched"} = 1;

	# Set basename
	($basename = basename($in_file)) =~ s/\.[^.]+$//;
	($extension = $in_file) =~ s/.*\.([^.]+$)/$1/;

	# Generate output filenames and file handles
	foreach my $ident (keys %barcodes) {
		my $new_filename;
		if ($ident eq "unmatched") {
			$new_filename = $out_dir_unmatched . "/" . $basename . "_" . $suffix . "_" . $ident . "." . $extension;
		} 
		else {
			$new_filename = $out_dir . "/" . $basename . "_" . $suffix . "_" . $ident . "." . $extension;
		} 
		$filenames{$ident} = $new_filename;
		open my $file, ">$new_filename" or die "[ERROR] Failed to create output file '$new_filename'!\n"; 
		$files{$ident} = $file;
	}
}

sub match_sequences {
	my %barcodes = map { $_->[0] => 1 } @barcodes; #generate a uniq list of barcode identifiers;
	$barcodes{"unmatched"} = 1;

	#reset counters
	foreach my $ident ( keys %barcodes ) {
		$counts{$ident} = 0;
	}

	create_output_files;

	# Read file FASTQ file
	# split according to barcodes
	while ( read_record ) {
		chomp $seq_bases;

		print STDERR "sequence $seq_bases: \n" if $debug;

		my $best_barcode_mismatches_count = $barcodes_length;
		my $best_barcode_ident = undef;

		#Try all barcodes, find the one with the lowest mismatch count
		foreach my $barcoderef (@barcodes) {
			my ($ident, $barcode) = @{$barcoderef};

			# Get DNA fragment (in the length of the barcodes)
			# The barcode will be tested only against this fragment
			# (no point in testing the barcode against the whole sequence)
			my $sequence_fragment;
			if ($barcodes_at_bol) {
				$sequence_fragment = substr $seq_bases, 0, $barcodes_length;
			} else {
				$sequence_fragment = substr $seq_bases, - $barcodes_length;
			}

			my $mm = mismatch_count($sequence_fragment, $barcode) ; 

			# if this is a partial match, add the non-overlap as a mismatch
			# (partial barcodes are shorter than the length of the original barcodes)
			$mm += ($barcodes_length - length($barcode)); 

			if ( $mm < $best_barcode_mismatches_count ) {
				$best_barcode_mismatches_count = $mm ;
				$best_barcode_ident = $ident ;
			}
		}

		$best_barcode_ident = 'unmatched' 
			if ( (!defined $best_barcode_ident) || $best_barcode_mismatches_count>$allowed_mismatches) ;

		print STDERR "sequence $seq_bases matched barcode: $best_barcode_ident\n" if $debug;

		$counts{$best_barcode_ident}++;

		#get the file associated with the matched barcode.
		#(note: there's also a file associated with 'unmatched' barcode)
		my $file = $files{$best_barcode_ident};

		write_record($file, $best_barcode_ident);
	}
	
	close $fh;
	
}

#Quickly calculate hamming distance between two strings
#
#NOTE: Strings must be same length.
#      returns number of different characters.
#see  http://www.perlmonks.org/?node_id=500235
sub mismatch_count($$) { length( $_[ 0 ] ) - ( ( $_[ 0 ] ^ $_[ 1 ] ) =~ tr[\0][\0] ) }



sub print_results
{
	print "Barcode\tCount\tLocation\n";
	my $total = 0 ;
	foreach my $ident (sort keys %counts) {
		print $ident, "\t", $counts{$ident},"\t",$filenames{$ident},"\n";
		$total += $counts{$ident};
	}
	print "total\t",$total,"\n";
}


sub read_record
{
	$seq_name = $input_file_io->getline();

	return undef unless defined $seq_name; # End of file?

	$seq_bases = $input_file_io->getline();
	die "Error: bad input file, expecting line with sequences\n" unless defined $seq_bases;

	# If using FASTQ format, read two more lines
	if ($fastq_format) {
		$seq_name2  = $input_file_io->getline();
		die "Error: bad input file, expecting line with sequence name2\n" unless defined $seq_name2;

		$seq_qualities = $input_file_io->getline();
		die "Error: bad input file, expecting line with quality scores\n" unless defined $seq_qualities;
		chomp($seq_qualities);
	}
	return 1;
}

sub write_record($$)
{
	# Pass parameter
	my $file = shift;
	my $ident = shift;

	# Catch file handle problem
	croak "Bad file handle" unless defined $file;

	# Trim sequences
	if ($trim && $ident ne "unmatched" && length $seq_bases > $barcodes_length) {
		if ($barcodes_at_bol) {
			$seq_bases = substr $seq_bases, $barcodes_length if length $seq_bases > $barcodes_length;
		} else {
			$seq_bases = substr $seq_bases, -0, -$barcodes_length if length $seq_bases > $barcodes_length;
		}
	}
	
	# Print first two lines (FASTA & FASTQ)
	print $file $seq_name;
	print $file $seq_bases,"\n";

	# If FASTQ...
	if ($fastq_format) {
	
		# Trim quality scores
		if ($trim && $ident ne "unmatched" && length $seq_qualities > $barcodes_length) {
			if ($barcodes_at_bol) {
				$seq_qualities = substr $seq_qualities, $barcodes_length;
			} else {
				$seq_qualities = substr $seq_qualities, -0, -$barcodes_length;
			}
		}
		
		# Print two more lines (FASTQ only!)
		print $file $seq_name2;
		print $file $seq_qualities,"\n";
	}
}

sub open_and_detect_input_format
{
	$input_file_io  = new IO::Handle;
	open $fh, $in_file or die "$0: open: $!";;
	die "Failed to open STDIN " unless $input_file_io->fdopen(fileno($fh),"r");

	# Get the first characeter, and push it back
	my $first_char = $input_file_io->getc();
	$input_file_io->ungetc(ord $first_char);

	if ($first_char eq '>') {
		# FASTA format
		$fastq_format = 0 ;
		print STDERR "Detected FASTA format\n" if $debug;
	} elsif ($first_char eq '@') {
		# FASTQ format
		$fastq_format = 1;
		print STDERR "Detected FASTQ format\n" if $debug;
	} else {
		die "Error: unknown file format. First character = '$first_char' (expecting > or \@)\n";
	}
}

sub usage()
{
print<<EOF;
Barcode Splitter, by Assaf Gordon (gordon\@cshl.edu), 11sep2008;
Modified by Alexander Kanitz, 10-JUL-13 

This program splits a FASTA or FASTQ (format auto-detected) file into several files 
based on barcode matching. Output files will be writen to disk, a summary printed 
to the standard output.

usage: $0 --infile FILE --bc STRING --outdir PATH [--suffix STRING] [--eol] [--trim]
         [--mismatches POS. INT] [--partial POS. INT] [--quiet] [--debug] [--help]

Arguments:

--infile [FILE]		Required: FASTA or FASTQ input file
--bc [STRING]		Required: String of pipe-separated ("|") barcodes (details below)
--outdir [PATH]		Required: Directory in which output files are saved
--outdir_um [PATH] 	Indicate a different output path for unmatched sequences
--suffix [STRING]	Used to distinguish output from input file (details below)
--eol			Match barcodes at the 3'-end of sequences [DEFAULT: match at 5'-end]
--trim			Remove barcode sequences when writing output files (details below)
--mismatches [POS. INT]	Maximum number of allowed mismatches [DEFAULT: details below]
--partial [POS. INT]	Allow partial barcode overlap (explanation below) [DEFAULT: 0]
--quiet			Don't print count summary after finishing [DEFAULT: print summary]
--debug			Print debug information to STDERR
--help			This helpful help screen

Example (assuming 'test.fq' is a FASTQ file):

	perl ./barcode_split.pl --infile test.fq --bc 'AC|GA|CT' --outdir /path/to/output


Barcode format (--bc)
--------------
Barcodes need to be of equal length and are specified as a string, with each 
entry being separated by the pipe ("|") character. Place the string between 
(single or double) quotation marks to avoid misinterpretation of the pipe 
character, like so:

	'ACG|GAT'
	'ACG|GAT|CTA'

Single barcodes are supported and may be used to trim e.g. leading nucleotides 
that were added during sequencing library preparation. The program is, however, 
not a robust tool for adapter removal!


Writing output (--suffix; --trim)
--------------
For each barcode, a new FASTA or FASTQ file will be created in the indicated 
output directory. The file type, basename and extension of the input file 
will be kept, but between basename and extension, the --suffix and the barcode 
in question will be inserted (separated by underscores). Anothe file will be 
created that contains in its filename the string "unmatched" instead of the 
barcode. This file will contain all the sequences that could not be matched 
to any of the specified barcodes.

The default value for --suffix is "barcode", but can be replaced by any 
alphanumeric string (and underscore), by indicating '--suffix mystring' when 
calling the program.

Running the example above (with the default --suffix value) will create the 
following files:

	/path/to/output/example_barcode_AC.fq
	/path/to/output/example_barcode_GA.fq
	/path/to/output/example_barcode_CT.fg
	/path/to/output/example_barcode_unmatched.fq
	
If --trim is specified, the number of nucleotides equal to the barcode length 
are removed from the 5'-end (default) or 3'-end (with --eol) when written to
the corresponding output file (this also applies to the quality control string 
in a FASTQ file). Specify this option if no further analysis of the barcodes 
is to be expected after splitting to e.g. prevent unnecessary mismatches in 
downstream alignment applications. The use of --trim and --partial is not 
supported as only a fixed number of nucleotides is resected, thus potentially 
leading to the resection of '--partial N' non-barcode nucleotides. Use at your 
own risk!


Allowed mismatches (--mismatches)
------------------
The default value for mismatches depends on the length of the supplied 
barcodes. Starting with 0, it gets increased by one for every three nucleotides,
i.e. for a barcode of length 1-2 '--mismatches' is set to 0, for length 3-5 it 
is set to 1, for length 6-8 to 2 mismatches and so on. Overwriting by specifying 
'--mismatches POS. INT' will set the allowed mismatches to a fixed number.

For exact matches, regardless of barcode length, set '--mismatches 0'.


Barcode matching (--eol; --mismatches; --partial)
----------------

** Without partial matching:

Count mismatches between the FASTA/Q sequences and the barcodes.
The barcode which matched with the lowest mismatches count (providing the
count is small or equal to '--mismatches POS. INT') 'gets' the sequences.

Example:

Input Sequence:
    GATTTACTATGTAAAGATAGAAGGAATAAGGTGAAG

Barcodes:
	'GATCT|ATCGT|GTGAT|TGTCT'
	BC1	GATCT
    BC2	ATCGT
    BC3	GTGAT
    BC4	TGTCT

Matching with '--mismatches 1':
	GATTTACTATGTAAAGATAGAAGGAATAAGGTGAAG
	GATCT (1 mismatch, BC1)
	ATCGT (4 mismatches, BC2)
	GTGAT (3 mismatches, BC3)
	TGTCT (3 mismatches, BC4)

This sequence will be classified as 'BC1' (it has the lowest mismatch count).
If '--mismatches 0' was specified, this sequence would be classified as 
'unmatched' (because, although BC1 had the lowest mismatch count, it is above 
the maximum allowed mismatches).

Matching with '--eol' (end of line) does the same, but from the other side
of the sequence.

** With partial matching (very similar to indels):

Same as above, with the following addition: barcodes are also checked for
partial overlap (number allowed non-overlapping bases: '--partial POS. INT').

Example:

Input sequence (as above except for missing 'G' at the 5'-end):
	ATTTACTATGTAAAGATAGAAGGAATAAGGTGAAG
	
Barcode:
	BC1 from above

Matching WITHOUT partial overlapping:
	ATTTACTATGTAAAGATAGAAGGAATAAGGTGAAG
	GATCT (4 mismatches)

Matching WITH partial overlapping will also try this match:
	-ATTTACTATGTAAAGATAGAAGGAATAAGGTGAAG
	GATCT (1 mismatch)

Notes: 

Scoring counts a missing base as a mismatch, so the final
mismatch count is 2 (1 'real' mismatch, 1 'missing base' mismatch).
If running with '--mismatches 2' (meaning allowing up to 2 mismatches), this 
seqeunce will be classified as BC1.

The use of --partial with --trim is not supported. Use at your own risk!
See above for details.

EOF

exit 1;
}

### ISSUES
# Solve more elegantly what will happen if the sequence consists of the barcode only! (particularly with --trim)
# Partial overlap with trim not really supported
