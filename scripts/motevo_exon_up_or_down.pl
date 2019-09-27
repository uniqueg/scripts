#!/usr/bin/perl

#==================#
#   HEADER START   #
#==================#
### Name: motevo_exon_up_or_down.pl
### Created: Nov 28, 2013
### Modified: Nov 28, 2013
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
### Version: v0.1
### Adapted from: n/a
### Requirements: 
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
my $motevo = '';
my $exon = '';
my $up = '';
my $down = '';
my $options_result = GetOptions (
	'h|usage|help' => \$usage,
	'q|quiet' => \$quiet,
	#-----------------------#
	'm|motevo=s' => \$motevo,
	'e|exon=s' => \$exon,
	'u|up=s' => \$up,
	'd|down=s' => \$down
);
die $usage_info if $usage || !$options_result;
die $usage_info if !$motevo || !$exon || !$up || !$down; 

#---> GLOBAL VARIABLES <---# 

#==========================#
#    PRE-REQUISITES END    #
#==========================#


#================#
#   MAIN START   #
#================#
#---> STATUS MESSAGE <---#
print "Starting '$0'...\n\n" unless $quiet;

#---> MAIN VARIABLES <---#
my ($exon_hash_ref, $up_hash_ref, $down_hash_ref);

#---> BODY <---#
$exon_hash_ref = &bed_to_short_coord($exon);
$up_hash_ref = &bed_to_short_coord($up);
$down_hash_ref = &bed_to_short_coord($down);
&append_to_motevo($motevo, $exon_hash_ref, $up_hash_ref, $down_hash_ref);

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
'Usage: perl ./motevo_exon_up_or_down.pl [OPTIONS] --motevo [MOTEVO_OUT] --exon [BED] --up [BED] --down [BED]

Description: Adds a column to a MotEvo output file that contains either "UP", "DOWN", "EXON" or "N/A"; output written to STDOUT.

==================================================
Required arguments:
--motevo	MotEvo output file
--exon	exons BED file
--up	upstream of exons BED file
--down	downstream of exons BED file
==================================================
Optional arguments:
--usage|help	Show this information and die
--quiet	Shut up!
';
}
#-----------------------#
sub bed_to_short_coord {
### Function: Combines the coordinate values of a BED file into a short form: {chromosome}:{start}-{stop}{strand}, e.g. chr1:10000-20000+ and loads into a hash
### Accepts: Input filename
### Returns: Hash of short coordinates
### Dependencies: n/a
### Type: Specialized

	#---> PASS ARGUMENTS ---#
	my $file = shift;

	#---> STATUS MESSAGE <---#
	print STDERR "Processing file '$file'..." . "\n" unless $quiet;

	#---> SUBROUTINE VARIABLES <---#
	my %hash;

	#---> BODY <---#

		#---> Open file <---#
		open FILE, "<", $file;

		#---> Push line to array <---#
		while (<FILE>) {
			my ($chr, $start, $stop, $name, $score, $strand) = split /\t/;
			my $short_coord = $chr . ":" . $start . "-" . $stop . $strand;
			$hash{$short_coord} = "";
		}

		#---> Close file <---#
		close FILE;

	#---> STATUS MESSAGE <---#
	print STDERR "File '$file' processed." . "\n" unless $quiet;
	
	#---> RETURN VALUE <---#
	return \%hash;
	
}
#-----------------------#
sub append_to_motevo {
### Function: Adds a column to a MotEvo output file that contains either "UP", "DOWN", "EXON" or "N/A"; output written to STDOUT.
### Accepts: MotEvo output filename; 2. reference to hash containing short coordinates of exons ; 3. reference to hash containing short coordinates of upstream exon regions; 4. reference to hash containing short coordinates of downstream exon regions
### Returns: n/a
### Dependencies: n/a
### Type: Specialized
	#---> PASS ARGUMENTS ---#
	my ($motevo, $exon_hash_ref, $up_hash_ref, $down_hash_ref) = @_;
	
	#---> STATUS MESSAGE <---#
	print STDERR "Processing file '$motevo'..." . "\n" unless $quiet;
	
	#---> BODY <---#

		#---> Open file <---#
		open FILE, "<", $motevo;

		#---> Add column value based on the short coordinates <---#
		while (my $line = <FILE>) {

			#---> Chop newline character <---#
			chomp $line;

			#---> Extract short coordinates <---#
			my $short_coord = (split /\|/, (split /\t/, $line)[3])[0];

			#---> Print original line plus TAB <---#
			print $line . "\t";

			#---> Print "EXON" if short coordinates in exon hash <---#
			if ( defined $exon_hash_ref->{$short_coord} ) {
				print "EXON";
			}
			#---> Print "UP" if short coordinates in up hash <---#
			elsif ( defined $up_hash_ref->{$short_coord} ) {
				print "UP";
			}
			#---> Print "DOWN" if short coordinates in down hash <---#
			elsif ( defined $down_hash_ref->{$short_coord} ) {
				print "DOWN";
			}
			#---> Else print "N/A" <---#			
			else {
				print "N/A";
			}

			#---> Add newline character <---#
			print "\n";

		}

		#---> Close file <---#
		close FILE;
	
	#---> STATUS MESSAGE <---#
	print STDERR "File '$motevo' processed." . "\n" unless $quiet;

}
#=======================#
#    SUBROUTINES END    #
#=======================#

