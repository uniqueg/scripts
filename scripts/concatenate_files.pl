#!/usr/bin/perl

### Author: Alexander Kanitz
### Created: 14-JAN-2013/ Modified: 15-JAN-2013
### Description: Concatenates files with the indicated extension in a given folder into the indicated output file 
### Usage: perl /path/to/wrapper_pat_search.pl </path/to/files> <ext> <path/to/output_file>

### A. Pre-requisites
## Pragmas
use strict;
use warnings;
## Command-line arguments / initialization
my $dir = shift;
my $ext = shift;
my $out_file = shift;
###

### B. Main program
# Read files
my $files_array_ref = &getFilesWithExtFromFolder($dir, $ext);
# Concatenate files
&concatFiles($files_array_ref, $out_file);
# Exit
exit;
###

### C. Subroutines

sub getFilesWithExtFromFolder {
### Returns array of files with indicated extension in indicated folder
	## Pass arguments
	my $dir = shift;
	my $ext = shift;
	# Open directory
	opendir DIR, $dir;
	# Read fasta files from directory
	my @files = map {m/\.[^.]+$/;$_} grep {/${ext}$/} readdir DIR;
	# Close directory
	closedir DIR;
	# Screen output
	print "The following files were found:\n";
	## Add directory to each file
	foreach (@files) {
		$_ = $dir . $_;
		# Screen output
		print "$_\n";
	}
	# Return array reference
	return \@files;
}

sub concatFiles {
### Concatenates files in array and writes to output file
	# Pass arguments
	my @files = @{shift()};
	my $out_file = shift;
	# Open output file
	open OUT, ">$out_file";
	# Traverse through @files array
	foreach (@files) {
		# Open input bed file
		open BED, $_;
		# Read input bed file line by line
		while (<BED>) {
			print OUT;
		}
		# Close input bed file
		close BED;
	}
	# Close output file
	close OUT;
	# Screen output
	print "The files were concatenated into file \"$out_file\".\n";
}