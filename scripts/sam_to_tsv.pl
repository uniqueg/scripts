#!/usr/bin/perl

#######
### AUTHOR INFO & CREATION DATE:
### ----------------------------
### Author:		Alexander Kanitz
### Created: 	17-DEC-2012
### Modified:	18-DEC-2012
#######

#######
### FUNCTION:
### ---------
### The script converts one or more SAM alignment format (.sam) files into tab-separated value files that can be analyzed with 'compare_tsv.pl'.
#######

#######
### ARGUMENTS:
### ----------
### 1. A folder containing the bam alignment files. Note that NO '/' is accepted at the end of the folder.
### 2. A pattern (mind Perl usage of patterns!) indicating which files to process. Note that the files that are to be processed MUST have the file extension '.tsv' which is automatically added and therefore should NOT be included in the pattern. If all '.tsv' files in the folder are to be processed (preferred, as the pattern matching was not tested extensively!), indicate '.*'.
### General Notes: Compare the examples in USAGE for further information.
#######

#######
### OUTPUT:
### -------
### For each line in the SAM alignment file, one line per alignment with the following entries is written to <FILE>.tsv (in the same directory as the input file): 1. sequence ID, 2. chromosome, 3. strand, 4. start position, 5. end position, 6. sequence.
#######

#######
### USAGE:
### ------
### perl /path/to/compare_tsv.pl </path/to/files> <pattern>
### Examples:
### 1. perl ~/compare_tsv.pl ~/sam_files .* (converts all files in folder '~/sam_files' with extension '.sam' to .tsv files with the same base names; script in home directory)
### 2. perl ~/compare_tsv.pl . ^MAPPER (converts all files starting with MAPPER and the extension .sam in the current folder to .tsv files with the same base names; script in home directory) 
#######

#######
### OVERVIEW:
### ---------
### A. Pre-requisites
### B. File handling
### C. Conversion
### D. Clean-up
### E. Subroutines
#######

### A. Pre-requisites
## Define pragmas, read arguments etc.
## Define used pragmas and import libraries
use warnings;
use strict;
# Define error message for wrong usage
my $usage = 'Usage: perl /path/to/compare_tsv.pl </path/to/files> <pattern>\n';
## Assign arguments (path, reference file, pattern) to variables
my $path = shift or die "$usage";
my $pat = shift or die "$usage";
# Append file extension '.tsv' to pattern
$pat = $pat."\.sam\$";
# Change working directory to indicated path 
chdir($path);
###

### B. File handling 
### Open indicated directory, traverse files and add to array of mapping data files if pattern matches
# Declare array variable for mapping data files
my @sam_files;
# Open indicated directory or die
opendir (DIR, $path) or die $!;
# Push each file in directory that matches the indicated pattern + '.tsv' in array
while (my $file = readdir(DIR)) {
	push (@sam_files, $file) if ($file =~ m/$pat/);
}
# Close indicated directory
closedir (DIR);
###

### C. Conversion
### Traverse through array of files, convert each line of each file and write results to output files  
foreach (@sam_files) { 
	my $sam_file = $_;
	
	## D1. Comparisons/Counts
	# Open file handle that processes input .sam file
	open (SAM, $sam_file);
	# Open file handle that is used to write output .tsv file
	open (TSV, ">$sam_file.tsv");
	# Traverse through mapping data file line by line
	while (<SAM>) {
		# Check if header line; if not, process
		if (substr($_,0,1) ne "@") {
			chomp($_);
			# Split line by tabs
			my @cols = split("\t");
			# Extract strand information from sam flag (position 5) using dec2bin subroutine
			my $str = substr(&dec2bin($cols[1]), 27, 1);
			if ($str == 0) {
				$str = "+";
			}
			else {
				$str = "-";
			}
			# Calculate length
			my $len = length($cols[9]);
			# Calculate end position
			my $end = $cols[3] + $len - 1;
			## If coordinates are available and complete, print full information  
			if ( ($cols[2] ne "") && ($cols[2] ne "*") && ($cols[3] ne "") ) {
				print TSV "$cols[0]\t$cols[2]\t$str\t$cols[3]\t$end\t$len\n";
				## Test if multiple alignments are present in the tags (bwa -> 'XA:Z:' tag)
				# Find 'XA:Z:' tag instances
				my @tags = grep { /XA:Z:/ } @cols;
				# Check whether tag was found
				if ($tags[0]) {
					# Trim leading tag name ('XA:Z:') and trailing semicolon
					$tags[0] = substr($tags[0],5,length($tags[0])-6);
					# Split multiple entries by separator (';')
					my @mm_xy = split(";", $tags[0]);
					# Traverse through each set of coordinates
					foreach (@mm_xy) {
						# Split individual entries by separator (',')
						my @mm = split(",");
						# Declare variables
						my ($mm_str, $mm_start, $mm_end);
						## Check for presence of '-' before start positon; if found, set strand to "-", substring starting position and derive end position 
						if (substr($mm[1],0,1) eq "-") {
							$mm_str = "-";
							$mm_start = substr($mm[1],1);
							$mm_end = $mm_start + $len - 1;	
							# Print alternative alignment if edit distance is '0'
							print TSV "$cols[0]\t$mm[0]\t$mm_str\t$mm_start\t$mm_end\t$len\n" if $mm[3] == 0;
						}
						## Check for presence of '-' before start positon; if found, set strand to "-", substring starting position and derive end position
						elsif (substr($mm[1],0,1) eq "+") {
							$mm_str = "+";
							$mm_start = substr($mm[1],1);
							$mm_end = $mm[1] + $len - 1;	
							# Print alternative alignment if edit distance is '0'
							print TSV "$cols[0]\t$mm[0]\t$mm_str\t$mm_start\t$mm_end\t$len\n" if $mm[3] == 0;
						}
					}
				}
			}
			## If no coordinates are available, print "NA"
			else {
				print TSV "$cols[0]\tNA\tNA\tNA\tNA\t$len\n"; 
			}
		}
	}
	# Close input and output file handles
	close (TSV);
	close (SAM);
	# Print status message to STDOUT
	print "File $sam_file converted to $sam_file.tsv.\n";
}
	
### D. Clean-up
exit;
###

### E. Subroutines

## E1. Convert integer to 32-digit binary number
sub dec2bin {
    my $str = unpack("B32", pack("N", shift));
    return $str;
}
###