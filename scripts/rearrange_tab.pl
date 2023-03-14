#!/usr/bin/perl

### Author: Alexander Kanitz
### Created: 08-FEB-2013
### Modified: 08-FEB-2013
### Adapted from: clipz_coordinates_to_bed.pl
### Description: Rearranges columns of a TAB file; constants may be introduced 
### Arguments: 1. input filename; 2. output filename; 3. default value for missing values; 4+. desired column index order in rearranged file (relative to input file); missing values will be replaced with the default value specified in 3.; if constant values shall be inserted at a specific location, use a string instead of a number instead; if the constant value to be inserted is a number (or for clarity), surround the replacement string/number by forward slashes (e.g. /3/) to distinguish from column indices; forward slashes are removed, so if forward slashes at the beginning and/or end are desired, use double slashes where applicable 
### Output: 1. Rearranged TAB file derived from input TAB file and - if applicable - manually assigned constants
### Usage: perl ./rearrrange_tab.pl input.tab output.tab N/A 5 4 /3/ constant 1

### A. Pre-requisites
## Pragmas
use strict;
use warnings;
## Command-line arguments / initialization
my $in_file = shift;
my $out_file = shift;
my $default_na = shift;
my $new_order_array_ref = \@ARGV;
###

### B. Main program
# 
rearrangeTab($in_file, $new_order_array_ref, $out_file);
# Print status message
print "Done.\n";
# Exit
exit;
###

### C. Subroutines
sub rearrangeTab {
### Accepts: 1. tab-delimited input file; 2. reference to array holding new order; 3. output file name; 4; default value for unavailable values (default: "N/A")
### Returns: 2. tab-delimited output file (rearranged)
### Type: Generic
	## Pass arguments
	my $in_file = shift;
	my $new_order_array_ref = shift;
	my $out_file = shift;
	my $default_NA = @_ ? shift : "N/A";
	# Open input file handle
	open IN, $in_file;
	# Open output file handle
	open OUT, ">$out_file";
	## Traverse through input file line by line
	while (<IN>) {
		# Remove trailing newline character
		chomp;
		# Split line by tab;
		my @in_line = split /\t/;
		# Declare variable for output line array
		my @out_line;
		## Traverse through each element in the array referenced by $new_order_array_ref
		for my $index (0..$#$new_order_array_ref) {
			## IF a non-numerical character or an entry surrounded by '/' for the value of the current index of the new order array, use as string (i.e. write found value, minus surrounding slashes if available, to @out_line element)
			if ($$new_order_array_ref[$index] =~ m:^0-9]: | $$new_order_array_ref[$index] =~ m:^/.*/$:) {
				# Assign value of the current index of the new order array (i.e. a string) to a helper variable
				my $string = $$new_order_array_ref[$index];
				# Remove surrounding slashes '/' if present
				$string =~ s:^/(.*)/$:$1:;
				# Assign string to value of the current index of the @out_line array
				$out_line[$index] = $string;
			}
			## ELSE... 
			else {
				## ...check IF value in column specified by new order array reference is present (mind Perl array index numbering!)
				if (defined $in_line[$$new_order_array_ref[$index] - 1]) {
					## If so, assign to @out_line array element with the indicated column index (mind Perl array index numbering!)
					$out_line[$index] = $in_line[$$new_order_array_ref[$index] - 1];
				}
				## ELSE print warning message and assign default value for missing values
				else {
					# Get column number (mind Perl array index numbering!)
					my $column = $index + 1;
					# Print warning message to STDOUT
					print "Warning: No value in row $., column $column, writing '$default_NA' instead.\n";
					# Assign default value for missing values to @out_line array element with the indicated column index
					$out_line[$index] = $default_NA;
				}
			}
		}
		# Join output line array elements by tab
		print OUT join("\t", @out_line) . "\n";
	}
	# Close output file handle
	close OUT;
	# Close input file handle
	close IN;
	## Status messages
	print "Processed file '$in_file'.\n";
	print "Output written to file '$out_file'.\n";
}
###