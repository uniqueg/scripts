#!/usr/bin/perl

#==================#
#   HEADER START   #
#==================#
### Name: merge_duplicates.pl
### Created: Mar 11, 2013
### Modified: Mar 11, 2013
### Author: Alexander Kanitz
### Company: Zavolan Group, Biozentrum, University of Basel
### Version: v0.1
### Adapted from: merge_duplicates_from_clipz_overlap_cl_sites.pl
### Requirements: n/a
#==================#
### Description: Merges all TAB files with filenames matching a specified pattern in a specified path based on matching values in the specified columns; rules specifying how to treat different values in columns that are not indicated, have to be specified in the subroutine 'merge_rules' 
### Output: One TAB file with merged entries (mind merging rules!) is produced for each input file
### Usage: perl ./merge_duplicates.pl for information on required and optional arguments
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
my $in_path = '';
my $pattern = '';
my $columns = '';
my $out_path = '';
my $prefix = 'merged_';
my $options_result = GetOptions (
	'usage|help' => \$usage,
	'quiet|no_verbose' => \$quiet,
	'column=i@' => \$columns,
	#-----------------------#
	'in-path=s' => \$in_path,
	'pattern=s' => \$pattern,
	'out-path=s' => \$out_path,
	'prefix=s' => \$prefix
);
die $usage_info if $usage || !$options_result;
die $usage_info if !$in_path || !$pattern || !$out_path; 
#==========================#
#    PRE-REQUISITES END    #
#==========================#


#================#
#   MAIN START   #
#================#
#---> STATUS MESSAGE <---#
print "Starting 'merge_duplicates.pl'...\n\n" unless $quiet;

#---> MAIN VARIABLES <---#
my $filenames_array_ref;

#---> BODY <---#
$filenames_array_ref = &files_regex($in_path, $pattern);
&merge_main($filenames_array_ref, $columns, $out_path);

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
### Type: Specialized
'Usage: perl ./merge_duplicates.pl [OPTIONS] --in-path [PATH] --pattern [STRING | REGEX] --column [INTEGER] --out-path [PATH]
==================================================
Required arguments:
--in-path	Path to input files
--pattern	Pattern for input file selection
--column	Index of column that is merged if identical values occur (1 means 1st column!); multiple columns can be specified, but require an additional "--column" each 
--out-path	Path to output files
==================================================
Optional arguments:
--usage|help	Show this information and die
--quiet		Shut up!
--prefix	String added at the beginning of output filename to indicate merging [STRING | DEFAULT = "merged_"]
';
}
#-----------------------#
sub files_regex {
### Function: Finds all files in a specified folder that contain a specified pattern in their filenames
### Accepts: 1. Path to folder; 2. Pattern [STRING]
### Returns: Reference to array of filenames
### Dependencies: n/a
### Type: Generic
	#---> PASS ARGUMENTS ---#
	my ($path, $pattern) = @_;
	
	#---> SUBROUTINE VARIABLES <---#
	my @filenames;
	
	#---> BODY <---#
	# Open directory handle DIR or die
	opendir (DIR, $path) or die $!;
	#-----------------------#
	## Traverse through files in DIR
	while (my $file = readdir(DIR)) {
 		# Push file to array @filenames if filename contains $pattern 
		push @filenames, $file if $file =~ m/$pattern/;
	}	
	#-----------------------#
	# Close directory handle DIR
	closedir (DIR);
	
	#---> STATUS MESSAGE <---#
	print "Directory '$path' processed. " . scalar @filenames . " file(s) matched the indicated pattern '$pattern'." . "\n" unless $quiet;
	
	#---> RETURN VALUE <---#
	return \@filenames;
}
#-----------------------#
sub merge_main {
### Function: Main routine coordinating file by file processing of input files 
### Accepts: 1. Reference to array of filenames; 2. Reference to array of columns to base the merging on; 3. Path to output files
### Returns: n/a
### Dependencies: n/a
### Type: Specialized
	#---> PASS ARGUMENTS ---#
	my ($filenames_array_ref, $columns_array_ref, $out_path) = @_;

	#---> STATUS MESSAGE <---#
	print "\nStarting to search for duplicate entries..." . "\n" unless $quiet;
	
	#---> BODY <---#
	## Traverse through each element $file of array @$filenames_array_ref
	foreach my $file ( @$filenames_array_ref ) {
		# Read file
		my $file_aoa_ref = tab_to_aoa($file);
	    # Merge file
		my $merged_file_aoa_ref = &aoa_merge($file_aoa_ref, $columns_array_ref);
		# Write to output file
		&aoa_to_file($merged_file_aoa_ref, "$out_path/$prefix$file");
		
		#---> STATUS MESSAGE <---#
		print "File '$file' processed." . "\n\n" unless $quiet;
	}
}
#-----------------------#
sub tab_to_aoa {
### Function: Reads a TAB file into an array of arrays
### Accepts: TAB file
### Returns: Reference to array of arrays (elements outer arrays = TAB file rows; elements inner arrays = values of single TAB file row)
### Dependencies: n/a
### Type: Generic
	#---> PASS ARGUMENTS ---#
	my $in_file = shift;	
	
	#---> SUBROUTINE VARIABLES <---#
	my @AoA;	
	
	#---> BODY <---#
	# Open file handle TAB
	open TAB, $in_file;
	#-----------------------#
	## Crawl through TAB file line by line
	while (<TAB>) {
		# Remove trailing new line characters
		chomp;
		# Split line into array @line by tabs
		my @line = split /\t/, $_;
		# Push line values array @line to respective hash value / outer array
		push @AoA, [ @line ];
	}		
	#-----------------------#
	# Close TAB file handle
	close TAB;
	#-----------------------#
		
	#---> STATUS MESSAGE <---#
	print "Read TAB file '$in_file'" . "\n" unless $quiet;
	
	#---> RETURN VALUE <---#
	return \@AoA;
}
#-----------------------#
sub aoa_merge {
### Function: Merges outer array elements of a specified array of arrays when a specified set of inner array elements is identical; rules how to deal with non-identical inner array elements (that were not specified) are defined in the non-generic subroutine 'merge_rules'
### Accepts: 1. Reference to array of arrays; 2. Reference to array of indices (starting with 1!), indicating the inner array elements which require identical values for merging of entries to occur
### Returns: Reference to array of arrays
### Dependencies: Subroutine 'merge_rules'
### Type: Generic (subrouting 'merge_rules' needs to be adapted!)
	#---> PASS ARGUMENTS ---#
	my ($aoa_ref, $columns_array_ref) = @_;
	
	#---> SUBROUTINE VARIABLES <---#
	# Declare temporary hash %hoa
	my %hoa;
	
	#---> BODY <---#
	## Traverse through each element $array_ref of array @$aoa_ref
	foreach my $array_ref ( @$aoa_ref ) {
	    # Declar variable for temporary array used for building the %hoa hash key
	    my @key;
		## Traverse through each element $column of array @$columns_array_ref
		foreach my $column ( @$columns_array_ref ) {
		    # Add value of column ($column - 1) in array @$array_ref to hash key array @key
		    push @key, @$array_ref[$column - 1];
		} 	      
		# Concatenate elements of hash key array @key
		my $key = join "_", @key;
		## Check whether hash %hoa entry with key $key already exists
		if (exists $hoa{$key}) {
			# If so, call subroutine defining merging rules
			$hoa{$key} = &merge_rules($hoa{$key}, $array_ref);
		} 		
		## Else...
		else {
			# Assign array reference $array_ref to hash %hoa key $key
			$hoa{$key} = $array_ref;
		}		
	}	
	#-----------------------#	
	# Assign values of hash %hoa back to $aoa_ref
	@$aoa_ref = sort values %hoa;
	
	#---> STATUS MESSAGE <---#
	print "Merged TAB file entries according to identical values in the following columns: @$columns_array_ref" . "\n" unless $quiet;	
	
	#---> RETURN VALUE <---#
	return $aoa_ref;
}
#-----------------------#
sub merge_rules {
### Function: Defines merging rules for columns whose values do not match (scores etc.) 
### Accepts: 1. Reference to array containing reference/target value; 2. Reference to array containing query value
### Returns: Reference to array with (updated) value determined by merging rules
### Dependencies: n/a
### Type: Specific
	#---> PASS ARGUMENTS ---#
	my ($target_array_ref, $query_array_ref) = @_;

	#---> BODY <---#
	# Sum values										###	REPLACE WITH
	$$target_array_ref[1] += $$query_array_ref[1];		### DESIRED OPERATION

	# Return array reference
	return $target_array_ref;	
}
#-----------------------#
sub aoa_to_file {
### Description: Writes the indicated, referenced array of arrays to the specified output file; the default settings produce a tab-separated file, but the intra- and interarray separators may be specified manually
### Accepts: 1. Reference to array of arrays; 2. Output file name; 3. String used to separate elements of inner arrays [DEFAULT = TAB "\t"]; 4. String used to separate outer array elements [DEFAULT = NEWLINE "\n"]; NOTE: If any one of the separators is set manually, always specify both!
### Returns: n/a 
### Dependencies: n/a
### Type: Generic
	#--- PASS ARGUMENTS ---#
	my ($aoa_ref, $out_file, $sep_intra, $sep_inter) = (@_, "\t", "\n");

	#--- BODY ---#
	# Open file handle OUT
	open OUT, ">$out_file";
	#-----------------------#
	## Loop through references to inner arrays (elements of outer array)
	foreach my $array_ref (@$aoa_ref) {
		# Print current inner array to output file; array elements are separated by $sep_intra (e.g. tab), different arrays by $sep_inter (e.g. newline)
		print OUT join($sep_intra,@$array_ref) . $sep_inter;
	}
	#-----------------------#
	# Close file handle OUT
	close OUT;

	#---> STATUS MESSAGE <---#
	print "Output written to file '$out_file'.\n" unless $quiet;	
}
#=======================#
#    SUBROUTINES END    #
#=======================#