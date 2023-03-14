#!/usr/bin/perl

### Author: Alexander Kanitz
### Created: 18-JAN-2013
### Modified: 21-JAN-2013


### Pre-requisites

## Pragmas
use strict; 
use warnings;
## Used Packages
use LWP::UserAgent;
use File::Path;
## Initialize global variables
my $url_base = "http://lyrics.chambermusik.com/";
my $url = "main.html";
my $out_dir = "/home/kanitz/Dropbox/Private/Coding/Rap_Lyrics/Chambermusik/";
my $log = "/home/kanitz/Dropbox/Private/Coding/Rap_Lyrics/Chambermusik/log";
## Declare global variables
my $main;
my $hrefs_array_ref;


### Main

# Initialize log file
&initializeLog($log);
# Load root html file to variable
$main = &returnFileFromWWW("$url_base$url", $log);
# Extract hyperlinks from root html file
$hrefs_array_ref = &extractHyperLinksFromString($main);
# Add URL basename to hyperlink array 
&addPrefixToArrayElementsInPlace($hrefs_array_ref, $url_base); 
# Download files, keep original directory structure
&saveFilesFromWWWKeepDirStructure($hrefs_array_ref, $url_base, $out_dir, $log);
# Print status message
print "Done.";
# Exit
exit;


### Sub-routines

sub initializeLog {
	# Pass argument
	my $log = shift;
	# Open log file handle (new entries will be appended!)
	open LOG, ">>$log";
	# Get date/time stamp
	my $time = &dateTimeStamp;
	## Check whether log file is empty 
	if (-z $log ) {
    	# If empty, print date/time stamp only
    	print LOG $time . "\n";
	}
	else {
		# Else, add separator
		print LOG "\n------------\n\n" . $time . "\n";	
	}
    # Close log file handle
    close LOG;
}
sub dateTimeStamp {
### Code modified from http://www.netadmintools.com/art212.html (18-JAN-2013)
	# Obtain date and time through localtime function
	my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time);
	# Print date/time stamp to variable $stamp
	my $stamp = sprintf "%4d-%02d-%02d %02d:%02d:%02d", $year+1900, $mon+1, $mday, $hour, $min, $sec;
	# Return $stamp
	return $stamp;
}
sub returnFileFromWWW {
### Code modified from http://www.caveofprogramming.com/perl/perl-wget-retrieving-html-pages-in-perl/ (18-JAN-2013)
	## Pass arguments
	my $file = shift;
	my $log = shift;
	# Create fake browser (user agent)
	my $ua = LWP::UserAgent->new();
	# Get file
	my $response = $ua->get($file);
	# Open log file handle
	open LOG, ">>$log";
	## Catch file error
	unless($response->is_success) {
		# Print error message to log file
		print LOG "Error: " . $response->status_line;
	}
    # Save result to variable $return
    my $return = $response->decoded_content;
	# Print status message to log file
    print LOG "Processed file '$file' (" . length($response->decoded_content) . " bytes).\n";
    # Close log file handle
    close LOG;
    # Return $return
    return $return;
}
sub extractHyperLinksFromString {
	# Pass argument
	my $html = shift;
	# Find hyperlinks using regular expression and put into array
	my @hrefs = $html =~ /href="(.*?html?)"/g;
	# Remove duplicates
	&rmArrayDuplicatesInPlace(\@hrefs);
	# Return array reference
	return \@hrefs;
}
sub rmArrayDuplicatesInPlace {
	# Pass argument
	my $array_ref = shift;
	# Map array to hash
	my %hash = map { $_, 1 } @{$array_ref};
    # Re-load array with hash keys
    @{$array_ref} = keys %hash;
}
sub addPrefixToArrayElementsInPlace {
	## Pass arguments
	my $array_ref = shift;
	my $prefix = shift;
	## Crawl through each array element
	foreach (@{$array_ref}) {
		# Add prefix
		$_ = $prefix . $_;
	}
}
sub saveFilesFromWWWKeepDirStructure {
### Code modified from http://www.caveofprogramming.com/perl/perl-wget-retrieving-html-pages-in-perl/ (18-JAN-2013)
	## Pass arguments
	my $files_array_ref = shift;
	my $base = shift;
	my $out_dir = shift;
	my $log = shift;
	## Crawl through all URLs in array
	foreach (@{$files_array_ref}) {
		# Get file
		my $href = &returnFileFromWWW($_, $log);
		# Extract hyperlinks
		my $txt_array_ref = &extractTextFilesFromString($href);
		## Crawl through all hyperlinks in URL
		foreach my $file (@{$txt_array_ref}) {
			# Split file			
			my @path = split /\//, $file;
			# Remove last element
			pop @path;
			# Re-join @path array
			my $path = join '/', @path;
			# Create folders
			mkpath("$out_dir$path", 0, 0755);
			# Create fake browser (user agent)
			my $ua = LWP::UserAgent->new();
			# Get file
			my $response = $ua->get("$base$file");
			# Open log file handle
			open LOG, ">>$log";
			## Catch file error
			unless($response->is_success) {
				# Print error message to log file
				print LOG "Error: " . $response->status_line;
			}
			# Open output file handle
			open OUT, ">$out_dir$file";
			# Avoid 'wide characters in print' warning
			binmode OUT, ":utf8";
		    # Print result to output file
		    print OUT $response->decoded_content;
			# Close output file handle 
		    close OUT;
			# Print status message to log file
		    print LOG "Downloaded file '$base$file' (" . length($response->decoded_content) . " bytes)\nLocation: '$out_dir$file'\n";			
			# Close log file handle
   		 	close LOG;
		}
	}
}
sub extractTextFilesFromString {
	# Pass argument
	my $html = shift;
	# Find hyperlinks using regular expression and put into array
	my @txt = $html =~ /href="(.*?txt)"/g;
	# Remove duplicates
	&rmArrayDuplicatesInPlace(\@txt);
	# Return array reference
	return \@txt;
}