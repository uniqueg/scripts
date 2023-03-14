#!/usr/bin/perl

### Author: Alexander Kanitz
### Created: 19-JAN-2013
### Modified: 21-JAN-2013


### Pre-requisites

## Pragmas
use warnings;
use strict;
# Initialize directory variable $dir
my $dir = '/media/MUSIC, PHOTOS, BOOKS & PERSONAL/DOCUMENTS/Dropbox/Private/Coding/Rap_Lyrics/Chambermusik/';
## Initialize output file locations
my $out = '/media/MUSIC, PHOTOS, BOOKS & PERSONAL/DOCUMENTS/Dropbox/Private/Coding/Rap_Lyrics/Chambermusik/output';
## Declare variables
my %HoH;
my $HoH_ref = \%HoH;

### Main

# Get file list array reference
my $files_array_ref = &recurGetTxtFilesFromRootDir($dir);
# Load hash of hashes of hashes of arrays reference (artist, album, song, verses) from file list array reference
my $HoHoHoA_ref = &song2ArtistAlbumSongVerseHoHoHoA($files_array_ref);
# Calculate general stats: total verses, bars, words, characters
$HoH_ref = &generalStats($HoHoHoA_ref, $HoH_ref);
# Calculate word usage: guns
$HoH_ref = &wordUsageGuns($HoHoHoA_ref, $HoH_ref);
# Calculate different words per word ratio for each artist; total
#$HoH_ref = &wordUsageTotal($HoHoHoA_ref, $HoH_ref);
# Calculate different words per word ratio for each artist; per verse
#$HoH_ref = &wordUsageVerse($HoHoHoA_ref, $HoH_ref);
# Print results to file
&printHoHToTab($HoH_ref, $out);
# Print status
print "Done.\n";
# Exit
exit;

# Verses
# Bars
# Words
# Unique words
# Characters


### Subroutines

sub recurGetTxtFilesFromRootDir {
### Modified from http://www.perlmonks.org/?node_id=136482 (19-JAN-2013)
### Accepts 1: full path to directory
### Returns: list of '.txt' files in directory and all subdirectories
	# Pass argument
	my $path = shift;
	# Open directory handle
	opendir DIR, $path
		# Catch error
		or die "Unable to open $path: $!";
	# Obtain files from directory
	my @files =
		# Third: Prepend the full path
		map { $path . '/' . $_ }
		# Second: take out '.' and '..'
		grep { !/^\.{1,2}$/ }
		# First: get all files
		readdir (DIR);
	# Close directory handle
	closedir DIR;
	## Crawl through all files 
	for (@files) {
		# Add all new files from directory and subdirectories
		push @files, @{&recurGetTxtFilesFromRootDir($_)} if -d $_;
	}
	# Keep only files with extension .txt that are not directories	
	@files= grep { (/\.txt$/) && (! -l $_) } @files;
	# Return files array reference  
	return \@files;
}
sub song2ArtistAlbumSongVerseHoHoHoA {
	# Print status
	print "Creating lookup table...\n";
	## Pass arguments
	my @files = @{shift()};
	# Declare variables
	my ($album, $song, %HoHoHoA); 												# default to basename $txt!!!!!!!!!!!!!!
	## Crawl through all files in array @files
	foreach my $txt (@files) {
		# Open text file handle
		open TXT, $txt;
		## Crawl through txt file line by line
		while (<TXT>) {
			# Remove trailing newline character
			chomp;
			# Assign line to dedicated variable
			my $line = $_;
			# Extract album information
			$album = $line if $line =~ s/^album:\s+//i;
			# Extract song title and concatenate to album
			$song = $album . ":" . $line if $line =~ s/^song:\s+//i;
			# Start verse recording IF artist info line, i.e. line
			if	(
				$line =~ /^\[.*\]$/											&&		# starts with '[' and ends with ']' 				AND
				$line =~ /[A-Z]/											&&		# contains at least one uppercase character			AND
				$line !~ /chorus|hook|interlude|intro|outro|sample|skit/i	&&		# does not contain any of the mentioned keywords	AND
				$line !~ /break|bridge|chrous|chours|intrlude/i				&&		# does not contain any of the mentioned keywords	AND
				$line !~ /\(|\)|\{|\}|\?|,|\].*\[|\*.*\*|<.*>]/						# does not contain (, ), {, }, ?, comma, ]..[, *..* or <..>
				#!# incorporate  
				) {
					# Assign artist line to dedicated variable
					my $artist = $line;
					# Declar verse variable
					my $verse = "";
					# Extract artist (disregards further artists given in round brackets)
					$artist =~ s/^\[(.*)(\s\(.*\))?\]$/$1/;
					## Crawl through verse line by line
					while (<TXT>) {
						# Assign verse line to dedicated variable
						my $verse_line = $_;
						# Exit loop if current line starts with a square bracket '['
						last if $verse_line =~ /^\[/;
						## Skip line IF
						if  (
							$verse_line eq "\n" 					||		# empty												OR
							$verse_line =~ /^\{.*\}$/				||		# starts and ends with curly brackets				OR
							$verse_line =~ /->/								# contains an arrow '->'
							) {
							next;	
						}
						# Add current line to $verse
						$verse .= $verse_line; 
					}
	    			# Push verse to HoHoHoA (artist, album, song, verses) if $verse is >= than 250 characters
	    			push @{ $HoHoHoA{$artist}{$album}{$song} }, $verse if length $verse >= 250;
					}
		}
		# Close bed file handle
		close TXT;
	}
	# Return bed file hash of hashes of arrays
	return \%HoHoHoA;
}
sub generalStats {
	# Print status
	print "Calculating general statistics...\n";
	# Pass arguments
	my %HoHoHoA = %{shift()};
	my %HoH = %{shift()};
	## Crawl through each artist
	foreach my $artist (sort keys %HoHoHoA) {
		# Declare variables
		my (@verses, @words, @bars, @chars);
		## Crawl through each album		
		foreach my $album (keys %{$HoHoHoA{$artist}}) {
			## Crawl through each song
			foreach my $song (keys %{$HoHoHoA{$artist}{$album}}) {
				## Crawl through each verse
				foreach my $verse (@{$HoHoHoA{$artist}{$album}{$song}}) {
					# Add individual verses to dedicated array
					push @verses, $verse;
					# Add individual bars to dedicated array
					push @bars, split(/\n/, $verse);
					# Remove punctuation characters
					$verse = &rmPunctuation($verse);
					# Add individual words to dedicated array
					push @words, split(/\s+/, $verse);
					# Remove whitespace
					$verse =~ s/\s+//g;
					# Add individual characters to dedicated array
					push @chars, split(//, $verse);
				}
			}
		}
		# Put data to hash
		$HoH{$artist}{'total_verses'} = scalar(@verses);
		$HoH{$artist}{'total_bars'} = scalar(@bars);
		$HoH{$artist}{'total_words'} = scalar(@words);
		$HoH{$artist}{'total_characters'} = scalar(@chars);
		$HoH{$artist}{'chars/word'} = scalar(@chars)/scalar(@words);
		$HoH{$artist}{'chars/bar'} = scalar(@chars)/scalar(@bars);
		$HoH{$artist}{'chars/verse'} = scalar(@chars)/scalar(@verses);
		$HoH{$artist}{'words/bar'} = scalar(@words)/scalar(@bars);
		$HoH{$artist}{'words/verse'} = scalar(@words)/scalar(@verses);
		$HoH{$artist}{'bars/verse'} = scalar(@bars)/scalar(@verses);
	}
	# Return
	return \%HoH;
}
sub wordUsageTotal {
	# Pass arguments
	my %HoHoHoA = %{shift()};
	my %HoH = %{shift()};
	## Crawl through each artist
	foreach my $artist (sort keys %HoHoHoA) {
		## Declare variables
		my (@total_words, @unique_words);
		my ($ratio, $total_words, $unique_words);
		## Crawl through each album		
		foreach my $album (keys %{$HoHoHoA{$artist}}) {
			## Crawl through each song
			foreach my $song (keys %{$HoHoHoA{$artist}{$album}}) {
				## Crawl through each verse
				foreach my $verse (@{$HoHoHoA{$artist}{$album}{$song}}) {
					# Remove punctuation characters
					$verse = &rmPunctuation($verse);
					# Add individual words to total words array
					push @total_words, split(/\s+/, $verse);
				}
			}
		}
		# Compute unique words
		@unique_words = @{&rmArrayDuplicates(\@total_words)};
		# Calculate number of total words
		$total_words = @total_words;
		# Calculate number of unique words
		$unique_words = @unique_words;
		# Calculate ratio of unique to total words		
		$ratio = $unique_words / $total_words;
		# Put data to hash
		$HoH{$artist}{'total_words'} = $total_words;
		$HoH{$artist}{'unique_words'} = $unique_words;
		$HoH{$artist}{'unique/total_words'} = $total_words/$unique_words;
	}
	# Return
	return \%HoH;
}
sub wordUsageVerse {
	# Pass arguments
	my %HoHoHoA = %{shift()};
	my %HoH = %{shift()};
	## Crawl through each artist
	foreach my $artist (sort keys %HoHoHoA) {
		## Declare variables
		my @ratios;
		my $mean_ratio;
		## Crawl through each album		
		foreach my $album (keys %{$HoHoHoA{$artist}}) {
			## Crawl through each song
			foreach my $song (keys %{$HoHoHoA{$artist}{$album}}) {
				## Crawl through each verse
				foreach my $verse (@{$HoHoHoA{$artist}{$album}{$song}}) {
					## Declare variables
					my (@total_words, @unique_words);
					my $ratio;
					# Remove punctuation characters
					$verse = &rmPunctuation($verse);
					# Add individual words to total words array
					push @total_words, split(/\s+/, $verse);
					# Compute unique words
					@unique_words = @{&rmArrayDuplicates(\@total_words)};
					# Calculate ratio of unique to total words		
					$ratio = @unique_words / @total_words;
					# Add ratio to array
					push @ratios, $ratio;					
				}
			}
		}
		# Calculate mean
		$mean_ratio = &meanArray(\@ratios);
		# Print output
		$HoH{$artist}{'unique/total_words_per_verse'} = $mean_ratio;
	}
	# Return
	return \%HoH;
}
sub meanWordLength {
	# Pass arguments
	my %HoHoHoA = %{shift()};
	my %HoH = %{shift()};
	## Crawl through each artist
	foreach my $artist (sort keys %HoHoHoA) {
		## Declare variables
		my @total_words;
		my $mean_length;
		## Crawl through each album		
		foreach my $album (keys %{$HoHoHoA{$artist}}) {
			## Crawl through each song
			foreach my $song (keys %{$HoHoHoA{$artist}{$album}}) {
				## Crawl through each verse
				foreach my $verse (@{$HoHoHoA{$artist}{$album}{$song}}) {
					# Remove punctuation characters
					$verse = &rmPunctuation($verse);
					# Add individual words to total words array
					push @total_words, split(/\s+/, $verse);
				}
			}
		}
		## Replace words in array with word lengths
		foreach (@total_words) {
			$_ = length($_);
		}
		# Calculate mean
		$mean_length = &meanArray(\@total_words);
		# Print output
		$HoH{$artist}{'mean_word_length'} = $mean_length;
	}
	# Return
	return \%HoH;
}
sub wordUsageGuns {
	# Print status
	print "Calculating trigger happy index...\n";
	# Pass arguments
	my %HoHoHoA = %{shift()};
	my %HoH = %{shift()};
	# Initialize word list regex
	my $regex = '(gun|firearm|tool|ratchet|blast|pop|popp|culture power|44|45|culture|slug|bullet|tech|tec|ruger|wesson|wessun|heat|AK|spray|355|click|walther|burner|chrome|clip|weapon|gat|heater|nine|strapped|9|four five|four four|bang|millimeter|uzi|47|22|shoot|shot|fire)(|\'s|s|in|ed|ing)';
	## Crawl through each artist
	foreach my $artist (sort keys %HoHoHoA) {
		# Declare variables
		my @hits_total;
		## Crawl through each album		
		foreach my $album (keys %{$HoHoHoA{$artist}}) {
			## Crawl through each song
			foreach my $song (keys %{$HoHoHoA{$artist}{$album}}) {
				## Crawl through each verse
				foreach my $verse (@{$HoHoHoA{$artist}{$album}{$song}}) {
					# Remove punctuation characters
					$verse = &rmPunctuation($verse);
					# Add individual words to dedicated array
					my @words = split(/\s+/, $verse);
					# Search for all occurrences of either of the words in $regex in @words
					my @hits = grep(/$regex/, @words);
					# Add to @hits_total
					@hits_total = (@hits_total, @hits);
				}
			}
		}
		# Put data to hash
		$HoH{$artist}{'trigger_happy_index'} = scalar(@hits_total);
	}
	# Return
	return \%HoH;
}	
sub printHoHToTab {
	## Pass arguments
	my %HoH = %{shift()};
	my $out = shift;
	# Open output file handle
	open OUT, ">$out";
	## PRINT HEADER
	# Print first column header
	print OUT "artist\t";
	# Get keys of outer hash
	my @keys = keys %HoH;
	# Get random key
	my $random_key = $keys[rand @keys];
	## Crawl through inner hash
	foreach my $inner (keys %{$HoH{$random_key}}) {
			# Print key
			print OUT "$inner\t";
		}
	# Print newline
	print OUT "\n";	
	## PRINT DATA
	## Crawl through outer hash
	foreach my $outer (keys %HoH) {
		# Print outer hash key
		print OUT "$outer\t";
		## Crawl through inner hash
		foreach my $inner (keys %{$HoH{$outer}}) {
			print OUT "$HoH{$outer}{$inner}\t";
		}
		print OUT "\n";
	# Close output file handle
	} 
	close OUT;
}
sub rmArrayDuplicates {
	# Pass argument
	my @array = @{shift()};
	# Map array to hash
	my %hash = map { $_, 1 } @array;
    # Re-load array with hash keys
    @array = keys %hash;
    # Return array reference
    return \@array;
}
sub rmPunctuation {
	# Pass argument
	my $string = shift;
	# Remove punctuation
	$string =~ s/[\!\?\"\'\;\:\,\.\{\}\(\)\<\>\[\]\-]//g;
	# Return string
	return $string;
}
sub meanArray {
	# Pass argument
	my @array = @{shift()};
	# Declar variable
	my $sum = 0;
	## Sum lengths of words in word array
	foreach (@array) {
		$sum += $_;
	}
	# Divide by number of words
	my $mean = $sum / @array;
	# Return mean
	return $mean;
}