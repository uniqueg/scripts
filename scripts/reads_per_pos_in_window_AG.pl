

# perl /import/bc2/home/zavolan/grubera/Upf1/scripts/ag-make-profile-around-a-site.pl /import/bc2/home/zavolan/grubera/Upf1/PEAK /import/bc2/home/zavolan/grubera/Upf1/1553_WT_H.mappers.transcriptome 7 200 > tmp

# argument 0: /import/bc2/home/zavolan/grubera/Upf1/PEAK
# argument 1: /import/bc2/home/zavolan/grubera/Upf1/1553_WT_H.mappers.transcriptome
# argument 2: 7
# argument 3: 200
# out: > tmp

### Pre-requisites
## Pragmas
use strict;
use warnings;
## Set variables
my $peaks = "HEAD";								#$ARGV[0]; # Peak file
my $reads = "1553_WT_H.mappers.transcriptome";	#$ARGV[1]; # Read file
my $ext = 7;									#$ARGV[2]; # Set $extent to 200
my $span = 200;									#$ARGV[3]; # Set $span to 200

### Populate AoA @sites with data from '/import/bc2/home/zavolan/grubera/Upf1/PEAK'; outer array: one entry per input line; inner array: ID and site? entry in each line 
# Initialize empty AoA @sites
my @sites = ();
# Open input filehandle: '/import/bc2/home/zavolan/grubera/Upf1/PEAK'
open( F, $peaks );
## Crawl through input file
while (<F>) {
  # Remove trailing newline character
  chomp;
  # Split line by tab into array @F
  my @F = split(/\t/);
  # Push values of column 1 (ID string) and column 2 (integer) to AoA @sites
  push @sites, [ $F[0], $F[1] ];
}
# Close input filehandle
close(F);
# Print status

### Set all values of hash %profile (keys = integers between -$span and +$span (+200) to 0
# Inititate empty hash %profile
my %profile = ();
## Loop through integers $i between -$span and +$span
foreach my $i ( -$span .. $span ) {
  # Set value of each key $i of hash %profile to 0 
  $profile{$i} = 0;
}

### Populate HoA %CLIP with data from '/import/bc2/home/zavolan/grubera/Upf1/1553_WT_H.mappers.transcriptome'; keys: IDs, values: array of  start position, end position and read count (possibly multiple of 3!)
my %CLIP = %{ _readSample( $reads, \@sites, $ext ) };

###
## Crawl through each array reference $site of AoA @sites
foreach my $site (@sites) {

  # Dereference array reference and load into ID and position variables
  my ( $chromname, $position ) = @{$site};
  
  ## Set hash %map that maps the absolute to relative coordinates
  ## Initialize hash %localprofile that maps the number of reads per position to relative coordinates
  # Initialize empty hash %localprofile 
  my %localprofile = ();
  my %map          = ();
  # Initialize counter to -$span - 1
  my $cnt          = -$span - 1;			# initializes to -201??
  ## Loop through integers $i between ($position -$span) and ($position +$span) 
  foreach my $i ( ( $position - $span ) .. ( $position + $span ) ) {
    # Increase counter $cnt by 1
    $cnt++;
    # Set hash %localprofile (keys: integers between -$span and +$span; values: 0)
    $localprofile{$cnt} = 0;
    # Set hash %map (keys: integers between ($position - $span) and ($position + $span); values: integers between -$span and +$span)
    $map{$i}            = $cnt;
  }

  ## Dereference reference to array saved in %HoA %CLIP values for the current ID  
  my @reads = @{ $CLIP{$chromname} };

  # Set variable $total to 0
  my $total = 0;
  # for each element of 
  foreach my $read (@reads) {
    foreach my $p ( $read->[0] .. $read->[1] ) {
      if ( defined $map{$p} ) {
        $localprofile{ $map{$p} } += $read->[2];
        $total += $read->[2];
      }
    }
  }

  if ( $total > 0 ) {
    foreach my $key ( keys %localprofile ) {
      $localprofile{$key} = $localprofile{$key} / $total;
    }
  }


  foreach my $key ( keys %localprofile ) {
    $profile{$key} += $localprofile{$key};
  }
foreach my $key ( sort { $a <=> $b } keys %localprofile ) {
  print "$key\t$localprofile{$key}\n";
}  
}
#foreach my $key ( sort { $a <=> $b } keys %profile ) {
#  print "$key\t$profile{$key}\n";
#}


sub _readSample {

# Pass input file '/import/bc2/home/zavolan/grubera/Upf1/1553_WT_H.mappers.transcriptome' (ARGV[1])
  my $file   = $_[0];
# Pass reference to AoA @sites
  my $sites  = $_[1];
# Pass extend parameter ARGV[1] (currently 7)
  my $extend = $_[2];

# Inititate empty HoA %sample
  my %sample = ();
## Loop through each element $w of AoA @sites
  foreach my $w ( @{$sites} ) {
    # Set keys $w of hash %sample to IDs stored in the inner array of AoA @sites (element 0) and set value to empty array 
    $sample{ $w->[0] } = [];
  }

  # Open input filehandle: 'import/bc2/home/zavolan/grubera/Upf1/1553_WT_H.mappers.transcriptome'
  open( F, $file );
  ## Crawl through input file
  while (<F>) {
    # Remove trailing newline character
    chomp;
    # Split line by tab into array @F
    my @F = split(/\t/);
    # Skip line if column 2 (ID value) of the current line of $file does not exist in HoA %sample
    next if ( not defined $sample{ $F[1] } );								# BETTER: exists???
    # Push array of values of start position (column 3 - $extend), end position (column 4) and read counts (column 5) to ID key (from column 1; matches ID of AoA @sites) of HoA %sample 
    push @{ $sample{ $F[1] } }, [ $F[2] - $extend, $F[3], $F[4] ];
  }
  # Close input filehandle
  close(F);
  # return reference to HoA %sample
  return \%sample;
}
  
# my $R = '
#dataWT <- read.table("1553_WT_H.mappers.transcriptome.peak.200");
#dataTR <- read.table("1554_Puro_H.mappers.transcriptome.peak.200");
#dataEJC <- read.table("1822_HeLa_eIF4AIII_CLIP-seq2.mappers.transcriptome.peak.200");
#datamRNA <- read.table("271_EWSR1_siCTRL-A_HeLa_SetA.mappers.transcriptome.peak");
#
##pdf("test1a.pdf");
#plot(dataWT$V1,dataWT$V2/5000,type="l",ylim=c(0,0.012),col="blue",lwd=2);
#lines(dataTR$V1,dataTR$V2/5000,col="green",lwd=2);
#lines(dataEJC$V1,dataEJC$V2/5000,col="black",lwd=2);
##lines(datamRNA$V1,datamRNA$V2/5000,col="red",lwd=2);
#dev.off();
#  
#  '; 
