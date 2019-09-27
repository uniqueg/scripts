use strict;
use warnings;
use Getopt::Long;
use File::Temp qw/ tempfile tempdir /;

my $aln2seq_IOeff;
my $blastscores;
my $removeAdaptor;
my $adapter;
my $help;
my $bin = 5000000;
my $minReadLength = 15;
my $outdir = ".";
my $tempDir;

my $result = GetOptions(
  "aln2seq=s"       => \$aln2seq_IOeff,
  "blastscores=s"   => \$blastscores,
  "removalscript=s" => \$removeAdaptor,
  "a|adapter=s"       => \$adapter,
  "b|binsize=i"       => \$bin,
  "minreadlength=i"   => \$minReadLength,
  "tempdir=s"         => \$tempDir,
  "h|usage|help"            => \$help,
);

my $usage_info = 'Usage: perl ' . $0 . ' [OPTIONS] --aln2seq [FILE] --blastscores [FILE] --removalscript [FILE] --adapter [STRING|DNA SEQUENCE] [FILE|FASTA]

Description: The indicated adapter sequence and fragements thereof are removed from the 3`-ends of each of the reads in the supplied FASTA input file.

==================================================
Required arguments:
--aln2seq [FILE]	This file is required by the actual removal script supplied by --removalscript.
--blastscores [FILE]	This file is required by the actual removal script supplied by --removalscript.
--removalscript [FILE]	This is the actual removal script called by this wrapper.
--adapter [STRING]	Sequence of the adapter to be removed (DNA alphabet).
==================================================
Optional arguments:
--binsize [INT]		The input file is split up into chunks of INT entries prior to adapter removal (Default: 5000000).
--mindreadlength [INT]	Reads with a length of less than INT AFTER removal will be discarded (Default: 15).
--tempdir [PATH]	Directory in which temporary files are to be stored (Default: If available, the folder in which the input FASTA file is located, else the current working directory).
--usage|help		Show this information and die!

Comments:
- This script calls the script supplied in --removalscript to do the actual removal.

Version 1.0 (2014-02-11)
Adapted from "ag-remove-adapter" written by Andreas R. Gruber
Adapted by Alexander Kanitz, starting 2014-02-11
';
die $usage_info if $help || !$result || ( not defined $ARGV[0] ) || ( not defined $adapter ) || ( not defined $aln2seq_IOeff ) || ( not defined $blastscores ) || ( not defined $removeAdaptor );

if ( not -e $ARGV[0] ) {
  print STDERR "[ERROR] Input file '$ARGV[0]' not found.\n";
  exit 1;
}

if ( not -e $aln2seq_IOeff ) {
  print STDERR "[ERROR] Please check path to 'aln2seq_IOeff' executable.\n";
  exit 1;
}

if ( not -e $blastscores ) {
  print STDERR "[ERROR] Please check path to 'blast.scores' file.\n";
  exit 1;
}

if ( not -e $removeAdaptor ) {
  print STDERR "[ERROR] Please check path to 'removeAdaptor.pl' script.\n";
  exit 1;
}

# make a temp file
if (defined $tempDir) {
	$outdir = $tempDir;
} 
elsif ( $ARGV[0] =~ m/(.*\/).*/ ) {
	$outdir = $1;
}

# write adapter to disk
my @files = ( File::Temp->new( DIR => $outdir, SUFFIX => '.adapter' ) );
$files[0]->print("$adapter\n");
$files[0]->eof();

# convert fasta file to sol file
push @files, File::Temp->new( DIR => $outdir, SUFFIX => '.sol' );
my $cnt   = 0;
my $tmpfh = $files[$#files];
open( IN, "$ARGV[0]" ) || die "can't open pipe to $ARGV[0]";
while (<IN>) {
  chomp;
  next if ( $_ =~ m/^>/ );
  $cnt++;
  $tmpfh->print("$_\t1\n");
  if ( $cnt % $bin == 0 ) {
    $tmpfh->eof();
    push @files, File::Temp->new( DIR => $outdir, SUFFIX => ".sol" );
    $tmpfh = $files[$#files];
  }
}
close(IN);
$tmpfh->eof();

foreach my $file (@files) {
  my $in = $file->filename;
  next if ( $in !~ m/.sol$/ );
  my $out = File::Temp->new( DIR => $outdir, SUFFIX => ".sol.aln" );
  push @files, $out;
  my $outfn     = $out->filename;
  my $adapterfn = $files[0]->filename;
  `$aln2seq_IOeff $in $adapterfn $blastscores > $outfn`;
}

# remove adaptor
$cnt = 0;
foreach my $file (@files) {
  my $in = $file->filename;
  next if ( $in !~ m/.sol.aln$/ );
  open( IN2, "$removeAdaptor $in |" ) || die "can't open pipe";
  while (<IN2>) {
    chomp;
    my @F = split(/\t/);
    next if ( length( $F[0] ) < $minReadLength );
    $cnt++;
    print ">$cnt\n$F[0]\n";
  }
  close(IN2);
}

exit 0;
