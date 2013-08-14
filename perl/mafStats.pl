#!/usr/bin/perl -w
# Last Time-stamp: <2013-08-14 14:14:16 at>
# date-o-birth: <2013-07-22 13:36:59 at>
# mafStats.pl
# does some basic statistics on maf alignments
# write N maf blocks per file (default: 1)

use strict;
use RNAz;
use Getopt::Long;
use Pod::Usage;

my $version       = 0;
my $help          = 0;
my $man           = 0;

my $index_start   = 1;
my $no_blocks     = 1;
my $pre           = "mafsplit";

GetOptions(
  'i|index:i'    => \$index_start,
  'n:i'          => \$no_blocks,
  'o|out:s'      => \$pre,
  'h|help'       => \$help,
  'man'          => \$man,
  'v|version'    => \$version
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -verbose => 2 ) if $man;

if ($version) {
  print "\nrnazWindow.pl is part of RNAz $RNAz::rnazVersion\n\n";
  print "http://www.tbi.univie.ac.at/~wash/RNAz\n\n";
  exit(0);
}


my $fileName = shift @ARGV;
my $fh;

if ( !defined $fileName ) {
  $fh = *STDIN;
} else {
  open( $fh, "<$fileName" ) || die("Could not open file $fileName ($!)");
}

my $alnFormat = checkFormat($fh);

my $alnCounter = 0;
my $filecount  = $index_start;
my @tmp_aln    = ();
my %stats = ();

while ( my $alnString = getNextAln( $alnFormat, $fh ) ) {

  $alnCounter++;
  # name, start, length, strand, fullLength, seq, org, chrom
  my $fullAln = parseAln( $alnString, $alnFormat );

  my %counts  = ();
  my %counts2 = ();
  
  for my $i (0..$#$fullAln){

    # multiple occurences
    foreach my $j ("org"){
      $counts{$j}{$fullAln->[$i]->{$j}}="";
    }
  }

  # Block length of current alignment
  my $l = length($fullAln->[0]->{seq});
  $stats{blockLength}{$l}++;

  # mean pairwise ID
  my $meanPairID = meanPairID($fullAln) * 100;
  $stats{meanPairID}{$meanPairID}++;
  
  # count the occurences
  foreach my $k (keys %counts){
    my $c = scalar(keys %{$counts{$k}});
    $stats{$k}{$c}++;
  }
}

# print counts
print "\# totalBlocks $alnCounter\n";

foreach my $k (keys %stats){
  printf("%-15s",$k);
  foreach my $i (keys %{$stats{$k}}){
    printf(" %6s %6s", $i, $stats{$k}{$i});
  }
  print "\n";
}

  
  
#  $alnCounter++;
#  push(@tmp_aln, $alnString);
#  
#  if (($alnCounter == $no_blocks) || eof()){
#
#    my $out = join(".", $pre,$filecount,lc($alnFormat));
#    open(OUT,">$out") || die "could not open rnazMafSplit outfile $out: $!\n";
#
#    print OUT "a score=0\n$_\n" foreach @tmp_aln;
#    close(OUT);
#
#    $filecount++;
#
#    # Clean up
#    @tmp_aln    = ();
#    $alnCounter = 0;
#  }


__END__

=head1 NAME

C<rnazMafSplit.pl> - Split Maf alignemnt file into individual
maf blocks and write each block into a separate file.

Output files are named as follows:

  mafsplit.INDEX.ALNFORMAT

=head1 SYNOPSIS

 rnazMafSplit.pl [options] [file]

=head1 OPTIONS

=over 8

=item B<-n> N

number of maf blocks per file (Default: B<1>)

=item B<-i, --index> N

start indexing blocks at index i (Default: B<1>)

=item B<-o, --out>

prefix of output files. (Default: B<mafsplit>)

=item B<-v, --version>

Prints version information and exits.

=item B<-h, --help>

Prints a short help message and exits.

=back

=head1 DESCRIPTION

To run RNAz on large data sets, e.g. genome wide alignments, it is
useful to split the input data, in this case genome wide alignemnst in maf format,
and run the analysis on a cluster machine.
    
=head1 EXAMPLES

 # rnazMafSplit.pl some_alignment.maf

Split the alignment some_alignment.maf and write each maf block into
a aseparate file. Files will be named by a running index,
which corresponds to the number of the maf block in the maf file.

=head1 AUTHORS

Andrea Tanzer, Stefan Washietl <at@tbi.univie.ac.at>

=cut


