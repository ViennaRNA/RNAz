#!/usr/bin/perl -w
# Last Time-stamp: <2013-08-14 20:22:27 at>
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

my $print_line;

GetOptions(
  'i|index:i'    => \$index_start,
  'n:i'          => \$no_blocks,
  'o|out:s'      => \$pre,
  'line|l'       => \$print_line,
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
my $c_hs_only = 0;
my $c_all     = 0;
my $c_onecol  = 0;

# print right away
if ($print_line){
  printf "%-6s %6s %6s %6s %10s\n",
  "blockNr", "blockLength", "meanPairID", "nrOrganisms", "orgSet";
}


while ( my $alnString = getNextAln( $alnFormat, $fh ) ) {
  $c_all++;
  # name, start, length, strand, fullLength, seq, org, chrom
  my $fullAln = parseAln( $alnString, $alnFormat );
  $c_hs_only++, next if scalar(@$fullAln) <2;
      
  my %counts  = ();

  # Block length of current alignment
  my $l = length($fullAln->[0]->{seq});
  $c_onecol++, next if $l <= 1;
  $stats{blockLength}{$l}++;

  # mean pairwise ID
  my $meanPairID = meanPairID($fullAln) * 100;
  $stats{meanPairID}{$meanPairID}++;
  
  # multiple occurences
  for my $i (0..$#$fullAln){

    foreach my $j ("org"){
      $counts{$j}{$fullAln->[$i]->{$j}} = "";
    }
  }
  # count the occurences
  foreach my $k (keys %counts){
    my $c = scalar(keys %{$counts{$k}});
    $stats{$k}{$c}++;
  }

  # organisms as string
  my $org_label = join("_", (sort keys %{$counts{org}}));
  $alnCounter++;

  # print right away
  if ($print_line){
    printf "%-6i %6i %6.2f %2i %30s\n",
    $alnCounter, $l, $meanPairID, (scalar(keys %{$counts{org}})), $org_label;
  }
}

unless ($print_line){
  # print counts
  print "\# totalBlocks $c_all\n";
  print "\# hsOnly      $c_hs_only\n";
  print "\# oneCol      $c_onecol\n";
  print "\# BlocksUsed  $alnCounter\n";
  
  foreach my $k (keys %stats){
    printf("%-15s",$k);
    foreach my $i (keys %{$stats{$k}}){
      printf(" %6s %6s", $i, $stats{$k}{$i});
    }
    print "\n";
  }
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

C<mafStats.pl> - Gather some basic statistics about maf alignement.

=head1 AUTHORS

Andrea Tanzer, Stefan Washietl <at@tbi.univie.ac.at>

=cut


