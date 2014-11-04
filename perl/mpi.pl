#!/usr/bin/perl -w

use strict;
use FindBin;
use lib $FindBin::Bin;
use RNAz;
use Getopt::Long;
use Pod::Usage;
use List::Util qw(sum);

my $fileName = shift @ARGV;
my $fh;

sub mean {
    my @arr=@{$_[0]};
    return sum(@arr)/@arr;
}

if ( !defined $fileName ) {
  $fh = *STDIN;
} else {
  open( $fh, "<$fileName" ) || die("Could not open file $fileName ($!)");
}

my $alnFormat = checkFormat($fh);
# print "format: $alnFormat\n";

while ( my $alnString = getNextAln( $alnFormat, $fh ) ) {
  print "======================================\n";
  my $alntmp  = parseAln( $alnString, $alnFormat );  
  my @aln=@{$alntmp};  
  my $N = @aln;
  for my $i (0..$N-1){
    for my $j ($i..$N-1){
       next if $i==$j;
      my $pi=meanPairID( [$aln[$i],$aln[$j]] );
      $aln[$i]->{name}=~/(.*?)\..*/;   
      my $name1=$1;
      $aln[$j]->{name}=~/(.*?)\..*/;   
      my $name2=$1;
      print "Pairwise identity between $name1 and $name2 is $pi\n";
    }
  }
}


