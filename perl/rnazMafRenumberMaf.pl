#!/usr/bin/perl
# Last Time-stamp: <2013-08-28 18:11:27 at>
# date-o-birth: <2013-08-28 11:18:45 at>, Vienna
# renumber maf block numbering

use strict;
use warnings;





# --------------------------------------------------

my $max = 19;

#for my $i (10..$max){
for my $i (X,Y,M){

  my $c_file = "chr".$i.".maf.gtf";
  my $p_file = "chr".($i-1).".maf.gtf";
  my $o_file = "chr".$i.".maf.gtf";
  
  my $prev = `tail -1  $p_file`;

  my $c = $1 if $prev =~ /blockNr \"(\S+)\"/;
  print "blockNr $c\n";

  open(IN,  "<../$c_file") || die "could not open $c_file: $!\n";
  open(OUT, ">$o_file") || die "could not open $o_file: $!\n";
      
  while(<IN>){
    next if /^\s+$/ || /^\#/;
    chomp;
    
    my $cc = $1 if /blockNr \"(\S+)\"/;
#    print OUT "blockNr $c $cc ";
    $cc += $c;
    s/(blockNr) \"(\S+)\"/$1 \"$cc\"/;
#    print OUT "$cc";
    print OUT "$_\n";
  }
  close(OUT);
  close(IN);
}

