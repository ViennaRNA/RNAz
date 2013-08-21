#!/usr/bin/perl
# Last Time-stamp: <2013-08-16 18:41:51 at>
# date-o-birth: <2013-08-16 16:03:00 at>, vienna
# plot data from input alignemnts as ggplot2 R pplots/scripts

use strict;
use warnings;
use Getopt::Long;

my $tab_in;
my $bin_length=10;
my $bin_id = 5;
my $bin;

GetOptions('i:s'       => \$tab_in,
	   'bin|b'     => \$bin,
	   'binL|b:i'  => \$bin_length,
	   'binID|b:i' => \$bin_id);

my %data = ();
my @head = ();

while(<>){
  next if (/^\s+$/ || /^\#/);
  chomp;

  my @s = split(" ",$_);

  @head = @s,next unless /^\d/;
#  next unless /^\d+/;
  my ($blockLength, $meanPairID, $nrOrganisms) = @s[1,2,3];

  if ($bin){
    my $i;
    for ($i=1; $i<$blockLength; $i+=$bin_length){
      if ($blockLength >= $i && $blockLength < ($i+$bin_length)){
	
	my $j;
	for ($j=1; $j<$meanPairID; $j+=$bin_id){
	  if ($meanPairID >= $j && $meanPairID < ($j+$bin_id)){
	    $data{$nrOrganisms}{$i}{$j}++;
	  }
	}
      }
    }
  }
  else{
    print "else";
    $data{$nrOrganisms}{$blockLength}{$meanPairID}++;
  }
}




print("nrOrganisms blockLength meanPairID counts\n");
foreach my $i (sort {$a<=>$b} keys %data){
  foreach my $j (sort {$a<=>$b} keys %{$data{$i}}){
    print("$i $j $_ $data{$i}{$j}{$_}\n") foreach (sort {$a<=>$b} keys %{$data{$i}{$j}});
  }
}
