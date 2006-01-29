#!/usr/local/bin/perl -w

# Counts items and coverage of a standard UCSC bed file

use strict;

my $lineCount=0;
my $total=0;
my $bases=0;

my $prevEnd=0;
my $prevSeqID='';

my $max=0;
my $min=99999999999;

while (<>){

  next if (/track/);
  next if (/^$/);

  my @fields=split(/\t/,$_);

  my $seqID=$fields[0];
  my $start=$fields[1];
  my $end=$fields[2];

  chomp($end);

  my $currRegion=$end-$start;

  $lineCount++;

  $total+=$currRegion;

  if (($lineCount!=1) and
	  ($prevSeqID eq $seqID) and
	  ($prevEnd>$start)){
	$bases+=($end-$prevEnd);
  } else {
	$bases+=$currRegion;
  }

  $max=($currRegion>$max ? $currRegion : $max);
  $min=($currRegion<$min ? $currRegion : $min);

  $prevEnd=$end;
  $prevSeqID=$seqID;
}

my $average=int(($total/$lineCount)+0.5);

print "Items: $lineCount\n";
print "Item bases: $bases\n";
print "Item total: $total\n";
print "Average item: $average\n";
print "Maximum item: $max\n";
print "Minimum item: $min\n";
print "\n";
