#!/usr/bin/perl -w

use strict;
use FindBin;
use lib $FindBin::Bin;
use RNAz;
use Getopt::Long;
use Pod::Usage;

my $gff=0;
my $bed=0;
my $clusters=0;
my $hits=0;
my $help='';
my $ucsc='';

GetOptions('gff' => \$gff,
		   'g' => \$gff,
		   'bed' => \$bed,
		   'b' => \$bed,
		   'clusters' => \$clusters,
		   'c' => \$clusters,
		   'hits' => \$hits,
		   'i' => \$hits,
		   'ucsc'=>\$ucsc,
		   'help'=>\$help,
		   'h'=>\$help
		  );

($clusters)=1 if (!$clusters and !$hits);
$bed=1 if (!$bed and !$gff);

my ($seqID, $name, $start, $end, $P, $strand);

$,="\t";

my $isCluster=0;
my $isHit=0;

while (my $line=<>){

  next if $line=~/^\s?\#/;
  next if $line=~/^\s+$/;

  my @fields=split(/\t/,$line);
  chomp(@fields);

  if ($line=~/^\s?locus\d+/){
	$isCluster=1;$isHit=0;
	($seqID, $start, $end, $name, $P)=
	  ($fields[1],$fields[2],$fields[3],$fields[0],$fields[6]);
	
  } elsif ($line=~/^\s?window\d+/){
	$isHit=1;$isCluster=0;
	($seqID, $start, $end, $name, $P,$strand)=
	  ($fields[2],$fields[3],$fields[4],$fields[0],$fields[16],$fields[5]);
  }
	
  if ($ucsc){
	if ($seqID=~/^(.*)\.(.*)$/){
	  $seqID=$2;
	}
	if ($P<0.5){
	  $P=0;
	} elsif (($P>=0.5) and ($P<0.9)){
	  $P=500;
	} elsif ($P>=0.9){
	  $P=1000;
	}
  }

  if ($bed){
	if ($isCluster and $clusters){
	  print $seqID,$start,$end,$name,$P,"\n";
	}
	if ($isHit and $hits){
	  print $seqID,$start,$end,$name,$P,$strand,"\n";
	}
  }

  if ($gff){
	if ($isCluster and $clusters){
	  print $seqID,"RNAz",$start+1,$end,$P,".",".","id \"$name\"","\n";
	}
	if ($isHit and $hits){
	  print $seqID,"RNAz",$start+1,$end,$P,$strand,".","id \"$name\"","\n";
	}
  }
}
