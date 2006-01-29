#!/usr/bin/perl -w

use strict;
use FindBin;
use lib $FindBin::Bin;
use RNAz;
use Getopt::Long;
use Pod::Usage;

my $BEDfile='';
my $help='';

GetOptions('file=s' => \$BEDfile,
		   'f=s' => \$BEDfile,
		   'help'=>\$help,
		   'h'=>\$help
		  );

open(BED,"<$BEDfile")||die("Could not read $BEDfile ($!)");

my %track=();

while (<BED>){
  next if /track/;
  (my $chrom, my $start, my $end, my $name)=split;

  if (!exists($track{$chrom})){
	$track{$chrom}=[{start=>$start,end=>$end,name=>$name}];
  } else {
	push @{$track{$chrom}},{start=>$start,end=>$end,name=>$name};
  }	
}

foreach my $key (keys %track){
  $track{$key}=[sort {$a->{start}<=>$b->{start}} @{$track{$key}}];
}

while (my $line=<>){

  # Only consider "cluster" entries for annotation, simply print "hits"
  if (!($line=~/\s?^locus/)){
	print $line;
	next;
  }

  (my $clusterID,my $findChrom,my $findStart,my $findEnd)=split(/\t/,$line);

  # In BED files usually only chromosome identifier are stored
  # (e.g. chr6), while in MAFs and rnazCluster.pl output you find
  # "hg17.chr6".  If the sequence idenitfiers in the original MAF are
  # of the form x.y, only y is used for comparison with the BED
  if ($findChrom =~ /^(.*)\.(.*)$/){
	$findChrom=$2;
  }

  next if (!exists($track{$findChrom}));

  # Look for two neighbouring BED entries using intervallschachtelung
  my $n1=0;
  my $n2=@{$track{$findChrom}}-1;

  while (($n2-$n1)>1){
	my $divide=$n1+int(($n2-$n1)/2);
	if (($track{$findChrom}->[$divide]->{start})<=$findStart){
	  $n1=$divide
	} else {
	  $n2=$divide;
	}
  }

  my $hit=undef;

  # if overlaps two BED entries the first is taken
  if (overlaps($findStart,$findEnd,$track{$findChrom}->[$n1]->{start},$track{$findChrom}->[$n1]->{end})){
	$hit=$n1;
  } elsif (overlaps($findStart,$findEnd,$track{$findChrom}->[$n2]->{start},$track{$findChrom}->[$n2]->{end})){
	$hit=$n2;
  }

  # If there is overlap, add the name from bed in "" as last field
  # to the cluster line
  if (defined $hit){
	my $ann=$track{$findChrom}->[$hit]->{name};
	$ann=~s/\t/ /g;
	chomp($line);
	print "$line\t\"$ann\"\n";
  } else {
	chomp($line);
	print "$line\t-\n";
  }
}

sub overlaps{
  (my $queryStart,my $queryEnd,my $subjectStart,my $subjectEnd)=@_;
  return 0 if (($queryEnd<$subjectStart) or ($queryStart>$subjectEnd));
  return 1;
}
