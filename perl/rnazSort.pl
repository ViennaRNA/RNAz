#!/usr/bin/perl -w

use strict;
use FindBin;
use lib $FindBin::Bin;
use RNAz;
use Getopt::Long;
use Pod::Usage;


my @fieldsList=qw(hitID clusterID seqID start end strand N columns
				  identity meanMFE consensusMFE energyTerm covarianceTerm combPerPair z SCI decValue P );
my $reverse=0;
my $noClusters=0;
my $help=0;

GetOptions('reverse' => \$reverse,
		   'r' => \$reverse,
		   'no-clusters'=>\$noClusters,
		   'n'=>\$noClusters,
		   'help'=>\$help,
		   'h'=>\$help
		  );

my $sortKeyInput= shift @ARGV;
my $sortKey='';

foreach my $key (@fieldsList){
  if (lc($key) eq lc($sortKeyInput)){
	$sortKey=$key;
	last;
  }
}

if (!$sortKey){
  print STDERR "Unknown sort key.\n";
  exit(1);
}

my %clustersHeader=();
my %clusters=();
my %hits=();

while (my $line=<>){

  next if $line=~/^\s?\#/;
  next if $line=~/^\s+$/;

  if ($line=~/^(locus\d+)/){
	$clustersHeader{$1}=$line;
  }

  if ($line=~/^(window\d+)\s+(locus\d+)/){
	$clusters{$2}=[] if !defined $clusters{$2};
	push @{$clusters{$2}},$1;
	my %fields;
	@fields{@fieldsList}=split(/\s+/,$line);
	$hits{$1}={%fields};
  }
}

if (!$noClusters){

  foreach my $hitID (sort sortFunction keys %hits){
	
	next if not defined $clusters{$hits{$hitID}->{clusterID}};
	
	print $clustersHeader{$hits{$hitID}->{clusterID}};
	
	foreach my $hitInCluster (sort sortFunction @{$clusters{$hits{$hitID}->{clusterID}}}){
	  printLine($hitInCluster);
	}
	
	$clusters{$hits{$hitID}->{clusterID}}=undef;
	
  }
} else {
  foreach my $hitID (sort sortFunction keys %hits){
	printLine($hitID);
  }
}


sub sortFunction{

  if ($sortKey eq "hitID" or
	  $sortKey eq "clusterID" or
	  $sortKey eq "seqID" or
	  $sortKey eq "strand"){
	if (!$reverse){
	  return $hits{$a}->{$sortKey} cmp $hits{$b}->{$sortKey};
	} else {
	  return $hits{$b}->{$sortKey} cmp $hits{$a}->{$sortKey};
	}
  }

  if ($sortKey eq "start" or
	  $sortKey eq "end" or
	  $sortKey eq "z" or
	  $sortKey eq "meanMFE" or
	  $sortKey eq "consensusMFE" or
	  $sortKey eq "energyTerm" or
	  $sortKey eq "covarianceTerm"){
	if (!$reverse){
	  return $hits{$a}->{$sortKey} <=> $hits{$b}->{$sortKey};
	} else {
	  return $hits{$b}->{$sortKey} <=> $hits{$a}->{$sortKey};
	}
  }

  if ($sortKey eq "N" or
	  $sortKey eq "columns" or
	  $sortKey eq "identity" or
	  $sortKey eq "SCI" or
	  $sortKey eq "combPerPair" or
	  $sortKey eq "decValue" or
	  $sortKey eq "P"){
	if (!$reverse){
	  return $hits{$b}->{$sortKey} <=> $hits{$a}->{$sortKey};
	} else {
	  return $hits{$a}->{$sortKey} <=> $hits{$b}->{$sortKey};
	}
  }
  return $hits{$a}->{$sortKey}<=>$hits{$b}->{$sortKey};
}

sub printLine{
  my $hitID=shift;
  my @tmp;
  foreach my $key (@fieldsList){
	push @tmp, $hits{$hitID}->{$key};
  }
  print join("\t",@tmp);
  print "\n";
}

