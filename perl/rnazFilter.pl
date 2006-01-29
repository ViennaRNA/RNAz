#!/usr/bin/perl -w

use strict;
use FindBin;
use lib $FindBin::Bin;
use RNAz;
use Getopt::Long;
use Pod::Usage;

my $count=0;
my $help=0;


GetOptions('count'=>\$count,
		   'help'=>\$help,
		   'h'=>\$help );


my @fieldsList=qw(hitID clusterID seqID start end strand N columns
				  identity meanMFE consensusMFE energyTerm covarianceTerm combPerPair z SCI decValue P );

my %fields;

my $filter=shift @ARGV;

foreach my $field (@fieldsList){
  $filter=~s/($field)/\$fields\{$1\}/g;
}

my $nCluster=0;
my $nWindow=0;

my $currClusterLine='';

foreach my $line (<>){

  next if $line=~/^\s?\#/;
  next if $line=~/^\s+$/;

  if ($line=~/^\s?locus\d+/){
	$currClusterLine=$line;
	next;
  }

  @fields{@fieldsList}=split(/\s+/,$line);

  local $SIG{__WARN__} =
	sub {
	  print STDERR ("You have an error in your filter statement.\n");
	  exit(1);
	};
  local $SIG{__DIE__}=$SIG{__WARN__};

  my $flag=eval $filter;

  if ($flag){
	if ($currClusterLine){
	  print $currClusterLine if (!$count);
	  $currClusterLine='';
	  $nCluster++;
	}
	print $line if (!$count);
	$nWindow++;
  }
}

if ($count){
  print "Loci: $nCluster\n";
  print "Windows: $nWindow\n";
}
