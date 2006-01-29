#!/usr/bin/perl -w

use strict;
use FindBin;
use lib $FindBin::Bin;
use RNAz;
use Getopt::Long;
use Pod::Usage;

my $cutoff=0.5;
my $printHits=0;
my $printClusters=0;
my $printHeader=0;
my $html=0;
my $help=0;

GetOptions('cutoff:f' => \$cutoff,
		   'c:f' => \$cutoff,
		   'hits'=>\$printHits,
		   'clusters'=>\$printClusters,
		   'i'=>\$printHits,
		   'l'=>\$printClusters,
		   'header'=>\$printHeader,
		   'd'=>\$printHeader,
		   'html'=>\$html,
		   'help'=>\$help,
		   'h'=>\$help
		  );

if (!$printHits and !$printClusters){
  $printHits=$printClusters=1;
}

if ($printHeader){

  if ($printClusters){
	print "# clusterID\tsequenceID\tstart\tend\tmaxN\tmaxIdentity\tmaxP\tminZ\n";
  }

  if ($printHits){
	print "# hitID\tclusterID\tsequenceID\tstart\tend\tstrand\tN\tcolumns\tidentity\tmeanMFE\tconsensusMFE\t'energy'\t'covariance'\tz\tSCI\tdecValue\tP\n";
  }
}

if ($html){
  mkdir "results" || die("Could not create directory for HTML files ($!)");
}


my ($currName,$currStart,$currEnd,$currStrand);
my ($prevName,$prevStart,$prevEnd);

$prevName=''; $prevStart=0; $prevEnd=0;

my ($maxP, $maxN, $maxID, $maxZ);
$maxN=0;$maxID=0;$maxZ=100;$maxP=0;
my $minStart=99000000;
my $maxEnd=0;

my @hits=();

my $clusterID=1;
my $hitID=1;

my $fileName=shift @ARGV;
my $fh;

if (!defined $fileName){
  $fh=*STDIN;
} else {
  open($fh,"<$fileName") || die("Could not open file $fileName ($!)");
}

while (my $rnazString=getNextRNAz($fh)){

  my $results=parseRNAz($rnazString);

  $currStart=$results->{aln}->[0]->{start};
  $currEnd=$results->{aln}->[0]->{end};
  $currName=$results->{aln}->[0]->{name};
  $currStrand=$results->{aln}->[0]->{strand};

  if ($results->{P}>$cutoff){
	
	# if there is NO overlap start new cluster
	if (!(($currName eq $prevName) and ($currStart <= $prevEnd))){

	  if ($maxEnd!=0){ # omit first one
		if ($printClusters){
		  print "locus$clusterID\t$prevName\t$minStart\t$maxEnd\t$maxN\t$maxID\t$maxP\t$maxZ\n";
		}

		if ($printHits){
		  foreach (@hits){
			s/CLUSTER/locus$clusterID/;
			print;
		  }
		}
		@hits=();
		$maxN=0;
		$maxID=0;
		$maxZ=100;
		$maxP=0;
		$minStart=999000000;
		$maxEnd=0;

		$clusterID++;
	  }
	}

	$maxN=$results->{N} if ($results->{N}>$maxN);
	$maxP=$results->{P} if ($results->{P}>$maxP);
	$maxID=$results->{identity} if ($results->{identity}>$maxID);
	$maxZ=$results->{z} if ($results->{z}<$maxZ);
	$minStart=$currStart if ($currStart<$minStart);
	$maxEnd=$currEnd if ($currEnd>$maxEnd);
	
	my $tmp="window$hitID\tCLUSTER\t$currName\t$currStart\t$currEnd\t$currStrand";
	foreach my $key (qw(N columns identity meanMFE consensusMFE energy covariance combPerPair z sci decValue P)){
	  $tmp.="\t$results->{$key}";
	}
	push @hits,"$tmp\n";

	
	
	$hitID++;
	
	($prevName,$prevStart,$prevEnd)=($currName,$currStart,$currEnd);

  }
}

# don't forget last cluster...

if (@hits){
  if ($printClusters){
	print "locus$clusterID\t$prevName\t$minStart\t$maxEnd\t$maxN\t$maxID\t$maxP\t$maxZ\n";
  }

  if ($printHits){
	foreach (@hits){
	  s/CLUSTER/locus$clusterID/;
	  print;
	}
  }
}




