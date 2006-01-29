#!/usr/bin/perl -w

use strict;
use FindBin;
use lib $FindBin::Bin;
use RNAz;
use Getopt::Long;
use Pod::Usage;


my $maxSeqs=6;
my $minSeqs=2;
my $nSamples=1;
my $maxID=99;
my $optID=80;
my $noReference=0;
my $help=0;

GetOptions('num-seqs:i' => \$maxSeqs,
		   'n:i' => \$maxSeqs,
		   'min-seqs:i' => \$minSeqs,
		   'num-samples:i' => \$nSamples,
		   'a:i' => \$nSamples,
		   'max-id:i' => \$maxID,
		   'opt-id:i' => \$optID,
		   'no-reference' => \$noReference,
		   'help'=>\$help,
		   'h'=>\$help
		  );

my $fileName=shift @ARGV;
my $fh;

if (!defined $fileName){
  $fh=*STDIN;
} else {
  open($fh,"<$fileName") || die("Could not open file $fileName ($!)");
}

my $alnFormat=checkFormat($fh);

while (my $alnString=getNextAln($alnFormat,$fh)){

  my $aln=parseAln($alnString,$alnFormat);
  #print formatAln($aln,$alnFormat);

  my $keepFirst=1;
  $keepFirst=0 if ($noReference);

  my $samples=pruneAln(alnRef=>$aln,
					   maxN=>$maxSeqs,
					   minN=>$minSeqs,
					   optSim=>$optID/100,
					   maxID=>$maxID/100,
					   keepfirst=>$keepFirst,
					   numAln=>$nSamples);

  foreach $aln (@$samples){

	print formatAln($aln,$alnFormat);
	
  }

}
