#!/usr/bin/perl -w

use strict;
use FindBin;
use lib $FindBin::Bin;
use RNAz;
use Getopt::Long;
use Pod::Usage;

my $seqDir='';
my $blastDir='';
my $dbName='rfam';
my $cutoff='1e-06';
my $help='';

GetOptions('database=s' => \$dbName,
		   'd=s' => \$dbName,
		   'seq-dir:s' => \$seqDir,
		   's:s' => \$seqDir,
		   'blast-dir=s' => \$blastDir,
		   'b=s' => \$blastDir,
		   'help'=>\$help,
		   'h'=>\$help
		  );

if ($seqDir){
  if (not -d $seqDir){
	die("Data directory '$seqDir' not a valid ($!).");
  }
} else {
  $seqDir='.';
}

if ($blastDir){
  if (not -d $blastDir){
	die("Blast directory '$blastDir' not a valid ($!).");
  }
} else {

  $blastDir=$ENV{BLASTDB};

  if (not -d $blastDir){
	die("No blast directory found. Set the BLASTDB variable or use the XXX-command.");
  }
}
#print getSeq('/scratch2/wash/genomes/sacCer1/chr3.fa',625461,625567,'+');

#exit;

while (my $line=<>){

  # Only consider "locus" entries for annotation, simply print "windows"
  if (!($line=~/^locus/)){
	print $line;
	next;
  }

  chomp($line);

  (my $clusterID,my $chrom,my $start,my $end)=split(/\t/,$line);

  my @guess=($chrom);

  if ($chrom =~ /^(.*)\.(.*)$/){
	push @guess,$2;
  }

  my $tmp=$#guess;
  foreach my $i (0..$tmp) {
	push @guess,$guess[$i].".fa";
	push @guess,$guess[$i].".fasta";
  }

  my $realFile='';

  foreach my $file (@guess){
	if (-e "$seqDir/$file"){
	  $realFile="$seqDir/$file";
	  last;
	}
  }

  if (not $realFile){
	warn("No sequence data for '$chrom' found.\n");
	print "$line\t?\n";
  } else {
	my $seq=getSeq($realFile, $start, $end, '+');

	my @result=blastSeq($blastDir, $dbName, $cutoff,$seq);

	if (not @result){
	  print "$line\t-\n";
	} else {

	  my $name=$result[0]->{subjectID};
	  my $eValue=$result[0]->{e};
	  print "$line\t$name\n";
	}
  }
}
