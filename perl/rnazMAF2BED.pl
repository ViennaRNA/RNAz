#!/usr/bin/perl -w


#TODO: handle UCSC style a la: hg17.chr17

use strict;
use Getopt::Long;
use Pod::Usage;


my $help='';
my $seqID='';


GetOptions('--seq-id=s'=>\$seqID,
		   '-s=s'=>\$seqID,
		   'help'=>\$help,
		   'h'=>\$help
		  );



while (my $line=<>){

  chomp($line);

  if ($line=~/^s\s+(.+?)\s+(\d+)\s+(\d+).*/){

	my ($id,$start,$length)=($1,$2,$3);

	my $end=$start+$length;
	
	if ($id=~/$seqID/){
	  print "$id\t$start\t$end\n";
	}
  }
}

