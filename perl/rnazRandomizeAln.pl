#!/usr/bin/perl -w

use strict;
use FindBin;
use lib $FindBin::Bin;
use RNAz;
use Getopt::Long;
use Pod::Usage;


my $level=2;
my $window=120;
my $slide=120;
my $help='';

GetOptions('window:i' => \$window,
		   'w:i' => \$window,
		   'slide:i' => \$slide,
		   's:i' => \$slide,
		   'level:i' => \$level,
		   'l:i' => \$level,
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

  my $fullAln=parseAln($alnString,$alnFormat);

  my @shuffledWindows=();
	
  my $length=length($fullAln->[0]->{seq});
  my $sliceStart=0;
  my $sliceEnd=0;
  while ($sliceStart<$length){
	$sliceEnd=$sliceStart+$window;
	$sliceEnd=$length if ($sliceEnd>$length);
	my $slice=sliceAlnByColumn($fullAln,$sliceStart,$sliceEnd);
	$slice=shuffleAln($slice,$level);
	push @shuffledWindows,$slice;
	$sliceStart+=$slide;
	last if ($sliceEnd==$length);
  }

  foreach my $entry (@{$fullAln}){
	$entry->{seq}='';
  }
  foreach my $w (@shuffledWindows){
	foreach my $i (0..@{$fullAln}-1){
	  $fullAln->[$i]->{seq}.=$w->[$i]->{seq};
	}
  }
  print formatAln($fullAln,$alnFormat);
}

