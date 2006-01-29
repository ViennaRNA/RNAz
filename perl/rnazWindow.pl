#!/usr/bin/perl -w

use strict;
use FindBin;
use lib $FindBin::Bin;
use RNAz;
use Getopt::Long;
use Pod::Usage;


my $window=120;
my $slide=40;
my $minLength=50;
my $maxGap=0.25;
my $maxMasked=0.1;
my $minN=2;
my $maxN=6;
my $refSeq=1;
my $noReference=0;
my $minID=50;
my $optID=80;
my $maxID=100;
my $forwardStrand=0;
my $reverseStrand=0;
my $bothStrands=0;
my $help=0;

GetOptions('window:i' => \$window,
		   'w:i' => \$window,
		   'slide:i' => \$slide,
		   's:i' => \$slide,
		   'min-length:i' => \$minLength,
		   'max-gap:i' => \$maxGap,
		   'max-masked:i' => \$maxMasked,
		   'max-seqs:i' => \$maxN,
		   'min-seqs:i' => \$minN,
		   'min-id:i' => \$minID,
		   'max-id:i' => \$maxID,
		   'opt-id:i' => \$optID,
		   'both-strands'=>\$bothStrands,
		   'forward'=>\$forwardStrand,
		   'reverse'=>\$reverseStrand,
		   'no-reference' => \$noReference,
		   'help'=>\$help,
		   'h'=>\$help
		  );

# If no strand is specified, default is forward
$forwardStrand=1 if (!$forwardStrand and
					 !$reverseStrand and
					 !$bothStrands);


$refSeq=0 if ($noReference);

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

  my $sliceStart=0;
  my $sliceEnd=0;
  my $length=length($fullAln->[0]->{seq});

  my $refName=$fullAln->[0]->{name};

  while ($sliceStart<$length){
	$sliceEnd=$sliceStart+$window;
	if ($sliceEnd>$length){
	  $sliceEnd=$length;
	  $sliceStart=$length-$window;
	  $sliceStart=0 if ($sliceStart<0);
	}

	my $slice=sliceAlnByColumn($fullAln,$sliceStart,$sliceEnd);

	my $sliceLength=$sliceEnd-$sliceStart;
	
	#print "BEFORE:\n";
	#print(formatAln($slice,"CLUSTAL"));

	if ($refSeq){

	  my $numGaps=($slice->[0]->{seq}=~tr/-./-/);

	  if ($numGaps/$sliceLength>$maxGap){
		$slice->[0]=undef;
		#print "Removing seq 0: too much gaps\n";
	  } else {

		for my $i (1..@$slice-1){

		  my @tmpAln=({seq=>$slice->[0]->{seq}},
					  {seq=>$slice->[$i]->{seq}});

		  removeCommonGaps(\@tmpAln);

		  my $numGaps0=($tmpAln[0]->{seq}=~tr/-./-/);
		  my $numGaps1=($tmpAln[1]->{seq}=~tr/-./-/);
		
		  my $tmpLength=length($tmpAln[0]->{seq});

		  if (($numGaps0+$numGaps1)/$tmpLength>$maxGap){
			$slice->[$i]=undef;
			#print "Removing seq $i: too much gaps\n";
		  }
		}
	  }
	} else {
	  for my $i (0..@$slice-1){
		my $numGaps=($slice->[$i]->{seq}=~tr/-./-/);
		if ($numGaps/$sliceLength>$maxGap){
		  $slice->[$i]=undef;
		  #print "Removing seq $i: too much gaps\n";
		}
	  }
	}

	for my $i (0..@$slice-1){
	  next if not defined $slice->[$i];
	  my $numMasked=($slice->[$i]->{seq}=~tr/a-z/a-z/);
	  if ($numMasked/$sliceLength>$maxMasked){
		$slice->[$i]=undef;
		#print "Removing seq $i: too much masked letters\n";
	  }
	}

	for my $i (0..@$slice-1){
	  next if not defined $slice->[$i];
	  my $tmpSeq=$slice->[$i]->{seq};
	  #print $tmpSeq, ":",rangeWarn([{seq=>$tmpSeq}]),"\n";
	  if (rangeWarn([{seq=>$tmpSeq}])){
		$slice->[$i]=undef;
		#print "Removing seq $i: out of range/too short\n";
	  }
	}
	
	my @tmp;
	foreach (@$slice){
	  next if (!defined $_);
	  push @tmp,$_;
	}
	$slice=\@tmp;
	
	# Nothing left
	goto SKIP if (!@$slice);

	# Reference sequence discarded
	goto SKIP if (($refSeq) and ($slice->[0]->{name} ne $refName));

	# Too few sequences
	goto SKIP if (@$slice<$minN);

	removeCommonGaps($slice);

	
	if ($alnFormat eq "CLUSTAL"){
		
	  for my $i (0..@$slice-1){
		$slice->[$i]->{start}=$sliceStart;
		$slice->[$i]->{end}=$sliceEnd;
		$slice->[$i]->{strand}='+';
	  }
	}

	if (@$slice > $maxN){

	  $slice=pruneAln(alnRef=>$slice,
					  maxN=>$maxN,
					  minN=>2,
					  optSim=>$optID/100,
					  maxID=>$maxID/100,
					  keepfirst=>$refSeq)->[0];
	}


	goto SKIP if length($slice->[0]->{seq})<$minLength;
	
	goto SKIP if meanPairID($slice)*100<$minID;


	my @strands=();
	
	push @strands,'+' if $forwardStrand;
	push @strands,'-' if $reverseStrand;
	
	@strands=('+','-') if $bothStrands;
	
	foreach my $strand (@strands){

	  if ($strand eq '-'){
		$slice=revAln($slice);
	  }

	  #print "\n\nAFTER:\n\n";

	  print formatAln($slice,$alnFormat);
	}
	
  SKIP:
	$sliceStart+=$slide;
	last if ($sliceEnd==$length);
  }
}



