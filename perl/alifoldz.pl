#!/usr/bin/perl -w

######################################################################
#
# alifoldz.pl --- Assessing a multiple sequence alignment for the
#                 existence of an unusual stable and conserved RNA
#                 secondary structure.
#
# Stefan Washietl <wash@tbi.univie.ac.at>
#
# Time-stamp: <04/03/10 11:31:21 wash>
######################################################################

use File::Temp;
use strict "vars";
#use Math::NumberCruncher;
use Getopt::Long;


# IMPORTANT:
#
# If RNAalifold or RNAfold is not in your path of executables,
# edit the following two variables to point to your custom locations:

my $RNAalifold="RNAalifold";
my $RNAfold="RNAfold";


# Program wide parameters which are set per command line-options

my $window=0;       # window size
my $slide=0;        # sliding offset
my $sampleN=100;    # number of random samples
my $forward;        # true if forward strand is to be scored
my $reverse;        # true if reverse strand is to be scored
my $both;           # true if both strands are to be scored
my $single='';      # if true, take first sequence in alignment and calculate
                    # calculate z-scores with RNAfold
my $threshold=-3.0; # only score segments with native MFEs lower than this
my $foldpars='';    # parameters for RNAalifold/RNAfold
my $help='';        # only display man-page.

my $random;
my $realign;
my $minID=60;

######################################################################
#
# void usage();
#
# Print manual page and exit.
#
######################################################################

sub usage{

print <<EOF;

NAME

  alifoldz.pl - Assessing a multiple sequence alignment for the
  existence of an unusual stable and conserved RNA secondary
  structure.

USAGE

   alifoldz.pl [OPTIONS] < alignment.aln

OPTIONS

  --n, -n      Number of random samples for the z-score calculation.
               Default: 100.

  --single     Score a single sequence (the first of an alignment, or
               given as a FASTA file) using RNAfold.

  --forward    Score the foreward, the reverse or both strands.
  --reverse    Default is both.
  --both

  --foldpars   The parameters for RNAalifold or RNAfold. Refer to the
               documention of RNAalifold and RNAfold for details.
               Default are default parameters.
               IMPORTANT: use quotes like this --foldpars "-T 25 -nc 5"!

  --window,-w  Score alignment using a sliding window of a specified
  --slide,-x   window size and step-size. Default is complete
               alignment.

  --threshold  Score only windows which have a native MFE below this
    -t         value. Default: -3.0

  --random     Also shuffles native sequence before calculating
    -r         the z-score. Negative control for estimating false
               positives.

  --help, -h   Display this help message.

AUTHOR

   Stefan Washietl <wash\@tbi.univie.ac.at>

EOF
exit;
}

######################################################################
#
# @shuffledAln shuffleAln(\@aln)
#
# Shuffles randomly columns of the alignment that have the same
# gap pattern and local conservation pattern.
#
# Takes reference to a array of array (rows/columns) which holds
# the alignment and returns the shuffled array.
#
######################################################################

sub shuffleAln{

  my @aln=@{$_[0]};
  my $maxRow=$#aln;
  my $maxCol=@{$aln[0]}-1;

  my @list=(); # stores mask for each column
  my %hash=(); # stores columns for each mask as hash of array
               # with mask as key

  # creates characteristic mask for each column and
  # writes @list and %hash
  foreach my $currCol (0..$maxCol){
	my %seen=();
	my $mask='';
	my $counter=0;
	foreach my $currRow (0..$maxRow){
	  my $currNt=$aln[$currRow][$currCol];
	  if ($currNt eq '-'){
		$mask.='-';
	  } else {
		if ($seen{$currNt}){
		  $mask.=$seen{$currNt};
		} else {
		  $counter++;
		  $seen{$currNt}=$counter;
		  $mask.=$counter;
		}
	  }
	}
	push @list, $mask;
	if (!exists $hash{$mask}){
	  $hash{$mask}=[$currCol];
	} else {
	  push @{$hash{$mask}},$currCol;
	}
  }

  # each list of columns with the same mask
  # are shuffled (Fisher-Yates from perlfaq)
  foreach my $arrayRef (values %hash){
	my $i;
	for ($i = @$arrayRef; --$i; ) {
	  my $j = int rand ($i+1);
	  @$arrayRef[$i,$j] = @$arrayRef[$j,$i];
	}
  }

  # columns are reassembled to a shuffled alignment
  my @shuffledAln;
  foreach my $currCol (0..$maxCol){
	my $randomCol=shift @{$hash{$list[$currCol]}};
	foreach  my $currRow (0..$maxRow){
	  $shuffledAln[$currRow][$currCol]=$aln[$currRow][$randomCol];
	}
  }
  return @shuffledAln;
}


######################################################################
#
# void writeAln(\@aln, [$filehandle])
#
# Writes alignment given as reference to array of array in ClustalW
# format to file given by $filehandle. If no $filehandle is specified
# writes to STDOUT.
#
######################################################################

sub writeAln{
  my @aln=@{$_[0]};
  my $fileHandle=$_[1];

  if (!$fileHandle){
	$fileHandle='STDOUT';
  }

  my $maxRow=$#aln;
  my $maxCol=@{$aln[0]}-1;
  my $maxName=0;
  my @names=();
  push @names, "seq$_" foreach (0..$maxRow);
  foreach my $name (@names){
	$maxName=($maxName<length($name))?length($name):$maxName;
  }
  for my $i (0..$#names){
	my $buffer=" "x(($maxName+6)-length($names[$i]));
	$names[$i].=$buffer;
  }
  my $columnWidth=60;
  my $currPos=0;
  my @seqs=();
  push @seqs, join('',@$_) foreach (@aln);
  my $length=length($seqs[0]);
  print $fileHandle "CLUSTAL W(1.81) multiple sequence alignment\n\n\n";
  while ($currPos<$length){
	for my $i (0..$#names){
	  print $fileHandle $names[$i];
	  print $fileHandle substr($seqs[$i],$currPos,$columnWidth);
	  print $fileHandle "\n";
	}
	print $fileHandle "\n\n";
	$currPos+=$columnWidth;
  }
}

######################################################################
#
# $mfe alifold(\@aln)
#
# Returns the alifold pseudo MFE for \@aln. Calls RNAalifold (given
# in the global variable $RNAalifold)
#
######################################################################

sub alifold{
  my @aln=@{$_[0]};
  (my $fh, my $fileName)=File::Temp::tempfile();
  writeAln(\@aln,$fh);
  my @result=`$RNAalifold $foldpars <$fileName 2>/dev/null`;

  if (not @result){
	print "\nERROR in RNAalifold. Maybe parameters are not correct: \"$foldpars\"\n\n";
	print "Valid options:\n";
	print "[-cv float] [-nc float]\n";
	print "[-p[0]] [-C] [-T temp] [-4] [-d[2|3]] [-noGU] [-noCloseGU]\n";
	print "[-noLP] [-e e_set] [-P paramfile] [-nsp pairs] [-S scale] [-noconv]\n\n";
	print "See documentation for RNAalifold for details.\n\n";
	exit;
  }
  my $mfe=undef;
  if ($result[1]=~/\(\s*(-?\s*\d+\.\d+)\s*=\s*-?\d+\.\d+\s*\+\s*(-?\d+\.\d+)\s*\)/){
	$mfe=$1;
  }
  unlink($fileName);
  return $mfe;
}

######################################################################
#
# $mfe rnafold(\@seq)
#
# Returns the MFE for \@seq. In this case \@seq is a reference to an
# array of the nucleotides. Calls RNAfold (given in the global
# variable $RNAalifold)
#
######################################################################

sub rnafold{
  my @aln=@{$_[0]};
  (my $fh, my $fileName)=File::Temp::tempfile();
  print $fh ">dummy\n",join('',@aln);

  my @result=`$RNAfold $foldpars<$fileName 2>/dev/null`;

  if (not @result){
	print "\nERROR in RNAfold. Maybe parameters are not correct: \"$foldpars\"\n\n";
	print "Valid options:\n";
	print "[-p[0]] [-C] [-T temp] [-4] [-d[2|3]] [-noGU] [-noCloseGU]\n";
	print "[-noLP] [-e e_set] [-P paramfile] [-nsp pairs] [-S scale] [-noconv]\n\n";
	print "See documentation for RNAfold for details.\n\n";
	exit;
  }

  my $mfe=undef;
  if ($result[2]=~/\(\s*(-?\s*\d+\.\d+)\s*\)/){
	$mfe=$1;
  }
  unlink($fileName);
  return $mfe;
}


######################################################################
#
# @slice sliceAln(\@aln, $from, $to)
#
# Returns the region $from-$to of \@aln. Use (non-biological) index 0!
#
######################################################################

sub sliceAln{
  my @aln=@{$_[0]};
  shift;
  (my $from, my $to)=@_;
  my $maxRow=$#aln;
  my @slice;
  for my $i (0..$maxRow){
	my @tmp=@{$aln[$i]};
	$slice[$i]=[@tmp[$from..$to]];
  }
  return @slice;
}

######################################################################
#
# @slice sliceSeq(\@seq, $from, $to)
#
# Returns the subsequence $from-$to of \@seq. Use (non-biological)
# index 0!
#
######################################################################

sub sliceSeq{
  my @seq=@{$_[0]};
  shift;
  (my $from, my $to)=@_;
  return @seq[$from..$to];

}

######################################################################
#
# @shuffledSeq shuffleSeq(\@seq)
#
# Shuffles the nucleotides of sequence given by \@seq. Thus, generates
# mononucleotide shuffled random sequences.
#
######################################################################

sub shuffleSeq{
  my @seq=@{$_[0]};
  for (my $i=@seq; --$i; ) {
	my $j = int rand ($i+1);
	@seq[$i,$j] = @seq[$j,$i];
  }
  return @seq;
}

######################################################################
#
# @revComplement revAln(\@aln)
#
# Returns the reverse complement of the alignment \@aln.
#
######################################################################

sub revAln{
  my @aln=@{$_[0]};
  my @rev=();
  my $maxRow=$#aln;
  for my $i (0..$maxRow){
	my $tmp=reverse @{$aln[$i]};
	$tmp=~tr/AGCTU/TCGAA/;
	$rev[$i]=[split(//,$tmp)];
  }
  return @rev;
}

######################################################################
#
# @revComplement revSeq(\@seq)
#
# Returns the reverse complement of the sequence \@seq.
#
######################################################################

sub revSeq{
  my @aln=@{$_[0]};
  my @rev=();
  my $tmp=reverse @aln;
  $tmp=~tr/AGCTU/TCGAA/;
  @rev=split(//,$tmp);
}

sub readClustal{

  my $fileName=shift;

  open(FILE,"<$fileName") or die("Could not open $fileName ($!).");

  my @aln=();
  my %order;
  my $order;
  my %alignments;

  while(<FILE>) {
	next if ( /^\s+$/ );	
	my ($seqname, $aln_line) = ('', '');	
	if( /^\s*(\S+)\s*\/\s*(\d+)-(\d+)\s+(\S+)\s*$/ ) {
	  # clustal 1.4 format
	  ($seqname,$aln_line) = ("$1/$2-$3",$4);
	} elsif( /^(\S+)\s+([A-Z\-.]+)\s*$/ ) {
	  ($seqname,$aln_line) = ($1,$2);
	} else {
	  next;
	}
	if( !exists $order{$seqname} ) {
	  $order{$seqname} = $order++;
	}
	$alignments{$seqname}.= $aln_line;
  }
  foreach my $name ( sort {$order{$a} <=> $order{$b}} keys %alignments){
	push @aln, [split(//,$alignments{$name})];
  }

  return @aln;
}




sub realignSlice{

  my @aln=@{$_[0]};

  my $length=@{$aln[0]};

  my @newAln=();


  # reference sequence is first

  for my $i (1..$#aln){
	my $match=0;
	for my $j (0..($length-1)){
	  $match++ if ($aln[$i][$j] eq $aln[0][$j]);
	}
	if (($match/$length)*100>=$minID){
	  push @newAln,$aln[$i];
	}
	
  }

	#my $seq=join('',@$ref);

	#my $numGaps=($seq=~tr/-/-/);
	#$seq=~s/-//g;

	#if (($numGaps/$length*100)<=$maxGaps){
	#  push @newAln,[split(//,$seq)];
	#}
	#print(($noGaps/$length)*100);
	#print "$seq\n";
    #}

  for my $i (0..$#newAln){
	my $seq=join('',@{$newAln[$i]});
	$seq=~s/-//g;
	$newAln[$i]=[split(//,$seq)];
  }

  #foreach my $ref (@newAln){
	#my $seq=join('',@$ref);
	#print "$seq\n";
  #}

  (my $fh, my $fileName)=File::Temp::tempfile();


  if (@newAln>=2){
	
	#writeAln(\@newAln,$fh);
	for my $i (0..$#newAln){
	  print $fh ">dummy$i\n";
	  print $fh join('',@{$newAln[$i]});
	  print $fh "\n";
	}

	system("clustalw -outorder=input -infile=$fileName >/dev/null");
	@newAln=readClustal("$fileName.aln");
	#system("cat $fileName.aln");
	unlink("$fileName.aln");
	unlink("$fileName.dnd");
  } else {
	@newAln=();
  }

  unlink("$fileName");

  return @newAln;
}



######################################################################
# MAIN PROGRAM
######################################################################

# Get user options

# handle --foldpars parameter without Getopt
my $killPos=undef;
foreach my $i (0..$#ARGV){
  #print help message also for '-help' and '?'
  if (($ARGV[$i] eq '-help') or ($ARGV[$i] eq '?')){
	usage();
  }
  if ($ARGV[$i] eq "--foldpars"){
	$foldpars=$ARGV[$i+1];
	$killPos=$i;
  }
}
# remove --foldpar "-nc 10 ..." from @ARGV ...
splice(@ARGV,$killPos,2) if (defined $killPos);

# ... and let do Getopt the rest
GetOptions ('w:i' => \$window,
			'window:i' => \$window,
			'x:i' => \$slide,
			'slide:i' => \$slide,
			'n:i'  => \$sampleN,
			'forward' => \$forward,
			'reverse' => \$reverse,
			'both' => \$both,
			'single' => \$single,
			'threshold:f' => \$threshold,
			't:f' => \$threshold,
			'random'=>\$random,
			'realign'=>\$realign,
			'minid:f'=>\$minID,
			"help"=>\$help,
			"h"=>\$help
		   );

usage() if ($help);

# interpret --forward, --reverse, --both
if ($forward and $reverse){
  $forward=1;
  $reverse=1;
} elsif ($forward){
  $forward=1;
  $reverse=0;
} elsif ($reverse){
  $forward=0;
  $reverse=1;
}
if ($both){
  $forward=1;
  $reverse=1;
}
# nothing specified -> default is both
if (!$forward and !$reverse and !$both){
  $forward=1;
  $reverse=1;
}

# both window and slide have to be specified, or nothing
if ($window and !$slide){
  print "You have specified a window but no slide. Use both --window and --slide.\n";
  exit;
}
if (!$window and $slide){
  print "You have specified a slide but no window. Use both --window and --slide.\n";
  exit;
}

# test if RNAalifold or RNAfold exist and work

if (!$single){
  # 2>&1 redirect STDERR to STDOUT, hope that works everywhere...
  my @test=`$RNAalifold -h 2>&1`; 
  if (!@test){
	warn "Could not find a working RNAalifold. Please make sure that's in your PATH,\n or specify a location by editing the source of the script $0\n";
  }
} else {
  my @test=`$RNAfold -h 2>&1`;
  if ((!@test)){
	warn "Could not find a working RNAfold. Please make sure that's in your PATH,\n or specify a location by editing the source of the script $0\n\n";
  }
}


# read input from STDIN

my @aln=();

while (my $line=<>){
  my $currSeq;
  next if ($line=~/^\s+$/);
  # starting with '> seqname", obviously fasta formatted
  if ($line=~/\s*>\s*(.*$)/){
	while ($line=<>){
	  if ($line=~/\s*>\s*(.*$)/){
		push @aln,[split(//,$currSeq)] if ($currSeq);
		$currSeq="";
		next;
	  }
	  chomp($line);
	  print $line;
	  $currSeq.=$line;
	}
	push @aln,[split(//,$currSeq)];
	last;
  }
  # a line with "Clustal" indicates a clustal format
  my %order;
  my $order;
  my %alignments;

  if ($line=~/clustal/i){
	while(<>) {
	  next if ( /^\s+$/ );	
	  my ($seqname, $aln_line) = ('', '');	
	  if( /^\s*(\S+)\s*\/\s*(\d+)-(\d+)\s+(\S+)\s*$/ ) {
		# clustal 1.4 format
		($seqname,$aln_line) = ("$1/$2-$3",$4);
	  } elsif( /^(\S+)\s+([A-Z\-.]+)\s*$/ ) {
		($seqname,$aln_line) = ($1,$2);
	  } else {
		next;
	  }
	  if( !exists $order{$seqname} ) {
		$order{$seqname} = $order++;
	  }
	  $alignments{$seqname}.= $aln_line;
	}
	foreach my $name ( sort {$order{$a} <=> $order{$b}} keys %alignments){
	  push @aln, [split(//,$alignments{$name})];
	}
  }
}

# if --single is given take first in alignment (or the only one in the fasta file)
# and remove gaps

my @single;
if ($single){
  @single=@{$aln[0]};
  my $tmp=join('',@single);
  $tmp=~s/(-|\.)//g;
  @single=split(//,$tmp);
}

# print status message

my $maxRow=$#aln;
my $maxCol=@{$aln[0]}-1;

if ($single){
  $maxCol=$#single;
}

my $from=0;
my $to;

print "\n###################################################################\n";
print "# alifoldz.pl \n#\n";

if (!$single){
print "#           Input: ", $maxRow+1," sequences of ",$maxCol+1," columns\n";
} else {
print "#           Input: 1 sequence of ",$maxCol+1," length\n";
}
print "#   Sample Number: $sampleN\n";

if ($window) {
print "#          Window: $window\n";
} else {
print "#          Window: full sequence\n";
}

if ($slide) {
print "#           Slide: $slide\n";
} else {
print "#           Slide: full sequence\n";
}

print "#          Strand: forward and reverse\n" if ($forward and $reverse);
print "#          Strand: forward\n" if ($forward and !$reverse);
print "#          Strand: reverse\n" if (!$forward and $reverse);

print "#   MFE threshold: $threshold\n";

if ($realign){
print "#    Re-alignment: exclude sequences with pairwise ID < $minID\%\n";
} else {
print "#    Re-alignment: OFF";
}
if ($random){
print "#  Random control: ON\n";
} else {
print "#  Random control: OFF\n";
}


my $program=($single?$RNAfold:$RNAalifold);
print "#    Program call: $program $foldpars\n#\n";
print "###################################################################\n\n";




# define output variables and output format

my ($native,$mean,$stdv,$z,$x,$y,$strand);

print "  From      To    Strand    Native MFE    Mean MFE     STDV      Z\n";
print " ------------------------------------------------------------------\n";

format STDOUT=
@>>>>>  @>>>>>       @     @#####.##   ^#####.##     ^##.##   ^##.#
$x,     $y,     $strand,   $native,     $mean,     $stdv,     $z
.

my $minZ=9999; #memorize the minimum Z-score

# Scan input in windows or as a whole
while ($from<$maxCol){
  $to=$from+$window-1;

  # if no --window is given scan alignment/sequence as a whole
  $to=$maxCol if ($from+$window>$maxCol or $window==0);

  $x=$from+1;
  $y=$to+1;

  my @slice;

  if (!$single){

	@slice=sliceAln(\@aln,$from,$to);

	if ($random){
	  @slice=shuffleAln(\@slice);
	}
	
	if ($realign){
	  @slice=realignSlice(\@slice);
	  if (@slice==0){
		$from+=$slide;
		next;
	  }
	}
  } else {
	@slice=sliceSeq(\@single,$from,$to);
	if ($random){
	  @slice=shuffleSeq(\@slice);
	}
  }


  if ($forward){
	if (!$single){
	  $native=alifold(\@slice);
	} else {
	  $native=rnafold(\@slice);
	}

	my @samples=();
	if ($native<=$threshold){
	  for my $i (1..$sampleN){
		if (!$single){
		  my @shuffled=shuffleAln(\@slice);
		  push @samples, alifold(\@shuffled);
		} else {
		  my @shuffled=shuffleSeq(\@slice);
		  push @samples, rnafold(\@shuffled);
		}
	  }
		
	  $mean=calc_mean(\@samples);
	  $stdv=calc_stdv(\@samples);
	  $z=($native-$mean)/$stdv;

	  $mean=sprintf("%.2f",$mean);
	  $stdv=sprintf("%.2f",$stdv);
	  $z=sprintf("%.1f",$z);

	} else {
	  $mean=undef;
	  $stdv=undef;
	  $z=undef;
	}

	if (defined $z){
	  $minZ=$z if ($z<$minZ);
	}	
	$strand="+";
	write;
  }
  if ($reverse){
	my @revSlice;

	if (!$single){
	  @revSlice=revAln(\@slice);
	} else {
	  @revSlice=revSeq(\@slice);
	}

	if (!$single){
	  $native=alifold(\@revSlice);
	} else {
	  $native=rnafold(\@revSlice);
	}

	if ($native<=$threshold){	
	  my @samples=();
	  for my $i (1..$sampleN){
		if (!$single){
		  my @shuffled=shuffleAln(\@revSlice);
		  push @samples, alifold(\@shuffled);
		} else {
		  my @shuffled=shuffleSeq(\@revSlice);
		  push @samples, rnafold(\@shuffled);
		}
	  }
	  $mean=calc_mean(\@samples);
	  $stdv=calc_stdv(\@samples);
	  $z=($native-$mean)/$stdv;

	  $mean=sprintf("%.2f",$mean);
	  $stdv=sprintf("%.2f",$stdv);
	  $z=sprintf("%.1f",$z);

	} else {
	  $mean=undef;
	  $stdv=undef;
	  $z=undef;
	 }

	if (defined $z){
	  $minZ=$z if ($z<$minZ);
	}	
	
	$strand="-";
	write;
  }
  $from+=$slide;
  last if ($to==$maxCol);
}


sub calc_mean{

  my @list=@{$_[0]};
  my $sum=0;
  $sum+=$_ foreach (@list);
  return $sum/(@list);

}

sub calc_stdv{

  my @list=@{$_[0]};

  my $mean=calc_mean(\@list);
  my $sum=0;

  foreach (@list){
        my $diff=$_-$mean;
        $sum+=($diff*$diff);
  }

  #not N-1 in denominator consistent with Number Cruncher but not
  #quite correct...
  return sqrt($sum/(@list));

}





print "\n";
print "$minZ\n";
