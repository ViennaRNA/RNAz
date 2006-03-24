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
my $version=0;
my $help=0;
my $man=0;

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
		   'man'=>\$man,
		   'h'=>\$help,
		   'version'=>\$version,
		   'v'=>\$version
		  ) or pod2usage(2);


pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

if ($version){
  print "\nrnazWindow.pl is part of RNAz $RNAz::rnazVersion\n\n";
  print "http://www.tbi.univie.ac.at/~wash/RNAz\n\n";
  exit(0);
}

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

__END__

=head1 NAME

C<rnazWindow.pl> - Slice alignments in overlapping windows and
process/filter alignment windows in various ways.

=head1 SYNOPSIS

 rnazWindow.pl [options] [file]

=head1 OPTIONS

=over 8

=item B<-w, --window>

Size of the window (Default: B<120>)

=item B<-s, --slide>

Step size (Default: B<40>)

=item B<--max-gap>=X

Maximum fraction of gaps. If a reference sequence is used
(i.e. C<--no-reference> is not set), each sequence is compared to the
reference sequence and if in the pairwise comparison the fraction of
columns with gaps is higher than X the sequence is discarded. If no
reference sequence is used, all sequences with a fraction of gaps
higher than X are discarded. (Default: B<0.25>)

=item B<--max-masked>=X

Maximum fraction of masked (=lowercase letters) in a sequence. All
sequences with a fraction of more than X lowercase letters are
discarded. This is usually used for excluding repeat sequences marked
by C<RepeatMasker> but any other information can be encoded by using
lowercase letters. (Default: B<0.1>)

=item B<--min-id>=X

Discard alignment windows with an overall mean pairwise identity
smaller than X%. (Default: B<50>)

=item B<--min-seqs>=N

Minimum number of sequences in an alignment. Discard any windows with
less than N sequences (Default:B<2>).

=item B<--max-seqs>=N

Maximum number of sequences in an alignment. If the number of
sequences in a window is higher than N, a subset of sequences is used
with exactly N sequences. The greedy algorithm of the program
C<rnazSelectSeqs.pl> is used which optimizes for a user specified mean
pairwise identity (see C<--opt-id>). (Default: B<6>)

=item B<--min-length>=N

Minimum number of columns of an alignment slice. After removing
sequences from the alignment, ``all-gap" columns are removed. If the
resulting alignment has fewer than N columns, the complete alignment
is discarded.

=item B<--opt-id>=X

If the number of sequences has to be reduced (see C<--max-seqs>) a
subset of sequences is chosen which is optimized for this value of
mean pairwise identity. (In percent, default: B<80>)

=item B<--max-id>=X

One sequence from pairs with pairwise identity higher than X % this is
removed (default: B<99>, i.e. only almost identical sequences are
removed) B<NOT IMPLEMENTED>

=item B<--forward>

=item B<--reverse>

=item B<--both-strands>

Output forward, reverse complement or both of the sequences in the
windows. Please note: C<RNAz> has the same options, so if you use
C<rnazWindow.pl> for an RNAz screen, we recommend to set the option
directly in C<RNAz> and leave the default here. (Default:
-B<--forward>)

=item B<--no-reference>

By default the first sequence is interpreted as reference
sequence. This means, for example, that if the reference sequence is
removed during filtering steps the complete alignment is
discarded. Also, if there are too many sequences in the alignment, the
reference sequence is never removed when choosing an appropriate
subset. Having a reference sequence is crucial if you are doing
screens of genomic regions. For some other applications it might not
be necessary and in such cases you can change the default behaviour by
setting this option.


=item B<-v, --version>

Prints version information and exits.

=item B<-h, --help>

Prints a short help message and exits.

=item B<--man>

Prints a detailed manual page and exits.

=back

=head1 DESCRIPTION

In many cases it is necessary to slice, pre-process and filter
alignments to get the optimal input for RNAz. This can be a tedious
task if you have a large number of alignments to analyze. This program
performs the most common pre-processing and filtering steps.

Basically it slices the input alignments (C<CLUSTAL W> or C<MAF>
format) in overlapping windows. The resulting alignments windows are
further processed and only ``reasonable" alignment windows are finally
printed out, i.e. not too much gaps/repeats, not too few or too many
sequences...

=head1 EXAMPLES

 # rnazWindow.pl --min-seqs=4 some.aln

Slices the alignment -C<some.aln> in overlapping windows of size 120,
slide 40 and filters the windows for an optimal input to RNAz
(=default behaviour). Only alignments with at least four sequences
are printed.

=head1 AUTHORS

Stefan Washietl <wash@tbi.univie.ac.at>

=cut



