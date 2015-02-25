#!/usr/bin/perl -w

use strict;
use FindBin;
use lib $FindBin::Bin;
use RNAz;
use Getopt::Long;
use Pod::Usage;

my $window        = 120;
my $slide         = 40;
my $maxLength     = undef;
my $minLength     = 50;
my $maxGap        = 0.25;
my $maxMasked     = 0.1;
my $minN          = 2;
my $maxN          = 6;
my $numSamples    = 1;
my $refSeq        = 1;
my $noReference   = 0;
my $minID         = 50;
my $optID         = 80;
my $maxID         = 100;
my $forwardStrand = 0;
my $reverseStrand = 0;
my $bothStrands   = 0;
my $verbose       = 0;
my $version       = 0;
my $help          = 0;
my $man           = 0;
my $noRangeWarn   = 0;

GetOptions(
  'window:i'      => \$window,
  'w:i'           => \$window,
  'slide:i'       => \$slide,
  's:i'           => \$slide,
  'm:i'           => \$maxLength,
  'max-length:i'  => \$maxLength,
  'min-length:i'  => \$minLength,
  'max-gap:f'     => \$maxGap,
  'max-masked:f'  => \$maxMasked,
  'max-seqs:i'    => \$maxN,
  'num-samples:i' => \$numSamples,
  'min-seqs:i'    => \$minN,
  'min-id:f'      => \$minID,
  'max-id:f'      => \$maxID,
  'opt-id:f'      => \$optID,
  'both-strands'  => \$bothStrands,
  'forward'       => \$forwardStrand,
  'reverse'       => \$reverseStrand,
  'no-reference'  => \$noReference,
  'verbose'       => \$verbose,
  'help'          => \$help,
  'man'           => \$man,
  'h'             => \$help,
  'version'       => \$version,
  'v'             => \$version,
  'no-rangecheck'  => \$noRangeWarn
) or pod2usage(2);

pod2usage(1) if $help;
pod2usage( -verbose => 2 ) if $man;

if ($version) {
  print "\nrnazWindow.pl is part of RNAz $RNAz::rnazVersion\n\n";
  print "http://www.tbi.univie.ac.at/~wash/RNAz\n\n";
  exit(0);
}

# If no strand is specified, default is forward
$forwardStrand = 1 if ( !$forwardStrand
  and !$reverseStrand
  and !$bothStrands );

$refSeq = 0 if ($noReference);

$maxLength = $window if not defined $maxLength;

my $originalWindow = $window; # Save command line option because
                              # $window is changed for small
                              # alignments dynamically

my $fileName = shift @ARGV;
my $fh;

if ( !defined $fileName ) {
  $fh = *STDIN;
} else {
  open( $fh, "<$fileName" ) || die("Could not open file $fileName ($!)");
}

my $alnFormat = checkFormat($fh);
# print "format: $alnFormat\n";
my $alnCounter = 0;

while ( my $alnString = getNextAln( $alnFormat, $fh ) ) {

  $window = $originalWindow; # reset $window size in case it has been altered 

  $alnCounter++;

  #print STDERR "Processing aln $alnCounter\n";

  my $fullAln = parseAln( $alnString, $alnFormat );

  my @tmp = ();
  foreach (@$fullAln) {
    push @tmp, { %{$_} };
  }
  my $shrinkingAln = [@tmp];

  my $sliceStart     = 0;
  my $prevSliceStart = 0;
  my $sliceEnd       = 0;
  my $length         = length( $fullAln->[0]->{seq} );

  if ( $length <= $maxLength ) {
    $window = $length;
  }

  my $refName = $fullAln->[0]->{name};

  my $windowCounter = 0;

  while ( $sliceStart < $length ) {

    $windowCounter++;

    #print STDERR "Processing window $windowCounter\n";

    $sliceEnd = $sliceStart + $window;
    if ( $sliceEnd > $length ) {
      $sliceEnd   = $length;
      $sliceStart = $length - $window;
      $sliceStart = 0 if ( $sliceStart < 0 );
    }

    #my $slice=sliceAlnByColumn($fullAln,$sliceStart,$sliceEnd);

    # correct ends without warning if outside of valid range
    $sliceStart = 0 if ( $sliceStart < 0 );
    $sliceEnd = length( $fullAln->[0]->{seq} ) if ( $sliceEnd > length( $fullAln->[0]->{seq} ) );

    # make deep copy of list of hash
    my @newAln = ();
    foreach (@$fullAln) {
      push @newAln, { %{$_} };
    }

    #print "prev: $prevSliceStart, $sliceStart\n";

    foreach my $i ( 0 .. $#newAln ) {

      if ( ( defined $newAln[$i]->{start} ) and ( defined $newAln[$i]->{start} ) ) {
        my $oldStart = $newAln[$i]->{start};
        my $oldEnd   = $newAln[$i]->{end};

        #$newAln[$i]->{start}=alnCol2genomePos($newAln[$i]->{seq},$oldStart,$sliceStart);
        #$newAln[$i]->{end}=alnCol2genomePos($newAln[$i]->{seq},$oldStart,$sliceEnd-1)+1;

        $newAln[$i]->{start} = alnCol2genomePos(
          $shrinkingAln->[$i]->{seq},
          $shrinkingAln->[$i]->{start},
          $sliceStart - $prevSliceStart, 'after'
        );

        $newAln[$i]->{end} = alnCol2genomePos(
          $shrinkingAln->[$i]->{seq},
          $shrinkingAln->[$i]->{start},
          $sliceEnd - $prevSliceStart - 1, 'before'
        ) + 1;

      }

      $newAln[$i]->{seq} = substr( $newAln[$i]->{seq}, $sliceStart, $sliceEnd - $sliceStart );

    }

    my $slice = [@newAln];

    my $sliceLength = $sliceEnd - $sliceStart;

    foreach my $i ( 0 .. @{$shrinkingAln} - 1 ) {

      $shrinkingAln->[$i]->{seq} =
        substr( $fullAln->[$i]->{seq}, $sliceStart, $length - $sliceStart );
      $shrinkingAln->[$i]->{start} = $slice->[$i]->{start};
    }

    #	print "BEFORE:\n";
    #	print(formatAln($slice,"CLUSTAL"));
    
    if ($refSeq) {
      my $numGaps = ( $slice->[0]->{seq} =~ tr/-./-/ );

      if ( $numGaps==$sliceLength or $numGaps / $sliceLength > $maxGap ) {
        $slice->[0] = undef;

        if ($verbose) {
          print STDERR
            "Alignment $alnCounter, window $windowCounter: Removing seq 1: too many gaps.\n";
        }
      } else {

        for my $i ( 1 .. @$slice - 1 ) {

          my @tmpAln = ( { seq => $slice->[0]->{seq} }, { seq => $slice->[$i]->{seq} } );

          removeCommonGaps( \@tmpAln );

          my $numGaps0 = ( $tmpAln[0]->{seq} =~ tr/-./-/ );
          my $numGaps1 = ( $tmpAln[1]->{seq} =~ tr/-./-/ );

          my $tmpLength = length( $tmpAln[0]->{seq} );

          if ($tmpLength==0 or ( $numGaps0 + $numGaps1 ) / $tmpLength > $maxGap ) {
            $slice->[$i] = undef;
            if ($verbose) {
              my $ii = $i + 1;
              print STDERR
                "Alignment $alnCounter, window $windowCounter: Removing seq $ii: too many gaps.\n";
            }
          }
        }
      }
    } else {
      for my $i ( 0 .. @$slice - 1 ) {
        my $numGaps = ( $slice->[$i]->{seq} =~ tr/-./-/ );
        if ( $numGaps / $sliceLength > $maxGap ) {
          $slice->[$i] = undef;
          if ($verbose) {
            my $ii = $i + 1;
            print STDERR
              "Alignment $alnCounter, window $windowCounter: Removing seq $ii: too many gaps.\n";
          }
        }
      }
    }

    for my $i ( 0 .. @$slice - 1 ) {
      next if not defined $slice->[$i];
      my $numMasked = ( $slice->[$i]->{seq} =~ tr/a-z/a-z/ );
      if ( $numMasked / $sliceLength > $maxMasked ) {
        $slice->[$i] = undef;
        if ($verbose) {
          my $ii = $i + 1;
          print STDERR
            "Alignment $alnCounter, window $windowCounter: Removing seq $ii: too many masked letters.\n";
        }
      }
    }

    for my $i ( 0 .. @$slice - 1 ) {
      next if not defined $slice->[$i];
      my $tmpSeq = $slice->[$i]->{seq};

      #print $tmpSeq, ":",rangeWarn([{seq=>$tmpSeq}]),"\n";
      my $warning = rangeWarn( [ { seq => $tmpSeq } ] );
      if ($warning) {
        $slice->[$i] = undef unless $noRangeWarn;
        if ($verbose) {
          my $ii = $i + 1;
          if ( $warning == 1 ) {
            print STDERR
              "Outside training range: Alignment $alnCounter, window $windowCounter: Seq $ii: too short or too long.\n";
          }
          if ( $warning == 2 ) {
            print STDERR
              "Outside training range: Alignment $alnCounter, window $windowCounter: Seq $ii: base composition out of range.\n";
          }
          if ( $warning == 3 ) {
            print STDERR
              "Outside training range: Alignment $alnCounter, window $windowCounter: Seq $ii: base composition out of range/too short or long.\n";
          }
          print STDERR "Removing Seq $ii\n" unless $noRangeWarn;
        }
      }
    }

    my @tmp;
    foreach (@$slice) {
      next if ( !defined $_ );
      push @tmp, $_;
    }
    $slice = \@tmp;

    # Nothing left
    if ( !@$slice ) {
      if ($verbose) {
        print STDERR "Alignment $alnCounter discarded: No sequences left.\n";
      }
      goto SKIP;
    }

    # Reference sequence discarded
    if ( ($refSeq) and ( $slice->[0]->{name} ne $refName ) ) {
      if ($verbose) {
        print STDERR
          "Alignment $alnCounter discarded: Reference sequence was discarded in previous filter steps.\n";
      }
      goto SKIP;
    }

    # Too few sequences
    if ( @$slice < $minN ) {
      if ($verbose) {
        print STDERR "Alignment $alnCounter discarded: Too few sequences left.\n";
      }
      goto SKIP;
    }

    removeCommonGaps($slice);

    if ( $alnFormat eq "CLUSTAL" ) {

      for my $i ( 0 .. @$slice - 1 ) {
        $slice->[$i]->{start}  = $sliceStart;
        $slice->[$i]->{end}    = $sliceEnd;
        $slice->[$i]->{strand} = '+';
      }
    }

    my $slices;

    if ( @$slice > $maxN ) {

      $slices = pruneAln(
        alnRef    => $slice,
        maxN      => $maxN,
        minN      => 2,
        optSim    => $optID / 100,
        maxID     => $maxID / 100,
        numAln    => $numSamples,
        keepfirst => $refSeq
      );
    } else {

      $slices = [$slice];

    }

    foreach my $slice (@$slices) {

      if ( length( $slice->[0]->{seq} ) < $minLength ) {
        if ($verbose) {
          print STDERR "Alignment $alnCounter discarded: Too short.\n";
        }
        next;
      }

      if ( meanPairID($slice) * 100 < $minID ) {
        if ($verbose) {
          print STDERR "Alignment $alnCounter discarded: Mean pairwise identity out of range.\n";
        }
        next;
      }

      my @strands = ();

      push @strands, '+' if $forwardStrand;
      push @strands, '-' if $reverseStrand;

      @strands = ( '+', '-' ) if $bothStrands;

      foreach my $strand (@strands) {

        if ( $strand eq '-' ) {
          $slice = revAln($slice);
        }

        #print "\n\nAFTER:\n\n";

        print formatAln( $slice, $alnFormat );
      }
    }

  SKIP:
    $prevSliceStart = $sliceStart;
    $sliceStart += $slide;
    last if ( $sliceEnd == $length );
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

=item B<-w, --window>=N

Size of the window (Default: B<120>)

=item B<-s, --slide>=N

Step size (Default: B<120>)

=item B<-m, --max-length>

Slice only alignments longer than N columns. This means blocks longer
than the window size given by B<--window> but shorter than N are kept
intact and not sliced. Per default this length is set to the window
size given by B<--window> (or 120 by default).

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

=item B<--num-samples>=N

Number of different subsets of sequences that is sampled if there are
more sequences in the alignment than C<--max-seqs>. (Default: B<1>)

=item B<--min-length>=N

Minimum number of columns of an alignment slice. After removing
sequences from the alignment, ``all-gap" columns are removed. If the
resulting alignment has fewer than N columns, the complete alignment
is discarded. Default: 50

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

=item B<--no-rangecheck>

By default, all sequences of all windows are discarded, if the length 
or base composition is outside the training range of RNAz, independent 
of the window-size commandline parameter. 
However with the flag --no-rangecheck such sequences outside the training 
range are not discarded. As of version 2.0, RNAz can cope with sequences 
outside this traning range. However the same quality of the RNAz results 
cannot be guaranteed if sequences outside the default range are present.

=item B<--verbose>

Verbose output on STDERR, describing all performed filtering steps.

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


