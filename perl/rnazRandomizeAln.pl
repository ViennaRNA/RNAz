#!/usr/bin/perl -w

# $Id: rnazRandomizeAln.pl,v 1.2 2006/03/24 15:43:14 wash Exp $

use strict;
use FindBin;
use lib $FindBin::Bin;
use RNAz;
use Getopt::Long;
use Pod::Usage;


my $level=2;
my $window=120;
my $slide=120;
my $version=0;
my $help='';
my $man=0;

GetOptions('window:i' => \$window,
		   'w:i' => \$window,
		   'slide:i' => \$slide,
		   's:i' => \$slide,
		   'level:i' => \$level,
		   'l:i' => \$level,
		   'version'=>\$version,
		   'v'=>\$version,
		   'help'=>\$help,
		   'man'=>\$man,
		   'h'=>\$help
		  ) or pod2usage(1);

pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

if ($version){
  print "\nrnazRandomizeAln.pl is part of RNAz $RNAz::rnazVersion\n\n";
  print "http://www.tbi.univie.ac.at/~wash/RNAz\n\n";
  exit(0);
}


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

__END__

=head1 NAME

C<rnazRandomizeAln.pl> - Randomize alignments by shuffling the columns.

=head1 SYNOPSIS

 rnazRandomizeAln.pl [options] [file]

=head1 OPTIONS

=over 8

=item B<-w> N, B<--window>=N

=item B<-s> N, B<--slide>=N

Long alignment blocks should be shuffled locally in order to maintain
local characteristics of the alignment. Therefore alignments can be
shuffled in windows. You can specify here the size of a window and the
offset. Defaults are window=120 and slide=120, i.e. the alignments are
shuffled in non-overlapping windows of 120 columns.

=item B<-l> N, B<--level>=N

The shuffling algorithm tries to mantain local conservation patterns,
i.e. it shuffles only columns of the same degree of conservation. This
becomes limiting if you have many sequences in your
alignment. Therfore you can choose the level of ``coarse graining" with
this option.

To decide which columns have the same degree of conservation, the mean
pairwise identity (MPI) of each column is calculated and finally only
columns of the same value are shuffled. You can adjust the rounding of
the MPI and thus the ``coarse graining" level with this option. If you
have two columns with say 0.52 and 0.48 MPI you get:

level 0: 1 and 0

level 1: 50 and 50

level 2: 52 and 48

So on level 0 you only have ``conserved" (MPI > 0.5) and ``non-conserved"
(MPI E<lt> 0.5) columns while on level 2 you need almost exactly the same
MPI to shuffle two columns.

Default value is 2.

=item B<-v, --version>

Prints version information and exits.

=item B<-h --help>

Prints a brief help message and exits.

=item B<--man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

C<rnazRandomizeAln.pl> reads a multiple sequence alignment in Clustal
W or MAF format and returns a randomized version in the same
format. The program uses the algorithm described in Washietl &
Hofacker, J. Mol. Biol. 342(1):19 (2004). It generates alignments of
the same length, the same base composition, the same gap pattern, the
same overall conservation and the same local conservation patterns
(see also option B<--level>).

=head1 EXAMPLES

 # rnazRandomizeAln.pl -l 1 some.maf > random.maf

Randomizes the file C<some.maf> using a less stringent parameter for
maintaining conservation patterns.

=head1 AUTHORS

Stefan Washietl <wash@tbi.univie.ac.at>

=cut
