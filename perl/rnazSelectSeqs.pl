#!/usr/bin/perl -w

# $Id: rnazSelectSeqs.pl,v 1.3 2006/03/24 15:43:14 wash Exp $

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
my $version=0;
my $man=0;


GetOptions('num-seqs:i' => \$maxSeqs,
		   'n:i' => \$maxSeqs,
		   'min-seqs:i' => \$minSeqs,
		   'num-samples:i' => \$nSamples,
		   'a:i' => \$nSamples,
		   'max-id:i' => \$maxID,
		   'opt-id:i' => \$optID,
		   'i:i' => \$optID,
		   'no-reference' => \$noReference,
		   'x' => \$noReference,
		   'help'=>\$help,
		   'man'=>\$man,
		   'h'=>\$help,
		   'version'=>\$version,
		   'v'=>\$version,
		  ) or pod2usage(1);

pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

if ($version){
  print "\nrnazSelectSeqs.pl is part of RNAz $RNAz::rnazVersion\n\n";
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
__END__

=head1 NAME

C<rnazSelectSeqs.pl> - Select subsets of sequences from an alignment.

=head1 SYNOPSIS

 rnazSelectSeqs.pl [options] [file]

=head1 OPTIONS

=over 8

=item B<-n> N, B<--num-seqs>=N

Number of sequences in the output alignment(s). (Default:B<6>)

=item B<-a> N, B<--num-samples>=N

Number of output alignments (Default: B<1>)

=item B<-i> X, B<--opt-id>=X

The resulting alignment(s) is (are) optimized for this value of mean
pairwise identity (in percent, default: B<80>)

=item B<--max-id>=X

Sequences from pairs with pairwise identity higher than X% are removed
(default: B<99>, i.e. only almost identical sequences are removed)

=item B<-x, --no-reference>

By default the first sequence (=reference sequence) is always present
in the output alignment(s). If you do not care having it removed, set
this flag.

=item B<-v, --version>

Prints version information and exits.

=item B<-h, --help>

Prints a short help message and exits.

=item B<--man>

Prints a detailed manual page and exits.

=back

=head1 DESCRIPTION

C<rnazSelectSeqs.pl> reads a multiple sequence alignment in C<Clustal
W> or C<MAF> format and returns an alignment in the same format with a
user specified number of sequences. The subset is greedily optimized
for a user specified mean pairwise identity. There are options to
removes sequences which are too similar. It is also possible to sample
more than one alignment. The program uses a simple heuristic to
accomplish that.

=head1 EXAMPLES

 # rnazSelectSeqs.pl -n 4 -a 3 miRNA.maf

Samples three subsets of four sequences from the alignment C<miRNA.maf>.

 # rnazSelectSeqs.pl -n 5 -i 70 miRNA.maf

Selects a subset of five sequences optimized to a mean pairwise identity of 70%.

=head1 AUTHORS

Stefan Washietl <wash@tbi.univie.ac.at>

Ivo Hofacker <ivo@tbi.univie.ac.at>

=cut
