#!/usr/local/bin/perl -w

# $Id: rnazBEDstats.pl,v 1.2 2006-03-24 15:43:13 wash Exp $

use strict;
use FindBin;
use lib $FindBin::Bin;
use RNAz;
use Getopt::Long;
use Pod::Usage;
use strict;

my $help=0;
my $man=0;
my $version=0;

GetOptions('help'=>\$help,
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

my $lineCount=0;
my $total=0;
my $bases=0;

my $prevEnd=0;
my $prevSeqID='';

my $max=0;
my $min=99999999999;

while (<>){

  next if (/track/);
  next if (/^$/);

  my @fields=split(/\t/,$_);

  my $seqID=$fields[0];
  my $start=$fields[1];
  my $end=$fields[2];

  chomp($end);

  my $currRegion=$end-$start;

  $lineCount++;

  $total+=$currRegion;

  if (($lineCount!=1) and
	  ($prevSeqID eq $seqID) and
	  ($prevEnd>$start)){
	$bases+=($end-$prevEnd);
  } else {
	$bases+=$currRegion;
  }

  $max=($currRegion>$max ? $currRegion : $max);
  $min=($currRegion<$min ? $currRegion : $min);

  $prevEnd=$end;
  $prevSeqID=$seqID;
}

my $average=int(($total/$lineCount)+0.5);

print "Items: $lineCount\n";
print "Item bases: $bases\n";
print "Item total: $total\n";
print "Average item: $average\n";
print "Maximum item: $max\n";
print "Minimum item: $min\n";
print "\n";


__END__

=head1 NAME

C<rnazBEDstats.pl> - Reports some statistics on a BED annotation file.

=head1 SYNOPSIS

 rnazBEDstats.pl [options] [file]

=head1 OPTIONS

=over 8

=item B<-v, --version>

Prints version information and exits.

=item B<-h, --help>

Prints a short help message and exits.

=item B<--man>

Prints a detailed manual page and exits.

=back

=head1 DESCRIPTION

C<rnazBEDstats.pl> reads a BED formatted annotation file and prints
some basic statistics. It counts the single annotations (``items") but
also the bases covered by these items. ``Item bases" means the number
of bases that are covered by the items (overlapping regions are not
counted). ``Item total" is simply the sum of all items (overlapping
regions are counted). Important: The BED file B<must be sorted> for
this program to work. You can use C<rnazBEDsort.pl> for this task.

=head1 EXAMPLES

 # rnazBEDstats.pl some.bed

Sorts the file C<some.bed> and prints statistics for it.

=head1 AUTHORS

Stefan Washietl <wash@tbi.univie.ac.at>

=cut

