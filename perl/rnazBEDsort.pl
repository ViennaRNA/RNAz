#!/usr/local/bin/perl -w

# $Id: rnazBEDsort.pl,v 1.2 2006-03-24 15:43:13 wash Exp $

#missing: re-introduce track line, more intelligent sort of seq-ids (chr9 before chr10)

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


my %lines=();

while (my $line=<>){

  chomp($line);

  next if ($line=~/track/);
  next if ($line=~/^$/);

  my @fields=split(/\t/,$line);

  my $seqID=$fields[0];
  my $start=$fields[1];
  my $end=$fields[2];

  if (!exists($lines{$seqID})){
	$lines{$seqID}=[[@fields]];
  } else {
	push @{$lines{$seqID}},[@fields];
  }
}


foreach my $key (sort keys %lines){

  my @tmp=sort {$a->[1] <=> $b->[1]} @{$lines{$key}};

  foreach my $line (@tmp){
	#print $line->[1];
	print join("\t",@$line),"\n";
  }
}

__END__

=head1 NAME

C<rnazBEDsort.pl> - Sorts a BED annotation file.

=head1 SYNOPSIS

 rnazBEDsort.pl [options] [file]

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

C<rnazBEDsort.pl> reads a BED formatted annotation file and sorts the
lines by sequence identifier and genomic location. Note: this simple
script is not very memory efficient so you could run into problems if
you try to sort really large BEDs.

=head1 EXAMPLES

 # rnazBEDsort.pl some.bed

Sorts the file C<some.bed> and prints the results to standard out.

=head1 AUTHORS

Stefan Washietl <wash@tbi.univie.ac.at>

=cut



