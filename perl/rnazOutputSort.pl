#!/usr/bin/perl -w

use strict;
use FindBin;
use lib $FindBin::Bin;
use RNAz;
use Getopt::Long;
use Pod::Usage;

my $version = 0;
my $man     = 0;
my $help    = 0;

GetOptions(
  'version' => \$version,
  'man'     => \$man,
  'v'       => \$version,
  'help'    => \$help,
  'h'       => \$help
) or pod2usage(1);

pod2usage(1) if $help;
pod2usage( -verbose => 2 ) if $man;

if ($version) {
  print "\nrnazOutputSort.pl is part of RNAz $RNAz::rnazVersion\n\n";
  print "http://www.tbi.univie.ac.at/~wash/RNAz\n\n";
  exit(0);
}

my $fileName = shift @ARGV;
my $fh;

if ( !defined $fileName ) {
  $fh = *STDIN;
} else {
  open( $fh, "<$fileName" ) || die("Could not open file $fileName ($!)");
}

my %data = ();

while ( my $rnazString = getNextRNAz($fh) ) {

  my $results = parseRNAz($rnazString);

  my ( $currStart, $currEnd, $currName, $currStrand );

  $currStart  = $results->{refSeqStart};
  $currEnd    = $results->{refSeqEnd};
  $currName   = $results->{refSeqName};
  $currStrand = $results->{refSeqStrand};

  if ( not exists $data{$currName} ) {
    $data{$currName} = [ { start => $currStart, rnazString => $rnazString } ];
  } else {
    push @{ $data{$currName} }, { start => $currStart, rnazString => $rnazString };
  }
}


foreach my $key ( keys %data ) {
  foreach my $item ( sort { $a->{start} <=> $b->{start} } @{ $data{$key} } ) {
    print "\n\n############################  RNAz ".$RNAz::rnazVersion."  ##############################\n\n";
    print $item->{rnazString};
  }
}

__END__

=head1 NAME

C<rnazOutputSort.pl> - Sorts output of RNAz by genomic coordinates
(only needed if input MAFs are unsorted)

=head1 SYNOPSIS

 rnazOutputSort.pl [options] [file]

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

C<rnazOutputSort.pl> - Sorts output of RNAz by genomic coordinates
(only needed if input MAFs are unsorted). Reads output from RNAz from
STDIN or from file given and writes the sorted output to STDOUT.

=head1 EXAMPLES

 # rnazOutputSort.pl rnaz.out

=head1 AUTHORS

Stefan Washietl <wash@mit.edu>

=cut
