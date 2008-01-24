#!/usr/bin/perl -w

# $Id: rnazMAF2BED.pl,v 1.3 2008-01-24 10:26:45 wash Exp $


#TODO: handle UCSC style a la: hg17.chr17

use strict;
use FindBin;
use lib $FindBin::Bin;
use RNAz;
use Getopt::Long;
use Pod::Usage;


my $help='';
my $seqID='';
my $version=0;
my $cluster=0;
my $man='';

GetOptions('--seq-id=s'=>\$seqID,
		   '-s=s'=>\$seqID,
		   '--cluster'=>\$cluster,
		   '-c'=>\$cluster,
		   'version'=>\$version,
		   'v'=>\$version,
		   'help'=>\$help,
		   'man'=>\$man,
		   'h'=>\$help
		  ) or pod2usage(1);

pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

if ($version){
  print "\nrnazMAF2BED.pl is part of RNAz $RNAz::rnazVersion\n\n";
  print "http://www.tbi.univie.ac.at/~wash/RNAz\n\n";
  exit(0);
}

if (!$cluster){

  while (my $line=<>){

	chomp($line);

	if ($line=~/^s\s+(.+?)\s+(\d+)\s+(\d+).*/){

	  my ($id,$start,$length)=($1,$2,$3);

	  # Use first ID in first alignment as reference if no --seq-id is
	  # given
	  if (not $seqID){
		if ($id=~/(.*)\..*/){
		  $seqID=$1;
		} else {
		  $seqID=$id;
		}
	  }

	  my $end=$start+$length;

	  if ($id=~/$seqID/){
		print "$id\t$start\t$end\n";
	  }
	}
  }
} else {

  my ($id,$start,$length);
  my ($currName,$currStart,$currEnd,$currStrand);
  my ($prevName,$prevStart,$prevEnd);
  my $minStart=99000000;
  my $maxEnd=0;
  $prevStart=0; $prevEnd=0; $prevName='';

  while (my $line=<>){

	chomp($line);

	if ($line=~/^s\s+(.+?)\s+(\d+)\s+(\d+).*/){

	  ($id,$start,$length)=($1,$2,$3);

	  # Use first ID in first alignment as reference if no --seq-id is
	  # given
	  if (not $seqID){
		if ($id=~/(.*)\..*/){
		  $seqID=$1;
		} else {
		  $seqID=$id;
		}
	  }

	  next if !($id=~/$seqID/);

	  $currName=$id;
	  $currStart=$start;
	  $currEnd=$start+$length;

	} else {

	  next;

	}

	if (!(($currName eq $prevName) and ($currStart <= $prevEnd))){

	  if ($maxEnd!=0){
		
		print "$id\t$minStart\t$maxEnd\n";
		
		$minStart=999000000;
		$maxEnd=0;

	  }
	}
	$minStart=$currStart if ($currStart<$minStart);
	$maxEnd=$currEnd if ($currEnd>$maxEnd);
	
	($prevName,$prevStart,$prevEnd)=($currName,$currStart,$currEnd);
  }
}
	

__END__

=head1 NAME

C<rnazMAF2BED.pl> - Convert sequence information from MAF formatted
multiple sequence alignment to a BED style annotation format.

=head1 SYNOPSIS

 rnazMAF2BED.pl [options] [file]

=head1 OPTIONS

=over 8

=item B<-s, --seq-id>

Specify the sequence identifier of the sequence which should be used
as a reference to create the output. Use for example C<hg17> if you
want to get all sequences containing C<hg17> in the idenitfier
(e.g. C<hg17.chr10>, C<hg17.chr22>,...). If this option is omitted,
the first sequence identifier of the first sequence in the first
alignment block is used.

=item B<-c, --cluster>

Combine overlapping alignments and report non-overlapping
regions in the BED output. 

=item B<-v, --version>

Prints version information and exits.

=item B<-h, --help>

Prints a short help message and exits.

=item B<--man>

Prints a detailed manual page and exits.

=back

=head1 DESCRIPTION

This simple programs extracts the position information for a given
sequence out of a MAF alignment and outputs it in a BED style
annotation format.

=head1 EXAMPLES

 # rnazMAF2BED.pl -s hg17 some.maf

Get the regions of the hg17 sequences in the alignment C<some.maf>.

=head1 AUTHORS

Stefan Washietl <wash@tbi.univie.ac.at>

=cut

