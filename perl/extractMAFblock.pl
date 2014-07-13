#!/usr/bin/perl -w
# Last Time-stamp: <2014-07-13 18:32:51 at>
# date-o-birth: <2013-08-22 12:11:45 at>, vienna
# parent: rnazMAF2BED.pl
# extract maf blocks that are overlapped by some feature
# reads bedtools intersect output
# both annotation files must be gtf format

# $Id: extractMAF.pl,v 1 2014-07-13 10:26:45 at Exp $


use strict;
use FindBin;
use lib $FindBin::Bin;
use Getopt::Long;
use Pod::Usage;

use AnnotUtil;
use RNAz;



my $help    = 0;
my $man     = 0;
my $version;

my $out_format = "CLUSTAL";
my $file_size = 100;
my $prefix    = "lfold";
my $pre_mafin;


GetOptions('outformat|f:s'=> \$out_format,
	   'mafin:s'      => \$pre_mafin,
	   'out|o:s'      => \$prefix,
	   'version|v'    => \$version,
	   'help|h'       => \$help,
	   'man'          => \$man
    ) or pod2usage(1);

# Known output formats
my %out_formats = ("CLUSTAL" => undef,
		   "MAF"     => undef);

# Error message: no output format or not exsisting format specified
die("Please specify a valid output format:\n", join(", ", (sort keys %out_formats)),"\n") unless (defined(uc($out_format)) || exists($out_formats{uc($out_format)}));

# set the format
#$out_formats{$out} = 1;

pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

if ($version){
  print "\convertMAF.pl is part of RNAz $RNAz::rnazVersion\n\n";
  print "http://www.tbi.univie.ac.at/~wash/RNAz\n\n";
  exit(0);
}

# --------------------------------------------------
# Read intersect bed output
# and print gtf of gquad entry with blockNr
# --------------------------------------------------

my %data = ();

while(<>){
  next if /^\#/ || /^\s/;
  chomp;
  my ($fields1,$attribs1,$fields2,$attribs2) = splitIntersectBedGtf($_);

  my $gqID    = $attribs1->{transcript_id}->[0];
  my $blockNr = $attribs2->{blockNr}->[0];
  s/["]//g foreach ($blockNr, $gqID);

  # calculate the file in which the maf block resides
  # $file_size is the number of maf blocks per file
  # calculation is based on convertion of alignment number
  # within a file to block nr sithin a chromosome
  # see processMafGtf.pl
  # my $shift = $file_size * ($off-1);
  # off is file name
  # example: file 1   blockNr 1..100
  #               2           101-200
  # file = int((100+ (100-1))/100 ) = 1
  # file pos = 100 - 100*(1-1) = 200
  
  my $fileNr  = int(($blockNr+($file_size-1))/$file_size);
  my $filepos = $blockNr - $file_size*($fileNr-1);

  $data{$fileNr}{$filepos} = $blockNr;
#  push(@{$data{$fileNr}{$filepos}[1]}, $gqID);

  # print gtf file

#  $attribs1->{blockNr}->[0] = $blockNr;
  
}



# --------------------------------------------------
# Print alignments
# --------------------------------------------------

foreach my $fileNr (keys %data ){

  my @filepos = (sort {$a <=> $b} keys %{$data{$fileNr}});

# --------------------------------------------------
# Open the determined maf file
# --------------------------------------------------

  my $fileName = join(".",$pre_mafin,$fileNr,"maf");
  my $fh;
  
  if ( !defined $fileName ) {
    $fh = *STDIN;
  }
  else{
    if ( $fileName =~ /\.gz$/ ) {
      open( $fh, "zcat $fileName |" ) || die("Could not open file $fileName ($!)");
    }
    else{
      open( $fh, "<$fileName" ) || die("Could not open file $fileName ($!)");
    }
  }
    

  # Some Variables
  my $alnCounter = 0;
  
  my $alnFormat = checkFormat($fh); # works only with MAF for now!
  #print "format: $alnFormat\n";
  
  while ( my $alnString = getNextAln( $alnFormat, $fh ) ) {
    
    # name, start, length, strand, fullLength, seq, org, chrom
    my ($fullAln, $orgs) = parseAln( $alnString, $alnFormat, "1" );
    
    # only reference species in alignment
    next if scalar(@$fullAln) < 2;
    
    $alnCounter++;
    print "$alnCounter\n";
    
    # Print alignment in clustal format
    if ($filepos[0] == $alnCounter){

      my $out = join(".",$prefix,$data{$fileNr}{$filepos[0]},"aln");
      open(BLOCK, ">$out") || die "could not open output file $out $!\n";

      print BLOCK (formatAln($fullAln,uc($out_format)));

      close(BLOCK);

      shift @filepos;
    }
    last unless defined($filepos[0]);
  }
  close($fh);
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

