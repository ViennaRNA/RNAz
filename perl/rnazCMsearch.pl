#!/usr/bin/perl -w

# $Id: rnazInfernal.pl,v 1.2 2006/03/24 15:43:13 wash Exp $

use strict;
use FindBin;
use lib $FindBin::Bin;
use RNAz;
use Getopt::Long;
use Pod::Usage;
use File::Spec qw/path catfile/;

# Set the path to your cmsearch executable here
my $cmsearchExecutable='cmsearch';

my $seqDir='';
my $cmDir='';
my $cmOpts='';
my $cmName='';
my $cutoff='1e-06';
my $version=0;
my $man=0;
my $help='';

GetOptions('seq-dir:s' => \$seqDir,
	   's:s' => \$seqDir,
	   'cm-dir=s' => \$cmDir,
	   'c=s' => \$cmDir,
	   'cmsearch-opts=s' => \$cmOpts,
	   'version'=>\$version,
	   'man'=>\$man,
	   'v'=>\$version,
	   'help'=>\$help,
	   'h'=>\$help
	   ) or pod2usage(1);

pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

if ($version){
  print "\nrnazCmsearch.pl is part of RNAz $RNAz::rnazVersion\n\n";
  print "http://www.tbi.univie.ac.at/~wash/RNAz\n\n";
  exit(0);
}


if ($seqDir){
  if (not -d $seqDir){
	die("Data directory '$seqDir' not a valid directory ($!).");
  }
} else {
  $seqDir='.';
}

if ($cmDir){
  if (not -d $cmDir){
	die("Cmsearch model directory '$cmDir' not a valid directory ($!).");
  }
} else {

#  $cmsearchDir=$ENV{CMSEARCHDB};
#
#  if (not -d $cmsearchDir){
#	die("No cmsearch directory found. Set the CMSEARCHDB variable or use the --cmsearch-dir option.");
#  }
}

my @PATH=File::Spec->path();

my $found=0;

foreach my $dir ("",@PATH){
  $found=1 if (-x File::Spec->catfile($dir, $cmsearchExecutable));
}

if (!$found){
  print STDERR "Could not find cmsearch executable. Ensure that 'cmsearch' is in your PATH or set the path directly in the perl program.\n";
  exit(2);
}

while (my $line=<>){

  # Only consider "locus" entries for annotation, simply print "windows"
  if (!($line=~/^locus/)){
	print $line;
	next;
  }

  chomp($line);

  (my $clusterID,my $chrom,my $start,my $end)=split(/\t/,$line);

  my @guess=($chrom);

  if ($chrom =~ /^(.*)\.(.*)$/){
	push @guess,$2;
  }

  my $tmp=$#guess;
  foreach my $i (0..$tmp) {
	push @guess,$guess[$i].".fa";
	push @guess,$guess[$i].".fasta";
  }

  my $realFile='';

  foreach my $file (@guess){
	if (-e "$seqDir/$file"){
	  $realFile="$seqDir/$file";
	  last;
	}
  }

  if (not $realFile){
      warn("No sequence data for '$chrom' found.\n");
      print "$line\t?\n";
  } else {
      my $seq=getSeq($realFile, $start, $end, '+');
      
      my @result=cmsearchSeq($cmDir, $cmOpts, $seq, $cmsearchExecutable);
      
      if (not @result){
	  print "$line\t-\n";
      } else {
	  print "$line\t\"", join(",", @result),"\"\n";
      }
  }
}

sub cmsearchSeq {
    my($cmDir, $cmOpts, $seq, $cmsearchExecutable) = @_;
    open(TMP,">/tmp/cmsearch$$.fa");
    print TMP ">dummy\n$seq\n";
    close(TMP);
    my $cutoff = log(1+2*length($seq))/log(2);
    my @output;
    
    my $opts = $cmOpts;
    $opts .= "-T $cutoff" unless $cmOpts =~ /-T/;
	
    foreach my $model (<$cmDir/*.cm>) {
	open(CMSEARCH, "$cmsearchExecutable $opts $model /tmp/cmsearch$$.fa|")  or die "can't run cmsearch";
	#print STDERR $model;
	my $m = $model;
	$m =~ s/\.cm$//;
	$m =~ s/$cmDir\///;
	while (<CMSEARCH>) {
	    if (/Score = (\S+),/) {
		push @output, "$m|$1" ;
	    }
	#    print STDERR $m,$_;
	}
	close(CMSEARCH);
    }
    
    return (@output);
}
__END__

=head1 NAME

C<rnazCmsearch.pl> - Compares predicted loci from data files as generated
by C<rnazCluster.pl> to a sequence database using CMSEARCH.

=head1 SYNOPSIS

 rnazCmsearch.pl [options] [file]

=head1 OPTIONS

=over 8

=item B<-b> name, B<--cmsearch-dir>=name

The directory with the covairance models for CMSEARCH. Required option.

=item B<-s> name, B<--seq-dir>=name

Directory with sequence files. For each sequence identifier in your
input file you need to have a corresponding FASTA formatted file. The
files should be named with the sequence identifier and the extension
C<.fa> or C<.fasta>. If your identifier in your input file is for
example C<contig100> then you should have a file named
C<contig100.fa>. (If your identifier is of the form
``assembly.chromosome" as for example used by UCSC alignments, it is
also possible to name the file C<chr22.fa> for a sequence identifier
C<hg17.chr22>).

=item B<--cmsearch-opts>=string

You can add additional options for cmsearch here. E.g. use
--cmsearch-opts="-T 40" to increase the score threshold to 40.
By default a score threshold of log_2(2*length(seq)) is used.

=item B<-v, --version>

Prints version information and exits.

=item B<-h --help>

Prints a brief help message and exits.

=item B<--man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

C<rnazCMsearch.pl> is a simple program to compare your hits to a sequence
database using CMSEARCH. To use it you need 
(i) a directory with covariance models (e.g. those for Rfam families)
(ii) the sequence files to which the coordinates in your results file refer
(iii) The C<cmsearch> program from the Infernal package

B<Beware that this search can take a very long time!>

Make sure that you have the sequence files available and named
correctly (see notes for the C<--seq-dir> option). In this example we
assume that the files are in the subdirectory C<seq>

You can run the following command to compare each locus in the file
C<results.dat> with each of the covariance models in the directory C<rfam>):

 # rnazCMsearch.pl --seq-dir=seq --cm-dir=rfam \
                  results.dat > annotated.dat

Any cmsearch hit the name of the matching model and the score is added in 
double quotes as additional field to the locus line.

=head1 AUTHORS

Ivo Hofackerd <ivo@tbi.univie.ac.at>

=cut

