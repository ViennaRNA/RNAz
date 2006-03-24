#!/usr/bin/perl -w

# $Id: rnazBlast.pl,v 1.2 2006-03-24 15:43:13 wash Exp $

use strict;
use FindBin;
use lib $FindBin::Bin;
use RNAz;
use Getopt::Long;
use Pod::Usage;
use File::Spec qw/path catfile/;

# Set the path to your blast executable here
my $blastExecutable='blastall';

my $seqDir='';
my $blastDir='';
my $dbName='';
my $cutoff='1e-06';
my $version=0;
my $man=0;
my $help='';

GetOptions('database=s' => \$dbName,
		   'd=s' => \$dbName,
		   'seq-dir:s' => \$seqDir,
		   's:s' => \$seqDir,
		   'blast-dir=s' => \$blastDir,
		   'b=s' => \$blastDir,
		   '--e-value:f' => \$cutoff,
		   'e:f' => \$cutoff,
		   'version'=>\$version,
		   'man'=>\$man,
		   'v'=>\$version,
		   'help'=>\$help,
		   'h'=>\$help
		  ) or pod2usage(1);

pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

if ($version){
  print "\nrnazBlast.pl is part of RNAz $RNAz::rnazVersion\n\n";
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

if ($blastDir){
  if (not -d $blastDir){
	die("Blast directory '$blastDir' not a valid directory ($!).");
  }
} else {

  $blastDir=$ENV{BLASTDB};

  if (not -d $blastDir){
	die("No blast directory found. Set the BLASTDB variable or use the --blast-dir option.");
  }
}

my @PATH=File::Spec->path();

my $found=0;

foreach my $dir ("",@PATH){
  $found=1 if (-x File::Spec->catfile($dir, $blastExecutable));
}

if (!$found){
  print STDERR "Could not find blast executable. Ensure that 'blastn' is in your PATH or set the path directly in the perl program.\n";
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

	my @result=blastSeq($blastDir, $dbName, $cutoff, $seq,$blastExecutable);

	if (not @result){
	  print "$line\t-\n";
	} else {

	  my $name=$result[0]->{subjectID};
	  my $eValue=$result[0]->{e};
	  print "$line\t\"$name|$eValue\"\n";
	}
  }
}

__END__

=head1 NAME

C<rnazBlast.pl> - Compares predicted loci from data files as generated
by C<rnazCluster.pl> to a sequence database using BLAST.

=head1 SYNOPSIS

 rnazBlast.pl [options] [file]

=head1 OPTIONS

=over 8

=item B<-b> name, B<--blast-dir>=name

The directory with your BLAST database. If not set, the value from the
C<BLASTDB> environment variable is used.

=item B<-d> name, B<--database>=name

Name of the BLAST database to compare with. Must exist in the
directory set with C<--blast-dir> or in the directory set by
C<BLASTDB>.

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

=item B<-e> X, B<--e-value>=X

E-value cutoff. All hits with E < X are reported. (Default: 1e-06)

=item B<-v, --version>

Prints version information and exits.

=item B<-h --help>

Prints a brief help message and exits.

=item B<--man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

C<rnazBlast.pl> is a simple program to compare your hits to a sequence
database using BLAST. To use it you need (i) a sequence database (ii)
the sequence files to which the coordinates in your results file refer
(iii) a NCBI BLAST installation, i.e. a C<blastall> executable
somewhere.

First you have to create a BLAST index file for your sequence
database. You should have a FASTA formatted file of your
database. Assume for example that the file C<rfam> contains all
sequences of the Rfam database. Run the following command 

 # formatdb -t rfam -i rfam -p F

Make sure that you have the sequence files available and named
correctly (see notes for the C<--seq-dir> option). In this example we
assume that the files are in the subdirectory C<seq>

You can run the following command to compare each locus in the file
C<results.dat> with the newly created C<rfam> database (which is in
the subdirectory C<rfam>):

 # rnazBlast.pl --database=rfam --seq-dir=seq \
                --blast-dir=rfam --e-value=1e-06 \
                  results.dat > annotated.dat

If there is a hit better than E=1e-06 the name of the matching
sequence and the E-value is added in double quotes as additional field
to the locus line.

=head1 AUTHORS

Stefan Washietl <wash@tbi.univie.ac.at>

=cut

