#!/usr/bin/perl -w

use File::Basename;

# creates the man.tex file which is included in tutorial.tex

my $outFile='man.tex';

my @pods=("../man/RNAz.pod",
		  "../perl/rnazWindow.pl",
		  "../perl/rnazCluster.pl",
		  "../perl/rnazSelectSeqs.pl",
		  "../perl/rnazFilter.pl",
		  "../perl/rnazSort.pl",
		  "../perl/rnazAnnotate.pl",
		  "../perl/rnazBlast.pl",
		  "../perl/rnazIndex.pl",
		  "../perl/rnazBEDsort.pl",
		  "../perl/rnazBEDstats.pl",
		  "../perl/rnazMAF2BED.pl",
		  "../perl/rnazRandomizeAln.pl");

open(OUT,">$outFile");

foreach my $pod (@pods){

  (my $baseName,my $dir,my $ext) = fileparse($pod, qr/\..*/);

  system("pod2latex -modify -h1level 2 $pod");
  system("pod2html --noindex --title=$baseName --infile $pod --outfile html/$baseName.html");
  open(IN,"<$baseName.tex") || die("$baseName.tex not found ($!)");

  my $in='';

  while (<IN>){
	$in.=$_;
  }

  print OUT $in;

  print OUT '\newpage';

  unlink "$baseName.tex";


}
