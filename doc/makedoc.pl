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


if (!-e 'html'){
  mkdir "html";
}

open(OUT,">$outFile");

print "Converting PODs...\n";

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


# add link to stylesheet in each html file

my $css='<link rel="stylesheet" href="../style.css" type="text/css"> ';

foreach my $file (glob("html/*.html")){

  open(IN,"<$file");

  my @lines=<IN>;

  close(IN);

  open(OUT,">$file");

  foreach my $line (@lines){
	$line=~s/(<\/head>)/$css \n$1/;
	print OUT $line;
  }
}


print "TeXing manual...\n";

`latex manual`;
`bibtex manual`;
`latex manual`;
`latex manual`;
`dvipdf manual.dvi`;
`cp manual.pdf ..`;


print "Creating manpage...\n";
chdir "../man";
`pod2man RNAz.pod > RNAz.1`
