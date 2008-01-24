#!/usr/bin/perl -w

# $Id: rnazIndex.pl,v 1.3 2008-01-24 10:26:45 wash Exp $

use strict;
use FindBin;
use lib $FindBin::Bin;
use RNAz;
use Getopt::Long;
use Pod::Usage;

my $gff=0;
my $bed=0;
my $fasta=0;
my $clusters=0;
my $hits=0;
my $help='';
my $ucsc='';
my $version=0;
my $man=0;
my $html='';
my $seqDir='';
my $forward=0;
my $reverse=0;

GetOptions('gff' => \$gff,
		   'g' => \$gff,
		   'bed' => \$bed,
		   'fasta'=>\$fasta,
		   'f'=>\$fasta,
		   'b' => \$bed,
		   'loci' => \$clusters,
		   'forward'=>\$forward,
		   'reverse'=>\$reverse,
		   'l' => \$clusters,
		   'windows' => \$hits,
		   'w' => \$hits,
		   'ucsc'=>\$ucsc,
		   'html'=>\$html,
		   'seq-dir:s' => \$seqDir,
		   'version'=>\$version,
		   'v'=>\$version,
		   'help'=>\$help,
		   'man'=>\$man,
		   'h'=>\$help
		  ) or pod2usage(1);

pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

if ($version){
  print "\nrnazIndex.pl is part of RNAz $RNAz::rnazVersion\n\n";
  print "http://www.tbi.univie.ac.at/~wash/RNAz\n\n";
  exit(0);
}


($clusters)=1 if (!$clusters and !$hits);
$forward=1 if (!$forward and !$reverse);
$bed=1 if (!$bed and !$gff and !$fasta);

my ($seqID, $name, $start, $end, $P, $strand);

$,="\t";

my $isCluster=0;
my $isHit=0;

if (($fasta) and (not -d $seqDir)){

  print STDERR "ERROR: Sequence directory '$seqDir' does not exist.\n";
  exit(1);
}


if ($html){

print
'<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
<html>
  <head>
    <title>Index</title>
    <meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">
    <STYLE type="text/css">
    TABLE {	border-spacing:1px;	background-color:#000000;}
    TD {font-size:10pt;padding:7px;font-family:helvetica,arial,impact,sans-serif;text-align:center}
    TH {font-size:10pt;background-color:#333399;color:#FFFFFF;padding:7px;font-family:helvetica,arial,impact,sans-serif;}
    .darkcell { background-color:#d4d9ff;}
    .lightcell { background-color:#e5e8ff;}
    P {font-size:11pt;font-family:helvetica,arial,impact,sans-serif;}
    H1 {font-size: 20pt;color:#333399;text-align: center;font-family: helvetica,arial,impact,sans-serif;padding: 10px;}
    H2 {color:#FFFFFF; background-color:#333399;width: auto; padding:3px; border-width:normal;font-size: 14pt; font-family:helvetica,arial,impact,sans-serif;}
    IMG {border-style:solid; border-width:thin; border-color:black;}
    PRE {font-size:10pt;font-family:courier,monospace;border-style:solid; border-width:thin;padding:5px;background-color:#e5e8ff;}
</STYLE>
  </head>
  <body>
<h1>Index</h1>
<table><tr><th>Locus ID</th><th>Seq. ID</th><th>Location</th><th>Custom annotation</th><th>Window ID</th><th>Location</th><th>Strand</th><th>N</th><th>Length</th><th>Mean pairwise ID</th><th>z</th><th>SCI</th><th>P</th></tr>
';

  my @clusters=();
  my %hits=();

  my $currCluster='';

  while (my $line=<>){
	if ($line=~/^\s?(locus\d+)/){
	  $currCluster=$1;
	  push @clusters, $line;
	}
	elsif ($line=~/^\s?window\d+/){
	  if (exists $hits{$currCluster}){
		push @{$hits{$currCluster}},$line;
	  } else {
		$hits{$currCluster}=[$line];
	  }
	}
  }

  my $cellColor='darkcell';

  foreach my $cluster (@clusters){
	#print $cluster;
	my @fields=split(/\t/,$cluster);
	($seqID, $start, $end, my $name, $P)=
	  ($fields[1],$fields[2],$fields[3],$fields[0],$fields[6]);
	
	my @customAnn=();

	my $i=8;
	
	while (defined $fields[$i]){
	  push @customAnn, $fields[$i];
	  $i++;
	}
	my $customAnn="&ndash;";
	
	if (@customAnn){
	  $customAnn=join("<br>",@customAnn);
	}
	
	my $clusterSize=@{$hits{$name}};

	my $clusterName=$name;

	$clusterName=~s/locus//;

	my $location=niceNumber($start)." &ndash; ".niceNumber($end);
	
	print "<tr><td rowspan=$clusterSize class=\"$cellColor\"><a href=\"locus$clusterName/index.html\">$clusterName</a></td>";
	print "<td rowspan=$clusterSize class=\"$cellColor\">$seqID</td>";
	print "<td rowspan=$clusterSize class=\"$cellColor\">$location</td>";
	print "<td rowspan=$clusterSize class=\"$cellColor\">$customAnn</td>";

	my $isFirst=1; # do not print <tr> if first entry
	foreach my $hit (@{$hits{$name}}){
	  my @fields=split(/\t/,$hit);
	  ($start, $end, $name, $strand, my $N, my $L, my $identity, my $z, my $SCI, my $P)=
		($fields[3],$fields[4],$fields[0],$fields[5],$fields[6],$fields[7],$fields[8],$fields[14],$fields[15],$fields[17]);


	  $name=~s/window//;

	  $location=niceNumber($start)." &ndash; ".niceNumber($end);
	
	  print "<tr>" if not $isFirst;
	  $isFirst=0;
	  print "<td class=\"$cellColor\"><a href=\"locus$clusterName/index.html\#window$name\">$name</a></td><td class=\"$cellColor\">$location</td><td class=\"$cellColor\">$strand</td>
             <td class=\"$cellColor\">$N</td><td class=\"$cellColor\">$L</td><td class=\"$cellColor\">$identity</td>
             <td class=\"$cellColor\">$z</td><td class=\"$cellColor\">$SCI</td><td class=\"$cellColor\">$P</td></tr>\n";
	}
	
	if ($cellColor eq 'darkcell'){
	  $cellColor='lightcell';
	} else {
	  $cellColor='darkcell';
	}
  }

my $timestamp=localtime;

print "</table><hr><p>Generated by <tt>rnazIndex.pl</tt> (part of <tt><a href=\"http://www.tbi.univie.ac.at/~wash/RNAz\">RNAz</a> $RNAz::rnazVersion</tt>) on $timestamp</p></body></html>";

}
else {


while (my $line=<>){

  next if $line=~/^\s?\#/;
  next if $line=~/^\s+$/;

  my @fields=split(/\t/,$line);
  chomp(@fields);

  if ($line=~/^\s?locus\d+/){
	$isCluster=1;$isHit=0;
	($seqID, $start, $end, $name, $P)=
	  ($fields[1],$fields[2],$fields[3],$fields[0],$fields[6]);
	
  } elsif ($line=~/^\s?window\d+/){
	$isHit=1;$isCluster=0;
	($seqID, $start, $end, $name, $P,$strand)=
	  ($fields[2],$fields[3],$fields[4],$fields[0],$fields[16],$fields[5]);
  }
	
  if ($ucsc){
	if ($seqID=~/^(.*)\.(.*)$/){
	  $seqID=$2;
	}
	#if ($P<0.5){
	#  $P=0;
	#} elsif (($P>=0.5) and ($P<0.9)){
	#  $P=500;
	#} elsif ($P>=0.9){
	#  $P=1000;
	#}
	$P=int($P*1000);
  }

  if ($bed){
	if ($isCluster and $clusters){
	  print $seqID,$start,$end,$name,$P,"\n";
	}
	if ($isHit and $hits){
	  print $seqID,$start,$end,$name,$P,$strand,"\n";
	}
  }

  if ($gff){
	if ($isCluster and $clusters){
	  print $seqID,"RNAz","structuredRNA",$start+1,$end,$P,".","id \"$name\"","\n";
	}
	if ($isHit and $hits){
	  print $seqID,"RNAz","structuredRNA",$start+1,$end,$P,$strand,"id \"$name\"","\n";
	}
  }

  if ($fasta){
	if ($isCluster and $clusters){

	  if ($forward) {
		$strand='+';
	  } else {
		$strand='-';
	  }

	  print ">$name ($seqID:$start..$end,$strand)\n";
	  my $seq=extractSeq($seqID,$start,$end);
	  $seq=~s/(.{60})/$1\n/g;
	  chomp($seq);
	  print $seq,"\n";
	}
	if ($isHit and $hits){
	  print ">$name ($seqID:$start..$end,$strand)\n";
	  my $seq=extractSeq($seqID,$start,$end,$strand);
	  $seq=~s/(.{60})/$1\n/g;
	  chomp($seq);
	  print $seq,"\n";
	}
  }
}
}


sub extractSeq{

  (my $chrom, my $start, my $end, my $strand)=@_;

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

	my @string=();
	
	foreach my $entry (@guess){
	  push @string, "\"$seqDir/$entry\"";
	}
	my $tmp=join(' or ',@string);
	print STDERR "No sequence data for '$chrom' found.\n";
	print STDERR "I could not find files named: $tmp";
	return ("");
  } else {
	return getSeq($realFile, $start, $end, $strand);
  }

}

__END__

=head1 NAME

C<rnazIndex.pl> - Convert data files as generated by C<rnazCluster.pl>
to different formats.

=head1 SYNOPSIS

 rnazIndex.pl [options] [file]

=head1 OPTIONS

=over 8

=item B<-g, --gff>

Generate GFF formatted output.

=item B<-b, --bed>

Generate BED formatted output.

=item B<-f, --fasta>

Get sequences in FASTA format for loci or windows. See options
C<--seq-dir>, C<--forward>, C<--reverse>!

=item B<--seq-dir>

Directory with sequence files. You only need this for FASTA output
(see option C<--fasta>). The files should be named with the sequence
identifier and the extension C<.fa> or C<.fasta>. If your identifier
in your input file is for example C<contig100> then you should have a
file named C<contig100.fa>. (If your identifier is of the form
``assembly.chromosome" as for example used by UCSC alignments, it is
also possible to name the file C<chr22.fa> for a sequence identifier
C<hg17.chr22>).

=item B<--forward, --reverse>

Only relevant for FASTA output (see option C<--fasta>). You can set if
you want the forward or reverse complement of the sequence
corresponding to a locus. Since loci don't have strand information you
might consider both strands for further analysis. Windows have strand
information, so if you export windows as FASTA these options are
ignored.

=item B<--ucsc>

In UCSC MAF alignment files it is common to use sequence identifiers
like for example ``hg17.chr22". However, in BED are usually specific
for a given assembly and therefore only ``chr22" is used in the BED
files. With this option you change any identifier of the form ``X.Y"
into ``Y". Moreover, the scores are multiplied by 1000 and rounded to
integers since the UCSC genome browser expects scores between 0 and
1000.

=item B<-l, --loci>

Use the locus information to generate the lines for the GFF and BED
files. This is the default.

=item B<-w, --windows>

Print the "windows" and not the "loci". Probably, rarely used
function.

=item B<--html>

With this option you get a HTML table which links to the the HTML
pages which you can create by using the C<--html> option in
C<rnazCluster.pl>. Redirect the output to some file which resides in
the C<results> directory created by C<rnazCluster.pl> and open the
file with your favourite web-browser.

=item B<-h, --help>

Prints a short help message and exits.

=item B<--man>

Prints a detailed manual page and exits.

=back

=head1 DESCRIPTION

C<rnazIndex.pl> reads tab-delimited data files as generated by
C<rnazCluster.pl> and converts them to GFF, BED or HTML formatted
files.

GFF is the most widely used annotation file format and supported by
many programs and systems
(http://www.sanger.ac.uk/Software/formats/GFF).

BED is the native annotation file format used by the UCSC genome
browser (http://genome.ucsc.edu).

=head1 EXAMPLES

 # rnazIndex.pl --gff results.dat > results.gff

Converts the C<results.dat> file to GFF format.

 # rnazIndex.pl --ucsc --bed results.dat > results.bed

Create UCSC style BED format.

 # rnazIndex.pl --html results.dat > results/index.html

Generates HTML formatted table.

 # rnazIndex.pl --forward --fasta --seq-dir=seq results.dat

Exports sequences in FASTA format.

=head1 AUTHOR

Stefan Washietl <wash@tbi.univie.ac.at>

=cut
