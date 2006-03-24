#!/usr/bin/perl -w

# $Id: rnazCluster.pl,v 1.2 2006-03-24 15:43:13 wash Exp $

use strict;
use FindBin;
use lib $FindBin::Bin;
use RNAz;
use Getopt::Long;
use Pod::Usage;
use File::Temp qw/tempfile/;
use File::Spec qw/path catfile/;

my $cutoff=0.5;
my $printHits=0;
my $printClusters=0;
my $printHeader=0;
my $html=0;
my $version=0;
my $man=0;
my $help=0;

# Set custom locations of the following programs for the --html output:

my $alifoldProg='RNAalifold';
my $gsProg='gs';
my $colorrnaProg='colorrna.pl';
my $coloralnProg='coloraln.pl';




GetOptions('cutoff:f' => \$cutoff,
		   'c:f' => \$cutoff,
		   'windows'=>\$printHits,
		   'loci'=>\$printClusters,
		   'w'=>\$printHits,
		   'l'=>\$printClusters,
		   'header'=>\$printHeader,
		   'd'=>\$printHeader,
		   'html'=>\$html,
		   'version'=>\$version,
		   'man'=>\$man,
		   'v'=>\$version,
		   'help'=>\$help,
		   'h'=>\$help
		  ) or pod2usage(1);

pod2usage(1) if $help;
pod2usage(-verbose => 2) if $man;

if ($version){
  print "\nrnazCluster.pl is part of RNAz $RNAz::rnazVersion\n\n";
  print "http://www.tbi.univie.ac.at/~wash/RNAz\n\n";
  exit(0);
}


if (!$printHits and !$printClusters){
  $printHits=$printClusters=1;
}

if ($printHeader){

  if ($printClusters){
	print "# locusID\tseqID\tstart\tend\tmaxN\tmaxIdentity\tmaxP\tminZ\n";
  }

  if ($printHits){
	print "# windowID\tclusterID\tseqID\tstart\tend\tstrand\tN\tcolumns\tidentity\tmeanMFE\tconsensusMFE\tenergyTerm\tcovarianceTerm\tcombPerPair\tz\tSCI\tdecValue\tP\n";
  }
}

if ($html){
  mkdir "results" || die("Could not create directory for HTML files ($!)");

  my $alifoldFlag=0;
  my $colorrnaFlag=0;
  my $coloralnFlag=0;
  my $gsFlag=0;

  my @PATH=File::Spec->path();

  foreach my $dir ("",@PATH){
	
	$alifoldFlag=1 if (-x File::Spec->catfile($dir, $alifoldProg));
	$colorrnaFlag=1 if (-x File::Spec->catfile($dir, $colorrnaProg));
	$coloralnFlag=1 if (-x File::Spec->catfile($dir, $coloralnProg));
	$gsFlag=1 if (-x File::Spec->catfile($dir, $gsProg));
	
	last if $alifoldFlag && $colorrnaFlag && $coloralnFlag && $gsFlag;
  }

  if (!($alifoldFlag && $colorrnaFlag && $coloralnFlag && $gsFlag)){

	print STDERR <<END;

To use the --html option you need additional programs. Make sure that
the following executables are within your PATH of executables or
directly set the location of the programs in the source code of the
rnazCluster.pl script:
END

	print STDERR "'RNAalifold' not found.\n" if not $alifoldFlag; 
	print STDERR "'colorrna.pl' not found.\n" if not $colorrnaFlag;
	print STDERR "'coloraln.pl' not found.\n" if not $coloralnFlag;
	print STDERR "'gs' not found.\n" if not $gsFlag;

	exit(1);
  }
}


my ($currName,$currStart,$currEnd,$currStrand);
my ($prevName,$prevStart,$prevEnd);

$prevName=''; $prevStart=0; $prevEnd=0;

my ($maxP, $maxN, $maxID, $maxZ);
$maxN=0;$maxID=0;$maxZ=100;$maxP=0;
my $minStart=99000000;
my $maxEnd=0;

my @hits=();

my $clusterID=1;
my $hitID=1;

my $fileName=shift @ARGV;
my $fh;

if (!defined $fileName){
  $fh=*STDIN;
} else {
  open($fh,"<$fileName") || die("Could not open file $fileName ($!)");
}

while (my $rnazString=getNextRNAz($fh)){

  my $results=parseRNAz($rnazString);

  $currStart=$results->{refSeqStart};
  $currEnd=$results->{refSeqEnd};
  $currName=$results->{refSeqName};
  $currStrand=$results->{refSeqStrand};

  if ($results->{P}>$cutoff){
	
	# if there is NO overlap start new cluster
	if (!(($currName eq $prevName) and ($currStart <= $prevEnd))){

	  if ($maxEnd!=0){ # omit first one

		if ($printClusters){
		  print "locus$clusterID\t$prevName\t$minStart\t$maxEnd\t$maxN\t$maxID\t$maxP\t$maxZ\n";
		}

		if ($printHits){
		  foreach my $hit (@hits){
			$hit->{clusterID}="locus$clusterID";
			hitLine($hit);
		  }
		}

		if ($html){
		  createHTML(\@hits);
		}
		
		@hits=();
		$maxN=0;
		$maxID=0;
		$maxZ=100;
		$maxP=0;
		$minStart=999000000;
		$maxEnd=0;
		
		$clusterID++;
	  }
	}

	$maxN=$results->{N} if ($results->{N}>$maxN);
	$maxP=$results->{P} if ($results->{P}>$maxP);
	$maxID=$results->{identity} if ($results->{identity}>$maxID);
	$maxZ=$results->{z} if ($results->{z}<$maxZ);
	$minStart=$currStart if ($currStart<$minStart);
	$maxEnd=$currEnd if ($currEnd>$maxEnd);

	#my %currHit=(hitID=>$hitID,name=>$currName,start=>$currStart,end=>$currEnd,strand=>$currStrand);
	
	#my $tmp="window$hitID\tCLUSTER\t$currName\t$currStart\t$currEnd\t$currStrand";
	#foreach my $key (qw(N columns identity meanMFE consensusMFE energy covariance combPerPair z sci decValue P)){
	#  $tmp.="\t$results->{$key}";
	#  $currHit{$key}=$results->{$key};
	#}
	#push @hits,"$tmp\n";

	push @hits,{%$results,hitID=>"window$hitID"};

	$hitID++;
	
	($prevName,$prevStart,$prevEnd)=($currName,$currStart,$currEnd);

	#exit if $hitID>10;

  }
}

# don't forget last cluster...

if (@hits){
  if ($printClusters){
	print "locus$clusterID\t$prevName\t$minStart\t$maxEnd\t$maxN\t$maxID\t$maxP\t$maxZ\n";
  }

  if ($printHits){
	foreach my $hit (@hits){
	  $hit->{clusterID}="locus$clusterID";
	  hitLine($hit);
	}
  }
}

sub toPNG{

  my $file=shift;

  (my $baseName)=($file=~/(.*)\.ps/);

  system("gs -r72 -q -dNOPAUSE -dQUIET -dBATCH -dEPSCrop -dTextAlphaBits=4 -dGraphicsAlphaBits=2  -sDEVICE=png16m -sOutputFile=$baseName.png $file");

}



sub hitLine{

  my %data=%{$_[0]};
  my @tmp=();
  foreach my $key (qw(hitID clusterID refSeqName refSeqStart refSeqEnd refSeqStrand N columns identity meanMFE consensusMFE energy covariance combPerPair z sci decValue P)){
	push @tmp,$data{$key};
  }
  print join("\t",@tmp);
  print "\n";
}



sub createHTML{

  my @hits=@{$_[0]};

  my $locusTemplate='
<!DOCTYPE HTML PUBLIC "-//W3C//DTD HTML 4.01 Transitional//EN">
<html>
  <head>
    <title>Locus %clusterID%</title>
    <meta http-equiv="Content-Type" content="text/html; charset=iso-8859-1">
    <STYLE type="text/css">
    TABLE {	border-spacing:1px;	background-color:#000000;}
    TD {font-size:10pt;padding:7px;font-family:helvetica,arial,impact,sans-serif; background-color:#e5e8ff;font-weight: bold;}
    TH {font-size:10pt;background-color:#d4d9ff;padding:7px;font-family:helvetica,arial,impact,sans-serif; text-align:right;  }
    P {font-size:11pt;font-family:helvetica,arial,impact,sans-serif;}
    H1 {font-size: 20pt;color:#333399;text-align: center;font-family: helvetica,arial,impact,sans-serif;padding: 10px;}
    H2 {color:#FFFFFF; background-color:#333399;width: auto; padding:3px; border-width:normal;font-size: 14pt; font-family:helvetica,arial,impact,sans-serif;}
    IMG {border-style:solid; border-width:thin; border-color:black;}
    PRE {font-size:10pt;font-family:courier,monospace;border-style:solid; border-width:thin;padding:5px;background-color:#e5e8ff;}


</STYLE>
  </head>
  <body>
    <h1>Locus %clusterID%</h1>
	
	<table>
		<tr>
		  <th>Sequence ID</th>
		  <td>%seqID%</td>
		</tr>
		<tr>
		  <th>Location</th>
		  <td>%minStart% &ndash; %maxEnd%</td>
		</tr>
		<tr>
		  <th>Length</th>
		  <td>%length%</td>
		</tr>
		<tr>
		  <th>Max. P</th>
		  <td>%maxP%</td>
		</tr>
	</table>
<map name="map">
%imageMap%
</map>
<p><img src="map.png" usemap="#map" alt="overview"></p>
	%hits%
<hr>
<p>Generated by <tt>rnazCluster.pl</tt> (part of <tt><a href="http://www.tbi.univie.ac.at/~wash/RNAz">RNAz</a> $RNAz::rnazVersion</tt>) on %timestamp%</p>
  </body>
</html>';

  my $windowTemplate='
<h2><a name="%hitID%">Window %hitIDnumber%</a></h2>
	<table>
		<tr>
		  <th>Location</th>
		  <td>%start% &ndash; %end%</td>
		</tr>
		<tr>
		  <th>Length</th>
		  <td>%length%</td>
		</tr>
%rnazUpper%
	</table>
<p>Download alignment: <a href="%hitID%.aln">ClustalW</a> | <a href="%hitID%.maf">MAF</a></p>
<pre>%rnazLower%</pre>
<p>
<img src="%hitID%_aln.png" alt="alignment"><br>
</p>
<p><a href="%hitID%_aln.ps">Postscript</a></p>

<p><img src="%hitID%_alirna.png" alt="secondary structure"></p>
<p><a href="%hitID%_alirna.ps">Postscript</a></p>
<p><img src="%hitID%_alidot.png" alt="dotplot"></p>
<p><a href="%hitID%_alidot.ps">Postscript</a></p>

';

  chdir "results";
  mkdir "locus$clusterID";
  chdir "locus$clusterID";

  open(INDEX,">index.html") or die("Could not open index.html for writing ($!)");

  my $tmpLocusTemplate=$locusTemplate;

  my $tmpStart=niceNumber($minStart);
  my $tmpEnd=niceNumber($maxEnd);

  $tmpLocusTemplate=~s/%clusterID%/$clusterID/gs;
  $tmpLocusTemplate=~s/%seqID%/$prevName/gs;
  $tmpLocusTemplate=~s/%minStart%/$tmpStart/gs;
  $tmpLocusTemplate=~s/%maxEnd%/$tmpEnd/gs;
  $tmpLocusTemplate=~s/%maxP%/$maxP/gs;

  my $length=$maxEnd-$minStart;

  $tmpLocusTemplate=~s/%length%/$length/gs;

  my $tmpWindowTemplate=$windowTemplate;
  my $hitsOutput='';

  my $cellClass="darkcell";

  my $imageMap=createMap(\@hits);

  foreach my $hit (@hits){

	my $id=$hit->{hitID};

	(my $idNumber)=($id=~/.*(\d+)/);

	$tmpStart=niceNumber($hit->{refSeqStart});
	$tmpEnd=niceNumber($hit->{refSeqEnd});
	
	$tmpWindowTemplate=~s/%hitID%/$id/gs;
	$tmpWindowTemplate=~s/%hitIDnumber%/$idNumber/gs;
	$tmpWindowTemplate=~s/%start%/$tmpStart/gs;
	$tmpWindowTemplate=~s/%end%/$tmpEnd/gs;
	my $length=$hit->{refSeqEnd}-$hit->{refSeqStart};
	$tmpWindowTemplate=~s/%length%/$length/gs;

	(my $rnazUpper, my $rnazLower)=($hit->{rawOutput}=~/(.*\#$)(.*)$/ms);

	$rnazUpper=~s/#//g;
	$rnazUpper=~s/^(.*?):(.*?)$/<tr><th>$1<\/th><td>$2<\/td><\/tr>/msg;
			
	$tmpWindowTemplate=~s/%rnazUpper%/$rnazUpper/gs;
	$tmpWindowTemplate=~s/%rnazLower%/$rnazLower/gs;

	$hitsOutput.=$tmpWindowTemplate;
	$tmpWindowTemplate=$windowTemplate;
	open(ALN,">$id.maf");

	print ALN formatAln($hit->{aln},"MAF");

	close(ALN);
			
	open(ALN,">$id.aln");

	foreach my $seq (@{$hit->{aln}}){
	  $seq->{start}=undef;
	  $seq->{end}=undef;
	}
			
	print ALN formatAln($hit->{aln},"CLUSTAL");


	my $devnull=File::Spec->devnull();
	system("RNAalifold -p <$id.aln >$devnull 2>$devnull");
	system("coloraln.pl < $id.aln >$id\_aln.ps");
	system("colorrna.pl alirna.ps alidot.ps >$id\_alirna.ps");
	unlink("alirna.ps");
	rename("alidot.ps","$id\_alidot.ps");

	toPNG("$id\_alirna.ps");
	toPNG("$id\_alidot.ps");
	toPNG("$id\_aln.ps");
					
	close ALN;

  }

  $tmpLocusTemplate=~s/%imageMap%/$imageMap/gs;

  $tmpLocusTemplate=~s/%hits%/$hitsOutput/gs;

  my $now=localtime;
			
  $tmpLocusTemplate=~s/%timestamp%/$now/gs;

  print INDEX $tmpLocusTemplate;

  close INDEX;
  chdir "../../";
}

sub createMap{


  my @hits=@{$_[0]};

  my $psTemplate='%!PS-Adobe-3.0 EPSF-3.0
%%BoundingBox: 0 0 %width% %height%
%%EndComments

/string { % (Text) x y red green blue
  gsave
	setrgbcolor
	18 add
	moveto
	show
  grestore
} def


0 %height% translate
1 -1 scale
/Helvetica findfont
[12 0 0 -12 0 0] makefont setfont

%/Font /Helvetica findfont def
%/FontSize 12 def

  /rightarrow { % X Y length width red green blue
  gsave
	setrgbcolor
	newpath

	% X,Y
	3 index
	3 index
	moveto

	% X+length, Y
	3 index
	2 index
	add
	10 sub
	3 index
	lineto

	% Draw arrow
	3 index
	2 index
	add
	1 index
	2 div
	4 index
	add
	lineto

	% X+length, Y+width
	3 index
	2 index
	add
	10 sub
	3 index
	2 index
	add
	lineto

	% X, Y+width
	3 index
	3 index
	2 index
	add
	lineto
	
	clear
	closepath
	fill
  grestore
} def

/leftarrow { % X Y length width red green blue
  gsave
	setrgbcolor
	newpath

	% X,Y
	3 index
	10 add
	3 index
	moveto

	% X+length, Y
	3 index
	2 index
	add
	3 index
	lineto

	% X+length, Y+width
	3 index
	2 index
	add
	3 index
	2 index
	add
	lineto

	% X, Y+width
	3 index
	10 add
	3 index
	2 index
	add
	lineto
	
	% Draw arrow, 1/20 of head
	3 index
	1 index
	2 div
	4 index
	add
	lineto

	clear
	closepath
	fill
  grestore
} def
  


0 0 0 setrgbcolor

%100 100 200 7 0 1 0 rightarrow
%(This is a string) 100 100 1 0 0 string

  %data%

showpage';


  open(MAP,">map.ps") or die("Could not open map.ps for writing ($!)");

my $data='';
my $imageMap='';

# $maxEnd and $minStart are global variables!

my $realClusterLength=$maxEnd-$minStart;

my $clusterLength=400;

my $factor=$clusterLength/$realClusterLength;

my $lineStep=40;
my $padding=20;
my $arrowWidth=6;
my $rulerWidth=3;
my $color90='0 0 1';
my $color50='0.8 0.8 1';
my $fontColor='0 0 0';
my $rulerColor='0 0 0';
my $currY=$padding;


foreach my $hit (@hits){

my $realLength=$hit->{refSeqEnd}-$hit->{refSeqStart};
my $length=$realLength*$factor;
my $start=($hit->{refSeqStart}-$minStart)*$factor;

my $color=$color50;
my $command='rightarrow';

$command='leftarrow' if ($hit->{refSeqStrand} eq '-');

$color=$color90 if ($hit->{P}>=0.9);

$start+=$padding;

$data.="$start $currY $length $arrowWidth $color $command\n";

my $tmpX1=int($start);
my $tmpY1=int($currY);
my $tmpX2=int($start+$length);
my $tmpY2=int($currY+$arrowWidth);

$imageMap.="<area shape=\"rect\" coords=\"$tmpX1,$tmpY1,$tmpX2,$tmpY2\" href=\"\#$hit->{hitID}\" alt=\"$hit->{hitID}\">\n";


my $id=$hit->{hitID};
$id=~s/window//;

my $P=sprintf("%.2f",$hit->{P});

my $ann="Window $id \($P\)";

$data.="($ann) $start $currY $fontColor string\n";


$currY+=$lineStep;

}


$data.="$padding $currY  $clusterLength $rulerWidth $rulerColor rightarrow";

my $ann=$hits[0]->{refSeqName}.": $minStart - $maxEnd";

$data.="($ann) $padding $currY $fontColor string\n";

my $width=$clusterLength+2*$padding;

$currY+=$lineStep+$padding;

$psTemplate=~s/%height%/$currY/gs;
$psTemplate=~s/%width%/$width/gs;
$psTemplate=~s/%data%/$data/gs;

print MAP $psTemplate;

toPNG("map.ps");

return $imageMap;

}


__END__

=head1 NAME

C<rnazCluster.pl> - Cluster RNAz hits and print a summary of the results.

=head1 SYNOPSIS

 rnazCluster.pl [options] [file]

=head1 OPTIONS

=over 8

=item B<-c> X, B<--cutoff>=X

Only consider hits with RNAz class probablility P>X (Default:B<0.5>)

=item B<-w>, B<--windows>

=item B<-l>, B<--loci>

Set these flags to print information for ``windows" and/or ``loci" in
the output. By default, both single windows and combined loci are
printed.

=item B<-d>, B<--header>

Print a header explaining the fields of the output (see below for a
detailed description of the fields).

=item B<--html>

Generates HTML formatted output of the results in the subdirectory
C<results>. For this option to work you need to have installed
ghostscript and a few programs from the ViennaRNA package. More
precisely you need the following executables in your PATH: C<gs>,
C<RNAalifold>, C<colorrna.pl>, C<coloraln.pl>. Alternatively you can
adjust the locations of these programs directly in the
C<rnazCluster.pl> script. Please note that if you use this option the
program will get B<very slow> because the figures have to be
generated. It is also important that you have run RNAz with the
B<C<--show-gaps>> option!

=item B<-v, --version>

Prints version information and exits.

=item B<-h, --help>

Prints a short help message and exits.

=item B<--man>

Prints a detailed manual page and exits.

=back

=head1 DESCRIPTION

C<rnazCluster.pl> reads RNAz output files and combines hits in
overlapping windows to ``loci". It prints a summary of the windows
and/or loci as a tabulator delimited text to the standard output. An
explanation of the fields can be found below. See the user manual for
a more detailed meaning of these values.

To work properly, your RNAz output file needs to contain position
information. This means there must have been genomic locations in your
original alignments you scored with RNAz (i.e. MAF files with a
reference sequence). Moreover, the original input alignments have to
be B<ordered by the genomic location of the reference sequence>.

If you want HTML output please see the notes for the C<--html> option
above.

=head1 FIELDS

B<"Window" lines>

=over 8

=item 1. B<windowID>

Consecutive numbered ID for each window

=item 2. B<locusID>

The locus which this window belongs to

=item 3. B<sequenceID>

Identifier of the sequence (e.g. human.chr1 or contig42)

=item 4. B<start>

Start position of the reference sequence in the window

=item 5. B<end>

End position of the reference sequence in the window

=item 6. B<strand>

Indicates if the reference sequence is from the positive or
negative strand

=item 7. B<N>

Number of sequences in the alignment

=item 8. B<columns>

Number of columns in the alignment

=item 9. B<identity>

Mean pairwise identity of the alignment

=item 10. B<meanMFE>

Mean minimum free energy of the single sequences as
calculated by the RNAfold algorithm

=item 11. B<consensusMFE>

``consensus MFE" for the alignment as calculated by
the RNAalifold algorithm

=item 12. B<energyTerm>

Contribution to the consensus MFE which comes from the energy part of
the RNAalifold algorithm

=item 13. B<covarianceTerm>

Contribution to the consensus MFE which comes from the covariance part
of the RNAalifold algorithm

=item 14. B<combPerPair>

Number of different base combinations per predicted
pair in the consensus seconary structure

=item 15. B<z>

Mean z-score of the sequences in the alignment

=item 16. B<SCI>

Structure conservation index for the alignment

=item 17. B<decValue>

Support vector machine decision value

=item 18. B<P>

RNA class probability as calculated by the SVM

=back

B<"Loci" lines>

=over 8

=item 1. B<locusID>

Consecutive numbered ID for each locus

=item 3. B<sequenceID>

Identifier of the sequence (e.g. human.chr1 or contig42)

=item 4. B<start>

Start position of the reference sequence in the window

=item 5. B<end>

End position of the reference sequence in the window

=item 6. B<strand>

Indicates if the reference sequence is from the positive or
negative strand

=item 7. B<maxN>

Maximum number of sequences in the alignments of this locus

=item 9. B<maxIdentity>

Maximum mean pairwise indentity in the alignments of this locus

=item 9. B<maxP>

Maximum RNA class probability in the alignments of this locus

=item 9. B<minZ>

Minimum z-score in the alignments of this locus.


=back


=head1 EXAMPLES

 # rnazCluster.pl rnaz.out

Parses and clusters the hits in the file C<rnaz.out> and prints loci
and cluster information to the standard output.

 # rnazCluster.pl -c 0.9 --html rnaz.out > results90.out

Clusters all hits from the file C<rnaz.out> with P>0.9, writes the
tab-delimited output to the file C<results90.out> and, at the same
time, generates a website in a subdirectory called C<results>.

=head1 AUTHORS

Stefan Washietl <wash@tbi.univie.ac.at>

=cut
