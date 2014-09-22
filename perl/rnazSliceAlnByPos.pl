#!/usr/bin/perl -w

use strict;
use FindBin;
use lib $FindBin::Bin;
use RNAz;
use AnnotUtil;
use Getopt::Long;
use Pod::Usage;



# input
my $is_file;       # intersect file
my $pre_mafin;     # name maf in

# output
my $prefix;        # name outfile

# parameters
my $flank = 100;            # flanking region
my $alnFormat  = "MAF";     # input aln format
my $out_format = "CLUSTAL"; # output aln format
my $file_size  = 100;       # flanking regions [nt]

GetOptions(
  'is:s'	   => \$is_file,
  'mafin|maf:s'	   => \$pre_mafin,
  'out|o:s'	   => \$prefix,
  'flank|f:i'	   => \$flank,
  'outformat|of:s' => \$out_format
    );

# Known output formats
my %out_formats = ("CLUSTAL" => undef,
                   "MAF"     => undef);

# Error message: no output format or not exsisting format specified
#die("Please specify a valid output format:\n", join(", ", (sort keys %out_formats)),"\n") unless (defined(uc($out_format)) || exists($out_formats{uc($out_format)}));

# set the format
#$out_formats{$out} = 1;

my %files = ();
my %gqs   = ();

open(IS, "<$is_file") || die "could not open intersect file $is_file $!\n";

while(<IS>){
  next if /^\#/ || /^\s/;
  chomp;
  my ($fields1,$attribs1,$fields2,$attribs2) = splitIntersectBedGtf($_);

  my $gqID    = $attribs1->{transcript_id}->[0];
  my $blockNr = $attribs2->{blockNr}->[0];
  s/["]//g foreach ($blockNr, $gqID);

  # calculate the file in which the maf block resides
  # $file_size is the number of maf blocks per file
  # calculation is based on convertion of alignment number
  # within a file to block nr within a chromosome
  # see processMafGtf.pl
  # my $shift = $file_size * ($off-1);
  # off is file name
  # example: file 1   blockNr 1..100
  #               2           101-200
  #
  # blockNr   = 100
  # file1     = int((100+ (100-1))/100 ) = 1
  # file pos1 = 100 - 100*(1-1) = 100
  #
  # blockNr   = 200
  # file2     = int((200+ (100-1))/100 ) = 2
  # file pos2 = 200 - 100*(2-1) = 100
  
  my $fileNr  = int(($blockNr+($file_size-1))/$file_size);
  my $filepos = $blockNr - $file_size*($fileNr-1);

  # store file nr and file pos (= the nth alignment block in maf file)
  push(@{$files{$fileNr}{$filepos}}, $gqID);

  # store intersect data for this gquad
  @{$gqs{$gqID}} = ($fields1,$attribs1,$fields2,$attribs2);
}
close(IS);


# --------------------------------------------------
# Process maf files
# --------------------------------------------------

# open GTF output
my $sliceGQlocal = join(".",$prefix,"gquad_local","gtf");
open(GQGTF, ">$sliceGQlocal") || die "could not open gtf out file $sliceGQlocal $!\n";

my $sliceout = join(".",$prefix,"slice","gtf");
open(SGTF, ">$sliceout") || die "could not open gtf out file $sliceout $!\n";


my $sliceID = 0;

foreach my $fileNr (keys %files ){
  
  my @filepos = (sort {$a <=> $b} keys %{$files{$fileNr}});
  
  # --------------------------------------------------
  # Open the determined maf file
  # --------------------------------------------------
  
  my $fileName = join(".",$pre_mafin,$fileNr,"maf");
  my $fh;
  open( $fh, "<$fileName" ) || die("Could not open file $fileName ($!)");
  
  # read all alignment blocks in curent maf file
  # IMPORTANT!!!!!!! use same parames (if any) used in other studies
  # to filter aln blocks!!!
  my $alnCounter = 0;
  
  $alnFormat = checkFormat($fh);

  while ( my $alnString = getNextAln( $alnFormat, $fh ) ) {
    
    # Process curent alignment block
    # --------------------------------------------------
    
    my ($fullAln, $orgs) = parseAln( $alnString, $alnFormat, "1");
        
    # only reference species in alignment
    next if scalar(@$fullAln) < 2;

    $alnCounter++;      # no. of alignment block in file
    
    # Check if current aln block needs processing
    next unless ($filepos[0] == $alnCounter);
    
    # Make a slice for each gquad that lies within the current aln block
    foreach my $gqID (@{$files{$fileNr}{$filepos[0]}}){
      $sliceID++;

      # make a copy of $fullAln
      my @tmpAln = @$fullAln;
      
      # Description of slices
      # --------------------------------------------------
      
      # get this from intersection file
      my $sliceStart  = ($gqs{$gqID}[0]->{start}) - $flank;
      my $sliceEnd    = $gqs{$gqID}[0]->{end}   + $flank;
      
      my $alnStart    = $gqs{$gqID}[2]->{start};
      my $alnEnd      = $gqs{$gqID}[2]->{end};
      
      # check if either side of slice is within range of aln block
      # if not, adjust to aln coords
      $sliceStart = $alnStart if ($alnStart > $sliceStart);
      $sliceEnd   = $alnEnd   if ($alnEnd   < $sliceEnd);
      
      
      my $slice = sliceAlnByPos(\@tmpAln,"0",($sliceStart-1),$sliceEnd);
      
      
      # remove emty lines and common gaps in slice
      my $gqLength   = $gqs{$gqID}[0]->{end}   - $gqs{$gqID}[0]->{start} + 1;
      $slice = &removeEmptyAlnLines($slice, $gqLength);
#    removeCommonGaps( $slice );
      
      # slice has only RefSeq or only one seq
      $sliceID--, next if ( scalar(@$slice) < 2 );
      
      # mean pairwise ID
      my $meanPairID = meanPairID($slice) * 100;
      
      # make reverse complement of the alignment if gquad was found on minus strand
      $slice = revAln($slice) if $gqs{$gqID}[0]->{strand} eq "-";
      
      # Print slice
      my $sliceOut = join(".", $prefix, $sliceID, "aln");
      open(SLICE, ">$sliceOut") || die "could not open output file $sliceOut $!\n";
      print SLICE (formatAln( $slice, $out_format ));
      close(SLICE);
      
      
      # Print annotation of gquad - local gap free coordinates
      # = local coordinates of gquad within reference sequence
      # in bed format = zero-based
      
      
      my $localStart = $gqs{$gqID}[0]->{start} - $sliceStart   +1;
      my $localEnd   = $localStart + $gqLength -1;
      
      my %gqattribs = ("gene_id"       => ["\"".$sliceID."\""],
		       "transcript_id" => ["\"".$gqID."\""],
		       "blockNr"       => [$gqs{$gqID}[3]->{blockNr}->[0]]
	  );
      
      my %gqfields = ("chr"        => $sliceID,
		      "source"     => "gquad",
		      "type"       => "exon",
		      "start"      => $localStart,
		      "end"        => $localEnd,
		      "score"      => ".",
		      "strand"     => $gqs{$gqID}[0]->{strand},
		      "phase"      => ".",
		      "attributes" => ""
	  );
      printGtfFieldsAndAttributes( \%gqfields, \%gqattribs, *GQGTF);
      
      
      # Print slice annotation
      my $sliceLength = $sliceEnd - $sliceStart +1;
    
      my @orgs = ();
      foreach my $k (sort keys %$orgs){
	push(@orgs, "\"".$k."\"");
      }
      
      my %sliceattribs = ("gene_id",       => ["\"".$sliceID."\""],
			  "transcript_id"  => ["\"".$sliceID."\""],
			  "meanPairID"	 => ["\"".$meanPairID."\""],
			  "slice_length"	 => ["\"".$sliceLength."\""],
			  "org_block"	 => \@orgs,
			  "nr_org"	 => ["\"".(scalar(keys %$orgs))."\""],
			  "blockNr"        => [$gqs{$gqID}[3]->{blockNr}->[0]]
	  );
      
      my %slicefields = ("chr"        => $gqs{$gqID}[0]->{chr},
			 "source"     => "alnslice",
			 "type"       => "exon",
			 "start"      => $sliceStart,
			 "end"        => $sliceEnd,
			 "score"      => ".",
			 "strand"     => $gqs{$gqID}[0]->{strand},
			 "phase"      => ".",
			 "attributes" => ""
	  );
      printGtfFieldsAndAttributes( \%slicefields, \%sliceattribs, *SGTF);
      
    } # done with all gqIDs in current maf block
    
    # done with this aln block -> remove it from list
    shift @filepos;
    last unless defined($filepos[0]); # list if relevant maf blocks is empty. done!
  } # I have red all blocks within current file
  close($fh);
}

close(GQGTF);
close(SGTF);



# --------------------------------------------------
# Subfunctions
# --------------------------------------------------


sub removeEmptyAlnLines{
  my ($slice,$gqLength) = @_;

  my $sliceLength = length($slice->[0]->{seq});

  my $maxgaps = $sliceLength - $gqLength;
  
  for (my $i=0; $i<=$#$slice; $i++) {
    my $numGaps = ( $slice->[$i]->{seq} =~ tr/-./-/ );
    splice(@$slice,$i,1), $i-- if ( $numGaps > $maxgaps);
  }
  return($slice);
}
