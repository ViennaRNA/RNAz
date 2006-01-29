#!/usr/local/bin/perl -w


#missing: re-introduce track line, more intelligent sort of seq-ids (chr9 before chr10)


use strict;

my %lines=();

while (my $line=<>){

  chomp($line);

  next if ($line=~/track/);
  next if ($line=~/^$/);

  my @fields=split(/\t/,$line);

  my $seqID=$fields[0];
  my $start=$fields[1];
  my $end=$fields[2];

  if (!exists($lines{$seqID})){
	$lines{$seqID}=[[@fields]];
  } else {
	push @{$lines{$seqID}},[@fields];
  }
}


foreach my $key (sort keys %lines){

  my @tmp=sort {$a->[1] <=> $b->[1]} @{$lines{$key}};

  foreach my $line (@tmp){
	#print $line->[1];
	print join("\t",@$line),"\n";
  }

}
	
