#!/usr/bin/perl

use strict;
use warnings;

my $model = $ARGV[0];


my @lines = `cat $model`;
open (FILE, ">$model.inc");
print FILE "char* $model\_string=\n";
my $k = 0;
foreach my $line (@lines)
{
    chomp $line;
    if ($k == $#lines)
    {
	    print FILE '"',$line,'\n";',"\n";
    }
    else
    {
	print FILE '"',$line,'\n"',"\n";
    }
    $k++;
}
close FILE;

