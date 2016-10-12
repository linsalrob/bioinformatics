#!/usr/bin/perl -w

# transpose a file containing tab seperated columns

use strict;
my $f=shift || die "file to transpose";

open(IN, $f) || die "Can't open $f";
my $linecount=0;
my $newfile;
my $maxr=0;
while (<IN>)
{
    my $rowcount=0;
    chomp;
    my @a=split /\t/;
    while (@a)
    {
        $newfile->[$rowcount]->[$linecount]=shift(@a);
        $rowcount++;
    }
    $linecount++;
    ($rowcount > $maxr) ? ($maxr=$rowcount) : 1;
}

for (my $i=0; $i<$maxr; $i++) {print join("\t", @{$newfile->[$i]}), "\n"}

