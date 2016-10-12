#!/usr/bin/perl -w

# renumber the first column before we run xipe
use strict;
my $c=1;
my %count;

while (<>){
 my @a=split /\t/;
 unless ($count{$a[0]}) {$count{$a[0]}=$c; $c++}
 $a[0]=$count{$a[0]};
 print join "\t", @a;
}
