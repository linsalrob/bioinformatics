#!/usr/bin/perl -w

# convert data to to the three column xipe format

use strict;
my $keyf=shift || die "$0 <key file> <xipe input file>";
my $xf  =shift || die "$0 <key file> <xipe input file>"; 

open(IN, $keyf) || die "Can't open $keyf";
my $max=0;
my %count;
while (<IN>)
{
    chomp;
    my @a=split /\t/;
    $count{$a[0]}=$a[1];
    ($a[1] > $max) ? ($max = $a[1]) : 1;
}
$max++;

open(IN, $xf) || die "can't open $xf";
open(OUT, ">$xf.xipe") || die "Can't write to $xf.xipe";
while (<IN>)
{
    chomp;
    next unless ($_);
    my @a=split /\t/;
    unless ($a[0]) {print STDERR "No a0 in $_\n"; next}
    unless ($count{$a[0]}) {$count{$a[0]}=$max; $max++}
    print OUT join("\t", $count{$a[0]}, "0", $a[1]), "\n";
}
close IN;
close OUT;

open(OUT, ">$keyf") || die "Can't open $keyf";
print OUT map {"$_\t$count{$_}\n"} keys %count;
close OUT;

