#!/usr/bin/perl -w

use strict;
my $genomefile=shift || die "phage database file";
my %genome;
open(IN, $genomefile) || die "CAn't open $genomefile";
while (<IN>) {
 my @a=split /\t/;
 $genome{$a[1]}=$a[0];
}

die "other file(s) to correct?" unless ($ARGV[0]);

foreach my $f (@ARGV) {
 open(IN, $f) || die "Can't opnbe $f";
 open(OUT, ">/home2/rob/temp.txt") || die "Can't open /home2/rob/temp.txt";
 while (<IN>)
 {
  my @a=split /\t/;
  if ($genome{$a[6]}) {$a[7]=$genome{$a[6]}}
  else {print STDERR "No genome for $a[6]\n"}
  print OUT join("\t", @a);
 }
 
 `mv $f $f.old`;
 `cp -i /home2/rob/temp.txt $f`;
}
