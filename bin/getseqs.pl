#!/usr/bin/perl -w

# get some sequences from a fasta file

use strict;
my $inf=shift || die "file of sequences we need?";
open(IN, $inf) || die "Can't open $inf";
my %need;
while (<IN>) {chomp; $need{$_}=1}
close IN;

my $fa=shift || die "fasta file?";
open(IN, $fa) || die "Can't iope  $fa";
my $print;
while (<IN>) {
 if (/^>(\S+)/) { 
  my $t=$1;
  if ($need{$t}) {$print = 1}
  else {undef $print}
 }
 print if ($print);
}

