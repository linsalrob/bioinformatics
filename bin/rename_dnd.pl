#!/usr/bin/perl -w

use strict;

my ($dnd, $corr)=@ARGV;
unless ($dnd && $corr) {die "$0 <dnd file> <corresp. file>"}

open(IN, $corr) || die "Cam't open $corr";
my %f;
while (<IN>) {
 chomp;
 my @a=split /\t/;
 $f{$a[0]}=$a[1];
}

open(IN, $dnd) || die "dnd?";
while (<IN>) {
 while (/(\w+)\:/) {
  my $b=$1; 
  if ($f{$b}) {s/$b\:/$f{$b}\@/}
  else {s/$b\:/$b\@/}
 }
 s/\@/\:/g;
 print;
}
