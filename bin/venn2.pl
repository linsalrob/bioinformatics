#!/usr/bin/perl -w

# create some more venn diagram data. But better.

use strict;

my @files=@ARGV;
die "$0 <list of files>" unless ($ARGV[0]);

my $count; my $revcount;
foreach my $f (@files) {
 open(IN, $f) || die "Can't open $f";
 while (<IN>) {
  my @a=split /\t/;
  $count->{$a[1]}->{$f}++;
  $revcount->{$f}->{$a[1]}++;
 }
}


while (@files) {
 my $t=shift(@files);
 print "There are ", scalar(keys %{$revcount->{$t}}), " subsystems in $t\n";
 foreach my $f (@files) {
  my $compare=0;
  foreach my $ss (keys %{$revcount->{$t}}) {if ($count->{$ss}->{$f}) {$compare++}}
  print "There are $compare subsystems in $t and $f\n";
 }
}


