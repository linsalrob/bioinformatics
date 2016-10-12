#!/usr/bin/perl -w


# get information from some files


#/home/rob/bin/phredphrap/bin/phrap pedulla.fasta.screen -minmatch 14 -maxmatch 30 -bandwidth 14 -minscore 30 -vector_bound 99 -bypasslevel 1 -maxgap 30 -repeat_stringency 0.95 -new_ace -view

use strict;

my @words=(qw[minmatch maxmatch bandwidth minscore vector_bound bypasslevel maxgap repeat_stringency]);

my %min; my %max;
foreach my $wd (@words) {$min{$wd}=[10000, 'unknown']; $max{$wd}=[0, 'unknown']}

my $pwd=`pwd`; chomp($pwd);
print STDERR "Working on $pwd\n";
opendir(DIR, $pwd) || die "Can't open .";
while (my $f=readdir(DIR)) {
 next unless ($f =~ /.phrap.out/);
 my $h= join " ", `head -n 1 $f`;
 foreach my $wd (@words) {
  $h =~ /$wd\s+(\S+)/;
  my $val=$1;
  next unless ($val);
  unless ($val =~ /^\d+$/ || $val =~ /^\d+\.\d+$/) {print STDERR "$val is not a number\n"; next}
  if ($val < $min{$wd}->[0]) {$min{$wd}=[$val, $f]}
  if ($val > $max{$wd}->[0]) {$max{$wd}=[$val, $f]}
 }
}

foreach my $wd (@words) {
 print join "\t", $wd, @{$min{$wd}}, @{$max{$wd}}, "\n";
}
