#!/usr/bin/perl -w

# take the roles list of the subsystems and convert them to a series of numbers
# we only want to end up with three columns: 
# subsystem	role	peg
use strict;

my $f=shift || die "$0 <extracted roles and subsystems>";

my $one=1; my %one; 
my $two=1; my %two;
my $inss;

if (-e "numbered_ss.txt") {die "numbered_ss.txt already exists. Please move or rename"}
if (-e "ss_key.txt") {die "ss_key.txt already exists. Please move or rename"}

open (OUT, ">numbered_ss.txt") || die "Can't open numbered_ss.txt";

open(IN, $f) || die "Can't open $f";
my $warn=1;
while (<IN>) {
 chomp;
 my @list=split /\t/;
 if ($#list==4) {
  if ($warn) {
   print STDERR "There are 5 columns in $f. We are guessing that the first two are the classification, and we are using the last three. Please check this\n";
   undef $warn;
  }
  my $trash=shift @list;
  $trash=shift @list;
 }
 if ($warn && $#list!=2) {
  print STDERR "There are $#list columns in $f and we are not sure that we are doing this right!\n";
  undef $warn;
 }
 
 
 unless ($one{$list[0]}) {$one{$list[0]}=$one; $one++}
 unless ($two{$list[1]}) {$two{$list[1]}=$two; $two++}
 print OUT $one{$list[0]}, "\t", $two{$list[1]}, "\t", $list[2], "\n";
 $inss->{$list[0]}->{$list[1]}=1;
}

close OUT;
open (OUT, ">ss_key.txt") || die "Can't open ss_key.txt";
foreach my $k (sort {$a cmp $b} keys %$inss) {
 foreach my $y (sort {$a cmp $b} keys %{$inss->{$k}}) {
  print OUT join "\t", $k, $one{$k}, $y, $two{$y};
  print OUT "\n";
 }
}
