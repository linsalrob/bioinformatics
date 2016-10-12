#!/usr/bin/perl -w

# convert a table that has ss/role/fig to numbers

use strict;


open(IN, $ARGV[0]) || die "Can't open $ARGV[0]";
my $ss;
while (<IN>) {
 chomp;
 my @a=split /\t/;
 push @{$ss->{$a[0]}->{$a[1]}}, $a[2];
}

my $sscount; my $c1=0;
my $rolecount; my $c2=0;
foreach my $sub (sort {$a cmp $b} keys %$ss) {
 unless (exists $sscount->{$sub}) {$c1++; $sscount->{$sub}=$c1}
 foreach my $role (sort {$a cmp $b} keys %{$ss->{$sub}}) {
  unless (exists $rolecount->{$role}) {$c2++; $rolecount->{$role}=$c2} 
  foreach my $peg (@{$ss->{$sub}->{$role}}) {
   print "$c1\t$c2\t$peg\n";
  }
 }
}

print STDERR "Last c1: $c1 Last c2: $c2\n";
foreach my $sub (sort {$a cmp $b} keys %$ss) {
 foreach my $role (sort {$a cmp $b} keys %{$ss->{$sub}}) {
  print "$sub\t", $sscount->{$sub}, "\t$role\t", $rolecount->{$role}, "\n";
 }
}

