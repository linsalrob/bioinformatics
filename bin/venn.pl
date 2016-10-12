#!/usr/bin/perl -w

# Generate some data for a venn diagram type picture

use strict;

my $count;
my %pos; my $i=1;
foreach my $f (@ARGV) {
 $pos{$f}=$i; $i++;
 open(IN, $f) || die "Can't open $f";
 while (<IN>) {
  my @a=split /\t/;
  $count->{$a[1]}->{$f}++;
 }
}

# set things to zero unless they exist
#foreach my $ss (keys %$count) {
# map { $count->{$ss}->{$_}=0 unless ($count->{$ss}->{$_}) } keys %pos;
#}

my $cross; my $byf;
foreach my $f (keys %pos) {$cross->{$f}=0}
foreach my $ss (keys %$count) {
 my $list=join " and ", sort {$pos{$a} <=> $pos{$b}} keys %{$count->{$ss}};
 $cross->{$list}++;
 push @{$byf->{$list}}, $ss;
}

print map {"$_\t" . $cross->{$_} . "\n"} sort {$a cmp $b} keys %$cross;

print "\n\n\n";
foreach my $f (keys %pos, "black and red", "red and black") {
 print "$f: ", join "\t", @{$byf->{$f}}, "\n" if ($byf->{$f});
}
