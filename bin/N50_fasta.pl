#!/usr/bin/perl -w

# Calculate the N50 for a fasta file.
# The N50 is the largest entity E such that at least half of the total size of the entities is contained in entities larger than E.

use strict;
use Rob;

my $faf=shift || die "fasta file?";
my $fa = Rob->read_fasta($faf);
my %length;
my $total;
foreach my $k (keys %$fa) {
	my $len = length($fa->{$k});
	$length{$len}++;
	$total += $len;
}

my @contigsizes = sort {$b <=> $a} keys %length;
my $currsize=0;
while ($currsize < int($total/2)) {
	my $l = shift @contigsizes;
	$currsize += ($l * $length{$l});
}
print "The N50 is $contigsizes[0]\n";
