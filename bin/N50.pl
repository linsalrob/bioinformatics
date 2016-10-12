#!/usr/bin/perl -w

# given a list of tuples of [contig length, number of contigs] calculate the N50 (the largest entity E such that at least half of the total size of the entities is contained in entities larger than E).

use strict;
my $f=shift || die "File of [contig length, number of contigs]";
open(IN, $f) || die "Can't open $f";
my %length;
my $total;
while (<IN>) {
	chomp;
	my ($len, $n)=split /\t/;
	if ($length{$len}) {$length{$len}+=$n}
	else {$length{$len}=$n}
	$total += ($len * $n);
}

my @contigsizes = sort {$b <=> $a} keys %length;
my $currsize=0;
while ($currsize < int($total/2)) {
	my $l = shift @contigsizes;
	$currsize += ($l * $length{$l});
}
print "The N50 is $contigsizes[0]\n";
