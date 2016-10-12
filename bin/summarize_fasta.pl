#!/usr/bin/perl -w

use strict;
use Rob;
my $rob=new Rob;

foreach my $f (@ARGV) {
	my $seq;
	if ($f =~ /fastq/) {
		my $fa = $rob->read_fastq($f);
		map {$seq->{$_} = $fa->{$_}->[0]} keys %$fa;
	}
	else {
		$seq = $rob->read_fasta($f);
	}

	my ($min, $max, $length, $n)=(["", 1e20], ["", 0], 0, 0);
	map {
		$n++;
		my $l = length($seq->{$_});
		($l > $max->[1]) ? ($max = [$_, $l]) : 1;
		($l < $min->[1]) ? ($min = [$_, $l]) : 1;
		$length+=$l;
	} keys %$seq;

	my $al = sprintf("%.3f", ($length/$n));
	if ($#ARGV > 0) {print "\nFile:    $f\n"}
	print <<EOF;
Total length: 	$length
Number of seqs: $n
Average:	$al

Longest:	$max->[0] ($max->[1] bp)
Shortest:	$min->[0] ($min->[1] bp) 
EOF
}

