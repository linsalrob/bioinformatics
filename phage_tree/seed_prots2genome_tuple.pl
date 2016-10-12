#!/usr/bin/env perl 


use strict;
use Rob;
my $faf = shift || die "fasta file?";
my $fa=Rob->read_fasta($faf);
map {
	m/fig\|(\d+\.\d+)\./;
	print join("\t", $_, $1, length($fa->{$_})), "\n";
} sort {uc($a) cmp uc($b)} keys %$fa;

