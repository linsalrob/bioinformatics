#!/usr/bin/perl -w

use strict;


open(IN, $ARGV[0]) || die "Can't open $ARGV[0]";
while (<IN>) {
	chomp;
	my ($id1,$id2,$iden,$ali_ln,$mismatches,$gaps,$b1,$e1,$b2,$e2,$psc,$bsc,$ln1,$ln2) = split /\t/;
	print join("\t", $id2,$id1,$iden,$ali_ln,$mismatches,$gaps,$b2,$e2,$b1,$e1,$psc,$bsc), "\n";
}

