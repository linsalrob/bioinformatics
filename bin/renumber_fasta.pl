#!/usr/bin/perl 
#
# renumber fasta sequences and write out an id map

use strict;

my $usage =<<EOF;
$0 <fasta file> <output file> <id map file> <number to start at (optional, default=0)>
EOF

my $f=shift || die $usage;
my $o=shift || die $usage;
my $i=shift || die $usage;
my $c=shift || 0;

open(IN, "$f") || die $!;
open(OUT, ">$o") || die $!;
open(ID, ">$i") || die $!;
while(<IN>) {
	if (s/^>//) {
		print ID "$c\t$_";
		print OUT ">", $c++, "\n";
	}
	else {
		print OUT;
	}
}
close IN;
close OUT;
close ID;
