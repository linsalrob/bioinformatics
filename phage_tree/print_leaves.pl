#!/usr/bin/perl -w

use strict;
use lib '/home/redwards/bioinformatics/Modules';
use ParseTree;
use Getopt::Std;
my %opts;
getopts('f:', \%opts);
unless ($opts{f}) {
	die <<EOF;
$0
-f tree file

Print out all the leaves on the tree

EOF
}


my $pt = new ParseTree();

my $treein;

open(IN, $opts{f}) || die "Can't open $opts{f}";
while (<IN>) {chomp; $treein.=$_}
close IN;
my $tree = $pt->parse($treein);

printnodes($tree);


sub printnodes {
	my ($node)=@_;
	if (ref($node) eq "ARRAY") {
		printnodes($node->[0]);
		printnodes($node->[1]);
	}
	else {
		print "$node\n";
	}
}



