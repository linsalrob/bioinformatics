#!/usr/bin/perl -w

use strict;
use lib '/home/redwards/bioinformatics/Modules';
use ParseTree;
use constant { LEFT => 0, RIGHT => 1};
use Data::Dumper;
use Getopt::Std;
my %opts;
getopts('f:x:m:', \%opts);
unless ($opts{f} && $opts{m} && $opts{x}) {
	die <<EOF;
$0
-f tree file
-x maximum number of children in subtree
-m minimum number of children in subtree

We have to specify a range, because for example a subtree may only have 2 nodes, and be branching  at the root of a tree with ~400 nodes

EOF
}


my $pt = new ParseTree();

my $treein;

open(IN, $opts{f}) || die "Can't open $opts{f}";
while (<IN>) {chomp; $treein.=$_}
close IN;
my $tree = $pt->parse($treein);
# add the number of children to each node
$tree = $pt->count_nodes($tree);

# now we just need to iterate through and find a node with 100 children, and print out that list of children!!
# except we want the first node where n < 100. Bugger
my $subtree = find_n_children($tree, ["root", "", $opts{x}+1], 100);


sub find_n_children {
	my ($node, $parnt)=@_;
	return if (ref($node) ne "ARRAY");
	if ($parnt->[2] <= $opts{x} && $node->[2] >= $opts{m}) {
		&printnodes($parnt);	
		print "\n";
	}
	else {
		find_n_children($node->[0], $node);
		find_n_children($node->[1], $node);
	}
}

sub printnodes {
	my ($node)=@_;
	if (ref($node) eq "ARRAY") {
		printnodes($node->[0]);
		printnodes($node->[1]);
	}
	else {
		print "$node; ";
	}
}



