#!/usr/bin/perl -w

use strict;
use lib '/home/redwards/bioinformatics/Modules';
use ParseTree;
use constant { LEFT => 0, RIGHT => 1};
use Data::Dumper;

my $pt = new ParseTree();

my $file = shift || die "tree file?";
my $treein;

open(IN, $file) || die "Can't open $file";
while (<IN>) {chomp; $treein.=$_}
close IN;
my $tree = $pt->parse($treein);
#print Dumper($tree);
#print "Done\n";

# find a node that has 100 children:
print  "The tree is ", ref($tree), "\n";
my $ct=0;
($tree, $ct) = &count_children($tree, 0);

print "There are $ct children\n";

sub count_children {
	my ($ref, $count)=@_;
	if (ref($ref) eq "ARRAY") {
		($ref->[0], $count) = &count_children($ref->[0], $count);
		($ref->[1], $count) = &count_children($ref->[1], $count);
		return ($ref, $count);
	}
	else {
		$count++;
		print "$ref\t$count\n";
		return ($ref, $count);
	}
}
