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
my $count=0;
&count_children($tree);

print "There are $count children\n";

sub count_children {
	my ($ref)=@_;
	if (ref($ref) eq "ARRAY") {
		&count_children($ref->[0]);
		&count_children($ref->[1]);
	}
	else {
		$count++;
		print "$ref\t$count\n";
	}
}
