#!/usr/bin/perl -w

use strict;
use ParseTree;
my $parser = new ParseTree;
use Data::Dumper;

my $file = shift || die "Tree file (newick format) to parse";
open(IN, $file) || die "Can't open $file";
my $rawtree="";
while (<IN>) {chomp; $rawtree .= $_}
my $tree = $parser->parse($rawtree, 1);

#print STDERR Dumper($tree);
# exit;

my %dist;

parsenode($tree, 0);

map {print "$_\t$dist{$_}\n"} keys %dist;


sub parsenode {
	my ($node, $dist)=@_;
	my ($left, $right, $newdist)=@$node;

	if ($newdist) {$dist+=$newdist}
	
	if (ref($left) eq "ARRAY" && ref($right) eq "ARRAY") {
		parsenode($left, $dist);
		parsenode($right, $dist);
	}
	elsif ($left eq "node") {  
		$dist{$right}=$dist;
	}
	else {
		print STDERR "Uh oh : $left and $right\n";
	}
}





